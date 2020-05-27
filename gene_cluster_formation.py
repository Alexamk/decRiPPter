# License: GNU Affero General Public License v3 or later

import pickle
import os
import sys
import argparse
import numpy as np
import time
import copy
import shutil

from collections import OrderedDict

PYTHON_VERSION = sys.version_info[0]
if PYTHON_VERSION == 2:
    import ConfigParser as configparser
elif PYTHON_VERSION == 3:
    import configparser


#import antismash_wrapper_copy as antismash
from lib.lib import store_pickle,load_pickle,flip_dict,Container,Min,Max,Length,Contains,IsIn,DictValue,NOT,OR,AND,read_template
from lib.Genes import Operon, OperonCollection
from lib.log import return_logger, setup_loggers

import lib.parse_antismash as antismash
import lib.hmmer_operons as hmm
import lib.group_operons_precursors as group_precursors
import lib.group_operons_domains as group_jaccard
import lib.group_operons_EC as group_EC
import lib.group_operons_dom_prec as group_pj
import lib.SVG_operon as arrow
import lib.write_html as html
#import lib.RRE as RRE
import lib.write_operons as wgbk
import lib.write_operons_js as woj

#******* Operon formation functions
    
def find_all_operons(genome_dict, distance = 500):
    '''Finds all operons for the parsed hits with the simple method'''
    logger.debug("Function find_all_operons")
    all_operons = {}
    counter = 1
    operon_basename = 'genecluster_'
    for gene in genome_dict.itergenes(SVM_hit=True):
        operon_name = operon_basename + str(counter)
        genes = find_operon(gene,distance,operon_name)
        logger.debug('Gene cluster found with %i genes for %s' %(len(genes), gene))
        all_operons[operon_name] = genes
        counter += 1
    operon_dict = OperonCollection(all_operons,active=True)
    return operon_dict

def print_flanking_genes(gene,nr_flank=5):
    # Print the "nr_flank" flanking genes of the given gene
    flanks = [gene]
    logger.debug('left flank')
    for dir in ['left','right']:
        current_gene = gene
        for i in range(nr_flank):
            flank = getattr(current_gene,'%s_flank' %dir)
            if flank == None:
                break
            if dir == 'left':
                flanks.insert(0,flank)
            else:
                flanks.append(flank)
            current_gene = flank
        current_gene = gene
    for flank in flanks:
        flank.print_short()
        if flank == gene:
            logger.debug('right flank')
        
def find_operon(gene,max_dist,operon_name):
    genes = {gene.name:gene}
    logger.debug('Finding gene cluster simple method for gene %s' %(gene))
    for dir in ['left','right']:
        while True:
            flank = getattr(gene,'%s_flank' %dir)
            dist = getattr(gene,'%s_dist' %dir)
            if flank != None and dist <= max_dist and flank.strand == gene.strand:
                genes[flank.name] = flank
                gene = flank
            else:
                break
    operon = Operon(genes, set_genes=True, name=operon_name, active=True)
    return operon

def form_operons_island(genome_dict,settings):
    all_operons = {}
    operon_counter = 1
    operon_basename = 'operon_'
    all_islands = {}
    for hit in genome_dict.itergenes(SVM_hit=True):
        logger.debug('Finding gene cluster island method for gene %s' %(gene))
        operon_name = operon_basename+str(operon_counter)
        island_dict = find_islands(hit,island_distance=settings['island_dist'],gene_distance=settings['island_gene_dist'])
        base_island = island_dict.find_operon_by_gene(hit) # This island contains the hit 
        if len(base_island.genes) == 1:
            #SVM hit alone in an island - add it to the closest island
            island_dict = fuse_hit_to_closest(base_island,island_dict)
            base_island = island_dict.find_operon_by_gene(hit)        
        if len(island_dict) == 1:
            # One island --> base island is the operon
            operon = base_island
            operon.name = operon_name
        else:
            operon = fuse_islands(genome_dict,base_island,settings['island_cog_cutoff'],settings['std_factor'],operon_name)
            operon.name = operon_name
        all_operons[operon_name] = operon
        all_islands[operon_name] = island_dict
        operon_counter += 1
        
    operons_dict = OperonCollection(all_operons,active=True)
    logger.info('%i gene clusters formed' %(len(operons_dict)))
    return operons_dict, all_islands
    
def find_islands(hit,island_distance=500,gene_distance=50):
    counter = 1
    base_name = 'island_%i' %(counter)
    base_island = find_operon(hit,gene_distance,base_name)
    all_islands = {base_name:base_island}
    
    for dir in ['left','right']:
         edge_island = base_island
         while True:
            if dir == 'left':
                edge_gene = edge_island.genes[0]
            elif dir == 'right':
                edge_gene = edge_island.genes[-1]
            flank = getattr(edge_gene,'%s_flank' %dir)
            dist = getattr(edge_gene,'%s_dist' %dir)
            if flank != None and dist <= island_distance and flank.strand == edge_gene.strand:
                counter += 1
                island_name = 'island_%i' %(counter)
                logger.debug('Forming island %s' %island_name)
                edge_island = find_operon(flank,gene_distance,island_name)
                all_islands[island_name] = edge_island
            
            else:
                if dir == 'left':
                    logger.debug('Changing direction')
                elif dir == 'right':
                    logger.debug('Breaking island search')
                break
    island_dict = OperonCollection(all_islands)
    logger.debug('find_islands function found %i island(s).' %len(island_dict))
    return island_dict
    
def fuse_hit_to_closest(base_island,island_collection):
    logger.debug('Fusing hit to closest island')
    # First find out the closest island
    if not (base_island.left_flank or base_island.right_flank):
        logger.debug('No islands to fuse the hit to.')        
        return island_collection
    elif (not base_island.left_flank) and base_island.right_flank:
        dir = 'right'
        fuse_island = base_island.right_flank
    elif (not base_island.right_flank) and base_island.left_flank:
        fuse_island = base_island.left_flank
        dir = 'left'
    elif base_island.right_dist < base_island.left_dist:
        dir = 'right'
        fuse_island = base_island.right_flank
    else:
        fuse_island = base_island.left_flank
        dir = 'left'
    logger.debug('Fusing base island to %s flank (%s)' %(dir,fuse_island.name))
    new_operon = Operon.fuse(base_island,fuse_island)
    new_collection = dict([(operon.name,operon) for operon in island_collection if (operon != base_island and operon != fuse_island)])
    new_collection[new_operon.name] = new_operon
    return OperonCollection(new_collection)
        
def fuse_islands(genome_dict,start_island,COG_cutoff,std_factor,operon_name,max_outliers=2):
    # Islands are fused when their COG averages are within the cutoffs
    # A certain number of outliers may be added, up to max_outliers, if at least one more island is added behind the outlier without going over the cutoff
    logger.debug('Fusing islands to gene cluster')
    current_outliers = 0
    start_island.outlier=False
    current_islands = [start_island]
    edge_islands = {'left':start_island,'right':start_island}
    current_COGs = start_island.COG
    consecutive_failed_attempts = 0
    dirs = ['left','right']
    dir_index = 0
    while consecutive_failed_attempts < 2:
        dir = dirs[dir_index]
        dir_index = not dir_index
        logger.debug('Direction: %s. Consecutive failed attempts: %i' %(dir,consecutive_failed_attempts))
        while True:
            edge_island = edge_islands[dir]
            flank_island = getattr(edge_island, '%s_flank' %(dir))
            if len(current_COGs) < 2:
                if hasattr(start_island.genome,'COG_global_std'):
                    current_std = start_island.genome.COG_global_std
                else:
                    current_std = 0.05
            else:
                current_std = np.std(current_COGs)
            logger.debug('Current nr of islands: %i. Current COG score: %.3f +- %.3f' %(len(current_islands),np.average(current_COGs),current_std))
            if flank_island != None:
                logger.debug('Flanking island: %s' %(flank_island))
                if not flank_island.COG_avg:
                    # Always add for now, although this should not occur in normal runs
                    edge_islands[dir] = flank_island
                    current_islands.append(flank_island)
                    logger.debug('Adding island without COG score')
                else:
                    COG_difference = (abs(np.average(current_COGs) - flank_island.COG_avg) - current_std - flank_island.COG_std)
                    if COG_difference < COG_cutoff:
                        if edge_island.outlier:
                            # Now also add the outlier
                            current_outliers += len(edge_island.genes)
                            logger.debug('Adding outlier. Remaining outliers: %s.' %(max_outliers-current_outliers))
                            current_islands.append(edge_island)
                        flank_island.outlier=False
                        edge_islands[dir] = flank_island
                        current_islands.append(flank_island)
                        current_COGs.extend(flank_island.COG)
                        consecutive_failed_attempts = 0
                        logger.debug('Island added! COG difference: %.2f' %(COG_difference))
                    elif current_outliers < max_outliers:
                        # Set this island as the edge gene, but don't add it yet
                        flank_island.outlier=True
                        edge_islands[dir] = flank_island
                        logger.debug('Remembering outlier island. COG difference: %.2f.' %(COG_difference))
                    else:
                        consecutive_failed_attempts += 1
                        logger.debug('Island was below cutoff! COG_difference: %.2f' %(COG_difference))
                        break
            else:
                consecutive_failed_attempts += 1
                logger.debug('No more islands in this direction! Changing direction')
                break
    operon = Operon.fuse(*current_islands)
    logger.debug('fuse_islands function fused %i islands: %s' %(len(current_islands),', '.join([i.name for i in current_islands])))
    return operon
    
    
#***** Fusing the operons together
#      Also assign their final names here

def fuse_operons(all_operons,max_dist,namebase, rename=True):
    # Even faster version that sorts the first time more efficiently, by only going over the operons once
    
    logger.debug('Fusing overlapping gene clusters')
    counter = 0
    final_operons = dict([(operon.name,operon) for operon in all_operons])
    # Presort operons per scaffold,strand to prevent repetitive sorting
    all_scaffolds = set()
    for operon in all_operons:
        all_scaffolds.add(operon.scaffold.name)
    
    sorted_dict = {}
    for operon in all_operons.itersubset(sorted=True):
        key = (operon.scaffold.name,operon.strand)
        if key not in sorted_dict:
            sorted_dict[key] = []
        sorted_dict[key].append(operon)
    
    for operon in all_operons.itersubset(sorted=True):
        if operon.name not in final_operons:
            continue # To prevent operons already fused to be added again
        original_operon = operon
        found_same = False
        for operon2 in sorted_dict[(operon.scaffold.name,operon.strand)]:
            if operon2.name not in final_operons:
                continue
            if operon2 == original_operon:
                found_same = True
                continue
            if not found_same:
                continue
            if ((operon.start >= operon2.start and operon.start - operon2.end <= max_dist) or
                (operon2.start >= operon.start and operon2.start - operon.end <= max_dist)):
                # Fuse them 
                logger.debug('Fusing %s to %s' %(operon.name,operon2.name))
                logger.debug('%s coords: %i, %i, %s, %s' %(operon.name,operon.start,operon.end,operon.strand,operon.scaffold.name))
                logger.debug('%s coords: %i, %i, %s, %s' %(operon2.name,operon2.start,operon2.end,operon.strand,operon.scaffold.name))
                new_operon = Operon.fuse(operon,operon2)
                # Update the output operons
                final_operons[operon.name] = new_operon
                final_operons.pop(operon2.name)
                # Continue comparing operons using the new boundaries
                operon = new_operon
                counter += 1
            else:
                # Don't try the rest of the operons
                break
    logger.info('%i operons were fused' %counter)
    if rename:# Rename all the operons
        renamed_operons = {}
        counter = 1
        for operon_name in final_operons:
            operon = final_operons[operon_name]
            new_name = '%s_%i' %(namebase,counter)
            operon.name = new_name
            renamed_operons[new_name] = operon
            counter += 1
        all_operons_out = OperonCollection(renamed_operons,active=True)
    else:
        all_operons_out = OperonCollection(final_operons, active=True)
    return all_operons_out,sorted_dict

#***** Extending the operons

def extend_operons(operons,nr_genes_to_add = 5):
    # This function will add extra genes to the genes to each side.
    # From this point onward, the genes currently in the operon will be called 
    # the core genes, the added ones the extra_genes. All the genes will be
    # grouped under the genes list.
    logger.debug('Extending operons with %i genes on either side' %nr_genes_to_add)
    for operon in operons.itersubset():
        genes = [i for i in operon.itergenes(sorted=True)]
        for _ in range(nr_genes_to_add):
            left_gene = genes[0].left_flank
            right_gene = genes[-1].right_flank
            if left_gene != None:
                genes.insert(0,left_gene)
                operon[left_gene.name] = left_gene
            if right_gene != None:
                genes.append(right_gene)
                operon[right_gene.name] = right_gene
        operon.extra_start = genes[0].start
        operon.extra_end = genes[-1].end
    # Now mark all genes in a core part of an operon as such
    # (This still means a gene can be part of multiple operons,
    #  core in one, not core in others)
    for gene in operons.itergenes(core=True):
        gene.is_core = True
    for gene in operons.itergenes():
        if not hasattr(gene,'is_core'):
            gene.is_core = False
    
#****** Setting operon domains, filtering them and scoring them

def assign_gene_types(operons, mibig_domaindict, domain_types):
    for itr in operons, mibig_domaindict:
        for gene in itr.itergenes():
            genetypes = []
            for key in domain_types:
                for domain in gene.domains:
                    if domain in domain_types[key]:
                        genetypes.append(key)
                        break
            if genetypes == []:
                genetypes = ['other']
            gene.genetypes=genetypes

def set_operon_domains(operons, mibig_dict, domain_categories, genetypes):
    # The domains associated with the operon here are those of the core genes
    
    def find_domains(operon,domains,name,all_domains=set(),**kwargs):
        domains_found = set()
        for gene in operon.itergenes(**kwargs):
            for domain in gene.domains:
                all_domains.add(domain)
                if domain in domains:
                    domains_found.add(domain)
        setattr(operon,name,sorted(list(domains_found)))
        return(all_domains)
    
    def set_genes(operon,genetype,name,**kwargs):
        genes = []
        for gene in operon.itergenes(**kwargs):
            if genetype in gene.genetypes:
                genes.append(gene)
        setattr(operon,name,genes)
    
    counter = 0
    for itr in operons,mibig_dict:
        for operon in itr:
            counter += 1
            if counter % 1000 == 0:
                logger.debug('At operon %i' %counter)
            all_domains = set()
            for category in domain_categories:
                domains = domain_categories[category]
                if category == 'precursor':
                    all_domains = find_domains(operon,domains,category,all_domains,core=True)
                else:
                    all_domains = find_domains(operon,domains,category,all_domains,core=True,SVM_hit=False)
                _ = find_domains(operon,domains,'%s_all' %category)
            setattr(operon,'all_domains',sorted(list(all_domains)))
            for genetype in genetypes:
                if genetype == 'precursor':
                    set_genes(operon,genetype,'genes_%s' %genetype,core=True)
                else:
                    set_genes(operon,genetype,'genes_%s' %genetype,core=True,SVM_hit=False)
                set_genes(operon,genetype,'genes_%s_all' %genetype)
            
            
def set_operon_unique_domains(operons,mibig_dict):
    for itr in operons,mibig_dict:
        for operon in itr:
            operon.unique_domains = set(operon.all_domains)

def score_operons(operons,score_order,settings):
    for operon in operons:
        cog = operon.COG_avg if getattr(operon,'COG_avg') else 1
        precursors = [gene for gene in operon.itergenes(core=True,SVM_hit=True)]
        genes_added = set(precursors)
        base_score = 0
        for key in score_order:
            if hasattr(operon,'genes_%s' %key) and 'score_%s' %key in settings:
                genes = [gene for gene in getattr(operon,'genes_%s' %key) if gene not in genes_added]
                for gene in genes:
                    genes_added.add(gene)
                    base_score += settings['score_%s' %key]
                if len(genes) > 1 and 'bonus_%s' %key in settings:
                    base_score += settings['bonus_%s' %key]
        
        score =  base_score* len(precursors) ** settings['weight_precursors'] / \
                 (cog ** settings['weight_cog'])
        operon.score = score
    
#***** Reading/Writing functions

def initialize(settings,name):
    path = settings['paths']['pickles']
    EC_classes = load_pickle(os.path.join(path, "EC_classes.pkl"))
    
    mibig_domaindict = load_pickle(os.path.join(path, 'mibig_domains_1.4.pkl'))

    all_domains = get_domain_files(settings)
    domain2ec = load_pickle(os.path.join(path, 'domain2ec.pkl'))
    domain_descr = load_pickle(os.path.join(path, 'domain_descr.pkl'))

    genome_dict = load_pickle(os.path.join(path, 'genome_dict_%s.pkl' %name))
    logger.debug('Setting quick references for genome_dict')
    genome_dict.set_quick()
    logger.debug('Switching genome_dict to active')  
    genome_dict.switch_active()
    return(EC_classes,all_domains,domain_descr,domain2ec,mibig_domaindict,genome_dict)

def get_domain_files(settings):
    domain_path = settings['paths']['domain']
    domains_per_category = {}
    for f in os.listdir(domain_path):
        read_first = False
        with open(os.path.join(domain_path,f)) as infile:
            for l in infile:
                if l.startswith('#') or l == '\n':
                    continue
                if not read_first:
                    category = l.strip()
                    domains_per_category[category] = {}
                    read_first = True
                    continue
                l = l.strip()
                if '\t' in l:
                    domain,descr = l.split('\t')
                else:
                    domain,_,descr = l.partition(' ')
                    while descr.startswith(' '):
                        descr = descr[1:]
                domains_per_category[category][domain] = descr
    return domains_per_category

def load_operons(path):
    # Loads the most recently stored operon file
    all_operon_files = [i for i in os.listdir(path) if 'all_operons_fused' in i]
    logger.debug('Load gene_clusters files found: %s' %(','.join(all_operon_files)))
    if len(all_operon_files) == 1:
        res = load_pickle(os.path.join(path,all_operon_files[0]))
        return res
    else:
        all_times = [i.split('_')[-1] for i in all_operon_files]
        all_times = [i.split('.') for i in all_times]
        at_sorted = sorted(all_times,key = lambda x: (x[2],x[1],x[0]), reverse=True)
        most_recent = at_sorted[0]
        index = all_times.index(most_recent)
        f = all_operon_files[index]
        Logger.debug('Loading file: %s' %f)
        operons = load_pickle(path + f)
        return(operons)

def store_operons(operons,genome_dict,path,name):
    if operons.active:
        start_active = True
        operons.switch_inactive()
    else:
        start_active = False
    y,m,d = time.gmtime()[0:3]
    time_operons = '%s.%s.%s' %(d,m,y)
    pfile = open(os.path.join(path, '%s_%s.pkl' %(name,time_operons)),'w')
    pickle.dump(operons,pfile)
    pfile.close()  
    if start_active:
        operons.switch_active(genome_dict)

def set_operon_paths(operons, outpath):
    for operon in operons:
        operon_folder = os.path.join(outpath, operon.name)
        if not os.path.isdir(operon_folder):
            os.mkdir(operon_folder)
        operon.path = operon_folder
        operon.fasta_path = os.path.join(operon_folder,operon.name + '.fasta')
        operon.genbank_path = os.path.join(operon_folder,operon.name + '.gbk')
    
def set_collection_paths(collections_by_type,settings):
    out_path = settings['paths']['out']
    entry_path = os.path.join(out_path,'Subgroups')
    if not os.path.isdir(entry_path):
        os.mkdir(entry_path)
    for coll_type in collections_by_type:
        collections = collections_by_type[coll_type]
        type_path = os.path.join(entry_path,coll_type)
        if not os.path.isdir(type_path):
            os.mkdir(type_path)
        for name,coll in collections.items():
            coll.path = os.path.join(type_path,coll.name + '.html')
            
def make_collections_per_genomes(operons,genome_dict):
    collections_per_genome = {}
    for genome in genome_dict:
        operons_genome = {}
        for operon in operons.itersubset(genome=genome):
            operons_genome[operon.name] = operon
        operon_coll = OperonCollection(operons_genome,name='operons_genome_%s' %genome.name, descr='Grouped by genome', \
                                       collection_type='genome',realname=genome.descr,active=True)
        collections_per_genome[genome.name] = operon_coll
    return collections_per_genome
    
def prep_collections(collections_by_type):
    # Get all domains that are present in at least 50% of operons
    for collections in collections_by_type.values():
        for coll in collections.values():
            logger.debug('Setting common domains for collection %s' %coll)
            coll.set_common_domains()
            if not hasattr(coll,'mibig'):
                coll.calc_avg_COG()

def set_pairs(mibig_dict,keywords):
    for bgc in mibig_dict:
        pairs = Container()
        bgc.pairs = pairs
        for k in keywords:
            setattr(bgc.pairs,k,{})

def get_pair_data(operon,operons,keywords,extra_operons=False):
    pair_data = {}
    for keyword in keywords:
        pairs = getattr(operon.pairs,keyword)
        for operon2 in pairs:
            # Check that the operon pair also passes the filter
            if operon2 not in operons and not(extra_operons and operon2 in extra_operons):
                continue
            if operon2 not in pair_data:
                pair_data[operon2] = {}
            pair_data[operon2][keyword] = pairs[operon2]

    return pair_data
    
def write_edge_data(operon,COG_text,pair_data,keywords,handle):
    antismash_ripps = ['bacteriocin','cyanobactin','lantipeptide','lanthipeptide',\
                   'lassopeptide','linaridin','thiopeptide','sactipeptide',\
                   'proteusin','glycocin','bottromycin','microcin']
                   
    for operon2 in pair_data:
        handle.write('%s\t%s' %(operon.name,operon2))
        text = ''
        pair_types = []
        for keyword in keywords:
            text += '\t'
            if keyword in pair_data[operon2]:
                edge = pair_data[operon2][keyword]
                if keyword == 'precursor':
                    text += str(round(edge.weight/100.0,3))
                else:
                    text += str(round(edge.weight,2))
                pair_types.append(keyword)
        # Determine the edge type
        if 'precursor' in pair_types and 'jaccard' in pair_types:
            text = '\tprecursor_domains' + text
        elif 'precursor' in pair_types:
            text = '\tprecursor' + text
        elif 'jaccard' in pair_types:
            text = '\tjaccard' + text
        # Write some general properties of the first operon
        if hasattr(operon,'score'):
            text += '\t%s' %(str(round(operon.score,2)))
        else:
            text += '\t'
        text += '\t%s' %(COG_text)
        text += '\t%i' %(len(operon.genes_kripp))
        if settings['antismash_parse']:
            # First crude information (all gene cluster types found)
            # Then filtered information (only RiPP gene clusters; or "other" if none are found)
            # Only consider gene clusters when their overlap is larger than 0.5
            
            if hasattr(operon,'mibig'):
                # Get the subtypes; rewrite a bit so they resemble those from antismash
                # TODO: Fix the annotation in the mibig domaindict so that they just match (or add extra description)
                if operon.subtype == 'lasso peptide':
                    subtype = 'lassopeptide'
                elif operon.subtype == 'Two-component bacteriocin':
                    subtype = 'bacteriocin'
                elif ' / ' in operon.subtype:
                    subtype = operon.subtype.replace(' / ','-')
                else:
                    subtype = operon.subtype
                text += '\t%s\t%s' %(operon.subtype,subtype)
            else:
                # First get all clusters found overlapping
                # All clusters are fused together, even if they were part of two different gene clusters
                # This should usually make no difference, but sometimes will in the cases where two antismash
                # gene clusters are partially overlapping with the operon, without overlapping with each other
                cluster_types = []
                for cluster,cluster_data in operon.overlaps_antismash.items():
                    if cluster_data['overlap_type'] == 'total' or \
                       (cluster_data['overlap_type'] == 'partial' and cluster_data['overlap_fraction'] > settings['antismash_minimal_overlap']):
                        cluster_types.extend(cluster_data['cluster_type'])
                text += '\t' + '-'.join(sorted(cluster_types))
                filtered_cluster_types = [i for i in cluster_types if i in antismash_ripps]
                # Check to see if there are any antismash RiPPs detected
                if len(cluster_types) == 0:
                    # No antismash detection
                    text += '\t'
                elif len(filtered_cluster_types) > 0:
                    text += '\t' + '-'.join(filtered_cluster_types)
                else:
                    text += '\tother'
        if settings['run_rre']:
            if hasattr(operon,'mibig'):
                text += '\t'
            else:
                text += '%s\t' %str(int(operon.nr_RRE))
        if hasattr(operon,'mibig'):
            handle.write(text + '\t1\n' )
        else:
            handle.write(text + '\t-1\n' )


def write_cytoscape(operons,mibig_dict,outfile,settings,**kwargs):
    
    keywords = ['precursor','jaccard','precursor_jaccard']
    
    with open(outfile,'w') as handle:
        handle.write('Operon1\tOperon2\tEdge_type\tPrecursor_similarity\tDomain_similarity\tAverage_weight_prec_dom\tscore\tCOG_average\tKnown_RiPP_genes')
        if settings['antismash_parse']:
            handle.write('\tantismash_full\tantismash_filtered')
        if settings['run_rre']:
            handle.write('\tRRE')
        handle.write('\tMIBiG\n')
        for operon in operons:
            if not hasattr(operon,'COG_avg') or operon.COG_avg == None:
                COG_text = ''
            else:
                COG_text = str(round(operon.COG_avg,2))
            
            # First gather all the data per operon pair
            pair_data = get_pair_data(operon,operons,keywords,mibig_dict)
            # Now write all the lines; one per operon pair
            write_edge_data(operon,COG_text,pair_data,keywords,handle)
            
        # Write the pairs between the mibig entries
        for operon in mibig_dict:
            COG_text = ''
            pair_data = get_pair_data(operon,mibig_dict,keywords,operons)
            write_edge_data(operon,COG_text,pair_data,keywords,handle)

def write_operons(operons,settings,outfile,**kwargs):
    with open(outfile,'w') as handle:
        header = ['Name','Genome','Scaffold','Start','End','Score','Number of genes','Number of precursors','COG average','COG stdev',\
                  'Domain groups','Precursor groups','E.C groups']
        if settings['antismash_parse']:
            header.append('antiSMASH overlaps')
        if settings['run_rre']:
            header += ['RRE']
        header += ['Number of biosynthetic genes','Number of transporter genes','Number of peptidase genes',\
                   'Number of regulator','Number of genes with known RiPP domains','Number of domains','Biosynthetic domains',\
                   'Transporter domains','Peptidase domains','Regulator domains','Known RiPP domains']
        handle.write('\t'.join(header) + '\n')
        for operon in operons.itersubset(**kwargs):
            data = [operon.name,operon.genome.name,operon.scaffold.name,str(operon.start),str(operon.end),str(round(operon.score,2)),\
                    str(len(operon.core_genes)),str(len(operon.precursors))]
            if operon.COG_avg == None :
                data += ['N\\A']
            else:
                data += [str(round(operon.COG_avg,3))]
            if operon.COG_std == None:
                data += ['N\\A']
            else:
                data += [str(round(operon.COG_std,3))]
            if hasattr(operon,'jaccard_groups'):
                data.append(operon.jaccard_groups)
            else:
                data.append('')
            for v in ['precursor_groups','EC_groups']:
                if hasattr(operon,v):
                    data.append(','.join(sorted(getattr(operon,v))))
                else:
                    data.append('')
            if settings['antismash_parse']:
                antismash_text = str(len(operon.overlaps_antismash))
                if len(operon.overlaps_antismash) > 0:
                    for cluster in operon.overlaps_antismash:
                        antismash_data = operon.overlaps_antismash[cluster]
                        antismash_text += ',%s' %('-'.join(antismash_data['cluster_type']))
                data.append(antismash_text)
            if settings['run_rre']:
                data.append(str(operons.nr_RRE))
            nr_genes_per_type = {}
            data += [str(len(operon.genes_biosyn)),str(len(operon.genes_transporter)),str(len(operon.genes_peptidase)),\
                     str(len(operon.genes_regulator)),str(len(operon.genes_kripp)),str(len(operon.all_domains)),','.join(operon.biosyn),
                     ','.join(operon.transporter),','.join(operon.peptidase),','.join(operon.regulator),','.join(operon.kripp)]
            handle.write('\t'.join(data) + '\n')
            
def write_subset(operons,group_collections,settings,filtername,headers,domain_descr,genome_dict,attributes,mibig_domaindict,copy_operons=True,**kwargs):
    # Writes a subset of the operons. Sets filtered on the operons
    settings_temp = copy.deepcopy(settings)
    out_path = os.path.join(settings['paths']['base'],'Output_%s' %filtername) + os.sep
    old_path = settings['paths']['out']
    logger.debug('Old path: %s' %old_path)
    logger.debug('New path: %s' %out_path)
    # Overwrite the paths with new values
    for key in settings_temp['paths']:
        p = settings_temp['paths'][key]
        if old_path in p:
            settings_temp['paths'][key] = p.replace(old_path,out_path)
    # Make some folders
    logger.info('Making folders')
    operon_path = settings_temp['paths']['operons']
    entry_path = os.path.join(out_path,'Subgroups')
    filtered_path = os.path.join(operon_path,'Filtered')
    unfiltered_path = os.path.join(operon_path,'Unfiltered')
    to_write = [out_path,operon_path,filtered_path,unfiltered_path,entry_path]
    for p in ['precursor','jaccard','EC','genome','precursor_jaccard']:
        to_write.append(os.path.join(entry_path, p))
    for p in to_write:
        if not os.path.isdir(p):
            os.mkdir(p)
    
    # Filter collections and set their paths
    logger.info('Filtering collections')
    group_collections_filtered = filter_collections_change_paths(group_collections,old_path,out_path,**kwargs)
    prep_collections(group_collections_filtered)
    # Set new paths and copy operon html files; remember old ones
    logger.info('Changing operons')
    old_paths = {}
    new_operons = {}
    operon_counter = 0
    for operon in operons.itersubset(**kwargs):
        operon_counter += 1
        new_operons[operon.name] = operon
        old_paths[operon.name] = {}
        for p in ['path','fasta_path','genbank_path']:
            old_op_path = getattr(operon,p)
            old_paths[operon.name][p] = old_op_path
            new_path = old_op_path.replace(old_path,out_path)
            setattr(operon,p,new_path)
            if p == 'path':
                new_folder = old_op_path.replace(old_path,out_path)
                if not os.path.isdir(new_folder):
                    os.mkdir(new_folder)
                old_op_path = os.path.join(old_op_path,operon.name + '.html')
                new_path = old_op_path.replace(old_path,out_path)
            if not os.path.isfile(new_path):
                if copy_operons:
                    logger.debug('Copying %s to %s' %(old_op_path,new_path))
                    shutil.copy(old_op_path,new_path)
                    
    logger.info('New operon collection with %i operons' %operon_counter)
    new_operons = OperonCollection(new_operons,prep=False,active=True)
    if not copy_operons:
        svg_text_dict = arrow.main(new_operons,all_domains,domain_descr,color_dict,genetypes,settings['paths']['html'],group_collections_filtered,settings_temp)
        logger.debug('Writing operon genbank files')
        wgbk.write_genbank(new_operons)
        logger.debug('Writing operon fasta files')
        wgbk.write_fasta_all(new_operons)
    
    # Write the operon js files
    settings_temp['paths']['operons_js_file'] = operons_js_file = os.path.join(out_path,'Operons_dict.js') 
    settings_temp['paths']['operon_group_js_file'] = operon_group_js_file = os.path.join(out_path,'Groups_dict.js')
    woj.write_operon_dict(new_operons,mibig_domaindict,attributes,operons_js_file)
    woj.write_group_dict(group_collections_filtered,operon_group_js_file)    
    
    # Write the new entry and index files    
    logger.info('Writing new index/entry files')
    t = html.main(settings_temp,new_operons,group_collections_filtered,headers,domain_descr,genome_dict)
    
    # Write a summary
    logger.info('Writing summary')
    summary_file = os.path.join(out_path,'%s_overview.tsv' %(name))
    write_operons(new_operons,settings_temp,summary_file,**kwargs)
    
    # Write a cytoscape file
    cytoscape_file = os.path.join(out_path,'%s_cytoscape.ssn' %name)
    write_cytoscape(new_operons,mibig_domaindict,cytoscape_file,settings_temp,**kwargs)
    
    # Reset the paths
    logger.info('Resetting paths')
    for operon in new_operons:
        op = old_paths[operon.name]
        for p in op:
            setattr(operon,p,op[p])
    return(new_operons,old_paths,svg_text_dict)
    
def filter_collections_change_paths(group_collections,old_outfolder,new_outfolder,**kwargs):
    out = {}
    for category in group_collections:
        operon_collections = group_collections[category]
        out[category] = {}
        for coll_name,operon_coll in operon_collections.items():
            new_operons = {}
            for operon in operon_coll.itersubset(**kwargs):
                # If any operon is given here, the collection can be passed on
                new_operons[operon.name] = operon
            # Also get the MIBiG operons
            for operon in operon_coll.itersubset(mibig=True):
                new_operons[operon.name] = operon
            if len([i for i in new_operons.values() if not hasattr(i,'mibig')]) > 1 or (len(new_operons) > 0 and category == 'genome'):
                new_coll = OperonCollection(new_operons,prep=False,collection_type=operon_coll.collection_type,name=operon_coll.name,\
                                            descr=operon_coll.descr)
                if hasattr(operon_coll,'realname'):
                    new_coll.realname = operon_coll.realname
                if hasattr(operon_coll,'muscle_path'):
                    new_coll.muscle_path = operon_coll.muscle_path
                # Change the path
                new_path = operon_coll.path.replace(old_outfolder,new_outfolder)
                new_coll.path = new_path
                out[category][coll_name] = new_coll
    return(out)
            
def read_config(config_file):
    config = configparser.ConfigParser()
    config.read(config_file)
    settings = {}
    for section in config.sections():
        items = config.items(section)
        for item,value in items:
            try:
#                print(section,item,value)
                value = int(value)
            except ValueError:
                try:
                    value = float(value)
                except ValueError:
                    try:
                        value = config.getboolean(section,item)
                    except ValueError:
                        if ',' in value:
                            value = value.split(',')
            settings[item] = value
    return settings
            
def parse_arguments():
    # Find the config file in the same folder as the 
    script_path = os.path.realpath(sys.argv[0])
    folder = os.path.dirname(script_path)
    config_file = os.path.join(folder,'config.ini')
    if not os.path.isfile(config_file):
        raise ValueError('Error! config.ini not found in base directory')
    settings = read_config(config_file)
                
    parser = argparse.ArgumentParser()
    
    parser.add_argument('name',metavar='PROJECT NAME',type=str,help='A name for the project to run')            
    parser.add_argument('-o','--outputfolder',help='path which contains the results folder',default='')
    parser.add_argument('-c','--cores', help='number of cores to be used on this run',type=int)
    parser.add_argument('-v', '--verbosity', help='Verbosity (0-3), higher means more output', choices = [0, 1, 2, 3], type=int, default=2)
    
    parser.add_argument('--gene_cluster_formation_method',help='Gene cluster formation method',choices=['simple','island'])
    parser.add_argument('--load_operons',help='load the most recent pickle file of created operons',default=False,action='store_true')
    
    parser.add_argument('--skip-hmm',help=argparse.SUPPRESS,default=False,action='store_true')
    
    args = parser.parse_args()
    for key,value in args.__dict__.items():
        # Update only if the args value is not None
        if value != None:
            settings[key] = value
    if settings['skip_hmm']:
        settings['run_hmm'] = False
        
    return args,settings
    
def check_paths(settings):
    base_path = settings['paths']['base']
    if not os.path.isdir(base_path):
        logger.error('No results folder found from previous run. Please run the genome_prep.py script first')
        raise ValueError('No results folder found from previous run. Please run the genome_prep.py script first')
    
    hmm_paths = settings['hmm_paths']
    for database_folder, path in hmm_paths:
        if not os.path.isfile(path):
#            database_name = os.path.basename(database_folder)
            logger.error('Path to hmm database is not found - please define it in the config.ini file')
            raise ValueError('Path to hmm database is not found - please define it in the config.ini file')
    
def build_operons(genome_dict, settings, rename=True):
    if not settings['load_operons']:
        logger.info("Forming gene clusters...")
        if settings['gene_cluster_formation_method'] == "island":
            all_operons, island_dict = form_operons_island(genome_dict,settings)
        elif settings['gene_cluster_formation_method'] == "simple":
            all_operons = find_all_operons(genome_dict,distance = settings['simple_dist'])
        operons, sorted_dict = fuse_operons(all_operons,50,operon_namebase, rename)        
        operons.clear_operon_flanks()
        if settings['store_operons']:
            logger.info('Storing gene clusters')
            store_operons(operons,genome_dict,pickle_path,name)
    else:
        operons = load_operons(pickle_path)
    return operons
        
def filter_operons(operons, settings):
    filter_dict = create_filter_criteria(settings)
    new_operons = {}
    for operon in operons.itersubset(**filter_dict):
        new_operons[operon.name] = operon
    filtered_operons = OperonCollection(new_operons,prep=False,active=True)
    return filtered_operons
        
def create_antismash_criteria(settings):
    return dict(overlaps_antismash=AND(NOT(DictValue(None,DictValue('overlap_type','total'))),
                                                                NOT(DictValue(None,AND(DictValue('overlap_type','partial'),
                                                                                       DictValue('overlap_fraction',Min(settings['antismash_minimal_overlap'])))))))
    
def create_length_criteria(settings):
    return {'core_genes': Length(Min(settings['gene_cluster_min_length']), req_type=list)}
    
def create_COG_criteria(settings):
    return {'COG_avg': Max(settings['gene_cluster_max_cog'])}
    
def generate_setting_to_filter(setting_keyword, setting_value):
    parts = setting_keyword.split('_')
    if parts[-1] == 'all':
        genetype = parts[-2]
        flanking = True
    else:
        genetype = parts[-1]
        flanking = False
    keyword = 'genes_%s' %genetype
    if flanking:
        keyword += '_all'
        
    filterdict = {keyword: Length(Min(setting_value), req_type=list)}
    return filterdict

def create_filter_criteria(settings):
    filter_dict_all = {}
    if settings['filter_antismash']:
        filter_dict_all.update(create_antismash_criteria(settings))
    keywords = ['gene_cluster_domains_min_biosyn',
                'gene_cluster_domains_min_biosyn_all',
                'gene_cluster_domains_min_transporter',
                'gene_cluster_domains_min_transporter_all',
                'gene_cluster_domains_min_peptidase',
                'gene_cluster_domains_min_peptidase_all',
                'gene_cluster_domains_min_regulator',
                'gene_cluster_domains_min_regulator_all',
                'gene_cluster_domains_min_kripp',
                'gene_cluster_domains_min_kripp_all',
                ]
                
    for keyword in keywords:
        setting_value = settings[keyword]
        filter_dict_keyword = generate_setting_to_filter(keyword, setting_value)
        filter_dict_all.update(filter_dict_keyword)
        
    filter_dict_all.update(create_length_criteria(settings))
    filter_dict_all.update(create_COG_criteria(settings))
    
    return filter_dict_all

if __name__ == '__main__':
    args,settings = parse_arguments()
    logger = return_logger('gene_cluster_formation', False)
    setup_loggers(settings)
    
    # Temporary to disable RREFinder until it is reimplemented
    settings['run_rre'] = False
    
    name = settings['name']
    path = os.path.abspath(os.path.join(settings['outputfolder'],name)) + os.sep 
    operon_namebase = 'RiPP_gene_cluster'
    # Step 0: Preparation
    
    genome_path = path + 'Genomes/'
    path_out = path + 'Output/'
    data_path = path + 'data/'
    pickle_path = data_path + 'Pickles/'
    COG_path = path + 'COG/'
    html_path = data_path + 'html_files/'
    domain_path = data_path +'domains/'
    blast_muscle_path = path + 'precursor_comparison/'
    hmm_path = path + 'hmm/'
    jaccard_path = path + 'domain_comparison/'
    RRE_path = path + 'RRE/'
    jaccard_prec_path = path + 'domains_precursor_overlap/'
    
    operon_path = path_out + 'Gene_clusters/'
    entry_path = path_out + 'Entries/'
    cytoscape_file = path_out + 'all_gene_clusters.ssn'
    hmm_path_pfam = hmm_path + 'pfam/'
    hmm_path_tigrfam = hmm_path + 'tigrfam/'
    operons_js_file = os.path.join(html_path,'Operons_dict.js')
    operon_group_js_file = os.path.join(html_path,'Operons_groups.js')
    filter_scripts = os.path.join(html_path,'Filter_operons.js')
    index_file_all = os.path.join(path_out,'Index.html') 
    summary_file = os.path.join(path_out,'%s_overview.tsv' %(name))
    
    jquery_path = os.path.join(data_path,'jquery','jquery-3.4.1.min.js')
    tablesorter_path=os.path.join(data_path,'jquery','tablesorter-master','js','jquery.tablesorter.js')
    
    paths = dict(out=path_out,pickles=pickle_path,html=html_path,domain=domain_path,genome=genome_path,COG=COG_path,\
            precursor=blast_muscle_path,hmm_basic=hmm_path,jaccard=jaccard_path,operons=operon_path,entries=entry_path,\
            RRE=RRE_path,base=path,operons_js_file=operons_js_file,operon_group_js_file=operon_group_js_file,filter_scripts=filter_scripts,\
            jquery_path=jquery_path,tablesorter_path=tablesorter_path,jaccard_prec_path=jaccard_prec_path,index_file=index_file_all)
    settings['paths'] = paths
    settings['hmm_paths'] = ((hmm_path_pfam, settings['pfam_db_path']),(hmm_path_tigrfam, settings['tigrfam_db_path']))
    
    check_paths(settings)
    
    for p in blast_muscle_path,jaccard_path,jaccard_prec_path,hmm_path,hmm_path_pfam,hmm_path_tigrfam,path_out,operon_path,entry_path:
        if not os.path.isdir(p):
            os.mkdir(p)
            
    EC_classes,all_domains,domain_descr,domain2ec,mibig_domaindict,genome_dict = initialize(settings,name)
    
    # Temp fix for backwards compatibility
    for gene in genome_dict.itergenes():
        gene.operon = None
    
    genetypes = {'precursor':'Putative Precursors','kripp':'Known RiPP domains',
                 'peptidase':'Peptidase domains','transporter':'Transporter domains',
                 'regulator':'Regulatory domains','biosyn':'Biosynthetic domains',
                 'other':'Other domains'}
                 
    attributes = ['all_domains','length','length_core','COG_avg','COG','score','mibig','nr_antismash','filtered']
    for category in ['biosyn','regulator','kripp','peptidase','transporter']:
        attributes.append(category)
        attributes.append('%s_all' %category)
        attributes.append('genes_%s' %category)
        attributes.append('genes_%s_all' %category)
                 
    color_dict = OrderedDict([('precursor', '#ff0000'),
                              ('kripp', '#66ff33'), 
                              ('peptidase', '#ff8c00'), 
                              ('transporter', '#7d9ec0'),
                              ('regulator', '#cc00ff'), 
                              ('biosyn', '#008b00'), 
                              ('other', '#cccccc')])
                              
    score_order = ['kripp','transporter','regulator','peptidase','biosyn']
    
    # Step 1: Search or load operons
    operons = build_operons(genome_dict, settings)
    
    # Step 2: extend found operons
    logger.info('Extending gene clusters')
    extend_operons(operons, settings['genecluster_extension'])
    logger.info('Setting quick_dict')
    operons.set_quick()
    logger.info('Setting precursors')
    operons.set_operon_precursors()
    
    # Step 3: Find all pfam/tigrfam domains
    logger.info('Detecting protein from Pfam/TIGRFAM domains')
    _ = hmm.main(operons,settings)

    # Step 3.5 (Optional)
    # Check overlap with other tools
    # 1) antiSMASH; 2) RRE-finder; to be implemented
    
    if settings['antismash_parse']:
        logger.info('Parsing antiSMASH data')
        all_clusters = antismash.parse_main(operons,genome_dict,settings)

    # Step 4: Filter operons by
        # 4a) length
        # 4b) COG
        # 4c) domains
    logger.info('Filtering gene clusters')
    assign_gene_types(operons, mibig_domaindict, all_domains)
    set_operon_domains(operons,mibig_domaindict, all_domains, genetypes)
    set_operon_unique_domains(operons, mibig_domaindict)
    
    operons_nr_old = len(operons)
    operons = filter_operons(operons, settings)
    logger.info('Filtered %i gene clusters; %i remaining' %(operons_nr_old - len(operons), len(operons)))

    logger.info('Preparing gene clusters for next analyses')
    logger.debug('Scoring gene clusters')
    score_operons(operons,score_order, settings)
    logger.debug('Setting gene cluster paths')
    set_operon_paths(operons, operon_path)
    logger.debug('Setting gene cluster paths')
    operons.set_gene_to_operon()
    

#    if settings['run_rre']:
#        RRE_wrapper(path,RRE_path,genome_dict,operons,settings)
#    else:
#        for operon in operons:
#            operon.nr_RRE = 'N\\A'

    # Step 5: Organize operons
        # 5a) By BLAST similarity of precursors
        # 5b) By Jaccard index of core operon
        # 5c) By E.C. families
        # 5d) By genomes
        # 5e) Both precursor and domain similarity
        
    # Operoncollection 
    group_collections = {}
    
    if settings['calc_networks']:
        logger.info('Forming gene cluster families')
        # Set initial pairs info
        keywords_pairs = ['precursor','jaccard','precursor_jaccard']
        logger.debug('Setting initial pairs')
        operons.set_pairs(keywords_pairs)
        set_pairs(mibig_domaindict,keywords_pairs)
        # 5a)
        logger.debug('Forming families based on precursors')
        precursor_group_names,precursor_operon_collections,ops,original_pairs_prec,sp,sg = group_precursors.main(operons,mibig_domaindict,settings)
        group_collections['precursor'] = precursor_operon_collections
        
        # 5b)
        logger.debug('Forming families based on protein domains')
        jaccard_group_names,jaccard_group_collections,jd,fp = group_jaccard.main(operons,jaccard_path,mibig_domaindict,settings)
        group_collections['jaccard'] = jaccard_group_collections
        
        # 5c)
        logger.debug('Grouping gene clusters based on E.C. domains')
        EC_group_collections = group_EC.main(operons,domain2ec)
        group_collections['EC'] = EC_group_collections
        
        # 5d)
        logger.debug('Grouping gene clusters based on E.C. domains')
        genome_collections = make_collections_per_genomes(operons,genome_dict)
        group_collections['genome'] = genome_collections
        
        # 5e)
        logger.debug('Forming families based on precursors and protein domains')
        grouped_both = group_pj.main(operons,mibig_domaindict,original_pairs_prec,jd,settings)
        group_collections['precursor_jaccard'] = grouped_both
        set_collection_paths(group_collections,settings)
        
    prep_collections(group_collections)
    
    # Step 6: Write SVG files
        # One html page per operon
    logger.info('Writing output')
    logger.debug('Writing out gene cluster SVGs')
    arrow.main(operons,all_domains,domain_descr,color_dict,genetypes,html_path,group_collections,settings)
    
    # Step 7: Writing output files    
    
    # Write genbank files
    logger.debug('Writing operon genbank files')
    wgbk.write_genbank(operons)
    logger.debug('Writing operon fasta files')
    wgbk.write_fasta_all(operons)
    
    # Write ssn file for cytoscape
    logger.debug('Writing cytoscape file')
    write_cytoscape(operons, mibig_domaindict, cytoscape_file, settings)

    # Write the operon js files to go with the html pages
    logger.debug('Writing javascript files for html pages')
    woj.write_operon_dict(operons, mibig_domaindict, attributes, operons_js_file)
    woj.write_group_dict(group_collections,operon_group_js_file)    
    
    # Write the html pages for the index and groups
    logger.debug('Writing index html pages')
    headers = {'precursor':'Precursor-based clustering',
               'jaccard':'Domain-based clustering',
               'genome':'Operons by genome',
               'EC':'Operons by enzyme class',
               'precursor_jaccard':'Clustered by precursors and domains'}
    t = html.main(settings, operons, group_collections, headers, domain_descr, genome_dict)

    # Write a summary
    logger.debug('Writing summary')
    write_operons(operons, settings, summary_file)
    
    logger.info('Finished')
     
    exit()
    
    



                            
                        

        
                        
