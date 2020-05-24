# License: GNU Affero General Public License v3 or later

#COG script
#Original algorithm developed by Michaelis Hadjithomas
#Script heavily modified and adapted by Alexander Kloosterman


import os
import matplotlib.pyplot as plt
import subprocess
import random
import time
import numpy as np
import pylab as pyplt

#from scipy.stats import skewnorm

#from lmfit.models import SkewedGaussianModel
#from lmfit_wrapper import fit_truecogs
from lib import autovivify, sortdictkeysbylenvalues, fuse_dict_list, run_mcl
from log import return_logger

#from lmfit.models import SkewedGaussianModel
#from lib.lmfit_wrapper import fit_truecogs
#from lib.lib import autovivify, sortdictkeysbylenvalues,fuse_dict_list

logger = return_logger(__name__, False)

def generate_scaffold_lists(genome_dict):
    gene2scf = {}
    scf_genelist = {}
    for genome in genome_dict:
        for scf in genome:
            scf_genelist[scf.name] = []
            for gene in scf.itergenes(intergenic=False):
                gene2scf[gene.name] = scf.name
                scf_genelist[scf.name].append(gene.name)
    return gene2scf,scf_genelist
    

def read_allvall(allVallfile,genome_dict,genomes=False):
    '''Parses the allvall blast results from a single file'''
    #Returns a gen2tax dict and an expanded score dictionary for the bbh analysis    
    score_dict_exp = autovivify(levels=4,final=list)
    with open(allVallfile, 'r') as f:
        for line in f:
            if not '#' in line:
                line=line.rstrip()
                # in the BLAST, outfmt 6 is used with the following keywords
                # qseqid sseqid pident length mismatch gapopen qlen slen evalue bitscore
                # qseqid and sseqid are already split up in gene and taxon
                (a,b,pid,length,mismatch,gapopen,ql,sl,evalue,bits)=line.split('\t')
                t1 = genome_dict[a].genome.name
                t2 = genome_dict[b].genome.name
                pid = float(pid)
                ql = float(ql)
                length = int(length)
                # Hits are only added if the alignment length is at least 75% as long as the query
                if float(length) < 0.75*ql:
                    continue
                # If an iterable of genomes is added, these are considered the only genomes allowed
                # This is mostly done to allow analysis of smaller subgroups of genomes without having to rerun the BLAST
                if genomes and (t1 not in genomes or t2 not in genomes):
                    continue
                bits = float(bits)
                evalue = float(evalue)
                score_dict_exp[t1][t2][a][b] = [ql,sl,pid,bits,evalue]
    return(score_dict_exp)

def get_bbh(d_allvall,min_simil = 0.0):
    # Finds all bidirectional best hits from the previously parsed dictionary
    # Takes into account cases where multiple proteins are the 'best hits'; adds them all (if they're bbhs)
    pairs_per_genome = autovivify(levels=1,final=list)
    all_pairs = set()
    pair_dict = autovivify(levels=1,final=set)
    for genome_q in d_allvall:
        genomes_s = d_allvall[genome_q]
        for genome_s in genomes_s:
            if genome_s <= genome_q:
                continue
                # Only check once for each genome pair
            genes_q = genomes_s[genome_s]
            for gene_q in genes_q:
                best_score,best_ev,best_genes = get_best_scoring_genes(d_allvall,genome_q,genome_s,gene_q)
                for best_gene in best_genes:
                    # Get all the best hits
                    best_score_rev,best_ev,best_genes_rev = get_best_scoring_genes(d_allvall,genome_s,genome_q,best_gene)
                    if gene_q in best_genes_rev:
                        # This best hit is a bbh
                        sim_ab = d_allvall[genome_q][genome_s][gene_q][best_gene][2]
                        sim_ba = d_allvall[genome_s][genome_q][best_gene][gene_q][2]
                        if sim_ab > min_simil and sim_ba > min_simil:
                            # Add both ways here
                            all_pairs.add((best_gene,gene_q))
                            all_pairs.add((gene_q,best_gene))
                            pairs_per_genome[(genome_q,genome_s)].append((gene_q,best_gene))
                            pairs_per_genome[(genome_s,genome_q)].append((best_gene,gene_q))
                            pair_dict[gene_q].add(best_gene)
                            pair_dict[best_gene].add(gene_q)
    return(all_pairs,pairs_per_genome,pair_dict)

def get_best_scoring_genes(d_allvall,genome_a,genome_b,gene_a):
    # Finds the gene(s) with the highest score
    if genome_a in d_allvall and genome_b in d_allvall[genome_a] and gene_a in d_allvall[genome_a][genome_b]:
        genes_b = d_allvall[genome_a][genome_b][gene_a]
    else:
        genes_b = []
    best_score = 0
    best_gene = ""
    best_ev = 10
    best_genes = []
    for gene_b in genes_b:
        (sl,ql,score,bits,ev) = d_allvall[genome_a][genome_b][gene_a][gene_b]
        if bits > best_score or (bits == best_score and ev < best_ev):
            best_score = bits
            best_gene = gene_b
            best_ev = ev
    # Find all the genes with the same score:
    best_genes.append(best_gene)
    for gene_b in genes_b:
        (sl,ql,score,bits,ev) = d_allvall[genome_a][genome_b][gene_a][gene_b]
        if bits == best_score and ev == best_ev and gene_b != best_gene:
            # Multiple genes can give the same BLAST result. These are all considered as potential truecogs further on.
            best_genes.append(gene_b)
    return(best_score,best_ev,best_genes)


def get_truecog_pairs(pairs,all_pairs,genome_dict):
    # Figure out which of the pairs are truecog-pairs (flanked by other pairs)
    truecog_pairs_per_genome = {}
    all_truecog_pairs = set()
    true_pair_dict = autovivify(levels=1,final=set)
    for genome_a,genome_b in pairs:
        pairs_g = pairs[genome_a,genome_b]
        truecog_pairs_per_genome[genome_a,genome_b] = []
        for genename_a,genename_b in pairs_g:
            # Get the objects from the names
            gene_a = genome_dict[genename_a]
            gene_b = genome_dict[genename_b]
            if (gene_a.left_flank and gene_b.left_flank and gene_a.right_flank and gene_b.right_flank):
                if ((gene_a.left_flank.name,gene_b.left_flank.name) in all_pairs and (gene_a.right_flank.name,gene_b.right_flank.name) in all_pairs) or\
                   ((gene_a.left_flank.name,gene_b.right_flank.name) in all_pairs and (gene_a.right_flank.name,gene_b.left_flank.name) in all_pairs):
            # If both genes have flanking genes, and the flanking genes are bbhs between the genomes...
                    # ... add these genes as a truecog-pair
                    truecog_pairs_per_genome[genome_a,genome_b].append((genename_a,genename_b))
                    all_truecog_pairs.add((genename_a,genename_b))
                    true_pair_dict[genename_a].add(genename_b)
                    true_pair_dict[genename_b].add(genename_a)
    return(truecog_pairs_per_genome,all_truecog_pairs,true_pair_dict)
    
def add_genome_cog(base_group,genome,true_pair_dict,base_truecogs,genome_dict):
    # Adds a genome to a group of genomes, and returns which truecogs remain
    truecogs_missing = 0
    new_truecogs = []
    for truecogs in base_truecogs:
        # truecogs is a group of truecogs, at least one per base genome
        candidates = autovivify(levels=1,final=int)
        # candidates keeps track of candidate genes that might be a truecog in the new group
        # The number indicates the amount of genomes of the base_genomes that this gene is a truecog_pair with
        for base_genome in base_group:
            # Analyze all the genes of this truecog group belonging to this genome 
            genome_truecogs = [tc for tc in truecogs if genome_dict[tc].genome.name == base_genome] 
            basegenome_candidates = set()
            added = False
            for truecog in genome_truecogs:
                gene_true_pairs = true_pair_dict[truecog]
                for gene in gene_true_pairs:
                    if genome_dict[gene].genome.name == genome:
                        basegenome_candidates.add(gene)
                        added = True
                        # One of the genes of the new genome is a truecog_pair with one of the truecogs of thise base genome
                        # If this is true for all genomes, the truecog is conserved
            for cand in basegenome_candidates:
                candidates[cand] += 1
            if not added:
                # None of the genes belonging to this trueCOG and one of the base genomes are a true_pair with any gene of the genome to be added
                # The truecog is not conserved when this genome is added
                truecogs_missing += 1
                break
                
        else:    
            # Keep only candidates that are true_pairs with all other genomes 
            # This step prevents the case where multiple candidate truecogs exist, but none are conserved true_pairs between all genomes in the basegroup        
            to_add = []
            for c in candidates:
                if candidates[c] < len(base_group):
                    continue
                else:
                    to_add.append(c)
            if len(to_add) == 0:
                truecogs_missing += 1
            else:
                new_truecogs.append(truecogs + tuple(to_add))
    return truecogs_missing,new_truecogs
    
def add_best_genome(base_group,truecog_pairs,true_pairs_gene,genomes,criteria,genome_dict,base_truecogs = [],min_truecogs=10):
    if base_truecogs == []:
        if len(base_group) != 2:
            raise StandardError('Too many genomes without a specified truecog group')
        base_truecogs = truecog_pairs[base_group[0],base_group[1]]
    new_truecogs = {}
    best_all_cutoffs = []
    best_average_cutoff = 0
    best_min_cutoff = 0
    best_genome = None
    best_truecogs = []   
    for genome in genomes:
#        print(genome)
        if genome in base_group:
#            print('Genome in base group! Continuing...')
            continue
        truecogs_missing,new_truecogs_genome = add_genome_cog(base_group,genome,true_pairs_gene,base_truecogs,genome_dict)
        if criteria == 'nr_truecogs':
            if len(new_truecogs_genome) > len(best_truecogs):
                best_truecogs = new_truecogs_genome
                best_genome = genome
                best_average_cutoff,best_all_cutoffs = None,None
        
        elif criteria == 'mc':
            new_truecogs[genome] = new_truecogs_genome
            if len(new_truecogs_genome) >= min_truecogs:
                cutoff_dict,_,all_cutoffs,_ = calculate_taxon_pairwise(new_truecogs_genome,genome_dict,score_dict_exp,cutoff_mode,cutoff_factor)
                min_cutoff = min(all_cutoffs)
                average_cutoffs_per_genome = get_average_cutoffs(cutoff_dict)
                average_cutoff = average_cutoffs_per_genome[genome]
                # Get largest distance (lowest cutoff) between all genome pairs
                # If the distance here is lower than previous, use this genome as the best candidate
                if min_cutoff > best_min_cutoff:
                    best_genome = genome
                    best_average_cutoff = average_cutoff  
                    best_min_cutoff = min_cutoff                                       
                    best_all_cutoffs = all_cutoffs
                    best_truecogs = new_truecogs_genome
#                    print('Best genome replaced with %s! New best_cutoff: %s, new nr_truecogs: %i, min_cutoff: %.2f' \
#                        %(best_genome,best_average_cutoff,len(new_truecogs_genome),min_cutoff))
                else:
                    pass
#                    print('Best genome not replaced. Previous: %s, tested: %s, prev_cutoff: %s, tested_cutoff: %s, prev_tc: %i, tested_tc: %i, prev_mc: %.2f, tested_mc: %.2f'\
#                             %(best_genome,genome,best_average_cutoff,average_cutoff,len(best_truecogs),len(new_truecogs_genome),best_min_cutoff,min(all_cutoffs)))
    return(best_genome,best_average_cutoff,best_all_cutoffs,best_truecogs,new_truecogs)
    
def add_random_genome(base_group,truecog_pairs,true_pairs_gene,genomes,criteria,genome_dict,base_truecogs = [],min_truecogs=10,min_cutoff_req=50.0):
    # Randomly try adding genomes until one is found that fulfills the criteria
    genomes_left = [i for i in genomes if i not in base_group]
    while genomes_left != []:
        genome = random.choice(genomes_left)
        truecogs_missing,new_truecogs_genome = add_genome_cog(base_group,genome,true_pairs_gene,base_truecogs,genome_dict)
        
        if criteria == 'nr_truecogs':
            if len(new_truecogs_genome) >= min_truecogs:
                return(genome,None,None,new_truecogs_genome,{})
            else:
                genomes_left.remove(genome)
                
        elif criteria == 'mc':
            new_truecogs[genome] = new_truecogs_genome
            if len(new_truecogs_genome) >= min_truecogs:
                cutoff_dict,_,all_cutoffs,_ = calculate_taxon_pairwise(new_truecogs_genome,genome_dict,score_dict_exp,cutoff_mode,cutoff_factor)
                min_cutoff = min(all_cutoffs)
                average_cutoffs_per_genome = get_average_cutoffs(cutoff_dict)
                average_cutoff = average_cutoffs_per_genome[genome]
                # Get largest distance (lowest cutoff) between all genome pairs
                # If the distance here is lower than previous, use this genome as the best candidate
                if min_cutoff > min_cutoff_req:
                    return(genome,average_cutoff,all_cutoffs,new_truecogs_genome,{})
                else:
                    genomes_left.genome(genome)
    else:
        logger.debug('No more genomes left to sample - breaking function')
        return None,None,None,None,None
    
def get_average_cutoffs(cutoff_dict,verbose=False):
    all_cutoffs = {}
    for a in cutoff_dict:
        for b in cutoff_dict[a]:
            cutoff = cutoff_dict[a][b]
            if a not in all_cutoffs:
                all_cutoffs[a] = []
            if b not in all_cutoffs:
                all_cutoffs[b] = []
            all_cutoffs[a].append(cutoff)
            all_cutoffs[b].append(cutoff)
    average_cutoffs = {}
    for g in all_cutoffs:
        if verbose:
            logger.debug(g,np.average(all_cutoffs[g]))
        average_cutoffs[g] = np.average(all_cutoffs[g])
    return average_cutoffs
    
def make_cog_group(group,truecog_pairs,true_pair_dict,genome_dict):
    if len(group) < 2:
        return None
    base_group = [group[0],group[1]]
    base_truecogs = truecog_pairs[group[0],group[1]]
    for genome in group[2:]:
        truecogs_missing,new_truecogs = add_genome_cog(base_group,genome,true_pair_dict,base_truecogs,genome_dict)
        base_group += [genome]
        logger.debug('group size: %s\tnr of truecogs: %s\t' %(len(base_group),len(new_truecogs)))
        base_truecogs = new_truecogs
    return(base_truecogs)

def get_unique_truecogs(truecogs):
    seen = set()
    unique_truecogs = []
    for truecog_group in truecogs:
        for gene in truecog_group:
            if gene in seen:
                break
            else:
                seen.add(gene)
        else:
            unique_truecogs.append(truecog_group)
    return unique_truecogs
    
def calculate_taxon_pairwise(true_list,genome_dict,score_dict_exp,cutoff_mode,cutoff_factor,\
                                genomes_allowed=False,cutoff_sgm = 0.05,outf=False,bins='auto'):
    # Calculate the cutoff for each pair of genomes based on the truecogs given
    # First cleanup: get all unique truecogs (remove extra truecogs due to some truecogs pairing up with multiple genes)
    unique_truecogs = get_unique_truecogs(true_list)
    
    # Now gather all the similarity scores for each pair of genomes
    true_scores = autovivify(2, list)
    for cog_group in unique_truecogs:        
        for a in cog_group:
            for b in cog_group:
                t_a = genome_dict[a].genome.name
                t_b = genome_dict[b].genome.name
                if t_a >= t_b or (genomes_allowed and (t_a not in genomes_allowed and t_b not in genomes_allowed)):
                    continue                
                
                score_ab = score_dict_exp[t_a][t_b][a][b][2]
                score_ba = score_dict_exp[t_b][t_a][b][a][2]
                
                score = (score_ab+score_ba)/2.0 # Just in case the score from prot a vs b and b vs a is slightly different, take the average
                
                if t_a < t_b:
                    true_scores[t_a][t_b].append(score)
                else:
                    true_scores[t_b][t_a].append(score)
               
    if outf:  
        cutoff_table=open(outf,'wt')
    cutoff_dict = autovivify(2, float)
    all_cutoffs = []
    
    histogram_data = autovivify(2, list)
    
    for t1 in true_scores:
        for t2 in true_scores[t1]:
            scores = true_scores[t1][t2]
            scores_sorted = tuple(sorted(scores))
            if scores_sorted in scores_tested:
                # To save time when doing multiple group formations in combination with sgm
                pars,cutoff = scores_tested[scores_sorted]
            else:
                y,bounds = np.histogram(scores,bins=bins)
#                highest_point = max(y)
                if cutoff_mode == 'sgm':
                    result = fit_truecogs(scores,bins='auto',plot=False,model=SkewedGaussianModel)
                    if result:
                        pars = result.params
                        skew = skewnorm(pars['gamma'],pars['center'],pars['sigma'])
                        cutoff = skew.ppf(cutoff_sgm)
                # Sometimes the data is fit to just the tail of the sgm, and the cutoffs are outside the 0-100 range
                if cutoff_mode == 'gm' or not result or (result and (cutoff > 100.0 or cutoff < 0.0)):
                    mean = np.mean(scores)
                    std = np.std(scores)
                    cutoff = round(mean-cutoff_factor*(std),2)
                    pars = {'center':mean,'sigma':std,'amplitude':1,'gamma':0}
                
                scores_tested[scores_sorted] = pars,cutoff
            cutoff_dict[t1][t2]=cutoff
            all_cutoffs.append(cutoff)
            histogram_data[t1][t2] = [pars,cutoff]
            
            if outf:
                cutoff_table.write( "%s\t%s\tmean=\t%s\tcutoff=\t%s\n" % (t1,t2,np.mean(true_scores[t1][t2]),cutoff_dict[t1][t2]) )
    if outf:
        cutoff_table.close()
        
    return cutoff_dict,true_scores,all_cutoffs,histogram_data


def draw_histograms_pairwise(true_scores,histogram_data,f):
    logger.info('Draw histograms function')
    font = 10
    distinct_taxa = len(true_scores) + 1
    figsize = (2*distinct_taxa,2*distinct_taxa-1)
    sp1 = 0
    minX=100
    for t1 in true_scores:
        for t2 in true_scores[t1]:
            pars,cutoff = histogram_data[t1][t2]
            if sorted(true_scores[t1][t2])[0] < minX:
                minX = int(sorted(true_scores[t1][t2])[0])
            if cutoff < minX:
                minX = cutoff
    logger.debug('Minx: %s' %minX)
    lower_bound = int(minX/10)*10
    x_incr = (100-minX) / 10.0
    plt.rcParams["figure.figsize"] = figsize
    for t1 in sorted(true_scores.keys()):
        logger.debug(t1)
        sp2 = 0
        for t2 in sorted(true_scores[t1].keys()):
            pars,cutoff = histogram_data[t1][t2]
            logger.debug("%s\t%s\t%s\t%f\t%f" %(len(true_scores[t1][t2]), t1,t2,np.mean(true_scores[t1][t2]), np.std(true_scores[t1][t2])))
            # the histogram of the data
            pyplt.subplot2grid((distinct_taxa,distinct_taxa-1), (sp1,sp1+sp2))
            y,x,_ = plt.hist(true_scores[t1][t2], bins='auto', normed=1, facecolor='r', alpha=0.75)
#                print(y,x)
            if sp2 ==0:
                plt.ylabel(t1, fontsize = font)
            if sp1 == 0:
                plt.title(t2, fontsize = font)
            x_plot = np.linspace(lower_bound,100,1000)
            center = pars['center']
            sigma = pars['sigma']
            gamma = pars['gamma']
            ampl = pars['amplitude']
            
            skew = skewnorm(gamma,center,sigma)
            y_plot = ampl* skew.pdf(x_plot)
            plt.plot(x_plot,y_plot,color='b')
            
            maxY = 1.1 * max(y.max(),max(y_plot))
            
            plt.text(lower_bound+x_incr*0.5, maxY*0.9, r'mu=%.2f' %center, fontsize=font - 2)
            plt.text(lower_bound+x_incr*0.5, maxY*0.8, r'sigma=%.2f' %sigma, fontsize=font - 2 )
            plt.text(lower_bound+x_incr*0.5, maxY*0.7, r'gamma=%.2f' %gamma , fontsize = font - 2)
            plt.annotate('Cutoff\n%s' % round(cutoff,2),xy=(max(cutoff-2,minX+5),maxY*0.1),
                         xytext=(max(cutoff-7,minX+2),maxY*0.3),
                         fontsize = font - 2,
                         arrowprops=dict(facecolor='darkgreen', width = 2, headwidth=5, shrink=0.05, alpha=0.8),
                         )
            plt.axis([minX, 100, 0, maxY])
            plt.xticks(range(lower_bound,101,10), fontsize = font * 0.75)
#                print "%s %s %s" % (t1,t2,range(int(minX/10)*10,100,10))
            plt.tick_params(
                            axis='y',          # changes apply to the x-axis
                            which='both',      # both major and minor ticks are affected
                            left='off',      # ticks along the bottom edge are off
                            right='off',         # ticks along the top edge are off
                            labelleft='off') # labels along the bottom edge are off
            #plt.yticks(np.arange(0,maxY,maxY/3))
            plt.grid(False)
            sp2 += 1
#                if sp2+sp1 == 4:
#                    break
        sp1 += 1
#            if sp1 == 4:
#                break
    plt.savefig(f) 
    
def get_cutoffs_pairs(truecog_pairs,genomes_allowed):
    # Calculate for each genome pair the cutoff based on their shared truecogs
    cutoffs_pairs = {}
    for genome_a,genome_b in truecog_pairs:
        if genome_a not in genomes_allowed or genome_b not in genomes_allowed:
            continue
        base_truecogs = truecog_pairs[genome_a,genome_b]
        
        cutoff_dict,true_scores,all_cutoffs,_ = calculate_taxon_pairwise(base_truecogs,genome_dict,score_dict_exp,cutoff_mode,cutoff_factor)                                                        
        if len(all_cutoffs) != 1:
            logger.error('More than two genomes compared during pairwise comparison')
            raise ValueError('More than two genomes compared during pairwise comparison')
        cutoff = all_cutoffs[0]
        cutoffs_pairs[(genome_a,genome_b)] = cutoff
    return cutoffs_pairs

def form_basegroup(truecog_pairs,genomes_allowed=False,genome=False):
    # Find the two genomes that share the most truecog pairs
    best_pairs = 0
    best_genomes = tuple()
    for genome_a,genome_b in truecog_pairs:
        if genomes_allowed and (genome_a not in genomes_allowed or genome_b not in genomes_allowed):
            continue
        if genome and genome_a != genome and genome_b != genome:
            continue
        nr_pairs = len(truecog_pairs[genome_a,genome_b])
        if nr_pairs > best_pairs:
            best_pairs = nr_pairs
            best_genomes = (genome_a,genome_b)
    return(best_genomes,best_pairs)

def form_basegroup_mc(cutoffs_pairs=False,genomes=False,genome=False):
    # Out of a group of genomes, get the pair that has the highest average truecog similarity
    best_cutoff = 0
    best_genomes = tuple()
    for genome_a,genome_b in cutoffs_pairs:
        if genomes and (genome_a not in genomes or genome_b not in genomes):
            continue
        if genome and genome_a != genome and genome_b != genome:
            continue
        cutoff = cutoffs_pairs[(genome_a,genome_b)]
        if cutoff > best_cutoff:
            best_cutoff = cutoff
            best_genomes = (genome_a,genome_b)
    return(best_genomes,best_cutoff)
    

def groups_bottomup(truecog_pairs,true_pairs_gene,genomes,genome_dict,seedmethod='closest',seen=[],criteria='nr_truecogs',add_method="best",\
                        min_truecogs=10,min_cutoff_req=50.0,cutoffs_pairs=False):
    all_groups = []
    seeds = []
    singles = []
    break_after_first = False
    if criteria == 'mc' and cutoffs_pairs==False:
        cutoffs_pairs = get_cutoffs_pairs(truecog_pairs,genomes)
    while True:
        logger.debug('Starting group formation! Number of genomes left: %s' %len(genomes))
        if seedmethod == 'random':
            nr_pairs = 0
            while nr_pairs < min_truecogs:
                base_group = tuple(random.sample(genomes,2))
                base_pairs = truecog_pairs[base_group[0],base_group[1]]
                nr_pairs = len(base_pairs)
                seeds.append(base_pairs)
        else:
            if seedmethod != 'closest':
                seed_genome = seedmethod #Assumes the genome passed to this function
                break_after_first = True 
            else:
                seed_genome = False
            if criteria == 'nr_truecogs':  
                base_group,nr_pairs = form_basegroup(truecog_pairs,genomes,seed_genome)
            elif criteria == 'mc':
                base_group,best_cutoff = form_basegroup_mc(cutoffs_pairs,genomes,seed_genome)
                nr_pairs = len(truecog_pairs[base_group])
                current_mc = best_cutoff
            
        if len(base_group) < 2 or nr_pairs < min_truecogs:
            # No more suitable group can be formed
            logger.debug('No more suitable pairs can be formed with the remaining genomes (best group had %i genomes and %i pairs - breaking function' %(len(base_group),nr_pairs))
            singles.extend(genomes)
            break
        base_truecogs = truecog_pairs[base_group[0],base_group[1]]
        logger.debug('Base group: %s' %' '.join(base_group))
        logger.debug('Nr truecogs: %s' %(nr_pairs))
        if criteria=='mc':
            logger.debug('Cutoff: %s' %(current_mc))
        while True:
            if add_method == 'best':
                best_genome,best_average_cutoff,best_all_cutoffs,best_truecogs,new_truecogs = add_best_genome(base_group,truecog_pairs,true_pairs_gene,\
                                                                                    genomes,criteria,genome_dict,base_truecogs,min_truecogs)
            elif add_method == 'random':
                best_genome,best_average_cutoff,best_all_cutoffs,best_truecogs,new_truecogs = add_random_genome(base_group,truecog_pairs,true_pairs_gene,genomes,criteria,\
                                                                                genome_dict,base_truecogs,min_truecogs,min_cutoff_req=min_cutoff_req)
            if best_genome and len(best_truecogs) >= min_truecogs and\
               (criteria == 'nr_truecogs' or (criteria == 'mc' and (min(best_all_cutoffs) > min_cutoff_req or min(best_all_cutoffs) > current_mc))):
                base_group = base_group + (best_genome,)
                nr_old_truecogs = len(base_truecogs)
                base_truecogs = best_truecogs
                percent = 100 - float(len(base_truecogs))/nr_old_truecogs*100
                text = 'Added a new genome: %s! trueCOGs jump: %s to %s (lost %.2f percent); Group size: %i' %(best_genome,nr_old_truecogs,len(base_truecogs),percent,len(base_group))
                if criteria == 'mc':
                    old_cutoff = current_mc
                    current_mc = min(best_all_cutoffs)
                    text += '\nCutoff jump: %.2f to %.2f' %(old_cutoff,current_mc)
                logger.debug(text)
            else:
                logger.debug('Breaking after group formation - starting with the next!')
                text = 'Total group size: %s; Number of trueCOGs: %s' %(len(base_group),len(base_truecogs))
                if criteria == 'mc':
                    text += 'Lowest cutoff: %.2f' %current_mc
                    if current_mc < min_cutoff_req:
                        text += '\nGenomes still too far away. Adding to singles'
                        singles.extend(base_group)
                    else:
                        all_groups.append((base_group,base_truecogs))
                        genomes = [genome for genome in genomes if genome not in base_group and genome not in singles]
                elif criteria == 'nr_truecogs':
                    genomes = [genome for genome in genomes if genome not in base_group]                    
                    all_groups.append((base_group,base_truecogs))
                break
        if len(genomes) < 2 or break_after_first:
            logger.debug('Breaking function - %i genomes remaining' %(len(genomes)))
            if len(genomes) < 2:
                singles.extend(genomes)
            break
    return all_groups,singles

    
def bottomup_groups(genomes,truecog_pairs_per_genome,true_pair_dict,genome_dict,criterium,min_truecogs,min_group_size,settings):
    # Go through all the genomes, and use them as seeds
    # Make a lot of different groups with the same seed. 
    # Keep the ones that are large enough. 
    # Then determine for each genome the smallest amount 
    # of groups that cover it's connections to as many
    # genomes as possible.
    
    all_groups_found = []
    largest_group_per_genome = {}
    t0 = time.time()
    current_iter = 0
    add_method = settings['cog_addmethod']
    seedmethod = settings['cog_seedmethod']
    genomes_left_to_analyze = genomes[:]
        
    while genomes_left_to_analyze != []:
        current_iter += 1
        
        if seedmethod == 'total': # Try formation with all groups
            genome = genomes_left_to_analyze[0]
            genomes_left_to_analyze = genomes_left_to_analyze[1:]
        elif seedmethod == 'smart':
            # Prioritize formation with genomes for which no group was found yet
            if genomes_left_to_analyze == []:
                # Get the genome which has the smallest group
                logger.debug('No more genomes without a group found but still not passed minimum iterations')
                logger.debug('Selecting the genome with the smallest group')
                smallest = len(genomes)
                for genome in largest_group_per_genome:
                    if largest_group_per_genome[genome] < smallest:
                        smallest_group_genome = genome
                        smallest = largest_group_per_genome[genome]
                genome = smallest_group_genome
            else:
                genome = genomes_left_to_analyze[0]
        
        
        
        genomes_remain = genomes[:]
        logger.debug('Seeding with genome %s' %genome)
        logger.debug('At iteration %i' %current_iter)
        while True:
        
            all_groups,singles = groups_bottomup(truecog_pairs_per_genome,true_pair_dict,genomes_remain,genome_dict,\
                                                    seedmethod=genome,criteria=criterium,min_truecogs=min_truecogs,add_method=add_method)
            if len(all_groups) == 0:
                break
            group = all_groups[0][0] # First 0-index: first group; Second 0-index: genomes (Second 1-index: base_truecogs)
            truecogs = all_groups[0][1]
            logger.debug('Genome: %s, Group found: %i' %(genome,len(group)))
            genomes_remain = [i for i in genomes_remain if i not in group]
            genomes_remain.append(genome)
            if len(group) >= min_group_size and (sorted(group),truecogs) not in all_groups_found:
                all_groups_found.append((sorted(group),truecogs))
            if len(genomes_remain) < min_group_size:
                break
        
        # Get the largest group per genome
        improved = False
        for group,truecogs in all_groups_found:
            for genome_2 in group:
                if genome_2 not in largest_group_per_genome or len(group) > largest_group_per_genome[genome_2]:
                    largest_group_per_genome[genome_2] = len(group)
                    improved = True
        if not improved:
            t1 = time.time()
            if (settings['cog_stop_split_iterations'] and current_iter >= settings['cog_stop_split_iterations']) or\
               (settings['cog_stop_split_time'] and (t1-t0)/60.0 > settings['cog_stop_split_time']):
                # No more group finding
                logger.debug('Breaking group finding after no improvement for %i iterations and %.2f minutes' %(current_iter,(t1-t0)/60.0))
                break
        else:
            logger.debug('Improvement found!')
            t0 = time.time()
            current_iter = 0
        genomes_left_to_analyze = [i for i in genomes if i not in largest_group_per_genome]
    
    groups_dict = {}
    truecog_dict = {}
    counter = 1
    for group,truecogs in all_groups_found:
        groups_dict[counter] = group
        truecog_dict[counter] = truecogs
        counter += 1
    return groups_dict,truecog_dict

def get_groups_per_genome(all_groups,genomes):
    groups_per_genome = {}
    for genome in genomes:
        groups_per_genome[genome] = []
        for counter in all_groups:
            group = all_groups[counter]
            if genome in group:
                groups_per_genome[genome].append(counter)
    return groups_per_genome

def required_groups_per_genome(groups_per_genome,all_groups,genomes_required):
    req_groups_per_genome = {}
    genomes_covered_per_genome = {}
    all_req = []
    for genome in groups_per_genome:
        groups_genome = groups_per_genome[genome]
        # Each time, add the group that contributes the most genomes
        # In case of ties, go for the bigger group
        groups = dict([(i,all_groups[i]) for i in groups_genome])
        genomes_covered = []
        keys_req = []
        c = 0
        while sorted(genomes_covered) != sorted(genomes_required):
            c += 1
            best_coverage = []
            best_key = ''
            for key in groups:
                group = groups[key]
                genomes_new = [i for i in group if i not in genomes_covered]
                if len(genomes_new) > 0 and (len(genomes_new) > len(best_coverage) or\
                   (len(genomes_new) == len(best_coverage) and len(group) > len(groups[best_key]))):
                    best_key = key
                    best_coverage = genomes_new
            if best_key == '':
                break
            genomes_covered.extend(best_coverage)
            keys_req.append(best_key)
        req_groups_per_genome[genome] = keys_req
        all_req.extend([i for i in keys_req if i not in all_req])
        genomes_covered_per_genome[genome] = genomes_covered
    return req_groups_per_genome,all_req,genomes_covered_per_genome
    
def convert_truecogs_to_pairs(truecogs):
    pairs = set()
    for truecog_group in truecogs:
        for gene in truecog_group:
            for gene2 in truecog_group:
                pairs.add((gene,gene2))
    return pairs


def filter_cogs(cutoff_dict,score_dict_exp,base_truecogs,genomes,propagate_truecogs,truecog_pairs):
    passed=autovivify(levels=2,final=float)
    for t1 in score_dict_exp:
        t2s = score_dict_exp[t1]
        for t2 in t2s:
            if t2 == t1:
                continue
            if t1 not in genomes or t2 not in genomes:
                continue
            g1s = t2s[t2]
            for g1 in g1s:
                g2s = g1s[g1]
                for g2 in g2s:
                    t0 = time.time()
                    [ql,sl,score,bits,ev] = g2s[g2]
                    if t1 < t2:
                        cutoff = cutoff_dict[t1][t2]
                    else:
                        cutoff = cutoff_dict[t2][t1]
                    if score >= cutoff:
                        if g1 < g2:
                            passed[g1][g2]=score
                        else:
                            passed[g2][g1]=score 
                #Decide if true cogs should be forced as cogs even if below threshold
                    if propagate_truecogs and score < cutoff:
                        if (g1,g2) in truecog_pairs:
                            if g1 < g2:
                                passed[g1][g2] = 100
                            else:
                                passed[g2][g1] = 100
                                                       
    return(passed)


def scoreDB2file(scoreDB,filename):
    f = open(filename,'wt')    
    for a in scoreDB:
        for b in scoreDB[a]:
            if a > b:
                score = float(scoreDB[b][a])
            else:
                score = float(scoreDB[a][b])
            if score == 0.0:
                logger.debug("Something went wrong %s %s" %(a,b))
            f.write("%s\t%s\t%s\n" % (a,b,score))
    f.close()
    return(1)

def cog_synopsis(infile,genome_dict,group_statsf,gene_statsf,genomes_covered_per_genome):
    c = 0
    with open(infile) as f:
        group_stats = open(group_statsf,'w')
        gene_stats = open(gene_statsf,'w')
        count_cog = 0
        for line in f:
            count_cog += 1
            cog_group = line.rstrip().split('\t')
            count_tax = set()
            for gene in cog_group:
                count_tax.add(genome_dict[gene].genome.name)
            for gene in cog_group:
                gene_obj = genome_dict[gene]
                genomes_covered = genomes_covered_per_genome[gene_obj.genome.name]
                COG_total = len([i for i in count_tax if i in genomes_covered])
                COG_perc = COG_total / float(len(genomes_covered))
                               
                if not hasattr(gene_obj,'COG') or\
                   (hasattr(gene_obj,'COG') and len(gene_obj.COG_compared) < len(genomes_covered)):
                    gene_obj.update(COG=COG_perc,COG_total=COG_total,COG_compared=genomes_covered)
                    c += 1
                
                gene_stats.write( "%s\thomologs found in\t%s\tgenomes\n" % (gene,COG_total) )
            group_stats.write( "Group\t%s\thas representatives from\t%s\tgenomes\n" % (count_cog,len(count_tax)))
    logger.info('COG groups found: %i' %count_cog)
    logger.info('Assigned %s COG scores to genes' %c)
    return 1
    
def assign_COG_remaining_genes(genome_dict,genomes_allowed):
    # Give a COG score to all other genes
    c = 0
    for gene in genome_dict.itergenes(genomes=genomes_allowed,intergenic=False):
        if not hasattr(gene,'COG'):
            gene.update(COG=(1.0/len(genomes_allowed)),COG_total=1,COG_compared=genomes_allowed)
            c += 1
    logger.info('COG found %i unique genes' %c)
    return 1
    
def get_genomes_covered_per_genome(all_groups):
    genomes_covered_per_genome = {}
    for group,base_truecogs in all_groups:
        for genome in group:
            genomes_covered_per_genome[genome] = group
    return genomes_covered_per_genome
    
def clear_COG(genome_dict):
    # Removes COG scores (mostly for debugging purposes)
    for gene in genome_dict.itergenes():
        for attr in ['COG','COG_total','COG_compared']:
            if hasattr(gene,attr):
                delattr(gene,attr)
            
def main(path,allvallfile,name,settings,genome_dict,genomes_allowed=False):
    # Set a few filenames
    # All BBHs; input for mcl
    dokfile = path + 'bbh_pairs_%s.txt' %name
    # All groups of bbh as determined by mcl
    bbh_mcl = path + 'bbh_%s.mcl' %name
    
    min_simil_bbh = settings['min_simil_bbh']
    cutoff_mode = settings['cutoff_mode']
    distinct_taxa = len(genome_dict)
    nr_genomes = distinct_taxa
    min_taxa = distinct_taxa # Subtract a number here (typically 1-3) to allow trueCOG groups with genes not conserved in each genome
    # Determines whether or not trueCOGs are always considered COGs for the final output
    propagate_truecogs = settings['propagate_truecogs']
    
    if not genomes_allowed:
        genomes_allowed = genome_dict.keys()
    
#    gene2scf,scf_genelist = generate_scaffold_lists(genome_list)
    logger.info('Parsing allvall BLAST results')
    (score_dict_exp) = read_allvall(allvallfile,genome_dict,genomes=genomes_allowed)
    logger.info('Getting bidirectional best hits')
    all_pairs,pairs_per_genome,pair_dict = get_bbh(score_dict_exp,min_simil_bbh)
    logger.info('Getting pairwise truecogs')
    truecog_pairs_per_genome,all_truecog_pairs,true_pair_dict = get_truecog_pairs(pairs_per_genome,all_pairs,genome_dict)
    
    # Now comes the testing of the main group
    # If it doesn't match the criteria, the group will be split
    
    min_truecogs_req = settings['cog_min_truecogs_req']
    min_cutoff_req = settings['cog_min_cutoff_req']
    cutoff_factor = settings['cutoff_factor']
    min_group_size = settings['cog_min_group_size']
#    split_method = settings['cog_split_method']
    split_method = 'bottomup_total' 
    global scores_tested
    scores_tested = {}
    
    logger.debug('Testing main group')
    base_truecogs = make_cog_group(genomes_allowed,truecog_pairs_per_genome,true_pair_dict,genome_dict)
    unique_truecogs = get_unique_truecogs(base_truecogs)
    if len(unique_truecogs) < min_truecogs_req:
        passed = False
        logger.debug('Nr truecogs: %i' %(len(unique_truecogs)))
    else:
        cutoff_dict,true_scores,all_cutoffs,histogram_parameters = calculate_taxon_pairwise(base_truecogs,genome_dict,score_dict_exp,cutoff_mode,cutoff_factor) 
        logger.debug('Nr truecogs: %i\nAverage cutoff: %.3f\nLowest cutoff: %.3f' %(len(unique_truecogs),np.average(all_cutoffs),min(all_cutoffs)))
        if min(all_cutoffs) < min_cutoff_req:
            passed = False
        else:
            passed = True
    if not passed:
        # Proceed with splitting groups
        logger.info('All the genomes together do not form a group that passes COG criteria\nSplitting the group into multiple subgroups')
        
        # Starting with a random genome, genomes are added until the group no longer passes the criterium for minimum nr of trueCOGs
        # This process is repeated multiple times, so that multiple overlapping groups are created which together link up as many genomes as possible
        
        # Genomes can be added randomly (as long as the group fits the requirements), or the "best" genome can be added each time, i.e. the genome 
        # whose addition would result in the highest amount of trueCOGs in the remaining group.
        # The latter may take a long time to calculate if many genomes are analyzed at the same time. 
        
        criterium = 'nr_truecogs'
        add_method = 'random'
        
        all_groups_dict,truecog_dict = bottomup_groups(genomes_allowed,truecog_pairs_per_genome,true_pair_dict,genome_dict,criterium,min_truecogs_req,min_group_size,settings)
        groups_per_genome = get_groups_per_genome(all_groups_dict,genomes_allowed)
        req_groups_per_genome,all_req_keys,genomes_covered_per_genome = required_groups_per_genome(groups_per_genome,all_groups_dict,genomes_allowed)
        all_groups_req_dict = dict([(i,all_groups_dict[i]) for i in all_req_keys])
        keys_req_sorted = sortdictkeysbylenvalues(all_groups_req_dict,reverse=True)
        all_groups = [(all_groups_req_dict[key],truecog_dict[key]) for key in keys_req_sorted]
        
    else:
        all_groups = [[genomes_allowed,base_truecogs]]
        genomes_covered_per_genome = get_genomes_covered_per_genome(all_groups)
    
    # For each subgroup, calculate the pairs, and write out the results
    # Add all to the same COG dictionary
    
    # First remove all groups smaller than allowed
    tbr = []
    for group in all_groups:
        if len(group[0]) < min_group_size:
            tbr.append(group)
    for group in tbr:
        all_groups.remove(group)
        
    i = 0
    analyzed_genomes = []
    all_passed = {}
    logger.info('Further analysis per COG subgroup.')
    for group,base_truecogs in all_groups:
        analyzed_genomes.extend([g for g in group if g not in analyzed_genomes])
        if len(all_groups) == 1:
            j = ''
        else:
            i += 1
            j = '%i' %i
        logger.info('Analyzing group %i (length: %i)' %(i,len(group)))
        # Contains the orthologue cutoff value for each pair of taxons
        taxon_pair_file = os.path.join(path,'taxon_pairwise_stats_%s.%s.txt' %(name,j))
        logger.info('Calculating pairwise cutoffs')
        cutoff_dict,true_scores,all_cutoffs,histogram_data = calculate_taxon_pairwise(base_truecogs,genome_dict,score_dict_exp,cutoff_mode,cutoff_factor,\
                                                                            genomes_allowed=group,outf=taxon_pair_file)
        logger.info('Finding gene pairs based on cutoffs')
        # Convert the truecogs to a set containing all pairs to speed up lookup
        
        truecog_pairs = convert_truecogs_to_pairs(base_truecogs)
        
        # Draw histograms for each group if necessary
        if settings['cog_histograms']:
            for group,truecogs in all_groups:
                histogram_path = os.path.join(path,'COG_histograms_group_%s.png' %j)
                draw_histograms_pairwise(true_scores,histogram_data,histogram_path)

        
        passed = filter_cogs(cutoff_dict,score_dict_exp,base_truecogs,group,propagate_truecogs,truecog_pairs)
        # If two genes are a COG in any collection of genomes, they are considered COGs
        fuse_dict_list(all_passed,passed)
        if split_method != 'bottomup_total':
            # Each group is handled seperately
            # Contains all orthologues (pid above cutoff)
            cog_input = os.path.join(path,"mcl_%s.%s.input" % (name,j))
            # Contains all groups of orthologues as made by mcl
            cog_output = os.path.join(path,"cog_%s.%s.mcl" % (name,j))
            # Information of each group of orthologues
            group_stats_file = os.path.join(path,'group_stats_%s.%s.txt' % (name,j))
            # Information for each gene
            gene_stats_file = os.path.join(path,'gene_stats_%s.%s.txt' % (name,j))
            r1 = scoreDB2file(passed,cog_input)
            r2 = run_mcl(cog_input,cog_output,settings['cores'])
            # For now, get the data for each genome in each of the groups
            # In the case of bottomup_total option, two more options arise
            # a) Use the largest group of each genome, and compare with that group (i.e. process within this code block)
            # b) Use all the information - given a genome, compare all genes in that genome two other genomes, even if they are from different groups
            # b) is used here for now; optionally add an option to allow a user to choose.
            cog_synopsis(cog_output,genome_dict,group_stats_file,gene_stats_file,genomes_covered_per_genome)
    
    del score_dict_exp
    del all_pairs
    del pairs_per_genome
    del pair_dict
    del truecog_pairs_per_genome
    del all_truecog_pairs
    del true_pair_dict
    
    if split_method == 'bottomup_total':
        logger.info('Fusing the subgroups result for MCL')
        # The groups are fused together at this point
        j = 'fused'
        # Contains all orthologues (pid above cutoff)
        cog_input = os.path.join(path,"mcl_%s.%s.input" % (name,j))
        # Contains all groups of orthologues as made by mcl
        cog_output = os.path.join(path,"cog_%s.%s.mcl" % (name,j))
        # Information of each group of orthologues
        group_stats_file = os.path.join(path,'group_stats_%s.%s.txt' % (name,j))
        # Information for each gene
        gene_stats_file = os.path.join(path,'gene_stats_%s.%s.txt' % (name,j))
        r1 = scoreDB2file(all_passed,cog_input)
        r2 = run_mcl(cog_input,cog_output,settings['cores'])
        cog_synopsis(cog_output,genome_dict,group_stats_file,gene_stats_file,genomes_covered_per_genome)
        
    assign_COG_remaining_genes(genome_dict,analyzed_genomes)
    skipped_genomes = [i for i in genomes_allowed if i not in analyzed_genomes]
            
    return(skipped_genomes)



