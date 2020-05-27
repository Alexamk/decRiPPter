# License: GNU Affero General Public License v3 or later

#Genome pipeline script
#Part of the decRiPPter pipeline

#Written by Alexander Kloosterman

#This script does the following things:
#    - Download (if necessary) genomes and preprocesses them into genome dictionary, used in the rest of the pipeline.
#    - Find small peptides and score with the Support Vector Machine model, predicting putative precursors.
#    - BLAST all non-intergenic proteins against each other, and run the COG algorithm to determine the number of orthologues.  

#Further details can be found in the README

import time
import os
import sys
import shutil
import argparse
import ConfigParser as configparser
import traceback
import random
import math
    
from lib.lib import load_pickle,store_pickle,find_file,open_genbank,replace_pipe,rev_comp,fuse_dict_add
from lib.log import return_logger, setup_loggers
from lib.multiprocessing_wrappers import operator
from lib.Genes import Gene,Genome,Scaffold,GenomeCollection

import lib.parse_taxonomy as pt
import lib.download_genomes as dg
import lib.smorfs as sm
import lib.prodigal as pr
import lib.blast as blast
import lib.COG as COG
#import lib.split_genomes_mash as mash
    
logger = return_logger('genome_prep', False)
    
def load_cog(genome_dict,pickle_path,f='COG_data.pkl'):
    COG_data = load_pickle(os.path.join(pickle_path,f))
    counter = 0
    for gene_obj in genome_dict.itergenes(intergenic = False):
        data = COG_data.get(gene_obj.name, None)
        if data:
            gene_obj.COG = data['COG']
            gene_obj.COG_total = data['COG_total']
            counter += 1
    logger.debug('Loaded %i COG scores' %counter)

def store_cog(genome_dict,pickle_path,f='COG_data.pkl'):
    COG_dict = {}
    for gene in genome_dict.itergenes(intergenic=False):
        COG_dict[gene.name] = {'COG':gene.COG,'COG_total':gene.COG_total}
    outf = os.path.join(pickle_path,f)
    store_pickle(COG_dict,outf)
    
def store_file_parsed(genome_dict, pickle_path):
    f = 'files_parsed.pkl'
    files_pared = {}
    for genome in genome_dict:
        files_parsed[genome.name] = genome.file_parsed
    
def read_extra_genomes(path):
    f=open(path)
    out = []
    for l in f:
        tabs = l.strip().split('\t')
        out.extend(tabs)
    return(out)
    
def process_in_command(in_command, genome_folder):
    filenames_clean = []
    if os.path.isfile(in_command):
        filename_clean = process_infile(in_command, genome_folder)
        filenames_clean.append(filename_clean)
    elif os.path.isdir(in_command):
        for infile in os.listdir(in_command):
            file_path = os.path.join(in_command, infile)
            if os.path.isfile(file_path):
                filename_clean = process_infile(file_path, genome_folder)
                filenames_clean.append(filename_clean)
    return filenames_clean

def process_infile(infile, out_folder):
    infile_folder, infile_base = os.path.split(infile)
    clean_name, ext = os.path.splitext(infile_base)
    genome_folder = os.path.join(out_folder, clean_name)
    if not os.path.isdir(genome_folder):
        os.mkdir(genome_folder)
    dest_file = os.path.join(genome_folder, infile_base)
    if not os.path.isfile(dest_file):
        shutil.copy(infile, dest_file)
    return clean_name
    
def read_genome_names(path,subset=False):
    folder2name = {}
    for folder in os.listdir(path):
        if subset and folder not in subset:
            continue
        p = os.path.join(path,folder)
        if os.path.isfile(p):
            continue
        files = os.listdir(p)
        if 'genome_name' in files:
            with open(path + folder + os.sep + 'genome_name') as inf:
                name = inf.readline().strip()
            folder2name[folder] = name
        else:
            folder2name[folder] = folder
    return folder2name
    
def filter_accs_files(genome_path,taxids,ext,extra,negative=[]):
    remaining = []
    no_files = []
    for t in taxids:
        path_g = os.path.join(genome_path, t)
        f = find_file(path_g,ext,extra,negative)
        if f:
            remaining.append(t)
        else:
            no_files.append(t)
    return(remaining,no_files)
    
def prep_parse_operator(genomes_annotated, settings, genome_seqs, gene_names, scaffold_names, files_parsed):
    genome_path = settings['paths']['genome']
    cores = settings['cores']
    
    worker_args = [genome_path]
    results = operator(genomes_annotated, cores, prep_parse_worker, worker_args)
    
    for genome, genome_gene_names, genome_scaffold_names, gbk_file, seqs in results:
        files_parsed[genome] = gbk_file
        gene_names = fuse_dict_add(gene_names, genome_gene_names)
        scaffold_names = fuse_dict_add(scaffold_names, genome_scaffold_names)
        genome_seqs[genome] = seqs
        
    return(genome_seqs, gene_names, scaffold_names, files_parsed)
    
def prep_parse_worker(genome, genomes_path):
    path_g = os.path.join(genomes_path, genome)
    gbkfile = find_file(path_g,['.gbk','.gbff'],['genomic'])
    full_path = os.path.join(path_g, gbkfile)
    seqs = open_genbank(full_path)
    gene_names = {}
    scaffold_names = {}
    names_to_ignore = ['none']
    for contig in seqs:
        scaffold = contig.id
        if scaffold == 'unknown':
            scaffold = contig.name
        if scaffold not in scaffold_names:
            scaffold_names[scaffold] = 0
        scaffold_names[scaffold] += 1
        added = set()
        for feat in contig.features:
            if feat.type == 'gene' or feat.type == 'CDS':
                name = None
                if 'locus_tag' in feat.qualifiers:
                    name = feat.qualifiers['locus_tag'][0]
                elif 'gene' in feat.qualifiers:
                    name = feat.qualifiers['gene'][0]
                elif 'product' in feat.qualifiers:
                    name = feat.qualifers['product'][0]
                if name in names_to_ignore:
                    # Ignore these cases
                    continue
                if name == None:
                    nonames += 1
                    continue
                name = name.replace(' ','') # In case the line broke, BioPython adds a space inbetween
                if name not in gene_names:
                    gene_names[name] = 0
                start,end = int(feat.location.start),int(feat.location.end)
                if (start,end) not in added:
                    # Don't add the name twice if the gene was added on the exact same location.
                    # For example, when both a "gene" record and a "CDS" record describe the same gene.
                    # However, only considering CDS may skip older genbank which only contain "gene" records 
                    gene_names[name] += 1
                    added.add((start,end))
                
    return genome, gene_names, scaffold_names, full_path, seqs
           
def gbk_parser(settings, genome_seqs, gene_names, scaf_names, files_parsed, genomeacc2name):
    include_pseudo = settings['pseudo']
    smorf_name_base = settings['smorf_base_name']
    overlap = settings['maximum_overlap_smorf']
    genome_path = settings['paths']['genome']
    
    jobs = gather_parse_gbk_jobs(genome_seqs, files_parsed)
    parse_args = [genome_path, gene_names, scaf_names, include_pseudo, files_parsed, genomeacc2name, smorf_name_base, overlap]
    results = operator(jobs, settings['cores'], parse_gbk_worker, parse_args, quiet=False)
    genome_dict = {}
    newly_named_genes = 0
    for genome, nr_new_genes in results:
        genome_dict[genome.accession] = genome
        newly_named_genes += nr_new_genes
    Genomes = GenomeCollection(genome_dict,active=False)
    return Genomes,newly_named_genes
    
def gather_parse_gbk_jobs(genome_seqs, files_parsed):
    jobs = []
    for genome, file_parsed in files_parsed.items():
        if genome in genome_seqs:
            seqs = genome_seqs[genome]
        else:
            seqs = None
        job = (genome, file_parsed, seqs)
        jobs.append(job)
    return jobs
    
def parse_gbk_worker(genome_data,path,gene_names,scaf_names,include_pseudo,files_parsed,genomeacc2name,smorf_name_base,overlap):
    genome, file_parsed, seqs = genome_data
    genome_folder = os.path.join(path,genome) 
    if seqs == None:
        seqs = open_genbank(os.path.join(genome_folder,file_parsed))
    scafname2realname = {}
    new_smorfs = []
    scaffolds_genome = {}
    newgenenamebase = '%s_gene' %(genome)
    newscaffoldnamebase = '%s_scaffold' %(genome)
    newly_named_scaffolds = 0
    newly_named_genes = 0
    smorf_nr = 1
    records_allowed = ['gene','CDS','tRNA','rRNA'] # Between these records, the DNA is considered intergenic (for the smORFs)
    records_parsed = ['CDS'] # These records are added to the dictionaries
    names_to_ignore = ['none']
    namesreplaced = {}
    name2realname = {}
    
    if genome in genomeacc2name:
        genome_name = genomeacc2name[genome]
    else:
        genome_name = genome
    
    for contig in seqs:
        name_prev = 'lower_edge'
        end_prev = 0
        current_index = -1
        scaffold_name = contig.id
        if scaffold_name == 'unknown':
            scaffold_name = contig.name
        realname = scaffold_name
        if scaf_names[scaffold_name] > 1:
            newscaf = '%s_%s' %(newscaffoldnamebase,newly_named_scaffolds)
            scaffold_name = newscaf
            newly_named_scaffolds += 1
        
        genes_scaffold = {}
        scaffold_seq = str(contig.seq)
        
        for feat in contig.features:
                        
            current_index += 1
            if feat.type in records_allowed:
                is_pseudo = 'pseudo' in feat.qualifiers
                if (feat.type in records_parsed and not (include_pseudo == False and is_pseudo))\
                or  feat.type == 'gene' and is_pseudo and include_pseudo == True:
                    # not parsed unless from the appropriate type (CDS) and not a pseudogene, unless that is allowed as well 
                    # Exception for genes that are annotated as pseudo, since these sometimes do not have a CDS record
                    parse = True
                else:
                    parse = False
                
                start,end = int(feat.location.start),int(feat.location.end)
                
                name,name2realname,newly_named_genes,namesreplaced,renamed = name_gene(feat,names_to_ignore,gene_names,scaf_names,newly_named_genes,namesreplaced,\
                                                                                    newgenenamebase,name2realname,genome,scaffold_name)
                if name == None:
                    continue
                    
                # If it is a gene record and a pseudo gene, check that the next record is not a tRNA/rRNA with the same name. If so, don't parse.
                if include_pseudo and  is_pseudo and feat.type == 'gene':
                    # print('Pseudo gene found - Checking ahead for first allowed record')
#                        print('Current index: %s' %current_index)
                    found = False
                    temp_adj = 1
                    while not found:
#                            print 'Temp_adj: %s' %temp_adj
                        if temp_adj + current_index == len(contig.features):
#                                print('No matching record found - gene is parsed')
                            break
                        next_record = contig.features[current_index + temp_adj]
                        if next_record.type in records_allowed:
                            found = True
                            next_name,_,_,_,_ = name_gene(next_record,names_to_ignore,gene_names,scaf_names,newly_named_genes,namesreplaced,\
                                                                                    newgenenamebase,name2realname,genome,scaffold_name)
#                                print('Next gene name: %s\nNext gene type: %s' %(next_name,next_record.type))
                            if next_name == name and (next_record.type == 'tRNA' or next_record.type == 'rRNA'):
#                                    print('Setting parse to false: rRNA/tRNA pseudo gene')
                                parse = False
                            else:
                                pass
#                                    print('Leaving parse to true')
                            break
                        temp_adj += 1
                    
                if feat.strand == 1:
                    strand = '+'
                elif feat.strand == -1:
                    strand = '-'
                else:
                    logger.critical('Nonlegal strand %s in genome %s' %(strand, genome))
                    raise ValueError('Nonlegal strand %s in genome %s' %(strand, genome))
                    
                 # Get the intergenic region, in which all proteins < 100 aa are found
                 # Overlap refers to the amount of nt that is allowed to overlap with the genes
                     
                if name != name_prev:
                    # print('Analyzing intergenic region')
                    start_intergenic_seq = pos_adj = max(0,end_prev - overlap)
                    end_intergenic_seq = min(start+overlap,len(contig.seq))
                    intergenic_seq = str(contig.seq[start_intergenic_seq:end_intergenic_seq])
                    smorfs_intergenic = sm.smorfs_main(intergenic_seq,pos_adj=pos_adj)
                    
                    for smorf_data in smorfs_intergenic:
                        smorf_name = '%s_%s_%s' %(genome,smorf_name_base,str(smorf_nr).zfill(6))
                        smorf_nr += 1 
                        new_smorfs.append((genome,scaffold_name,smorf_name) + smorf_data)
                    
                if parse:
                    # Add the sequence to the protseq dict, name2accesion dict, and possibly to the other smorfs
                    # print('Parsing gene')
                    if 'translation' in feat.qualifiers:
                        protseq = feat.qualifiers['translation'][0]
                    else:
                        logger.debug('No translation found!')
                        protseq = 'M' + str(feat.extract(contig).seq.translate())[1:]
                    dnaseq_prot = scaffold_seq[start:end]
                    if strand == '-':
                        dnaseq_prot = rev_comp(dnaseq_prot)
                    gene = Gene(start=start,end=end,strand=strand,name=name,transl=protseq,\
                            pseudo=is_pseudo,SVM_candidate=(len(protseq)<100 and not is_pseudo),intergenic=False,SVM_hit=False,active=False) 
                    if 'protein_id' in feat.qualifiers:
                        gene.protein_id = feat.qualifiers['protein_id'][0]
                    if 'locus_tag' in feat.qualifiers:
                        gene.locus_tag = feat.qualifiers['locus_tag'][0]
                    if 'gene' in feat.qualifiers:
                        gene.gene = feat.qualifiers['gene'][0]
                    if 'gene_synonym' in feat.qualifiers:
                        gene.gene_synonym = feat.qualifiers['gene_synonym']
                        
                    # Set the realname (if the gene was renamed)
                    if renamed:
                        gene.realname = name2realname[name]
                        
                    # Pick a representative name for display later
                    if hasattr(gene,'locus_tag'):
                        gene.repr_name = gene.locus_tag
                    elif hasattr(gene,'gene'):
                        gene.repr_name = gene.gene
                    else:
                        gene.repr_name = gene.name
                        
                    if feat.type == 'CDS' or gene.name not in genes_scaffold:
                        # Overwrite the record if a CDS is parsed, as these usually hold better information about amino acid sequences (in the case of pseudo genes)
                        genes_scaffold[gene.name] = gene
                
                # Always update the end_prev, so that smorfs are not found overlapping with tRNAs/rRNAs
                # Also always update the name_prev to prevent double-taking on intergenic sequences
                name_prev = name
                end_prev = end
        scaffold = Scaffold(genes_scaffold,name=scaffold_name,genome=genome,seq=str(contig.seq),realname=realname,active=False)
        for gene in scaffold.itergenes():
            gene.scaffold = scaffold.name
        scaffolds_genome[scaffold.name] = scaffold
    genome_obj = Genome(scaffolds_genome,accession=genome,name=genome,descr=genome_name,path=genome_folder,active=False)
#    if genome in files_parsed:
#        genome_file = files_parsed[genome]
#        genome_obj.parsed_file = genome_file
    for gene in genome_obj.itergenes():
        gene.genome = genome_obj.accession
        if gene.SVM_candidate:
            new_smorfs.append((gene.genome,gene.scaffold,gene.name,gene.start,gene.end,gene.strand,len(gene.transl),gene.get_seq(reference_obj=genome_obj),gene.transl))
    smorfs_path = os.path.join(genome_folder,'smORFs.txt') 
    sm.write_smorfs(new_smorfs,smorfs_path)
    return genome_obj,newly_named_genes           
           
def split_genomes_randomly(genome_dict,settings,nr_proteins):
    # Split genomes evenly over estimated amount of groups
    avg_per_genome = float(nr_proteins)/len(genome_dict)
    max_prot = settings['max_proteins']
    nr_genomes_per_group = max_prot / avg_per_genome
    nr_groups = int(math.ceil(len(genome_dict) / nr_genomes_per_group))
    groups_out = [[] for _ in range(nr_groups)]
    genomes = genome_dict.keys()
    random.shuffle(genomes)
    i = 0
    for genome in genomes:
        groups_out[i].append(genome)
        i += 1
        if i == len(groups_out):
            i = 0
    return groups_out
    
def name_gene(feat,names_to_ignore,gene_names,scaffold_names,newly_named_genes,names_replaced,newgenenamebase,name2realname,genome,sc):
    name = 'Not_found'
    renamed = False
    start,end = int(feat.location.start),int(feat.location.end)
    if 'locus_tag' in feat.qualifiers:
        name = feat.qualifiers['locus_tag'][0]
    elif 'gene' in feat.qualifiers:
        name = feat.qualifiers['gene'][0]
    if name in names_to_ignore:
        # Ignore these cases
        name = None
    illegal = '| '
    if name != None:
        for c in illegal:
            name = name.replace(c,'')
    if gene_names[name] > 1 or name in scaffold_names or name == 'Not_found': # in the last case, no /gene entry or /locus_tag entry was found
        if 'gene_synonym' in feat.qualifiers:
            name = feat.qualifiers['gene_synonym'][0]
            # Split by ';' in case of multiple synonyms
            if ';' in name:
                name = name.split(';')[0]
        else:
            if name in names_replaced and (start,end,genome,sc) in names_replaced[name]:
                newname = names_replaced[name][(start,end,genome,sc)]
                logger.debug('Not replacing name %s, already replaced with %s before (same start/end)' %(name,newname))
                name = newname
            else:
                newname = '%s_%i' %(newgenenamebase,newly_named_genes)
                newly_named_genes += 1
                name2realname[newname] = name
                if name not in names_replaced:
                    names_replaced[name] = {}
                names_replaced[name][(start,end,genome,sc)] = newname
                logger.debug('Renaming %s to %s (Not a unique name, and no alternatives found)' %(name,newname))
                name = newname
                renamed = True
    # Remove spaces as they can break BLAST output
    if name != None:
        name = name.replace(' ','')
    
    return(name,name2realname,newly_named_genes,names_replaced,renamed)

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
    script_path = os.path.realpath(sys.argv[0])
    folder = os.path.dirname(script_path)
    config_file = os.path.join(folder,'config.ini')
    settings = read_config(config_file)
    
    parser = argparse.ArgumentParser()

    parser.add_argument('name',metavar='PROJECT NAME',type=str,help='A name for the project to run')

    base_group = parser.add_argument_group('Basic options')
    base_group.add_argument('-o','--outputfolder',help='path where the results folder will be created',default=os.getcwd())
    base_group.add_argument('-tf','--taxonomy_folder',help=argparse.SUPPRESS)
    base_group.add_argument('-c','--cores', help='number of cores to be used on this run',type=int)
    base_group.add_argument('-i', '--in', help='Input file or folder', default=False)
    base_group.add_argument('-v', '--verbosity', help='Verbosity (0-3), higher means more output', choices = [0, 1, 2, 3], type=int, default=2)

    download_group = parser.add_argument_group('NCBI download options')
    download_group.add_argument('-t','--taxid', help='specify a taxid from which to download genomes',default=0,type=int)
    download_group.add_argument('--assembly_level',help=argparse.SUPPRESS,action='append')
    download_group.add_argument('--refseq_category',help=argparse.SUPPRESS,action='append')
    
    # further specification of which genomes to use
#    filter_group = parser.add_mutually_exclusive_group()
    # 'Genome selection options','These genomes may be a subset of downloaded/found genomes or user-provided in the genomes folder'
    download_group.add_argument('-rg','--reuse_genomes', 
    help='Used to reuse genome files already in the Genomes folder. Specify a file that contains a list of genomes to use exclusively, or use "all" to use all genomes in the Genome folder',default=False)
    
    # parsing/prodigal arguments
    parse_group = parser.add_argument_group('Parsing options')
    parse_group.add_argument('--run_prodigal',help='set how and whether to run prodigal (default: tries to annotate .fna files if .gbk files are missing)',\
                             choices=['auto','never','always'],dest='run_prodigal')
                             
    parse_group.add_argument('-p','--pseudo',help='Toggles the detection of pseudo genes (default on)',action='store_false')
    parse_group.add_argument('--load_pickles', help='instead of parsing files, load previously parsed pickles (use if parsing was done before)',action='store_true',default=False)
    download_group.add_argument('--store_cog',help='Stores the COG data seperately',\
                                              default=False,action='store_true')
    download_group.add_argument('--load_cog',help='Load the COG data from a stored file',default=False,action='store_true')

    # To prepare for COG
    split_group = parser.add_argument_group('BLAST/COG preparation option')
    split_group.add_argument('--skip_blast',help=argparse.SUPPRESS,choices=[True,False,'auto'])
    split_group.add_argument('--max_proteins',help=argparse.SUPPRESS,type=int)
    split_group.add_argument('--min_group_size',help=argparse.SUPPRESS,type=int)
    split_group.add_argument('--split_genomes', help=argparse.SUPPRESS, choices=['random','load']) 

    args = parser.parse_args()
    for key,value in args.__dict__.items():
        # Update only if the args value is not None
        if value != None:
            settings[key] = value
            
    return settings,args
    
def check_args_compatibility(settings):
    # At least some input is required
    if not settings['in'] and not settings['taxid'] and not settings['reuse_genomes']:
        logger.error('No input was given in any way')
        raise ValueError('No input was given in any way')
    if settings['tax_file'] not in ['refseq', 'genbank']:
        logger.error('No valid parse_file parameter given')
        raise ValueError('No valid parse_file parameter given')
    
def check_data_folder(data_folder):
    required_files = [
                    'RiPP_dataset_xlsx',
                    ]
    check_files(data_folder, required_files)
                    
                    
def check_pickle_folder(pickle_folder):
    required_files = [
                    'SVM_rbf.pkl',
                    'SVM_poly_4.pkl',
                    'SVM_poly_3.pkl',
                     ]
    check_files(pickle_folder, required_files)
      
def check_files(folder, required_files):
    for required_file in required_files:
        if not os.path.isfile(os.path.join(pickle_folder, required_file)):
            logger.error('File %s not found in folder %s' %(required_file, folder))
            raise ValueError('File %s not found in folder %s' %(required_file, folder))
    
def check_enough_genomes(genomes_to_analyze, genomes_not_processed, settings):
    min_genomes = settings['cog_min_group_size']
    nr_analyzed = len(genomes_to_analyze) - len(genomes_not_processed)
    if nr_analyzed < min_genomes:
        logger.error('%i genomes going to be analyzed with a minimum group size of %i' %(nr_analyzed, min_genomes))
        raise ValueError('%i genomes going to be analyzed with a minimum group size of %i' %(nr_analyzed, min_genomes))
        
        
if __name__ == '__main__':
    t_start = time.time()
    
    settings,args = parse_arguments()
    
    name = settings['name']
    base_path = os.path.abspath(os.path.join(settings['outputfolder'], name) + os.sep)
    
    check_args_compatibility(settings)

    if not os.path.isdir(base_path):
        os.mkdir(base_path)
        
    setup_loggers(settings)
    
    tax_path = settings['taxonomy_folder']
    
    genome_path = os.path.join(base_path, 'Genomes')
    COG_path_base = os.path.join(base_path, 'COG')
    data_path = os.path.join(base_path, 'data')
    pickle_path = os.path.join(data_path, 'Pickles')
        
    script_path = os.path.realpath(sys.argv[0])
    script_folder = os.path.dirname(script_path)
    base_data_path = os.path.join(script_folder,'data')
    
    paths = dict(genome = genome_path,
                 COG = COG_path_base,
                 data = data_path,
                 pickle = pickle_path,
                 )
    settings['paths'] = paths
    
    # Preparing the folders if necessary
    present_folder = os.listdir(base_path)
    required_folders = ['Genomes','COG']
    
    logger.debug('Copying files and creating folders')

    for req in required_folders:
        if req not in present_folder:
            os.mkdir(os.path.join(base_path,req))
    if 'data' not in present_folder:
        shutil.copytree(base_data_path,data_path.rstrip(os.sep)) 
        #The data folder with requirements is in the same folder as the original script
    else:
        files = os.listdir(data_path)
        req_files = os.listdir(base_data_path)
        for req in req_files:
            if req not in files:
                from_path = os.path.join(base_data_path,req)
                to_path = os.path.join(pickle_path,req)
                logger.debug('Copying %s to %s' %(from_path, to_path))
                if os.path.isfile(os.path.join(base_data_path,req)):
                    shutil.copyfile(from_path, to_path)
                elif os.path.isdir(os.path.join(base_data_path,req)):
                    shutil.copytree(from_path, to_path)

    tax_file = settings['tax_file']
    skip_blast = settings['skip_blast']
    blast_type = settings['blast_type']
    
    run_prodigal = settings['run_prodigal']
    taxid=int(settings['taxid'])
    
    include_pseudo=settings['pseudo']
    
    if settings['load_pickles']:
        logger.info('Loading previously parsed genome_dict pickle')
        genome_dict = load_pickle(os.path.join(pickle_path,'genome_dict_%s.pkl' %name))
        genome_dict.switch_active()
    else:
        genomes_to_analyze = []
        genome_acc_to_name = {}
        if settings['in']:
            logger.info('Copying genomes from --in command')
            filenames_clean = process_in_command(settings['in'], genome_path)
            genomes_to_analyze.extend(filenames_clean)

        if taxid != 0:
            # Download genomes from the chosen taxid
            logger.info('Taxid given - downloading genomes')
            logger.debug("Parsing taxonomy nodes.dmp...")
            if not os.path.isfile(os.path.join(tax_path,'nodes.dmp')):
                logger.info('Taxonomy files not found. Downloading...')
                dg.download_taxonomy_files(tax_path)
                
            id2parent,parent2ids,id2rank = pt.parse_nodes(os.path.join(tax_path,'nodes.dmp'))
            logger.info("Parsing assembly file...")
            if tax_file == 'refseq':
                assembly_file = os.path.join(tax_path,'assembly_summary_refseq.txt')
            elif tax_file == 'genbank':
                assembly_file = os.path.join(tax_path,'assembly_summary_genbank.txt')

            if not os.path.isfile(assembly_file):
                dg.download_assembly_file(tax_path,type=tax_file)
            acc2data,name2acc,genomeacc2name_dl = pt.parse_assembly(assembly_file)
            logger.info("Scanning for full genome of relevant taxid...")
            genomes_to_download, others, notin_taxonomy = pt.find_taxid_in_parent(acc2data,taxid,id2parent,\
                                                    assembly_level_req=settings['assembly_level_req'],refseq_category_req=settings['refseq_category'])
            logger.info('Downloading %s genomes' %len(genomes_to_download))
            downloaded = dg.main(genome_path, run_prodigal, genomes_to_download, acc2data)
            logger.info('Downloaded %i genomes' %(len(downloaded)))
            if len(downloaded) < len(genomes_to_download):
                logger.warning('Failed to download %i genomes!' %(len(genomes_to_download) - len(downloaded)))
            genomes_to_analyze.extend(downloaded)
            genome_acc_to_name.update(genomeacc2name_dl)
            
        if settings['reuse_genomes']:
            # Uses genomes from a chosen list (doesn't download yet - mostly used for a second run or for user genomes)
            # Genomes must be in the Genomes folder
            # In case of user provided genomes, the name of the genome will be the folder name, unless a file is present
            # in that folder named "genome_name", that contains the genome name
            logger.info('Adding genomes from --reuse-genomes option: %s' %settings['reuse_genomes'])
            if settings['reuse_genomes'] == 'all':
                genomes_added = os.listdir(genome_path)
            else:
                genomes_added = read_extra_genomes(settings['reuse_genomes'])
            for genome in genomes_added:
                if genome not in genomes_to_analyze:
                    genomes_to_analyze.append(genome)

        logger.info('Selecting genomes with relevant files...')
        if run_prodigal == 'always':
            genomes_prodigal_processed, not_processed = filter_accs_files(genome_path,genomes_to_analyze,['_prodigal.gbk'],[])
            genomes_dna_files, genomes_not_processed = filter_accs_files(genome_path,not_processed,['fna'],['genomic'])
            logger.info('Files previously annotated with prodigal found for %i out of %i genomes' %(len(genomes_prodigal_processed), len(genomes_to_analyze)))
            logger.info('DNA fasta files found for %i out of %i remaining genomes' %(len(genomes_dna_files), len(not_processed)))
            prodigal_genomes = genomes_dna_files
            genomes_annotated = genomes_prodigal_processed
        else:
            genomes_annotated, genomes_no_gbk = filter_accs_files(genome_path,genomes_to_analyze, ['gbk','gbff'], ['genomic'], ['prodigal'])
            logger.info('Genbank files found for %i out of %i genomes' %(len(genomes_annotated), len(genomes_to_analyze)))
            if run_prodigal == 'auto':
                genomes_dna_files, genomes_not_processed = filter_accs_files(genome_path, genomes_no_gbk, ['fna'], ['genomic'])
                logger.info('DNA fasta files found for %i out of %i remaining genomes' %(len(genomes_dna_files), len(genomes_no_gbk)))
                prodigal_genomes = genomes_dna_files
            elif run_prodigal == 'never':
                prodigal_genomes = []
                genomes_not_processed = genomes_no_gbk
                
        check_enough_genomes(genomes_to_analyze, genomes_not_processed, settings)    
        
        genomes_with_files = genomes_annotated + prodigal_genomes
        logger.info("%i genomes are going to be parsed/processed" %len(genomes_with_files))
                
        if run_prodigal != 'never' and len(prodigal_genomes) != 0: 
            # Function to determine which genomes need prodigal
            logger.info('Annotating DNA fasta files with prodigal')
            genome_seqs_prod, gene_names_prod, scaf_names_prod, files_parsed_prod = pr.fasta_to_gbk_parser(prodigal_genomes, genome_path, settings)
            # Remaining genomes
            logger.debug('Prodigal genomes added')
        else:
            genome_seqs_prod = {}
            scaf_names_prod = {}
            files_parsed_prod = {}
            gene_names_prod = {}
        if genomes_annotated != []:
            # Prepwork only for genomes not just annotated with prodigal, those are already parsed
            t0 = time.time()
            logger.info('Parse preparation of %i genbank files to identify non-unique names' %(len(genomes_annotated)))
            genome_seqs, gene_names, scaf_names, files_parsed = prep_parse_operator(genomes_annotated, settings, genome_seqs_prod, gene_names_prod, scaf_names_prod, files_parsed_prod)
            t1 = time.time()
            logger.info('Time prep: %s' %(t1-t0))
        else:
            genome_seqs = genome_seqs_prod
            gene_names = gene_names_prod
            scaf_names = scaf_names_prod
            files_parsed = files_parsed_prod
        # Adding prodigals parsed genomes
            
        logger.info('Starting the parsing of %i genomes' %(len(genome_seqs)))
        genome_dict, newly_named_genes = gbk_parser(settings, genome_seqs, gene_names, scaf_names, files_parsed, genome_acc_to_name)
        
        logger.debug('Renamed %s genes\n' %(newly_named_genes))
        
        logger.info('Analysing smORFs and assigning SVM score\n')
        genome_dict,failed = sm.smorfs_operator(settings, genome_dict)
        t_sa_start = time.time()
        genome_dict.set_quick()
        genome_dict.switch_active()
        genome_dict.prep()
        t_sa_end = time.time()
        logger.info('Prepped genomes in %.2f seconds' %(t_sa_end - t_sa_start))
    
    # Determine the way to split up the groups, if at all
    nr_proteins = len([i for i in genome_dict.itergenes(intergenic=False)])
    logger.debug('Total non-intergenic proteins: %i' %nr_proteins)
    logger.debug('Maximum allowed proteins: %i' %settings['max_proteins'])
    
    # Only necessary when calculating cog scores
    if settings['load_cog']: 
        load_cog(genome_dict, pickle_path)
    else: 
        if nr_proteins > settings['max_proteins']:
            # Split up genomes
            logger.info('Total number of proteins exceeds max_proteins. Genomes will be split before BLASTing') 
            if settings['split_genomes']  == 'random':
                # Splitting randomly in a given amount of groups
                logger.debug('Splitting up genomes randomly')
                repr_groups = split_genomes_randomly(genome_dict,settings,nr_proteins)
                # Storing the mixup for any next attempts or later reference
                store_pickle(repr_groups,'repr_groups_split_rand.pkl',pickle_path)
            elif settings['split_genomes']  == 'load':
                logger.info('Loading random splitup from previous attempt')
                repr_groups = load_pickle(os.path.join(pickle_path, 'repr_groups_split_rand.pkl'))
        else:
            logger.info('Number of proteins below max_proteins. Genomes will not be split up')
            repr_groups = [genomes_with_files]
        
        if settings['largest_group_only']:
            logger.debug('Using only the largest group of genomes')
            longest = []
            for group in repr_groups:
                if len(group) > len(longest):
                    longest = group
            repr_groups = [longest]
    
        # Determine COG scores within each group
        i = 1
        all_skipped = []
        logger.info('Starting further analysis per group. Number of groups: %s' %len(repr_groups))
        for group in repr_groups:
            logger.info('Group number: %s, size: %s' %(i,len(group)))
            if len(repr_groups) > 1:
                appendix = str(i) + '_' + str(len(group))
                COG_folder = 'COG%s' %appendix
                COG_path = os.path.join(COG_path_base, COG_folder)
                if not os.path.isdir(COG_path):
                    os.mkdir(COG_path)
            else:
                COG_path = COG_path_base
            all_proteome_path = os.path.join(COG_path, 'All_proteins.fasta')
            allvall_name = os.path.join(COG_path, 'allVall_out_%s.tab' %name)
            
            logger.info('Starting allVall BLASTs')
            if not os.path.isfile(all_proteome_path):
                _ = genome_dict.write_fasta(filename=all_proteome_path,genomes=group,intergenic=False)

            maxtargets = max(500, 10*len(group))
            if blast_type == 'NCBI':
                ncbi_db = all_proteome_path.rpartition('.')[0]
                required_files = [ncbi_db + '.pin',ncbi_db + '.phr',ncbi_db + '.psq',allvall_name]
                if skip_blast == False or (skip_blast == 'auto' and not all([os.path.isfile(f) for f in required_files])):
                    logger.info('Making NCBI BLAST database')
                    ncbi_db = blast.makeblastdb(all_proteome_path)
                    logger.info('Running allvall NCBI BLAST')
                    res,outfiles = blast.blast_mp(COG_path,ncbi_db,all_proteome_path,settings['cores'],name,list,parse_function=replace_pipe,\
                                                  db_args=dict(num_alignments=maxtargets))
                    allvall_file = blast.merge_blast_files(res,allvall_name,remove=True)
                else:
                    logger.debug('Skipping BLAST')
            elif blast_type == 'Diamond':
                name_clean = all_proteome_path.rpartition('.')[0]
                dmnd_database = name_clean + '.dmnd'
                logger.info('Making DIAMOND database')
                blast.make_diamond_database(all_proteome_path,threads=settings['cores'])
                if skip_blast == False or not os.path.isfile(allvall_name):
                    logger.info('Running allvall DIAMOND')

                    blast.run_diamond(all_proteome_path, dmnd_database, allvall_name, 
                                    maxhits=maxtargets, moresens=True, threads=settings['cores'], evalue=1)
                else:
                    logger.debug('Skipping BLAST')

            logger.info('Running COG.main...')
            skipped_genomes = COG.main(COG_path, allvall_name, name, settings, genome_dict, group)
            all_skipped.extend(skipped_genomes)
            if len(skipped_genomes) > 0:
                logger.warning('Forced to skip %i out of %i genomes for this group' %(len(skipped_genomes), len(group)))
            i += 1
        if len(all_skipped) > 0:
            logger.warning('Total genomes skipped: %s' %len(all_skipped))
            if len(all_skipped) == len(genome_dict):
                logger.error('No viable genome groups found. Try again with more lenient trueCOG requirements')
                raise ValueError('No viable genome groups found. Try again with more lenient trueCOG requirements')
        genome_dict.remove(all_skipped)
        
    if settings['store_cog']:
        store_cog(genome_dict,pickle_path)
    store_pickle(files_parsed, 'gbk_files_for_antismash.pkl', pickle_path)
    genome_dict.switch_inactive()
    delattr(genome_dict,'quick_dict')
    store_pickle(genome_dict,'genome_dict_%s.pkl' %name,pickle_path)
    t1 = time.time()
    logger.info('Genome prep done. Total time: %s' %(t1-t_start))
    exit()
    
    
    
    
    
