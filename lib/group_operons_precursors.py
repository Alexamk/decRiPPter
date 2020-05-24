# License: GNU Affero General Public License v3 or later

# BLAST all precursors against each other
# Use them to form networks/groups

import os
import pickle
import time
import numpy as np

from blast import makeblastdb, blast_mp
from lib import store_pickle, load_pickle, run_mcl, read_mcl, pairs_to_groups, write_ssn, clean, muscle, Edge
from Genes import OperonCollection

from multiprocessing_wrappers import operator

from log import return_logger
logger = return_logger(__name__, False)


def main(operons,mibig_dict,settings):
    # Set the name for this collection
    collection_name = 'precursor_groups'
    # BLAST all precursors against each other
    path = settings['paths']['precursor']
    muscle_path = os.path.join(path,'fasta_groups')
    fasta_file_all = path + 'precursors.fasta'
    ssn_file = path+'smORFs_all.ssn'
    ssn_file2 = path+'smORFs_all_filtered.ssn'
    pfile = os.path.join(settings['paths']['pickles'],'smorfs_pairs.pkl') 
    if settings['precursor_blast_run']:
        text,fasta_dict = operons.write_fasta(filename=fasta_file_all,extra_dict=mibig_dict.precursors,SVM_hit=True,core=True)
        db_file = fasta_file_all.rpartition('.')[0]
        if not os.path.isfile(db_file):
            db_file = makeblastdb(fasta_file_all)
        f_args = [settings['prec_min_bitscore'],settings['prec_min_ev']]
        smorfs_pairs, outfiles = blast_mp(path,db_file,fasta_file_all,settings['cores'],'precursors',dict,parse_function=form_pairs,parse_function_args=f_args,time_pause=10) 
        if not os.path.isfile(ssn_file):
            write_ssn(smorfs_pairs,ssn_file,i=0)
        store_pickle(smorfs_pairs,pfile)
    else:
        smorfs_pairs = load_pickle(pfile)
    if settings['precursor_mcl']:
        mcl_file = path + 'smORFs_mcl.txt'
        if not os.path.isfile(mcl_file):
            r1 = run_mcl(ssn_file,mcl_file,settings['cores'])
        smorf_groups = read_mcl(mcl_file,min_group_size=settings['min_group_size'])
        filtered_pairs = filter_smorfs_pairs(smorf_groups,smorfs_pairs)
        if not os.path.isfile(ssn_file2):
            write_ssn(filtered_pairs,ssn_file2,i=0)
    else:
        smorf_groups = pairs_to_groups(smorfs_pairs,settings['min_group_size'])
        
#    clean(os.path.join(path,'split_files'))
    group_names,operons_per_group = assign_groups(operons,smorf_groups,mibig_dict,collection_name)
    operon_collections,_ = make_group_collections(operons,operons_per_group,collection_name,muscle_path)
    muscle_operator(operon_collections,path,settings['cores'])
    t0 = time.time()
    operon_pairs = get_pairing_precursors(operons,mibig_dict,smorfs_pairs)
    original_pairs = operon_pairs
    if settings['precursor_mcl']:
        operon_pairs = get_pairing_precursors(operons,mibig_dict,filtered_pairs)
    t1 = time.time()
    set_precursor_pairs(operons,mibig_dict,operon_pairs)
    t2 = time.time()
    logger.debug('Pairing precursors time: %.2f' %(t1-t0))
    logger.debug('set precursor_pairs time: %.2f' %(t2-t1))
    return group_names,operon_collections,operon_pairs,original_pairs,filtered_pairs,smorf_groups

def form_pairs(allvall_file,bitscore_cutoff,ev_cutoff):
    all_scores = {}
    f = open(allvall_file)
    total = 0
    allowed = 0
    for l in f:
        total += 1
        tabs = l.strip().split('\t')
        name1 = tabs[0]
        name2 = tabs[1]
        pid = float(tabs[2])
        bitscore = float(tabs[-1])
        ev = float(tabs[-2])
        
#        print 'E-value: %s (cutoff: %s)' %(ev,ev_cutoff)
        if name1 == name2 or bitscore <= bitscore_cutoff or ev >= ev_cutoff:
            continue
        else:
            allowed += 1
            if name1 not in all_scores:
                all_scores[name1] = {}
            all_scores[name1][name2] = bitscore,pid
#    print('allowed: %s out of %s' %(allowed,total))
    return(all_scores)
    
def blasthits2muscle(group_names,operons,path,clw=True):
    f = open(path,'w')
    base_cmd = ['muscle']
    out = {}
    if clw:
        base_cmd.append('-clw')
    i = 0
    for name in group_names:
        i += 1
        f.write('Entry %s\n' %i)
        text,smorf_dict = operons.write_fasta(precursor_group=name)
        p = Popen(base_cmd,stdin = PIPE,stdout = PIPE, stderr = PIPE)
        stdout = p.communicate(text)
        f.write(stdout[0] + '\n\n')
        out['entry_%i' %i] = stdout[0]
    f.close()       
    return(out)
    
def assign_groups(operons,smorf_groups,mibig_dict,collection_name):
    # Assign a number to each smorf group
    # Add a variable precursor_i to each smorf
    i = 0
    group_names = []
    operons_per_group = {}
    for group in smorf_groups:
        i += 1
        group_name = '%s_%i' %(collection_name,i)
        group_names.append(group_name)
        operons_per_group[group_name] = []
        for gene in group:
            gene_obj = operons[gene]
            if gene_obj == None:
                # Gene from the MIBiG dict
                bgc,_,prec = gene.partition('_')
                operon = mibig_dict[bgc]
                if operon == None:
                    raise KeyError('%s is not found in operon database or mibig dictionaries' %gene)
            else:
                # Assign the operon to this group
                operon = operons.find_by_precursor(gene_obj)
                if not hasattr(gene_obj,collection_name):
                    setattr(gene_obj,collection_name,group_name)
                else:
                    raise ValueError('%s is part of multiple groups. Likely error in pairs to group function' %gene)
            operons_per_group[group_name].append(operon)
    return group_names,operons_per_group
            
def make_group_collections(operons,operons_per_group,collection_name,muscle_path):
#    # Rework the operons per group to be operoncollection objects
    operons_per_group_coll = {}
    for group in operons_per_group:
        operons_group = operons_per_group[group]
        # Don't use the group if it only contains operons from MIBiG
        if all([hasattr(operon,'mibig') and operon.mibig for operon in operons_group]):
            continue
        operon_collection = OperonCollection(dict([(operon.name,operon) for operon in operons_group]),prep=False,collection_type=collection_name,name=group,\
                                             descr='Precursor BLAST',muscle_path=os.path.join(muscle_path,group+'.clw'))
        operon_collection.set_quick()
        operons_per_group_coll[group] = (operon_collection)
        for operon in operons_group:
            if not hasattr(operon,collection_name):
                setattr(operon,collection_name,[])
            getattr(operon,collection_name).append(group)
    return operons_per_group_coll,operons_per_group

def blasthit2muscle(operon_collections,path):
    path_fasta = os.path.join(path,'fasta_groups')
    path_muscle = os.path.join(path,'all_muscle.txt')
    if not os.path.isdir(path_fasta):
        os.mkdir(path_fasta)
    out_text = ''
    for name,coll in operon_collections.items():
        fasta_file = os.path.join(path_fasta,coll.name+'.fasta')
        muscle_file = os.path.join(path_fasta,coll.name+'.clw')
        if not os.path.isfile(fasta_file):
            coll.write_fasta(filename=fasta_file,precursor_groups=coll.name)
        if not os.path.isfile(muscle_file):
            muscle(fasta_file,muscle_file,clw=True)
        with open(muscle_file) as infile:
            text = infile.read()
        out_text += '%s\n\n%s' %(coll.name,text)
    with open(path_muscle,'w') as handle:
        handle.write(out_text)
#    clean(path_fasta)
    
def muscle_operator(operon_collections,path,nr_cores):
    path_fasta = os.path.join(path,'fasta_groups')
    path_muscle = os.path.join(path,'all_muscle.txt')
    if not os.path.isdir(path_fasta):
        os.mkdir(path_fasta)
    for coll in operon_collections.values():
        fasta_file = os.path.join(path_fasta,coll.name+'.fasta')
        if not os.path.isfile(fasta_file):
            coll.write_fasta(filename=fasta_file,precursor_groups=coll.name)
    
    jobs = [i.name for i in operon_collections.values()]
    worker_args = [path_fasta,path_muscle]
    results = operator(jobs,nr_cores,muscle_worker,worker_args,sleeptime=1)
    out_text = ''.join(results)
    with open(path_muscle,'w') as handle:
        handle.write(out_text)
#    clean(path_fasta)
    
def muscle_worker(name,path_fasta,path_muscle):
    fasta_file = os.path.join(path_fasta,name+'.fasta')
    muscle_file = os.path.join(path_fasta,name+'.clw')
    if not os.path.isfile(muscle_file):
        muscle(fasta_file,muscle_file,clw=True)
    with open(muscle_file) as infile:
        text = infile.read()
    return('%s\n\n%s' %(name,text))
    
def filter_smorfs_pairs(smorf_groups,smorf_pairs):
    filtered = {}
    for group in smorf_groups:
        allowed = set(group)
        for smorf in group:
            if smorf in smorf_pairs:
                for smorf_paired in smorf_pairs[smorf]:
                    if smorf_paired in allowed and smorf != smorf_paired:
                        if smorf not in filtered:
                            filtered[smorf] = {}
                        filtered[smorf][smorf_paired] = smorf_pairs[smorf][smorf_paired]
    return filtered
    
def get_pairing_precursors(operons,mibig_dict,smorf_pairs):
    operon_pairs = {}
    for smorf,paired in smorf_pairs.items():
        gene_obj = operons[smorf]
        if gene_obj == None:
            bgc,_,prec = smorf.partition('_')
            op1 = mibig_dict[bgc]
        else:
            op1 = operons.find_by_precursor(gene_obj)
        for sp in paired:
            gene_obj2 = operons[sp]
            if gene_obj2 == None:
                bgc,_,prec = sp.partition('_')
                op2 = mibig_dict[bgc]
            else:
                op2 = operons.find_by_precursor(gene_obj2)
            if op1 == op2:
                continue
            if op1.name not in operon_pairs:
                operon_pairs[op1.name] = {}
            if op2.name not in operon_pairs[op1.name]:
                operon_pairs[op1.name][op2.name] = {}
            operon_pairs[op1.name][op2.name][smorf,sp] = paired[sp]
    return operon_pairs
    
def set_precursor_pairs(operons,mibig_dict,operon_pairs):
    for op1 in operon_pairs:
        for op2 in operon_pairs[op1]:
            pairs = operon_pairs[op1][op2]
            pids = []
            for pair in pairs:
                data = pairs[pair]
                pids.append(data[1])
            avg_pid = np.average(pids)
            e = Edge(avg_pid,pairs)
            if op1 in operons:
                operons[op1].pairs.precursor[op2] = e
            else:
                mibig_dict[op1].pairs.precursor[op2] = e
            if op2 in operons:
                operons[op2].pairs.precursor[op1] = e
            else:
                mibig_dict[op2].pairs.precursor[op1] = e

