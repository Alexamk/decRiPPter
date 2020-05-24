# License: GNU Affero General Public License v3 or later

import time
import os

from subprocess import Popen, PIPE, call
from multiprocessing import Process, Queue

from lib import fuse_dict_list, clean, run_cmd

from log import return_logger
logger = return_logger(__name__, False)

def main(operons,settings): 
    # Split the sequences over fasta files
    nr_procs = settings['cores']
    run_hmm = settings['run_hmm']
    logger.debug('run_hmm = %s' %run_hmm)
    logger.debug('nr_procs = %i' %nr_procs)
    fasta_path = os.path.join(settings['paths']['hmm_basic'],'fasta')
    if not os.path.isdir(fasta_path):
        os.mkdir(fasta_path)
    if nr_procs == 1:
        seqfile = os.path.join(fasta_path,'seqs.fasta')
        operons.write_fasta(filename=seqfile)
    elif nr_procs > 1:
        seqfiles = []
        for n in range(nr_procs):
            seqfiles.append(os.path.join(fasta_path,'seqs%i.fasta' %n))
        operons.write_fasta(filename=seqfiles)
    
    for out_path,db_path in settings['hmm_paths']:
        final_file = out_path + 'result.domtbl'
        if run_hmm != True and os.path.isfile(final_file):
            logger.debug('Parsing results from a previous run')
            results = parse_hmm_domtbl_hmmsearch(final_file)
        elif run_hmm == False and not os.path.isfile(final_file):
            logger.error('No domains found from previous run and hmm is set not to run')
            raise ValueError('No domains found from previous run and hmm is set not to run')
        else:
            if nr_procs == 1:
                results,_,_ = hmmer_worker(seqfile,out_path,db_path,False,'')
            elif nr_procs > 1:
                logger.debug('Running hmmsearch on %i threads' %(nr_procs))
                # Split the files over the workers
                results,tbl_files,domtbl_files = hmmer_operator(seqfiles,out_path,db_path,run_hmm)
                # Concatenate split result files for easier parsing on a second run, and clean up all split files
                clean_n_concatenate(out_path)
        # Resolve overlap between found domains
        domain_dict,location_dict = resolve_overlap(results)    
        # Update the genes with the found domains
        update_operons(operons,domain_dict)
    
    return domain_dict

def hmmer_operator(seqfiles,out_path,db_path,run_hmm):
    processes = []
    q = Queue()
    if not run_hmm: 
        # Parse the result file
        return
    for i in range(len(seqfiles)):
        data = [seqfiles[i],out_path,db_path,q,i]
        p = Process(target = hmmer_worker, args = data)
        processes.append(p)
    for p in processes:
        p.start()
    res = {}
    time.sleep(1)
    domtbl_files = []
    tbl_files = []
    while True:
        while not q.empty():
            temp,tbl_file,domtbl_file = q.get()
            fuse_dict_list(res,temp)
            domtbl_files.append(domtbl_file)
            tbl_files.append(tbl_file)
        processrunning = False
        for i in processes:
            if i.is_alive():
                processrunning=True
        if processrunning:
            time.sleep(1)
        else:
            break
    for p in processes:
        p.join()
    while not q.empty():
        temp,tbl_file,domtbl_file = q.get()
        fuse_dict_list(res,temp)
        domtbl_files.append(domtbl_file)
        tbl_files.append(tbl_file)
    return res,domtbl_files,tbl_files

def hmmer_worker(seqfile,out_path,db_path,q,n):
    base = 'result%s' %n
    tbl_file = out_path + base + '.tbl'
    dom_file = out_path + base + '.domtbl'
    r = run_hmmsearch(seqfile,out_path,1,hmm_db_path=db_path,outfile_base=base)
    out = parse_hmm_domtbl_hmmsearch(dom_file)
    if q:
        q.put([out,tbl_file,dom_file])
    else:
        return([out,tbl_file,dom_file])
    
def run_hmmsearch(in_path,out_path,thr,hmm_db_path,outfile_base='domains'):
    cmds = ['hmmsearch', '-o', out_path+outfile_base + '.out', '--tblout', out_path+outfile_base+'.tbl', '--domtblout', out_path + outfile_base + '.domtbl', \
            '--acc', '--cut_ga', '--cpu', str(thr),hmm_db_path, in_path]
    _ = run_cmd(cmds)
    
def parse_hmm_domtbl_hmmsearch(p):
    # Parse per domain found
    # Sort out overlap later
    with open(p) as f:
        outd = {}
        for _ in range(3):
            l = f.readline()
        for l in f:
            if l.startswith('#'):
                continue
            tabs = l.strip().split(' ')
            tabs = [tab for tab in tabs if tab != '']
            protein_name = tabs[0]
            domain_found = tabs[4]
            if '.' in domain_found:
                domain_found = domain_found.rpartition('.')[0]
            domain_evalue = tabs[12]
            domain_evalue = float(domain_evalue)
            seq_start = int(tabs[17])
            seq_end = int(tabs[18])
            if seq_start > seq_end:
                pass
            if protein_name not in outd:
                outd[protein_name] = []
            outd[protein_name].append([domain_found,seq_start,seq_end,domain_evalue])
    return outd
    
def clean_n_concatenate(path):
    # Cleanup all partial files made by the individual processes
    files_org = os.listdir(path)
    exts = {'.out':False,'.domtbl':True,'.tbl':True}
    for ext in exts:
        concat = exts[ext]
        files = [os.path.join(path,f) for f in os.listdir(path) if f.endswith(ext) and any([i in f for i in '1234567890'])] # only remove the files with numbers in them
        if concat:
            f = files[0]
            f_clean,_,ext = f.rpartition('.')
            f_out = f_clean.rstrip('1234567890') + '.' + ext
            with open(f_out,'w') as handle:
                for f in files:
                    with open(f) as inf:
                        for _ in range(3):
                            l = inf.readline()
                            if files[0] == f:
                                handle.write(l) # Copy the header once
                        for l in inf:
                            if l.startswith('#'):
                                continue
                            else:
                                handle.write(l)
        for f in files:
            pass
            os.remove(f)    
    
def resolve_overlap(domain_dict):
    # Takes in a different dictionary parsed from domtbl hmm file
    # Adds potential pfams one at a time, starting from lowest evalue
    # Overlaps are allowed if the overlap is less than half the non-\
    # overlapping parts for each domain
    location_dict = {}
    prod2domain = {}
    for gene in domain_dict:
        domains = domain_dict[gene]
        # sort by e-value, lowest first
        domains.sort(key=lambda x: x[3]) 
        locations_added = []
        prod2domain[gene] = []
        location_dict[gene] = {}
        for domain,start,end,evalue in domains:
            add = True
            length = end-start+1
            for start_prev,end_prev in locations_added:
                if start >= end_prev or start_prev >= end:
                    # No overlap, add
                    continue
                else:
                    # Get length overlap for both
                    length_prev = end_prev - start_prev + 1 
                    if start >= start_prev:
                        overlap = end_prev - start + 1
                    elif start_prev >= start:
                        overlap = end - start_prev + 1 
                    non_overlap = length - overlap
                    non_overlap_prev = length_prev - overlap
                    if overlap > non_overlap*0.5 or overlap > non_overlap_prev*0.5:
                        add = False
            if add:
                locations_added.append([start,end])
                prod2domain[gene].append(domain)
                if gene not in location_dict:
                    location_dict[gene] = {}
                if domain not in location_dict[gene]:
                    location_dict[gene][domain] = []
                location_dict[gene][domain].append([start,end,evalue])
    return(prod2domain,location_dict)
    
def update_operons(operons,domain_dict):
    c = 0
    g = 0
    added = set()
    for gene in operons.itergenes():
        if not hasattr(gene,'domains'):
            gene.domains = []
        if gene.name in domain_dict and gene.name not in added:
            domains = domain_dict[gene.name]
            gene.domains.extend(domains)
            c += len(domains)
            g += 1
            added.add(gene.name)
    logger.info('Detected %i protein domains in %i gene products' %(c,g))
