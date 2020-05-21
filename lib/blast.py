from lib import fuse_dict_list, parse_fasta, dict_to_fasta, split_dict, run_cmd
from log import return_logger

from multiprocessing import Process, Queue
import os
import time

logger = return_logger(__name__, False)

# Wrapper for DIAMOND BLAST

def make_diamond_database(protein_file,dbfile=False,threads=False):
    if not dbfile:
        name_clean = protein_file[:-6]
        dbfile = '%s.dmnd' %name_clean
    commands = ['diamond','makedb','--in',protein_file,'-d',dbfile]
    if threads:
        commands += ['-p',str(threads)]
    logger.debug(' '.join(commands))
    _ = run_cmd(commands)
    return(dbfile)

def run_diamond(protein_file,database,outname,replacetabs=False,tmpdir='/dev/shm/',maxhits=100,sens=False,moresens=False,threads=False,evalue=False,quiet=False):
    commands = ['diamond','blastp','--query',protein_file,'--db',database,\
    '--max-target-seqs',str(maxhits),'--out',outname,'--tmpdir',tmpdir,'--outfmt',\
    '6', 'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qlen', 'slen', 'evalue', 'bitscore']
    if threads:
        commands += ['-p',str(threads)]
    if evalue:
        commands += ['--evalue',str(evalue)]
    if moresens:
        commands.append('--more-sensitive')
    elif sens:
        commands.append('--sensitive')

    _ = run_cmd(commands)    
    
    if replacetabs:
        f = open(outname)
        fout = open(out_clean + '.tab','w')
        for l in f:
            l_new = l.replace('|','\t')
            fout.write(l_new)
        f.close()
        fout.close()
        outname = out_clean+'.tab'


# Wrapper functions for NCBI BLAST
# Uses it's own multiprocessing by splitting up the queries into multiple files, one for each process

def makeblastdb(infile,dbtype='prot',title = False,outfile=False):
    # Makes NCBI BLAST databases
    current_path = os.getcwd()
    # Changing directory to prevent some errors
    os.chdir(infile.rpartition(os.sep)[0])
    if not outfile:
        outfile = infile.rpartition('.')[0]
    commands = ["makeblastdb","-in",'"%s"' %infile,"-dbtype",dbtype,"-out",outfile]
    if title:
        commands += ["-title",title]
    _ = run_cmd(commands)
    os.chdir(current_path)
    return(outfile)

def blast_worker(query,out,db_file,nr,queue,name,parse_function,parse_function_args,evalue,outfmt,kwargs):
    # Blast this process' part of the query
    logger.debug('Blasting process %s' %nr)
    commands = ['blastp','-query',query,'-db',db_file,'-out',out,'-evalue', str(evalue), '-outfmt',outfmt]
    if kwargs != {}:
        for key,value in kwargs.items():
            commands += ['-'+str(key),str(value)]
    if not os.path.isfile(out):
        _ = run_cmd(commands)
    else:
        logger.debug('%s already found. Not running BLAST' %out)
    if parse_function:
        logger.debug('Parsing function process: %s' %nr)
        if parse_function_args:
            parse_function_args = [out] + list(parse_function_args)
            res = parse_function(*parse_function_args)
        else:
            res = parse_function(out)
    else:
        res = out
    logger.debug('Process %s: blasting done!' %nr)
    queue.put(res)

def blast_mp(path,db_file,input_file,threads,name,outtype,time_pause=60,parse_function=False,parse_function_args=False,evalue=1,\
             outfmt='6 qseqid sseqid pident length mismatch gapopen qlen slen evalue bitscore',db_args={}):
    # Wrapper for the multiprocessing of the BLASTs
    if True:
        logger.debug('Running blast_mp')
    
    split_folder = 'split_files'
    if not split_folder in os.listdir(path):
        os.mkdir(path + split_folder)
    split_path = path + split_folder + os.sep
    
    outfiles = [split_path + '%s_result_part%s.txt' %(name,n) for n in range(threads)]
    queryfiles = [split_path + '%s_query_part%s.txt' %(name,n) for n in range(threads)]    
    # Split the data over the queryfiles if necessary
    seq_dict = parse_fasta(input_file)
    split_dicts = split_dict(seq_dict,threads)
    for i in range(threads):
        text = dict_to_fasta(split_dicts[i],queryfiles[i])
    
    q = Queue()
    processes = []
    for n in range(threads):
        query = queryfiles[n]
        outfile = outfiles[n]
        args = [query,outfile,db_file,n,q,name,parse_function,parse_function_args,evalue,outfmt,db_args]
        p = Process(target=blast_worker,args=args)
        processes.append(p)
    for p in processes:
        p.start()
    res = outtype()
    time.sleep(time_pause)
    while True:
        while not q.empty():
            temp = q.get()
            if outtype == dict:
                fuse_dict_list(res,temp)
            elif outtype == list:
                res.append(temp)
        processrunning = 0
        for i in processes:
            if i.is_alive():
                processrunning+=1
        if processrunning>0:
            time.sleep(time_pause)
        else:
            break
    for p in processes:
        p.join()
    while not q.empty():
        temp = q.get()
        if outtype == dict:
            fuse_dict_list(res,temp)
        elif outtype == list:
            res.append(temp)
    return(res,outfiles)
    
def merge_blast_files(infiles,outfile,remove=False):
    f_out = open(outfile,'w')
    for fil in infiles:
        f = open(fil)
        for l in f:
            f_out.write(l)
        f.close()
        if remove:
            os.remove(fil)
    f_out.close()
     
 
