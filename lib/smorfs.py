# License: GNU Affero General Public License v3 or later

# Functions to find smORFs in sequences and score them
from multiprocessing import Process, Queue
import time
import os
import numpy as np

from features import RiPP
from lib import rev_comp, split_data, load_pickle
from log import return_logger
from Genes import Gene
from multiprocessing_wrappers import operator

logger = return_logger(__name__, False)

table = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
    "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
    "TAT":"Y", "TAC":"Y", "TAA":"STOP", "TAG":"STOP",
    "TGT":"C", "TGC":"C", "TGA":"STOP", "TGG":"W",
    "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
    "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
    "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
    "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",}

def load_SVM_models(settings):
    pickle_path = settings['paths']['pickle']
    SVM3 = load_pickle(os.path.join(pickle_path , 'SVM_poly_3.pkl'))
    SVM4 = load_pickle(os.path.join(pickle_path , 'SVM_poly_4.pkl'))
    SVMr = load_pickle(os.path.join(pickle_path , 'SVM_rbf.pkl'))
    
    return SVM3, SVM4, SVMr
    
def smorfs_operator(settings, genome_dict):
    ripp_features = ['OrganismACC', 'Scaffold', '(Inter)Locustag', 'Start(nt)', 'End(nt)', 'Strand', 'Length', \
    'NtSequence', 'Sequence', 'Cys30', 'Cys20', 'Cys/Ser30', 'Cys/Ser20', \
    'Charge', 'AvgCharge', 'AvgHydrop', 'Length', 'Entropy', 'EntropyRatio', 'A', 'C', \
    'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W',\
    'Y', 'AVIL', 'CGP', 'DE', 'MFYW', 'RHK', 'STNQ', 'SVM_poly4', 'SVM_poly3', 'SVM_rbf', 'SVM_average']
    
    smorf_name = 'smORF'
    threads = settings['cores']
    genome_path = settings['paths']['genome']
    
    SVM3, SVM4, SVMr = load_SVM_models(settings)
    
    jobs = genome_dict.keys()
    target_args = [settings, genome_path, ripp_features, SVM3, SVM4, SVMr]
    results = operator(jobs, threads, smorfs_worker, target_args, quiet=False)
    total = 0
    errors = 0
    genomes_failed = []
    for genome, filtered_smorfs, error, nr_analyzed in results:
        logger.debug('Genome %s: Parsing %s hits' %(genome,len(filtered_smorfs)))
        try:
            genome_dict = parse_hits(filtered_smorfs,genome_dict,genome,smorf_name)
        except:
            genomes_failed.append(genome)
        errors += error
        total += nr_analyzed
    logger.info('Finished smORF SVM scoring. %i out %i gave errors' %(errors, total))
    return(genome_dict, genomes_failed)
    
def smorfs_worker(genome,settings,path,ripp_features,SVM3,SVM4,SVMr):
    genome_path = os.path.join(path, genome)
    # 1) Score the smorfs
    scanned_smorf_path = os.path.join(genome_path, 'SVM_smORFs.out')
    filtered_smorf_path = os.path.join(genome_path, 'SVM_smORFs_filtered.out')   
    
    if os.path.isfile(scanned_smorf_path):
        scanned_smorfs = read_smorfs(scanned_smorf_path,skipfirst=True)
        errors = 0
        total = len(scanned_smorfs)
    else:
        smorfs_genome = read_smorfs(os.path.join(genome_path,'smORFs.txt'))
        scanned_smorfs,errors,total= run_SVM(smorfs_genome,SVM3,SVM4,SVMr)
        write_smorfs(scanned_smorfs,scanned_smorf_path,header=ripp_features)
    
    if os.path.isfile(filtered_smorf_path):
        filtered_smorfs = read_smorfs(filtered_smorf_path,skipfirst=True,d=True,keyindex=2)
    else:
        filtered_smorfs = filter_smorfs_score(scanned_smorfs,settings['svm_cutoff'])
        filtered_smorfs = filter_smorfs_overlap(filtered_smorfs)
        write_smorfs(filtered_smorfs,filtered_smorf_path,header=ripp_features,mult_values=False)
    
    return(genome,filtered_smorfs,errors,total)
    
def filter_smorfs_score(smorfs,threshold=0.99,max_overlap=10):
    # Parses all the significant smORFs and sorts them per tag
    sign_smorfs = []
    for smorf in smorfs:
        avg_score = float(smorf[-1])
        if avg_score > threshold:
            sign_smorfs.append(smorf)
    return(sign_smorfs)

def filter_smorfs_overlap(sign_smorfs):
    '''remove overlapping smORFs on the same strand'''
    # Only remove overlaps if one is an extension of the other
    # I.e. they are on the same frame and one completely overlaps the other
    out_dict = {}
    removed = set()
    for smorf_data in sign_smorfs:
        start = int(smorf_data[3])
        end = int(smorf_data[4])
        for smorf_data2 in sign_smorfs:
            if smorf_data == smorf_data2 or\
               smorf_data[0] != smorf_data2[0] or\
               smorf_data[1] != smorf_data2[1] or\
               smorf_data[5] != smorf_data2[5] or\
               smorf_data[2] in removed:
                continue
            start2 = int(smorf_data2[3])
            end2 = int(smorf_data2[4])
            # Check if 1 completely overlaps 2
            if (start == start2 and end > end2) or (end == end2 and start < start2):
                removed.add(smorf_data2[2])
            # This keeps all the partially overlapping smORFs on different frames
    for smorf_data in sign_smorfs:
        if smorf_data[2] not in removed:
            out_dict[smorf_data[2]] = smorf_data
    return(out_dict)

def parse_hits(hit_smorfs,genome_dict,genome,smorf_name):
    # Parse the smORFs to integrate in the genome_dict
    
    for tag in hit_smorfs:
        smorf = hit_smorfs[tag]
        genome = smorf[0]
        scaffold = smorf[1]
        name = smorf[2]
        start = int(smorf[3])
        end = int(smorf[4])
        strand = smorf[5]
        dnaseq_prot = smorf[7]
        protseq = smorf[8]
        if smorf_name in name:
            smorf_gene = Gene(genome=genome_dict[genome].accession,scaffold=genome_dict[genome][scaffold].name,\
                            start=start,end=end,strand=strand,name=name,transl=protseq,\
                            pseudo=False,SVM_candidate=True,SVM_hit=True,intergenic=True,active=False) #transl=protseq,dna=dnaseq_prot,\
            genome_dict[genome][scaffold][name] = smorf_gene
        else:
            gene = genome_dict[genome][scaffold][name]
            gene.update(SVM_hit=True,intergenic=False)
            genome_dict[genome][scaffold][name] = gene
    return(genome_dict)
    
# Functions to find small CDS in (intergenic) sequences
def smorfs_main(seq, smorf_settings, pos_adj=0):
    #print('smORFs forward')
    smorfs_for = get_smorfs(seq, table, maxlen = smorf_settings['max_len'], minlen = smorf_settings['min_len'], strand = '+')
    #print('smORFs reverse')
    seq_rc = rev_comp(seq)
    smorfs_rev = get_smorfs(seq_rc, table, maxlen = smorf_settings['max_len'], minlen = smorf_settings['min_len'], strand = '-')
    #print('Mirror coordinates on rev. comp. smorfs')
    smorfs_rev_adj = mirror_coordinates(smorfs_rev,len(seq))
    allsmorfs = smorfs_for + smorfs_rev_adj
    allsmorfs_adj = adj_pos(allsmorfs,pos_adj)
    allsmorfs_adj.sort(key=lambda x: x[0])
    #print('SmORFs found: %s' %len(allsmorfs))
    return(allsmorfs_adj)

def get_startpoints(seq):
    start = ['ATG','GTG','CTG','TTG']
    startpoints = []
    for s in start:
        i = 0
        while i > -1:
            i = seq.find(s,i+1)
            if i > -1:
                startpoints.append(i)
    return(startpoints)

def get_smorfs(seq,table,maxlen = 100,minlen = 5,strand = '+'):
    all_aa = []
    startpoints = get_startpoints(seq)
    dna_maxlen = 3*maxlen
    for i in startpoints:
        aa_seq = 'M'
        if len(seq)-i <= dna_maxlen:
            upper = len(seq)-1
        else:
            upper = i+dna_maxlen
        for j in range(i+3,upper,3):
            triplet = seq[j:j+3]
            if len(triplet) == 3:
                if triplet in table:
                    aa = table[triplet]
                else:
                    break
                if aa == 'STOP':
                    if len(aa_seq) >= minlen:
                        aa_seq += '*'
                        dnaseq = seq[i:j+3]
                        all_aa.append((i,j+3,strand,len(aa_seq)-1,dnaseq,aa_seq))
                    break
                else:
                    aa_seq += aa
    return(all_aa)

def mirror_coordinates(all_aa,ls):
    all_aa_rev = []
    for start,end,strand,length,dna,aa in all_aa:
        start_mirror = ls - start
        end_mirror = ls - end
        all_aa_rev.append((end_mirror,start_mirror,strand,length,dna,aa))
    return(all_aa_rev)

def adj_pos(all_aa,adj):
    all_aa_adj = []
    for start,end,strand,length,dna,aa in all_aa:
        all_aa_adj.append((start+adj,end+adj,strand,length,dna,aa))
    return(all_aa_adj)
    
    
# Writing/reading functions
def write_smorfs(smorfs_tagged,outpath,header=False,mult_values=False):
    outfile = open(outpath,'w')  
    if header:
        outfile.write('\t'.join(header)+'\n')
    if type(smorfs_tagged) == dict:
        keys = sorted(smorfs_tagged.keys())
        for key in keys:
            smorfs = smorfs_tagged[key]
            if mult_values:
                for smorf in smorfs:
                    outfile.write('\t'.join([str(i) for i in smorf]) + '\n')
            else:
                outfile.write('\t'.join([str(i) for i in smorfs]) + '\n')
    elif type(smorfs_tagged) == list or type(smorfs_tagged) == tuple:
        for sm in smorfs_tagged:
            outfile.write('\t'.join([str(i) for i in sm]) + '\n')
    outfile.close()

def read_smorfs(path,skipfirst=False,d=False,keyindex=False):
    if d and type(keyindex) == int or type(keyindex) == list:
        out = {}
    else:
        out = []
    f = open(path)
    if skipfirst:
        f.readline()
    for l in f:
        tabs = l.strip().split('\t')
        if type(out) == dict:
            if type(keyindex) == list:
                key = tuple([tabs[keyi] for keyi in keyindex])
            else:
                key = tabs[keyindex]
            if key in out:
                logger.error('Double name found when parsing SVM hits.')
                raise ValueError('Double name found when parsing SVM hits.')
            out[key] = tabs
        else:           
            out.append(tuple(tabs))
    return(out)

# SVM running functions
def run_SVM(smorfs, SVM3, SVM4, SVMr):
    errors = 0
    total = 0
    smorfs_out = []
    for smorf in smorfs:
        seq = smorf[-1].rstrip('*')
        RiPP_object = RiPP(seq)
        RiPP_object.calculate_features()
        features = RiPP_object.get_features()
        m = np.array([np.array(features[1:])])
        # --- run ML
        try:
            AVGS = []
            avgs = []
            for svm in SVM4:
                y = SVM4[svm]['learner'].predict(m)
                avgs.append(y)
            pred = np.mean(avgs)
            AVGS.append(pred)

            avgs = []
            for svm in SVM3:
                y = SVM3[svm]['learner'].predict(m)
                avgs.append(y)
            pred = np.mean(avgs)
            AVGS.append(pred)

            avgs = []
            for svm in SVMr:
                y = SVMr[svm]['learner'].predict(m)
                avgs.append(y)
            pred = np.mean(avgs)
            AVGS.append(pred)
            results = tuple([str(j) for j in list(smorf[:-1])+features+AVGS+[np.mean(AVGS)]])
            smorfs_out.append(results)
        except ValueError: 
            errors += 1
        total += 1
    return(smorfs_out,errors,total)

