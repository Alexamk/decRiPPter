import os

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature,FeatureLocation
from Bio.Alphabet import IUPAC
from subprocess import call
from multiprocessing import Queue, Process
from multiprocessing_wrappers import operator

from lib import find_file,fuse_dict_add, run_prodigal_cmd
from log import return_logger
from multiprocessing_wrappers import operator

logger = return_logger(__name__, False)

def process_from_fasta(genomes,path,override_prodigal=False):
    # Run prodigal for each genome
    # Parse the fasta via biopython
    # Make a feature list based on prodigals output and expand the seqrecord with it
    # Pass it onto regular genbank parsing
    genome_seqs = {}
    all_names = set()
    gene_nr = 1
    
    for genome in genomes:
        genome_path = os.path.join(path, genome)
        sco_file = os.path.join(genome_path, genome + '.sco')
        fasta_file = find_file(genome_path,['fna','fasta'],['genomic'])
        fasta_file = os.path.join(genome_path, fasta_file)        
        if override_prodigal or not os.path.isfile(sco_file):
            # Always run if override_prodigal is set to do so, else only run if no sco file is found
            run_prodigal_cmd(fasta_file,sco_file)
        
        sr = SeqIO.parse(fasta_file,'fasta')
        features,gene_nr,names = parse_prodigal(sco_file,prev_nr=gene_nr)
        all_names = fuse_dict_add(all_names,names)
        all_seqs = []
        for record in sr:
            name = record.id
            scaffold_features = features[name]
            record.features = scaffold_features
            all_seqs.append(record)
        genome_seqs[genome] = all_seqs
    # Get all scaffold names and gene names so that they can be checked for double names laters
    scaffold_names = get_scaffold_names(genome_seqs)
    # Write the files to a genbank file; also adds translations to the features
    files_parsed = write_genbank(genome_seqs,path)
    logger.debug('Prodigal parsing: renamed %i genes' %(gene_nr-1))
    return(genome_seqs,all_names,scaffold_names,files_parsed)
    
def fasta_to_gbk_parser(genomes,path,settings):
    logger.debug('Starting prodigal annotation of %i genomes, split over %i processes' %(len(genomes),settings['cores']))
    worker_args = [settings,path]
    
    results = operator(genomes,settings['cores'],prodigal_worker,worker_args)
    files_parsed = {}
    genome_seqs = {}
    all_names = {}
    scaffold_names = {}
    for genome, seqs, file_parsed, gene_names, scaf_names in results:
        genome_seqs[genome] = seqs
        files_parsed[genome] = file_parsed
        all_names = fuse_dict_add(all_names,gene_names)
        scaffold_names = fuse_dict_add(scaffold_names,scaf_names)

    return(genome_seqs, all_names, scaffold_names, files_parsed)

def prodigal_worker(genome,settings,path):
    genome_path = os.path.join(path, genome)
    sco_file = os.path.join(genome_path, genome + '.sco')
    fasta_file = find_file(genome_path,['fna','fasta'],['genomic'])
    fasta_file = os.path.join(genome_path, fasta_file)   
    out_genbank = os.path.join(genome_path, genome + '_prodigal.gbk')        
    if settings['override_prodigal'] or not os.path.isfile(sco_file):
        # Always run if override_prodigal is set to do so, else only run if no sco file is found
        run_prodigal_cmd(fasta_file,sco_file)
    
    sr = SeqIO.parse(fasta_file, 'fasta')
    features,gene_nr,gene_names = parse_prodigal(sco_file, start_nr=1, prefix=genome+'_')
    all_seqs = []
    scaf_names_count = {}
    for record in sr:
        renamed = False
        scaffold_name = record.id
        if scaffold_name not in scaf_names_count:
            scaf_names_count[scaffold_name] = 0
        scaf_names_count[scaffold_name] += 1
        scaffold_features = features[scaffold_name]
        record.features = scaffold_features
        all_seqs.append(record)
    if not os.path.isfile(out_genbank):
        write_genbank_single(all_seqs,out_genbank)
    return(genome, all_seqs, out_genbank, gene_names, scaf_names_count)
    
def get_scaffold_names(genome_seqs):
    scaffold_names = set()
    for genome in genome_seqs:
        seqs = genome_seqs[genome]
        for contig in seqs:
            scaffold = contig.id
            if scaffold == 'unknown':
                scaffold = contig.name
            scaffold_names.add(scaffold)
    return scaffold_names
    
def check_scaffolds_names(genome_seqs):
    scaffold_names = set()
    scaf_names_double = set()
    for genome in genome_seqs:
        seqs = genome_seqs[genome]
        for contig in seqs:
            scaffold = contig.id
            if scaffold == 'unknown':
                scaffold = contig.name
            if scaffold in scaffold_names:
                scaf_names_double.add(scaffold)
            else:
                scaffold_names.add(scaffold)
    return scaffold_names,scaf_names_double

def write_genbank(genome_seqs,path):
    gbk_files = {}
    for genome in genome_seqs:
        genome_path = path + genome + os.sep
        out_genbank = genome_path + genome + '_prodigal.gbk'
        if os.path.isfile(out_genbank):
            continue
        seqs = genome_seqs[genome]
        for seq in seqs:
            str_sequence = str(seq.seq)
            new_sequence = Seq(str_sequence,IUPAC.ambiguous_dna)
            seq.seq = new_sequence
            for feat in seq.features:
                if feat.type == 'CDS':
                    dna_feature = feat.extract(seq).seq
                    translation = dna_feature.translate()
                    feat.qualifiers.update(dict(translation=translation))
        SeqIO.write(seqs,out_genbank,'genbank')
        gbk_files[genome] = out_genbank
    return gbk_files
    
def write_genbank_single(seqs,gbk_path):
    for seq in seqs:
        str_sequence = str(seq.seq)
        new_sequence = Seq(str_sequence,IUPAC.ambiguous_dna)
        seq.seq = new_sequence
        for feat in seq.features:
            if feat.type == 'CDS':
                dna_feature = feat.extract(seq).seq
                translation = dna_feature.translate()
                feat.qualifiers.update(dict(translation=translation))
    SeqIO.write(seqs,gbk_path,'genbank')
        

    
def parse_prodigal(infile,start_nr = 1,name_base = 'prodigal_',prefix=''):
    feature_dict = {}
    all_names = {}
    with open(infile) as f:
        text = f.read()
    
    tabs = text.split('\n#')
    groups = [(tabs[i],tabs[i+1]) for i in range(0,len(tabs),2)]
    for header,data in groups:
        features = []
        name_start = header.index('seqhdr=')
        scaffold_name = header[name_start+7:].strip('"')
        if ' ' in scaffold_name:
            scaffold_name = scaffold_name.split(' ')[0]
        data_lines = data.split('\n')
        if '' in data_lines:
            data_lines.remove('')
        for line in data_lines[1:]:
            line_data = line[1:].split('_')
            gene_nr = line_data[0]
            start = int(line_data[1]) - 1 # To adjust to pythonic index
            end = int(line_data[2])
            if line_data[3] == '+':
                strand = 1
            else:   
                strand = -1
            gene_name = '%s%s%i' %(prefix,name_base,start_nr)
            if gene_name not in all_names:
                all_names[gene_name] = 0
            all_names[gene_name] += 1
            start_nr += 1
            quals = {'locus_tag':[gene_name]}
            # Convert to SeqFeature
            fl = FeatureLocation(start,end,strand)
            sf = SeqFeature(fl,type='CDS',strand=strand)
            sf.qualifiers = quals
            features.append(sf)
        feature_dict[scaffold_name] = features
        
    return(feature_dict,start_nr,all_names)
    
    
