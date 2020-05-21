# Convert the .gbk files from MiBIG into Operon object for comparison with other found operons
# Method 1) Get all the genes in the genbank file
# Method 2) Try to restrict to the genes in operon with predicted precursors

import os
import shutil
import pickle
import sys
import argparse
import pandas as pd

PYTHON_VERSION = sys.version_info[0]
if PYTHON_VERSION == 2:
    import ConfigParser as configparser
elif PYTHON_VERSION == 3:
    import configparser

import Operon_formation as op
import lib.hmmer_operons as hmm

from genome_prep import parse_genbank_prep,parse_gbk_to_dict
from lib.smorfs import smorfs_analysis
from lib.lib import open_genbank,open_json,scour_dict_v,load_pickle,store_pickle,fuse_dict_list
from lib.Genes import OperonCollection, Operon

global debug
debug = True

def get_ripp_bgcs(path_gbk,path_json,path_out):
    errors = []
    json_files = os.listdir(path_json)
    gbk_files = os.listdir(path_gbk)
    ripp_gbks = []
    out_folders = os.listdir(path_out)
    for filename in json_files:
        if filename[0] == "." or filename[-5:] != ".json":
            continue
        try:
            json_p = open_json(path_json+filename)
        except:
            errors.append("json opening error "+filename)
            continue
        if filename[:-5]+".gbk" not in gbk_files:
            errors.append("no gbk " + filename)
            continue
        cl = scour_dict_v(json_p,u'biosyn_class')[-1]
        if 'RiPP' in cl:
            name_clean = filename[:-5]
            gbk_file = name_clean + '.gbk'
            ripp_gbks.append(name_clean)
            if name_clean not in out_folders:
                os.mkdir(path_out + name_clean)
            if not os.path.isfile(path_out + name_clean + os.sep + gbk_file):
                shutil.copy(path_gbk + gbk_file, path_out + name_clean + os.sep)
    return(ripp_gbks,errors)

def set_mibig_domains(mibig_dict):
    for bgc in mibig_dict:
        all_domains = []
        for gene in bgc.itergenes():
            if hasattr(gene,'domains'):
                all_domains.extend(gene.domains)
        bgc.all_domains = all_domains

def add_mibig_info(mibig_dict,f):
    data = pd.read_excel(f,sheet_name='MIBiG')
    for index,row in data.iterrows():
        name = row['MIBiG accession']
        mibig_bgc = mibig_dict[name]
        if not hasattr(mibig_bgc,'precursors'):
            mibig_bgc.precursors = {}
        mibig_bgc.organism = row['Organism']
        mibig_bgc.subtype = row['RiPP subtype']
        mibig_bgc.precursors[row['Precursor']] = row['RiPP sequence']
        mibig_bgc.product = row['Main product']
        mibig_bgc.mibig = True
    # Now set all precursors as an attribute of the mibig_dict
    all_prec = {}
    for mibig_bgc in mibig_dict:
        if hasattr(mibig_bgc,'precursors'):
            for prec,seq in mibig_bgc.precursors.items():
                new_prec = '%s_%s' %(mibig_bgc.name,prec)
                all_prec[new_prec] = seq
        else:
            print('No precursors found for %s' %mibig_bgc.name)
    mibig_dict.precursors = all_prec
    
def convert_to_operons(genome_dict):
    op_dict = {}
    for genome in genome_dict:
        gene_dict = {}
        for gene in genome.itergenes():
            gene_dict[gene.name] = gene
        operon = Operon(gene_dict,name=genome.name)
        op_dict[genome.name] = operon
    collection = OperonCollection(op_dict,name='MIBiG')
    return collection
    
    
#def read_bagel_dict(f):
#    data_c1 = pd.read_excel()
    
    

def main(jsonpath,gbkpath,path_out,settings):
    bgcs,errors = get_ripp_bgcs(gbkpath,jsonpath,path_out)
    print('Parsing genbank (%i files)' %len(bgcs))
    genome_seqs,names_double,scaf_names_double,file_parsed = parse_genbank_prep(bgcs,path_out)
    mibig_dict,smorfs,name2realname,newly_named_genes,namesreplaced,newly_named_scaffolds,scafname2realname =\
                                 parse_gbk_to_dict(path_out,genome_seqs,names_double,scaf_names_double,include_pseudo=True)
    return mibig_dict
    # Method 1: All genes in the genbank file
    # Method 2: First form operons around the precursors
    method = 1
    if method == 1:
        print('Running hmmer')    
        hmm.main(mibig_dict,settings)
    elif method == 2:
        pickle_path = settings['pickle_path']
        SVM3 = load_pickle(pickle_path + 'SVM_poly_3.pkl')
        SVM4 = load_pickle(pickle_path + 'SVM_poly_4.pkl')
        SVMr = load_pickle(pickle_path + 'SVM_rbf.pkl')
        mibig_dict = smorfs_analysis(smorfs,mibig_dict,path_out,SVM3,SVM4,SVMr,settings['cores'],'smORF')
        # Discontinued at this point
        # 25% of the operons do not have a predicted precursor
        # --> Use method 1 for now
        # Possible fix later: Assign the locations of each precursor to each mibig operon
#        operons = op.find_all_operons(mibig_dict,distance = settings['simple_dist'])
#        operons = op.fuse_operons(operons,50,'operon')
    set_mibig_domains(mibig_dict)
    return mibig_dict
    
def parse_arguments():
    config = configparser.ConfigParser()
    configfile = 'config.ini'
    if not os.path.isfile(configfile):
        configfile = '../config.ini'
    if not os.path.isfile(configfile):
        print('Warning: config file not found in cwd')
        settings = {}
    else:
        config.read(configfile)
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
                
    parser = argparse.ArgumentParser()
    parser.add_argument('-o',help='path to the mibig_databases',dest='mibig_folder')
    parser.add_argument('-v',help='MiBIG version',dest='version')
    parser.add_argument('--pfam',help='path to the PFAM databases',dest='pfam_db_path')
    parser.add_argument('--tigrfam',help='path to the TIGRFAM databases',dest='tigrfam_db_path')
    parser.add_argument('-c','--cores', help='number of cores to be used on this run',type=int)
    parser.add_argument('-bpp',help='Base data path, where the results will be stored',dest='base_data_path')
    
    args = parser.parse_args()
    for key,value in args.__dict__.items():
        # Update only if the args value is not None
        if value != None:
            settings[key] = value
    return settings
    
    
if __name__ == '__main__':
    settings = parse_arguments()
    mibig_path = settings['mibig_folder']
    version = settings['version']
    jsonpath = os.path.join(mibig_path, 'mibig_json_%s/' %version)
    gbkpath = os.path.join(mibig_path, 'mibig_gbk_%s/' %version)
    RiPP_path = os.path.join(mibig_path, 'RiPP/')
    hmmpath = os.path.join(RiPP_path,'hmm/')
    hmmpath_pfam = os.path.join(hmmpath,'pfam/')
    hmmpath_tigrfam = os.path.join(hmmpath, 'tigrfam/')
    path_out = os.path.join(RiPP_path + 'Clusters/')
    pfam_db_path = settings['pfam_db_path']
    tigrfam_db_path = settings['tigrfam_db_path']
    data_path = settings['base_data_path']
    settings['hmm_paths'] = ((hmmpath_pfam,pfam_db_path),(hmmpath_tigrfam,tigrfam_db_path))
    settings['paths'] = {}
    settings['paths']['hmm_basic'] = hmmpath
    for folder in [RiPP_path,hmmpath,hmmpath_pfam,hmmpath_tigrfam,path_out]:
        if not os.path.isdir(folder):
            os.mkdir(folder)
    
    mibig_dict = main(jsonpath,gbkpath,path_out,settings)
    exit()
    excel_file = '/vol/home/kloostermanam/Documents/Scripts/ripp_v2/data/RiPP_dataset.xlsx'
    add_mibig_info(mibig_dict,excel_file)
    # Remove the geneclusters that are misannotated as RiPPs
    tbr = ['BGC0001559','BGC0001734']
    for t in tbr:
        _ = mibig_dict.pop(t)
    storage_path = os.path.join(data_path,'Pickles','mibig_domains_%s.pkl' %version)
    store_pickle(mibig_dict,storage_path)
