# License: GNU Affero General Public License v3 or later

# Wrapper script to run antiSMASH in it's own environment, after step 1
# antiSMASH - requires a different python environment;
#             this script can be called standalone to run antismash

import os
import argparse
import pickle

from subprocess import run, PIPE

from lib.log import return_logger, setup_loggers

logger = return_logger('antiSMASH_wrapper', False)


def run_cmd(commands):
    logger.debug('Running system command "%s"' %(' '.join(commands)))
    if logger.getEffectiveLevel() <= 10:
        completed_process = run(commands, stdout=PIPE, stderr=PIPE)
    else:
        cp = run(commands)
#    if completed_process.stderr:
#        logger.debug(stderr)
#    if completed_process.stdout:
#        logger.debug(stdout)
#    return stdout, stderr
    
def load_pickle(f,t=dict):
    logger.debug('Loading pickle %s' %f)
    with open(f,'rb') as f_open:
        res = pickle.load(f_open)
    return res

def run_antismash_minimal(infile, outfolder, cores):
    # Since only the gene cluster boundaries are neccessary, antismash minimal 
    # is used by default.
    commands = ['antismash','--minimal','-c',str(cores),'--output-dir',outfolder,infile]
    run_cmd(commands)

def run_main(args):
    path = os.path.join(args.outputfolder, args.name)
    pickle_path = os.path.join(path, 'data', 'Pickles')
    pickle_file = os.path.join(pickle_path,'gbk_files_for_antismash.pkl' )
    genome_path = os.path.join(path, 'Genomes')
    gbk_files_parsed = load_pickle(pickle_file)
    run_all_antismash(gbk_files_parsed, args.cores)
        
def run_all_antismash(gbk_files_parsed, cores):
    logger.info('Running antiSMASH minimal for %i genomes on %i cores' %(len(gbk_files_parsed), cores))
    for genome, gbk_file in gbk_files_parsed.items():
        genome_path = os.path.dirname(gbk_file)
        outfolder = os.path.join(genome_path, 'antismash_minimal')
        if os.path.isdir(outfolder):
            logger.warning('antiSMASH results already found for genome %s. Continuing' %genome)
            continue
        logger.info('Running antiSMASH minimal for genome %s' %genome)
        run_antismash_minimal(gbk_file, outfolder, cores)

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('name',metavar='PROJECT NAME',type=str,help='The name of the project')
    parser.add_argument('-o','--outputfolder',help='path where the results folder is located',default=os.getcwd())
    parser.add_argument('-c','--cores',help='number of cores to use',default=1,type=int)
    parser.add_argument('-v', '--verbosity', help='Verbosity (0-3), higher means more output', choices = [0, 1, 2, 3], type=int, default=2)
    
    args = parser.parse_args()
    settings = {}
    for key,value in args.__dict__.items():
        # Update only if the args value is not None
        if value != None:
            settings[key] = value
    
    return args, settings

if __name__ == '__main__':
    # As a stand-alone function, can be called as a wrapper for antismash
    # Will run antismash minimal on all the genomes given an project folder
    # If a genome_dict pickled file is found, the genbank files parsed to create 
    # that object will be used for antiSMASH. 
    #TODO Allow more optional antismash arguments
    args, settings = parse_arguments()
    setup_loggers(settings)
    genome_dict = run_main(args)

    
    
    
