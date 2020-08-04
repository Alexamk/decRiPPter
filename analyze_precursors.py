# License: GNU Affero General Public License v3 or later

# Add-on for the decRiPPter pipeline (https://github.com/Alexamk/decRiPPter)
# Standalone script to just analyze precursor sequences
# Can read in a fasta file or raw sequences from stdin

import argparse
import os

try:
    import cPickle as pickle
except ImportError:
    import pickle


from lib.lib import parse_fasta
from lib.smorfs import run_SVM

def load_pickle(f,t=dict):
    try:
        with open(f) as f_open:
            res = pickle.load(f_open)
        return res
    except IOError:
        raise ValueError('Pickle %s not found' %f)

def get_sequences(settings):
    if settings['input_type'] == 'seq':
        sequences = [('STDIN_precursor', settings['input'])]
    elif settings['input_type'] == 'fasta':
        seq_dict = parse_fasta(settings['input'])
        sequences = [(header, seq) for header, seq in seq_dict.items()]
    return sequences
    
def write_output(settings, sequences, outfile):
    with open(outfile, 'w') as handle:
        for sequence in sequences:
            if settings['output_type'] == 'simple':
                to_write = [sequence[0], sequence[-1]]
            else:
                to_write = sequence
            handle.write('%s\n' %('\t'.join(to_write)))
                
def filter_sequences(settings, sequences):
    filtered = []
    for sequence in sequences:
        score = float(sequence[-1])
        if score >= settings['cutoff']:
            filtered.append(sequence)
    return filtered

def load_pickles(settings, SVM_files):
    SVM_models = []
    for SVM_file in SVM_files:
        model = load_pickle(os.path.join(settings['SVM_model_path'], SVM_file))
        SVM_models.append(model)
    return SVM_models

def check_args(settings, SVM_files):
    if not os.path.isdir(settings['outputfolder']):
        raise ValueError('Outputfolder not found. Please specify an existing folder with the -o flag')
    if not os.path.isdir(settings['SVM_model_path']):
        raise ValueError('SVM model path not found. Please point to it with the --SVM_model_path flag') 
    pickle_files = os.listdir(settings['SVM_model_path'])
    for r in SVM_files:
        if r not in pickle_files:
            raise ValueError('File %s not found in SVM model path.' %r)

def parse_arguments():
    
    settings = {}
    parser = argparse.ArgumentParser()
    
    parser.add_argument('-i', '--input', help='Input: either a fasta file or a single peptide sequence to analyze (specify which with -t)', required=True)
    parser.add_argument('-t', '--input_type', help='Type of input given with the -i flag. Either fasta or seq (default)', choices=['fasta', 'seq'], default='seq')
    parser.add_argument('-c', '--cutoff', help='The cutoff of the SVM score to use (0-1).', default=0.9, type=float)
    parser.add_argument('-o','--outputfolder',help='Folder in which the results will be written to',default='.')
    
    parser.add_argument('--output_type', choices=['simple', 'detailed'], default='simple',
                        help='Choose between a simple output (only headers with final score) or a detailed one (all protein features and individual scores of the three SVMs')
#    parser.add_argument('--keep_negatives', help='Also show the precursors in the output that did not make the cutoff', default=False, action='store_true')
    parser.add_argument('--SVM_model_path', help='Point to the folder containing the three SVM models (downloaded with decRiPPter)',
                        default='data/Pickles/')
    
    args = parser.parse_args()
    for key,value in args.__dict__.items():
        # Update only if the args value is not None
        if value != None:
            settings[key] = value
        
    return args,settings



# def main():
if __name__ == '__main__':
    args, settings = parse_arguments()
    SVM_files = ['SVM_poly_3.pkl', 'SVM_poly_4.pkl', 'SVM_rbf.pkl']
    if settings['input_type'] == 'seq':
        outfile = os.path.join(settings['outputfolder'], 'STDIN_decRiPPter_precursor_prediction.txt')
    else:
        basename = os.path.basename(settings['input']).rpartition('.')[2]
        outfile = os.path.join(settings['outputfolder'],'%s_decRiPPter_precursor_prediction.txt' %basename)
    check_args(settings, SVM_files)
    print('Loading models')
    SVM3, SVM4, SVMr = load_pickles(settings, SVM_files)
    print('Analyzing sequences')
    sequences = get_sequences(settings)
    sequences_analyzed, errors, total = run_SVM(sequences, SVM3, SVM4, SVMr)
    print('Analyzed %i sequences. %i gave errors.' %(total, errors))
    filtered_sequences = filter_sequences(settings, sequences_analyzed)
    print('%i sequences passed the filter.' %(len(filtered_sequences)))
    write_output(settings, filtered_sequences, outfile)
    
    
    
