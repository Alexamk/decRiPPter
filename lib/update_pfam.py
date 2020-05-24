# License: GNU Affero General Public License v3 or later

from lib import load_pickle,store_pickle
import os

def read_pfam(p):
    d = {}
    with open(p) as f:
        for l in f:
            tabs = l.strip().split('\t')
            pfam_acc,clan_acc,clan_ID,pfam_ID,pfam_descr = tabs
            d[pfam_acc] = pfam_descr
    
    return(d)
    
if __name__ == '__main__':
    import sys
    args = sys.argv
    if len(args) < 2:
        print('USAGE: python update_pfam.py /path/to/Pfam-A.clans.tsv')
        exit()
    else:
        path = args[1]
        d = read_pfam(path)
        data_path = '../data/Pickles/'
        descr_path = os.path.join(data_path,'domain_descr.pkl')
        old_d = load_pickle(descr_path)
        old_d.update(d)
        store_pickle(old_d,descr_path)
        
