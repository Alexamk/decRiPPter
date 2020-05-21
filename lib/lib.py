import os
import sys
import numpy as np
import json

try:
    import cPickle as pickle
except ImportError:
    import pickle

from collections import defaultdict
from Bio import SeqIO

from subprocess import call, Popen, PIPE

from log import return_logger
logger = return_logger(__name__, False)

PYTHON_VERSION = sys.version_info[0]

def autovivify(levels=1, final=dict):
    '''Helper function to easily make layered dictionaries'''
    return (defaultdict(final) if levels < 2 else
            defaultdict(lambda: autovivify(levels - 1, final)))    
 
def clean(folder):
    # Removes a given file, or all the files in a given folder
    # Not a recursive function so that it can't accidentally remove all your files
    if type(folder) == list:
        for file in folder:
            if os.path.isfile(file):
                os.remove(file)
    elif os.path.isfile(folder):
        os.remove(folder)
    elif os.path.isdir(folder):
        files = os.listdir(folder)
        for f in files:
            os.remove(os.path.join(folder,f))
        os.rmdir(folder)

def sortdictkeysbyvalues(dict, index = False, double_dict = False,reverse = False):
    if type(index) == int and double_dict:
        items = [(min([i[index] for i in value.values()]),key) for key,value in dict.items()]
    elif double_dict == True and type(index) == bool and index == False:
        items = [(min([i for i in value.values()]),key) for key,value in dict.items()]
    elif type(index) == int:
        items = [(value[index], key) for key, value in dict.items()]
    else:
        items = [(value, key) for key, value in dict.items()]
    items.sort(reverse=reverse)
    return [key for value, key in items]
    
def sortdictkeysbylenvalues(dict,reverse=False):
    items = [(len(value),key) for key,value in dict.items()]
    items.sort(reverse=reverse)
    return [key for value, key in items]
    
def flip_dict(d):
    d_out = {}
    for key,value in d.items():
        d_out[value] = key
    return d_out

# To more easily load/save pickles

if PYTHON_VERSION == 2:
    def load_pickle(f,t=dict):
        logger.debug('Loading pickle %s' %f)
        try:
            with open(f) as f_open:
                res = pickle.load(f_open)
            return res
        except IOError:
            logger.error('Pickle %s not found' %f)
            raise ValueError('Pickle %s not found' %f)
            
elif PYTHON_VERSION == 3:
    def load_pickle(f,t=dict):
        try:
            with open(f,'rb') as f_open:
                res = pickle.load(f_open,fix_imports)
            return res
        except IOError:
            logger.error('Pickle file %s not found' %(f))
            raise ValueError('Pickle file %s not found' %(f))
            
def store_pickle(data,path,basepath=''):
    with open(os.path.join(basepath, path),'w') as pfile:
        pickle.dump(data,pfile)
        
def find_file(genome_path,ext,extra_req,negative=[]):
    files = os.listdir(genome_path)
    files = [i for i in files if any([i.endswith(j) for j in ext])]
    if len(files) > 1:
        files2 = [i for i in files if any([(j in i) for j in extra_req])]
        if len(files2) > 1:
            files3 = [i for i in files2 if not any([(j in i) for j in negative])]
            if len(files3) > 1:
                logger.error('Multiple files fitting requirements found - genome: %s' %genome)
                raise ValueError('Multiple files fitting requirements found - genome: %s' %genome)
            else:
                f = files3[0]
        else:
            f = files2[0]
    elif files == []:
        return 
    else:
        f = files[0]
    return f
    
def muscle(inf,outf,clw=False,quiet=True):
    commands = ['muscle','-in',inf,'-out',outf]
    if clw:
        commands.append('-clw')
    if quiet:
        commands.append('-quiet')
    logger.debug(' '.join(commands))
    call(commands)
    return outf
        
def replace_pipe(infile):
    outfile = infile.rpartition('.')[0] + '.tab'
    f = open(infile)
    fout = open(outfile,'w')
    for l in f:
        l_new = l.replace('|','\t')
        fout.write(l_new)
    f.close()
    fout.close()
    clean(infile)
    return outfile
        
        
def split_data(l,nr):
    nr_entries = len(l)
    nr_each = int(np.ceil(nr_entries/float(nr)))
    out = []
    for i in range(nr):
        temp = l[i*nr_each:(i+1)*nr_each]
        out.append(temp)
    return(out)
    
def split_dict(d,nr):
    nr_entries = len(d)
    nr_each = int(np.ceil(nr_entries/float(nr)))
    out = [dict() for i in range(nr)]
    i = 0
    for key in d:
        out[i][key] = d[key]
        i += 1
        if i > nr - 1:
            i -= nr
    return out
    
def fuse_dict_list(d1,d2):
    '''Updates the first dictionary with the second one. Works iteratively on layered dictionaries with shared keys betweed d1 and d2. Lists as values are extended'''
    if (type(d1) == dict or dict in d1.__class__.__bases__) and (type(d2) == dict or dict in d2.__class__.__bases__):
        for key in d2:
            if key not in d1:
                d1[key] = d2[key]
            else:
                fuse_dict_list(d1[key],d2[key])
    if type(d1) == list and type(d2) == list:
        d1.extend(d2)
        
def fuse_dict_add(d1,d2):
    for key in d2:
        if key not in d1:
            d1[key] = d2[key]
        else:
            d1[key] += d2[key]
    return d1
        
        
def open_genbank(in_path):
    all_seqs = []
    for seq_record in SeqIO.parse(in_path, "genbank"):
        all_seqs.append(seq_record)
    return(all_seqs)
        
def open_json(path):
    json_file = "" 
    data = open(path,"r")
    for line in data:
        json_file = json_file + line.rstrip("\n")
    json_loaded = json.loads(json_file)
    data.close()
    return(json_loaded)
        
def rev_comp(s):
    d = {"a":"t","c":"g","g":"c","t":"a", "A": "T","C":"G","G":"C","T":"A", \
        "y": "r", "Y":"R", "r": "y", "R":"Y", "w":"w","W":"W","s":"s","S":"S", \
        "k":"m","K":"M","m":"k","M":"K","d":"h","D":"H","v":"b","V":"B", \
        "h":"d","H":"D","b":"v","B":"V","x":"x","X":"X","-":"-","n":"n","N":"N"}
    out = ""
    for c in range(len(s)-1,-1,-1):
        out += d[s[c]]
    return(out)
    
def dict_to_fasta(d,f=False,mode='w'):
    s = ''
    for key in d:
        seq = d[key]
        s += '>%s\n%s\n' %(key,seq)
    if f:
        with open(f,mode) as handle:
            handle.write(s)
    return s
    
def parse_fasta(path):
    infile = open(path)
    out = {}
    name = False
    for line in infile:
        line = line.strip()
        if ">" in line:
            if name and name not in out:
                out[name] = seq
            name_start = line.find(">")
            name = line[name_start+1:]
            seq = ""
        else:
            seq += line
    if name not in out:
        out[name] = seq
    return(out)

def scour_dict(d,k):
    any_true = False
    logger.debug(type(d))
        
    if type(d) == list:
        for element in d:
            res = scour_dict(element,k)
            if res[-1]:
                return(res)
                break
        else:
            return([False])
                
    elif type(d) == dict:
        for key in d:
            res =  scour_dict(d[key],k)
            if res[-1]:
                return([key] + res)
                break
        else:
            return([False])
                
    else:
        res = (k == d)
        return([res])
        
def scour_dict_v(d,k):
    if type(d) == list:
        for element in d:
            res = scour_dict_v(element,k)
            if res[-2]:
                return(res)
                break
        else:
            return([False,None])       
    elif type(d) == dict:
        for key in d:
            res = (key==k)
            if res:
                return([res,d[k]])
                break
            else:
                res = scour_dict_v(d[key],k)
                if res[-2]:
                   return([key] + res)
        else:
            return([False,None])         
    else:
        return([False,None])

# Functions to help with networks/mcl

def write_ssn(allvall,outfile,cutoff=False,i=None):
    f_out = open(outfile,'w')
    for key1 in allvall:
        for key2 in allvall[key1]:
            value = allvall[key1][key2]
            if i != None:
                value = value[i]
            if not(cutoff and value < cutoff):   
                f_out.write('%s\t%s\t%s\n' %(key1,key2,value))
    f_out.close() 

def read_mcl(mcl_file,min_group_size=0):
    all_groups = []
    with open(mcl_file,'r') as f:
        for l in f:
            group = l.strip().split('\t')
            if len(group) >= min_group_size:
                all_groups.append(group)
    return(all_groups)
    
# Functions to call commands

def run_mcl(inf,outf,threads,I=1.5):
    logger.info("MCL START")
    I = str(I)
    cmd = ["mcl",inf,"-o",outf,"--abc","-I", I, '-te', str(threads)]
    stdout, stderr = run_cmd(cmd)
    return(stdout)

def run_prodigal_cmd(infile,outfile):
    cmds = ['prodigal','-i',infile,'-o',outfile,'-f','sco']
    stdout, stderr = run_cmd(cmds)
    
def run_cmd(commands):
    logger.debug('Running system command %s' %(' '.join(commands)))
    p = Popen(commands, stdout=PIPE, stderr=PIPE)
    p.wait()
    stdout,stderr = p.stdout.read(),p.stderr.read()
    if stderr:
        logger.debug(stderr)
    else:
        logger.debug(stdout)
    return stdout, stderr


# Functions to process the found groups/clusters
def pairs_to_groups(score_dict,min_group_size):
    '''Forms groups from a dictionary containing scores'''
    groups = []
    for key in score_dict:
        matches = score_dict[key]
        for key2 in matches:
            for group in groups:
                if key in group or key2 in group:
                    group.add(key)
                    group.add(key2)
                    break
            else:
                newgroup = set((key,key2))
                groups.append(newgroup)
    # Now fuse overlapping groups
    added = 1
    while added != 0:
        skip = []
        added = 0
        new_groups = []
        for group in groups:
#            print('At group: %s' %group)
            if group in skip:
#                print('Skipping gr1')
                continue
            new_group = set(group)
            for group2 in groups:
#                print('Comparing to group %s' %group2)
                if group == group2:
#                    print('Skipping gr2')
                    continue
                for gene in new_group:
                    if gene in group2:
                        # Overlap --> add the groups
                        new_group = new_group | group2
#                        print('Overlap found. New group: %s' %new_group)
                        # Don't use this group as baseline for comparison again
                        skip.append(group2)
                        added += 1
#                        print('Total added: %s' %added)
                        break
#            print('Adding new group to newgroups: %s' %(new_group))
            new_groups.append(new_group)
        if added > 0:
            groups = new_groups
#            print('Replacing groups with new_groups')
    groups_filtered = [list(group) for group in groups if len(group) >= min_group_size]
    return(groups_filtered)

def assign_groups(operons,smorf_groups):
    # Assign a number to each smorf group
    # Add a variable precursor_group to each smorf
    i = 0
    group_names = []
    for group in smorfs_groups:
        i += 1
        group_name = 'precursor_group_%i' %i
        group_names.append(group_name)
        for gene in group:
            gene = operons[gene]
            if gene == None:
                # Gene from the extra dict, e.g. MIBiG
                continue
            if not hasattr(gene,'precursor_group'):
                setattr(gene,'precursor_group',group_name)
            else:
                raise ValueError('%s is part of multiple groups! Fix pairs to group function' %gene)
            
def assign_groups_operons(operons):
    for operon in operons:
        precursor_groups = []
        for precursor in operon.precursors:
            if hasattr(precursor,'precursor_group'):
                precursor_groups.append(gene.precursor_group)   
     
def assign_group_collection(operons):
    precursor_groups = {}
    for name in group_names:
        ops = [op for op in operons.itersubset(precursor_group=name)]
        precursor_groups[name] = ops
        
# Function to read HTML template

def read_template(f,p=False):
    html_pieces = []
    with open(f) as handle:
        current_part = ''
        waiting = False
        for line in handle:
            if line == '<!--BREAK-->\n':
                html_pieces.append(current_part)
                current_part = ''
                waiting = True
            elif waiting and line == '<!--CONTINUE-->\n':
                waiting = False
            elif not waiting:
                current_part += line
            else:
                continue
        html_pieces.append(current_part)
    return html_pieces
    
    
# Generic container object for various purposes
class Container():
    # Container for provided settings and to store info on genes
    def __init__(self):
        pass
    def __repr__(self):
        if hasattr(self,'name'):
            return('Container object (%s)' %self.name)
        else:
            return('Container object')
    def setattrs(self,**kwargs):
        for key,value in kwargs.items():
            setattr(self,key,value)
            
class Edge():
    def __init__(self,weight,data):
        self.weight = weight
        self.data = data
        
    def __repr__(self):
        return('Edge object (weight: %.2f; data: %s' %(self.weight,self.data))
# Classes to use in itergenes and itersubset functions
class Max():
    def __init__(self,a):
        self.a = a
    def __eq__(self,k):
        return((type(k) in [int,float] or any([i in [int,float] for i in k.__class__.__bases__]) ) and k<=self.a)
    def __ne__(self,k):
        return( not(type(k) in [int,float] or any([i in [int,float] for i in k.__class__.__bases__])) or k>self.a)
    def __repr__(self):
        return('Max(%s)' %(self.a.__repr__()))

class Min():
    def __init__(self,a):
        self.a = a
    def __eq__(self,k):
        return( (type(k) in [int,float] or any([i in [int,float] for i in k.__class__.__bases__]) ) and k>=self.a)
    def __ne__(self,k):
        return( not(type(k) in [int,float] or any([i in [int,float] for i in k.__class__.__bases__])) or k<self.a)
    def __repr__(self):
        return('Min(%s)' %self.a.__repr__())

class Contains():
    def __init__(self,*args):
        self.a = args
    def __eq__(self,k):
        for item in self.a:
            if item not in k:
                return(False)
        else:
            return(True)
    def __ne__(self,k):
        for item in self.a:
            if item not in k:
                return(True)
        else:
            return(False)
    
    def __repr__(self):
        return('Contains(%s)' %','.join(self.args))

class IsIn():
    def __init__(self,*args):
        self.a = args
        
    def __eq__(self,k):
        return(k in self.args)
        
    def __ne__(self,k):
        return(k not in self.args)
        
    def __repr__(self):
        return('IsIn(%s)' %(self.args.__repr()))

class DictValue():
    def __init__(self,key,value):
        self.key = key
        self.value = value
        
    def __eq__(self,k):
        if self.key != None:
            return((type(k) == dict or dict in k.__class__.__bases__) and self.key in k and k[self.key] == self.value)
        else:
            return((type(k) == dict or dict in k.__class__.__bases__) and self.value in k.values())
        
    def __ne__(self,k):
        if self.key != None:
            return(not(type(k) == dict or dict in k.__class__.__bases__) or self.key not in k or k[self.key] != self.value)
        else:
            return(not(type(k) == dict or dict in k.__class__.__bases__) or self.value not in k.values())
        
    def __repr__(self):
        return('DictValue(key=%s,value=%s)' %(self.key.__repr__(),self.value.__repr__()))
        


class Length():
    def __init__(self,a,req_type=None):
        self.a = a
        self.req_type = req_type
    def __eq__(self,k):
        if self.req_type:
            return(len(k) == self.a and type(k) == self.req_type)
        else:
            return(len(k) == self.a)
    def __ne__(self,k):
        if self.req_type:
            return(len(k) != self.a or type(k) != self.req_type)
        else:
            return(len(k) != self.a)
    
    def __repr__(self):
        return('Length(%s)' %self.a.__repr__())


class OR():
    def __init__(self,*conditions):
        self.conditions = conditions
    
    def __eq__(self,k):
        return (any(c==k for c in self.conditions))
    
    def __ne__(self,k):
        return (all(c!=k for c in self.conditions))
        
    def __repr__(self):
        return('OR(%s)' %(', '.join([i.__repr__() for i in self.conditions])))
        
class AND():
    def __init__(self,*conditions):
        self.conditions = conditions
    
    def __eq__(self,k):
        return (all(c==k for c in self.conditions))
    
    def __ne__(self,k):
        return (any(c!=k for c in self.conditions))
        
    def __repr__(self):
        return('AND(%s)' %(', '.join([i.__repr__() for i in self.conditions])))
        
class NOT():
    def __init__(self,a):
        self.a = a
        
    def __eq__(self,k):
        return(k != self.a)
        
    def __ne__(self,k):
        return(k == self.a)
    
    def __repr__(self):
        return('NOT(%s)' %self.a.__repr__())
        

        
