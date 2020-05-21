# Organization

# Each genome, scaffold and gene has a unique name/identifier
# One is assigned to it if double names exist

# A genecollection is a collection of genes. Genes are contained in a dictionary
# Operons/BGCs are genecollections
# Scaffolds are genecollections

# Each genome is an object, containing its scaffolds in a dictionary
# Each gene is an object

from collections import OrderedDict
import numpy as np

from lib import Container, Edge
from log import return_logger

logger = return_logger(__name__, False)

def sortx_start(x):
    return(x.start)
def sortx_name(x):
    return(x.name)

class D(dict):
    # Similar to a normal dictionary except
    # a) It can contain attrs
    # b) It iterates over values by default
    # c) If sortkey is defined, it iterates over values sorted by sortkey
    def __init__(self,*args,**kwargs):
        if len(args) > 0 and type(args[0]) == dict:
            args = [(k,v) for k,v in args[0].items()]
        dict.__init__(self,args)
        self.__dict__.update(kwargs)

    def __getitem__(self,k):
        if hasattr(self,'quick_dict'):
            try:
                return(self.quick_dict[k])
            except KeyError:
                return None
        else:
            try:
                return(dict.__getitem__(self,k))
            except:
                return
        
    def __setitem__(self,k,v):
        dict.__setitem__(self,k,v)
        if hasattr(self,'quick_dict') and hasattr(v,'name'):
            self.quick_dict[k] = v
            
    def __eq__(self,other):
        if isinstance(other, self.__class__):
            if dict.__eq__(self,other) and self.name == other.name:
                return True
        return False
        
    def __ne__(self,other):
        if dict.__ne__(self,other) or self.name != other.name or not isinstance(other, self.__class__):
            return True
        return False

    def __iter__(self):
        for item in self.values():
            yield item
            
    def __repr__(self):
        return('%s(' %self.__class__.__name__ + ','.join(['%s=%s' %(key,value) for key,value in sorted(self.__dict__.items()) if key not in ['sortkey','sortrev','subcollection_keyword'] or type(value) in [list,dict,set,tuple]]) + ')')
    
    def iter_sorted(self):
        if hasattr(self,'sortkey') and hasattr(self,'sortrev'):
            for v in sorted(self.values(),key=self.sortkey,reverse=self.sortrev):
                yield v
        else:
            for v in self.values():
                yield v
        
    def setattrs(self,**kwargs):
        self.__dict__.update(kwargs)
        
    def iter_all(self):
        yield self
        for item in self:
            if hasattr(item,'iter_all'):
                for item2 in item.iter_all():
                    yield item2
            else:
                yield item
        
    def write_fasta(self,outtype='prot',filename=False,extra_dict=False,max_len=75,**kwargs):
        # All the kwargs are passed on to the itergenes function to select the genes
        # All the dna sequences of the genes are written to a fasta file
        
        def return_entry(header,seq,max_len):
            out = ''
            i = 0
            while i * 75 < len(seq):
                out += seq[i*max_len:(i+1)*max_len] + '\n'
                i += 1
            return '>%s\n%s' %(header,out)
        
        written = set()
        out_dict = OrderedDict()
        out = ''
        
        if outtype == 'dna':
            for item in self.itergenecoll(**kwargs):
                out_dict[item.name] = item.seq
                if item.name not in written:
                    out += return_entry(item.name,item.seq,max_len)
                    written.add(item.name)
        elif outtype == 'prot':
            for gene in self.itergenes(**kwargs):
                out_dict[gene.name] = gene.transl
                if gene.name not in written:
                    out += return_entry(gene.name,gene.transl,max_len)
                    written.add(gene.name)
        if extra_dict:
            out_dict.update(extra_dict)
            for key in extra_dict:
                out += '>%s\n%s\n' %(key,extra_dict[key])
        if filename:
            if type(filename) == str:
                with open(filename,'w') as f:
                    f.write(out)
            elif type(filename) == list:
                # First split the dicts
                out_dicts = [OrderedDict() for _ in range(len(filename))]
                n = 0
                for key in out_dict:
                    out_dicts[n][key] = out_dict[key]
                    n += 1
                    n = n % (len(filename))
                # Now write the files
                for n in range(len(filename)):
                    with open(filename[n],'w') as handle:
                        for key in out_dicts[n]:
                            text = return_entry(key,out_dicts[n][key],max_len)
                            handle.write(text)
        return out,out_dict
        
        
class GeneCollection(D):
    def __init__(self,*args,**kwargs):
        if 'set_genes' in kwargs and kwargs['set_genes'] == True:
            set_genes = kwargs.pop('set_genes')
        else:
            set_genes = False
        D.__init__(self,*args,**kwargs)
        self.sortkey=sortx_start
        self.sortrev=False
        self.subcollection_keyword = 'genes'
        if set_genes:
            self.set_genes()
        self.flanks_set=False
        
    def __getitem__(self,k):
        if hasattr(self,'quick_dict'):
            try:
                return(self.quick_dict[k])
            except KeyError:
                return None
        else:
            try:
                return(dict.__getitem__(self,k))
            except:
                return 
        
    def itergenes(self,core=False,sorted=False,**kwargs):
        if core and hasattr(self,'core_genes') and sorted:
            itr = sorted(self.core_genes,key=self.sortkey)
        elif core and hasattr(self,'core_genes'):
            itr = self.core_genes
        elif sorted:
            itr = self.iter_sorted()
        else:
            itr = self
        if self.subcollection_keyword in kwargs:
            subset = kwargs.pop(self.subcollection_keyword)
        else:
            subset = False
        for gene in itr:
            if subset and sub_item.name not in subset:
                continue
            for key,value in kwargs.items():
                if not (hasattr(gene,key) and getattr(gene,key) == value):
                    break
            else:
                yield gene

    def itergenecoll(self):
        yield self
                
    def set_genes(self):
        self.genes = sorted(self.values(),key=self.sortkey,reverse=self.sortrev)
            
    def setflanks(self):
        if not self.flanks_set:
            prev_gene = None
            intergenics = []
            for gene in self.itergenes(sorted=True):
                if prev_gene == None:
                    dist_to_prev = None
                else:
                    dist_to_prev = gene.start - prev_gene.end
                gene.update(left_flank=prev_gene,left_dist=dist_to_prev)
                if gene.intergenic:
                    intergenics.append(gene)
                else:
                    for intergene in intergenics:
                        right_gene = gene
                        right_dist = gene.start - intergene.end
                        intergene.update(right_flank=right_gene,right_dist=right_dist)
                    if prev_gene != None:
                        prev_gene.update(right_flank=gene,right_dist=dist_to_prev)
                    prev_gene=gene
                    intergenics=[]
            if prev_gene:
                prev_gene.update(right_flank=None,right_dist= None)
            for intergene in intergenics:
                intergene.update(right_flank=None,right_dist= None)
            self.flanks_set=True
        
    def switch_inactive(self):
        if self.active:
            for item in ['genome','scaffold']:
                if hasattr(self,item):
                    name = getattr(self,item).name
                    setattr(self,item,name)
            self.active=False
        
    def switch_active(self,collection):
        if not self.active:
            for item in ['genome','scaffold']:
                if hasattr(self,item):
                    name = getattr(self,item)
                    obj = collection[name]
                    setattr(self,item,obj)
            self.active=True

class Operon(GeneCollection):
    def __init__(self,*args,**kwargs):
        GeneCollection.__init__(self,*args,**kwargs)
        genes = [gene for gene in self.iter_sorted()]
        self.core_genes = genes
        gene = genes[0]
        for attr in ['scaffold','genome']:
            if not hasattr(self,attr):
                setattr(self,attr,getattr(gene,attr))
        gene_end = genes[-1]
        self.start = gene.start
        self.end = gene_end.end
#        if hasattr(self.scaffold,'seq'):
#            self.seq = self.scaffold.seq[self.start:self.end]
        self.strand = gene.strand
        self.calc_cog()
        
    def __repr__(self):
        if getattr(self,'COG_avg'):
            return('Operon(name=%s,COG=%.3f+-%.3f)' %(self.name,self.COG_avg,self.COG_std))
        else:
            return('Operon(name=%s,COG=None)' %self.name)

    def calc_cog(self):
        all_cog = []
        for gene in self.core_genes:
            if hasattr(gene,'COG'):
                all_cog.append(gene.COG)
        self.COG = all_cog
        if all_cog == []:
            self.COG_avg = None
            self.COG_std = None
        elif len(all_cog) == 1:
            self.COG_avg = all_cog[0]
            if hasattr(self.values()[0].genome,'COG_global_std'):
                self.COG_std = self.values()[0].genome.COG_global_std
            else:
                self.COG_std = 0.05
        else:
            self.COG_avg = np.average(all_cog)
            self.COG_std = np.std(all_cog)
            
    def fuse(*args,**kwargs):
        '''Returns a new Operon fused from Operons in the args'''
        new_gene_dict = {}
        for arg in args:
            for gene in arg:
                new_gene_dict[gene.name] = gene
        if 'name' in kwargs:
            name = kwargs['name']
        else:
            name = args[0].name
        new_operon = Operon(new_gene_dict,name=name,active=args[0].active)
        return new_operon
        
    def get_seq(self,extra=False):
        if extra:
            return(self.scaffold.seq[self.extra_start:self.extra_end])
        else:
            return(self.scaffold.seq[self.start:self:end])


class CollectionCollection(D):
    
    def __init__(self,*args,**kwargs):
        D.__init__(self,*args,**kwargs)
        self.subcollection_keyword = False
        
    def __getitem__(self,k):
        if hasattr(self,'quick_dict'):
            try:
                return(self.quick_dict[k])
            except KeyError:
                return None
        else:
            try:
                return(dict.__getitem__(self,k))
            except:
                for subitem in self:
                    res = subitem[k]
                    if res != None:
                        return res

    def set_quick(self):
        self.quick_dict = {}
        for item in self.iter_all():
            if hasattr(item,'name'):
                self.quick_dict[item.name] = item


    def itergenes(self,**kwargs):
#        print('itergenes on type: %s' %type(self))
        if self.subcollection_keyword in kwargs:
            subset = kwargs.pop(self.subcollection_keyword)
        else:
            subset = False
#        print('Subset = %s' %subset)
#        print('remaining kwargs = %s' %kwargs)
        for sub_item in self:
            if subset and sub_item.name not in subset:
                continue
            for gene in sub_item.itergenes(**kwargs):
                yield gene
                
    def itergenecoll(self,**kwargs):
        # Return all GeneCollection subitems
        if self.subcollection_keyword in kwargs:
            subset = kwargs.pop(self.subcollection_keyword)
        else:
            subset = False
        for sub_item in self:
            if subset and sub_item.name not in subset:
                continue
            for item in sub_item.itergenecoll():
                if GeneCollection in type(item).__bases__:
                    for key,value in kwargs.items():
                        if not (hasattr(item,key) and getattr(item,key) == value):
                            break
                    else:
                        yield item
                else:
                    yield item
                    
    def itersubset(self,sorted=False,**kwargs):
        if sorted:
            itr = self.iter_sorted()
        else:
            itr = self
        for item in itr:
            for key,value in kwargs.items():
                if hasattr(item,key):
                    item_value = getattr(item,key)
                    if value != getattr(item,key) and not(type(item_value) == list and value in item_value):
                        break
                else:
                    break
            else:
                yield item
                    
    def del_all(self,attr):
        for item in self:
            if hasattr(self,attr):
                delattr(self,attr)
    
    # The following functions 'activate'/'inactivate' the collection
    # In it's inactive form, it can be pickled (references are names, no infinite internal loops)
    # In it's active form, it can be used more easily (references are objects), but large enough dictionaries can't be pickled (recursion depth limit)     
    def switch_inactive(self):
        for gene in self.itergenes():
            gene.switch_inactive()
        for item in self:
            item.switch_inactive()
        self.active=False
    
    def switch_active(self,collection=False):
        if hasattr(self,'name'):
            logger.debug('%s switching to active' %self.name)
        if not collection:
            collection = self
        for gene in self.itergenes():
            gene.switch_active(collection)
        for item in self:
            item.switch_active(collection)
        self.active=True

    def find_by_precursor(self,prec_query):
        if hasattr(self,'precursor_dict'):
            try:
                return self.precursor_dict[prec_query.name]
            except KeyError:
                return 
        else:
            for item in self:
                if hasattr(item,'precursors') and prec_query in item.precursors:
                    return item

class GenomeCollection(CollectionCollection):
    def __init__(self,*args,**kwargs):
        D.__init__(self,*args,**kwargs)
        self.sortkey=sortx_name
        self.sortrev=False
        self.subcollection_keyword = 'genomes'
    
    def remove(self,genomes):
        for genome in genomes:
            self.pop(genome)
    
    def prep(self):
        self.set_quick()
        start_active = self.active
        if not start_active:
            self.switch_active(self)
        for genome in self:
            logger.debug('Setting flanks for genome %s' %genome.name)
            for scaffold in genome:
                scaffold.setflanks()
        if not start_active:
            self.switch_inactive()
            
class OperonCollection(CollectionCollection):
    def __init__(self,*args,**kwargs):
        if 'prep' in kwargs:
            prep = kwargs.pop('prep')
        else:
            prep = True
        D.__init__(self,*args,**kwargs)
        self.sortkey=sortx_start
        self.sortrev=False
        self.subcollection_keyword = 'operons'
        if prep:
            self.prep()
        
    def __repr__(self):
        return('OperonCollection(%i operons)' %len(self))

#    def set_operons(self):
#        self.operons = sorted(self.values(),key=self.sortkey,reverse=self.sortrev)
    
    def find_operon_by_gene(self,gene_query):
        # Return the first operon that contains a given gene name or gene object
        for operon in self:
            for gene in operon.itergenes(core=True):
                if (type(gene_query) == str and gene.name == gene_query) or (type(gene_query) == Gene and gene.name == gene_query.name):
                    return operon
    
    def prep(self):
        self.set_flanks()
        self.calc_avg_COG()
    
    def set_flanks(self):
        # For each operon, find it's flanking operon (on the same strand)
        # Go through the operons twice, once for each strand
        # Only works for operons that are already on the same strand (operon formation; island method)
        for strand in ['+','-']:
            prev_operon = None
            for operon in self.itersubset(sorted=True):
                if prev_operon: 
                    if operon.strand == strand and prev_operon.strand == strand:
                        dist = operon.start - prev_operon.end
                        operon.left_flank = prev_operon
                        prev_operon.right_flank = operon
                        operon.left_dist = dist
                        prev_operon.right_dist = dist
                else:
                    operon.left_flank = None
                    operon.left_flank = None
                prev_operon = operon
            operon.right_flank = None
            operon.right_dist = None
            
    def clear_operon_flanks(self):
    # Remove any remaining left/right flank references
        for operon in self:
            for attr in ['left_flank','right_flank']:
                if hasattr(operon,attr):
                    delattr(operon,attr)
            
    def calc_avg_COG(self):
        # calculate and set the average COG over all cores genes within the collection
        all_cog = []
        for operon in self:
            if hasattr(operon,'COG'):
                all_cog.extend(operon.COG)
        self.all_cog = all_cog
        if all_cog == []:
            self.COG_avg = None
            self.COG_std = None
        elif len(all_cog) == 1:
            self.COG_avg = all_cog[0]
            self.COG_std = None
        else:
            self.COG_avg = np.average(all_cog)
            self.COG_std = np.std(all_cog)
            
    def add_descr(self,descr_dict):
        # Add the relevant descriptions, relevant for subclusters
        descr = descr_dict[self.collection_type]
        if descr == 'genome':
            descr = 'Genome %s' %(operons_group[0].genome.name)
        elif descr == 'enzyme':
            descr = 'Common enzyme class %s' %group_name
            
    def set_common_domains(self,min_occurrence=0.5):   
        domains = {}
        for operon in self:
            for domain in set(operon.all_domains): # a set is used so that each domain is only added once per operon
                if domain not in domains:
                    domains[domain] = 0
                domains[domain] += 1
        self.all_domains = domains
        rel_domains = {}
        for domain in domains:  
            if float(domains[domain])/len(self) >= min_occurrence:
                rel_domains[domain] = domains[domain]
        self.rel_domains = rel_domains
                
    def set_operon_precursors(self):
        # Store the precursors seperately
        prec_dict = {}
        for operon in self:
            precursors = []
            for gene in operon.itergenes(SVM_hit=True,core=True):
                precursors.append(gene)
                if gene.name in prec_dict:
                    logger.warning('During setting precursors with gene name %s, gene was found twice' %(gene.name))
                prec_dict[gene.name] = operon
            setattr(operon,'precursors',precursors)
        self.precursor_dict = prec_dict

    def set_gene_to_operon(self):
        # Set a link for each core gene pointing to the operon it's part of
        for operon in self.itersubset():
            for gene in operon.itergenes(core=True):
                if hasattr(gene, 'operon') and gene.operon != None:
                    raise ValueError('Double link to operon found! Operon: %s, Gene: %s' %(operon,gene))
                gene.operon = operon

    def set_pairs(self,keywords):
        for operon in self:
            pairs = Container()
            operon.pairs = pairs
            for k in keywords:
                setattr(operon.pairs,k,{})
                

class Scaffold(GeneCollection):
    def __repr__(self):
        return 'Scaffold(name=%s,genome=%s,nr_genes=%i,sequence_length=%i)' %(self.name,self.genome,len(self),len(self.seq))
    
class Genome(CollectionCollection):
    def __init__(self,*args,**kwargs):
        D.__init__(self,*args,**kwargs)
        self.sortkey = sortx_name
        self.sortrev = False
        self.subcollection_keyword = 'scaffolds'
        for item in self:
            if hasattr(self,'active') and self.active:
                item.genome = self
            else:
                item.genome = self.name
        
class Gene(object):
    def __init__(self,**kwargs):
        # The operon the gene is part of, if any
        self.operon = None
        # The genome the gene is part of
        self.genome = None
        # The scaffold the gene is on
        self.scaffold = None
        
        # Any keywords passed here
        self.__dict__.update(kwargs)
        
    def __len__(self):
        return int(self.end)-int(self.start)
        
    def __repr__(self):
        text = 'name=%s,start=%s,end=%s' %(self.name,self.start,self.end)
        if hasattr(self,'scaffold') and hasattr(self.scaffold,'name'):
            text += ',' + self.scaffold.name
        if hasattr(self,'genome') and hasattr(self.genome,'name'):
            text += ',' + self.genome.name
        return 'Gene(%s)' %(text)
    
    def update(self,**kwargs):
        self.__dict__.update(**kwargs)
            
    def print_short(self):
        out = '%s %s %s %s' %(self.name,self.start,self.end,self.strand)
        if hasattr(self,'COG'):
            out += ' %.3f' %self.COG
        logger.debug(out)
        
    def write_fasta(self,filename=False):
        text = '>%s\n%s' %(self.name,self.transl)
        if filename:
            with open(filename,'w') as handle:
                handle.write(text)
        return(text)

    def switch_inactive(self):
        if self.active:
            if hasattr(self,'scaffold'):
                scaffold = self.scaffold.name
                self.scaffold = scaffold
            if hasattr(self,'genome'):
                genome = self.genome.name
                self.genome = genome
            for flank in ['left_flank','right_flank']:
                if hasattr(self,flank) and getattr(self,flank) != None:
                    fl = getattr(self,flank)
                    setattr(self,flank,fl.name)
        self.active=False
    
    def switch_active(self,collection):
        if self.active == False:
            scaffold = collection[self.scaffold]
            genome = collection[self.genome]
            self.scaffold = scaffold
            self.genome = genome
            for flank in ['left_flank','right_flank']:
                if hasattr(self,flank):
                    flank_name = getattr(self,flank)
                    if flank_name != None:
                        flank_obj = collection[flank_name]
                        setattr(self,flank,flank_obj)
        self.active=True

    def get_seq(self,reference_obj=False,seq=False):
        if seq:
            return(seq[self.start:self.end])
        elif type(self.scaffold) == str and reference_obj:
            scaffold = reference_obj[self.scaffold]
        else:
            scaffold = self.scaffold
        return scaffold.seq[self.start:self.end]


