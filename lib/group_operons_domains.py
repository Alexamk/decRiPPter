# Calculate the jaccard index of found operons 
# Use it to form networks/groups

import os
import time
from lib import store_pickle, load_pickle, run_mcl, read_mcl, pairs_to_groups, write_ssn, Edge
from log import return_logger
from Genes import OperonCollection
from multiprocessing_wrappers import operator

logger = return_logger(__name__, False)

def main(operons,path,mibig_dict,settings):
    logger.debug('Calculating gene cluster pairwise distances')
    ssn_file = path + 'operon_pairs.ssn'
    ssn_file_filtered = path + 'operon_pairs_filtered.ssn'
    pfile = os.path.join(settings['paths']['pickles'],'jaccard_dict.pkl')
    if settings['calculate_jaccard']:
        t0 = time.time()
        organized_operons,set_dict = organize_operons_by_domains(operons,mibig_dict)
        t1 = time.time()
        jaccard_dict_domains = calculate_all_jaccard_prefiltered(organized_operons,set_dict,settings['jaccard_cutoff'])
        t2 = time.time()
        jaccard_dict = convert_domainscore_to_operonscore(organized_operons,jaccard_dict_domains)
        t3 = time.time()
        logger.debug('Jaccard prefilter calculation time: %.4f' %(t1-t0))
        logger.debug('Jaccard actual calculation time: %.4f' %(t2-t1))
        logger.debug('Jaccard conversion time: %.4f' %(t3-t2))
        if not os.path.isfile(ssn_file):
            write_ssn(jaccard_dict,ssn_file,i=0)
        store_pickle(jaccard_dict,pfile)
    else:
        jaccard_dict = load_pickle(pfile)
    if settings['jaccard_mcl']:
        mcl_file = path + 'jaccard_groups.mcl'
        if not os.path.isfile(mcl_file):
            r1 = run_mcl(ssn_file,mcl_file,settings['cores'])
        groups = read_mcl(mcl_file,settings['min_group_size'])
    else:
        groups = pairs_to_groups(jaccard_dict,settings['min_group_size'])
    collection_name = 'jaccard_groups'
    logger.debug('Setting gene cluster pairs')
    group_names,operon_collections = assign_groups_operons(operons,groups,mibig_dict,collection_name)
    filtered_pairs = filter_operon_pairs(jaccard_dict,groups)
    write_ssn(filtered_pairs,ssn_file_filtered,i=0)
    set_jaccard_operon_pairs(operons,mibig_dict,filtered_pairs)
    return group_names,operon_collections,jaccard_dict,filtered_pairs
   
def calculate_all_jaccard_prefiltered(operons_grouped,set_dict,cutoff):
    c = 0
    outd = {}
    logger.debug('Calculating pairwise jaccard score for %i identical domain sets' %len(operons_grouped))
    for domains1 in operons_grouped:
        c += 1
#        if c % 100 == 0:
#            print('At set %i out of %i' %(c,len(operons_grouped)))
#            print('Breaking')
#            return outd
        for domains2 in operons_grouped:
            if domains1 == domains2:
                continue
            jaccard,domains_overlap = calc_jaccard_from_domains(domains1,domains2,set_dict)
            if jaccard >= cutoff:
                if domains1 not in outd:
                    outd[domains1] = {}
                outd[domains1][domains2] = jaccard,list(domains_overlap)
    return outd
    
def calc_jaccard_from_domains(domains1,domains2,set_dict):
    nr_tag1 = len(domains1)
    nr_tag2 = len(domains2)
    if nr_tag1 == 0 or nr_tag2 == 0:
        return 0.0,None
    intersection = set_dict[domains1] & set_dict[domains2]
    overlap = len(intersection)
    jaccard = float(overlap) / (nr_tag1 + nr_tag2 - overlap)
    return jaccard,intersection

def calc_jaccard_from_domains2(domains1,domains2):
    nr_tag1 = len(domains1)
    nr_tag2 = len(domains2)
    if nr_tag1 == 0 or nr_tag2 == 0:
        return 0.0,None
    intersection = set(domains1) & set(domains2)
    overlap = len(intersection)
    jaccard = float(overlap) / (nr_tag1 + nr_tag2 - overlap)
    return jaccard,intersection
    
def calc_jaccard_from_domains3(domains1,domains2,set_dict,done_dict):
    if (domains1,domains2) in done_dict:
        return done_dict[domains1,domains2]
        
    nr_tag1 = len(domains1)
    nr_tag2 = len(domains2)
    if nr_tag1 == 0 or nr_tag2 == 0:
        return 0.0,None
    intersection = set_dict[domains1] & set_dict[domains2]
    overlap = len(intersection)
    jaccard = float(overlap) / (nr_tag1 + nr_tag2 - overlap)
    done_dict[domains1,domains2] = done_dict[domains2,domains1] = jaccard,intersection
    return jaccard,intersection

def organize_operons_by_domains(operons,mibig_operons):
    outd = {}
    set_dict = {}
    for itr in mibig_operons,operons:
        for operon in itr:
            domains = tuple(operon.unique_domains)
            if domains not in outd:
                outd[domains] = []
                set_dict[domains] = set(domains)
            outd[domains].append(operon.name)
    logger.debug('%i identical domain sets found out of %i gene clusters' %(len(outd),len(operons)))
    return outd,set_dict
    
def convert_domainscore_to_operonscore(operons_by_domains,jaccard_dict_domains):
    outd = {}
    for key in jaccard_dict_domains:
        vs = jaccard_dict_domains[key]
        for key2 in vs:
            jaccard,intersection = vs[key2]
            
            operons1 = operons_by_domains[key]
            operons2 = operons_by_domains[key2]
            
            for operon1 in operons1:
                for operon2 in operons2:
                    if operon1 < operon2:
                        if operon1 not in outd:
                            outd[operon1] = {}
                        outd[operon1][operon2] = jaccard,intersection
    # Now add pairs for identical domainsets
    for key in operons_by_domains:
        if key == tuple():
            continue
        operons = operons_by_domains[key]
        
        for operon1 in operons:
            for operon2 in operons:
                if operon1 < operon2:
                    if operon1 not in outd:
                        outd[operon1] = {}
                    outd[operon1][operon2] = 1.0,list(key)
    return outd

            
def assign_groups_operons(operons,groups,mibig_dict,collection_name):
    group_names = []
    operons_per_group = {}
    i = 1
    for group in groups:
        group_name = '%s_%i' %(collection_name,i)
        group_names.append(group_name)
        operons_per_group[group_name] = []
        for operon in group:
            operon_obj = operons[operon]
            if operon_obj == None:
                # Get it from MIBiG
                operon_obj = mibig_dict[operon]
                if operon_obj == None:
                    logger.error('%s is not found in gene cluster database or mibig dictionaries' %operon)
                    raise KeyError('%s is not found in gene cluster database or mibig dictionaries' %operon)
            # Now update the operon itself, and add it to the dictionary
            if not hasattr(operon_obj,collection_name):
                # Add the group name as an attribute to the operon
                setattr(operon_obj,collection_name,group_name)
            else:
                logger.error('Gene cluster %s part of multiple groups' %operon)
                raise ValueError('Gene cluster %s part of multiple groups' %operon)
            operons_per_group[group_name].append(operon_obj)
        if not all([hasattr(operon,'mibig') for operon in operons_per_group[group_name]]): # At least one operon from the run
            i += 1
        else:
            del(operons_per_group[group_name])
    # Rework the operons per group to be operoncollection objects
    operons_per_group_coll = {}
    for group in operons_per_group:
        operons_group = operons_per_group[group]
        operon_collection = OperonCollection(dict([(operon.name,operon) for operon in operons_group]),prep=False,collection_type=collection_name,name=group,descr='Jaccard index of protein domains')
        operons_per_group_coll[group] = operon_collection
    return group_names,operons_per_group_coll

def filter_operon_pairs(jaccard_dict,groups):
    filtered = {}
    for group in groups:
        allowed = set(group)
        for operon in group:
            if operon not in jaccard_dict:
                continue
            pairs = jaccard_dict[operon]
            for operon2 in pairs:
                if operon2 in allowed:
                    if operon not in filtered:
                        filtered[operon] = {}
                    filtered[operon][operon2] = jaccard_dict[operon][operon2]
    return filtered
            
def set_jaccard_operon_pairs(operons,mibig_dict,filtered_pairs):
    for operon in filtered_pairs:
        op_obj = operons[operon]
        if op_obj == None:
            op_obj = mibig_dict[operon]
        pairs = filtered_pairs[operon]
        for operon2 in pairs:
            op_obj2 = operons[operon2]
            if op_obj2 == None:
                op_obj2 = mibig_dict[operon2]
            e = Edge(pairs[operon2][0],pairs[operon2][1])
            op_obj.pairs.jaccard[operon2] = e
            op_obj2.pairs.jaccard[operon] = e
        
        
        
        
        
        
        

