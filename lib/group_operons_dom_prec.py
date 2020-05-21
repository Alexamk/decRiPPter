import os
import numpy as np
from lib import run_mcl, read_mcl, pairs_to_groups, write_ssn, Edge
from Genes import OperonCollection

from log import return_logger
logger = return_logger(__name__, False)

def main(operons,mibig_dict,prec_pairs,jaccard_pairs,settings):
    path = settings['paths']['jaccard_prec_path']
    if not os.path.isdir(path):
        os.mkdir(path)
    ssn_file = os.path.join(path,'pairs.ssn')
    mcl_file = os.path.join(path,'groups.mcl')
    pairs = find_overlap_pairs(prec_pairs,jaccard_pairs)
    if not os.path.isfile(ssn_file):
        write_ssn(pairs,ssn_file)
    if not os.path.isfile(mcl_file):
        run_mcl(ssn_file,mcl_file,settings['cores'])
    groups = read_mcl(mcl_file)
    filtered_pairs = filter_operon_pairs(pairs,groups)
    set_pairs(filtered_pairs,operons,mibig_dict)
    group_names,group_collections = assign_groups_operons(operons,groups,mibig_dict,'precursor_jaccard_groups')
    return group_collections
    
def find_overlap_pairs(prec_pairs,jaccard_pairs):
    operon_pairs = {}
    for operon in jaccard_pairs:
        pairs = jaccard_pairs[operon]
        for pair in pairs:
            if operon in prec_pairs and pair in prec_pairs[operon]:
                prec_weights = []
                for bitscore, pid in prec_pairs[operon][pair].values():
                    prec_weights.append(pid)
                prec_weight = np.average(prec_weights)/100.0
                jac_weight = pairs[pair][0]
                new_weight = (prec_weight + jac_weight) / 2

                if operon not in operon_pairs:
                    operon_pairs[operon] = {}
                operon_pairs[operon][pair] = new_weight
    return operon_pairs
    
def set_pairs(operon_pairs,operons,mibig_dict):
    for op1 in operon_pairs:
        for op2 in operon_pairs[op1]:
            weight = operon_pairs[op1][op2]
            e = Edge(weight,[])
            if op1 in operons:
                operons[op1].pairs.precursor_jaccard[op2] = e
            else:
                mibig_dict[op1].pairs.precursor_jaccard[op2] = e
            if op2 in operons:
                operons[op2].pairs.precursor_jaccard[op1] = e
            else:
                mibig_dict[op2].pairs.precursor_jaccard[op1] = e

def filter_operon_pairs(pair_dict,groups):
    filtered = {}
    for group in groups:
        allowed = set(group)
        for operon in group:
            if operon not in pair_dict:
                continue
            pairs = pair_dict[operon]
            for operon2 in pairs:
                if operon2 in allowed:
                    if operon not in filtered:
                        filtered[operon] = {}
                    filtered[operon][operon2] = pairs[operon2]
    return filtered
    
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
                    raise KeyError('%s is not found in operon database or mibig dictionaries' %operon)
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
        operon_collection = OperonCollection(dict([(operon.name,operon) for operon in operons_group]),prep=False,\
                                            collection_type=collection_name,name=group,descr='Overlap in jaccard and precursor pairs')
        operons_per_group_coll[group] = operon_collection
    return group_names,operons_per_group_coll
    
