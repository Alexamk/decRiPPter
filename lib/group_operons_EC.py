# License: GNU Affero General Public License v3 or later

# Find all operons belonging to specific subcategories based on pfam2ec and TIGRFAM2ec links

from Genes import OperonCollection, CollectionCollection

from log import return_logger
logger = return_logger(__name__, False)

def main(operons,domain2ec):
    categories_ec = {
    '2.1': 'Transferring_One-Carbon_Groups',
    '2.2': 'Transferring_Aldehyde_or_Ketonic_Groups',
    '2.3': 'Acyltransferases',
    '2.4': 'Glycosyltransferases',
    '2.5': 'Transferring_Alkyl_or_Aryl_Groups_Other_than_Methyl_Groups',
    '2.6': 'Transferring_Nitrogenous_Groups',
    '2.7': 'Transferring_Phosphorus-Containing_Groups',
    '2.8': 'Transferring_Sulfur-Containing_Groups'}
    
    categories_domains = {'Halogenase':['PF04820']}
    
    # First convert ec categories to domain categories - i.e. a dict containing a category name and all allowed domains.
    for ec_cl in categories_ec:
        ec_cl_full = 'EC:%s' %ec_cl
        category_name = categories_ec[ec_cl]
        domains = find_ec(domain2ec,ec_cl_full)
        categories_domains[category_name] = domains
    # Now assign the groups to the operons and the operoncollection
    operon_collections = assign_groups_to_operons(operons,categories_domains)
    return operon_collections
    
def find_ec(domain2ec,ec_q):
    out = set()
    for dom in domain2ec:
        ecs = domain2ec[dom]
        for ec in ecs:
            if ec_q in ec:
                out.add(dom)
    return(sorted(list(out)))
    
def assign_groups_to_operons(operons,categories_domains):
    operons_per_category = {}
    for operon in operons:
        categories = []
        for category in categories_domains:
            domains = categories_domains[category]
            if any([dom in operon.all_domains for dom in domains]):
                categories.append(category)
                if category not in operons_per_category:
                    operons_per_category[category] = []
                operons_per_category[category].append(operon)
        operon.EC_groups = categories
    collections_per_group = {}
    for category in operons_per_category:
        operons_group = operons_per_category[category]
        collection = OperonCollection(dict([(operon.name,operon) for operon in operons_group]),collection_type='EC',name=category,\
                                      descr='Operons containing a common enzyme class',prep=False)
        collections_per_group[category] = collection
    return collections_per_group
    


    
'''Suggested interesting enzymes belonging to EC classes
2.1 Transferring One-Carbon Groups
2.2 Transferring Aldehyde or Ketonic Groups
2.3 Acyltransferases
2.4 Glycosyltransferases
2.5 Transferring Alkyl or Aryl Groups, Other than Methyl Groups
2.6 Transferring Nitrogenous Groups
2.7 Transferring Phosphorus-Containing Groups
2.8 Transferring Sulfur-Containing Groups

Additionally, for halogen transfer

PF04820	Trp_halogenase	Tryptophan halogenase

Possibly also:
https://www.nature.com/articles/nchembio.1649
--> PF13640
--> PF05721
'''



