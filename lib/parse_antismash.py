import os

from lib import load_pickle, open_genbank
from log import return_logger

logger = return_logger(__name__, False)

def read_all_antismash(genomes_folder, genome_dict, antismash_version):
    all_clusters = {}
#    all_biosyn = {}
#    all_cluster_genes = {}
    for genome in os.listdir(genomes_folder):
        logger.debug('Parsing antiSMASH results for genome %s' %genome)
        genome_path = os.path.join(genomes_folder,genome)
        antismash_folder = os.path.join(genome_path,'antismash_minimal')
        result_file = get_antismash_results_file(antismash_folder, antismash_version, genome)
        if result_file == None:
            logger.warning('Not parsing %s antiSMASH results' %genome)
            all_clusters[genome] = {}
            continue
        all_seqs = open_genbank(result_file)
        clusters = find_gbk_clusters(all_seqs, antismash_version)
        all_clusters[genome] = clusters
        
    return all_clusters

def find_gbk_clusters(all_seqs, antismash_version):
    cluster_identifier = get_cluster_identifier(antismash_version)
    logger.debug('Finding BGCs with cluster identifier: %s' %cluster_identifier )

    clusters = {}
    for seq in all_seqs:
        contig_id = seq.id
        if seq.id == 'unknown':
            contig_id = seq.name
        clusters[contig_id] = {}
        for feature in seq.features:
            if feature.type == cluster_identifier:
                coords = int(feature.location.start),int(feature.location.end)
                cluster_nr = get_cluster_number(feature, antismash_version)
                cluster_product = get_cluster_product(feature, antismash_version)
                logger.debug('BGC found: contig: %s; protocluster_number: %s; product: %s' %(contig_id, cluster_nr, cluster_product))                
                clusters[contig_id][cluster_nr] = [cluster_product,coords,{}]

    return clusters
       
def get_cluster_identifier(antismash_version):
    logger.debug('antiSMASH version %s' %antismash_version)
    if  antismash_version == 4:
        return 'cluster'
    elif antismash_version == 5:
        return 'protocluster'
       
def get_cluster_number(feature, antismash_version):
    if antismash_version == 4:
        ID_holder = feature.qualifiers.get('note')
    elif antismash_version == 5:
        ID_holder = feature.qualifiers.get('protocluster_number')
    if not ID_holder:
        logger.error('No cluster ID found in feature {}'.format(feature))
        raise ValueError('No cluster ID found in feature {}'.format(feature))
    cluster_nr = ID_holder[0]
    return cluster_nr
    
def get_cluster_product(feature, antismash_version):
    product_holder = feature.qualifiers.get('product')
    if not product_holder:
        logger.error('No product found in feature {}'.format(feature))
        raise ValueError('No product found in feature {}'.format(feature))
    return product_holder[0]

def get_antismash_results_file(folder, antismash_version, genome_name):
    files = os.listdir(folder)
    if antismash_version == 4:
        eligible_files = [f for f in files if f.endswith('.final.gbk')]
    elif antismash_version == 5:
        eligible_files = [f for f in files if f == '%s.gbk' %genome_name]
    if len(eligible_files) > 1:
        logger.error('More than one antiSMASH final result file found in folder %s' %folder)
        raise ValueError()
    elif len(eligible_files) == 0:
        logger.warning('No antiSMASH results found in folder %s' %folder)
        return 
    return eligible_files[0]

def compare_antismash_operons(antismash_clusters, operons):
    # Update all the operons with any overlaps with antismash clusters
    for operon in operons:
        start = operon.start
        end = operon.end
        operon.overlaps_antismash = {}
        if operon.scaffold.name not in antismash_clusters[operon.genome.name] and\
           operon.scaffold.name.rpartition('.')[0] not in antismash_clusters[operon.genome.name]:
            continue
            # Too short contig, not analyzed by antiSMASH
        try:
            clusters = antismash_clusters[operon.genome.name][operon.scaffold.name]
        except KeyError:
            clusters = antismash_clusters[operon.genome.name][operon.scaffold.name.rpartition('.')[0]]
        for cluster in clusters:
            cluster_type,coords,found = clusters[cluster]
            start_cluster,end_cluster = coords
            cluster_types = cluster_type.split('-')
            # Determine overlap - total overlap or partial overlap?
            if start_cluster <= start and end_cluster >= end:
                #total overlap
                operon.overlaps_antismash[cluster] = dict(cluster_type=cluster_types,coords=coords,overlap_type='total')
                clusters[cluster][2][operon.name] = 1.0
            elif (start_cluster < start and end_cluster > start) or\
                 (start_cluster < end and end_cluster > end):
                 # partial overlap
                 # Calculate partial overlap; give as fraction of operon covered
                if (start_cluster < start and end_cluster > start):
                    nt_overlap = end_cluster - start
                elif (start_cluster < end and end_cluster > end):
                    nt_overlap = end - start_cluster
                part_overlap = float(nt_overlap) / (operon.end - operon.start)
                operon.overlaps_antismash[cluster] = dict(cluster_type=cluster_types, coords=coords, overlap_type='partial', overlap_fraction=part_overlap)
                clusters[cluster][2][operon.name] = part_overlap
                # If the entire antiSMASH cluster is covered, also indicate this
                # This is rare as antiSMASH clusters are typically larger than RiPPTIDE's
                if start_cluster >= start and end_cluster <= end:
                    logging.debug('Entire antiSMASH BGC covered by decRiPPTer gene cluster %s' %operon)
                    operon.overlaps_antismash['antismash_overlapped'] = True


def parse_main(operons,genome_dict,settings):
    genome_path = settings['paths']['genome']
    all_clusters = read_all_antismash(genome_path, genome_dict, settings['antismash_version'])
    compare_antismash_operons(all_clusters,operons)
    return all_clusters

