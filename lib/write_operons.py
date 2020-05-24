# License: GNU Affero General Public License v3 or later

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature,FeatureLocation

import os
import datetime

from log import return_logger
logger = return_logger(__name__, False)

def write_genbank(operons):
    '''Write a genbank file per operon'''
    date = datetime.datetime.now()
    descr = 'decRiPPter predicted RiPP gene cluster'
    for operon in operons:
        logger.debug('Writing gbk file for gene cluster %s' %operon.name)
        annotation = {'keywords':['RiPP','decRiPPter']} 
        annotation['source'] = operon.genome.name
        annotation['organism'] = operon.genome.name
        annotation['date'] = date.strftime('%d-%b-%Y').upper()
        annotation['data_file_division'] = 'BCT'
        features = []
        for gene in operon.itergenes(sorted=True):
            loc = FeatureLocation(gene.start-operon.extra_start + 1,gene.end-operon.extra_start) # +1 to adjust for pythonic 0-index 
            quals = {'SVM_predicted_precursor':gene.SVM_hit, 'protein_domains': gene.domains, 'translation':gene.transl,'pseudo':gene.pseudo}
            
            if hasattr(gene,'protein_id'):
                quals['protein_id'] = gene.protein_id
            if hasattr(gene,'locus_tag'):
                quals['locus_tag'] = gene.locus_tag
            if hasattr(gene,'gene'):
                quals['gene'] = gene.gene
            if hasattr(gene,'gene_synonym'):
                quals['gene_synonym'] = gene.gene_synonym
            if gene.strand == '+':
                strand = 1
            else:
                strand = -1
            feat = SeqFeature(type='CDS',location=loc,strand=strand,id=gene,qualifiers=quals)
            
            features.append(feat)
        dna_operon = Seq(operon.get_seq(extra=True),IUPAC.unambiguous_dna)
        # Check name length: locus identifier + length of nr of bases can't be longer than 27 chars
        newname = operon.name
        current_len = len(newname) + len(str(len(dna_operon)))
        while current_len > 27:
            newname = newname[1:]
            current_len = len(newname) + len(str(len(dna_operon)))
        if newname != operon.name:
            logger.debug('Changing gbk name %s to %s' %(operon.name,newname))
            pass
        sr = SeqRecord(dna_operon,id=operon.name,name=newname,features=features,description=descr,annotations=annotation)
        
        SeqIO.write(sr,operon.genbank_path,'genbank')
        
def write_fasta_all(operons):
    for operon in operons:
        _ = operon.write_fasta(filename=operon.fasta_path)
