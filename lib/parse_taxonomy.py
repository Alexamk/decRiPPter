# License: GNU Affero General Public License v3 or later

def parse_nodes(path):
    f = open(path)
    id2parent = {}
    parent2ids = {}
    id2rank = {}
    for l in f:
        tabs = l.strip().split('\t')
        taxid = tabs[0]
        parent_taxid = tabs[2]
        rank = tabs[4]
        id2parent[taxid] = parent_taxid
        id2rank[taxid] = rank
        if not parent_taxid in parent2ids:
            parent2ids[parent_taxid] = [taxid]
        else:
            parent2ids[parent_taxid].append(taxid)
    f.close()    
    return(id2parent,parent2ids,id2rank)
    

def parse_assembly(path='assembly_summary_genbank.txt'):
    f = open(path)
    acc2data = {}   
    name2acc = {} 
    acc2name = {}
    f.readline()
    f.readline()
    for l in f:
        tabs = l.strip().split('\t')
        assembly_accession = tabs[0]
        name = tabs[7]
        taxid = str(tabs[6])
        strain = tabs[8]
        assembly_level = tabs[11]
        ftp = tabs[19][6:]
        refseq_category = tabs[4]
        acc2data[assembly_accession] = (taxid,assembly_level,strain,ftp,refseq_category)
        acc2name[assembly_accession] = name
        if name not in name2acc:
            name2acc[name] = [assembly_accession]
        else:
            name2acc[name] += [assembly_accession]
    f.close()
    return(acc2data,name2acc,acc2name)
    
def find_taxid_in_parent(acc2data,taxq,id2parent,assembly_level_req=["complete genome"],refseq_category_req=False):
    have = []
    donthave = []
    notintaxonomy = []
    taxq = str(taxq)
    for acc in acc2data:
        data = acc2data[acc]
        assembly_level = data[1]
        if assembly_level_req and assembly_level.lower() not in assembly_level_req : 
            continue
        refseq_category = data[-1]
        if refseq_category_req and refseq_category.lower() not in refseq_category_req:
            continue
        taxid = data[0]
        parentid = taxid        
        while parentid != '1' and parentid != taxq:
            if parentid in id2parent:
                parentid = id2parent[parentid]
            else:
                notintaxonomy.append((parentid,acc))
                break
        if parentid == 1:
            donthave.append(acc)
        if parentid == taxq:
            have.append(acc)
    return(have,donthave,notintaxonomy)
