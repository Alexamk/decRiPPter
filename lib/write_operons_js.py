# License: GNU Affero General Public License v3 or later

from log import return_logger
logger = return_logger(__name__, False)

def write_operon_dict(operons,mibig_dict,attributes,path,**kwargs):

    def write_operon(operon,handle):
        handle.write('  "%s": { \n' %operon.name)
        attr_text = ''
        for attr in attributes:
            if hasattr(operon,attr):
                value = getattr(operon,attr)
            elif attr == 'length':
                value = len(operon)
            elif attr == 'length_core':
                if hasattr(operon,'core_genes'):
                    value = len(operon.core_genes)
                else:
                    value = len(operon)
            elif attr == 'nr_antismash' and hasattr(operon,'overlaps_antismash'):
                value = len(operon.overlaps_antismash)
            else:
                value = ''
            if 'genes' in attr:
                value = [g.name for g in value]
            if type(value) == str:
                value = '"%s"' %value
            elif value == None:
                value = 'null'
            elif type(value) == bool:
                value = str(value).lower()
            if value not in [[],'',0]: # These are considered the default version
                attr_text += '    "%s": %s,\n' %(attr,value)
        attr_text.rstrip(',')
        handle.write(attr_text)
        handle.write('  },\n')
        
    with open(path,'w') as handle:
        handle.write('var Operons = {\n')
        for operon in operons.itersubset(**kwargs):
            write_operon(operon,handle)
        for operon in mibig_dict.itersubset(**kwargs):
            write_operon(operon,handle)
        handle.write('}')
                
    
def write_group_dict(group_collections_by_type,path,**kwargs):
    with open(path,'w') as handle:
        handle.write('var Groups = {\n')
        for category in group_collections_by_type:
            collections = group_collections_by_type[category]
            for coll_name,coll in collections.items():
                handle.write('  "%s": ' %coll_name)
                operon_names = []
                for operon in coll.itersubset(**kwargs):
                    operon_names.append(operon.name)
                handle.write('[%s],\n' %', '.join(['"%s"' %i for i in operon_names]))
        handle.write('}')
            


