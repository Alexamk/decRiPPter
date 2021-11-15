# License: GNU Affero General Public License v3 or later

# Script to write HTML pages for decRiPPter

import numpy as np
import os

from lib import read_template, sortdictkeysbyvalues

from log import return_logger
logger = return_logger(__name__, False)

def main(settings,operons,collections_by_type,headers,domain_descr,genome_dict):
    logger.debug('Writing index file')
    index_file = settings['paths']['index_file']
    index_text = write_index(settings,operons,collections_by_type,headers,domain_descr,genome_dict,index_file)
    logger.debug('Writing entry files')
    write_all_entries(settings,collections_by_type,domain_descr,index_file)

def write_index(settings,operons,collections_by_type,headers,domain_descr,genome_dict,out_path):
    html_path = os.path.abspath(settings['paths']['html'])
    outfolder = settings['paths']['out']
    style_sheet_path = os.path.relpath(os.path.join(html_path,'html_index.css'),start=outfolder)
    template_path = os.path.join(html_path,'Index_template.html')
    style_sheet_text = ' '*4 + '<link rel="stylesheet" href="%s">\n' %style_sheet_path
    script_text = '<script src="%s"></script>\n' %os.path.relpath(settings['paths']['jquery_path'],start=outfolder)
    script_text += '<script src="%s"></script>\n' %os.path.relpath(settings['paths']['tablesorter_path'],start=outfolder)
    script_text += '<script src="%s"></script>\n' %os.path.relpath(settings['paths']['operons_js_file'],start=outfolder)
    script_text += '<script src="%s"></script>\n' %os.path.relpath(settings['paths']['operon_group_js_file'],start=outfolder)
    script_text += '<script src="%s"></script>\n' %os.path.relpath(settings['paths']['filter_scripts'],start=outfolder)
    
    script_path = os.path.relpath(os.path.join(html_path,'Index.js'),start=outfolder)
    script_text += '<script src="%s"></script>\n' %script_path
    table_ids = []
    
    template = read_template(template_path)
    
    menu_text = ''
    content_text = ''
    # Get the menu headers
    for group_type in collections_by_type:
        logger.debug('Writing index section for group %s' %(group_type))
        collections = collections_by_type[group_type]
        header = headers[group_type]
        
        menu_text += '''<button class="tablinks" onclick="opentab(event, '%s')">%s</button>\n''' %(group_type,header)
        coll_text = ''
        # Get the content per group type; each collection is an operon
        for coll_name,collection in collections.items():
            logger.debug('Writing table row for collection %s' %coll_name)
            if group_type == 'EC':
                coll_name_table = coll_name.replace('_',' ')
            elif group_type == 'genome':
                coll_name_table = collection.realname
            else:   
                coll_name_table = coll_name
            domain_text,nr_domains = get_common_domain_text(collection,domain_descr,mode='text')
            nr_operons_passed_filter = len([op for op in collection if (hasattr(op,'filtered') and op.filtered == False)])
            nr_mibig_operons = len([op for op in collection if hasattr(op,'mibig') and op.mibig])
            nr_antismash = len([op for op in collection if hasattr(op,'overlaps_antismash') and op.overlaps_antismash != {}])
            nr_kripp = len([op for op in collection if len(op.kripp) > 0])
            scores = [op.score for op in collection if hasattr(op,'score')]
            if len(scores) > 0:
                score_text = '%.2f +- %.2f' %(np.average(scores),np.std(scores))
            else:
                score_text = 'N\A'
            
            if domain_text != '<p>No common domains found</p>':
                domain_text = '''%i<br><button class="collapsible">Show domains</button>
            <div class="content">
            <p> 
            %s
            </p>
            </div>''' %(nr_domains,domain_text)
            COG_text = get_COG_text(collection)
            coll_text += '''        <tr id="%s">
            <td><a href='%s'>%s</a></td>
            <td>%s</td>
            <td>%i</td>
            <td>%i</td>
            <td>%i</td>
            <td>%i</td>
            <td>%s</td>
            <td>%i</td>
            <td>%s</td>
        </tr>''' %(coll_name,os.path.relpath(collection.path,start=outfolder),coll_name_table,score_text,len(collection),nr_mibig_operons,\
                   nr_antismash,nr_kripp,COG_text,nr_operons_passed_filter,domain_text)
        
        table_text = '''<div class="tabcontent" id='%s'>
    <table class="tablesorter" id="%s_table">
        <thead>
        <tr>
            <th>Entry name</th>
            <th>Average score</th>
            <th>Number of operons</th>
            <th>Number of MIBiG operons</th>
            <th>Number of operons with antismash overlap</th>
            <th>Nr of operons with known RiPP domains</th>
            <th>Average COG score</th>
            <th>Nr of operons passed filter</th>
            <th>Common domains</th>
        </tr>
        </thead>
    <tbody>
    %s
    </tbody>
    </table>
    </div>\n''' %(group_type,group_type,coll_text)
        content_text += table_text
        
        table_ids.append('%s_table' %group_type)
        
    info_text = '''<b>Project name:</b> {0} <br>
    <b>Number of genomes:</b> {1}  <br>
    <b>Taxid given:</b> <br>
    '''.format(settings['name'],len(genome_dict))
    
    info_text = '''<aside>
%s
</aside>''' %info_text
    
    script_text  += write_script_load_tablesorter(table_ids)
    
    out_text = template[0] + style_sheet_text + template[1] + menu_text + template[2] + info_text + template[3] + content_text +  template[4] + script_text + template[5]
    with open(out_path,'w') as outfile:
        outfile.write(out_text)
    return out_text
        
def write_all_entries(settings,collections_by_type,domain_descr,index_file_all):
    # Write an entry page for each cluster/group found
    html_path = os.path.abspath(settings['paths']['html'])
    out_path = settings['paths']['entries']
    template_path = os.path.join(html_path,'Entry_template.html')
    table_ids = ["main_table"]
    tablesorter_script = write_script_load_tablesorter(table_ids)
    template = read_template(template_path)
    for group_type in collections_by_type:
        collections = collections_by_type[group_type]
        type_path = os.path.join(out_path,group_type)
        back_to_index_path = '<a href="%s">Back to index</a>' %os.path.relpath(index_file_all,start=type_path)
        style_sheet_path = os.path.relpath(os.path.join(html_path,'html_entry.css'),start=type_path)
        script_path = os.path.relpath(os.path.join(html_path,'Entry.js'),start=type_path)
        style_sheet_text = ' '*4 + '<link rel="stylesheet" href="%s">\n' %style_sheet_path
        script_text = '<script src="%s"></script>\n' %os.path.relpath(settings['paths']['jquery_path'],start=type_path)
        script_text += '<script src="%s"></script>\n' %os.path.relpath(settings['paths']['tablesorter_path'],start=type_path)
        script_text += '<script src="%s"></script>\n' %os.path.relpath(settings['paths']['operons_js_file'],start=type_path)
        script_text += '<script src="%s"></script>\n' %os.path.relpath(settings['paths']['operon_group_js_file'],start=type_path)
        script_text += '<script src="%s"></script>\n' %os.path.relpath(settings['paths']['filter_scripts'],start=type_path)
        script_text += '<script src="%s"></script>\n' %script_path

        for name_coll,collection in collections.items():
            page_id_text = '    <html_id>%s</html_id>\n' %name_coll
            table_ids = []
            info_text = write_entry_info(collection,type_path)
            table_text = ''
            domain_text,nr_domains = get_common_domain_text(collection,domain_descr,mode='table')
            for operon in collection:
                table_text += write_operon_td(operon,os.path.dirname(collection.path),settings)
            # Assemble the text
            out_text = template[0] + style_sheet_text + page_id_text + template[1] + collection.name + '\n' + \
                       template[2] + back_to_index_path + template[3] + info_text + template[4] + domain_text + \
                       template[5] + table_text + template[6] + script_text + tablesorter_script + template[7]
           
            with open(collection.path,'w') as handle:
                handle.write(out_text)
         
def get_COG_text(collection):
    if collection.COG_avg == None:
        COG_text = 'N\\A'
    elif collection.COG_std == None:
        COG_text = '%.2f' %collection.COG_avg
    else:
        COG_text = '%.2f +- %.2f' %(collection.COG_avg,collection.COG_std)
    return COG_text
         
def write_entry_info(collection,type_path):
    COG_text = get_COG_text(collection)
    
    text = '''        <b>Entry name: </b>%s<br>
        <b>Common denominator: </b>%s<br>
        <b>Number of operons: </b><span id="nr_operons">%s</span><br>
        <b>Average COG: </b>%s<br>\n''' %(collection.name,collection.descr,len(collection),COG_text)
    if collection.collection_type == 'precursor_groups':
        text += '<br>\n<a href="%s">Precursor alignment</a>\n<br>\n' %os.path.relpath(collection.muscle_path,start=type_path)
    return text

def write_script_load_tablesorter(ids):
    out = ''
    for id in ids:
        out += '''$("#%s").tablesorter();\n''' %id
    out = '''    <script>
    $( document ).ready(function() {
        $(function() {
            %s
        });
    });
    </script>\n''' %out
    return out

def write_operon_td(operon,coll_path,settings):
    domain_text = write_domain_text_operon(operon)
    if hasattr(operon,'mibig'):
        operon_path = 'https://mibig.secondarymetabolites.org/repository/%s/index.html' %operon.name
        operon_score = 'N\\A'
        passed_filter = 'N\\A'
    else:
        operon_path = os.path.relpath(os.path.join(operon.path,operon.name+'.html'),coll_path)
        operon_score = operon.score
        passed_filter = True
    operon_kripps = len(operon.kripp)
    if settings['antismash_parse']:
        if not hasattr(operon,'overlaps_antismash') or operon.overlaps_antismash == {}:
            antismash_text = 'None found'
        else:
            longest_overlap = 0
            best_cluster = ''
            for name,data in operon.overlaps_antismash.items():
                if data['overlap_type'] == 'total':
                    best_cluster = name
                    break
                else:
                    start,end = data['coords'] 
                    if operon.start < start and operon.end > end:
                        overlap = end - start
                    elif operon.start > start and operon.end > end:
                        overlap = end - operon.start
                    elif operon.start < start and operon.end < end:
                        overlap = operon.end - start
                    if overlap > longest_overlap:
                        longest_overlap = overlap
                        best_cluster = name
            data = operon.overlaps_antismash[name]
            antismash_text = '%s (%s overlap)' %('-'.join(data['cluster_type']),data['overlap_type'])
    else:
        antismash_text = 'N\\A'
    
    COG_text = ''
    if hasattr(operon,'COG_avg') and operon.COG_avg:
        COG_text = '%.2f +- %.2f' %(operon.COG_avg,operon.COG_std)
    else:
        COG_text = 'N\\A'
    
    text='''    <tr id="%s">
        <td><a href="%s">%s</a></td>
        <td>%s</td>
        <td>%s</td>
        <td>%s</td>
        <td>%s</td>
        <td>%s</td>
        <td>%s</td>
        <td>%s</td>
    </tr>\n''' %(operon.name,operon_path,operon.name,operon_score,len(operon),\
                 antismash_text,COG_text,passed_filter,domain_text,operon_kripps)
    return text

def write_domain_text_operon(operon):
    domain_wrapper = '\n' + ' '*12 + '''%i<br><button class="collapsible" >Show domains</button>
        <div class="content">
        <p> 
        REPLACEME
        </p>
        </div>\n''' %(len(operon.unique_domains))
    
    domain_text = ''
    for domain in sorted(operon.all_domains):
        domain_text += domain + '<br>'
    if domain_text == '':
        domain_text = 'No domains found'
    else:
        domain_text = domain_text[:-4] # Strip the rightmost <br>
        domain_text = domain_wrapper.replace('REPLACEME',domain_text)
    return domain_text
    
    

def get_common_domain_text(collection,domain_descr,mode='table'):
    # Common domains from a single entry (e.g. a cluster of operons)
    rel_domains = collection.rel_domains
    if rel_domains == {}:
        return '''<p>No common domains found</p>''',0
    elif mode == 'table':
        domain_text,nr_domains = write_common_domain_table(rel_domains,len(collection),domain_descr)
    elif mode == 'text':   
        domain_text,nr_domains = write_common_domain_text(rel_domains,len(collection),domain_descr)
    return domain_text,nr_domains
    
def write_common_domain_text(rel_domains,total,domain_descr):
    doms_sorted = sortdictkeysbyvalues(rel_domains,reverse=True)
    nr_domains = len(doms_sorted)
    text = '<p>\n'
    for domain in doms_sorted:
        count = rel_domains[domain]
        try:
            descr = domain_descr[domain]
        except KeyError:
            logger.warning('No domain description found for domain %s - please update the domaindescr file' %domain)
            descr = 'Unknown'
        if domain.startswith('PF'):
            link = 'https://pfam.xfam.org/family/%s' %domain
        else:
            link = 'http://tigrfams.jcvi.org/cgi-bin/HmmReportPage.cgi?acc=%s' %domain
        text += 4*' ' + '''<a href='%s'>%s</a> (%i/%i) %s<br>\n''' %(link,domain,count,total,descr)
    return text,nr_domains
    
def write_common_domain_table(rel_domains,total,domain_descr):
    base_text ='''<table>
    <caption>Common protein domains</caption>
    <tr>
        <th>Domain ID</th>
        <th>Domain name</th>
        <th>Times found</th>
    </tr>\n'''
    
    doms_sorted = sortdictkeysbyvalues(rel_domains,reverse=True)
    nr_domains = len(doms_sorted)
    for domain in doms_sorted:
        count = rel_domains[domain]
        try:
            descr = domain_descr[domain]
        except KeyError:
            logger.warning('No domain description found for domain %s - please update the domaindescr file' %domain)
            descr = 'Unknown'
        if domain.startswith('PF'):
            link = 'https://pfam.xfam.org/family/%s' %domain
        else:
            link = 'http://tigrfams.jcvi.org/cgi-bin/HmmReportPage.cgi?acc=%s' %domain
        base_text += '''    <tr>
        <td><a href="%s">%s</a></td>
        <td>%s</td>
        <td>%i/%i</td>
    </tr>\n''' %(link,domain,descr,count,total)
    
    base_text += '</table>\n'
    return base_text,nr_domains

