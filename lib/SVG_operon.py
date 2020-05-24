# License: GNU Affero General Public License v3 or later

# Script based of script by Peter Cimermancic
# Adapted by Alexander Kloosterman

# Making arrows

import os
    
from log import return_logger

logger = return_logger(__name__, False)
    
def SVG_operon(operon,all_domains,color_dict,domain_descr,html_path,settings, \
                ArrowHeight=20,HeadEdge=8,HeadLength=10,marginX=100,marginY=30,scaling=10.0,F=14):
                
    '''
    Create the main SVG document:
        - start from operon dictionary
        - find genes, start and stop positions, and strands
        - write the SVG files
    '''
    ALL_TEXT = ''
    
    # Line passing through all the genes
    l = line(marginX,marginY+ArrowHeight/2,(operon.extra_end-operon.extra_start)/scaling)
    ALL_TEXT += l + '\n'
    
    id = 0
    core = False
    boxes = ''
    lines = ''
    all_subs = ''
    descr_dict = {}
    prev_subtitles_printed = {}
    domain_firstline_prevgene = 100
    prev_x = 0
    max_y = 0
    genes = [gene for gene in operon.itergenes(sorted=True)]
    for gene in genes:
        # print(gene)
        id += 1
        start,stop = (gene.start-operon.extra_start + 1)/scaling,(gene.end-operon.extra_start)/scaling # +1 to adjust for pythonic 0-index 
        color = determine_gene_color(gene,all_domains,color_dict)
        extra_dict = {}
        domains = []
        for domain in gene.domains:
            try:
                descr = domain_descr[domain]
            except KeyError:
                logger.warning('No domain description found for domain %s - please update the domaindescr file' %domain)
                descr = 'Unknown'
            descr_dict[domain] = descr
            # Figure out which subclass this gene belongs to
            for key,value in all_domains.items():
                if domain in value:
                    if key not in extra_dict:
                        extra_dict[key] = domain
                    else:
                        extra_dict[key] += ' %s' %(domain)
                    break
            else:
                if 'other' not in extra_dict:
                    extra_dict['other'] = domain
                else:
                    extra_dict['other'] = '%s %s' %(extra_dict['other'],domain)
        for attr in ['COG','gene','name','locus_tag','protein_id','SVM_hit','intergenic','transl']:
            if hasattr(gene,attr):
                v = getattr(gene,attr)
                if type(v) == float:
                    v = round(v,3)
                extra_dict[attr] = str(v)
        # Get the DNA seq
            extra_dict['dna'] = gene.get_seq()
        # Set RRE data in gene
        if settings['run_rre']:
            if gene.RRE_hit == None:
                extra_dict['RRE'] = 'Not analyzed'
            else:            
                extra_dict['RRE'] = gene.RRE_hit
            # Add data for the best hit
            if gene.RRE_hit:
                extra_dict['rre_hit'] = gene.best_hit.split('|')[0]
                extra_dict['rre_loc'] = gene.RRE_data[gene.best_hit][6]
                extra_dict['rre_prob'] = gene.RRE_data[gene.best_hit][0]
            else:
                extra_dict['rre_hit'] = 'n/a'
                extra_dict['rre_loc'] = 'n/a'
                extra_dict['rre_prob'] = 'n/a'
        
        # Determine gene class
        if settings['run_rre'] and gene.RRE_hit and settings['pseudo'] and gene.pseudo:
            geneclass = 'RRE_pseudo'
        elif settings['run_rre'] and gene.RRE_hit:
            geneclass = 'RRE'
        elif settings['pseudo'] and gene.pseudo:
            geneclass = 'pseudo'
        else:
            geneclass = False
            
        extra_dict['coordinates'] = '%s -- %s' %(gene.start + 1,gene.end) # +1 to adjust for pythonic 0-index 
        
        # Add a dashed line between the core genes and the other genes
        if (not core and gene in operon.core_genes) or (core and gene not in operon.core_genes):
            core = not core
            if gene == genes[-1]:
                x_line = marginX + stop
            elif gene == genes[0]:
                x_line = marginX
            else:
                x_line = marginX + (start - prev_stop)/2 + prev_stop
                
            bound_length = HeadEdge + ArrowHeight
            y_line_start = marginY - 0.5*bound_length
            y_line_end =  marginY + 1.5*bound_length
            ALL_TEXT += line_exp(x_line,x_line,y_line_start,y_line_end,dashed=True,width=3) + '\n'
    
        polygon,box_outer = arrow(id,marginX+start,marginY,stop-start,ArrowHeight,gene.strand,HeadEdge,HeadLength,\
                                    color,svg_class=geneclass,title=gene.name,extra_dict=extra_dict,with_box=True)
        ALL_TEXT += polygon
        boxes += box_outer
        prev_stop = stop
        
        # Add COG scores and domains directly under the genes
        current_y = base_y = marginY + ArrowHeight * 2 + 16
        linenr = 1
        text_start = marginX+start+2
        # print(prev_subtitles_printed)
        if hasattr(gene,'COG'):
            text_cog,est_length = text(text_start,current_y,str(round(gene.COG,2)),extra_dict={'style':'display:none;'},text_class="gene_sub_cog")
            text_cog += '\n'
            prev_subtitles_printed[linenr] = text_start + est_length
#            lines += line_exp(marginX+start,marginX+start,current_y-20,current_y,extra_dict={'class':'gene_sub_cog','style':'display:none;'}) + '\n'
        else:
            text_cog = ''
        current_y += 15
        linenr += 1
        text_domains = ''
        # Make sure that there is enough space to print all the domains after each other
        while linenr in prev_subtitles_printed and prev_subtitles_printed[linenr] > text_start + 3 and len(gene.domains) > domain_firstline_prevgene - 2:
            linenr += 1
            current_y += 15
            # print('incrementing')
        # print('Starting at line %i' %linenr)
        domain_firstline_prevgene = linenr
        max_length_text = 0
        for domain in gene.domains:
            t,est_length = text(text_start,current_y,domain,extra_dict={'style':'display:none;'},text_class="gene_sub_domain")
            text_domains +=  t + '\n'
            current_y += 15
            if est_length > max_length_text:
                max_length_text = est_length
            for i in range(linenr):
                if i+1 not in prev_subtitles_printed or text_start + est_length > prev_subtitles_printed[i+1]:
                    prev_subtitles_printed[i+1] = text_start + est_length
            linenr += 1
        if len(gene.domains) > 0:
            lines += line_exp(marginX+start,marginX+start,base_y-20,(linenr-2)*15+base_y,\
                                    extra_dict={'class':'gene_sub_domain','style':'display:none;'}) + '\n'
            background = box(marginX+start+1,base_y + (domain_firstline_prevgene-2)*15 + 1, max_length_text * 0.8, len(gene.domains) * 15,cl='text_background')
        else:
            background = ''
        if current_y > max_y:
            max_y = current_y
        all_subs += text_cog + background + text_domains
        
    max_width = prev_stop + marginX*2
    ALL_TEXT += boxes
    ALL_TEXT += lines
    ALL_TEXT += all_subs
    
    svg_header = '''<svg width="%s" height="%s" xmlns="http://www.w3.org/2000/svg">
<style>
    @import url(%s)
</style>
''' %(max_width,max_y,os.path.join(os.path.relpath(html_path,operon.path),'svg_operon.css'))
    
    ALL_TEXT = svg_header + ALL_TEXT + '</svg>\n'
    
    return ALL_TEXT,descr_dict

def determine_gene_color(gene,all_domains,color_dict):
    if gene.SVM_hit == True:
        return(color_dict['precursor'])
    color = color_dict['other']
    for key in color_dict:
        if key in all_domains:
            domains_key = all_domains[key]
            if any([dom in domains_key for dom in gene.domains]):
                return(color_dict[key])
    return color
    
def write_operon_information_html(operon,group_collections,settings):
    if operon.COG:
        cog = '%s +- %s' %(round(operon.COG_avg,3),round(operon.COG_std,3))
    else:
        cog = 'N/A'
    info = '''<h2>Operon information</h2>
    <p>
    <b>Genome:</b> <a href="https://www.ncbi.nlm.nih.gov/assembly/{0}" target='_blank'>{1}</a><br>
    <b>Scaffold:</b> <a href="https://www.ncbi.nlm.nih.gov/nuccore/{2}" target='_blank'>{2}</a><br>
    <b>Coordinates:</b> {3} -- {4}<br>
    <b>Overall score:</b> {5}<br>
    <b>COG:</b> {6} <br>
    <b>Passed filter:</b> {9} <br>
    </p>
    <p><a href={7}><b>Genbank file</b></a><br>
    <b><a href={8}>Fasta file</b></a></p>\n'''.format(operon.genome.accession,operon.genome.descr,\
                                             operon.scaffold.realname,operon.start + 1,operon.end,\
                                             round(operon.score,2),cog,os.path.relpath(operon.genbank_path,start=operon.path),
                                             os.path.relpath(operon.fasta_path,start=operon.path), True) # +1 in start to adjust for pythonic index
     
    if settings['antismash_parse']:
        antismash_text = '''    <div class="info" id="antismash_overlap">
    <b>antiSMASH info</b><br>'''
        if operon.overlaps_antismash == {}:
            antismash_text += '''No overlap found.<br><br></div>\n'''
        else:
            for cluster in sorted(operon.overlaps_antismash.keys()):
                cluster_info = operon.overlaps_antismash[cluster]
                
                antismash_text +='''        <b>Cluster name: </b>{0}<br>
        <b>Cluster type: </b>{1}<br>
        <b>Cluster coordinates: </b>{2}<br>
        <b>Overlap type: </b>{3}<br>
        <hr>'''.format(cluster,'-'.join(cluster_info['cluster_type']),cluster_info['coords'],\
                             cluster_info['overlap_type'])
            # Remove the lastmost <hr>
            antismash_text = antismash_text[:-4]
            antismash_text +=  '<br></div>\n'
    else:
        antismash_text = ''
        
    if settings['run_rre']:
        genes_RRE = [gene for gene in operon.itergenes(RRE_hit=True)]
        if len(genes_RRE) == 0:
            RRE_text = '''    <div class="info" id="RRE_operon">No RRE hits found in the core operon.'''
        else:
            RRE_text = '''    <div class="info" id="RRE_operon">
    <table>
        <caption>RRE hits</caption>
        <tr>
            <th>Gene name</th>
            <th>Best hit</th>
            <th>Probability</th>
        </tr>\n'''
        
            for gene in genes_RRE:
                # Get the best RRE hit
                best_prob = gene.RRE_data[gene.best_hit][0]
                RRE_text += '''        <tr>
            <td>{0}</td>
            <td>{1}</td>
            <td>{2}</td>
        </tr>\n'''.format(gene.repr_name,gene.best_hit.split('|')[0],best_prob)
        
        RRE_text += '    </table></div>'
    else:
        RRE_text = ''
        
    # Add the links to operon groups
    operon_group_text ='''
    <table>
        <caption>Part of operon groups</caption>
    <tr>
        <th>Group type</th>
        <th>Group name</th>
    </tr>\n'''
    keywords = ['genome','jaccard','precursor','EC','precursor_jaccard']
    
    for keyword in keywords:
        gr_keyword = '%s_groups' %keyword
        if hasattr(operon,gr_keyword):
            groups = getattr(operon,gr_keyword)
        elif keyword == 'genome':
            groups = operon.genome.descr
        else:
            continue
        if groups == []:
            continue
        if type(groups) == list:
            groups = sorted(list(set(groups)))
            for group in groups:
                if group in group_collections[keyword]:
                    collection = group_collections[keyword][group]
                    operon_group_text += '''    <tr>
            <td>{0}</td>
            <td><a href="{2}">{1}</a></td>
        </tr>\n'''.format(keyword.capitalize(), group, os.path.relpath(collection.path,operon.path))
        elif type(groups) == str:
            if groups in group_collections[keyword]:
                collection = group_collections[keyword][groups]
                operon_group_text += '''    <tr>
            <td>{0}</td>
            <td><a href="{2}">{1}</a></td>
        </tr>\n'''.format(keyword.capitalize(), groups, os.path.relpath(collection.path,operon.path))
    operon_group_text += '</table>\n'
    info = info + antismash_text + RRE_text + operon_group_text 
    name = '<h2>{0}</h2>'.format(operon.name)
    
    return name,info                     

def write_domaindescr_html(domain_descr):
    text = '<div style="display:none;" id="descriptions"\n' 
    for dom in domain_descr:
        try:
            descr = domain_descr[dom]
        except KeyError:
            logger.warning('No domain description found for domain %s - please update the domaindescr file' %dom)
            descr = 'Unknown'
        text += '%s="%s"\n' %(dom,descr)
    text += '></div>'
    return text

def write_html_references(path):
    scr = '<script src=%s></script>' %os.path.join(path,'Operon.js')
    style = '<link rel="stylesheet" href="%s">' %os.path.join(path,'html_operon.css')
    return(scr,style)
    
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

def make_legend(color_dict,genetypes,html_path,settings,\
                ArrowHeight=20,HeadEdge=8,HeadLength=10,marginX=5,marginY=20,scaling=100.0,F1=20,F2=16):

    svg_header = '''<svg class=legend xmlns="http://www.w3.org/2000/svg">
    <style>
        @import url(%s/svg_operon.css)
    </style>
''' %(html_path)
    all_text = ''
    for genetype in color_dict:
        color = color_dict[genetype]
        svg = arrow('',marginX,marginY,50,ArrowHeight,'+',HeadEdge,HeadLength,color,title=genetypes[genetype],with_box=False)
        svg = svg_header + svg + '</svg>\n'
        
        descr = genetypes[genetype]
        html_text = '<div class=legend>' + svg + '<p>%s</p></div>\n' %(descr)
        all_text += html_text
    # For the pseudogenes
    if settings['pseudo']:
        svg = arrow('',marginX,marginY,50,ArrowHeight,'+',HeadEdge,HeadLength,color_dict['other'],title='Pseudogenes',svg_class='pseudo',with_box=False)
        svg = svg_header + svg + '</svg>\n'
        html_text = '<div class=legend>' + svg + '<p>Pseudogenes</p></div>\n'
        all_text += html_text
    
    # For the RRE marked genes
    if settings['run_rre']:
        svg = arrow('',marginX,marginY,50,ArrowHeight,'+',HeadEdge,HeadLength,color_dict['other'],svg_class='RRE',title='Genes with RRE hits',with_box=False)
        svg = svg_header + svg + '</svg>\n'
        html_text = '<div class=legend>' + svg + '<p>Genes with RRE hits</p></div>\n'
        all_text += html_text
        
    # For genes both pseudo and RRE
    if settings['pseudo'] and settings['run_rre']:
        svg = arrow('',marginX,marginY,50,ArrowHeight,'+',HeadEdge,HeadLength,color_dict['other'],svg_class='RRE_pseudo',title='Pseudogenes with RRE hits',with_box=False)
        svg = svg_header + svg + '</svg>\n'
        html_text = '<div class=legend>' + svg + '<p>Pseudogenes with RRE hits</p></div>\n'
        all_text += html_text
    
    return all_text

def make_RRE_box(settings):
    if settings['run_rre']:
        return '''    <div class="info" id="rre.info">
        <b>RRE element found:</b> <span id="gene.rre">yes/no</span><br>
        <u><b>RRE data highest hit</u></b><br>
        Hit name: <span id="gene.rre_hit">abc</span><br>
        Location: <span id="gene.rre_loc">1-100</span><br>
        hhpred probability: <span id="gene.rre_prob">50</span><br>
        </div>\n'''
    else:
        return ''

def arrow(id,X,Y,L,H,strand,h,l,color='#cccccc',svg_class=False,title='',extra_dict={},with_box=True):
    '''
    SVG code for arrow:
        - (X,Y) ... upper left (+) or right (-) corner of the arrow
        - L ... arrow length
        - H ... arrow height
        - strand
        - h ... arrow head edge width
        - l ... arrow head length
        - color
        - strand
    the edges are ABCDEFG starting from (X,Y)     
    '''
    class_text = ''
    if svg_class:
        class_text = 'class="%s" ' %svg_class
    outer_arrow = ''
    if strand == '+':
        
        A = [X,Y]
        B = [X+L-l,Y]
        C = [X+L-l,Y-h]
        D = [X+L,Y+H/2]
        E = [X+L-l,Y+H+h]
        F = [X+L-l,Y+H]
        G = [X,Y+H]
        
        if L < l:
            # squeeze arrow if length shorter than head length
            B = [X,Y]
            C = [X,Y-h]
            D = [X+L,Y+H/2]
            E = [X,Y+H+h]
            F = [X,Y+H]
    
    elif strand == '-':
        
        A = [X+L,Y]
        B = [X+l,Y]
        C = [X+l,Y-h]
        D = [X,Y+H/2]
        E = [X+l,Y+H+h]
        F = [X+l,Y+H]
        G = [X+L,Y+H]
        
        if L < l:
            # squeeze arrow if length shorter than head length
            B = [X+L,Y]
            C = [X+L,Y-h]
            D = [X,Y+H/2]
            E = [X+L,Y+H+h]
            F = [X+L,Y+H]
    func_text = ''
    if with_box:
        func_text = '''onclick="Toggle('%s')" ''' %id
    line = """ <polygon %s%sid="%s" 
                   points="%i,%i %i,%i %i,%i %i,%i %i,%i %i,%i %i,%i"
                   style="fill:%s;" """  % (class_text,func_text,id,A[0],A[1],B[0],B[1],C[0],C[1],D[0],D[1],E[0],E[1],F[0],F[1],G[0],G[1],color)
    for key, value in extra_dict.items():
        line += '%s="%s" ' %(key,value)
    line += '>\n'
    if title != '':
        line += '<title>%s</title>' %title
    line += '</polygon>\n'
    if with_box:
        box_pad = bp = 8
        extra_style = {'visibility':'hidden','color':'None'}
        box_outer = box(X-bp,Y-h-bp,L+2*bp,H+2*h+2*bp,id='%s.2'%id,cl='outer',extra_style=extra_style)
        return line,box_outer
    else:
        return line
    

def box(X,Y,width,height,id='bla',cl=False,extra_style={}):
    '''
    Draws a box (rectangle)
    X,Y = coordinates of topleft corner
    width,height
    '''
    class_text = ''
    if cl:
        class_text = 'class="%s" ' %cl
    if extra_style != {}:
        style_text = 'style="'
        for key,value in extra_style.items():
            style_text += '%s:%s;' %(key,value)
        style_text += '"'
    else:
        style_text = ''
    box = '''<rect %sx="%i" y="%i" width="%i" height="%i"
                id="%s" %s>
             </rect>''' %(class_text,X,Y,width,height,id,style_text)
    
    return box


def line(X,Y,L,width=2):
    '''
    Draw a line below genes
    '''
    
    line = """<line x1="%i" y1="%i" x2="%i" y2="%i"
              style="stroke-width:%s"/>""" % (X,Y,X+L,Y,width)
    return line
    
def line_exp(X1,X2,Y1,Y2,dashed=False,width=2,extra_dict={}):
    '''
    Draw any line
    '''
    
    extra_text = ''
    for key,value in extra_dict.items():
        extra_text += '%s="%s" ' %(key,value)
    
    line = '''<line %s x1="%i" x2="%i" y1="%i" y2="%i"\n''' %(extra_text,X1,X2,Y1,Y2)
    style = '''style="stroke-width:%s"/>''' %width
    dash = ''' stroke-dasharray="5, 5"\n''' 
    if dashed:
        return(line+dash+style)
    else:
        return(line+style)
        
def text(X,Y,s,extra_dict={},text_class=False,background=False):
    
    if text_class:
        classtext = 'class="%s" ' %text_class
    else:
        classtext = ''
    extra_text = ''
    for key,value in extra_dict.items():
        extra_text += '%s="%s" ' %(key,value)
        
    t = '''<text {0}{1} x="{3}" y="{4}">{2}</text>'''.format(classtext,extra_text,s,X,Y)
    if background:
        b = box(X,Y-15,len(s)*8,16,cl="text_background")
        t = b + '\n' + t
    return t,len(s)*10
    
def main(operons,all_domains,domain_descr,color_dict,genetypes,html_path,group_collections,settings):
    if not html_path.endswith(os.sep):
        html_path += os.sep
    html = read_template(os.path.join(html_path, 'Operon_template.html'))
    index_file = settings['paths']['index_file']
    svg_text_dict = {}
    for operon in operons:
        outfolder = operon.path
        back_to_index_path = '<a href="%s">Back to index</a>' %os.path.relpath(index_file,start=outfolder)
        script,style_header = write_html_references(os.path.relpath(html_path,start=operon.path))
        legend = make_legend(color_dict,genetypes,os.path.relpath(html_path,start=operon.path),settings)
        operon_header,operon_info = write_operon_information_html(operon,group_collections,settings)
        svg,domains_used = SVG_operon(operon,all_domains,color_dict,domain_descr,html_path,settings)
        domain_text = write_domaindescr_html(domains_used)
        RRE_box = make_RRE_box(settings)
        
        total_text = html[0] + style_header + html[1] + back_to_index_path + html[2] + operon_header + html[3] + \
        svg + html[4] + legend + html[5] + operon_info + html[6] + RRE_box + html[7] + domain_text + script + html[8]
        outfile = os.path.join(operon.path,'%s.html' %(operon.name)) 
        with open(outfile,'w') as handle:
            handle.write(total_text)
        
        svg_text_dict[operon.name] = svg
    return svg_text_dict
    
        
            




     



