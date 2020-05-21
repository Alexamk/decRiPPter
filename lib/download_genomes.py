from ftplib import FTP
from lib import find_file, clean
from log import return_logger

import gzip
import tarfile
import os
import traceback

logger = return_logger(__name__, False)

def handleDownload(f,block):
    f.write(block)

def ftp_connect():
    logger.debug('Connecting to ncbi through ftp')
    ftp = FTP('ftp.ncbi.nlm.nih.gov')
    ftp.login()
    return ftp
    
def ftp_disconnect(ftp):
    logger.debug('Disconnecting from ftp')
    try:
        ftp.quit()
    except:
        ftp.close()

def ftp_reconnect(ftp):
    ftp_disconnect(ftp)
    ftp = ftp_connect()
    return ftp

def gunzip(f):
    logger.debug('Gunzip on ' + f)
    clean_file = f.rpartition('.gz')[0]
    with open(clean_file,'w') as handle:
        with gzip.open(f) as inf:
            handle.write(inf.read())
    clean(f)
    
def unpack(f,move_files_up=False):
    folder = os.path.dirname(f)
    logger.debug('Unpacking %s to %s' %(f,folder))
    if f.endswith('.tar.gz'):
        tar = tarfile.open(f,'r:gz')
        tar.extractall(folder)
        file_clean = f[:-7]
    elif f.endswith('.tar'):
        tar = tarfile.open(f)
        tar.extractall(folder)
        file_clean = f[:-4]
    elif f.endswith('.gz'):
        file_clean = f.rpartition('.gz')[0]
        with open(file_clean,'w') as handle:
            with gzip.open(f) as inf:
                handle.write(inf.read())
    if os.path.isfile(f):
        clean(f)
    if move_files_up and os.path.isdir(file_clean):
        logger.debug(file_clean)

def download_genomes_ftp(ftp,prodigal,accessions,acc2data,dest_path,downloaded=[],ext=False,forbidden_words=False):
    logger.debug('Starting function "download_genomes_ftp"')
    for acc in accessions:
        logger.debug('Accession: %s' %acc)
        try:
            if acc not in downloaded:
                if acc not in acc2data:
                    # No data - no download
                    downloaded.append(acc)
                    continue
                folder = acc2data[acc][3][20:]
                logger.debug('cwd to ftp folder: %s' %folder)
                ftp.cwd(folder)
                lines = []
                ftp.retrlines('LIST',lines.append)
                logger.debug('lines found: ')
                for i in lines:
                    logger.debug(i)
                if ext:
                    # Force the script to download different extensions. Not used in the program right now, might be later
                    rel_files = [f for tabs in [line.split(' ') for line in lines] for f in tabs if ext in f]
                    if forbidden_words:
                        logger.debug('Filtering files containing the following words: ' + str(forbidden_words))
                        rel_files = [f for f in rel_files if not any([word in f for word in forbidden_words])]
                else:
                    if prodigal == 'never' or prodigal == 'auto':
                        rel_files = [f for tabs in [line.split(' ') for line in lines] for f in tabs if '.gbff.gz' in f and '_from_' not in f and 'wgsmaster' not in f]
                        if rel_files == [] and prodigal == 'auto':
                            rel_files = [f for tabs in [line.split(' ') for line in lines] for f in tabs if '.fna.gz' in f and '_from_' not in f]
                    else:
                        rel_files = [f for tabs in [line.split(' ') for line in lines] for f in tabs if '.fna.gz' in f and '_from_' not in f]
                    if len(rel_files) > 1:
                        logger.warning('Multiple files found for accession: %s\nPlease inspect ftp: %s' %(acc,folder))
                        logger.debug(rel_files)
                        continue
                logger.debug('%i files found: ' %len(rel_files))
                for i in rel_files:
                    logger.debug(i)
                if len(rel_files) == 0:
                    logger.warning('Warning! No files found for genome %s. Continuing...' %acc)
                    continue
                f = rel_files[0] # Should only be one file
                genome_folder = os.path.join(dest_path, acc)
                if not os.path.isdir(genome_folder):
                    os.mkdir(genome_folder)
                with open(os.path.join(genome_folder, f),"wb") as download_file:
                    logger.info('Downloading file: %s' %f)
                    ftp.retrbinary("RETR " + f, lambda x: handleDownload(download_file,x))
                unpack(os.path.join(genome_folder,f))
                downloaded.append(acc)
        except:
            error = traceback.format_exc()
            logger.debug(error)
            logger.info('Downloading of %s failed. Reconnecting...' %acc)
            return(downloaded)
            
    return(downloaded)

def find_downloaded_files(genome_path,prodigal,ext=False):
    logger.debug('Starting function "find_downloaded_files"')
    already_dl = []
    for folder in os.listdir(genome_path):
        if ext:
            exts = [ext] + [ext.rpartition('.gz')[0]]
        else:
            if prodigal == 'never':
                exts = ['.gbk','.gbff','.gbk.gz','.gbff.gz']
            elif prodigal == 'always':
                exts = ['.fna','.fna.gz']
            else:
                exts = ['.gbk','.gbk.gz','.gbff','.gbff.gz','.fna','.fna.gz']
        
        f = find_file(os.path.join(genome_path,folder),exts,[])
        if f:
            already_dl.append(folder)
    return already_dl
    
def main(genome_path, prodigal, accessions, acc2data, ext=False, forbidden=False):
    #TODO also download md5 and use it to verify downloaded file integrity
    
    logger.debug('Starting function "main"')
    downloaded = find_downloaded_files(genome_path,prodigal,ext)
    while any([i not in downloaded for i in accessions]):
        ftp = ftp_connect()
        downloaded = download_genomes_ftp(ftp,prodigal,accessions,acc2data,genome_path,downloaded,ext=ext,forbidden_words=forbidden)
        ftp_disconnect(ftp)
    check_gzip(genome_path)
    
def check_gzip(genome_path):
    # One final check to see all downloaded files are unpacked
    logger.debug('Starting function "check_gzip"')
    for folder in os.listdir(genome_path):
        f = find_file(os.path.join(genome_path,folder),['.gz'],[])
        if f:
            unpack(f)
                                
def download_assembly_file(path,type='refseq'):
    # Alternative type: genbank
    if type == 'refseq':
        logger.debug('Using refseq summary file')
        ftp_path = 'genomes/refseq/assembly_summary_refseq.txt'
        download_file = os.path.join(path,'assembly_summary_refseq.txt')
    elif type == 'genbank':
        ftp_path = 'genomes/genbank/assembly_summary_genbank.txt'
        download_file = os.path.join(path,'assembly_summary_genbank.txt')
    ftp = ftp_connect()
    logger.debug('Downloading assembly file to %s' %(download_file))
    with open(download_file,'w') as handle:
        ftp.retrbinary("RETR " + ftp_path, lambda x: handleDownload(handle,x))
    ftp_disconnect(ftp)
    
def download_taxonomy_files(path):
    ftp_path = 'pub/taxonomy/'
    req_file = 'taxdump.tar.gz'
    download_file = os.path.join(path, req_file)
    ftp = ftp_connect()
    logger.debug('Downloading taxonomy files to %s' %(download_file))
    with open(download_file,'w') as handle:
        ftp.retrbinary("RETR " + ftp_path+req_file, lambda x: handleDownload(handle,x))
    unpack(download_file,move_files_up=False)

def download_genomes(accessions,acc2data,dest_path,already_dl=[],gunzip=True):
    for acc in accessions:
        if acc not in already_dl:
            ftp = acc2data[acc][3]
            commands = ['rsync', '--copy-links', '--recursive', \
             '--times', '--verbose', 'rsync://' + ftp, dest_path + acc + '/']
            logger.info("Downloading %s..." %acc)            
            call(commands)
        else:
            logger.debug('Already downloaded: %s' %acc)
            continue
        folder_extra = os.listdir(os.path.join(dest_path,acc))[0]
        final_path = os.path.join(dest_path,acc,folder_extra)            
        for f in os.listdir(final_path):
            os.rename(os.path.join(final_path,f),os.path.join(dest_path,acc,f))
        os.rmdir(final_path.rstrip('/'))
        if gunzip:
            logger.info("Unpacking files...")
            filenames = [i for i in os.listdir(dest_path+acc+'/') if i[-3:] == '.gz']
            filenames_full = [os.path.join(dest_path,acc,i) for i in filenames]
            commands = ['gunzip', '-f'] + filenames_full
            logger.debug(' '.join(commands))
            call(commands)


