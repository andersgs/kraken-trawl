from __future__ import print_function
'''
Expects a environmental variable called ASPERA_KEY
'''

import click
import os
import subprocess
import pandas as pd
import re
import shlex
import gzip
from Bio import SeqIO
import sys
import pkg_resources
import logging

# SOME CONSTANTS #############################################################

ncbi_host = 'ftp-private.ncbi.nlm.nih.gov'
ncbi_user = 'anonftp'

taxonomy_files = [
    '/pub/taxonomy/gi_taxid_nucl.dmp.gz',
    '/pub/taxonomy/taxdump.tar.gz',
]

assembly_files = {
    'general': '/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt',
    'refseq': '/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt',
    'archea': '/genomes/refseq/archaea/assembly_summary.txt',
    'bacteria': '/genomes/refseq/bacteria/assembly_summary.txt',
    'fungi': '/genomes/refseq/fungi/assembly_summary.txt'
}

plamids = {
    'plasmid1': '/refseq/release/plasmid/plasmid.1.genomic.gbff.gz',
    'plasmid2': '/refseq/release/plasmid/plasmid.2.genomic.gbff.gz'
}

# SOME BASIC CHECK FUNCTIONS  ################################################


def check_aspera_exists():
    '''
    Check if aspera and aspera key are available. If so, return a string
    with the command.
    '''
    cmd = 'ascp'
    try:
        subprocess.check_output([cmd, '--help'])
        logging.info("Found ascp...")
    except:
        logging.critical("Could not find ascp.")
        raise IOError
    try:
        aspera_key = os.environ['ASPERA_KEY']
        logging.info("Found the ASPERA key")
    except:
        logging.critical("Could not find the ASPERA_KEY Environment variable.")
        raise IOError
    cmd += ' -d --overwrite=diff -k2 -i '\
        + aspera_key\
        + ' -Q -k1 -T -l500m --host=ftp-private.ncbi.nlm.nih.gov '\
        '--user=anonftp --mode=recv --file-manifest=text '\
        '--file-manifest-path={} --file-pair-list={} {}'
    return(cmd)


def check_kraken_exists():
    '''
    Check if kraken exists, and if so return a string with the command.
    '''
    cmd = 'kraken'
    try:
        subprocess.check_output([cmd, '--help'], stderr=subprocess.PIPE)
        logging.info("Found kraken...")
    except:
        logging.critical("Please install kraken...")
        raise IOError
    return(cmd)

# SOME SETTING THE SCENE FUNCTIONS ###########################################


def parse_aspera_manifest_file(path, files, new_name):
    '''
    Check the Aspera manifest file for any updates
    '''
    manifest_file = [f for f in os.listdir(path)
                     if re.match("aspera-transfer", f)][0]
    was_updated = {}
    try:
        fi = open(manifest_file, 'r')
        file_log = [line.strip() for line in fi if re.search('^\"', line)]
        fi.close()
        for f in files:
            log_line = [line for line in file_log if re.search(f, line)][0]
            if re.search('0B', log_line):
                was_updated[f] = False
            else:
                was_updated[f] = True
    except:
        logging.critical("Something went wrong when parsing "
                         "the aspera transfer manifest file")
        raise IOError
    os.rename(manifest_file, new_name)
    return was_updated


def create_kraken_db_folder(path, db_name, cmd, ncbi_taxonomy_files):
    '''
    A function to create the folder structure for a new kraken database
    '''
    db_path = os.path.join(os.path.abspath(path), db_name)
    folders_to_make = [
        db_path,
        os.path.join(db_path, 'library'),
        os.path.join(db_path, 'library', 'added'),
        os.path.join(db_path, 'taxonomy')
    ]
    try:
        for folder in folders_to_make:
            os.mkdir(folder)
    except:
        logging.warning("Necessary folders for kraken DB {} already exist. "
                        "Proceeding will overwrite a few things."
                        .format(db_name))
    fo = open("aspera_taxonomy_src_dest.txt", 'w')
    files_to_unzip = []
    for fi in ncbi_taxonomy_files:
        fi_name = os.path.basename(fi)
        fi_dest = fi_name
        fo.write('{}\n{}\n'.format(fi, fi_dest))
        files_to_unzip.append(fi_dest)
    fo.close()
    cmd = cmd.format(path, "aspera_taxonomy_src_dest.txt",
                     os.path.join(db_name, 'taxonomy'))
    logging.info("Running {}".format(cmd))
    p = subprocess.Popen(shlex.split(cmd))
    p.communicate()
    was_updated = parse_aspera_manifest_file(path,
                                             files_to_unzip,
                                             "aspera_taxonomy_manifest.txt")
    for fi in files_to_unzip:
        path_to_file = os.path.join(db_name, 'taxonomy', fi)
        if re.search('tar', fi):
            if was_updated[fi]:
                path_to_folder = os.path.join(db_name, 'taxonomy')
                cmd = 'tar xzvf {} -C {}'.format(path_to_file, path_to_folder)
                logging.info("Running {}".format(cmd))
                p = subprocess.Popen(shlex.split(cmd))
                p.communicate()
            else:
                logging.info("File {} was not updated in this "
                             "transfer. Nothing to do.".format(fi))
        else:
            if was_updated[fi]:
                outfilename = fi.rstrip('.gz')
                path_to_outfile = os.path.join(db_name,
                                               'taxonomy',
                                               outfilename)
                cmd = 'gunzip -c {} > {}'.format(path_to_file, path_to_outfile)
                logging.info("Running {}".format(cmd))
                p = subprocess.Popen(cmd, shell=True)
                p.communicate()
            else:
                logging.info("File {} was not updated in this transfer. "
                             "Nothing to do.".format(fi))
    return db_path

# SOME FILE LOADING FUNCTIONS ################################################


def load_assembly_table(infile):
    '''
    This function at the moment assumes that the file has been independently
    downloaded, but could be made to download the most recent version from
    Genbank.

    ## TODO ##
    Add functionality to download file from Genbank
    '''
    logging.info("Loading assembly table")
    try:
        return pd.read_csv(infile, skiprows=1, sep='\t')
    except:
        logging.critical("Could load assembly table.")
        raise IOError


def load_taxon_list(taxon_list):
    '''
    Load a list of taxons to include in the database.

    One taxon per line.  It should be a Genus species pair per line, or
    just a Genus.

    These names will be searched in the `organism_name` column in the
    `assembly_summary_genbank.txt` file.

    ## NOTE ##
    This file must necessarily be provided by the user
    '''
    try:
        fi = open(taxon_list, 'r')
        taxon_list = [taxon.strip() for taxon in fi]
        fi.close()
        # If searching for genus name, we want to make sure it matches exactly
        # the genus name, so must add a trailling white space
        genus_list = [taxon + ' ' for taxon in taxon_list
                      if len(taxon.split()) == 1]
        species_list = [taxon for taxon in taxon_list
                        if len(taxon.split()) > 1]
        genus_list.extend(species_list)
        taxon_list = genus_list
    except:
        logging.critical("Could not load the taxon list!")
        raise IOError
    return taxon_list


# LET US GET SOME GENOMES ####################################################


# RUN THIS FUNCTION IF A TAXON_LIST IS PROVIDED ##############################
def filter_assemblies(taxon_list, assembl_tab, filter_opt='strict'):
    '''
    Takes two inputs: the taxon list and the assembly table.
    It will output a filterd assembly table.

    The current approach takes a cascading style of filtering. If the most
    strict approach does not find anything, then go down to the next level.

    Filtering options:
        * 'strict' --- take only reference genomes, if available. This is
                        indicated by a non-NA in the  `refseq_category` column
                        (current other options are 'reference genome' or
                        'representative genome')
        * 'moderate' --- take reference genomes, if available. If not, accept
                        'Chromosome' or 'Complete Genome' in the assembly_level
                        column

        * 'liberal' --- if all other filters fail, then accept 'Scaffold' and
                        'Contig' assemblies in the assembly_level column
        * 'all' --- gives all the genomes associated with that taxon
    '''

    org_name = assembl_tab.organism_name
    filtered = []
    missing = []
    # this addresses the issue of assemblies that do not have a specified
    # path in the file
    ix_missing_ftp_path = ~assembl_tab.ftp_path.isin(['na'])
    assembl_tab = assembl_tab[ix_missing_ftp_path]
    logging.info("Total assemblies missing an ftp path was {}."
                 .format(sum(ix_missing_ftp_path)))
    for taxon in taxon_list:
        ix = org_name.str.match(taxon)
        # make sure there is at least one hit
        n_hits = sum(ix)
        if(n_hits == 0):
            missing.append(taxon)
            logging.warning("Did *not* find any genomes for {}"
                            .format(taxon))
            logging.warning("If you believe it should be there, "
                            "please check the spelling!".format(taxon))
            logging.warning("Skipping to next taxon...".format(taxon))
            continue
        else:
            print("Found {} hits to {}".format(n_hits, taxon))
            taxon_tab = assembl_tab[ix]
        if filter_opt == 'all':
            logging.warning("You chose to export **all** {} assemblies for {}."
                            .format(n_hits, taxon))
            filtered.append(taxon_tab)
            continue
        else:
            # because this category is always included, it should always be run
            refseq_cat = taxon_tab.refseq_category
            ix = refseq_cat.isin(['reference genome', 'representative genome'])
            n_hits = sum(ix)
            if n_hits > 0:
                logging.info("Found {} hits for {} in refseq."
                             .format(n_hits, taxon))
                refseq_tab = taxon_tab[ix]
                filtered.append(refseq_tab)
                continue
            elif filter_opt in ['moderate', 'liberal']:
                logging.warning("Did *not* find any reference genomes for {}."
                                .format(taxon))
                logging.warning("As per your request, "
                                "I am moving in to unknown territory!")
                asm_level = taxon_tab.assembly_level
                ix = asm_level.isin(['Chromosome', 'Complete Genome'])
                n_hits = sum(ix)
                if n_hits > 0:
                    logging.info("Found {} Chromosome/Complete Genome "
                                 "hits for {}".format(n_hits, taxon))
                    complete_tab = taxon_tab[ix]
                    filtered.append(complete_tab)
                elif filter_opt == 'liberal':
                    ix = asm_level.isin(['Scaffold', 'Contig'])
                    n_hits = sum(ix)
                    if n_hits > 0:
                        logging.info("Found {} Scaffold/Contig hits for {}."
                                     .format(n_hits, taxon))
                        scaffolds_tab = taxon_tab[ix]
                        filtered.append(scaffolds_tab)
                    else:
                        logging.critical("We seem to have a problem! "
                                         "Even with the most liberal of "
                                         "searchers, I was still unable to "
                                         "find a genome for {}. "
                                         "IF YOU SEE THIS MESSAGE, "
                                         "PLEASE REPORT IT!".format(taxon))
                else:
                    logging.warning("I was unable to find any assemblies "
                                    "for {} that fitted your filtering "
                                    "options!".format(taxon))
                    logging.warning("Moving on to the next taxon...")
            else:
                logging.warning("I could *not* find any hits for {} in "
                                "refseq. If you are feeling adventurous, "
                                "perhaps loosen your search criteria."
                                .format(taxon))
    out_tab = pd.concat(filtered)
    logging.info("For {} taxa, I found {} assemblies."
                 .format(len(taxon_list), out_tab.shape[0]))
    if len(missing) > 0:
        fo = open("missing_taxon.txt", 'w')
        for m in missing:
            fo.write('{}\n'.format(m))
        fo.close()
    return(out_tab)

# RUN THIS FUNCTION ON A TABLE OF ASSEMBLIES #################################
def parse_assemb_rows(row):
    '''
    Parse an assembly_summary row
    '''
    logging.info("Prepping {} for download.".format(row.organism_name))
    assemb_dict = {}
    assemb_dict['organism'] = row.organism_name
    assemb_dict['path'] = re.findall(
                                    pattern='ftp://ftp.ncbi.nlm.nih.gov(/.*)$',
                                    string=row.ftp_path)[0]
    assemb_dict['asm_name'] = os.path.basename(assemb_dict['path'])
    assemb_dict['gbk_file'] = assemb_dict['asm_name'] + '_genomic.gbff.gz'
    assemb_dict['dest'] = os.path.join(assemb_dict['asm_name'],
                                       assemb_dict['gbk_file'])
    assemb_dict['source'] = '/'.join([assemb_dict['path'],
                                     assemb_dict['gbk_file']])
    assemb_dict['ref_status'] = row.refseq_category
    assemb_dict['asm_level'] = row.assembly_level
    return assemb_dict


def download_gbk(assemb_tab, cmd, outdir='.'):
    '''
    This should take a list of accessions and paths, and produce a
    source/destination pairs file.

    It is then possible to use the --file-pair-file to download all at once,
    and we can add the following flags to the ascp command:

    --overwrite=diff -k2

    the -k2 means files are compared with sparse checksums.
    '''

    fo = open("aspera_assemblies_src_dest.txt", 'w')
    assembs_dic_list = [parse_assemb_rows(row) for
                        row in assemb_tab.itertuples()]
    for assembly in assembs_dic_list:
        fo.write("{}\n{}\n".format(assembly['source'], assembly['dest']))
    fo.close()

    cmd = cmd.format(outdir, 'aspera_assemblies_src_dest.txt', outdir)
    logging.info("Running {}".format(cmd))
    p = subprocess.Popen(shlex.split(cmd))
    p.communicate()
    files_to_check = [a['gbk_file'] for a in assembs_dic_list]
    was_updated = parse_aspera_manifest_file(outdir,
                                             files_to_check,
                                             "aspera_assemblies_manifest.txt")
    logging.info("Finished downloading all genomes.")
    return assembs_dic_list


# LOADING GENOMES TO THE KRAKEN STAGGING AREA ################################
def kraken_add(db_name, fasta_file, clean=True):
    '''
    Add a FASTA sequence to the staging area in Kraken
    '''
    cmd = "kraken-build --add-to-library {} --db {}".format(fasta_file,
                                                            db_name)
    logging.info("Running {}".format(cmd))
    cmd = shlex.split(cmd)
    subprocess.check_output(cmd)
    if clean:
        os.remove(fasta_file)


def inject_adapters(db_name):
    '''
    Inject adapter sequences in to the Kraken DB
    '''
    adapter_fasta = pkg_resources.resource_filename(__name__,
                                                    os.path.join("data",
                                                                 "adapter.fasta"))
    logging.debug("Found adapter files here: {}.".format(adapter_fasta))
    kraken_add(db_name, adapter_fasta, False)
    logging.info("Added Illumina adapters to {}".format(db_name))


def add_isolates(dic_list, db_name):
    for assembly in dic_list:
        print("Adding {} to kraken stagging area.".format(assembly['organism']), file = sys.stderr)
        genbank_zip_file = assembly['dest']
        fi = gzip.open( genbank_zip_file, 'rt')
        seqs = list(SeqIO.parse( fi, 'genbank'))
        new_seqs = []
        for s in seqs:
            tmp = SeqIO.SeqRecord(s.seq)
            tmp.id = 'gi|{}'.format(s.annotations['gi'])
            tmp.description = s.description
            tmp.name = s.name
            new_seqs.append(tmp)
        fi.close()
        fa_file = os.path.join(os.getcwd(), os.path.basename(genbank_zip_file).strip('gbff.gz') + ".fa")
        tmpf = open(fa_file, 'wt')
        SeqIO.write(new_seqs, tmpf, 'fasta')
        tmpf.close()
        kraken_add(db_name, fa_file)
        # cmd = 'kraken-build --add-to-library {} --db {}'.format(fa_file, db_name)
        # print(cmd, file = sys.stderr)
        # cmd = shlex.split(cmd)
        # p = subprocess.check_output(cmd)
        # os.remove(fa_file)
    print("Added all {} assemblies to kraken stagging area. DB is ready to build".format(len(dic_list)), file = sys.stderr)

### ACTUALLY BUILD THE DATABASE ################################################

def kraken_build(db_name, threads = 16, kmer_len = 31, minz_len = 15, clean = True):
    cmd = 'kraken-build --build --db {} --threads {} --kmer-len {} --minimizer-len {}'.format(db_name, threads, kmer_len, minz_len)
    print(cmd, file = sys.stderr)
    p = subprocess.Popen( shlex.split(cmd))
    p.communicate()
    if clean:
        cmd = 'kraken-build --clean --db {}'.format(db_name)
        p = subprocess.Popen( shlex.split(cmd))
        p.communicate()
    print("Finished building the database.", file = sys.stderr)

### GENERATE A LOG OF SEQUENCES ADDED TO THE DATABASE ##########################

def generate_log(dic_list):
    print("Generating a log of added genomes.", file = sys.stderr)
    fo = open("log", 'w')
    header = ['Organism', 'Accession', 'RefSeq', 'Assembly Level', 'Source', 'Destination']
    fo.write('\t'.join(header) + '\n')
    for assembly in dic_list:
        row = '\t'.join( [assembly['organism'], assembly['asm_name'], assembly['ref_status'], assembly['asm_level'],  assembly['source'], assembly['dest'] ]) + '\n'
        fo.write(row)
    fo.close()

@click.command()
@click.option("--assemb_file", help = "assembly_summary file")
@click.option("--kraken_db", help = "name of kraken db")
@click.option("--create", "-c", help = "create a new database at path <give path>", default = None)
@click.option("--kraken_db_path", help = "path to kraken db", default = '.')
@click.option("--taxon_list", help = "give it a taxon list to filter the assembly_summary.", default = None)
@click.option("--include_human", help = "include the human reference genome", is_flag = True)
@click.option("--filter_opt", help = "Can be all, strict, moderate, or liberal.", default = 'strict')
@click.option("--outdir", help = "Where to place downloaded genomes.", default = ".")
@click.option("--no_log", help = 'Do NOT output a tab-delimited list of genomes added', is_flag = True)
@click.option("--do_not_inject_adapters", help = "Do **not** add Illumina adapter and primer sequences to DB", is_flag = True, default = False)
@click.option("--threads", help = 'How many threads to use when building the database', default = 16)
@click.option("--kmer_len", help = 'kmer length to use when building database', default = 31)
@click.option("--minz_len", help = 'minimizer length to use when building database', default = 15)
@click.option("--clean", help = 'clean db of unnecessary files after making the DB', default = True, is_flag = True)
def kraken_trawl(assemb_file, kraken_db, create, kraken_db_path, taxon_list, include_human, filter_opt, outdir, no_log, do_not_inject_adapters, threads, kmer_len, minz_len, clean ):
    aspera_cmd = check_aspera_exists()
    kraken_cmd = check_kraken_exists()
    if create != None:
        kraken_db_path = create_kraken_db_folder(create, kraken_db, aspera_cmd, taxonomy_files)
    if kraken_db_path == '.':
        kraken_db_path = os.getcwd()
    if not do_not_inject_adapters:
        inject_adapters(kraken_db)
    if outdir == '.':
        outdir = os.getcwd()
    assembs_raw = load_assembly_table(assemb_file)
    if taxon_list != None:
        taxon_list = load_taxon_list(taxon_list)
        assembs = filter_assemblies(taxon_list, assembs_raw, filter_opt = filter_opt)
    else:
        assembs = assembs_raw
    if include_human:
        ix_human = assembs_raw.species_taxid.isin([9606])
        ix_human_ref = assembs_raw.refseq_category.isin(['reference genome'])
        human = assembs_raw[ix_human & ix_human_ref]
        frames = [assembs, human]
        assembs = pd.concat(frames)
    n_assemblies = assembs.shape[0]
    print("Total assemblies found: {}".format(n_assemblies), file = sys.stderr)
    print("Total assemblies expected: {}".format(len(taxon_list)), file = sys.stderr)
    assemb_dic_list = download_gbk( assembs, aspera_cmd, outdir = outdir)
    add_isolates(assemb_dic_list, kraken_db)
    if not no_log:
        generate_log(assemb_dic_list)
    kraken_build(kraken_db, threads = threads, kmer_len = kmer_len, minz_len = minz_len, clean = clean )

if __name__ == '__main__':
    kraken_trawl()
