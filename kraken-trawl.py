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
from Bio import Entrez
import tempfile
import sys

### SOME BASIC CHECK FUNCTIONS  ################################################

def check_aspera_exists():
    cmd = 'ascp'
    try:
        p = subprocess.check_output( [cmd, '--help'] )
        print("Found ascp", file = sys.stderr)
    except:
        print("Could not find ascp.", file = sys.stderr)
        raise IOError
    try:
        aspera_key = os.environ['ASPERA_KEY']
        print("Found the ASPERA key", file = sys.stderr)
    except:
        print("Could not find the ASPERA_KEY Environment variable.", file = sys.stderr)
        raise IOError
    cmd = cmd + ' -i ' + aspera_key + ' -Q -k1 -T -l500m anonftp@ftp-private.ncbi.nlm.nih.gov:'
    return(cmd)

def check_kraken_exists():
    cmd = 'kraken'
    try:
        p = subprocess.check_output( [cmd, '--help'], stderr = subprocess.PIPE )
        print("Found kraken", file = sys.stderr)
    except:
        print("Please install kraken", file = sys.stderr)
        raise IOError
    return(cmd)

### SOME FILE LOADING FUNCTIONS ################################################

def load_assembly_table(infile):
    '''
    This function at the moment assumes that the file has been independently
    downloaded, but could be made to download the most recent version from
    Genbank.

    ## TODO ##
    Add functionality to download file from Genbank
    '''
    print("Loading assembly table", file = sys.stderr)
    try:
        return( pd.read_csv( infile, skiprows = 1, sep = '\t' ) )
    except:
        print("Could load assembly table.", file = sys.stderr)
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
        genus_list = [taxon + ' ' for taxon in taxon_list if len(taxon.split()) == 1]
        species_list = [taxon for taxon in taxon_list if len(taxon.split()) > 1]
        genus_list.extend(species_list)
        taxon_list = genus_list
    except:
        print("Could not load the taxon list!", file = sys.stderr)
        raise IOError
    return taxon_list

### LET US GET SOME GENOMES ####################################################


### RUN THIS FUNCTION IF A TAXON_LIST IS PROVIDED ##############################
def filter_assemblies(taxon_list, assembl_tab, filter_opt ='strict'):
    '''
    Takes two inputs: the taxon list and the assembly table.
    It will output a filterd assembly table.

    Filtering options:
        * 'strict' --- take only reference genomes, if available. This is
                        indicated by a non-NA in the  `refseq_category` column (
                        current other options are 'reference genome' or
                        'representative genome'
                        )
        * 'moderate' --- take reference genomes, if available. If not, accept
                        'Chromosome' or 'Complete Genome' in the assembly_level column

        * 'liberal' --- if all other filters fail, then accept 'Scaffold' and
                        'Contig' assemblies in the assembly_level column
        * 'all' --- gives all the genomes associated with that taxon
    '''

    org_name = assembl_tab.organism_name
    filtered = []
    for taxon in taxon_list:
        ix = org_name.str.match(taxon)
        # make sure there is at least one hit
        n_hits = sum(ix)
        if(n_hits == 0):
            print("Did **not** find any genomes for {}".format(taxon), file=sys.stderr)
            print("If you believe it should be there, please check the spelling!".format(taxon), file=sys.stderr)
            print("Skipping to next taxon...".format(taxon), file=sys.stderr)
            continue
        else:
            print("Found {} hits to {}".format(n_hits, taxon))
            taxon_tab = assembl_tab[ix]
        if filter_opt == 'all':
            print("You chose to export **all** {} assemblies for {}.".format(n_hits, taxon), file = sys.stderr)
            filtered.append(taxon_tab)
            continue
        else:
            # because this category is always included, it should always be run
            refseq_cat = taxon_tab.refseq_category
            ix = refseq_cat.isin(['reference genome', 'representative genome'])
            n_hits = sum(ix)
            if n_hits > 0:
                print("Found {} hits for {} in refseq.".format(n_hits, taxon))
                refseq_tab = taxon_tab[ix]
                filtered.append(refseq_tab)
                continue
            elif filter_opt in ['moderate', 'liberal']:
                print("Did **not** find any reference genomes for {}.".format(taxon), file = sys.stderr)
                print("As per your request, I am moving in to unknown territory!")
                asm_level= taxon_tab.assembly_level
                ix = asm_level.isin(['Chromosome', 'Complete Genome'])
                n_hits = sum(ix)
                if n_hits > 0:
                    print("Found {} Chromosome/Complete Genome hits for {}".format(n_hits, taxon), file = sys.stderr)
                    complete_tab = taxon_tab[ix]
                    filtered.append(complete_tab)
                elif filter_opt == 'liberal':
                    ix = asm_level.isin(['Scaffold','Contig'])
                    n_hits = sum(ix)
                    if n_hits > 0:
                        print("Found {} Scaffold/Contig hits for {}.".format(n_hits, taxon, file = sys.stderr))
                        scaffolds_tab = taxon_tab[ix]
                        filtered.append(scaffolds_tab)
                    else:
                        print("We seem to have a problem! Even with the most liberal of searchers, I was still unable to find a genome for {}. IF YOU SEE THIS MESSAGE, PLEASE REPORT IT!".format(taxon), file = sys.stderr)
                else:
                    print("I was unable to find any assemblies for {} that fitted your filtering options!".format(taxon), file = sys.stderr)
                    print("Moving on to the next taxon...", file = sys.stderr)
            else:
                print("I could **not** find any hits for {} in refseq. If you are feeling adventurous, perhaps loosen your search criteria.".format(taxon), file = sys.stderr)
    out_tab = pd.concat(filtered)
    print("For {} taxa, I found {} assemblies.".format(len(taxon_list), out_tab.shape[0]), file = sys.stderr)
    return(out_tab)

### RUN THIS FUNCTION ON A TABLE OF ASSEMBLIES #################################
def download_gbk(assemb_tab, cmd, kraken_db):
    '''
    This should take a list of accessions and paths, and produce a
    source/destination pairs file.

    It is then possible to use the --file-pair-file to download all at once,
    and we can add the following flags to the ascp command:

    --overwrite=diff -k2

    the -k2 means files are compared with sparse checksums.
    '''
    path = re.findall(pattern='ftp://ftp.ncbi.nlm.nih.gov(/.*)$', string=path)[0]
    asm_name = os.path.basename(path)
    gbk_file = asm_name + '_genomic.gbff.gz'
    outpath = os.path.join( asm_name, gbk_file)
    if os.path.exists(outpath):
        print("Already downloaded {}".format(asm_name), file = sys.stderr)
    else:
        print("Downloading {}".format(asm_name), file = sys.stderr)
        if (not os.path.exists( asm_name )):
            os.mkdir( asm_name )
        cmd = cmd + path + '/' + gbk_file + ' ' + outpath
        print(cmd, file = sys.stderr)
        p = subprocess.check_output( shlex.split( cmd ) )
        add_isolate(outpath, kraken_db)
    print("Finished {}".format(asm_name), file = sys.stderr)

def check_missing_taxid(assemb_table_new, assemb_table_old, taxid):
    miss_id = []
    for t in taxid:
        if t not in assemb_table_new.species_taxid.tolist():
            print("{}".format(t))
            miss_id.append(t)
    #import pdb; pdb.set_trace()
    filtered_ix = assemb_table_old.species_taxid.isin(miss_id)
    filtered = assemb_table_old.loc[filtered_ix]
    filtered_list = filtered.species_taxid.unique().tolist()
    n_filtered = len(filtered_list)
    n_miss = len(miss_id) - n_filtered
    missing_list = [m for m in miss_id if m not in filtered_list]
    print("Missing IDs: {}".format(missing_list), file=sys.stderr)
    print("Filtered IDs: {}".format(filtered_list), file=sys.stderr)
    return(filtered_list)

### LOADING GENOMES TO THE KRAKEN STAGGING AREA ################################

def add_isolate(genbank_zip_file, db_name):
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
    cmd = 'kraken-build --add-to-library {} --db {}'.format(fa_file, db_name)
    print(cmd, file = sys.stderr)
    cmd = shlex.split(cmd)
    p = subprocess.check_output(cmd)
    os.remove(fa_file)
    return(0)

@click.command()
@click.option("--assemb_file", help = "assembly_summary file")
@click.option("--kraken_db", help = "name of kraken db")
@click.option("--taxon_list", help = "give it a taxon list to filter the assembly_summary.", default = None)
@click.option("--include_human", help = "give it a species list to filter the assembly_summary.", is_flag = True)
@click.option("--filter_opt", help = "Can be all, strict, moderate, or liberal.", default = 'strict')
def kraken_trawl(assemb_file, kraken_db, taxon_list, include_human, filter_opt):
    aspera_cmd = check_aspera_exists()
    kraken_cmd = check_kraken_exists()
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
    # print("Are the numbers of genomes found equal to number of taxa: {}".format(n_assemblies == len(taxid)), file = sys.stderr)
    # total_taxid = sum(assembs['species_taxid'].isin(taxid))
    # print("Are the numbers of taxid in genomes found equal to number of taxa searched: {}".format(n_assemblies == len(taxid)), file = sys.stderr)
    # for r in assembs.itertuples():
    #     download_gbk(r.asm_name, r.ftp_path, aspera_cmd, kraken_db)

if __name__ == '__main__':
    kraken_trawl()
