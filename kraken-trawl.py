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

def load_assembly_table(infile):
    print("Loading assembly table", file = sys.stderr)
    try:
        return( pd.read_csv( infile, skiprows = 1, sep = '\t' ) )
    except:
        print("Could load assembly table.", file = sys.stderr)
        raise IOError

def download_gbk(accession, path, cmd, kraken_db):
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

def get_tax_id(species_list):
    try:
        email = os.environ['EMAIL']
    except:
        print("To use species_list, you need to set EMAIL environment variable to a valid email address. For instance: export EMAIL=john_doe@mailserver.com.", file = sys.stderr)
        raise IOError
    Entrez.email = email
    fi = open(species_list, 'r')
    species_taxid = []
    for species in fi:
        species = species.strip() + '[Name Tokens]'
        print("Searching for: {}.".format(species), file = sys.stderr)
        search = Entrez.esearch(term = species, db = "taxonomy", retmode = "xml")
        rec = Entrez.read(search)
        res_ids = rec['IdList']
        if len(res_ids) > 0:
            for tid in res_ids:
                species_taxid.append(int(tid))
        else:
            print("Did not find a taxon id for {}.".format(species), file = sys.stderr)
            print(rec, file = sys.stderr)
            #raise IOError
        print("Done.", file = sys.stderr)
    fi.close()
    return species_taxid

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

@click.command()
@click.option("--input", help = "assembly_summary file")
@click.option("--kraken_db", help = "name of kraken db")
@click.option("--species_list", help = "give it a species list to filter the assembly_summary.", default = None)
@click.option("--include_human", help = "give it a species list to filter the assembly_summary.", is_flag = True)
def kraken-trawl(input, kraken_db, species_list, include_human):
    aspera_cmd = check_aspera_exists()
    kraken_cmd = check_kraken_exists()
    assembs_raw = load_assembly_table(input)
    if species_list != None:
        taxid = get_tax_id(species_list)
        #import pdb; pdb.set_trace()
        ix = assembs_raw.species_taxid.isin(taxid)
        assembs = assembs_raw.loc[ix]
        ix_ref_cat = assembs.refseq_category.isin(['reference genome', 'representative genome'])
        ix_asm_level = assembs.assembly_level.isin(['Chromosome', 'Complete Genome'])
        assembs_ref = assembs[ix_ref_cat | ix_asm_level]
        filtered = check_missing_taxid(assembs_ref, assembs, taxid)
        if len(filtered) > 0:
            ix = assembs_raw.species_taxid.isin(filtered)
            assembs_filtered = assembs_raw.loc[ix]
            frames = [assembs, assembs_filtered]
            assembs = pd.concat(frames)
        else:
            assembs = assembs_ref
        if include_human:
            ix_human = assembs_raw.species_taxid.isin([9606])
            ix_human_ref = assembs_raw.refseq_category.isin(['reference genome'])
            human = assembs_raw[ix_human & ix_human_ref]
            frames = [assembs, human]
            assembs = pd.concat(frames)
        n_assemblies = assembs.shape[0]
        print("Total assemblies found: {}".format(n_assemblies), file = sys.stderr)
        print("Total assemblies expected: {}".format(len(taxid)), file = sys.stderr)
        print("Are the numbers of genomes found equal to number of taxa: {}".format(n_assemblies == len(taxid)), file = sys.stderr)
        total_taxid = sum(assembs['species_taxid'].isin(taxid))
        print("Are the numbers of taxid in genomes found equal to number of taxa searched: {}".format(n_assemblies == len(taxid)), file = sys.stderr)
    else:
        assembs = assembs_raw
    for r in assembs.itertuples():
        download_gbk(r.asm_name, r.ftp_path, aspera_cmd, kraken_db)

if __name__ == '__main__':
    kraken-trawl()
