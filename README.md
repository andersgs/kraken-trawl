# `kraken-trawl`: A tool to help you build custom `kraken` databases

## Use case

`kraken` is a tool for rapid identification of species from short read data ([Wood and Salzberg 2014](#wood)). It
uses a `kmer` database, with `kmers` identified by their `NCBI` Taxon ID. It comes with a
few databases preconfigured, but sometimes they are not enough, and they are also significantly
out of date.

While `kraken` provides for a limited API to build your own, custom, database, it leaves it
up to the user to find the sequences of interest, and rename them in an appropriate way so
that the `kraken-build` tool can interpret the sequences correctly.

What `kraken-trawl` aims to do is to automate as much as possible the creation of
custom `kraken` databases. The user only needs a list of species names of interest.
`kraken-trawl` will then download the genomes, if available, transform them to the
appropriate name, and place them in a folder structure appropriate for `kraken-build`.
The user can create a new DB from scratch, or add to an existing DB.

## Dependencies

### Kraken

`kraken-trawl` depends on `kraken`. Download and installation instructions can
be found [here](http://ccb.jhu.edu/software/kraken/MANUAL.htm).

### Aspera

The current version assumes you have `aspera` `ascp` installed and in the path.
It also assumes that you have an environmental variable called `ASPERA_KEY` that
points to the aspera key for downloading from `NCBI`. For a guide on obtaining
`ascp` and a key from `NCBI` go [here](https://www.ncbi.nlm.nih.gov/books/NBK242625/).

### Python

`kraken-trawl` also depends has the following `python` dependencies:

* `python 2.7`
* `BioPython`
* `pandas >=0.19`
* `click >=5`

## Installation

The easiest way of installing `kraken-trawl` is using `pip`:

`pip install git+https://github.com/andersgs/kraken-trawl.git`

Use the `--user` option to install locally:

`pip install --user git+https://github.com/andersgs/kraken-trawl.git`

Use the `--install-option` to install the script in a particular location:

`pip install --install-option="--install-scripts=$HOME/bin" --user git+hhttps://github.com/andersgs/kraken-trawl.git`

Once installed type the following:

`kraken-trawl --help`


## Configuration

The current version relies on an environmental variable to work properly:

* `ASPERA_KEY`

Either type the following before using `kraken-trawl`, or paste this line in to your
 `bash_profile`:

`export ASPERA_KEY=/path/to/asperakey/aspera_key.openssh`

You may ask, why am I using environment variables for configuration? The reason is that this
is currently considered best practice in tool development. More here: [12factor app](https://12factor.net/config).

## Options

```
Usage: kraken-trawl [OPTIONS]

Options:
  --assemb_file TEXT        assembly_summary file
  --kraken_db TEXT          name of kraken db                                                                                                         -c, --create TEXT         create a new database at path <give path>
  --kraken_db_path TEXT     path to kraken db                                                                                                         --taxon_list TEXT         give it a taxon list to filter the
                            assembly_summary.                                                                                                         --include_human           include the human reference genome
  --filter_opt TEXT         Can be all, strict, moderate, or liberal.                                                                                 --outdir TEXT             Where to place downloaded genomes.
  --no_log                  Do NOT output a tab-delimited list of genomes                                                                                                       added
  --do_not_inject_adapters  Do **not** add Illumina adapter and primer
                            sequences to DB
  --threads INTEGER         How many threads to use when building the database
  --kmer_len INTEGER        kmer length to use when building database
  --minz_len INTEGER        minimizer length to use when building database
  --clean                   clean db of unnecessary files after making the DB
  --help                    Show this message and exit.
```

## Output

1. A series of folders each containing a `*.genomic.gbff.gz` file. Because we
    are using `aspera` if the folder structure is maintained, only novel genomes
    will be downloaded subsequently.
2. an `aspera_src_dest.txt` file containing the `SRC` and `DEST` for each genome.
    This file is used to for the `aspera` download .
3. individual genome `*.fna` files in the `<kraken_db>/libraries/added` folder.
4. a `log` file, with a `tab-delimited` information on each downloaded file:
    * Organism --- the organism name, including any strain names
    * Accession --- the accession id
    * RefSeq --- the RefSeq category, if any
    * Assembly Level --- the assembly level (i.e., `Complete Genome`, `Chromosome`, `Scaffold` , or `Contig`)
    * Source  --- the location of the `genbank` file on the `NCBI` servers
    * Destination --- the local folder with the `genbank` file
5. a `missing_taxon.txt` --- a list of taxons for which no genomes were found

## Examples

## History

## Missing Features

## References

<a href='wood'>Wood</a> DE, Salzberg SL: Kraken: ultrafast metagenomic sequence classification using exact alignments. Genome Biology 2014, 15:R46 [URL](http://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-3-r46)
