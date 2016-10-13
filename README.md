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

## Configuration

The current version relies on two environmental variables to work properly:

* `ASPERA_KEY`
* `EMAIL`

Either type the following before using `kraken-trawl`, or paste this line in to your
 `bash_profile`:

`export ASPERA_KEY=/path/to/asperakey/aspera_key.openssh`

If you wish to use the option of downloading from a species list, you must also
specify an EMAIL environmental variable for use in querying `NCBI` with `Entrez`.
Just type the following in your command line before running `kraken-trawl`, or
add this line to your `bash profile`:

`export EMAIL=my.valid.email@mailserver.com`

You may ask, why am I using environment for configuration? The reason is that this
is currently considered best practice in tool development. More here: [12factor app](https://12factor.net/config).

## Options


## Examples

## History

## Missing Features

## References

<a href='wood'>Wood</a> DE, Salzberg SL: Kraken: ultrafast metagenomic sequence classification using exact alignments. Genome Biology 2014, 15:R46 [URL](http://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-3-r46)
