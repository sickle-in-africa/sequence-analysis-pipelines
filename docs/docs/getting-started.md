# Getting Started

## Downloading

the source code for the pipeline library can be found on [github](https://github.com/gkm-software-dev/Sequence-Analysis-Pipelines/tree/experiment-jack). Simply clone the repository onto your local machine. 

## Installing...

### ...the library
Installation of the library itself is straightforward. Simply:

1. clone the repository as mentioned above to somewhere on your local machine where you intend to work, for example `~/Projects/Genomics/Sequence-Analysis-Pipelines`;

2. open the file `Sequence-Analysis-Pipelines/includes/locations.sh` and modify the line `pro_dir=/home/jack/Local/GeneMap/Sequence-Analysis-Pipelines` to point to where you have saved the library. 

### ...the library dependencies
The full pipeline library makes use of many external tools, each of which can be installed in the `tools/` folder inside the library's top directory, or somewhere else on your machine. Each tool needs to be discoverable by the library functions, so for each tool you install, the path of the tool executable must be added to the `includes/locations.sh` file. 

The Broad institute have published a tutorial for installing all the tools required to follow *gatk best practices*, which can be found [here](https://gatk.broadinstitute.org/hc/en-us/articles/360041320571--How-to-Install-all-software-packages-required-to-follow-the-GATK-Best-Practices).

The external dependencies for the whole pipeline library are:

1. **bwa mem** - Burrows-Wheeler aligner for aligning reads to a reference sequence. Source code and a manual can be found [here](http://bio-bwa.sourceforge.net/).

2. **gatk 4** - genome analysis toolkit for performing many operations on aligned and non aligned reads, including variant calling. Executable binary can be found on the Broad institute github pages [here](https://github.com/broadinstitute/gatk/releases). The gatk requires Java. 

3. **samtools** - a library for reading, writing, and manipulating fastq, sam, bam, and vcf files. Source code and manual can be found [here](http://www.htslib.org/).

## Running a pipeline