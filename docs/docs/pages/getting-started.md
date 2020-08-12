# Getting Started

## Downloading

the source code for the pipeline library can be found on [github](https://github.com/gkm-software-dev/Sequence-Analysis-Pipelines/tree/experiment-jack). Simply clone the repository onto your local machine. 

## Installing...

### ...the library
Installation of the library itself is straightforward. Simply:

1. clone the repository as mentioned above to somewhere on your local machine where you intend to work, for example `~/Projects/Genomics/Sequence-Analysis-Pipelines`;

2. open the file `Sequence-Analysis-Pipelines/includes/locations.sh` and modify the line `pro_dir=/home/jack/Local/GeneMap/Sequence-Analysis-Pipelines` to point to where you have saved the library. 

3. format the `data/` directory by running the following script: `<pipeline library path>/pipelines/format-data-folder.sh`. You may need to make it executable first. If successfull, you should notice that many empty sub directories in `data/` have been created, for example `data/bam`.

### ...the library dependencies
The full pipeline library makes use of many external tools, each of which can be installed in the `tools/` folder inside the library's top directory, or somewhere else on your machine. Each tool needs to be discoverable by the library functions, so for each tool you install, the path of the tool executable must be added to the `includes/locations.sh` file. 

The Broad institute have published a tutorial for installing all the tools required to follow *gatk best practices*, which can be found [here](https://gatk.broadinstitute.org/hc/en-us/articles/360041320571--How-to-Install-all-software-packages-required-to-follow-the-GATK-Best-Practices).

The external dependencies for the whole pipeline library are:

1. **bwa mem** - Burrows-Wheeler aligner for aligning reads to a reference sequence. Source code and a manual can be found [here](http://bio-bwa.sourceforge.net/).

2. **gatk 4** - genome analysis toolkit for performing many operations on aligned and non aligned reads, including variant calling. Executable binary can be found on the Broad institute github pages [here](https://github.com/broadinstitute/gatk/releases). The gatk requires Java. 

3. **samtools** - a library for reading, writing, and manipulating sam and bam files. Source code and manual can be found [here](http://www.htslib.org/).

4. **bcftools** - a library for reading, writing, and manipulating vcf files. Source and documentation is found [here](http://samtools.github.io/bcftools/bcftools.html).

5. **fastqc** - a tool for summarising the quality of input sequence reads (fastq) files. Source and documentation is [here](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).

6. **trimmomatic** - a tool for trimming raw fastq reads and removing reads smaller than some minimum length. Source and documentation is [here](http://www.usadellab.org/cms/?page=trimmomatic).

## Running a pipeline

Once you have set up all the required tools, you can test your library installation. In `<SAP path>/test-suites/` you will find some scripts that test different aspects of the pipeline library. Here we will run one called `simple-cohort.sh` that simulates a cohort of read groups from a genomic sequence mutated from an input reference, `<SAP path>/data/references/lambda-virus.fa`, aligns the read groups individually to the reference, calls the variants, and then compares the final vcf file of the cohort to the true vcf file from the simulation of the read groups. 

First, navigate to the pipelines folder with:

`cd <SAP path>/pipelines/` .

In this directory, open the file: `simulate-cohort.basic.json` and navigate to the line `"threads": 8,`. Change this value to the number of processors your machine has. You can see this visually on a linux machine by running the command `htop` from the command line, and then pressing `q` to exit (the number of bars at the top of the terminal display will be the number of local processors). 

Next, open the file `gatkall.basic.json`, naviate to `"threads": 8,` and again change this value to the number of processors your machine has. 

Finally, navigate to the test-suite folder with 

`cd <SAP path>/test-suites/`

and run:

`./simple-cohort.sh simple-cohort.basic.json` .

You may have to make this script executable first: simply type `chmod +x ./simple-cohort.sh`. The run may take from one to several minutes, depending on your computer hardware and parameter choices. 

If successfull, you will see, at the end of the terminal output stream, something like:

```
computing jaccard index...
  computing jaccard index......done
...done.
  Jaccard value: 0.314814814815
  pipeline runtime: 88.687579383 seconds
```

The Jaccard value is a [measure of how simular two sets are](https://en.wikipedia.org/wiki/Jaccard_index). Here we are measuring the difference between the truth vcf file, and the estimated one, where each vcf file is a set of variants. For more details, naviagate to the `<SAP path>/data/vcf/c_1.isec` directory. Saved here is the full comparison between truth and estimated vcf files for your simulation. 

*If the test-suite fails with errors*, the details of the errors are not printed to the terminal screen. To find the full log out put, check the directory `<SAP path>/data/logs` for information. The job you have just run has a random ID code, which was printed to the terminal. Use this to find the logs associated with the failed run to diagnose the problem. 