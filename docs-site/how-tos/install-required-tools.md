---
layout: default
title: installing required tools
parent: how-to's
---

How do I install the SAP library's required tools?
==================================================

Installation instructions for these tools can all be found on their project source web pages, but crucially for this library to find the tools, it is required that you install each tool in the `<sap-path>/tools` library, and that you add the relative path to the file `tool-list.json`. Specific install instructions for each tool are given below. 

No pipeline requires every tool, and you may only need to install the ones that will be used for your project. See the pipeline library reference manual for details on which tools are required by which pipeline. 

# bcftools

BCFtools is a set of utilities that manipulate variant calls in the Variant Call Format (VCF) and its binary counterpart BCF. A stable version of the tool can be downloaded from [here](http://www.htslib.org/download/).

1. Download the most recent version of `bcftools-X.X.X`, where `X.X.X` will be some version number, for example `1.10.2`, and extract the contents of the tar archive to the `<sap-path>/tools` directory with:
	```
	cd <your-downoads-path>
	tar -C <sap-path>/tools -xvf bcftools-X.X.X.tar.bz2
	```

2. Install the package to `<sap-path>/tools/bcftools-X-X-X` using the `--prefix` option in the configure step:
	```
	cd <sap-path>/tools/bcftools-X.X.X
	./configure --prefix=$(pwd)
	make
	make install
	```
	There should now be build files inside the `bcftools-X.X.X` folder. 

3. Test the tool installed correctly by typing:```
	./bcftools --help
	``` and you should see the help manual, which starts with something like this:```
	Program: bcftools (Tools for variant calling and manipulating VCFs and BCFs)
	Version: 1.10.2 (using htslib 1.10.2)
	...
	```

4. Tell the SAP library where to find the bcftools executable, `bcftools`. For this, you simply add the path of this file, *relative to the SAP library tools directory* to the file `<sap-path>/tools/tool-list.json` as in the following example:
```
"bcftools": "bcftools-1.10.2/bcftools",
```
Now, when the setup script is run, this tool will be accessible to the library. Remember that the quotes and comma at the end of the line in the above example are needed to ensure that `tool-list.json` is a valid json file. 

# bowtie2

Bowtie 2 is an ultrafast and memory-efficient tool for aligning sequencing reads to long reference sequences. It is particularly good at aligning reads of about 50 up to 100s or 1,000s of characters, and particularly good at aligning to relatively long (e.g. mammalian) genomes.

A stable version of bowtie2 can be downloaded via sourceforge from [here](https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.4.1/). 

## tool requirements
*If building bowtie2 from source:* 
* GCC, GNU Make and other basics
* Threading Building Blocks library (TBB)
* [zlib](https://www.zlib.net)

## steps
1. start by downloading the most recent version of bowtie2 from [sourceforge](https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.4.1/). If you are using linux we advise you obtain the binaries: `bowtie2-X.X.X-linux-x86_64` where X.X.X is the version number, for example 2.4.1.

2. Unzip the archive to the `<sap-path>/tools` directory. 

3. Test that the build works by changing to the bowtie2 directory and accessing the help menu:
```
cd <sap-path>/tools/bowtie2-X.X.X
./bowtie2 --help
```

4. Tell the SAP library where to find the bowtie2 executable, `bowtie2`. For this, you simply add the path of this file, *relative to the SAP library tools directory* to the file `<sap-path>/tools/tool-list.json` as in the following example:
```
"bcftools": "bowtie2-X.X.X/bowtie2",
```
Now, when the setup script is run, this tool will be accessible to the library. Remember that the quotes and comma at the end of the line in the above example are needed to ensure that `tool-list.json` is a valid json file. 

# bwa

BWA is a software package for mapping low-divergent sequences against a large reference genome, such as the human genome. It consists of three algorithms: BWA-backtrack, BWA-SW and BWA-MEM. The first algorithm is designed for Illumina sequence reads up to 100bp, while the rest two for longer sequences ranged from 70bp to 1Mbp.

A stable version of the package can be found, along with documentation, via [sourceforge](http://bio-bwa.sourceforge.net/)

## tool requirements

## steps
1. Download the latest version of the package from [sourceforge](https://sourceforge.net/projects/bio-bwa/files/) and extract the contents of the tar archive to the `<sap-path>/tools` directory with:
```
cd <your-downloads-path>
tar -C <sap-path>/tools -xvf bwa-X.X.X
```

# fastqc

# freebayes

# gatk v4

# gatk v3

# gmap

# samtools

Samtools is a suite of programs for interacting with high-throughput sequencing data, in particular: reading/writing/editing/indexing/viewing SAM/BAM/CRAM format files.

1. A stable version of the tool can be downloaded from [here](http://www.htslib.org/download/). Download the most recent version of `samtools-X.X.X`, where `X.X.X` will be some version number, for example `1.10.0`, and extract the contents of the tar archive to the `<sap-path>/tools` directory with:
	```
	cd <your-downoads-path>
	tar -C <sap-path>/tools -xvf samtools-X.X.X.tar.bz2
	```

2. Install the package to `<sap-path>/tools/samtools-X-X-X` using the `--prefix` option in the configure step:
	```
	cd <sap-path>/tools/samtools-X.X.X
	./configure --prefix=$(pwd)
	make
	make install
	```
There should now be build files inside the `samtools-X.X.X` folder. 

3. test the tool installed correctly by typing:
	```
	./samtools --help
	```
and you should see the help manual, which starts with something like this:
	```
	Program: samtools (Tools for alignments in the SAM format)
	Version: 1.10 (using htslib 1.10)
	...
	```

4. Tell the SAP library where to find the samtools executable, `samtools`. for this, you simply add the path of this file, *relative to the SAP library tools directory* to the file `<sap-path>/tools/tool-list.json` as in the following example:
	```
	  "samtools": "samtools-1.10/samtools",
	```
Now, when the setup script is run, this tool will be accessible to the library. Remember that the quotes and comma at the end are needed to ensure that `tool-list.json` is a valid json file. 

# stampy

# trimmomatic

# for developers: jekyll (to edit the docs)

