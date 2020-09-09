---
layout: default
title: installing required tools
parent: how-to's
---

How do I install the SAP library's required tools?
==================================================
{: .no_toc }

Installation instructions for these tools can all be found on their project source web pages, but crucially for this library to find the tools, it is required that you install each tool in the `<sap-path>/tools` library, and that you add the relative path to the file `<sap-path>/tools/tool-list.json`. Specific install instructions for each tool are given below. 

No pipeline requires every tool, and you may only need to install the ones that will be used for your project. See the pipeline library reference manual for details on which tools are required by which pipeline. 

The instructions below apply only to a 64-bit linux OS, the linux x86_64 platform. The command line versions only are given. 

1. TOC
{:toc}

## bcftools

BCFtools is a set of utilities that manipulate variant calls in the Variant Call Format (VCF) and its binary counterpart BCF. A stable version of the tool can be downloaded from [here](http://www.htslib.org/download/).

### tools required
{: .no_toc }
* cmake

### steps
{: .no_toc }
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

## bowtie2

Bowtie 2 is an ultrafast and memory-efficient tool for aligning sequencing reads to long reference sequences. It is particularly good at aligning reads of about 50 up to 100s or 1,000s of characters, and particularly good at aligning to relatively long (e.g. mammalian) genomes.

A stable version of bowtie2 can be downloaded via sourceforge from [here](https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.4.1/). 

### tool requirements
{: .no_toc }
*If building bowtie2 from source:* 
* GCC, GNU Make and other basics
* Threading Building Blocks library (TBB)
* [zlib](https://www.zlib.net)

### steps
{: .no_toc }
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
Now, when the setup script is run, this tool will be accessible to the library. Remember that the quotes and comma at the end of the line in the above example are needed to ensure that `tool-list.json` is a valid json file. If this is the last entry in the tool list, the final comma should be ommited. 

## bwa
BWA is a software package for mapping low-divergent sequences against a large reference genome, such as the human genome. It consists of three algorithms: BWA-backtrack, BWA-SW and BWA-MEM. The first algorithm is designed for Illumina sequence reads up to 100bp, while the rest two for longer sequences ranged from 70bp to 1Mbp.

A stable version of the package can be found, along with documentation, via [sourceforge](http://bio-bwa.sourceforge.net/)

### tool requirements
{: .no_toc }
* GNU make

### steps
{: .no_toc }
1. Download the latest version of the package from [sourceforge](https://sourceforge.net/projects/bio-bwa/files/) and extract the contents of the tar archive to the `<sap-path>/tools` directory with:
```
cd <your-downloads-path>
tar -C <sap-path>/tools -xvf bwa-X.X.X
```

2. Compile the package in `<sap-path>/tools/bwa-X-X-X` using GNU make:
```
cd <sap-path>/tools/bwa-X.X.X
make
```
There should now be build files inside the `bcftools-X.X.X` folder, and the bwa main executable file `bwa`. 

3. Test the build was successful by changing to the bwa directory and accessing the help menu:
```
cd <sap-path>/tools/bwa-X.X.X
./bwa
```

4. Tell the SAP library where to find the bwa executable, `<sap-path>/tools/bwa-X.X.X/bwa`. For this, you simply add the path of this file, *relative to the SAP library tools directory* to the file `<sap-path>/tools/tool-list.json` as in the following example:
```
"bcftools": "bwa-X.X.X/bwa",
```
Now, when the setup script is run, this tool will be accessible to the library. Remember that the quotes and comma at the end of the line in the above example are needed to ensure that `tool-list.json` is a valid json file. If this is the last entry in the tool list, the final comma should be ommited. 

## fastqc

FastQC aims to provide a simple way to do some quality control checks on raw sequence data coming from high throughput sequencing pipelines. It provides a modular set of analyses which you can use to give a quick impression of whether your data has any problems of which you should be aware before doing any further analysis.

A stable version of the package, along with specific help and documentation, can be found [here](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).

### tool requirements
{: .no_toc }
* java

### steps
{: .no_toc }
1. Download the latest version of fastqc from [here](https://www.bioinformatics.babraham.ac.uk/projects/download.html#fastqc) and extract the archive to the SAP library tools directory `<sap-path>/tools` with:
```
cd <your-downloads-path>
unzip fastqc_vX.X.X -d <sap-path>/tools
```
You should see a folder called `FastQC` in your tools folder with some java files (.jar files) and other things. These .jar files are already built, and so no compilation is necessary. 

2. Make the fastqc wrapper script `<sap-path>/tools/FastQC/fastqc` executable:
```
cd <sap-path>/tools/FastQC
chmod 755 ./fastqc
```

3. check that your java runtime enviroment is working by checking the verison:
```
java --version
```
and check that fastqc works by accessing the help:
```
./fastqc --help
```

4. Tell the SAP library where to find the fastqc executable, `<sap-path>/tools/FastQC/fastqc`by adding the path of this file, *relative to the SAP library tools directory* to the file `<sap-path>/tools/tool-list.json` as in the following example:
```
"fastqc" : "FastQC/fastqc",
```
Now, when the setup script is run, this tool will be accessible to the library. Remember that the quotes and comma at the end of the line in the above example are needed to ensure that `tool-list.json` is a valid json file. If this is the last entry in the tool list, the final comma should be ommited. 

## freebayes
freebayes is a Bayesian genetic variant detector designed to find small polymorphisms, specifically SNPs (single-nucleotide polymorphisms), indels (insertions and deletions), MNPs (multi-nucleotide polymorphisms), and complex events (composite insertion and substitution events) smaller than the length of a short-read sequencing alignment.

See the project github page [here](https://github.com/ekg/freebayes) for downloads and the manual. 

### tools required
{: .no_toc }
* git
* g++
* cmake
* the standard C and C++ development libraries
* liblzma
* pthread
* libbzip2.

### steps
{: .no_toc }
1. Download the latest version of the code, via git, to the SAP tools directory:
```
cd <sap-path>/tools
git clone https://github.com/ekg/freebayes
```

2. compile the package with cmake:
```
cd <sap-path>/tools/freebayes
make -j4
```
where 4 is the number of processors cmake will use to build the code (this value can be changed to suit your own machine).

3. Test that the build has worked by typing:
```
./freebayes --help
```
in the freebayes directory.

4. Tell the SAP library where to find the freebayes executable by adding the path of this file *relative to the SAP library tools directory* to the tool list json file `<sap-path>/tools/tool-list.json` as in the below example:
```
"freebayes" : "freebayes/bin/freebayes",
```
Note that the quotes, colon, and comma are needed to ensure that the tool list remains a valid json file. If this entry is the last in the tool list, the final comma must be ommitted. 

## gatk-4
Developed in the Data Sciences Platform at the Broad Institute, the genome analysis toolkit (gatk) offers a wide variety of tools with a primary focus on variant discovery and genotyping.

A stable version of the toolkit, as well as tutorials, manuals, best practice examples and further reading, can be found [here](https://gatk.broadinstitute.org/hc/en-us).

### tools required
{: .no_toc }
* java

### steps
{: .no_toc }
1. download the latest version of the toolkit from the Broad Institute github [here](https://github.com/broadinstitute/gatk/releases). Download the zip archive, `gatk-4.X.X.X.zip`, not the docker image. Unzip the archive to the SAP library tools directory `<sap-path>/tools/`:
```
cd <your-downloads-path>
unzip gatk-4.X.X.X -d <sap-path>/tools/
```
You should see some java (.jar) files in the uncompressed directory, and other files. These .jar files have already been built, and do not need compiling.

2. Test that your java enviroment and the toolkit is working by checking the java version and accessing the gatk help menu respectively, with:
```
java --version
```
and 
```
cd <sap-path>/tools/gatk-4.X.X.X
./gatk --help
```
3. Tell the SAP library where to find the main gatk executable by adding it to the tool list json file `<sap-path>/tools/tool-list.json` by adding a line like the following example:
```
"gatk": "gatk-4.X.X.X/gatk",
```
Note that the quotes, colon, and comma are needed to ensure that the tool list remains a valid json file. If this entry is the last in the tool list, the final comma must be ommitted. 

## gatk-3
*Still in development*

## gmap/gsnap
 gmap is a standalone program for mapping and aligning cDNA sequences to a genome. The program maps and aligns a single sequence with minimal startup time and memory requirements, and provides fast batch processing of large sequence sets.

A stable version of the program can be found [here](http://research-pub.gene.com/gmap/). 

### tool requirements
{: .no_toc }
* cmake

### steps
{: .no_toc }
1. download the latest version of gmap (which includes gsnap) from [here](http://research-pub.gene.com/gmap/) and extract the archive to the SAP library tools directory `<sap-path>/tools`:
```
cd <your-downloads-path>
tar -C <sap-path>/tools gmap-gsnap-20XX.XX.XX.tar.gz
```
where `20XX.XX.XX` is the date of the release, for example `2020-06-30`.  

2. Compile the program files using cmake. The object files and the final executable must be saved in the gsnap directory `<sap-path>/tools/gmap-20XX.XX.XX`, and you can specify this in the configure step with `$(pwd)`:
```
cd <sap-path>/tools/gmap-20XX.XX.XX
./configure --prefix=$(pwd)
make
make install
```
The gmap and gsnap executables as well as the other object files, will be in the bin directory `<sap-path>/tools/gmap-20XX.XX.XX/bin`. 

3. Test that the program installed correctly by accessing the help menu, with:
```
cd <sap-path>/tools/gmap-20XX.XX.XX/bin/gmap --help
```

4. Tell the SAP library where to find the gmap/gsnap executables by adding the path of each file to the tool list `<sap-path>/tools/tool-list.json` as in the following example:
```
"gmap": "gmap-20XX-XX-XX/bin/gmap",
"gmap_build": "gmap-20XX-XX-XX/bin/gmap_build",
"gsnap": "gmap-20XX-XX-XX/bin/gsnap",
```
Note that the quotes, colon, and comma are needed to ensure that the tool list remains a valid json file. If this entry is the last in the tool list, the final comma must be ommitted.

## samtools
Samtools is a suite of programs for interacting with high-throughput sequencing data, in particular: reading/writing/editing/indexing/viewing SAM/BAM/CRAM format files.

### tools required
{: .no_toc }
* cmake

### steps
{: .no_toc }
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
Now, when the setup script is run, this tool will be accessible to the library. Remember that the quotes and comma at the end are needed to ensure that `tool-list.json` is a valid json file. If this entry is the last in the tool list, the comma must be omitted. 

## stampy
Stampy is a package for the mapping of short reads from illumina sequencing machines onto a reference genome. It's recommended for most workflows, including those for genomic resequencing, RNA-Seq and Chip-seq. Stampy excels in the mapping of reads containing that contain sequence variation relative to the reference, in particular for those containing insertions or deletions.

Stable versions and documentation can be found [here](https://www.well.ox.ac.uk/research/research-groups/lunter-group/lunter-group/stampy). 

### tools required
{: .no_toc }
* python 2.6 or 2.7
* cmake

check your python version by typing:
```
python --version
```
and make sure that you are *not* in a custom python enviroment, like [Anaconda](https://www.anaconda.com/open-source) or one of your own. 

### steps
{: .no_toc }
1. Downoad the latest version [here](https://www.well.ox.ac.uk/research/research-groups/lunter-group/lunter-group/stampy) (you may be required to register before downloading) and unpack the archive to the SAP library tools directory `<sap-path>/tools` with:
```
cd <your-downloads-path>
tar -C <sap-path>/tools stampy-latest.tgz
```

2. change to the stampy directory, `<sap-path>/tools/stampy-X.X.X`, and compile with cmake:
```
cd <sap-path>/tools/stampy-X.X.X
make
```

3. Test that stampy is working on your machine by accessing the help menu:
```
cd <sap-path/tools/stampy-X.X.X
python stampy.py --help
```
4. Tell the SAP library where to find the samtools executable, `samtools`. for this, you simply add the path of the stampy python file, *relative to the SAP library tools directory* to the file `<sap-path>/tools/tool-list.json` as in the following example:
```
 "stampy" : "stampy-1.0.32/stampy.py",
```
Now, when the setup script is run, this tool will be accessible to the library. Remember that the quotes and comma at the end are needed to ensure that `tool-list.json` is a valid json file. If this entry is the last in the tool list, the comma must be omitted. 

## trimmomatic
Trimmomatic performs a variety of useful trimming tasks for illumina paired-end and single ended data.The selection of trimming steps and their associated parameters are supplied on the command line.

Stable versions of the package, as well as documentation, can be found [here](http://www.usadellab.org/cms/?page=trimmomatic).

### tools required
{: .no_toc }
* java

### steps
{: .no_toc }
1. Download the latest *binary* version of trimmomatic from [here](http://www.usadellab.org/cms/?page=trimmomatic), and extract the archive to the SAP Library tools directory with:
```
cd <your-downloads-path>
unzip Trimmomatic-X.X -d <sap-path>/tools
```
You should see, in the unpacked trimmomatic folder, a java .jar file and some other files. Since this is a .jar file it is already compiled and does not need to be built. 

2. Test that your version of java and the trimmomatic package are working by accessing your java version:
```
java --version
```
and the trimmomatic help menu
```
cd <sap-path>/tools/Trimmomatic-X.X
java -jar trimmomatic-X.X --help
```
respectively. 

4. Tell the SAP library where to find this .jar file, by adding its path *relative to the SAP library tool directory* to the tool list `<sap-path>/tools/tool-list.json` as in the following example:
```
"trimmomatic" : "trimmomatic-0.39/Trimmomatic-0.39.jar",
```
Now, when the setup script `<sap-path>/setup.sh` is run, this tool will be accessible to the library. Remember that the quotes and comma at the end are needed to ensure that `tool-list.json` is a valid json file. If this entry is the last in the tool list, the comma must be omitted. 

## for developers: jekyll
Jekyll is a static web site generator, often used for github project pages and static blogs. Pages can be written in simple markdown, and jekyll will convert them to html. Jekyll is a ruby gem, and needs a ruby enviroment to run. Jekyll is not necessary for any of the SAP pipelines, but it is helpful if you are a developer and want to contribute to the SAP library documentation. For more information see their [home page](https://jekyllrb.com/). 

### tools required
{: .no_toc }
* ruby

### steps
{: .no_toc }
1. check that a full ruby development enviroment has been installed on your computer, by checking the version:
```
ruby -v
```
(If no such enviroment exists this will throw an error). 

2. install jekyll bundler, with:
```
gem install jekyll bundler
```

3. test that you have correctly installed jekyll and bundler with:
```
jekyll --help
bundler --help
```

