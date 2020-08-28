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

3. Test the tool installed correctly by typing:
	```
	./bcftools --help
	```
and you should see the help manual, which starts with something like this:
	```
	Program: bcftools (Tools for variant calling and manipulating VCFs and BCFs)
	Version: 1.10.2 (using htslib 1.10.2)
	...
	```

4. Tell the SAP library where to find the bcftools executable, `bcftools`. for this, you simply add the path of this file, *relative to the SAP library tools directoty* to the file `<sap-path>/tools/tool-list.json` as in the following example:
	```
	  "bcftools": "bcftools-1.10.2/bcftools",
	```
Now, when the setup script is run, this tool will be accessible to the library. Remember that the quotes and comma at the end of the line in the above example are needed to ensure that `tool-list.json` is a valid json file. 

# bowtie2



# bwa

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