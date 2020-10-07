---
layout: default
title: setting up
parent: a lambda virus project
---

### Navigation
* [home](../../)
* [tutorials](../../tutorials/)
* [how-tos](../../how-tos/)

***

Setup the SAP library for a lambda-virus project
------------------------------------------------

In this first part of the tutorial we will see how to set up an instance of the SAP library --- format the data directory, connect installed tools, and build indices --- for a given project. Our project will be a whole genome sequencing project, and our reference sequence will be the lambda virus sequence, so these are both inputs for the SAP library setup script.


### steps

1. We first create a project directory in our home directory:
```
cd ~
mkdir my-seq-projects
```

1. Then, we clone the SAP library - change to `my-seq-projects` and clone the repo with git:
```
cd ~/my-seq-projects
git clone https://github.com/sickle-in-africa/sequence-analysis-pipelines
```
For the rest of the tutorial, we will use the following shorthand:
```
	<sap-path> : ~/my-seq-projects/sequence-analysis-pipelines
```

1. Now, we want to perform whole genome sequencing on our lambda virus reads that we will be simulating. First, we need to chose which pipeline to use. Lets go with `gatkall`. This pipeline requires:
	* gatk-4
	* samtools
	* bcftools
	* fastqc
	* trimmomatic

	so we need to install these tools. For each tool on this list, we need to install it to the instance of the SAP library we obtained in step 2 of this tutorial by carefully following [these](../../how-tos/install-required-tools.md) instructions. 

1. Navigate to the references directory, `<sap-path>/data/references`. The reference sequence we will be using is already saved in our library, namely the file:
```
	<sap-path>/data/references/lambda-virus.fa
```
Inspect it in an ordinary text editor. What does the comment line (the line starting with the > character) at the top say about the sequence? 

1. Change to the root directory, `<sap-path>` and run the setup script with the following options:
	* *study_type:* wgs
	* *reference basename:* lambda-virus.fa

	the call should look like this:
	```
	cd <sap-path>
	./setup wgs lambda-virus.fa
	```

1. Once this has run successfully, the data directory should be formatted. In the next section we will generate read files, which we will save in the new `reads`folder:
	```
	<sap-path>/data/reads
	```

[Back](index.md)

[Next](2_simulate.md)

### tutorial steps
1. [setting up the library](1_setup.md)
1. [simulating reads](2_simulate.md)
1. [calling variants](3_call-variants.md)
