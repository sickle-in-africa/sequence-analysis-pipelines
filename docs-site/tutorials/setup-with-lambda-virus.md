---
layout: default
title: a lambda virus project
parent: tutorials
---

A lambda virus project
=====================================================
{: .no_toc}

In this quick tutorial we will set up the SAP library for a sequencing project that uses the lambda virus genome saved in the references folder `<sap-path>/data/references/lambda-virus.fa`. We will learn how to:
* install the SAP library using the setup script
* simulate fake reads for an imaginary cohort with the `simulate` tool
* joint-call the variants from this simulated data set

1. We create a test project directory in our home directory:
```
cd ~
mkdir my-seq-projects
```

2. Then, we need to clone the SAP library, so we change to `my-seq-projects` and clone the repo:
```
cd ~/my-seq-projects
git clone https://github.com/sickle-in-africa/sequence-analysis-pipelines
```

3. Now, we want to perform whole genome sequencing on our lambda virus reads that we have received yesterday from Illumina. First, we need to chose which pipeline we will be using. Lets go with `gatkall`. This pipeline requires:
* gatk-4
* samtools
* bcftools
* fastqc
* trimmomatic
so we need to install these tools. For each tool on this list, we need to install it to the clone of the SAP library we obtained in step 2 of this example by carefully following [these]({% link how-tos/install-required-tools.md %}) instructions. 

4. Since the reference sequence we will already be using is already saved in our library, namely the file:
```~/my-seq-projects/sequence-analysis-pipelines/data/references/lambda-virus.fa
```
We can skip this step and move on to running the setup script.

5. change to the root directory, `~/my-seq-projects/sequence-analysis-pipelines` and run the setup script with options. Since it is a whole genome sequencing project we are doing, and we are aligning our reads to the sequence saved in `lambda-virus.fa`, we run:
```
cd ~/my-seq-projects/sequence-analysis-pipelines
./setup wgs lambda-virus.fa
```

6. Once this has run successfully, the data directory should be formatted. Your read files that you recieved yesterday from Illumina should be saved in the reads folder:
```
~/my-seq-projects/sequence-analysis-pipelines/data/reads
```

