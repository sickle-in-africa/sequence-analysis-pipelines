---
layout: default
title: wgs analysis
parent: how-to's
---

### Navigation
* [home](../)
* [tutorials](../tutorials/)
* [how-tos](../how-tos/)

***

How do I perform joint calling on a cohort of sample reads?
============================================================

A major use of the SAP library is to perform joint calling of 

## tools required

the required tools for whole genome sequence variant calling depend on your choice of aligner and variant caller. The required tools for any workflow are:
* gatk v4
* samtools

The for the optional extra tools, see the [pipelines reference manual]({% link references/pipelines/pipelines.md %}) and chose a pipeline. The required tools for each pipeline are given there. 

## steps

1. make sure that all the necessary external tools are installed, as descibed in [this guide](install-required-tools.md)

1. make sure the SAP library has been properly set up for your chosen study type, whole genome sequencing, and for your projects chosen reference sequence, by following [this guide](install-sap-library.md). 

1. save the reference seqence in the references folder, `<sap-path>/data/references`, and make sure the basename matches the name that was supplied the SAP setup script `<sap-path>/setup.sh`. 

1. chose a suitable cohort id to track your downstream data files, for example:
```
c_001-illumina-cameroon
```
and save your read files in the reads directory `<sap-path>/data/reads` in the following format:
```
	<cohort-id>.<sample-id>.raw_1P.fq.gz
	<cohort-id>.<sample-id>.raw_2P.fq.gz
```
where "1P" and "2P" are the forward and reverse read files from the paried end reading. In the same folder, save a sample list file:
```
	<cohort-id>.samples.list
```
with the following records:
```
# sample file for cohort: c_1
s_1	c_1.s_1.raw_1P.fq.gz	c_1.s_1.raw_2P.fq.gz	ILLUMINA
s_2	c_1.s_2.raw_1P.fq.gz	c_1.s_2.raw_2P.fq.gz	ILLUMINA
s_3	c_1.s_3.raw_1P.fq.gz	c_1.s_3.raw_2P.fq.gz	ILLUMINA
s_4	c_1.s_4.raw_1P.fq.gz	c_1.s_4.raw_2P.fq.gz	ILLUMINA
...
```
where the first column are the sample ids (they do not have to be of the form shown but they must uniquely identify each sample) and the last column is the platform. This is just for annotating the downstream files, so the exact form is not important. The first line is a comment, and can be anything you like: it is to help you identify the contents of the file (you can have as many comment lines as you need). 

1. create an inputs file for your pipeline run. In the inputs wgs folder, `<sap-path>/inputs/wgs`, create a file with the name:
```
	<pipeline-id>.<inputs id>.json
```
where `<pipeline-id>` must match one of the pipeline ids for the pipelines saved in the library, for example:
```
gatkall, bcfall, ngspipeline
```
The inputs id can be anything you like (no spaces please!) but keep it short because it will be used when you call the pipeline in the following step. Use this file to set all the parameters you will need to run the pipeline. This includes the cohort id, as well as the number of threads. If you are using ngspipeline, then you must also chose your caller and aligner. Make sure your choices have been installed!

1. Change back to the root directory and run your pipeline with the sap script:
```
./sap wgs <pipeline-id> <inputs-id>
```

## comments
* the inputs id should be short and easy to identify. You may have multiple input sets for a single pipeline.