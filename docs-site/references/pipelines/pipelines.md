---
layout: default
title: pipelines
parent: references
---

# pipelines
{: .no_toc}

1. TOC
{:toc}

## gatkall

A pipeline for whole genome sequencing written by Kevin Esoh. The pipeline:
* performs fastqc read quality inspection
* trims reads for length and quality,
* aligns reads to a reference with *bwa mem*
* calls variants with the *gatk (v4) haplotype caller*
* runs in parallel across the samples in the input cohort

To use the pipeline, the following options must be supplied to the sap script:
* *study type:* wgs
* *pipeline id:* gatkall

A basic call of this pipeline is:
```
./sap.sh wgs gatkall <inputs id>
```

### tools used
gatkall uses the following external tools:
* gatk-4
* samtools
* bcftools
* fastqc
* trimmomatic

### input parameters
the input parameters for this pipeline, along with some example input values, are:
```
"cohort_id": "c_1",
"ref_base": "lambda-virus.fa",
"threads": 8,
"trim_minlen": 40,
"trim_av_qual_min": 1
```

## bcfall

A pipeline for whole genome sequencing written by Kevin Esoh. The pipeline:
* performs fastqc read quality inspection
* trims reads for length and quality,
* aligns reads to a reference with *bwa mem*
* calls variants with *bcftools*
* runs in parallel across the samples in the input cohort

To use the pipeline, the following options must be supplied to the sap script:
* *study type:* wgs
* *pipeline id:* bcfall

A basic call of this pipeline is:
```
./sap.sh wgs bcfall <inputs id>
```

### tools used
gatkall uses the following external tools:
* bcftools
* samtools
* bcftools
* fastqc
* trimmomatic

### input parameters
the input parameters for this pipeline, along with some example input values, are:
```
"cohort_id": "c_1",
"ref_base": "lambda-virus.fa",
"threads": 8,
"trim_minlen": 40,
"trim_av_qual_min": 1
```

## ngspipeline

A pipeline for whole genome sequencing. Pipeline adapted from a paper in Nature: *Comparing whole genome sequence variant callers*. Full citation: Hwang, K., Lee, I., Li, H. et al. Comparative analysis of whole-genome sequencing pipelines to minimize false negative findings. Sci Rep 9, 3219 (2019). https://doi.org/10.1038/s41598-019-39108-2.

* aligns reads to a reference with a choice of aligners:
	* bwa
	* stampy
	* bowtie2
	* gsnap
* calls variants with a choice of callers:
	* freebayes
	* samtools (bcftools)
* runs in parallel across the samples in the input cohort (with the exception of the *freebayes* aligner).

To use the pipeline, the following options must be supplied to the sap script:
* *study type:* wgs
* *pipeline id:* ngspipeline

A basic call of this pipeline is:
```
./sap.sh wgs ngspipeline <inputs id>
```

### tools used
ngspipeline uses the following external tools:
* gatk-4
* samtools

as well as, depending on the user choice of aligner:
* gsnap (optional, "aligner_id": "gsnap")
* bowtie2 (optional, "aligner_id": "bowtie2")
* bwa (optional, "aligner_id": "bwa")
* stampy (optional, also requires bwa to be installed, "aligner_id": "stampy")

and choice of aligner:
* freebayes (optional, "caller_id": "freebayes")
* bcftools (optional, "caller_id": "samtools")

Note that at least one aligner and one caller must be installed. 

### input parameters
the input parameters for this pipeline, along with some example input values, are:
```
"aligner_id": "stampy",
"caller_id": "samtools",
"ref_base": "lambda-virus.fa",
"cohort_id": "c_1",
"threads": 8,
"maxmem": "1G",
"recal_realign_on": "no"
```

Note: the "recal_realign_on" must be "no" at the moment, there is no other woking value. 

