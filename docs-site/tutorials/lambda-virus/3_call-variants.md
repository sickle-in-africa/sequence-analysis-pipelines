---
layout: default
title: calling variants
parent: a lambda virus project
---

Joint-calling variants with the gatk
------------------------------------

In this final part of the tutorial, we will align the reads for the lambda virus samples in our cohort to the reference sequence, and then jointly call variants (SNPs and small indels) for the cohort. We will compare the estimated set of variants for the cohort against the true one from the simulation, and finish by considering different ways to improve the fidelity between the true and estimated variant sets. 

We will be using the pipeline `gatkall` in this pipeline, but the other wgs pipelines work similarly. Make sure you have all the correct tools installed for this pipeline. 

### steps

1. Move back to the inputs folder, `<sap-path>/inputs/wgs`. Locate the *basic* inputs file for the gatkall pipeline, `<sap-path>/inputs/wgs/gatkall.basic.json`, and open it in any text editor.

1. Navigate to the line:
	```
	  "threads": 8
	```
	in the file `<sap-path>/inputs/wgs/gatkall.basic.json` and change the value to match the number of cores you have access to.

1. Move back to the root directory, and run the main sap script `<sap-path>/sap.sh` with the following options:
	* *study type:* wgs
	* *pipeline id:* gatkall
	* *inputs id:* basic

	the call should look like this:
	```
	ch <sap-path>
	./sap.sh wgs gatkall basic
	```
	and you should see the following appear at the top of the output:
	```
checking a pipeline input json file was provided...
...done.
initializing pipeline input parameter values...
  setting default parameter values......done
  updating with arguments from input json file......done.
  checking that parameter values make sense......done
  checking that the necessary tools have been installed......done
==> output logs will be saved to <sap-path>/data/logs/f2d8cf8c/
==> temp. files will be saved to <sap-path>/data/temp/f2d8cf8c_
...done.
checking read file quality...
...done.
trimming read files...
	```

1. If the previous step ran successfully (i.e. you see "done" after each command and no "failed!" flags) then you can inspect the products of the pipeline run. Move to the fastqc directory, `<sap-path>/data/fastqc`, and inspect the .html files. These are your quality reports for the raw reads data. Select the quality report for the sample `s_1` and open it. In the summary table at the top, what was the mininum read lenth?

1. Move to the logs directory `<sap-path>/data/logs` and find the code from your run (you will see it at the top of the pipeline ouptut, quoted above). Check the log file for sample `s_2`. You will see here the output from the different pipeline tools. Locate the line:
```
TrimmomaticPE: Started with arguments:
```
This is the start of the output from the trimmomatic tool. Can you see, in the output from trimmmomatic, what percentage of input read pairs for sample `s_2` were dropped by the trimmomatic tool? 

1. Move to the vcf folder. Locate and open the file:
```
c_1.bwa.gatk.lambda-virus.genotyped.g.vcf
```
This is the cohort joint-called vcf file from the pipeline run. Read off the position of:
	* the first snp
	* the first insertion
	* the first deletion

	Now open the simulation's truth vcf file `c_1.truth.sorted.g.vcf` and examine the first few records.
	* Was the first snp in the joint-called file a true snp?
	* How about the first insertion?
	* The first deletion?

### possible errors

1. You may, during the pipeline run, the validate-bam-file step fails, i.e. you see the following appear in the terminal standard output:
```
  validating bam files......failed!
```
If you have run the same pipeline previously, before this run, then there is a chance that the bam files being validated are from a previous run, not the present one. It is helpful to clean the downstream data before repeating pipeline runs to eliminate this confusion. The clean-data pipeline can be used for this reason, but please *use it with caution* as the results are irreversible. For a how-to on downstream cleaning, see [this how-to]({% link how-tos/clean-data.md %}).

[Back]({% link tutorials/lambda-virus/2_simulate.md %})