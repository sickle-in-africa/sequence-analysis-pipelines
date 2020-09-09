---
layout: default
title: simulating reads
parent: a lambda virus project
---

Simulating reads for a cohort of lambda viruses
-----------------------------------------------

We have just set up the SAP library for a lambda-virus project. The data directory is preoperly formatted, and the indices for the lambda virus genome are built. Now, we need a cohort of sample reads. In this section we will simulate these reads. 

Before going further, it will be necessary to know the number of cores you have access to for analysis during this tutorial. On a linux computer we can simply type:
```
htop
```
to see the status of our processors. At the top of the output for this command you will see a number of bars:
```
1 [|||||||||||||||||||]
2 [|||||||||||||||||||]
...
```
each bar is a processor, and the levels show the activity on that processor. The number of these bars is the number of processors. 

### steps

1. determine the number of cores you have access to during this tutorial. For concreteness we will assume for the rest of the tutorial that you have 8, though it will in reality be different. 

2. Change to the pipeline inputs directory, `<sap-path>/inputs/wgs`. Locate the *basic* inputs file for the simulation pipeline, `<sap-path>/inputs/wgs/simulate-cohort.basic.json` and open in any text editor.

3. Navigate to the line:
```
  "threads": 8,
```
in file `<sap-path>/inputs/wgs/simulate-cohort.basic.json` and change the value to match the number of cores you have access to. Check that
```
  "n_samples"
```
is set to 5. 

4. navigate back to the root directory, `<sap-path>` and run the main sap script `<sap-path>/sap.sh` with the following options:
	* *study type:* wgs
	* *pipeline id:* simulate-cohort
	* *inputs id:* basic

	the call should look like this:
	```
	cd <sap-path>
	./sap.sh wgs simulate-cohort basic
	```
	and you should see the following output appear at the start of the pipeline run:
```
checking simulation input json file was provided...
...done.
initializing simulation input parameter values...
  setting default parameter values......done
  updating with arguments from input json file......done
  checking that parameter values make sense......done
==> output logs will be saved to <sap-path>/data/logs/397eebda/
==> temp. files will be saved to <sap-path>/data/temp/397eebda_
...done.
simulating reads for the cohort...
  mutating the reference......done
  simulating reads...
```

5. Go to the logs directory `<sap-path>/data/logs`. You will see a new directory has been created with a name that is a random 8 digit code. This is where all the output logs from the pipeline run go, and if there were any errors, they will be printed in detail here. What is the first line of the log file for sample `s_1` for this pipeline run?

6. Change to the reads directory, `<sap-path>/data/reads`. You will see here some new files have been created, for example:
```
c_1.s_1.raw_1P.fq.gz
c_1.s_1.raw_2P.fq.gz
c_1.s_2.raw_1P.fq.gz
c_1.s_2.raw_2P.fq.gz
...
```
these are the simulated read files. They are zipped, so to inspect them, we need to use the tool `zcat`. Run the followig command in the reads directory:
```
zcat c_1.s_1.raw_1P.fq.gz | head -20
```
how many individual reads have been printed (in whole or in part) to the terminal screen?

7. Finally, move to the vcf directory, `<sap-path>/data/vcf` and open the file:
```
c_1.truth.sorted.g.vcf
```
This is the true vcf file from the simulation, i.e. it describes the set of variants between the reference sequence, `<sap-path>/data/references/lambda-virus.fa` and the sequences of the simulated samples (there are no variants *between* the samples, they all represent the same mutated sequence but broken into non-equivalent read sets). In this truth vcf file, read off the position of:
	* the first snp
	* the first insertion
	* the first deletion


[Back]({% link tutorials/lambda-virus/1_setup.md %})

[Next]({% link tutorials/lambda-virus/3_call-variants.md %})