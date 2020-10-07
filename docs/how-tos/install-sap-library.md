---
layout: default
title: installing SAP library
parent: how-to's
---

### Navigation
* [home](../)
* [tutorials](../tutorials/)
* [how-tos](../how-tos/)

***

How do I install the SAP library?
=================================
{: .no_toc}
The SAP library can be downloaded from github [here](https://github.com/sickle-in-africa/sequence-analysis-pipelines) and installed by simply executing the setup script `<sap-library>/setup.sh` with options. The library has many external dependencies, like the [gatk](https://gatk.broadinstitute.org/hc/en-us) or [samtools](http://www.htslib.org/) and these must be installed *after* downloading the SAP library source code, but *before* running the setup script. Whenever a new external tool is added to the library, you can re-run the setup script, to connect it to the library. 

1. TOC
{: toc}

## tools required
See the how-to [*installing required tools*](install-required-tools.md) for a full list. The SAP library does *not* need every tool on this list installed to operate, and not every tool is required for a given pipeline; you should consult the references to see which tools are required by your chosen pipeline.

## steps
1. Download the library from [github](https://github.com/sickle-in-africa/sequence-analysis-pipelines) and save to somewhere you intend to do your analysis, for example `~/jack/projects/seq-project-1` using the following commands:
```
cd <your-genomics-project-path>
git clone https://github.com/sickle-in-africa/sequence-analysis-pipelines
```
Your sap path, `<sap-path>` will be the root directory of the git repository you have just cloned. In the above example, this would be: `~/Work/Projects/seq-project-1/sequence-analysis-pipelines`

2. Install the external tools you will need for your analysis. See [here](install-required-tools.md) for the full list, as well as installation instructions specific for each tool. If you are unsure which tools you may need, consult the references pages under your chosen analysis pipelines.


3. Save your project reference sequence, the species-specific sequence you will be aligning whole genome reads to for example, in the references directory `<sap-path>/data/references` with a sort name that uses no spaces or unusual characters (stick with letters, numbers, underscores and hypens to be safe). You will find a lambda virus genome sequence already there for training purposes: `lambda-virus.fa`.  

4. Run the setup script, `<sap-path>/setup.sh` with options, depending on your study. The basic usage of this script is:
```
cd <sap-path>
./setup.sh <study-type> <reference_basename>
```
where `<reference_basename>` is the local filename of the reference sequence file *without the extension*. So for example, if you are working on a whole genome sequencing project and will be using the lambda virus sequence already saved for training, then you would run:
```
cd <sap-path>
./setup.sh wgs lambda-virus
```
The setup script connects the library to any tools installed to the SAP library tools directory, builds all indexes for these tools from the input reference sequence---the lambda virus sequence in the above example, and formats the data directory `<sap-path>/data` for analysis. 

## comments
* the setup script can be repeatedly run, whenever a new tool is installed. It does not wipe or overwrite anything except the locations file `<sap-path>/includes/locations.sh` and so it will not delete your work.


