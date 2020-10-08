Setting up SAP library on CHPC
==============================

1. move to somewhere in SADaCC/sequence-analysis-projects
```
$ cd /mnt/lustre/groups/CBBI1243/SADaCC/sequence-analysis-projects/jmorrice/
```
2. clone the git repo
```
$ git clone https://github.com/sickle-in-africa/sequence-analysis-pipelines.git
Cloning into 'sequence-analysis-pipelines'...
remote: Enumerating objects: 98, done.
remote: Counting objects: 100% (98/98), done.
remote: Compressing objects: 100% (78/78), done.
remote: Total 721 (delta 33), reused 56 (delta 18), pack-reused 623
Receiving objects: 100% (721/721), 44.55 MiB | 10.39 MiB/s, done.
Resolving deltas: 100% (257/257), done.

```
3. rename it something useful and move into the folder
```
$ mv sequence-analysis-pipelines/ wgs-project
$ cd wgs-project
```
4. create a data directory, and references directory
```
$ mkdir data
$ mkdir data/references
$ ls
data  docs  includes  LICENSE  pipelines  README.md  reports  sap.sh  setup.sh  test-suites  tools
```
5. create a soft link in `data/references` to the lambda-virus folder in `reference-sequence-data`:
```
$ cd data/references/
$ ln -s /mnt/lustre/groups/CBBI1243/SADaCC/reference-sequence-data/lambda-virus lambda-virus
$ ls
lambda-virus
```
6. move now to the project root folder and run the setup script
```
$ cd ../..
$ ./setup.sh wgs lambda-virus
./setup.sh: line 24: includes/config.sh: No such file or directory

==> any output logs will be saved to 97fd9059_setup.log
setting the project sub-directory locations hash...
...done.
setting the tool path hash...
cat: /mnt/lustre/groups/CBBI1243/SADaCC/sequence-analysis-projects/jmorrice/wgs-project/tools/tool-list.json: No such file or directory
cat: /mnt/lustre/groups/CBBI1243/SADaCC/sequence-analysis-projects/jmorrice/wgs-project/tools/tool-list.json: No such file or directory
cat: /mnt/lustre/groups/CBBI1243/SADaCC/sequence-analysis-projects/jmorrice/wgs-project/tools/tool-list.json: No such file or directory
cat: /mnt/lustre/groups/CBBI1243/SADaCC/sequence-analysis-projects/jmorrice/wgs-project/tools/tool-list.json: No such file or directory
cat: /mnt/lustre/groups/CBBI1243/SADaCC/sequence-analysis-projects/jmorrice/wgs-project/tools/tool-list.json: No such file or directory
cat: /mnt/lustre/groups/CBBI1243/SADaCC/sequence-analysis-projects/jmorrice/wgs-project/tools/tool-list.json: No such file or directory
cat: /mnt/lustre/groups/CBBI1243/SADaCC/sequence-analysis-projects/jmorrice/wgs-project/tools/tool-list.json: No such file or directory
cat: /mnt/lustre/groups/CBBI1243/SADaCC/sequence-analysis-projects/jmorrice/wgs-project/tools/tool-list.json: No such file or directory
cat: /mnt/lustre/groups/CBBI1243/SADaCC/sequence-analysis-projects/jmorrice/wgs-project/tools/tool-list.json: No such file or directory
cat: /mnt/lustre/groups/CBBI1243/SADaCC/sequence-analysis-projects/jmorrice/wgs-project/tools/tool-list.json: No such file or directory
cat: /mnt/lustre/groups/CBBI1243/SADaCC/sequence-analysis-projects/jmorrice/wgs-project/tools/tool-list.json: No such file or directory
cat: /mnt/lustre/groups/CBBI1243/SADaCC/sequence-analysis-projects/jmorrice/wgs-project/tools/tool-list.json: No such file or directory
cat: /mnt/lustre/groups/CBBI1243/SADaCC/sequence-analysis-projects/jmorrice/wgs-project/tools/tool-list.json: No such file or directory
...done.
formatting data directory for analysis...
  created /mnt/lustre/groups/CBBI1243/SADaCC/sequence-analysis-projects/jmorrice/wgs-project/data/bam
  created /mnt/lustre/groups/CBBI1243/SADaCC/sequence-analysis-projects/jmorrice/wgs-project/data/fastqc
  created /mnt/lustre/groups/CBBI1243/SADaCC/sequence-analysis-projects/jmorrice/wgs-project/data/logs
  created /mnt/lustre/groups/CBBI1243/SADaCC/sequence-analysis-projects/jmorrice/wgs-project/data/media
  created /mnt/lustre/groups/CBBI1243/SADaCC/sequence-analysis-projects/jmorrice/wgs-project/data/reads
  created /mnt/lustre/groups/CBBI1243/SADaCC/sequence-analysis-projects/jmorrice/wgs-project/data/sam
  created /mnt/lustre/groups/CBBI1243/SADaCC/sequence-analysis-projects/jmorrice/wgs-project/data/temp
  created /mnt/lustre/groups/CBBI1243/SADaCC/sequence-analysis-projects/jmorrice/wgs-project/data/vcf
  created /mnt/lustre/groups/CBBI1243/SADaCC/sequence-analysis-projects/jmorrice/wgs-project/inputs
  created /mnt/lustre/groups/CBBI1243/SADaCC/sequence-analysis-projects/jmorrice/wgs-project/inputs/wgs
...done.
checking which tools have been installed...
  checking bowtie2......not installed!
  checking samtools......not installed!
  checking bcftools......not installed!
  checking bgzip......not installed!
  checking jq......installed
  checking bwa......not installed!
  checking gsnap......not installed!
  checking stampy......not installed!
  checking gatk......not installed!
  checking fastqc......not installed!
  checking trimmomatic......not installed!
  checking freebayes......not installed!
...failed!
```
I get lots of errors at this stage. Thats ok.

7. run the setup script again in the same way
```
$ ./setup.sh wgs lambda-virus
./setup.sh: line 24: includes/config.sh: No such file or directory

==> any output logs will be saved to 2949e0c7_setup.log
setting the project sub-directory locations hash...
...done.
setting the tool path hash...
...done.
formatting data directory for analysis...
...done.
checking which tools have been installed...
  checking jq......installed
...done.
 writing includes/locations.sh
...done.
building indices for all installed tools...
...done.
```
8. we get another another error. We need to create a config.sh script:
```
$ cd includes
$ echo '''
> # load CHPC modules:
> module load chpc/perl/5.28.0
> module load chpc/gnu/parallel-20180622
> module load chpc/java/11.0.6
> ''' > config.sh
```
9. run the setup script again, the error should be gone
```
$ cd ../
$ ./setup.sh wgs lambda-virus
```
10. now we need to add some tools. they are already installed so now we just need to create soft links to them.
```
$ cd tools
$ ls
jaccard  jq-1.6  tool-list.json
```
11. from within the `wgs-project/tools` folder, create links to the tools
```
$ ln -s /mnt/lustre/groups/CBBI1243/SADaCC/sequence-analysis-tools/gatk-4.1.7.0 gatk-4.1.7.0
$ ln -s /mnt/lustre/groups/CBBI1243/SADaCC/sequence-analysis-tools/bcftools-1.11 bcftools-1.11
$ ln -s /mnt/lustre/groups/CBBI1243/SADaCC/sequence-analysis-tools/samtools-1.11 samtools-1.11
$ ln -s /mnt/lustre/groups/CBBI1243/SADaCC/sequence-analysis-tools/htslib-1.11 htslib-1.11
$ ln -s /mnt/lustre/groups/CBBI1243/SADaCC/sequence-analysis-tools/bwa bwa
$ ln -s /mnt/lustre/groups/CBBI1243/SADaCC/sequence-analysis-tools/sequence-simulator sequence-simulator
$ ln -s /mnt/lustre/groups/CBBI1243/SADaCC/sequence-analysis-tools/FastQC FastQC
$ ln -s /mnt/lustre/groups/CBBI1243/SADaCC/sequence-analysis-tools/Trimmomatic-0.39 Trimmomatic-0.39
```
12. add these tools to the tool list
```
$ nano tool-list.json
```
the file `tool-list.json` should look like this:
```
{
  "#": "external tools to be installed:",
  "simulate": "sequence-simulator/simulate.pl",
  "samtools": "samtools-1.11/bin/samtools",
  "bcftools": "bcftools-1.11/bin/bcftools",
  "bgzip": "htslib-1.11/bgzip",
  "fastqc": "FastQC/fastqc",
  "trimmomatic": "Trimmomatic-0.39/trimmomatic-0.39.jar",
  "bwa": "bwa/bwa",
  "gatk": "gatk-4.1.7.0/gatk"
}
```
13. now finally run the setup script again:
```
$ cd ../
$ ./setup.sh wgs lambda-virus

==> any output logs will be saved to e3086e3a_setup.log
setting the project sub-directory locations hash...
...done.
setting the tool path hash...
...done.
formatting data directory for analysis...
...done.
checking which tools have been installed...
  checking samtools......installed
  checking bcftools......installed
  checking bgzip......installed
  checking jq......installed
  checking bwa......installed
  checking gatk......installed
  checking fastqc......installed
  checking trimmomatic......installed
...done.
 writing includes/locations.sh
...done.
building indices for all installed tools...
  building bwa index...skipping as index exists...done
  building samtools index......done
  building gatk dict...skipping as dictionary exists...done
...done.

```
14. now create some input files for the pipelines:
```
$ cd inputs/wgs/
nano simulate-cohort.basic.json
```
add these contents to the `simulate-cohort.basic.json` file:
```
{
  "cohort_id": "c_1",
  "ref_label": "lambda-virus",
  "read_type": "pe",

  "threads": 8,
  "n_samples": 5,

  "has_snps": 1,
  "psnp": 0.001,
  "has_indels": 1,
  "findel": 0.0002,
  "has_rearr": 0,
  "nrearr_av": 3,
  "nreads": 10000,
  "rdlen_av": 75,
  "rdlen_min": 20
}
```
and create a pipeline inputs file too:
```
$ nano gatkall.basic.json
```
and add these contents to the file:
```
{
  "cohort_id": "c_1",
  "ref_label": "lambda-virus",
  "threads": 8,
  "trim_minlen": 40,
  "trim_av_qual_min": 1
}
```
and create another identical file:
```
$ cp gatkall.basic.json bcfall.basic.json
```
Now we are ready to run the pipelines. 

15. simulate the cohort reads:
```
$ ./sap.sh wgs simulate-cohort basic
checking simulation input json file was provided...
...done.
initializing simulation input parameter values...
  setting default parameter values......done
  updating with arguments from input json file......done
  checking that parameter values make sense......done
==> output logs will be saved to /mnt/lustre/groups/CBBI1243/SADaCC/sequence-analysis-projects/jmorrice/wgs-project/data/logs/b50095d7/
==> temp. files will be saved to /mnt/lustre/groups/CBBI1243/SADaCC/sequence-analysis-projects/jmorrice/wgs-project/data/temp/b50095d7_
...done.
simulating reads for the cohort...
  mutating the reference......done
  simulating reads......done
  sorting truth gvcf files......done
...done.
```
16. call variants with bcfall
```
$ ./sap.sh wgs bcfall basic
checking a pipeline input json file was provided...
...done.
initializing pipeline input parameter values...
  setting default parameter values......done
  updating with arguments from input json file......done.
  checking that parameter values make sense......done
  checking that the necessary tools have been installed......done
==> output logs will be saved to /mnt/lustre/groups/CBBI1243/SADaCC/sequence-analysis-projects/jmorrice/wgs-project/data/logs/dc6948fc/
==> temp. files will be saved to /mnt/lustre/groups/CBBI1243/SADaCC/sequence-analysis-projects/jmorrice/wgs-project/data/temp/dc6948fc_
...done.
checking read file quality with FASTQC...
...done.
trimming read files...
...done.
mapping reads to the reference with bwa...
  preparing reads list......done
  aligning reads to reference with bwa......done
  converting mapped sam files to bam......done
  sorting bam files......done
  adding read group information to bam files......done
  marking duplicates......done
  validating bam files......done
...done.
calling variants with bcftools...
  preparing list of bam files......done
  piling up sample gvcf files before joint calling......done
  indexing piledup cohort gvcf file......done
  joint-calling variants for cohort with bcftools......done
  indexing the joint-called cohort gvcf file......done
  unzipping the joint-called cohort gvcf file......done
...done.
```
17. call variants with gatkall
```
$ ./sap.sh wgs gatkall basic
checking a pipeline input json file was provided...
...done.
initializing pipeline input parameter values...
  setting default parameter values......done
  updating with arguments from input json file......done.
  checking that parameter values make sense......done
  checking that the necessary tools have been installed......done
==> output logs will be saved to /mnt/lustre/groups/CBBI1243/SADaCC/sequence-analysis-projects/jmorrice/wgs-project/data/logs/68e137f8/
==> temp. files will be saved to /mnt/lustre/groups/CBBI1243/SADaCC/sequence-analysis-projects/jmorrice/wgs-project/data/temp/68e137f8_
...done.
checking read file quality...
...done.
trimming read files...
...done.
performing Mapping/Alignment with GATKv4 and BWA ...
  preparing reads list......done
  converting fastq files to unmapped bam......done
  aligning reads to reference with bwa......done
  converting mapped sam files to bam......done
  sorting bam files......done
  merging aligned and unaligned bam files......done
  marking duplicates......done
  validating bam files......done
...done.
calling variants with gatk...
  preparing list of bam files......done
  indexing bam files......done
  calling sample variants with gatk HaplotypeCaller......done
  preparing gvcf list......done
  combining sample gvcf files......done
  genotyping cohort combined gvcf file......done
...done.
```
18. now check the final vcf files overlap a bit
```
$ cd data/vcf/
$ ls
c_1.bwa.gatk.lambda-virus.combined.g.vcf       c_1.bwa.gatk.lambda-virus.s_2.raw.g.vcf.idx  c_1.bwa.samtools.lambda-virus.genotyped.g.vcf
c_1.bwa.gatk.lambda-virus.combined.g.vcf.idx   c_1.bwa.gatk.lambda-virus.s_3.raw.g.vcf      c_1.bwa.samtools.lambda-virus.genotyped.g.vcf.gz
c_1.bwa.gatk.lambda-virus.genotyped.g.vcf      c_1.bwa.gatk.lambda-virus.s_3.raw.g.vcf.idx  c_1.bwa.samtools.lambda-virus.genotyped.g.vcf.gz.tbi
c_1.bwa.gatk.lambda-virus.genotyped.g.vcf.idx  c_1.bwa.gatk.lambda-virus.s_4.raw.g.vcf      c_1.bwa.samtools.lambda-virus.piledup.g.vcf.gz
c_1.bwa.gatk.lambda-virus.s_1.raw.g.vcf        c_1.bwa.gatk.lambda-virus.s_4.raw.g.vcf.idx  c_1.bwa.samtools.lambda-virus.piledup.g.vcf.gz.tbi
c_1.bwa.gatk.lambda-virus.s_1.raw.g.vcf.idx    c_1.bwa.gatk.lambda-virus.s_5.raw.g.vcf      c_1.truth.g.vcf
c_1.bwa.gatk.lambda-virus.s_2.raw.g.vcf        c_1.bwa.gatk.lambda-virus.s_5.raw.g.vcf.idx  c_1.truth.sorted.g.vcf

$ less -S c_1.truth.sorted.g.vcf
$ less -S c_1.bwa.gatk.lambda-virus.genotyped.g.vcf
$ less -S c_1.bwa.samtools.lambda-virus.genotyped.g.vcf
$ 
```
find a variant at the same position with the same REF/ALT in all three files to verify that it has worked. For example I found the following SNP:
```
NC_001416.1     1485    .       G       C 
```
in all three vcf files. 




