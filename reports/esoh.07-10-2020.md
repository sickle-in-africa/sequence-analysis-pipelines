Technical Report
================

* Kevin Esoh
* 07-10-2020
* installing external tools on CHPC


## Libraries required for bcftools installation on a new linux computer

```
$ bunzip2 bcftools.1.11.X.tar.bz2
$ tar -xvf bcftools.1.11.X.tar
$ sudo apt install cmake
$ sudo apt-get install zlib1g-dev 
$ sudo apt-get install libbz2-dev 
$ sudo apt-get install -y liblzma-dev
$ sudo apt-get install libcurl4-openssl-dev
make
sudo make install
```
(worked on CHPC, without ```sudo``` of course)

## Samtools installation

```$ sudo apt-get install libncurses5 libncurses5-dev``` (could not get on the CHPC)
```
$ bunzip2 samtools.1.11.X.tar.bz2
$ tar -xvf samtools.1.11.X.tar
$ make all all-htslib
$ make install install-htslib
```
## freebayes: contains several submodules that need to be cloned as well
A simple ```git clone git@github.com:ekg/freebayes.git``` did not work for me

```git clone --recurse-submodules -j8 git@github.com:ekg/freebayes.git
```
works (-j option not avalable for git version on the CHPC)

```make -j4
```
(no cmake on CHPC)

## GATK 
Find GATK archives here https://console.cloud.google.com/storage/browser/gatk-software/package-archive/gatk/

### Install GATK3
Download latest version (GATK v3.8-1 from archive site)
```
$ bunzip2 package-archive_gatk_GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2
$ tar -xvf package-archive_gatk_GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar
```

Test GATK3
```
$ cd GenomeAnalysisTK-3.8-1-0-gf15c1c3ef
$ java -jar GenomeAnalysisTK.jar
```

Find GATK3 official docker here: https://hub.docker.com/r/broadinstitute/gatk3/tags

### GATK4
For a new Ubuntu 20 installation, python3 is the default, and python2 may not be installed
In this case 
```./gatk --help``` will fail
```sudo apt-get install python``` will resolve the issue

## Stampy
Get python-dev
```
$ sudo apt-get install python-dev
```
Now in the stampy directory, run
```
$ make 
```
This will build maptools.so required for stampy to run

### Stampy on CHPC
First activate python2 environment
```
$ module add chpc/BIOMODULES
$ module load python/2.7.15
```
Now in the stampy directory, run
```
make -j24
```
Commands run on the login node or lengau node directly without requesting for allocation of resources usually get killed.
``` -j24 ``` option helps speed up the process so it does not get killed before it completes.

### Important last installation
```
$ sudo apt install jq  
```
Else test fails with ```'includes/utilities.sh: line 102: jq: command not found'```

## Test Pipeline
```
$ ./setup.sh wgs lambda-virus
```

### Tutuorial - simulate reads
```
$ ./sap.sh wgs simulate-cohort basic
```
Requires perl modules
```
sudo cpan Math::Random.pm
sudo cpan JSON::Parse.pm
```
failed without superuser permision

### Finally, install parallel
```
$ sudo apt install parallel
```

### Runtime 
DELL intel® Core™ i7-8665U CPU @ 1.90GHz × 8 (512GB SSD), memory = 16GB (RAM)
```
- simulate-cohort basic: threads = 8, n_samples = 5
real	0m18.533s
user	1m16.816s
sys	0m0.136s

- simulate-cohort basic: threads = 8, n_samples = 15
real	0m51.903s
user	6m12.774s
sys	0m0.424s
```

gatkall
```
- quality control and variant calling: threads = 8, trim_minlen = 40
real	3m0.228s
user	21m26.388s
sys	0m41.100s

- quality control and variant calling: threads = 5, trim_minlen = 40
real	2m55.342s
user	20m18.777s
sys	0m38.120s

- quality control and variant calling: threads = 3, trim_minlen = 40
real	2m58.661s
user	19m16.933s
sys	0m34.663s

- quality control and variant calling: threads = 1, trim_minlen = 40
real	4m45.002s
user	12m15.260s
sys	0m24.858s
```
