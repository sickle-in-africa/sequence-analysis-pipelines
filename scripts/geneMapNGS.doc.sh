#!/usr/bin/env bash

function usage() {
echo -e """
\e[38;5;43mGeneMAP NGS Pipeline - Paired End\e[0m

A 'p' in front of a function means to run in parallel

Usage: geneMapNGS <command> [options] [h|help]

Commands:

-- Quality Processing
    fastqc/pfastqc
    trim/ptrim

 -- Mapping/Alignment
    gatkmap/pgatkmap
    bcfmap/pbcfmap

 -- Post-Mapping [Quality] Processing
    indelrealign
    bqsr/pbqsr

 -- Variant Calling
    bcfcall/pbcfcall
    gatkcall/pgatkcall

    -h,--help

    Enter help after 'geneMapNGS' to display full usage
    Example:
    geneMapNGS help

    Enter h or help after a command to show help for that command
    Example:
    geneMapNGS fastqc help
    geneMapNGS pgatkcall h
"""
}

function fqhelp() {
echo -e """
\e[38;5;43mGeneMAP NGS Pipeline - Paired End\e[0m [\e[38;5;43mQuality Processing - Fastqc\e[0m]

Usage: geneMapNGS [p]fastqc [options]

Options:
    -p,--path 	     		   : Path to fastq files. NB: Make sure all the files are paired i.e. forward
   				     and reverse. Enter '.' for current directory
                     		     [default: \e[38;5;2m"$(readlink -f `pwd`)"/\e[0m]

    -T,--threads     		   : Number of files to be processed simultaneously [default: 1]

				   : Examples:
				     geneMapNGS fastqc -p /some_directory/ -T 5
				     geneMapNGS pfastqc --path /some_directory/ --threads 10
"""
}

function trimhelp() {
echo -e """
\e[38;5;43mGeneMAP NGS Pipeline - Paired End\e[0m [\e[38;5;43mQuality Processing - Trimmomatic\e[0m]

Usage: geneMapNGS [p]trim [options]

Options:
    -p,--path        		   : Path to fastq files. NB: Make sure all the files are paired i.e. forward
                                     and reverse. Enter '.' for current directory
                     		     [default: \e[38;5;2m"$(readlink -f `pwd`)"/\e[0m]

    -a,--adapter     		   : trim with adapters by entering one of the following [default: none]
                     		     NB: Make sure you have the corresponding adapter file in your current 
				     directory

                     		     NP   --> NexteraPE-PE.fa
                     		     T3U  --> TruSeq3-PE-2.fa [Illumina universal]
                     		     T2P  --> TruSeq2-PE.fa
                     		     T2S  --> TruSeq2-SE.fa
                     		     T3P  --> TruSeq3-PE.fa
                     		     T3S  --> TruSeq3-SE.fa

    -l,--leadx       		   : Number of bases to clip from the 5' end [default: 0]
    -t,--trailx      		   : Number of bases to clip from the 3' end [default: 0]
    -T,--threads     		   : Number of files to be processed simultaneously [default: 1]

				   : Examples:
				     geneMapNGS trim -p /some_directory/ -T 4 -l 3 -t 30	(this will trim without adapter)
				    
				     geneMapNGS ptrim --path /some_directory/ --threads 10 --leadx 5 --trailx 15 --adapter NP
				    
				     This will attempt to trim with nextera adapter using the file 'NexteraPE-PE.fa'
				     If the file is not found in the current directory, it will attempt trimming without adapter
				     You can stop the pipeline by pressing on your keyboard 'CTRL+Z'
"""
}

function gmaphelp() {
echo -e """
\e[38;5;43mGeneMAP NGS Pipeline - Paired End\e[0m [\e[38;5;43mAlignment/Mapping - GATK-BWA\e[0m]

Usage: geneMapNGS [p]gatkmap [options]

Options:
    -s,--sample_list               : Sample/metadata file for FASTQ file. Contains four columns; 
                                     FR[forward read] RR[reverse read] SAM[sample name] PL[sequencing platform]
                                     Example:
                                             XXX_L001_R1_001.fastq.gz XXX_L001_R2_001.fastq.gz SAM1 ILLUMINA
                                             ERRXXXXXXXXXX_1.fastq.gz ERRXXXXXXXXXX_2.fastq.gz SAM2 ILLUMINA

    -r,--ref         		   : Reference FASTA file [required]

    -p,--path        		   : Path to fastq files. NB: Make sure all the files are paired i.e. forward 
                                     and reverse. Enter '.' for current directory
				     [default: \e[38;5;2m"$(readlink -f `pwd`)"/\e[0m]

    -T,--threads    		   : Number of files to be processed simultaneously [default: 1]
                                   : Examples:
                                     geneMapNGS gatkmap -r /path_to_reference/ref.fasta -p /path_to_gvcf_files/ -T 5
                                     geneMapNGS pgatkmap --ref /path_to_reference/ref.fasta --path /path_to_gvcf_files/ --threads 5
"""
}

function bmaphelp() {
echo -e """
\e[38;5;43mGeneMAP NGS Pipeline - Paired End\e[0m [\e[38;5;43mAlignment/Mapping - BWA-MEM\e[0m]

Usage: geneMapNGS [p]bcfmap [options]

Options:
    -r,--ref                       : Reference FASTA file [required]

    -p,--path                      : Path to fastq files. NB: Make sure all the files are paired i.e. forward
                                     and reverse. Enter '.' for current directory
                                     [default: \e[38;5;2m"$(readlink -f `pwd`)"/\e[0m]

    -T,--threads                   : Number of files to be processed simultaneously [default: 1]
                                   : Examples:
                                     geneMapNGS bcfmap -r /path_to_reference/ref.fasta -p /path_to_fastq_files/ -T 5
                                     geneMapNGS pbcfmap -ref /path_to_reference/ref.fasta --path /path_to_fastq_files/ --threads 5
"""
}

function realhelp() {
echo -e """
\e[38;5;43mGeneMAP NGS Pipeline - Paired End\e[0m [\e[38;5;43mInDel Realignment - GATK3\e[0m]

Usage: geneMapNGS [p]indelrealign [options]

Options:
    -b,--bam_list      		   : List containing BAM files if no FASTQ files are present. One BAM file per line
    -r,--ref           		   : Reference FASTA file [required]
    -d,--dbsnp         		   : A VCF file containig dbsnp rs IDs [optional]
    -k,--known_sites   		   : Reference VCF/dbSNP files containing sites or SNPs. You may provide multiple and
   				     separate with comma ','
                       		     E.g. dbsnp_138.hg19.vcf,1000G_phase1.indels.hg19.sites.vcf.gz
                       		     [default: NULL]

    -p,--path          		   : Path to BAM files. Enter '.' for current directory
                       		     [default: \e[38;5;2m"$(readlink -f `pwd`)"/\e[0m]

    -T,--threads       		   : Number of files to be processed simultaneously [default: 1]
				   : Examples:
				     geneMapNGS indelrealign -b bam.list -r ref.fasta -d dbsnpxxx.vcf -k dbsnp_138.hg19.vcf,1000G_phase1.indels.hg19.sites.vcf.gz -p /some_directory/ -T 5 
				      
				     \e[38;5;1mNB: Requires GATK3\e[0m
"""
}

function bqhelp() {
echo -e """
\e[38;5;43mGeneMAP NGS Pipeline - Paired End\e[0m [\e[38;5;43mBase Qulaity Score Recalibration\e[0m]

Usage: geneMapNGS [p]bqsr [options]

Options:
    -b,--bam_list                  : List containing BAM files if no FASTQ files are present. One BAM file per line
    -r,--ref                       : Reference FASTA file [required]
    -k,--known_sites               : Reference VCF/dbSNP files containing sites or SNPs. You may provide multiple and
                                     separate with comma ','
                                     E.g. dbsnp_138.hg19.vcf,1000G_phase1.indels.hg19.sites.vcf.gz
                                     [default: NULL]

    -p,--path                      : Path to BAM files. Enter '.' for current directory
                                     [default: \e[38;5;2m"$(readlink -f `pwd`)"/\e[0m]

    -T,--threads                   : Number of files to be processed simultaneously [default: 1]
                                   : Examples:
                                     geneMapNGS bqsr -b bam.list -r ref.fasta -d dbsnpxxx.vcf -k dbsnp_138.hg19.vcf,1000G_phase1.indels.hg19.sites.vcf.gz -p /some_directory/ -T 5 
"""
}

function gcallhelp() {
echo -e """
\e[38;5;43mGeneMAP NGS Pipeline - Paired End\e[0m [\e[38;5;43mVariant Calling - GATK\e[0m]

Usage: geneMapNGS [p]gatkcall [options]

Options:
    -b,--bam_list                  : List containing BAM files if no FASTQ files are present. One per line
    -g,--gvcf_list                 : List containing GVCF files if no FASTQ files are present. One per line
    -r,--ref                       : Reference FASTA file [required]
    -P,--ped                       : PED file containing family information [optional]. If provided, it should
    				     contain six columns
                                     FamilyID   IndividualID   PaternalID   MaternalID   Sex   Phenotype
                                     Example:
                                             FAM001   FAM001IND1   FAM001P1   FAM001P2   1   2
                                             FAM001   FAM001IND2   FAM001P1   FAM001P2   1   2
                                             FAM002   FAM002IND1   FAM002P1   FAM002P2   1   2
                                             FAM003   FAM003IND1   FAM003P1   FAM003P2   1   2

				     NB: Sex (1=male, 2=female, other=unknown)
				         Phenotype (1=unaffected, 2=affected, 0=missing, -9=missing)

    -d,--dbsnp                     : A VCF file containig dbsnp rs IDs [optional]
    -k,--known_sites               : Reference VCF/dbSNP files containing sites or SNPs. You may provide multiple
                                     and separate with comma ','
                                     E.g. dbsnp_138.hg19.vcf,1000G_phase1.indels.hg19.sites.vcf.gz
                                     [default: NULL]

    -p,--path                      : Path to BAM files. NB: Make sure all the files are paired i.e. forward and
                                     reverse
                                     NB: Enter '.' for current directory
                                     [default: \e[38;5;2m"$(readlink -f `pwd`)"/\e[0m]

    -T,--threads                   : Number of files to be processed simultaneously [default: 1]
    -o,--out			   : Output file name prefix [default: ngs]

                                   : Examples:
                                     geneMapNGS gatkcall -g gvcf.list -r /path_to_reference/ref.fasta -p /path_to_gvcf_files/ -T 5 -d /path_to_dbsnp_file/dbsnpxxx.vcf
                                     geneMapNGS pgatkcall -b bam.list -r /path_to_reference/ref.fasta -p /path_to_bam_files/ -T 5 -d /path_to_dbsnp_file/dbsnpxxx.vcf
"""
}

function bcallhelp() {
echo -e """
\e[38;5;43mGeneMAP NGS Pipeline - Paired End\e[0m [\e[38;5;43mVariant Calling - BCFTOOLS\e[0m]

Usage: geneMapNGS [p]bcfcall [options]

Options:
    -b,--bam_list                  : List containing BAM files if no FASTQ files are present. One per line
    -r,--ref                       : Reference FASTA file [required]
    -p,--path                      : Path to BAM files. NB: Make sure all the files are paired i.e. forward and
                                     reverse
                                     NB: Enter '.' for current directory
                                     [default: \e[38;5;2m"$(readlink -f `pwd`)"/\e[0m]

    -T,--threads                   : Number of files to be processed simultaneously [default: 1]
    -o,--out                       : Output file name prefix [default: ngs]

                                   : Examples:
                                     geneMapNGS bcfcall -b bam.list -r /path_to_reference/ref.fasta -p /path_to_bam_files/ -T 5
"""
}

function hlp() {
echo -e """
\e[38;5;43mGeneMAP NGS Pipeline - Paired End\e[0m

A 'p' in front of a function means to run in parallel

Usage: geneMapNGS <command> [options]

Commands:

 -- Quality Processing
    fastqc/pfastqc
    trim/ptrim

 -- Mapping/Alignment
    gatkmap/pgatkmap
    bcfmap/pbcfmap

 -- Post-Mapping [Quality] Processing
    indelrealign
    bqsr/pbqsr

 -- Variant Calling
    bcfcall/pbcfcall
    gatkcall/pgatkcall

\e[38;5;2mFastqc Options:\e[0m
--------------
    -p,--path                      : Path to fastq files. NB: Make sure all the files are paired i.e. forward
                                     and reverse. Enter '.' for current directory
                                     [default: \e[38;5;4m"$(readlink -f `pwd`)"/\e[0m]

    -T,--threads                   : Number of files to be processed simultaneously [default: 1]
                                   : Examples:
                                     geneMapNGS fastqc -p /some_directory/ -T 5
                                     geneMapNGS pfastqc --path /some_directory/ --threads 10

\e[38;5;2mTrimmomatic Options:\e[0m
-------------------
    -p,--path                      : Path to fastq files. NB: Make sure all the files are paired i.e. forward
                                     and reverse. Enter '.' for current directory
                                     [default: \e[38;5;4m"$(readlink -f `pwd`)"/\e[0m]

    -a,--adapter                   : trim with adapters by entering one of the following [default: none]
                                     NB: Make sure you have the corresponding adapter file in your current
                                     directory

                                     NP   --> NexteraPE-PE.fa
                                     T3U  --> TruSeq3-PE-2.fa [Illumina universal]
                                     T2P  --> TruSeq2-PE.fa
                                     T2S  --> TruSeq2-SE.fa
                                     T3P  --> TruSeq3-PE.fa
                                     T3S  --> TruSeq3-SE.fa

    -l,--leadx                     : Number of bases to clip from the 5' end [default: 0]
    -t,--trailx                    : Number of bases to clip from the 3' end [default: 0]
    -T,--threads                   : Number of files to be processed simultaneously [default: 1]
                                   : Examples:
                                     geneMapNGS trim -p /some_directory/ -T 4 -l 3 -t 30 (this will trim without adapter)
                                     geneMapNGS ptrim --path /some_directory/ --threads 10 --leadx 5 --trailx 15 --adapter NP

                                     This will attempt to trim with nextera adapter using the file 'NexteraPE-PE.fa'
                                     If the file is not found in the current directory, it will attempt trimming without adapter
                                     You can stop the pipeline by pressing on your keyboard 'CTRL+Z'

\e[38;5;2mGATK - Mapping/Alignment:\e[0m
-----------------------
    -s,--sample_list               : Sample/metadata file for FASTQ file. Contains four columns;
                                     FR[forward read] RR[reverse read] SAM[sample name] PL[sequencing platform]
                                     Example:
                                             XXX_L001_R1_001.fastq.gz XXX_L001_R2_001.fastq.gz SAM1 ILLUMINA
                                             ERRXXXXXXXXXX_1.fastq.gz ERRXXXXXXXXXX_2.fastq.gz SAM2 ILLUMINA

    -r,--ref                       : Reference FASTA file [required]

    -p,--path                      : Path to fastq files. NB: Make sure all the files are paired i.e. forward
                                     and reverse. Enter '.' for current directory
                                     [default: \e[38;5;4m"$(readlink -f `pwd`)"/\e[0m]

    -T,--threads                   : Number of files to be processed simultaneously [default: 1]
                                   : Examples:
                                     geneMapNGS gatkmap -r /path_to_reference/ref.fasta -p /path_to_gvcf_files/ -T 5
                                     geneMapNGS pgatkmap --ref /path_to_reference/ref.fasta --path /path_to_gvcf_files/ --threads 5

\e[38;5;2mBWA-MEM - Mapping/Alignment:\e[0m
---------------------------
    -r,--ref                       : Reference FASTA file [required]

    -p,--path                      : Path to fastq files. NB: Make sure all the files are paired i.e. forward
                                     and reverse. Enter '.' for current directory
                                     [default: \e[38;5;4m"$(readlink -f `pwd`)"/\e[0m]

    -T,--threads                   : Number of files to be processed simultaneously [default: 1]
                                   : Examples:
                                     geneMapNGS bcfmap -r /path_to_reference/ref.fasta -p /path_to_fastq_files/ -T 5
                                     geneMapNGS pbcfmap -ref /path_to_reference/ref.fasta --path /path_to_fastq_files/ --threads 5

\e[38;5;2mInDel Realignment Options:\e[0m
-------------------------
    -b,--bam_list                  : List containing BAM files if no FASTQ files are present. One BAM file per line
    -r,--ref                       : Reference FASTA file [required]
    -d,--dbsnp                     : A VCF file containig dbsnp rs IDs [optional]
    -k,--known_sites               : Reference VCF/dbSNP files containing sites or SNPs. You may provide multiple and
                                     separate with comma ','
                                     E.g. dbsnp_138.hg19.vcf,1000G_phase1.indels.hg19.sites.vcf.gz
                                     [default: NULL]

    -p,--path                      : Path to BAM files. Enter '.' for current directory
                                     [default: \e[38;5;4m"$(readlink -f `pwd`)"/\e[0m]

    -T,--threads                   : Number of files to be processed simultaneously [default: 1]
                                   : Examples:
                                     geneMapNGS indelrealign -b bam.list -r ref.fasta -d dbsnpxxx.vcf -k dbsnp_138.hg19.vcf,1000G_phase1.indels.hg19.sites.vcf.gz -p /some_directory/ -T 5

                                     \e[38;5;1mNB: Requires GATK3\e[0m

\e[38;5;2mBQSR Options:\e[0m
------------
    -b,--bam_list                  : List containing BAM files if no FASTQ files are present. One BAM file per line
    -r,--ref                       : Reference FASTA file [required]
    -d,--dbsnp                     : A VCF file containig dbsnp rs IDs [optional]
    -k,--known_sites               : Reference VCF/dbSNP files containing sites or SNPs. You may provide multiple and
                                     separate with comma ','
                                     E.g. dbsnp_138.hg19.vcf,1000G_phase1.indels.hg19.sites.vcf.gz
                                     [default: NULL]

    -p,--path                      : Path to BAM files. Enter '.' for current directory
                                     [default: \e[38;5;4m"$(readlink -f `pwd`)"/\e[0m]

    -T,--threads                   : Number of files to be processed simultaneously [default: 1]
                                   : Examples:
                                     geneMapNGS bqsr -b bam.list -r ref.fasta -d dbsnpxxx.vcf -k dbsnp_138.hg19.vcf,1000G_phase1.indels.hg19.sites.vcf.gz -p /some_directory/ -T 5

\e[38;5;2mGATK - Variant Call Options:\e[0m
---------------------------
    -b,--bam_list                  : List containing BAM files if no FASTQ files are present. One per line
    -g,--gvcf_list                 : List containing GVCF files if no FASTQ files are present. One per line
    -r,--ref                       : Reference FASTA file [required]
    -P,--ped                       : PED file containing family information [optional]. If provided, it should
                                     contain six columns
                                     FamilyID   IndividualID   PaternalID   MaternalID   Sex   Phenotype
                                     Example:
                                             FAM001   FAM001IND1   FAM001P1   FAM001P2   1   2
                                             FAM001   FAM001IND2   FAM001P1   FAM001P2   1   2
                                             FAM002   FAM002IND1   FAM002P1   FAM002P2   1   2
                                             FAM003   FAM003IND1   FAM003P1   FAM003P2   1   2

                                     NB: Sex (1=male, 2=female, other=unknown)
                                         Phenotype (1=unaffected, 2=affected, 0=missing, -9=missing)

    -d,--dbsnp                     : A VCF file containig dbsnp rs IDs [optional]
    -k,--known_sites               : Reference VCF/dbSNP files containing sites or SNPs. You may provide multiple
                                     and separate with comma ','
                                     E.g. dbsnp_138.hg19.vcf,1000G_phase1.indels.hg19.sites.vcf.gz
                                     [default: NULL]

    -p,--path                      : Path to BAM files. NB: Make sure all the files are paired i.e. forward and
                                     reverse
                                     NB: Enter '.' for current directory
                                     [default: \e[38;5;4m"$(readlink -f `pwd`)"/\e[0m]

    -T,--threads                   : Number of files to be processed simultaneously [default: 1]
    -o,--out                       : Output file name prefix [default: ngs]

                                   : Examples:
                                     geneMapNGS gatkcall -g gvcf.list -r /path_to_reference/ref.fasta -p /path_to_gvcf_files/ -T 5 -d /path_to_dbsnp_file/dbsnpxxx.vcf
                                     geneMapNGS pgatkcall -b bam.list -r /path_to_reference/ref.fasta -p /path_to_bam_files/ -T 5 -d /path_to_dbsnp_file/dbsnpxxx.vcf

\e[38;5;2mBCFTOOLS - Variant Call Options:\e[0m
-------------------------------
    -b,--bam_list                  : List containing BAM files if no FASTQ files are present. One per line
    -r,--ref                       : Reference FASTA file [required]
    -p,--path                      : Path to BAM files. NB: Make sure all the files are paired i.e. forward and
                                     reverse
                                     NB: Enter '.' for current directory
                                     [default: \e[38;5;4m"$(readlink -f `pwd`)"/\e[0m]

    -T,--threads                   : Number of files to be processed simultaneously [default: 1]
    -o,--out                       : Output file name prefix [default: ngs]

                                   : Examples:
                                     geneMapNGS bcfcall -b bam.list -r /path_to_reference/ref.fasta -p /path_to_bam_files/ -T 5

    -h,--help

 << Running One or More Steps >>
    -------------------------
    - Enter a command to run it
    - Enter two or more commands to run only those steps
    - Enter 'all' or 'pall' to run all steps in serial/parallel, one after another

    Examples:
    --------
    geneMapNGS fastqc [options]          : either FastQC in serial.
    geneMapNGS ptrim [options]           : trimmomatic in parallel
    geneMapNGS pmap [options]            : alignment/mapping with bwa-mem using GATKv4.x in parallel
    geneMapNGS bcfcall [options]         : variant calling with bcftools in seril
    geneMapNGS pgatkcall [options]       : joint (varaint) calling with GATKv4.x in parallel
    geneMapNGS pfastqc trim [options]    : fastqc in parallel, then trimmomatic in serial
    geneMapNGS pall [options]            : all the pipeline in parallel
"""
}
