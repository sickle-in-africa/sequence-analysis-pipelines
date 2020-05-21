Running the GeneMAP NGS Pipeline
---

Get minimal usage
```
./geneMapNGS.sh [-h|--help]
```

Get full usage
```
./geneMapNGS.sh [h|help]
```

Get usage for a particular function
```
./geneMapNGS.sh fastqc [h|help]
```

FastQC: Quality Check
---
```
./geneMapNGS.sh fastqc --path /path_to_fastq_file/ --threads 5
```
> If the fastq files are in your current directory, then you can either choose to leave out --path, or enter ```--path .```

```
./geneMapNGS.sh pfastqc -p . -T 5
```
> This runs fastq in parallel, using the current directory as path to the fastq files

Trimmomatic
---
Get full usage
```
./geneMapNGS.sh trim help
```

Trim without adapter: in serial
```
./geneMapNGS.sh trim --path /path_to_fastq_file/ --threads 10 --leadx 3 --trailx 30
```
> This runs in serial, trimming 3 bases from the 5' end and 30 bases from the 3' end without any adapter

Trim with adapter: in parallel
```
./geneMapNGS.sh ptrim --path /path_to_fastq_file/ --threads 10 --leadx 3 --trailx 30 --adapter NP
```
> This runs in parallel (10 jobs at a time), trimming 3 bases from the 5' end and 30 bases from the 3' end, using the NexteraPE adapter

Alignment/Mapping: BWA-MEM
---
Get full usage
```
./geneMapNGS.sh bcfmap help
```

Run BWA MEM
```
./geneMapNGS.sh pbcfmap --path /path_to_fastq_file/ --ref reference.fasta --threads 8
```
> This runs BWA MEM in parallel (8 jobs at a time)

Alignment/Mapping: BWA-MEM-GATK
---
NB: This requires GATKv4.x

Get full usage
```
./geneMapNGS.sh gatkmap help
```

Run BWA-MEM-GATK
```
./geneMapNGS.sh pgatkmap --sample_list samples.txt --path /path_to_fastq_file/ --ref reference.fasta --threads 8

./geneMapNGS.sh gatkmap -r /path_to_reference/ref.fasta -p /path_to_gvcf_files/ -T 5
```
> See the usage for how the sample file should look like

Indel Realignment: GATKv3.x
---
NB: This requires GATKv3.x

Get full usage
```
./geneMapNGS.sh indelrealign h
```
Run
```
./geneMapNGS.sh indelrealign -b bam.list -r ref.fasta -d dbsnpxxx.vcf -k dbsnp_138.hg19.vcf,1000G_phase1.indels.hg19.sites.vcf.gz -p /some_directory/ -T 5
```

Base Quality Score Recalibration: BQSR - GATKv4.x
---
Get full usage
```
./geneMapNGS.sh pbqsr h
```

Run BQSR
```
./geneMapNGS.sh pbqsr -b bam.list -r ref.fasta -d dbsnpxxx.vcf -k dbsnp_138.hg19.vcf,1000G_phase1.indels.hg19.sites.vcf.gz -p /some_directory/ -T 5
```

Variant Calling: BCFTOOLS
---
Get full usage
```
./geneMapNGS.sh bcfcall h
```

Run
```
./geneMapNGS.sh bcfcall -b bam.list -r /path_to_reference/ref.fasta -p /path_to_bam_files/ -T 5
```


Variant Calling: GATKv4.x
---
Get full usage
```
./geneMapNGS.sh gatkcall h
```

Run with a list of BAM files
```
./geneMapNGS.sh pgatkcall -b bam.list -r /path_to_reference/ref.fasta -p /path_to_bam_files/ -T 5 -d /path_to_dbsnp_file/dbsnpxxx.vcf
```

Run with a list of GVCF files
```
./geneMapNGS.sh pgatkcall -g gvcf.list -r /path_to_reference/ref.fasta -p /path_to_gvcf_files/ -T 5 -d /path_to_dbsnp_file/dbsnpxxx.vcf
```


