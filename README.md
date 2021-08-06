# ABBA_BABA_2021

I received seperate VCF files for each of the chromosomes of Xenopus laevis. Going to filter them and do ABBA BABA analysis

Testing codes

Create a reference genome folder and download reference genome in that folder
```
wget http://ftp.xenbase.org/pub/Genomics/JGI/Xenla9.2/XENLA_9.2_genome.fa.gz

```
Unzip and index reference genome
```bash
gunzip XENLA_9.2_genome.fa

module load StdEnv/2020 samtools/1.12
samtools faidx XENLA_9.2_genome.fa
```

create dictionary file for reference genome

```bash
module load nixpkgs/16.09
module load gatk/4.1.2.0
gatk --java-options "-Xmx2G" CreateSequenceDictionary -R   XENLA_9.2_genome.fa
```
Create a saved scripts folder and save following script









