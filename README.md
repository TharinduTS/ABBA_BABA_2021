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
To increase efficiency, I renamed chr 9_10 as chr 9 so I can submit job arrays
```bash
find . -type f -name '*9_10*' | while read FILE ; do     newfile="$(echo ${FILE} |sed -e 's/9_10/9/g')" ;     mv "${FILE}" "${newfile}" ; done
```

Submit array of jobs for different chromosomes to create depth tables - make sure memory yoy reserve is at least 4 times the memory in '-Xmx'
Create a saved scripts folder and save following script as create_depth_table.sh
```bash
#!/bin/sh
#SBATCH --job-name=bwa_505
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=4:00:00
#SBATCH --mem=4gb
#SBATCH --output=bwa505.%J.out
#SBATCH --error=bwa505.%J.err
#SBATCH --account=def-ben
#SBATCH --array=1-9

#SBATCH --mail-user=premacht@mcmaster.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=ALL

# paste the same code twice , seperately for L and S

module load nixpkgs/16.09
module load gatk/3.8


# For L

for i in ./../XL_vcf_files/DB_new_chr${SLURM_ARRAY_TASK_ID}L_out_updated.vcf; do java -Xmx1G -jar $EBROOTGATK/GenomeAnalysisTK.jar -T VariantsToTable -R ./../reference_genome/XENLA_9.2_genome.fa -V $i -F CHROM -F POS -GF DP -o ${i#./../XL_vcf_files/}_GVCF_DP.table ; done


# For S

for i in ./../XL_vcf_files/DB_new_chr${SLURM_ARRAY_TASK_ID}S_out_updated.vcf; do java -Xmx1G -jar $EBROOTGATK/GenomeAnalysisTK.jar -T VariantsToTable -R ./../reference_genome/XENLA_9.2_genome.fa -V $i -F CHROM -F POS -GF DP -o ${i#./../XL_vcf_files/}_GVCF_DP.table ; done

```
Move dp tables to a seperate directory
```bash
mkdir ../dp_tables
mv ./*DP.table ../dp_tables/
```
Run the aove script for all the files

```bash
#!/bin/sh
#SBATCH --job-name=bwa_505
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --mem=32gb
#SBATCH --output=bwa505.%J.out
#SBATCH --error=bwa505.%J.err
#SBATCH --account=def-ben

#SBATCH --mail-user=premacht@mcmaster.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=ALL

module load nixpkgs/16.09
module load gcc/7.3.0
module load r

cp /scratch/premacht/laevis_GBS_2020/saved_scripts/cal_moving_dp_and_find_excludes.r .
Rscript cal_moving_dp_and_find_excludes.r 
```








