# ABBA_BABA_2021

I received seperate VCF files for each of the chromosomes of Xenopus laevis. Going to filter them and do ABBA BABA analysis

working on cedar
```
/scratch/premacht/ABBA_BABA_updated
```

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

Submit array of jobs for different chromosomes to create depth tables - make sure memory you reserve is at least 4 times the memory in '-Xmx'
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
Save the following R script in saved scripts folder- This 'Cal moving average depath, cal cutoffs , plot mean_over_cuttoff and output a file with sites to exclude, ready to be used by VCFtools in the next step'

used loop to do this for all the chromosomes seperately

```r
# set working directory to current script directory 
# **uncomment this if you use this in local computer. KEEP COMMENTED OUT IF YOU ARE ON COMPUTECANADA) ***
# setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library (ggplot2)

all_files<-list.files(path = ".", pattern = "GVCF_DP.table")

for (chr_file in 1:length(all_files)) {
  


# get all the files with the site data
files <- all_files[chr_file]
my_data<-c()
temp<-c()

# read in the data and name the df based on the file name
for(f in 1:length(files)) {
  temp <- read.table(files[f], header = T)
  my_data <- rbind(my_data,temp)
  temp<-c()
}  

# my_data now has coverage per site info for each sample for all sites, genomewide
dim(my_data)
colnamez <- colnames(my_data)

# here is a function to calculate moving averages
# https://stackoverflow.com/questions/743812/calculating-moving-average
moving_fun <- function(x, w, FUN, ...) {
  # x: a double vector
  # w: the length of the window, i.e., the section of the vector selected to apply FUN
  # FUN: a function that takes a vector and return a summarize value, e.g., mean, sum, etc.
  # Given a double type vector apply a FUN over a moving window from left to the right, 
  #    when a window boundary is not a legal section, i.e. lower_bound and i (upper bound) 
  #    are not contained in the length of the vector, return a NA_real_
  if (w < 1) {
    stop("The length of the window 'w' must be greater than 0")
  }
  output <- x
  for (i in 1:length(x)) {
    # plus 1 because the index is inclusive with the upper_bound 'i'
    lower_bound <- i - w + 1
    if (lower_bound < 1) {
      output[i] <- NA_real_
    } else {
      output[i] <- FUN(x[lower_bound:i, ...])
    }
  }
  output
}

# example
# v <- seq(1:10)

# compute a MA(2)
# moving_fun(v, 2, mean)


# make a new dataframe that has the moving average for each sample throughout the 
# genome.  No worries if the window goes across chromosomes
mv_ave_df2 <- c()
for (i in 3:length(my_data)){  
  # calculate the moving average for each column
  moving_average <- moving_fun(my_data[,i], 50, mean)
  # add this to a new dataframe
  mv_ave_df2 <- cbind(mv_ave_df2,moving_average)
  # rename the column to match the sample
  colnames(mv_ave_df2)[i-2] = colnamez[i]
}  

# add chromosome and position data to mv_ave_df2
mv_ave_df3 <- data.frame(mv_ave_df2,my_data$CHROM,my_data$POS)


colnames(mv_ave_df3)[ncol(mv_ave_df3)-1] <- "CHROM"
colnames(mv_ave_df3)[ncol(mv_ave_df3)] <- "POS"


# calculate mean depth and sd per sample
mean_depth<- c()
sd_depth<- c()

for (i in 3:length(my_data)){
  x<- mean(my_data[,i],na.rm=T)
  mean_depth<-append(mean_depth,x, after = length(mean_depth))
  y<-sd(my_data[,i],na.rm=T)
  sd_depth<-append(sd_depth,y, after = length(sd_depth))
}  

cutoff_vector <- c()


# identify chr and positions in any sample that are >4 sd above that samples mean coverage
# based on the rolling average
# first make a vector with cutoff values for each sample
for (i in 3:ncol(mv_ave_df3)){
  cutoff <- mean_depth[i-2] + 4*sd_depth[i-2]
  cutoff_vector <- append(cutoff_vector,cutoff, after=length(cutoff_vector))
}

mean_over_cuttoff <- mean_depth/cutoff_vector
cutoff_data<-data.frame(1:length(mean_over_cuttoff),mean_over_cuttoff)
p<-ggplot(cutoff_data)+
  geom_point(aes(x=X1.length.mean_over_cuttoff.,y=mean_over_cuttoff))+
  theme_bw()

ggsave(filename = gsub("GVCF_DP.table","mean_over_cuttoff_plot.pdf",files),plot = p,width = 30,height = 10)
# give the elements in cutoff_vector some names

#manual
#names(cutoff_vector) <- c('bru_PF707.DP','download.DP')


# automated *******************************************
#get sample names from mv_ave_df3
xx<-names(mv_ave_df3)

# select just sample names
yy<-xx[1:(length(xx)-2)]

# take sample names for vector names
names(cutoff_vector)<-yy
# zz<-sub("_cuttrim_sorted_final_l_only.bam.DP","",yy)

# ****************************************************

# now cycle through each column and identify chr and pos of the bad ones
# manual
#subtest <- NULL
#subtest <- subset(mv_ave_df3, 
#              JM_no_label1_Draken_CCACGT_cuttrim_sorted_final_l_only.bam.DP > cutoff_vector["JM_no_label1_Draken_CCACGT_cuttrim_sorted_final_l_only.bam.DP"]  | 
#                JM_no_label2_Draken_TTCAGA_cuttrim_sorted_final_l_only.bam.DP > cutoff_vector["JM_no_label2_Draken_TTCAGA_cuttrim_sorted_final_l_only.bam.DP"] 
# )

# automated - Does the same thing as above**********
sub <- NULL

sub<-subset(mv_ave_df3,
            get(yy)>min(cutoff_vector),
            )


# **************************************************

# automated - only collects sites with a higher depth than that samples own cutoff*****************************************
#sub <- NULL

#for (i in yy) {
  
#  sub <- subset(mv_ave_df3,
#               get(i) >cutoff_vector[i] 
#  )
#}
# **************************************************

dim(sub)
dim(mv_ave_df3)


to_file <- data.frame(sub$CHROM, sub$POS)

# write to file
write.table(to_file, gsub("GVCF_DP.table","positions_to_exclude.txt",files), append = FALSE, sep = "\t", dec = ".",
            row.names = F, col.names = F,quote = FALSE)

}
```

Run the above script for all the files- use this in the folder with dp table files

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

cp ./../saved_scripts/cal_moving_dp_and_find_excludes.r .
Rscript cal_moving_dp_and_find_excludes.r 
```

put plots and positions to exclude in seperate directories
```bash
mkdir cutoff_plots
mkdir positions_to_exclude

mv dp_tables/*.pdf cutoff_plots/
mv dp_tables/*exclude.txt positions_to_exclude/
```
make a directory for positions excluded vcfs

```bash
mkdir positions_excluded
```

Filter selected sites submiting an array of jobs - writing this in saved scripts

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
module load intel/2018.3
module load vcftools/0.1.16


# For L

vcftools --vcf ./../XL_vcf_files/DB_new_chr${SLURM_ARRAY_TASK_ID}L_out_updated.vcf --out ./../positions_excluded/DB_new_chr${SLURM_ARRAY_TASK_ID}L_out_updated_positions_excluded.vcf --exclude-positions ./../positions_to_exclude/DB_new_chr${SLURM_ARRAY_TASK_ID}L_out_updated.vcf_positions_to_exclude.txt --recode


# For S

vcftools --vcf ./../XL_vcf_files/DB_new_chr${SLURM_ARRAY_TASK_ID}S_out_updated.vcf --out ./../positions_excluded/DB_new_chr${SLURM_ARRAY_TASK_ID}S_out_updated_positions_excluded.vcf --exclude-positions ./../positions_to_exclude/DB_new_chr${SLURM_ARRAY_TASK_ID}S_out_updated.vcf_positions_to_exclude.txt --recode

```
## Admixture

Create directory for admixture

```bash
mkdir ADMIXTURE
cd ADMIXTURE
```
copy filtered vcfs here
```bash
cp -r ../positions_excluded/ .
```
making directories for all chromosome data, l only, s only
```bash
mkdir all
mkdir l_only
mkdir s_only
```
inside all,
make directories
```bash 
mkdir filtered_thinned_VCFs
mkdir filtered_again_VCFs
mkdir filtered_VCFs
mkdir scripts
mkdir combined_files
```
then copy all filtered vcfs in filtered_VCFs

try filtering the data to remove positions that have >50% missing data.  This might decrease the size of the data substantially.  If the file sizes are way smaller (e.g. half as large) then no need to thin any more- run this in scripts folder
```bash
#!/bin/sh
#SBATCH --job-name=bwa_505
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=12:00:00
#SBATCH --mem=32gb
#SBATCH --output=bwa505.%J.out
#SBATCH --error=bwa505.%J.err
#SBATCH --account=def-ben

module load vcftools/0.1.16

for i in ../filtered_VCFs/*.vcf;do
j=${i#../filtered_VCFs/}

vcftools --vcf ${i} --max-missing 0.5 --out ../filtered_again_VCFs/${j%.vcf.recode.vcf}_missing_filtered.vcf --recode ;done
```

now thin the files to reduce the time needed for calculations. - If the size of the filtered files are not small enough

save this in scripts folder and run it

```bash
#!/bin/sh
#SBATCH --job-name=bwa_505
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=1:00:00
#SBATCH --mem=32gb
#SBATCH --output=bwa505.%J.out
#SBATCH --error=bwa505.%J.err
#SBATCH --account=def-ben

module load vcftools/0.1.16

for i in ../filtered_VCFs/*.vcf;do
j=${i#../filtered_VCFs/}

vcftools --vcf ${i} --thin 1000 --out ../filtered_thinned_VCFs/${j%.vcf.recode.vcf}_thinned.vcf --recode ;done
```

use bcftools to combine chromosomes into one file to feed into plink(inside filtered_thinned_VCFs)
you can get the list of files seperated by space by,
```bash
ls *.vcf | tr "\n" " "
```
then,

```bash
module load StdEnv/2020  gcc/9.3.0 bcftools/1.10.2
bcftools concat -o ../combined_files/autosomes.vcf DB_new_chr1L_out_updated_positions_excluded_missing_filtered.vcf.recode.vcf DB_new_chr1S_out_updated_positions_excluded_missing_filtered.vcf.recode.vcf DB_new_chr2L_out_updated_positions_excluded_missing_filtered.vcf.recode.vcf DB_new_chr2S_out_updated_positions_excluded_missing_filtered.vcf.recode.vcf DB_new_chr3L_out_updated_positions_excluded_missing_filtered.vcf.recode.vcf DB_new_chr3S_out_updated_positions_excluded_missing_filtered.vcf.recode.vcf DB_new_chr4L_out_updated_positions_excluded_missing_filtered.vcf.recode.vcf DB_new_chr4S_out_updated_positions_excluded_missing_filtered.vcf.recode.vcf DB_new_chr5L_out_updated_positions_excluded_missing_filtered.vcf.recode.vcf DB_new_chr5S_out_updated_positions_excluded_missing_filtered.vcf.recode.vcf DB_new_chr6L_out_updated_positions_excluded_missing_filtered.vcf.recode.vcf DB_new_chr6S_out_updated_positions_excluded_missing_filtered.vcf.recode.vcf DB_new_chr7L_out_updated_positions_excluded_missing_filtered.vcf.recode.vcf DB_new_chr7S_out_updated_positions_excluded_missing_filtered.vcf.recode.vcf DB_new_chr8L_out_updated_positions_excluded_missing_filtered.vcf.recode.vcf DB_new_chr8S_out_updated_positions_excluded_missing_filtered.vcf.recode.vcf DB_new_chr9L_out_updated_positions_excluded_missing_filtered.vcf.recode.vcf DB_new_chr9S_out_updated_positions_excluded_missing_filtered.vcf.recode.vcf

```
do the same for l only and s only 

compress and index(inside combined_files)
```bash
bgzip -c autosomes.vcf > autosomes.vcf.gz
tabix -p vcf autosomes.vcf.gz
```

convert to geno format using plink , make a bed file and remove any SNP with no data for autosomes:

```bash
module load nixpkgs/16.09  intel/2016.4 plink/1.9b_5.2-x86_64

plink --vcf ./autosomes.vcf.gz --make-bed --geno 0.999 --out ./autosomes --allow-extra-chr --const-fid
```
do all these for l_only and s_only seperately

we need to change the chr names in the .bim file because these cause problems for admixture:

```
awk -v OFS='\t' '{$1=0;print $0}' autosomes.bim > autosomes.bim.tmp
mv autosomes.bim.tmp autosomes.bim
```
create a directory to collect outputs from next step
```bash
mkdir outs_array
```
now run two job arrays to cal admixture(2:6 and 7:12-run array num,bers acccordingly)
```bash
PLINK v1.90b4.6 64-bit (15 Aug 2017)           www.cog-genomics.org/plink/1.9/
(C) 2005-2017 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to ./autosomes.log.
Options in effect:
  --allow-extra-chr
  --const-fid
  --geno 0.999
  --make-bed
  --out ./autosomes
  --vcf ./autosomes.vcf.gz

257860 MB RAM detected; reserving 128930 MB for main workspace.
Allocated 7259 MB successfully, after larger attempt(s) failed.
--vcf: ./autosomes-temporary.bed + ./autosomes-temporary.bim +
./autosomes-temporary.fam written.
2823721 variants loaded from .bim file.
96 people (0 males, 0 females, 96 ambiguous) loaded from .fam.
Ambiguous sex IDs written to ./autosomes.nosex .
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 96 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate is 0.821702.
0 variants removed due to missing genotype data (--geno).
2823721 variants and 96 people pass filters and QC.
Note: No phenotypes present.
--make-bed to ./autosomes.bed + ./autosomes.bim + ./autosomes.fam ... done.
[premacht@cedar5 combined_files]$ clear

[premacht@cedar5 combined_files]$ ls
autosomes.bed  autosomes.log    autosomes.vcf.gz               cal_admixture_array_7_to_12.sh
autosomes.bim  autosomes.nosex  autosomes.vcf.gz.tbi           outs_array
autosomes.fam  autosomes.vcf    cal_admixture_array_2_to_6.sh
[premacht@cedar5 combined_files]$ awk -v OFS='\t' '{$1=0;print $0}' autosomes.bim > autosomes.bim.tmp
[premacht@cedar5 combined_files]$ mv autosomes.bim.tmp autosomes.bim
[premacht@cedar5 combined_files]$ ls
autosomes.bed  autosomes.log    autosomes.vcf.gz               cal_admixture_array_7_to_12.sh
autosomes.bim  autosomes.nosex  autosomes.vcf.gz.tbi           outs_array
autosomes.fam  autosomes.vcf    cal_admixture_array_2_to_6.sh
[premacht@cedar5 combined_files]$ less cal_admixture_array_2_to_6.sh 
[premacht@cedar5 combined_files]$ sbatch cal_admixture_array_2_to_6.sh 
Submitted batch job 14392663
[premacht@cedar5 combined_files]$ sbatch cal_admixture_array_7_to_12.sh 
Submitted batch job 14392666
[premacht@cedar5 combined_files]$ ls
autosomes.bed         bwa505.14392663.out  bwa505.14393673.err  bwa505.14393749.out
autosomes.bim         bwa505.14392666.err  bwa505.14393673.out  bwa505.14393750.err
autosomes.fam         bwa505.14392666.out  bwa505.14393684.err  bwa505.14393750.out
autosomes.log         bwa505.14393648.err  bwa505.14393684.out  cal_admixture_array_2_to_6.sh
autosomes.nosex       bwa505.14393648.out  bwa505.14393719.err  cal_admixture_array_7_to_12.sh
autosomes.vcf         bwa505.14393649.err  bwa505.14393719.out  outs_array
autosomes.vcf.gz      bwa505.14393649.out  bwa505.14393725.err
autosomes.vcf.gz.tbi  bwa505.14393671.err  bwa505.14393725.out
bwa505.14392663.err   bwa505.14393671.out  bwa505.14393749.err
[premacht@cedar5 combined_files]$ less bwa505.14392666.err
[premacht@cedar5 combined_files]$ less cal_admixture_array_2_to_6.sh 
[premacht@cedar5 combined_files]$ module load StdEnv/2020 admixture/1.3.0

Lmod is automatically replacing "nixpkgs/16.09" with "gentoo/2020".

Lmod has detected the following error:  The following module(s) are unknown: "gcccore/.9.3.0"

Please check the spelling or version number. Also try "module spider ..."
It is also possible your cache file is out-of-date; it may help to try:
  $ module --ignore-cache load "gcccore/.9.3.0"

Also make sure that all modulefiles written in TCL start with the string #%Module



[premacht@cedar5 combined_files]$ module --ignore-cache load "gcccore/.9.3.0"
Lmod has detected the following error:  The following module(s) are unknown: "gcccore/.9.3.0"

Please check the spelling or version number. Also try "module spider ..."
It is also possible your cache file is out-of-date; it may help to try:
  $ module --ignore-cache load "gcccore/.9.3.0"

Also make sure that all modulefiles written in TCL start with the string #%Module



[premacht@cedar5 combined_files]$ module spider admixture

-----------------------------------------------------------------------------------------------------------
  admixture: admixture/1.3.0
-----------------------------------------------------------------------------------------------------------
    Description:
      ADMIXTURE is a software tool for maximum likelihood estimation of individual ancestries from
      multilocus SNP genotype datasets. It uses the same statistical model as STRUCTURE but calculates
      estimates much more rapidly using a fast numerical optimization algorithm.

    Properties:
      Bioinformatic libraries/apps / Logiciels de bioinformatique

    You will need to load all module(s) on any one of the lines below before the "admixture/1.3.0" module is available to load.

      StdEnv/2020
      nixpkgs/16.09
 
    Help:
      
      Description
      ===========
      ADMIXTURE is a software tool for maximum likelihood estimation of individual ancestries from 
      multilocus SNP genotype datasets. It uses the same statistical model as STRUCTURE but calculates estimates much 
      more rapidly using a fast numerical optimization algorithm.
      
      
      More information
      ================
       - Homepage: https://www.genetics.ucla.edu/software/admixture
      


 

[premacht@cedar5 combined_files]$ nixpkgs/16.09 admixture/1.3.0
-bash: nixpkgs/16.09: No such file or directory
[premacht@cedar5 combined_files]$ module load nixpkgs/16.09 admixture/1.3.0
[premacht@cedar5 combined_files]$ module list

Currently Loaded Modules:
  1) nixpkgs/16.09  (S)   3) icc/.2016.4.258   (H)   5) intel/2016.4 (t)      7) plink/1.9b_5.2-x86_64 (bio)
  2) gcccore/.5.4.0 (H)   4) ifort/.2016.4.258 (H)   6) gsl/2.6      (math)   8) admixture/1.3.0       (bio)

  Where:
   S:     Module is Sticky, requires --force to unload or purge
   bio:   Bioinformatic libraries/apps / Logiciels de bioinformatique
   math:  Mathematical libraries / Bibliothèques mathématiques
   t:     Tools for development / Outils de développement
   H:                Hidden Module

Inactive Modules:
  1) htslib/1.10.2   2) bcftools/1.10.2

 

[premacht@cedar5 combined_files]$ vi cal_admixture_array_2_to_6.sh 
[premacht@cedar5 combined_files]$ vi cal_admixture_array_7_to_12.sh 
[premacht@cedar5 combined_files]$ less cal_admixture_array_2_to_6.sh 

#!/bin/sh
#SBATCH --job-name=bwa_505
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=48:00:00
#SBATCH --mem=16gb
#SBATCH --output=bwa505.%J.out
#SBATCH --error=bwa505.%J.err
#SBATCH --account=def-ben
#SBATCH --array=2-6

#SBATCH --mail-user=premacht@mcmaster.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=ALL

module load nixpkgs/16.09 admixture/1.3.0

# submitting array
 admixture --cv autosomes.bed ${SLURM_ARRAY_TASK_ID} > ./outs_array/chrs_log${SLURM_ARRAY_TASK_ID}.out
```
