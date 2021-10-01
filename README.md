# ABBA_BABA_2021

I received seperate VCF files for each of the chromosomes of Xenopus laevis. Going to filter them and do ABBA BABA analysis

working on cedar
```
/scratch/premacht/ABBA_BABA_final
```
made three directories for three data sets and copied vcf files there

```bash
mkdir filtered_XLXGXM_radseq
mkdir nonfiltered_XLXGXM_radseq
mkdir XL
```
then copied the folder witgh VCF files in respective directories

To increase efficiency, I renamed chr 9_10 as chr 9 so I can submit job arrays
```bash
find . -type f -name '*9_10*' | while read FILE ; do     newfile="$(echo ${FILE} |sed -e 's/9_10/9/g')" ;     mv "${FILE}" "${newfile}" ; done
```

In the directory for each dataset (showing XL here), created directories
```bash
mkdir sample_filtered_VCF
mkdir scripts
```
Load vcftools and list samples
```bash 
module load StdEnv/2020  gcc/9.3.0
module load bcftools/1.11
bcftools query -l DB_new_chr9_10S_out.vcf_filtered.vcf.gz_filtered_removed.vcf
```

inside scripts

writing this to remove messy sample

```
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

# for S only
for i in ../XL_vcf_files/DB_new_chr${SLURM_ARRAY_TASK_ID}S_out.vcf_filtered.vcf.gz_filtered_removed.vcf; do vcftools --remove-indv BJE3632_Niewou_CATAAGT_cuttrim_sorted.bam --vcf ${i} --out ../sample_filtered_VCF/${i#../filtered_VCFs/}_s_rmved --recode;done

#for L only
for i in ../XL_vcf_files/DB_new_chr${SLURM_ARRAY_TASK_ID}L_out.vcf_filtered.vcf.gz_filtered_removed.vcf; do vcftools --remove-indv BJE3632_Niewou_CATAAGT_cuttrim_sorted.bam --vcf ${i} --out ../sample_filtered_VCF/${i#../filtered_VCFs/}_s_rmved --recode;done

```


Create a reference genome folder and download reference genome in that folder
```
mkdir reference_genome
wget http://ftp.xenbase.org/pub/Genomics/JGI/Xenla9.2/XENLA_9.2_genome.fa.gz
mv XENLA_9.2_genome.fa.gz reference_genome/
```
Unzip and index reference genome
```bash
cd reference_genome
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

make directories
```bash 
mkdir scripts
mkdir combined_files
```
and this only if other filtering steps are not enough ***********************************************
```bash
mkdir filtered_again_VCFs
mkdir filtered_VCFs
mkdir filtered_thinned_VCFs
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

module load StdEnv/2020
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
extra filtering ends here ***************************************************************************

Run this on admixture folder with final filtered vcf files to make seperate combined files for seperate genomes
making directories for all chromosome data, l only, s only

MAKE SURE TO CHANGE 'combined_files' TO WHATEVER THE DIRECTORY WITH FINAL FILTERED VCFS
```bash
cd combined_files
module load StdEnv/2020  intel/2020.1.217 bcftools/1.11
mkdir all
bcftools concat -o all/autosomes.vcf $(ls ../positions_excluded/*.vcf | tr "\n" " ")
mkdir l_only
bcftools concat -o l_only/autosomes.vcf $(ls ../positions_excluded/*L*.vcf | tr "\n" " ")
mkdir s_only
bcftools concat -o s_only/autosomes.vcf $(ls ../positions_excluded/*S*.vcf | tr "\n" " ")
```


compress and index(inside combined_files)
convert to geno format using plink , make a bed file and remove any SNP with no data for autosomes
do all these for l_only and s_only seperately
and we need to change the chr names in the .bim file because these cause problems for admixture:
create a directory to collect outputs from next step
```bash
for i in all l_only s_only; do
module load StdEnv/2020  intel/2020.1.217 bcftools/1.11
cd ${i}
bgzip -c autosomes.vcf > autosomes.vcf.gz
tabix -p vcf autosomes.vcf.gz
module load nixpkgs/16.09  intel/2016.4 plink/1.9b_5.2-x86_64
plink --vcf ./autosomes.vcf.gz --make-bed --geno 0.999 --out ./autosomes --allow-extra-chr --const-fid
awk -v OFS='\t' '{$1=0;print $0}' autosomes.bim > autosomes.bim.tmp
mv autosomes.bim.tmp autosomes.bim
cd .. ; done
```
now save two job arrays to cal admixture(2:6 and 7:12-run array numbers acccordingly) in scripts folder changing array numbers

as
cal_admix_2to6.sh


```bash
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
 admixture --cv autosomes.bed ${SLURM_ARRAY_TASK_ID} > autosomeslog${SLURM_ARRAY_TASK_ID}.out
```
and
cal_admix_7to13.sh
```bash
#!/bin/sh
#SBATCH --job-name=bwa_505
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=48:00:00
#SBATCH --mem=16gb
#SBATCH --output=bwa505.%J.out
#SBATCH --error=bwa505.%J.err
#SBATCH --account=def-ben
#SBATCH --array=7-13

#SBATCH --mail-user=premacht@mcmaster.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=ALL

module load nixpkgs/16.09 admixture/1.3.0

# submitting array
 admixture --cv autosomes.bed ${SLURM_ARRAY_TASK_ID} > autosomeslog${SLURM_ARRAY_TASK_ID}.out
 ```
now run both of them in all three genome types(all,l and s) by pasting this inside directory 'combined_files'

```bash
for i in all l_only s_only; do  
cd ${i} 
cp ../../scripts/cal_admix_2to6.sh . 
cp ../../scripts/cal_admix_7to13.sh . 
sbatch cal_admix_2to6.sh 
sbatch cal_admix_7to13.sh
cd ..; done

```

now with all outputs in the same directory,

1) Download the plotting Rscript - NOT DOING THIS HERE. USED CUSTOMIZED PLOTS INSTEAD AS IN NEXT STEP
```bash
wget https://github.com/speciationgenomics/scripts/raw/master/plotADMIXTURE.r
chmod +x plotADMIXTURE.r
```

2) Create a sample list assigning all the samples into populations like following(here I assigned all of them to tropicalis)

 you can create a tab seperated txt file with sample name in first column and species or population in second column using excel.

```txt
2014_Inhaca_10_Inhaca_ATATGT_cuttrim_sorted.bam tropicalis
2014_Inhaca_150_Inhaca_ATCGTA_cuttrim_sorted.bam        tropicalis
2014_Inhaca_152_Inhaca_CATCGT_cuttrim_sorted.bam        tropicalis
2014_Inhaca_24_Inhaca_CGCGGT_cuttrim_sorted.bam tropicalis
2014_Inhaca_38_Inhaca_CTATTA_cuttrim_sorted.bam tropicalis
2014_Inhaca_52_Inhaca_GCCAGT_cuttrim_sorted.bam tropicalis
2014_Inhaca_65_Inhaca_GGAAGA_cuttrim_sorted.bam tropicalis
946_Draken_TCGTT_cuttrim_sorted.bam     tropicalis
993_Draken_GGTTGT_cuttrim_sorted.bam    tropicalis
BJE2897_Lendu_TTCCTGGA_cuttrim_sorted.bam       tropicalis
BJE3252_Cameroon_TATCGGGA_cuttrim_sorted.bam    tropicalis
BJE3508_DeDorn_ATTGA_cuttrim_sorted.bam tropicalis
BJE3509_DeDorn_CATCT_cuttrim_sorted.bam tropicalis
BJE3510_DeDorn_CCTAG_cuttrim_sorted.bam tropicalis
BJE3511_DeDorn_GAGGA_cuttrim_sorted.bam tropicalis
BJE3512_DeDorn_GGAAG_cuttrim_sorted.bam tropicalis
BJE3513_DeDorn_GTCAA_cuttrim_sorted.bam tropicalis
BJE3514_DeDorn_TAATA_cuttrim_sorted.bam tropicalis
BJE3515_DeDorn_TACAT_cuttrim_sorted.bam tropicalis
BJE3525_Laigns_GAATTCA_cuttrim_sorted.bam       tropicalis
BJE3526_Laigns_GAACTTG_cuttrim_sorted.bam       tropicalis
BJE3527_Laigns_GGACCTA_cuttrim_sorted.bam       tropicalis
BJE3528_Laigns_GTCGATT_cuttrim_sorted.bam       tropicalis
BJE3529_Laigns_AACGCCT_cuttrim_sorted.bam       tropicalis
BJE3530_Laigns_AATATGG_cuttrim_sorted.bam       tropicalis
BJE3531_Laigns_ACGTGTT_cuttrim_sorted.bam       tropicalis
BJE3532_Laigns_ATTAATT_cuttrim_sorted.bam       tropicalis
BJE3533_Laigns_ATTGGAT_cuttrim_sorted.bam       tropicalis
BJE3534_BW_CTCG_cuttrim_sorted.bam      tropicalis
BJE3535_BW_TGCA_cuttrim_sorted.bam      tropicalis
BJE3536_BW_ACTA_cuttrim_sorted.bam      tropicalis
BJE3537_BW_CAGA_cuttrim_sorted.bam      tropicalis
BJE3538_BW_AACT_cuttrim_sorted.bam      tropicalis
BJE3539_BW_GCGT_cuttrim_sorted.bam      tropicalis
BJE3541_BW_CGAT_cuttrim_sorted.bam      tropicalis
BJE3542_BW_GTAA_cuttrim_sorted.bam      tropicalis
BJE3543_BW_AGCG_cuttrim_sorted.bam      tropicalis
BJE3544_BW_GATG_cuttrim_sorted.bam      tropicalis
BJE3545_BW_TCAG_cuttrim_sorted.bam      tropicalis
BJE3546_BW_TGCGA_cuttrim_sorted.bam     tropicalis
BJE3547_GRNP_TAGGAA_cuttrim_sorted.bam  tropicalis
BJE3548_GRNP_GCTCTA_cuttrim_sorted.bam  tropicalis
BJE3549_GRNP_CCACAA_cuttrim_sorted.bam  tropicalis
BJE3550_GRNP_CTTCCA_cuttrim_sorted.bam  tropicalis
BJE3551_GRNP_GAGATA_cuttrim_sorted.bam  tropicalis
BJE3552_GRNP_ATGCCT_cuttrim_sorted.bam  tropicalis
BJE3553_GRNP_AGTGGA_cuttrim_sorted.bam  tropicalis
BJE3554_GRNP_ACCTAA_cuttrim_sorted.bam  tropicalis
BJE3573_VicW_CGCGGAGA_cuttrim_sorted.bam        tropicalis
BJE3574_VicW_CGTGTGGT_cuttrim_sorted.bam        tropicalis
BJE3575_Kimber_GTACTT_cuttrim_sorted.bam        tropicalis
BJE3576_Kimber_GTTGAA_cuttrim_sorted.bam        tropicalis
BJE3577_Kimber_TAACGA_cuttrim_sorted.bam        tropicalis
BJE3578_Kimber_TGGCTA_cuttrim_sorted.bam        tropicalis
BJE3579_Kimber_TATTTTT_cuttrim_sorted.bam       tropicalis
BJE3580_Kimber_CTTGCTT_cuttrim_sorted.bam       tropicalis
BJE3581_Kimber_ATGAAAG_cuttrim_sorted.bam       tropicalis
BJE3582_Kimber_AAAAGTT_cuttrim_sorted.bam       tropicalis
BJE3632_Niewou_CATAAGT_cuttrim_sorted.bam       tropicalis
BJE3633_Niewou_CGCTGAT_cuttrim_sorted.bam       tropicalis
BJE3640_Niewou_CGGTAGA_cuttrim_sorted.bam       tropicalis
BJE3641_Niewou_CTACGGA_cuttrim_sorted.bam       tropicalis
BJE3642_Niewou_GCGGAAT_cuttrim_sorted.bam       tropicalis
BJE3644_Niewou_TAGCGGA_cuttrim_sorted.bam       tropicalis
BJE3645_Niewou_TCGAAGA_cuttrim_sorted.bam       tropicalis
BJE3647_Niewou_TCTGTGA_cuttrim_sorted.bam       tropicalis
BJE3654_ThreeSis_TGCTGGA_cuttrim_sorted.bam     tropicalis
BJE3655_ThreeSis_ACGACTAG_cuttrim_sorted.bam    tropicalis
BJE3656_ThreeSis_TAGCATGG_cuttrim_sorted.bam    tropicalis
BJE3657_ThreeSis_TAGGCCAT_cuttrim_sorted.bam    tropicalis
BJE3658_ThreeSis_TGCAAGGA_cuttrim_sorted.bam    tropicalis
BJE3659_ThreeSis_TGGTACGT_cuttrim_sorted.bam    tropicalis
BJE3660_ThreeSis_TCTCAGTG_cuttrim_sorted.bam    tropicalis
BJE3661_ThreeSis_CGCGATAT_cuttrim_sorted.bam    tropicalis
BJE3662_ThreeSis_CGCCTTAT_cuttrim_sorted.bam    tropicalis
BJE3663_ThreeSis_AACCGAGA_cuttrim_sorted.bam    tropicalis
BJE3664_ThreeSis_ACAGGGA_cuttrim_sorted.bam     tropicalis
BJE3665_ThreeSis_ACGTGGTA_cuttrim_sorted.bam    tropicalis
BJE3666_ThreeSis_CCATGGGT_cuttrim_sorted.bam    tropicalis
BJE3667_Citrus_CGCTT_cuttrim_sorted.bam tropicalis
BJE3668_Citrus_TCACG_cuttrim_sorted.bam tropicalis
BJE3669_Citrus_CTAGG_cuttrim_sorted.bam tropicalis
BJE3670_Citrus_ACAAA_cuttrim_sorted.bam tropicalis
BJE3671_Citrus_TTCTG_cuttrim_sorted.bam tropicalis
BJE3672_Citrus_AGCCG_cuttrim_sorted.bam tropicalis
BJE3673_Citrus_GTATT_cuttrim_sorted.bam tropicalis
BJE3674_Citrus_CTGTA_cuttrim_sorted.bam tropicalis
BJE3675_Citrus_ACCGT_cuttrim_sorted.bam tropicalis
BJE3676_Citrus_GCTTA_cuttrim_sorted.bam tropicalis
BJE3677_Citrus_GGTGT_cuttrim_sorted.bam tropicalis
BJE3678_Citrus_AGGAT_cuttrim_sorted.bam tropicalis
JM_no_label1_Draken_CCACGT_cuttrim_sorted.bam   tropicalis
JM_no_label2_Draken_TTCAGA_cuttrim_sorted.bam   tropicalis
RT5_Botsw_GGATTGGT_cuttrim_sorted.bam   tropicalis
Vred_8_Vred_GCTGTGGA_cuttrim_sorted.bam tropicalis
amnh17260_Nigeria_GTGAGGGT_cuttrim_sorted.bam   tropicalis

```
Now you can plot admixture for different k values with following command after loading the r module

```bash
module load gcc/9.3.0 r/4.1.0
```
run
```bash
Rscript plotADMIXTURE.r -p autosomes -i samples.list -k 11 -l tropicalis
```
If this asks you to install r packages, do that creating a personal library after opening R

here,
'autosomes' is the prefix/common name among files
'samples.list' is the file we created above with sample names and populations
'11' is the top most k value to plot(min value is 2 by default. you can change it with flag -m XX)
after '-l' you should list all populations in the text file you created seperated by commas (eg: tropicalis,laevis,pop2,pop3)

This will output plots in tiff format
# Customized script

I wanted a bunch of customization in the plots given by the default R script. So I wrote the following to be used in local machine R studio.

I wrote this in a way you can just copy this script to your working folder and run. This creates a directory named plot_outs and save the resulting plots in that.

For this to work,

1)put all outputs from previous sections including
   all files ending with P or Q(eg; autosomes.2.p)
   all files created from previous scripts with different formats like .bed .bim .fam ....etc.(eg; autosomes.bed)
   .log and all .out files
   
   Inside a directory named 'all_outputs' which is in your working directory.
   
2) sample list created with assigned populations as shown below(should be named as samples.list and the r script should be in the current working directory(not inside all_outputs)

samples.list
```txt
2014_Inhaca_10_Inhaca_ATATGT_cuttrim_sorted.bam	Inhaca
2014_Inhaca_150_Inhaca_ATCGTA_cuttrim_sorted.bam	Inhaca
2014_Inhaca_152_Inhaca_CATCGT_cuttrim_sorted.bam	Inhaca
2014_Inhaca_24_Inhaca_CGCGGT_cuttrim_sorted.bam	Inhaca
2014_Inhaca_38_Inhaca_CTATTA_cuttrim_sorted.bam	Inhaca
2014_Inhaca_52_Inhaca_GCCAGT_cuttrim_sorted.bam	Inhaca
2014_Inhaca_65_Inhaca_GGAAGA_cuttrim_sorted.bam	Inhaca
946_Draken_TCGTT_cuttrim_sorted.bam	Draken
993_Draken_GGTTGT_cuttrim_sorted.bam	Draken
BJE2897_Lendu_TTCCTGGA_cuttrim_sorted.bam	Other
BJE3252_Cameroon_TATCGGGA_cuttrim_sorted.bam	Other
BJE3508_DeDorn_ATTGA_cuttrim_sorted.bam	DeDorn
BJE3509_DeDorn_CATCT_cuttrim_sorted.bam	DeDorn
BJE3510_DeDorn_CCTAG_cuttrim_sorted.bam	DeDorn
BJE3511_DeDorn_GAGGA_cuttrim_sorted.bam	DeDorn
BJE3512_DeDorn_GGAAG_cuttrim_sorted.bam	DeDorn
BJE3513_DeDorn_GTCAA_cuttrim_sorted.bam	DeDorn
BJE3514_DeDorn_TAATA_cuttrim_sorted.bam	DeDorn
BJE3515_DeDorn_TACAT_cuttrim_sorted.bam	DeDorn
BJE3525_Laigns_GAATTCA_cuttrim_sorted.bam	Laigns
BJE3526_Laigns_GAACTTG_cuttrim_sorted.bam	Laigns
BJE3527_Laigns_GGACCTA_cuttrim_sorted.bam	Laigns
BJE3528_Laigns_GTCGATT_cuttrim_sorted.bam	Laigns
BJE3529_Laigns_AACGCCT_cuttrim_sorted.bam	Laigns
BJE3530_Laigns_AATATGG_cuttrim_sorted.bam	Laigns
BJE3531_Laigns_ACGTGTT_cuttrim_sorted.bam	Laigns
BJE3532_Laigns_ATTAATT_cuttrim_sorted.bam	Laigns
BJE3533_Laigns_ATTGGAT_cuttrim_sorted.bam	Laigns
BJE3534_BW_CTCG_cuttrim_sorted.bam	BW
BJE3535_BW_TGCA_cuttrim_sorted.bam	BW
BJE3536_BW_ACTA_cuttrim_sorted.bam	BW
BJE3537_BW_CAGA_cuttrim_sorted.bam	BW
BJE3538_BW_AACT_cuttrim_sorted.bam	BW
BJE3539_BW_GCGT_cuttrim_sorted.bam	BW
BJE3541_BW_CGAT_cuttrim_sorted.bam	BW
BJE3542_BW_GTAA_cuttrim_sorted.bam	BW
BJE3543_BW_AGCG_cuttrim_sorted.bam	BW
BJE3544_BW_GATG_cuttrim_sorted.bam	BW
BJE3545_BW_TCAG_cuttrim_sorted.bam	BW
BJE3546_BW_TGCGA_cuttrim_sorted.bam	BW
BJE3547_GRNP_TAGGAA_cuttrim_sorted.bam	GRNP
BJE3548_GRNP_GCTCTA_cuttrim_sorted.bam	GRNP
BJE3549_GRNP_CCACAA_cuttrim_sorted.bam	GRNP
BJE3550_GRNP_CTTCCA_cuttrim_sorted.bam	GRNP
BJE3551_GRNP_GAGATA_cuttrim_sorted.bam	GRNP
BJE3552_GRNP_ATGCCT_cuttrim_sorted.bam	GRNP
BJE3553_GRNP_AGTGGA_cuttrim_sorted.bam	GRNP
BJE3554_GRNP_ACCTAA_cuttrim_sorted.bam	GRNP
BJE3573_VicW_CGCGGAGA_cuttrim_sorted.bam	VicW
BJE3574_VicW_CGTGTGGT_cuttrim_sorted.bam	VicW
BJE3575_Kimber_GTACTT_cuttrim_sorted.bam	Kimber
BJE3576_Kimber_GTTGAA_cuttrim_sorted.bam	Kimber
BJE3577_Kimber_TAACGA_cuttrim_sorted.bam	Kimber
BJE3578_Kimber_TGGCTA_cuttrim_sorted.bam	Kimber
BJE3579_Kimber_TATTTTT_cuttrim_sorted.bam	Kimber
BJE3580_Kimber_CTTGCTT_cuttrim_sorted.bam	Kimber
BJE3581_Kimber_ATGAAAG_cuttrim_sorted.bam	Kimber
BJE3582_Kimber_AAAAGTT_cuttrim_sorted.bam	Kimber
BJE3632_Niewou_CATAAGT_cuttrim_sorted.bam	Niewou
BJE3633_Niewou_CGCTGAT_cuttrim_sorted.bam	Niewou
BJE3640_Niewou_CGGTAGA_cuttrim_sorted.bam	Niewou
BJE3641_Niewou_CTACGGA_cuttrim_sorted.bam	Niewou
BJE3642_Niewou_GCGGAAT_cuttrim_sorted.bam	Niewou
BJE3644_Niewou_TAGCGGA_cuttrim_sorted.bam	Niewou
BJE3645_Niewou_TCGAAGA_cuttrim_sorted.bam	Niewou
BJE3647_Niewou_TCTGTGA_cuttrim_sorted.bam	Niewou
BJE3654_ThreeSis_TGCTGGA_cuttrim_sorted.bam	Threesis
BJE3655_ThreeSis_ACGACTAG_cuttrim_sorted.bam	Threesis
BJE3656_ThreeSis_TAGCATGG_cuttrim_sorted.bam	Threesis
BJE3657_ThreeSis_TAGGCCAT_cuttrim_sorted.bam	Threesis
BJE3658_ThreeSis_TGCAAGGA_cuttrim_sorted.bam	Threesis
BJE3659_ThreeSis_TGGTACGT_cuttrim_sorted.bam	Threesis
BJE3660_ThreeSis_TCTCAGTG_cuttrim_sorted.bam	Threesis
BJE3661_ThreeSis_CGCGATAT_cuttrim_sorted.bam	Threesis
BJE3662_ThreeSis_CGCCTTAT_cuttrim_sorted.bam	Threesis
BJE3663_ThreeSis_AACCGAGA_cuttrim_sorted.bam	Threesis
BJE3664_ThreeSis_ACAGGGA_cuttrim_sorted.bam	Threesis
BJE3665_ThreeSis_ACGTGGTA_cuttrim_sorted.bam	Threesis
BJE3666_ThreeSis_CCATGGGT_cuttrim_sorted.bam	Threesis
BJE3667_Citrus_CGCTT_cuttrim_sorted.bam	Citrus
BJE3668_Citrus_TCACG_cuttrim_sorted.bam	Citrus
BJE3669_Citrus_CTAGG_cuttrim_sorted.bam	Citrus
BJE3670_Citrus_ACAAA_cuttrim_sorted.bam	Citrus
BJE3671_Citrus_TTCTG_cuttrim_sorted.bam	Citrus
BJE3672_Citrus_AGCCG_cuttrim_sorted.bam	Citrus
BJE3673_Citrus_GTATT_cuttrim_sorted.bam	Citrus
BJE3674_Citrus_CTGTA_cuttrim_sorted.bam	Citrus
BJE3675_Citrus_ACCGT_cuttrim_sorted.bam	Citrus
BJE3676_Citrus_GCTTA_cuttrim_sorted.bam	Citrus
BJE3677_Citrus_GGTGT_cuttrim_sorted.bam	Citrus
BJE3678_Citrus_AGGAT_cuttrim_sorted.bam	Citrus
JM_no_label1_Draken_CCACGT_cuttrim_sorted.bam	Draken
JM_no_label2_Draken_TTCAGA_cuttrim_sorted.bam	Draken
RT5_Botsw_GGATTGGT_cuttrim_sorted.bam	Other
Vred_8_Vred_GCTGTGGA_cuttrim_sorted.bam	Vred
amnh17260_Nigeria_GTGAGGGT_cuttrim_sorted.bam	Other

```

plotADMIXTURE.r

```r
#!/usr/bin/Rscript

# Usage: plotADMIXTURE.r -p <prefix> -i <info file, 2-column file with ind name and population/species name> 
#                        -k <max K value> -l <comma-separated list of populations/species in the order to be plotted>
# This R script makes barplots for K=2 and all other K values until max K (specified with -k). It labels the individuals 
# and splits them into populations or species according to the individual and population/species names in the 2-column file specified with -i.
# The order of populations/species follows the list of populations/species given with -l.
# Usage example: plotADMIXTURE.r -p fileXY -i file.ind.pop.txt -k 4 -pop pop1,pop2,pop3
# In this example, the script would use the files fileXY.2.Q, fileXY.3.Q, fileXY.4.Q to make barplots for the three populations.
# file.ind.pop.txt should contain one line for each individual in the same order as in the admixture files e.g.
# ind1 pop1
# ind2 pop1
# ind3 pop2
# ind4 pop3

# Author: Joana Meier, September 2019

# I am customizing this script. Making it easy to use with R studio-TS
# ***keeping this script in current directory for clarity. But setting the working directory to the directory with all output files***


#set working directory to all_outputs inside the current path
setwd(paste(dirname(rstudioapi::getSourceEditorContext()$path),"/all_outputs",sep=""))

#set default values here so it can run on R studio without parsing options
#These will come into action if you do not define these options like you do in bash

#change prefix here
p_input<-"autosomes"
#change sample list file here
i_input<-paste(dirname(rstudioapi::getSourceEditorContext()$path),"/samples.list",sep="")
# change maximum k value to plot here
k_input<-11
# change minimum k value to plot here
m_input<-2
#add a list of populations seperated by commas here. This should be exactly similiar to the populations in your sample list file. plots will be created according to this population order
# you will have to change this every time you edit sample list
l_input<-"GRNP,DeDorn,Citrus,Vred,Laigns,BW,Threesis,VicW,Kimber,Draken,Inhaca,Niewou,Other"
#add the location and file name for the plots here. I am setting this to the directory I am creating in the next line
o_input<-paste(dirname(rstudioapi::getSourceEditorContext()$path),"/plot_outs/plot",sep="")
  
#create a directory for plot output if it doesn't already exist in the directory with the script
dir.create(file.path(dirname(rstudioapi::getSourceEditorContext()$path), "plot_outs"))

# Read in the arguments
library("optparse")
option_list = list(
  make_option(c("-p", "--prefix"), type="character", default=p_input, 
              help="prefix name (with path if not in the current directory)", metavar="character"),
  make_option(c("-i", "--infofile"), type="character", default=i_input, 
              help="info text file containing for each individual the population/species information", metavar="character"),
  make_option(c("-k", "--maxK"), type="integer", default=k_input, 
              help="maximum K value", metavar="integer"),
  make_option(c("-m", "--minK"), type="integer", default=m_input, 
              help="minimum K value", metavar="integer"),
  make_option(c("-l", "--populations"), type="character", default=l_input, 
              help="comma-separated list of populations/species in the order to be plotted", metavar="character"),
  make_option(c("-o", "--outPrefix"), type="character", default=o_input, 
              help="output prefix (default: name provided with prefix)", metavar="character")
) 
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Check that all required arguments are provided
if (is.null(opt$prefix)){
  print_help(opt_parser)
  stop("Please provide the prefix", call.=FALSE)
}else if (is.null(opt$infofile)){
  print_help(opt_parser)
  stop("Please provide the info file", call.=FALSE)
}else if (is.null(opt$maxK)){
  print_help(opt_parser)
  stop("Please provide the maximum K value to plot", call.=FALSE)
}else if (is.null(opt$populations)){
  print_help(opt_parser)
  stop("Please provide a comma-separated list of populations/species", call.=FALSE)
}

# If no output prefix is given, use the input prefix
if(opt$outPrefix=="default") opt$outPrefix=opt$prefix

# Assign the first argument to prefix
prefix=opt$prefix

# Get individual names in the correct order
labels<-read.table(opt$infofile)

# Name the columns
names(labels)<-c("ind","pop")

# Add a column with population indices to order the barplots
# Use the order of populations provided as the fourth argument (list separated by commas)
labels$n<-factor(labels$pop,levels=unlist(strsplit(opt$populations,",")))
levels(labels$n)<-c(1:length(levels(labels$n)))
labels$n<-as.integer(as.character(labels$n))

# read in the different admixture output files
minK=opt$minK
maxK=opt$maxK
tbl<-lapply(minK:maxK, function(x) read.table(paste0(prefix,".",x,".Q")))

# Prepare spaces to separate the populations/species
rep<-as.vector(table(labels$n))
spaces<-0
for(i in 1:length(rep)){spaces=c(spaces,rep(0,rep[i]-1),0.15)}
spaces<-spaces[-length(spaces)]

# Plot the cluster assignments as a single bar for each individual for each K as a separate row
tiff(file=paste0(opt$outPrefix,".tiff"),width = 2000, height = 1300,res=200)
 par(mfrow=c(maxK-1,1),mar=c(0,1,0,0),oma=c(2,1,9,1),mgp=c(0,0.2,0),xaxs="i",cex.lab=1.2,cex.axis=0.8)
 
#Create a list for plots
 plot_list<-list()
 #assign colors to populations at once
 
 col1<-"red"
 col2<-"green"
 col3<-"orange"
 col4<-"skyblue"
 col5<-"darkblue"
 col6<-"yellow"
 col7<-"pink"
 col8<-"purple"
 col9<-"grey"
 col10<-"black"
 col11<-"forestgreen"
 col12<-"brown"
 
 #create color palettes for each k value
 col_palette_k2<-c(col1,col2)
 col_palette_k3<-c(col3,col2,col1)
 col_palette_k4<-c(col1,col2,col3,col4)
 col_palette_k5<-c(col1,col3,col2,col5,col4)
 col_palette_k6<-c(col1,col2,col6,col5,col4,col3)
 col_palette_k7<-c(col1,col7,col5,col6,col4,col3,col2)
 col_palette_k8<-c(col5,col8,col1,col7,col6,col3,col4,col2)
 col_palette_k9<-c(col8,col6,col7,col2,col5,col4,col3,col1,col9)
 col_palette_k10<-c(col8,col7,col4,col10,col2,col9,col6,col1,col3,col5)
 col_palette_k11<-c(col6,col10,col2,col4,col1,col7,col9,col8,col3,col11,col5)
 col_palette_k12<-c(col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11,col12)
 
 # make a variable with full sample names 
 full_sample_names<-labels$ind[order(labels$n)]
 
 # get simpler names  for those complex sample names
 library(plyr)
 #first renaming the sample names which does not follow the specific format manually**DO NOT USE"_" HERE AS SPLIT WILL USE THIS IN NEXT LINE
 simple_sample_list_changing<-mapvalues(full_sample_names, from = c("Vred_8_Vred_GCTGTGGA_cuttrim_sorted.bam","JM_no_label1_Draken_CCACGT_cuttrim_sorted.bam","JM_no_label2_Draken_TTCAGA_cuttrim_sorted.bam","946_Draken_TCGTT_cuttrim_sorted.bam","993_Draken_GGTTGT_cuttrim_sorted.bam","2014_Inhaca_10_Inhaca_ATATGT_cuttrim_sorted.bam","2014_Inhaca_150_Inhaca_ATCGTA_cuttrim_sorted.bam","2014_Inhaca_152_Inhaca_CATCGT_cuttrim_sorted.bam","2014_Inhaca_24_Inhaca_CGCGGT_cuttrim_sorted.bam","2014_Inhaca_38_Inhaca_CTATTA_cuttrim_sorted.bam","2014_Inhaca_52_Inhaca_GCCAGT_cuttrim_sorted.bam","2014_Inhaca_65_Inhaca_GGAAGA_cuttrim_sorted.bam"), 
                                                             to = c("Vred8","JM1","JM2","Draken946","Draken993","Inhaca10","Inhaca150","Inhaca152","Inhaca24","Inhaca38","Inhaca52","Inhaca65"))
 
 #convert facter list into chrs to rename
 s_list_chr<-as.character(simple_sample_list_changing)
 #then remove the parts after the first"_" from other samples
 shortened_sample_list<-sapply(strsplit(s_list_chr,split = "_"),`[`, 1)
 
 
 # following plots are written in a way you can just copy paste only changing k_val for the different number of 'k's
 #paste whats inside * marks for different k values and then change k_val
 #**********
 # Plot k=2
 # change only here
 k_val<-2
 

 
 
 bp<-barplot(t(as.matrix(tbl[[k_val-1]][order(labels$n),])), col=get(paste("col_palette_k",k_val,sep="")),xaxt="n", border=NA,ylab=paste0("K=",k_val,sep=""),yaxt="n",space=spaces)
 axis(3,at=bp,labels=shortened_sample_list,las=2,tick=F,cex=0.6)
 

 
 
 #***********
 
 #commwnt out this section to automate coloring. Keep commented to manually assign colours for each pop with the next section( you can do it once by selecting lines and ctrl+shift+c)
 # Plot higher K values
#  if(maxK>minK)lapply(2:(maxK-1), function(x) barplot(t(as.matrix(tbl[[x]][order(labels$n),])), col=rainbow(n=x+1),xaxt="n", border=NA,ylab=paste0("K=",x+1),yaxt="n",space=spaces))
#  axis(1,at=c(which(spaces==0.15),bp[length(bp)])-diff(c(1,which(spaces==0.15),bp[length(bp)]))/2,
#      labels=unlist(strsplit(opt$populations,",")))
# dev.off()
 
 #**********
 # Plot k
 # change only here
 k_val<-3
 

 
 
 bp<-barplot(t(as.matrix(tbl[[k_val-1]][order(labels$n),])), col=get(paste("col_palette_k",k_val,sep="")),xaxt="n", border=NA,ylab=paste0("K=",k_val,sep=""),yaxt="n",space=spaces)
 axis(1,at=c(which(spaces==0.15),bp[length(bp)])-diff(c(1,which(spaces==0.15),bp[length(bp)]))/2,
      labels = FALSE)
 
 
 
 
 #***********
 
 #**********
 # Plot k
 # change only here
 k_val<-4
 

 
 
 bp<-barplot(t(as.matrix(tbl[[k_val-1]][order(labels$n),])), col=get(paste("col_palette_k",k_val,sep="")),xaxt="n", border=NA,ylab=paste0("K=",k_val,sep=""),yaxt="n",space=spaces)
 axis(1,at=c(which(spaces==0.15),bp[length(bp)])-diff(c(1,which(spaces==0.15),bp[length(bp)]))/2,
      llabels = FALSE)
 
 
 
 #***********
 
 #**********
 # Plot k
 # change only here
 k_val<-5
 

 
 
 bp<-barplot(t(as.matrix(tbl[[k_val-1]][order(labels$n),])), col=get(paste("col_palette_k",k_val,sep="")),xaxt="n", border=NA,ylab=paste0("K=",k_val,sep=""),yaxt="n",space=spaces)
 axis(1,at=c(which(spaces==0.15),bp[length(bp)])-diff(c(1,which(spaces==0.15),bp[length(bp)]))/2,
      labels = FALSE)
 
 
 
 
 #***********
 #**********
 # Plot k
 # change only here
 k_val<-6
 
 
 
 
 bp<-barplot(t(as.matrix(tbl[[k_val-1]][order(labels$n),])), col=get(paste("col_palette_k",k_val,sep="")),xaxt="n", border=NA,ylab=paste0("K=",k_val,sep=""),yaxt="n",space=spaces)
 axis(1,at=c(which(spaces==0.15),bp[length(bp)])-diff(c(1,which(spaces==0.15),bp[length(bp)]))/2,
      labels = FALSE)
 
 
 
 #***********
 
 #**********
 # Plot k
 # change only here
 k_val<-7
 

 
 
 bp<-barplot(t(as.matrix(tbl[[k_val-1]][order(labels$n),])), col=get(paste("col_palette_k",k_val,sep="")),xaxt="n", border=NA,ylab=paste0("K=",k_val,sep=""),yaxt="n",space=spaces)
 axis(1,at=c(which(spaces==0.15),bp[length(bp)])-diff(c(1,which(spaces==0.15),bp[length(bp)]))/2,
      labels = FALSE)
 
 
 
 #***********
 #**********
 # Plot k
 # change only here
 k_val<-8
 
 
 
 
 bp<-barplot(t(as.matrix(tbl[[k_val-1]][order(labels$n),])), col=get(paste("col_palette_k",k_val,sep="")),xaxt="n", border=NA,ylab=paste0("K=",k_val,sep=""),yaxt="n",space=spaces)
 axis(1,at=c(which(spaces==0.15),bp[length(bp)])-diff(c(1,which(spaces==0.15),bp[length(bp)]))/2,
      labels = FALSE)
 
 
 
 #***********
 #**********
 # Plot k
 # change only here
 k_val<-9
 

 
 
 bp<-barplot(t(as.matrix(tbl[[k_val-1]][order(labels$n),])), col=get(paste("col_palette_k",k_val,sep="")),xaxt="n", border=NA,ylab=paste0("K=",k_val,sep=""),yaxt="n",space=spaces)
 axis(1,at=c(which(spaces==0.15),bp[length(bp)])-diff(c(1,which(spaces==0.15),bp[length(bp)]))/2,
      labels = FALSE)
 
 
 
 #***********
 #**********
 # Plot k
 # change only here
 k_val<-10
 

 
 
 bp<-barplot(t(as.matrix(tbl[[k_val-1]][order(labels$n),])), col=get(paste("col_palette_k",k_val,sep="")),xaxt="n", border=NA,ylab=paste0("K=",k_val,sep=""),yaxt="n",space=spaces)
 axis(1,at=c(which(spaces==0.15),bp[length(bp)])-diff(c(1,which(spaces==0.15),bp[length(bp)]))/2,
      labels = FALSE)
 
 
 
 
 #***********
 #**********
 # Plot k
 # change only here
 k_val<-11
 

 
 
 bp<-barplot(t(as.matrix(tbl[[k_val-1]][order(labels$n),])), col=get(paste("col_palette_k",k_val,sep="")),xaxt="n", border=NA,ylab=paste0("K=",k_val,sep=""),yaxt="n",space=spaces)
 axis(1,at=c(which(spaces==0.15),bp[length(bp)])-diff(c(1,which(spaces==0.15),bp[length(bp)]))/2,
      labels=unlist(strsplit(opt$populations,",")))
 
 
 
 # #***********
 # #**********
 # # Plot k
 # # change only here
 # k_val<-12
 # 
 # 
 # 
 # 
 # bp<-barplot(t(as.matrix(tbl[[k_val-1]][order(labels$n),])), col=get(paste("col_palette_k",k_val,sep="")),xaxt="n", border=NA,ylab=paste0("K=",k_val,sep=""),yaxt="n",space=spaces)
 # axis(1,at=c(which(spaces==0.15),bp[length(bp)])-diff(c(1,which(spaces==0.15),bp[length(bp)]))/2,
 #      labels=unlist(strsplit(opt$populations,",")))
 
 dev.off()
 
 ```



   
   
