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
now with all outputs in the same directory,

1) Download the plotting Rscript 
```bash
wget https://github.com/speciationgenomics/scripts/raw/master/plotADMIXTURE.r
chmod +x plotADMIXTURE.r
```

2) Create a sample list assigning all the samples into populations like following(here I assigned all of them to tropicalis)
This is a tab seperated txt file with sample name in first column and species or population in second column.

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
here,
'autosomes' is the prefix/common name among files
'samples.list' is the file we created above with sample names and populations
'11' is the top most k value to plot(min value is 2 by default. you can change it with flag -m XX)
after '-l' you should list all populations in the text file you created seperated by commas (eg: tropicalis,laevis,pop2,pop3)

This will output plots in tiff format
