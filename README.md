# ABBA_BABA_2021

I received seperate VCF files for each of the chromosomes of Xenopus laevis. Going to filter them and do ABBA BABA analysis

Testing codes










******************** Not going to use here *********************************************
```txt
# Select few samples to make a smaller VCF file and make it easier to process for the test run

view samples
```bash
module load bcftools
bcftools query -l DB_chr7S_out.vcf
```
Select few samples and make files out of them
First, create a txt file with the sample list to extract using vi
then, (this may give an error trying to read white space as a file. No need to worry)
```bash
while read i ; do bcftools view -s "$i" DB_chr7S_out.vcf > "$i".vcf ; done < sample_list_to_extract.txt
```
```
