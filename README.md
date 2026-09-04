# MAMA2_GP2

## Step 1 - run Regenie in each genotools population 

Use the notebook that is attached to the repository, just change the hardcoded configurations at the top cell. 

## Step 2 - harmonize summary statistics 

Install gwaslab (I used version 4.1.6 without major problems), using the output from regenie do:

```
python harmonize_sumstats.py \
    --input FULL_AAC_GWAS_PD.regenie.gz \
    --output-prefix AAC_harmonized \
    --plot AAC_harmonized.png \
    --title "AAC Regenie GWAS" \
    --threads 4
```
or 
```
python harmonize_sumstats.py --help
```
To see more customizable options.

## Step 3 - subset the sumstats to the hap map 3 variants universe 

This will help to do the MR-MEGA prep easier and faster. 
Since the sumstats could get really big, using polars, the subsetting gets easier

Both the script and the hap map 3 vars are attached to the repository:

```
for i in AAC AFR AJ EAS SAS; do
    python subset_sumstats_by_snp_polars.py --sumstats /config/workspace/ws_files/MAMA_2/jf_release12/*/output/Regenie_GWAS/${i}_R12_imputed/${i}_harmonized.gwaslab.tsv.gz  --variants hpm3snplist.bed --out ${i}_R12_imputed_hpm3snplist.tsv
done

for i in AMR CAS FIN MDE; do
    python subset_sumstats_by_snp_polars.py --sumstats /config/workspace/ws_files/MAMA_2/Release12_PR/output/Regenie_GWAS/${i}_R12_imputed/${i}_harmonized.gwaslab.tsv.gz  --variants hpm3snplist.bed --out ${i}_R12_imputed_hpm3snplist.tsv
done
```

## Step 3 - run the first round of MR-MEGA with the max amount of PCs available 

Important note: for MR-MEGA, the more "granular" you have your data the better. Example, if you have EUR sumstas from GP2, 23adnME, UKB, etc. it is better to input them individually. 

In their paper they say that the max amount of PCs available is equal or less than the number of studies minus two. 
Aka, if we run the thing for 9 pops (the non-EUR ones), the max number of PCs is 7. 

In reality, based on previous experiences with the tool, this rule of thumb is a bit different when you have a lower number of studies. For example, when I run it with 4 studies, the maximum number of PCs to include for the tool to run is actually 1. Since no PCs at all is just basically a fixed effects meta-analysis. I emailed the authors of the tool and they confirmed as well that the minimum amount of studies to run MR-MEGA is 4. Otherwise just stick with fixed and or random effects models. 

The command line to run MR-MEGA is the following:

```
MR-MEGA -i MR-MEGA_input.txt --pc 7 -o MAMA_gp2Only_nonEUR --name_marker SNPID --name_chr CHR --name_pos POS
```
Where the file mR-MEGA_input.txt looks like this (the paths to where your harmonized sumstats are located):
```
cat MR-MEGA_input.txt
/home/jupyter/workspace/ws_files/MAMA_2/jf_release12/aac/output/Regenie_GWAS/AAC_R12_imputed/AAC_harmonized.gwaslab.tsv
/home/jupyter/workspace/ws_files/MAMA_2/jf_release12/afr/output/Regenie_GWAS/AFR_R12_imputed/AFR_harmonized.gwaslab.tsv
/home/jupyter/workspace/ws_files/MAMA_2/jf_release12/aj/output/Regenie_GWAS/AJ_R12_imputed/AJ_harmonized.gwaslab.tsv
/home/jupyter/workspace/ws_files/MAMA_2/Release12_PR/output/Regenie_GWAS/AMR_R12_imputed/AMR_harmonized.gwaslab.tsv
/home/jupyter/workspace/ws_files/MAMA_2/Release12_PR/output/Regenie_GWAS/CAS_R12_imputed/CAS_harmonized.gwaslab.tsv
/home/jupyter/workspace/ws_files/MAMA_2/jf_release12/eas/output/Regenie_GWAS/EAS_R12_imputed/EAS_harmonized.gwaslab.tsv
/home/jupyter/workspace/ws_files/MAMA_2/Release12_PR/output/Regenie_GWAS/FIN_R12_imputed/FIN_harmonized.gwaslab.tsv
/home/jupyter/workspace/ws_files/MAMA_2/Release12_PR/output/Regenie_GWAS/MDE_R12_imputed/MDE_harmonized.gwaslab.tsv
/home/jupyter/workspace/ws_files/MAMA_2/jf_release12/sas/output/Regenie_GWAS/SAS_R12_imputed/SAS_harmonized.gwaslab.tsv
```

## Step 4 - plot the PCs of the first round of MR-MEGA and decide how many PCs to include

After the run is done you can parse the log file to generate the input to the plotting notebook:

```
awk '
BEGIN {OFS="\t"}
/^PCs[[:space:]]+PC0/ {in_pc=1; next}
in_pc && /\.tsv/ {
    sub(/^.*\//, "", $1)
    $1=$1
    print
}
' MAMA_gp2Only_nonEUR.log > cs_MAMA_2026.tsv
```
After this, in the same directory run the plotPCs_2026.ipynb notebook to get the plot and decide how many PCs to include

<img width="5972" height="1465" alt="PCs_all (1)" src="https://github.com/user-attachments/assets/3d2591bb-5ea4-40b6-ae68-08426ce4726a" />

PC1 + PC2 separates AFR, AAC, and EAS from others. AAC is kind of "mid-point" between AFR and others which is a good sign.
PC3 + PC4 - AMR + SAS
PC5 + PC6 - MDE, AJ, FIN, and to a lesser degree CAS.
PC7 really pulls CAS away from the rest

So, for the moment, and waiting for the EUR pop to run we are moving froward with 6

## Step 5 - run the actual MR-MEGA with the number of PCs needed 

```
MR-MEGA -i MR-MEGA_input.txt --pc 6 -o MAMA_gp2Only_nonEUR --name_marker SNPID --name_chr CHR --name_pos POS 2>&1 | tee MR-MEGA_gp2Only_nonEUR.logFile
```
To parse and plot the sumstats use again gwaslab:

```
python harmonize_mrmega.py \
    --input MAMA_gp2Only_nonEUR.result \
    --output-prefix mama_gp2Only_non-eur \
    --plot mama_gp2Only_non-eur.png \
    --title "MR-MEGA GP2 only non-EUR" \
    --threads 4
```
To check if there are any novel loci (the table with the known loci is attached to this repo, i got it from the significant associations in GWAS catalog for PD and the same from NDKP portal):

```
python get_novel_loci.py \
    --input mama_gp2Only_non-eur.gwaslab.tsv.gz \
    --known known_lociPD.tsv \
    --output-prefix mama_gp2Only_non-eur
```

## Step 6 - run random effects 

The logic is pretty much the same as for MR-MEGA but now using GWAMA (from the same creators of MR-MEGA), that allows you to run random and fixed effects meta-analyses. 

```
GWAMA --filelist gwama_input_list.txt \
  --name_marker SNPID --name_ea EA --name_nea NEA --name_eaf EAF \
  --name_or OR --name_or_95l OR_95L --name_or_95u OR_95U \
  --output gwama_random \
  --random
```
the gwama_input_list.txt is the same file we used before with the paths for the input sumstats like in MR-MEGA but with a different name here. 

To parse and plot the sumstats use gwaslab (the script is a bit different since some column names change compared to MR-MEGA):

```
python harmonize_gwama.py \
    --input gwama_random.out \
    --output-prefix gwmam_random_gp2Only_non-eur \
    --plot manhattan.png \
    --title "gp2 non-eur random effects meta" \
    --threads 4
```
And to find the potential "novel loci" use the same script as for MR-MEGA

## Step 7 - generate the miami plot with the novel loci included

For this the script doesn't have argparse (my bad), but is relatively easy to hard code the configurations for it (miami_plot.py). 

