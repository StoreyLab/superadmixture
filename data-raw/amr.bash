#!/bin/bash
SOFTWARE_PATH="../../../../reproducible-analysis/utils"
zstd="${SOFTWARE_PATH}/zstd/programs/zstd"
PLINK2="${SOFTWARE_PATH}/PLINK2/plink2"

## clean up the amr folder
if [  -d "amr-rawdata" ]; then
  rm -rf amr-rawdata 
fi 
mkdir -p amr-rawdata && cd $_

## download TGP 
if ! [ -f "all_phases3.pgen" ]; then 
  wget "https://www.dropbox.com/s/j72j6uciq5zuzii/all_hg38.pgen.zst?dl=1"
  mv "all_hg38.pgen.zst?dl=1" "all_phases3.pgen.zst"
  ./$zstd -d all_phases3.pgen.zst
  rm all_phases3.pgen.zst
fi

if ! [ -f "all_phases3.pvar.zst" ]; then
  wget "https://www.dropbox.com/s/vx09262b4k1kszy/all_hg38.pvar.zst?dl=1"
  mv "all_hg38.pvar.zst?dl=1" "all_phases3.pvar.zst"
fi 

if ! [ -f "all_phases3.psam" ]; then
  wget "https://www.dropbox.com/s/2e87z6nc4qexjjm/hg38_corrected.psam?dl=1"
  mv "hg38_corrected.psam?dl=1" "all_phases3.psam"
fi

## assign unique variant IDs to each locus
./$PLINK2 \
    --pfile all_phases3 vzs \
    --set-missing-var-ids '@:#' \
    --allow-extra-chr \
    --make-just-pvar zs \
    --out all_phases3_uniq

mv all_phases3.pvar.zst all_phases3.orig.pvar.zst
mv all_phases3_uniq.pvar.zst all_phases3.pvar.zst
rm all_phases3_uniq.log
rm all_phases3.orig.pvar.zst

## preserve loci that:
##  1. are autosomal, biallelic SNPs (model assumption)
##  2. are variant in the Yoruba samples (code YRI), and
##  3. have unique IDs
grep -E '\#|YRI' all_phases3.psam > all_phases3_YRI.psam

./$PLINK2 \
    --pfile all_phases3 vzs \
    --rm-dup exclude-all \
    --allow-extra-chr \
    --write-snplist zs \
    --out nodups

./$PLINK2 \
    --pfile all_phases3 vzs \
    --keep all_phases3_YRI.psam \
    --extract nodups.snplist.zst \
    --autosome \
    --allow-extra-chr \
    --snps-only just-acgt \
    --max-alleles 2 \
    --keep-founders \
    --mac 1 \
    --write-snplist zs \
    --out YRI

rm YRI.log
rm all_phases3_YRI.psam
rm nodups.log
rm nodups.snplist.zst

## preserve individuals marked as "AMR"
grep -E '\#|AMR' all_phases3.psam > all_phases3_AMR.psam

./$PLINK2 \
    --pfile all_phases3 vzs \
    --keep all_phases3_AMR.psam \
    --extract YRI.snplist.zst \
    --keep-founders \
    --allow-extra-chr \
    --mac 1 \
    --make-bed \
    --out amr
    
    
rm YRI.snplist.zst
rm all_phases3_AMR.psam
rm amr.log

## add subpopulation labels to FAM files.
awk -F" " '{
  if (NR == FNR) population[$1]=$6;
  else print population[$2],$2,$3,$4,$5,$6;
}' all_phases3.psam amr.fam > amr.fam.NEW

mv amr.fam.NEW amr.fam

## preserve loci that:
## 1. have MAF >= 0.01
## 2. are in approximate linkage equilibrium with each other
./$PLINK2 \
    --bfile amr \
    --make-bed \
    --out amr_maf_0.01 \
    --maf 0.01
    
## this command determines the loci to keep or exclude
./$PLINK2 --bfile amr_maf_0.01 --indep-pairwise 1000kb 0.3 --out amr_maf_0.01
  
## this actually filters the data
./$PLINK2 \
    --bfile amr_maf_0.01 \
    --extract amr_maf_0.01.prune.in \
    --make-bed \
    --out amr_final
    
rm amr.{bed,bim,fam}
rm amr_maf_0.01.{bed,bim,fam,log,prune.in,prune.out}
rm amr_final.log
