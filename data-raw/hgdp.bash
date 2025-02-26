#!/bin/bash
SOFTWARE_PATH="../../../../reproducible-analysis/utils"
zstd="${SOFTWARE_PATH}/zstd/programs/zstd"
PLINK2="${SOFTWARE_PATH}/PLINK2/plink2"

## clean up the amr folder
if [  -d "hgdp-rawdata" ]; then
  rm -rf hgdp-rawdata
fi
mkdir -p hgdp-rawdata && cd $_

if ! [ -f "hgdp_all.pgen" ]; then
  wget "https://www.dropbox.com/s/hppj1g1gzygcocq/hgdp_all.pgen.zst?dl=1"
  mv "hgdp_all.pgen.zst?dl=1" "hgdp_all.pgen.zst"
  ./$zstd -d "hgdp_all.pgen.zst"
  rm "hgdp_all.pgen.zst"
fi

if ! [ -f "hgdp_all.pvar.zst" ]; then
  wget "https://www.dropbox.com/s/1mmkq0bd9ax8rng/hgdp_all.pvar.zst?dl=1"
  mv "hgdp_all.pvar.zst?dl=1" "hgdp_all.pvar.zst"
fi

if ! [ -f "hgdp_all.psam" ]; then
  wget "https://www.dropbox.com/s/0zg57558fqpj3w1/hgdp.psam?dl=1"
  mv "hgdp.psam?dl=1" "hgdp_all.psam"
fi

if ! [ -f "hgdp_wgs.20190516.metadata.txt" ]; then
  wget "ftp://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/metadata/hgdp_wgs.20190516.metadata.txt"
fi


## set missing IDs to unique values to avoid these being detected as repeated IDs
./$PLINK2 --pfile hgdp_all vzs \
          --set-missing-var-ids '@:#' \
          --make-just-pvar zs \
          --out hgdp_all_uniq
          
# replace data
mv hgdp_all_uniq.pvar.zst hgdp_all.pvar.zst

# trash
rm hgdp_all_uniq.log

## convert files from `zvs` format to PLINK binary format
./$PLINK2 \
    --pfile hgdp_all vzs \
    --var-filter \
    --snps-only just-acgt \
    --max-alleles 2 \
    --make-bed \
    --out hgdp_wgs_autosomes

## add subpopulation labels to FAM files
awk -F" " '{
  if (NR == FNR) {
    population[$1]=$6;
    if ($10 == "M") {
      gender[$1]=1;
    } else if ($10 == "F") {
      gender[$1]=2;
    } else {
      gender[$1]=0;
    }
  } 
  else print population[$2],$2,$3,$4,gender[$2],$6;
}' hgdp_wgs.20190516.metadata.txt hgdp_wgs_autosomes.fam > hgdp_wgs_autosomes.fam.NEW
mv hgdp_wgs_autosomes.fam.NEW hgdp_wgs_autosomes.fam

## preserve loci that:
#  1. have MAF >= 0.01
## 2. are in approximate linkage equilibrium with each other
./$PLINK2 \
    --bfile hgdp_wgs_autosomes \
    --make-bed \
    --out hgdp_wgs_autosomes_maf_0.01 \
    --maf 0.01
    
## this command determines the loci to keep or exclude
./$PLINK2 --bfile hgdp_wgs_autosomes_maf_0.01 --indep-pairwise 1000kb 0.3 --out hgdp_wgs_autosomes_maf_0.01
  
## this actually filters the data
./$PLINK2 \
    --bfile hgdp_wgs_autosomes_maf_0.01 \
    --extract hgdp_wgs_autosomes_maf_0.01.prune.in \
    --make-bed \
    --out hgdp_wgs_autosomes_final

rm hgdp_wgs_autosomes.{bed,bim,fam,log}
rm hgdp_wgs_autosomes_maf_0.01.{bed,bim,fam,log,prune.in,prune.out}
rm hgdp_wgs_autosomes_final.log
