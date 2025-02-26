#' Genotypes of AMR subset of 1000 Genomes dataset
#'
#' @format A 10,000 \eqn{\times} 253 matrix of genotypes. Each element is an integer  of 0's, 1's and 2's.
#' \describe{
#'   Since the number of loci of the AMR subset of 1000 Genomes dataset is too large for a quick analysis,
#'   we created a subset of this dataset by first applying allele frequency filters and LD-pruning to the
#'   AMR dataset. We then randomly selected 10,000 SNPs out of LD-pruned SNP sets. We also re-ordered individuals
#'   according to their pairwise kinship level. This subset is available in `data/X_amr.rda`. This data has
#'   353 individuals and 10,000 loci. These data can be reproduced by scripts `data-raw/{amr.bash,amr.R}`.
#' }
"X_amr"

#' PLINK FAM file of AMR subset of 1000 Genomes dataset
#'
#' @format A dataframe of PLINK FAM file of AMR subset of 1000 Genomes dataset.
#' \describe{
#'  \item{fam}{The population code. MXL	abbreviates for Mexican Ancestry from Los Angeles USA.
#'      PUR	abbreviates for Puerto Rican from Puerto Rica. CLM abbreviates for Colombian from Medellian, Colombia.
#'      PEL abbreviates for	Peruvian from Lima, Peru.}
#'  \item{id}{The individual ID of 1000 Genomes Project.}
#'  \item{pat}{The father's ID of 1000 Genomes Project.}
#'  \item{mat}{The mother's ID of 1000 Genomes Project.}
#'  \item{sex}{The gender code. "1" denotes males. "2" denotes female. }
#'  \item{pheno}{The case-control status. All are -9.}
#' }
"fam_amr"

#' Admixture proportions of AMR subset of 1000 Genomes dataset
#'
#' @format A 3 \eqn{\times} 353 matrix of admixture proportions. Each column is the admixture proportions of the
#' corresponding individuals.
#'
"admix_props_amr"

#' Coancestry among populations of AMR subset of 1000 Genomes dataset
#'
#' @format A 3 \eqn{\times} 3 matrix of coancestry among populations. Each entry is coancestry coefficient for
#' the corresponding pair of the antecedent populations.
#'
"coanc_pops_amr"

#' Ancestral allele frequencies of AMR subset of 1000 Genomes dataset
#'
#' @format A length 10,000 vector of allele frequencies.
#'
"p_anc_amr"

#' Genotypes of HGDP dataset
#'
#' @format A 10,000 \eqn{\times} 929 matrix of genotypes. Each element is an integer  of 0's, 1's and 2's.
#' \describe{
#'   we created a subset of this dataset by first applying allele frequency filters and LD-pruning to the
#'   AMR dataset. We then randomly selected 10,000 SNPs out of LD-pruned SNP sets. This subset is available in
#'   `data/X_hgdp.rda`. This data has 929 individuals and 10,000 loci. These data can be reproduced by scripts
#'   `data-raw/{hgdp.bash,hgdp.R}`.
#' }
"X_hgdp"

#'

#' PLINK FAM file of HGDP dataset
#'
#' @format A dataframe of PLINK FAM file of HGDP dataset.
#' \describe{
#'  \item{fam}{The population code.}
#'  \item{id}{The individual ID of HGDP.}
#'  \item{pat}{The father's ID of HGDP.}
#'  \item{mat}{The mother's ID of HGDP.}
#'  \item{sex}{The gender code. "1" denotes males. "2" denotes female. }
#'  \item{pheno}{The phenotypic value}
#'  \item{population}{The population code.}
#'  \item{latitude}{}
#'  \item{longitude}{}
#'  \item{subpop}{The super-population label.}
#' }
"fam_hgdp"

#' Admixture proportions of HGDP dataset
#'
#' @format A 3 \eqn{\times} 929 matrix of admixture proportions. Each column is the admixture proportions of the
#' corresponding individuals.
#'
"admix_props_hgdp"

#' Coancestry among populations of HGDP dataset
#'
#' @format A 7 \eqn{\times} 7 matrix of coancestry among populations. Each entry is coancestry coefficient for
#' the corresponding pair of the antecedent populations.
#'
"coanc_pops_hgdp"

#' Ancestral allele frequencies of HGDP dataset
#'
#' @format A length 10,000 vector of allele frequencies.
#'
"p_anc_hgdp"
