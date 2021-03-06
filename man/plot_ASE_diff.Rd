% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ASEP functions.R
\name{plot_ASE_diff}
\alias{plot_ASE_diff}
\title{Boxplot of the estimated SNP-level ASE difference between two conditions}
\usage{
plot_ASE_diff(dat, phased = FALSE, minu_ref = TRUE)
}
\arguments{
\item{dat:}{bulk RNA-seq data of a given gene. Must contain variables: \cr
\itemize{
   \item `gene`: character, gene name;
   \item `id`: character, individual identifier;
   \item `snp`: character, the name/chromosome location of the heterzygous genetic variants;
   \item `ref`: numeric, the snp-level read counts for the reference allele if the haplotype phase of the data is unknown, and the snp-level read counts for allele aligned on the same paternal/maternal haplotype for both conditions if haplotype phase is known;
   \item `total`: numeric, snp-level total read counts for both alleles;
   \item `group`: character, the condition each RNA-seq sample is obtained from (i.e., pre- vs post-treatment);
   \item `ref_condition`: character, the condition used as the reference for pseudo haplotype phasing;
}}

\item{phased:}{a logical value indicates whether the haplotype phase of the data is known or not. Default is FALSE}

\item{minu_ref:}{a logical value indicates when calculating the difference, whether the "ref_condition" should be treated as the minuend, i.e. the difference is calculated as estimated ASE for minuend_condition minus that for the other condition. Default is TRUE}
}
\value{
The boxplot (ggplot2::geom_boxplot) of the estimated SNP-level ASE difference between two conditions across individuals
}
\description{
This function is used to showed the estimated SNP-level ASE difference, i.e. major allele proportion difference,
between two condition samples across all individuals and SNPs after haplotype phasing.
The individuals are aligned based on an increasing trend of their median ASE differences across all snps for a given gene.
}
