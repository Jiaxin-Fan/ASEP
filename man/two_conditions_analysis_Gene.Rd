% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ASEP functions.R
\name{two_conditions_analysis_Gene}
\alias{two_conditions_analysis_Gene}
\title{Perform differential ASE analysis in the population for a given gene}
\usage{
two_conditions_analysis_Gene(
  dat,
  phased = FALSE,
  varList = NULL,
  adaptive = TRUE,
  n_resample = 10^6,
  parallel = FALSE,
  n_core = 1
)
}
\arguments{
\item{dat:}{bulk RNA-seq data of a given gene. Must contain variables \cr
\itemize{
    \item `id`: character, individual identifier;
    \item `snp`: character, the name/chromosome location of the heterzygous genetic variants;
    \item `ref`: numeric, the snp-level read counts for the reference allele if the haplotype phase of the data is unknown, and the snp-level read counts for allele aligned on the same paternal/maternal haplotype for both conditions if haplotype phase is known;
    \item `total`: numeric, snp-level total read counts for both alleles;
    \item `group`: character, the condition each RNA-seq sample is obtained from (i.e., pre- vs post-treatment);
    \item `ref_condition`: character, the condition used as the reference for pseudo haplotype phasing;
}}

\item{phased:}{a logical value indicates whether the haplotype phase of the data is known or not. Default is FALSE}

\item{varList:}{a character string specifies fomula of covariates that users want to adjusted in the model. An example could be "`var1`+`var2`". Default is NULL}

\item{n_resample:}{a numeric value indicates the maximum number of resamplings performed to obtain estimated p-value. Default is 10^6}

\item{adaptive:}{a logical value indicates whether the resampling is done through an adaptive procedure or not. Only applies when n_resample >= 10^3. Default is TRUE \cr
By adaptive, it means first do 1000 resamplings, if the estimated p-value < 0.1, increase the number of resampling, by a factor of 10, to 10^4.
if then the estimated p-value < 0.01, increase the number of resampling again, by a factor of 10, to 10^5.
The procedure continuous until reaches the maximum number of resampling. \cr}

\item{parallel:}{a logical value indicates whether do parallel computing for the resampling precedure or not. Default is FALSE}

\item{n_core:}{a numeric value indicates number of clusters used for parallel computing when parameter "parallel" is set to TRUE. Default is 1}
}
\value{
A vector with two elements:
      \itemize{
         \item `LRT statistic`: numeric, the likelihood ratio test (LRT) statistics of differential ASE effect;
         \item `p-value`: the estimated p-value of the LRT statistic;
      }
}
\description{
This function is used to perform differential ASE detection for two conditions analysis of a given gene.
}
