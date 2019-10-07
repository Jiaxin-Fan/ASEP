# Utility Functions
#
# Author: Jiaxin Fan
###################################################

#' @title Haplotype pseudo alignment
#'
#' @description This function is used to obtain the pseudo haplotype phase of the RNA-seq data for a given gene, and align the major alleles across individuals.
#' @param dat: bulk RNA-seq dataset of a given gene. Must contain variables:
#'             \itemize{
#'                \item One condition analysis: \cr
#'                  - `id`: character, individual identifier; \cr
#'                  - `ref`: numeric, the snp-level read counts for the reference allele if the haplotype phase of the data is unknown, and the snp-level read counts for allele aligned on paternal/maternal haplotype if haplotype phase is known;\cr
#'                  - `total`: numeric, snp-level total read counts for both alleles;\cr
#'                \item Two conditions analysis:  \cr
#'                  - `id`: character, individual identifier; \cr
#'                  - `snp`: character, the name/chromosome location of the heterzygous genetic variants; \cr
#'                  - `ref`: numeric, the snp-level read counts for the reference allele if the haplotype phase of the data is unknown, and the snp-level read counts for allele aligned on the same paternal/maternal haplotype for both conditions if haplotype phase is known;\cr
#'                  - `total`: numeric, snp-level total read counts for both alleles;\cr
#'                  - `group`: character, the condition each RNA-seq sample is obtained from (i.e., pre- vs post-treatment); \cr
#'                  - `ref_condition`: character, the condition used as the reference for pseudo haplotype phasing; \cr
#'             }
#' @param phased: a logical value indicates whether the haplotype phase of the data is known or not. Default is FALSE
#' @param n_condition: a character string indicates whether the RNA-seq data contains data from only one condition or two conditions (i.e. normal vs diseased). Possible values are "one" or "two". Default is "one"
#' @return The psudo-phased RNA-seq data, with one more column "major" indicates the read counts for major alleles aligned across individuals
#' @export
#' @import stats

phasing<-function(dat, phased=FALSE, n_condition="one"){
  if (phased){
    if (n_condition == "one"){
      diff<-aggregate(ref-total/2 ~ id, dat, sum)
      # phase==1 means need to flip
      phasesub=data.frame(diff$id,ifelse(diff[,2]<0,1,0))
      colnames(phasesub)=c("id","phase")
      # merge the phase information with the raw data
      dat_phase=merge(dat,phasesub)
    }
    else if (n_condition == "two"){
      if("ref_condition" %in% colnames(dat)){
        ref_condition = unique(dat$ref_condition)
        if(length(ref_condition) == 1 & ref_condition %in% dat$group){
            dat_ref=dat[which(dat$group==ref_condition),]
            diff<-aggregate(ref-total/2 ~ id, dat_ref, sum)
            phasesub=data.frame(diff$id,ifelse(diff[,2]<0,1,0))
            colnames(phasesub)=c('id','phase')
            # merge the phase information with the raw data
            dat_phase=merge(dat,phasesub,by='id')
        }else{
            stop('Error: reference condition is not correctly specified.')
        }
      }else{
        stop('Error: input data should contain column "ref_condition".')
      }
    }
    else{
      stop('Error: "condition" only takes value "one" or "two".')
    }
  }
  else if (!phased){
    if (n_condition == "one"){
      phase=ifelse(dat$ref<dat$total/2,1,0)
      dat_phase<-data.frame(dat,phase)
    }
    else if (n_condition == "two"){
      if("ref_condition" %in% colnames(dat)){
        ref_condition = unique(dat$ref_condition)
        if(length(ref_condition) == 1 & ref_condition %in% dat$group){
          dat_ref<-dat[which(dat$group==ref_condition),]
          phasesub=data.frame(dat_ref$id,dat_ref$snp,ifelse(dat_ref$ref<dat_ref$total/2,1,0))
          colnames(phasesub)=c('id','snp','phase')
          # merge the phase information with the raw data, 'snp' must be differentiable
          dat_phase<- merge(dat, phasesub, by=c('id','snp'))
        }else{
          stop('Error: reference condition is not correctly specified.')
        }
      }else{
        stop('Error: input data should contain column "ref_condition".')
      }
    }
    else{
      stop('Error: "condition" only takes value "one" or "two".')
    }
  }
  else{
    stop('Error: "phased" must be a boolean.')
  }
  # 1 means need to flip
  dat_phase$major = dat_phase$ref
  dat_phase<-within(dat_phase, major[phase==1] <- total[phase==1]-ref[phase==1])
  dat_phase = dat_phase[ , !(names(dat_phase) %in% c("phase"))]
  return(dat_phase)
}

#' @title Fit generalized linear mixed-effect model
#'
#' @description This function is used to fit the generalized linear mixed-effect model through nonparametric maximum likelihood estimation for a given gene.
#' @param dat_phase: Psudo-phased bulk RNA-seq data of a given gene
#' @param n_condition: a character string indicates whether to perform one condition analysis or two conditions analysis. Possible values are "one" or "two". Default is "one"
#' @param resampled_data: a logical value indicates whether the data analyzed is obtained from resampling or not. Used for two conditions analysis only. Default is FALSE
#' @return One of: \cr
#'    \itemize{
#'        \item the test statistic for ASE detection - when modeling the original/resampled data under one condition analysis;
#'        \item the test statistic for differential ASE detection - when modeling the resampled data under two conditions analysis;
#'        \item the test statistic for differential ASE detection & the estimated ASE effect in the pooled sample - when modeling the original data under two conditions analysis;
#'    }
#' @export
#' @import npmlreg

modelFit<-function(dat_phase, n_condition="one", resampled_data=FALSE){
  mudiff<-NA
  if (n_condition == "one"){
    np=tryCatch({allvc(log(major/total)~1,random=~1|id,
                       family=gaussian, data=dat_phase, k=2, spike.protect = 1,
                       random.distribution='np',lambda=0.99, plot.opt = 0, verbose = FALSE)},
                error=function(e){FALSE})
    # check model convergence
    if((!is.vector(np))){
      if(np$EMconverged){
        mu1 = as.numeric(as.character(np$coefficients[1]))
        mu2 = as.numeric(as.character(np$coefficients[2]))
        sd1 = as.numeric(as.character(np$sdev$sdevk[1]))
        sd2 = as.numeric(as.character(np$sdev$sdevk[2]))
        mudiff = abs(mu1-mu2)/sqrt(sd1^2+sd2^2)
      }
    }
    return(mudiff)
  }
  else if (n_condition == "two"){
    np=tryCatch({allvc(log(major/total)~group,random=~group|id,
                       family=gaussian,lambda=0.99, spike.protect = 1,
                       data=dat_phase, k=2, random.distribution='np', plot.opt = 0, verbose = FALSE)},
                error=function(e){FALSE})
    # check model convergence
    if((!is.vector(np))){
      if(np$EMconverged){
        name_group = grep("group",names(np$coefficients),value=TRUE)[1]
        mass1diff=abs(exp(as.numeric(as.character(np$coefficients[names(np$coefficients) == "MASS1"]))+
                          as.numeric(as.character(np$coefficients[names(np$coefficients) == name_group]))+
                          as.numeric(as.character(np$coefficients[names(np$coefficients) == paste("MASS1:",name_group,sep='')])))-
                        exp(as.numeric(as.character(np$coefficients[names(np$coefficients) == "MASS1"]))))
        mass2diff=abs(exp(as.numeric(as.character(np$coefficients[names(np$coefficients) == "MASS2"]))+
                          as.numeric(as.character(np$coefficients[names(np$coefficients) == name_group])))-
                        exp(as.numeric(as.character(np$coefficients[names(np$coefficients) == "MASS2"]))))
        mudiff=unname(abs(mass1diff-mass2diff))
        # obtain the homozygous and heterozygous cis-rSNP group prediction
        # individual from homozygous cis group should have lower ASE, therefore, smaller mean than people from heterozygous group
        if(!resampled_data){
          if(np$coefficients[names(np$coefficients) == "MASS1"] <= np$coefficients[names(np$coefficients) == "MASS2"]){
            homo='1'
            heter='2'
          }
          else{
            homo='2'
            heter='1'
          }
          sumup<-aggregate(cbind(major,total) ~ id, dat_phase, sum)
          pool=data.frame(sumup$id, sumup$major/sumup$total)
          colnames(pool)=c("id","poolp")
          est.prob=data.frame(np$post.prob,dat_phase$id)
          colnames(est.prob)=c('1','2','id')
          est.prob=unique(est.prob)
          newprob=merge(est.prob,pool,by='id')
          newprob$poolp=as.numeric(unlist(newprob[homo]))*0.5+
            as.numeric(unlist(newprob[heter]))*newprob$poolp
          dat_allpoolp=merge(dat_phase, newprob, by="id")
          return(list(mudiff,dat_allpoolp))
        }
      }
    }
    return(mudiff)
  }
  else{
    stop('Error: "condition" only takes value "one" or "two".')
  }
}

#' @title Perform ASE analysis in the population for a given gene
#'
#' @description This function is used to perform ASE detection for one condition analysis of a given gene.
#' @param dat: bulk RNA-seq data of a given gene. Must contain variables:
#'             \itemize{
#'                  \item `id`: character, individual identifier;
#'                  \item `ref`: numeric, the snp-level read counts for the reference allele if the haplotype phase of the data is unknown, and the snp-level read counts for allele aligned on paternal/maternal haplotype if haplotype phase is known;
#'                  \item `total`: numeric, snp-level total read counts for both alleles;
#'             }
#' @param phased: a logical value indicates whether the haplotype phase of the data is known or not. Default is FALSE
#' @param n_resample: a numeric value indicates the maximum number of resamplings performed to obtain estimated p-value. Default is 10^6
#' @param adaptive: a logical value indicates whether the resampling is done through an adaptive procedure or not. Default is TRUE \cr
#'     By adaptive, it means first do 1000 resamplings, if the estimated p-value < 0.1, increase the number of resampling, by a factor of 10, to 10^4.
#'     if then the estimated p-value < 0.01, increase the number of resampling again, by a factor of 10, to 10^5.
#'     The procedure continuous until reaches the maximum number of resampling. \cr
#' @param parallel: a logical value indicates whether do parallel computing for the resampling precedure or not. Default is FALSE
#' @param n_core: a numeric value indicates number of clusters used for the parallel computing. Used only when parallel is set to TRUE. Default is 1
#' @return A vector with two elements:
#'       \itemize{
#'          \item `test statistic`: numeric, the test statistics of ASE effect;
#'          \item `p-value`: the estimated p-value of the test statistic;
#'       }
#' @import parallel

one_condition_analysis_Gene <- function(dat, phased=FALSE, n_resample=10^6, adaptive=TRUE, parallel=FALSE, n_core = 1){
  # check data
  dat=as.data.frame(dat)
  if (!(all(c("id","ref","total") %in% colnames(dat)))){
      stop('Error: Data must contain columns: "id","ref","total".')
  }
  # analysis
  dat_phase<-phasing(dat=dat, phased=phased, n_condition="one")
  dat_phase$total = as.numeric(as.character(dat_phase$total))
  mudiff = modelFit(dat_phase=dat_phase, n_condition="one")
  if(is.na(mudiff)){
    pvalue=NA
  }else{
    if(!parallel){
      if(!adaptive){
        xgene <- sapply(dat_phase$total,function(x){rbinom(n=n_resample,size=x,prob=0.5)})
        testnull <- apply(xgene,1,function(x){
          ref<-x
          dat_resam<-data.frame(ref,dat_phase[ , !(names(dat_phase) %in% c("ref","major"))])
          dat_resam_phase<-phasing(dat = dat_resam, phased=phased, n_condition="one")
          mudiffn = modelFit(dat_phase = dat_resam_phase, n_condition="one")
          return(mudiffn)
        })
        testnulldist<-na.omit(testnull)
        pvalue<-sum(ifelse(testnulldist>=mudiff,1,0))/length(testnulldist)
      }else{
        n_1 = log(n_resample,base=10)
        n_2 = floor(n_1)
        if(n_1==n_2){
          simulationList = c(10^3, diff(10^(3:n_2)))
        }else{
          simulationList = c(10^3, diff(10^(3:n_2)), n_resample-10^n_2)
        }
        pvalueList <- c(0.1^(1:(length(simulationList)-1)), 0.1^(length(simulationList)-1))[1:length(simulationList)]
        testnulldist = NULL
        for(i in 1:length(simulationList)){
          simulation <- simulationList[i]
          xgene<-sapply(dat_phase$total,function(x){rbinom(n=simulation,size=x,prob=0.5)})
          testnull<-apply(xgene,1,function(x){
            ref<-x
            dat_resam<-data.frame(ref,dat_phase[ , !(names(dat_phase) %in% c("ref","major"))])
            dat_resam_phase<-phasing(dat=dat_resam, phased=phased, n_condition="one")
            mudiffn = modelFit(dat_phase=dat_resam_phase, n_condition="one")
            return(mudiffn)
          })
          testnulldist<-c(testnulldist, na.omit(testnull))
          pvalue=sum(ifelse(testnulldist>=mudiff,1,0))/length(testnulldist)
          pvaluecutoff <- pvalueList[i]
          if(pvaluecutoff <= pvalue){
            break
          }
        }
      }
    }else{
      if(!adaptive){
        clus=makeCluster(n_core)
        xgene <- sapply(dat_phase$total,function(x){rbinom(n=n_resample,size=x,prob=0.5)})
        clusterExport(clus,list('dat_phase','phasing','modelFit'), envir=environment())
        testnull<-parRapply(clus,xgene,function(x){
          ref<-x
          dat_resam<-data.frame(ref,dat_phase[ , !(names(dat_phase) %in% c("ref","major"))])
          dat_resam_phase<-phasing(dat = dat_resam, phased=phased, n_condition="one")
          mudiffn = modelFit(dat_phase = dat_resam_phase, n_condition="one")
          return(mudiffn)
        })
        testnulldist<-na.omit(testnull)
        pvalue<-sum(ifelse(testnulldist>=mudiff,1,0))/length(testnulldist)
      }else{
        clus=makeCluster(n_core)
        n_1 = log(n_resample,base=10)
        n_2 = floor(n_1)
        if(n_1==n_2){
          simulationList = c(10^3, diff(10^(3:n_2)))
        }else{
          simulationList = c(10^3, diff(10^(3:n_2)), n_resample-10^n_2)
        }
        pvalueList <- c(0.1^(1:(length(simulationList)-1)), 0.1^(length(simulationList)-1))
        testnulldist = NULL
        for(i in 1:length(simulationList)){
          simulation <- simulationList[i]
          xgene<-sapply(dat_phase$total,function(x){rbinom(n=simulation,size=x,prob=0.5)})
          clusterExport(clus,list('dat_phase','phasing','modelFit'), envir=environment())
          testnull<-parRapply(clus,xgene,function(x){
            ref<-x
            dat_resam<-data.frame(ref,dat_phase[ , !(names(dat_phase) %in% c("ref","major"))])
            dat_resam_phase<-phasing(dat = dat_resam, phased=phased, n_condition="one")
            mudiffn = modelFit(dat_phase = dat_resam_phase, n_condition="one")
            return(mudiffn)
          })
          testnulldist<-c(testnulldist, na.omit(testnull))
          pvalue=sum(ifelse(testnulldist>=mudiff,1,0))/length(testnulldist)
          pvaluecutoff <- pvalueList[i]
          if(pvaluecutoff <= pvalue){
            break
          }
        }
      }
    }
  }
  return(c("test statistic" = mudiff, "p-value" = pvalue))
}

#' @title Perform differential ASE analysis in the population for a given gene
#'
#' @description This function is used to perform differential ASE detection for two conditions analysis of a given gene.
#' @param dat: bulk RNA-seq data of a given gene. Must contain variables \cr
#'       \itemize{
#'           \item `id`: character, individual identifier;
#'           \item `snp`: character, the name/chromosome location of the heterzygous genetic variants;
#'           \item `ref`: numeric, the snp-level read counts for the reference allele if the haplotype phase of the data is unknown, and the snp-level read counts for allele aligned on the same paternal/maternal haplotype for both conditions if haplotype phase is known;
#'           \item `total`: numeric, snp-level total read counts for both alleles;
#'           \item `group`: character, the condition each RNA-seq sample is obtained from (i.e., pre- vs post-treatment);
#'           \item `ref_condition`: character, the condition used as the reference for pseudo haplotype phasing;
#'       }
#' @param phased: a logical value indicates whether the haplotype phase of the data is known or not. Default is FALSE
#' @param n_resample: a numeric value indicates the maximum number of resamplings performed to obtain estimated p-value. Default is 10^6
#' @param adaptive: a logical value indicates whether the resampling is done through an adaptive procedure or not. Default is TRUE \cr
#'     By adaptive, it means first do 1000 resamplings, if the estimated p-value < 0.1, increase the number of resampling, by a factor of 10, to 10^4.
#'     if then the estimated p-value < 0.01, increase the number of resampling again, by a factor of 10, to 10^5.
#'     The procedure continuous until reaches the maximum number of resampling. \cr
#' @param parallel: a logical value indicates whether do parallel computing for the resampling precedure or not. Default is FALSE
#' @param n_core: a numeric value indicates number of clusters used for the parallel computing. Used only when parallel is set to TRUE. Default is 1
#' @return A vector with two elements:
#'       \itemize{
#'          \item `test statistic`: numeric, the test statistics of differential ASE effect;
#'          \item `p-value`: the estimated p-value of the test statistic;
#'       }
#' @import parallel

two_conditions_analysis_Gene <- function(dat, phased=FALSE, adaptive=TRUE, n_resample=10^6, parallel=FALSE, n_core = 1){
  # check data
  dat=as.data.frame(dat)
  if (all(c("ref","total","id","group","snp","ref_condition") %in% colnames(dat))){
    if (length(unique(dat$group))==2){
      # check if the data is paired
      if(phased){
        if(!(2 %in% aggregate(group~id, dat, function(x){length(unique(x))})[,"group"])){
          stop('Error: data needs to be paired.')
        }
      }else{
        if(!(2 %in% aggregate(group~id+snp, dat, function(x){length(unique(x))})[,"group"])){
          stop('Error: data needs to be paired.')
        }
      }
      dat=dat[order(dat$group,dat$id,dat$snp),]
    }else{
      stop('Error: "group" must have 2 levels.')
    }
  }else{
    stop('Error: Data must contain columns: "id","snp","ref","total","group","ref_condition".')
  }
  # analysis
  dat_phase<-phasing(dat=dat, phased=phased, n_condition="two")
  mod = modelFit(dat_phase=dat_phase, n_condition="two", resampled_data=FALSE)
  # model not converge
  if(length(mod)==1){
    mudiff = mod
    pvalue=NA
  }
  # converge
  if(length(mod)==2){
    mudiff = mod[[1]]
    dat_allpoolp = mod[[2]]
    gen<-cbind(dat_allpoolp$total,dat_allpoolp$poolp)

    if((!parallel) & (!adaptive)){
      xgene<-apply(gen,1,function(x){rbinom(n_resample,x[1],x[2])})
      testnull <- apply(xgene,1,function(x){
        ref<-x
        dat_resam<-data.frame(ref,dat_allpoolp[ , !(names(dat_allpoolp) %in% c("ref","major","poolp"))])
        dat_resam_phase<-phasing(dat=dat_resam, phased=phased, n_condition="two")
        mudiffn = modelFit(dat_phase=dat_resam_phase, n_condition="two", resampled_data=TRUE)
        return(mudiffn)
      })
      testnulldist<-na.omit(testnull)
      pvalue<-sum(ifelse(testnulldist>=mudiff,1,0))/length(testnulldist)
    }

    if(!parallel & adaptive){
      n_1 = log(n_resample,base=10)
      n_2 = floor(n_1)
      if(n_1==n_2){
        simulationList = c(10^3, diff(10^(3:n_2)))
      }else{
        simulationList = c(10^3, diff(10^(3:n_2)), n_resample-10^n_2)
      }
      pvalueList <- c(0.1^(1:(length(simulationList)-1)), 0.1^(length(simulationList)-1))
      testnulldist = NULL
      for(i in 1:length(simulationList)){
        simulation <- simulationList[i]
        xgene<-apply(gen,1,function(x){rbinom(simulation,x[1],x[2])})
        testnull<-apply(xgene,1,function(x){
          ref<-x
          dat_resam<-data.frame(ref,dat_allpoolp[ , !(names(dat_allpoolp) %in% c("ref","major","poolp"))])
          dat_resam_phase<-phasing(dat = dat_resam, phased=phased, n_condition="two")
          mudiffn = modelFit(dat_phase = dat_resam_phase, n_condition="two", resampled_data=TRUE)
          return(mudiffn)
        })
        testnulldist<-c(testnulldist, na.omit(testnull))
        pvalue=sum(ifelse(testnulldist>=mudiff,1,0))/length(testnulldist)
        pvaluecutoff <- pvalueList[i]
        if(pvaluecutoff <= pvalue){
          break
        }
      }
    }

    if (parallel & !adaptive){
      clus=makeCluster(n_core)
      xgene<-apply(gen,1,function(x){rbinom(n_resample,x[1],x[2])})
      clusterExport(clus,list('dat_allpoolp','phasing','modelFit'), envir=environment())
      testnull<-parRapply(clus,xgene,function(x){
        ref<-x
        dat_resam<-data.frame(ref,dat_allpoolp[ , !(names(dat_allpoolp) %in% c("ref","major","poolp"))])
        dat_resam_phase<-phasing(dat = dat_resam, phased=phased, n_condition="two")
        mudiffn = modelFit(dat_phase = dat_resam_phase, n_condition="two", resampled_data=TRUE)
        return(mudiffn)
      })
      testnulldist<-na.omit(testnull)
      pvalue<-sum(ifelse(testnulldist>=mudiff,1,0))/length(testnulldist)
    }

    if(parallel & adaptive){
      clus=makeCluster(n_core)
      n_1 = log(n_resample,base=10)
      n_2 = floor(n_1)
      if(n_1==n_2){
        simulationList = c(10^3, diff(10^(3:n_2)))
      }else{
        simulationList = c(10^3, diff(10^(3:n_2)), n_resample-10^n_2)
      }
      pvalueList <- c(0.1^(1:(length(simulationList)-1)), 0.1^(length(simulationList)-1))
      testnulldist = NULL
      for(i in 1:length(simulationList)){
        simulation <- simulationList[i]
        xgene<-apply(gen,1,function(x){rbinom(simulation,x[1],x[2])})
        clusterExport(clus,list('dat_allpoolp','phasing','modelFit'), envir=environment())
        testnull<-parRapply(clus,xgene,function(x){
          ref<-x
          dat_resam<-data.frame(ref,dat_allpoolp[ , !(names(dat_allpoolp) %in% c("ref","major","poolp"))])
          dat_resam_phase<-phasing(dat=dat_resam, phased=phased, n_condition="two")
          mudiffn = modelFit(dat_phase = dat_resam_phase, n_condition="two", resampled_data=TRUE)
          return(mudiffn)
        })
        testnulldist<-c(testnulldist, na.omit(testnull))
        pvalue=sum(ifelse(testnulldist>=mudiff,1,0))/length(testnulldist)
        pvaluecutoff <- pvalueList[i]
        if(pvaluecutoff <= pvalue){
          break
        }
      }
    }
  }
  return(c("test statistic" = unname(mudiff), "p-value" = pvalue))
}

#' @title Perform gene-level ASE analysis in the population across genes
#'
#' @description This function is used to detect ASE genes under one condition given bulk RNA-seq data that may contain multiple genes.
#' @param dat_all: bulk RNA-seq data. Must contain variables: \cr
#'       \itemize{
#'           \item `gene`: character, gene name;
#'           \item `id`: character, individual identifier;
#'           \item `ref`: numeric, the snp-level read counts for the reference allele if the haplotype phase of the data is unknown, and the snp-level read counts for allele aligned on paternal/maternal haplotype if haplotype phase is known;
#'           \item `total`: numeric, snp-level total read counts for both alleles;
#'       }
#' @param phased: a logical value indicates whether the haplotype phase of the data is known or not. Default is FALSE
#' @param n_resample: a numeric value indicates the maximum number of resamplings performed to obtain estimated p-value. Default is 10^6
#' @param adaptive: a logical value indicates whether the resampling is done through an adaptive procedure or not. Default is TRUE \cr
#'     By adaptive, it means first do 1000 resamplings, if the estimated p-value < 0.1, increase the number of resampling, by a factor of 10, to 10^4.
#'     if then the estimated p-value < 0.01, increase the number of resampling again, by a factor of 10, to 10^5.
#'     The procedure continuous until reaches the maximum number of resampling.
#' @param parallel: a logical value indicates whether do parallel computing for the resampling precedure or not. Default is FALSE
#' @param n_core: a numeric value indicates number of clusters used for the parallel computing. Used only when parallel is set to TRUE. Default is 1
#' @param save_out: a logical value indicates whether to write out the result for each gene stepwisely to a txt file. Default is FALSE
#' @param name_out: a character string indicates the output file name when save_out is set to TRUE, with the format of "XXX.txt"
#' @return A matrix with three columns:
#'       \itemize{
#'          \item `gene`: character, gene name;
#'          \item `test statistic`: numeric, the test statistics of ASE effect;
#'          \item `p-value`: the estimated p-value of the test statistic;
#'       }
#' @export
#' @import parallel

ASE_detection <- function(dat_all, phased=FALSE, adaptive=TRUE, n_resample=10^6, parallel=FALSE, n_core = 1, save_out=FALSE, name_out=NULL){
  dat_all = as.data.frame(dat_all)
  result_all=NULL
  if(save_out){
    file.create(name_out)
    re_out<-name_out
    write(c('gene','test statistic','p-value'), re_out, sep="\t",append=TRUE)
  }
  if ("gene" %in% colnames(dat_all)){
    dat_all$gene = as.character(dat_all$gene)
    for(gene in unique(dat_all$gene)){
      dat = dat_all[which(dat_all$gene==gene),]
      result = one_condition_analysis_Gene(dat=dat, phased=phased, adaptive=adaptive, n_resample=n_resample, parallel=parallel, n_core=n_core)
      result_all = rbind(result_all,c('gene'=gene,result))
      if(save_out){
        write(c('gene'=gene,result), re_out, sep="\t",append=TRUE)
      }
    }
  }else{
    stop('Error: Data must contain column: "gene".')
  }
  return(result_all)
}

#' @title Perform gene-level differential ASE analysis in the population across genes
#'
#' @description This function is used to detect differential ASE genes between two conditions given bulk RNA-seq data that may contain multiple genes.
#' @param dat_all: bulk RNA-seq data. Must contain variables: \cr
#'       \itemize{
#'          \item `gene`: character, gene name;
#'          \item `id`: character, individual identifier;
#'          \item `snp`: character, the name/chromosome location of the heterzygous genetic variants;
#'          \item `ref`: numeric, the snp-level read counts for the reference allele if the haplotype phase of the data is unknown, and the snp-level read counts for allele aligned on the same paternal/maternal haplotype for both conditions if haplotype phase is known;
#'          \item `total`: numeric, snp-level total read counts for both alleles;
#'          \item `group`: character, the condition each RNA-seq sample is obtained from (i.e., pre- vs post-treatment);
#'          \item `ref_condition`: character, the condition used as the reference for pseudo haplotype phasing;
#'       }
#' @param phased: a logical value indicates whether the haplotype phase of the data is known or not. Default is FALSE
#' @param n_resample: a numeric value indicates the maximum number of resamplings performed to obtain estimated p-value. Default is 10^6
#' @param adaptive: a logical value indicates whether the resampling is done through an adaptive procedure or not. Default is TRUE \cr
#'     By adaptive, it means first do 1000 resamplings, if the estimated p-value < 0.1, increase the number of resampling, by a factor of 10, to 10^4.
#'     if then the estimated p-value < 0.01, increase the number of resampling again, by a factor of 10, to 10^5.
#'     The procedure continuous until reaches the maximum number of resampling.
#' @param parallel: a logical value indicates whether do parallel computing for the resampling precedure or not. Default is FALSE
#' @param n_core: a numeric value indicates number of clusters used for the parallel computing. Used only when parallel is set to TRUE. Default is 1
#' @param save_out: a logical value indicates whether to write out the result for each gene stepwisely to a txt file. Default is FALSE
#' @param name_out: a character string indicates the output file name when save_out is set to TRUE, with the format of "XXX.txt"
#' @return  A matrix with three columns:
#'       \itemize{
#'          \item `gene`: character, gene name;
#'          \item `test statistic`: numeric, the test statistics of differential ASE effect;
#'          \item `p-value`: the estimated p-value of the test statistic;
#'       }
#' @export
#' @import parallel


differential_ASE_detection <- function(dat_all, phased=FALSE, adaptive=TRUE, n_resample=10^6, parallel=FALSE, n_core = 1, save_out=FALSE, name_out=NULL){
  dat_all = as.data.frame(dat_all)
  result_all=NULL
  if(save_out){
    file.create(name_out)
    re_out<-name_out
    write(c('gene','test statistic','p-value'), re_out, sep="\t",append=TRUE)
  }
  if ("gene" %in% colnames(dat_all)){
    dat_all$gene = as.character(dat_all$gene)
    for(gene in unique(dat_all$gene)){
      dat = dat_all[which(dat_all$gene==gene),]
      if( "ref_condition" %in% colnames(dat)){
        ref_condition = unique(dat$ref_condition)
        if(length(ref_condition)==1){
          result = two_conditions_analysis_Gene(dat=dat, phased=phased, adaptive=adaptive, n_resample=n_resample, parallel=parallel, n_core=n_core)
          result_all = rbind(result_all,c('gene'=gene,result))
          if(save_out){
            write(c('gene'=gene,result), re_out, sep="\t",append=TRUE)
          }
        }else{
          stop('Error: for each gene, there should be one and only one condition used as "ref_condition".')
        }
      }else{
        stop('Error: Data must contain column: "ref_condition".')
      }
    }
  }else{
    stop('Error: Data must contain column: "gene".')
  }
  return(result_all)
}

#' @title Boxplot of the estimated SNP-level ASE difference between two conditions
#'
#' @description This function is used to showed the estimated SNP-level ASE difference, i.e. major allele proportion difference,
#' between two condition samples across all individuals and SNPs after haplotype phasing.
#' The individuals are aligned based on an increasing trend of their median ASE differences across all snps for a given gene.
#' @param dat: bulk RNA-seq data of a given gene. Must contain variables: \cr
#'       \itemize{
#'          \item `gene`: character, gene name;
#'          \item `id`: character, individual identifier;
#'          \item `snp`: character, the name/chromosome location of the heterzygous genetic variants;
#'          \item `ref`: numeric, the snp-level read counts for the reference allele if the haplotype phase of the data is unknown, and the snp-level read counts for allele aligned on the same paternal/maternal haplotype for both conditions if haplotype phase is known;
#'          \item `total`: numeric, snp-level total read counts for both alleles;
#'          \item `group`: character, the condition each RNA-seq sample is obtained from (i.e., pre- vs post-treatment);
#'          \item `ref_condition`: character, the condition used as the reference for pseudo haplotype phasing;
#'       }
#' @param phased: a logical value indicates whether the haplotype phase of the data is known or not. Default is FALSE
#' @param minu_ref: a logical value indicates when calculating the difference, whether the "ref_condition" should be treated as the minuend, i.e. the difference is calculated as estimated ASE for minuend_condition minus that for the other condition. Default is TRUE
#' @return The boxplot (ggplot2::geom_boxplot) of the estimated SNP-level ASE difference between two conditions across individuals
#' @export
#' @import ggplot2

plot_ASE_diff <- function(dat, phased=FALSE, minu_ref=TRUE){
  if("gene" %in% colnames(dat)){
    gene = unique(dat$gene)
    if(length(gene)!=1){
      stop('Error: the input should contain data for only one gene.')
    }else{
      if("ref_condition" %in% colnames(dat)){
        ref_condition = unique(dat$ref_condition)
        if(length(ref_condition) != 1){
          stop('Error: for each gene, there should be one and only one condition used as "ref_condition".')
        }else{
          dat_phase = phasing(dat=dat, phased=phased, n_condition="two")
          dat_phase$MAF=dat_phase$major/dat_phase$total
          dat_phase$id = as.character(dat_phase$id)
          dat_phase = subset(dat_phase, select=c('id','snp','gene','total','group','MAF'))
          level = unique(dat_phase$group)
          if(length(level)==2){
            s1 = dat_phase[which(dat_phase$group!=ref_condition),]
            colnames(s1) = c('id','snp','gene','total','group',paste('MAF_',level[level != ref_condition], sep=''))
            s2 = dat_phase[which(dat_phase$group==ref_condition),]
            colnames(s2) = c('id','snp','gene','total','group',paste('MAF_',ref_condition, sep=''))
            s = merge(s1,s2,by=c('id','snp','gene'),all=T)
            s = na.omit(s)
            if(minu_ref){
              s$diff = unname(s[, names(s) %in% c(paste('MAF_',ref_condition, sep=''))]) -
                unname(s[, names(s) %in% c(paste('MAF_',level[level != ref_condition], sep=''))])
            }else if(!minu_ref){
              s$diff = unname(s[, names(s) %in% c(paste('MAF_',level[level != ref_condition], sep=''))])-
                unname(s[, names(s) %in% c(paste('MAF_',ref_condition, sep=''))])
            }else{
              stop('Error: "minu_ref" must be a boolen.')
            }
            dat_sort=aggregate(diff ~ id ,data=s, median)
            dat_sort=dat_sort[order(dat_sort$diff),]
            id_order=as.character(dat_sort$id)
            g = ggplot(data=s,aes(y = diff, x = id)) + geom_boxplot(aes(y = diff, x = id, group=id)) +
                geom_point(data=s,aes(y = diff, x = id, color=snp))+theme_classic()+
                geom_hline(yintercept = 0, linetype = "dashed",color='gray',size=1)+
                theme(plot.title = element_text(color="black", size=16, hjust=0.5),
                      axis.text.x = element_text(hjust = 0.5,size=12,angle=90,face="bold"),
                      axis.title.x = element_text(hjust = 0.5,size=14,face="bold"),
                      axis.text.y = element_text(size=12,face="bold"),
                      axis.title.y = element_text(size=14,face="bold"),
                      legend.title=element_text(size=14),
                      legend.text=element_text(size=12),
                      legend.position='right')+xlab('ID')+
                ylab('ASE difference')+scale_x_discrete(limits=id_order)+labs(colour="SNP")+
                ggtitle(unique(dat$gene))
              print(g)
          }else{
            stop('Error: Input must and only contain data from two conditions.')
          }
        }
      }else{
        stop('Error: Data must contain column: "ref_condition".')
      }
    }
  }else{
    stop('Error: Data must contain column: "gene".')
  }
}
