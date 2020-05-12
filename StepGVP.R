#!/usr/bin/env Rscript
# linux version
# ideas to test random factor stepwise procedure
arg  <- commandArgs(TRUE)
# this function is to read genetic relationship matrix of low-triangular format
# the output of this function is square format 
Read_GRM <- function(low3) {
   fun1 <- function(n) 1:n
   fun2 <- function(n) rep(n,n)
   len  <- length(low3)
   n    <- (sqrt(len*8+1)-1)/2
   u    <- unlist(sapply(1:n, fun1))
   v    <- unlist(sapply(1:n, fun2))
   uv   <- cbind(v, u)
   out  <- matrix(NA, n, n)
   out[uv[,1:2]] <- low3
   out[uv[,2:1]] <- low3
   return(out)
}
# function to read model file
Read_Model <- function(fmod){
   fphe   <- trait  <- fix_col <- cov_col <- NULL
   rand_0 <- rand_1 <- rand_2  <- rmat_1  <- rmat2 <- NULL
   for (i in 1:9) {
      mod <- read.table(fmod, skip=i-1,nrows=1,stringsAsFactors=FALSE)
      if (mod[1,1]=='Phe_File:') fphe    <- mod[1,2]
      if (mod[1,1]=='RanMat_1:') rmat_1  <- mod[1,2]
      if (mod[1,1]=='RanMat_2:') rmat_2  <- mod[1,2]       
      if (mod[1,1]=='Dep_Col:' ) trait   <- as.numeric(mod[1,-1])
      if (mod[1,1]=='Fix_Col:' ) fix_col <- as.numeric(mod[1,-1])
      if (mod[1,1]=='Cov_Col:' ) cov_col <- as.numeric(mod[1,-1])
      if (mod[1,1]=='Random_0:') rand_0  <- as.numeric(mod[1,-1])
      if (mod[1,1]=='Random_1:') rand_1  <- as.numeric(mod[1,-1])
      if (mod[1,1]=='Random_2:') rand_2  <- as.numeric(mod[1,-1])
   }
   if(is.null(fphe)  |is.null(trait) |is.null(fix_col)|is.null(cov_col)|
      is.null(rand_0)|is.null(rand_1)|is.null(rand_2) |is.null(rmat_1) |
      is.null(rmat_2)) {
      cat('\t An Error was found in model file !!!\n')
      q('no')   
   }
   list(fphe=fphe,trait=trait,fix_col=fix_col,cov_col=cov_col,rand_0=rand_0,
        rand_1=rand_1,rand_2=rand_2,rmat_1=rmat_1,rmat_2=rmat_2)
}
# function to print help information
Help_Function <- function(){
cat(
 "\t************************************************************************\n
  \t*          StepWise Genomic Variance Partition Software                *\n
  \t*                  Author: Hengde Li, Ph.D.                            *\n
  \t*                  Email: hengde@aliyun.com                            *\n
  \t*              Chinese Academy of Fishery Sciences                     *\n
  \t************************************************************************\n 
  \tThe Chinese Academy of Fishery Sciences owns the copyright of the source\n
  \tcode of this software.\n 
  \tThis Software is free of charge to use for ACADEMIC and RESEARCH purpose,\n
  \thowever, technical support can not be provided. \n 
  \tUse of the software should be acknowledged in publications by reference to\n
  \tthe the user guide, web or paper.\n 
  \tYou are not allowed to distribute modified versions of the software under\n
  \tthe same nor a different name without the approval from the author.\n  
  \tFor commercial use of this software, Please contact the author.\n\n        
  \tCommand:\n 
  \t     stepGVP.R -model model_file -out out_file -iter n_iteration -alpha\n 
  \t                alpha -varpct var_percentage -core n_core\n
  \tor simplified as:\n
  \t     stepGVP.R -m model_file -o out_file -i n_iteration -a alpha -v\n 
  \t               var_percentage -c n_core\n
  \tfor help:\n
  \t     stepGVP.R -help or stepGVP.R -h\n
  \tmodel_file: \n
  \t     The file describes phenotypic file and the statistical model. Please\n
  \t     see the details by opening model_file.\n
  \tout_file:\n
  \t     The prefix of out_file, and five files will be output,\n 
  \t     xx_log:  log file of iteration\n
  \t     xx_sol:  solutions for random factors\n
  \t     xx_llik: loglik test in the first iteration\n
  \t     xx_summ: summary files of final model\n
  \t     xx_beta: solutions for fixed factors\n
  \talpha: \n
  \t     The significance level, default value is 0.05\n
  \tvar_percentage: \n
  \t     Threshold of variance ratio to drop matrices, default value is 0.05.\n
  \tn_core:\n
  \t     The number of threads in Linux.\n
  \tIf any bug was found, please contact the author.\n")
}
# function to construct design matrix for a specified factor
Design_Matrix <- function(v) {
  n <- length(v)  # nobs
  l <- unique(v)  # level
  o <- matrix(0, nobs, length(l)) # design matrix
  for (i in 1:length(l)) o[v==l[i],i]<-1
  return(o) 
}
fmod <- NULL
fphe <- NULL
mat1 <- NULL
mat2 <- NULL
Help <- FALSE      # whether print help information
fout <- NULL       # output file prefix
nite   <- 10
alpha  <- 0.05
varpct <- 0.05
ncore  <- 1
for (i in 1:length(arg)) {
  if (arg[i]=='-model' | arg[i]=='-m') fmod   <- arg[i+1]
  if (arg[i]=='-out'   | arg[i]=='-o') fout   <- arg[i+1]
  if (arg[i]=='-iter'  | arg[i]=='-i') nite   <- as.numeric(arg[i+1])
  if (arg[i]=='-alpha' | arg[i]=='-a') alpha  <- as.numeric(arg[i+1])
  if (arg[i]=='-varpct'| arg[i]=='-v') varpct <- as.numeric(arg[i+1])
  if (arg[i]=='-core'  | arg[i]=='-c') ncore  <- as.numeric(arg[i+1])
  if (arg[i]=='-help'  | arg[i]=='-h') Help   <- TRUE
}
if(Help) {
   Help_Function()
   q('no')
}
if ((is.null(fmod)|is.null(fout))&!Help) {
   cat('\tPlease check the command !!!\n')
   q('no')
}
if (.Platform$OS.type=='windows' & ncore>1) {
   cat('\tThe OS is Windows, Only One processor can be used !\n')
}
library(parallel)
library(EMMREML)
if (detectCores()<ncore)  ncore <- detectCores()
model <- Read_Model(fmod)
flog  <- paste(fout, 'log',  sep='_')
fsol  <- paste(fout, 'sol',  sep='_')
flik  <- paste(fout, 'llik', sep='_')
fsum  <- paste(fout, 'summ', sep='_')
fbeta <- paste(fout, 'beta', sep='_')
fnull <- paste(fout, 'null', sep='_')
threshold <- qchisq(1.0-alpha,df=1)
phe    <- read.table(model$fphe)
y      <- phe[,model$trait]
covs   <- NULL
if (model$cov_col[1]>0) covs <- as.matrix(phe[,model$cov_col])
nobs   <- nrow(phe)
X <- c()
fix_code <- matrix(NA,0,2)

for (ifix in model$fix_col) {
  X    <- cbind(X, Design_Matrix(phe[,ifix]))
  fc   <- unique(phe[,ifix])
  fix_code <- rbind(fix_code, cbind(rep(which(model$fix_col==ifix),length(fc)),fc))
}
if (!is.null(covs)) {
  X <- cbind(X, covs)
  fix_code <- rbind(fix_code, cbind((1:ncol(covs))+nrow(fix_code),rep(1,ncol(covs))))
}
nran0  <- length(model$rand_0)
if (model$rand_0[1]==0) nran0 <- 0
nran1  <- length(model$rand_1)
if (model$rand_1[1]==0) nran1 <- 0
nran2  <- length(model$rand_2)
if (nran2>1 | model$rand_2[1]==0) {
   cat('\tAn Error was found in model file !!!\n')
   q('no')
}
mat1 <- NULL
if (nran1>0) {
  mat1  <- as.matrix(read.table(model$rmat_1,colClasses='character'))
}
if (nran1!=nrow(mat1)) {
  cat('\tThe CorMat number is NOT equal to the number of type 1 random factor!!!\n')
  q('no')
}
mat2  <- as.matrix(read.table(model$rmat_2,colClasses='character')) 
K     <- list()
Z     <- list()
ran_code <- matrix(NA,0,2)
if (nran0>0) {
  for (iz in 1:length(model$rand_0)) {
    Z[[iz]]  <- Design_Matrix(phe[,model$rand_0[iz]])
    lev      <- unique(phe[,model$rand_0[iz]])
    K[[iz]]  <- diag(length(lev))
    ran_code <- rbind(ran_code, cbind(rep(iz,length(lev)), lev))
  }
}
if (nran1>0) {
  for (iz in 1:length(model$rand_1)) {
    Z[[iz+nran0]] <- Design_Matrix(phe[,model$rand_1[iz]])
    K[[iz+nran0]] <- Read_GRM(scan(mat1[iz]))
    nlev          <- nrow(K[[iz+nran0]]) 
    ran_code <- rbind(ran_code, cbind(rep(iz+nran0,nlev),1:nlev))    
  }
}
# Ztest: the design matrix for the factor being tested
Ztest <- Design_Matrix(phe[,model$rand_2])
nm0    <- nran0 + nran1
write(NULL, flog)
if (nm0==0) {
  emm   <- glm(y~X)
  llik0 <- -2*logLik(emm)
  write(paste('-2*logLik:',llik0,sep='\t'),fnull)
}
rFactor <- NULL
if (nm0>0){
  K <- K[1:nm0]
  Z <- Z[1:nm0]
  if (nran0>0) rFactor <- c(rFactor, paste('Random_0',1:nran0,sep='_'))
  if (nran1>0) rFactor <- c(rFactor, mat1)

  emm    <- emmremlMultiKernel(y,X,Z,K)
  llik0  <- -2.0*emm$loglik
  write(paste('-2*logLik:',llik0,sep='\t'),fnull)
  write('\nVariance Components:',fnull, append=T)
  var_comp <- c(emm$weights*emm$Vu, emm$Ve)
  var_comp <- data.frame(c(rFactor,'VarE'),var_comp)
  write.table(var_comp,fnull,append=T,row.names=F,col.names=F,quote=F,sep='\t')
}
# Iteration start
hit <- NULL
MPause <- NULL
MDrop  <- NULL
for (iter in 1:nite) {
# the null model
  if (nm0>0){   
    K <- K[1:nm0]
    Z <- Z[1:nm0]  
  }
  # read the matrices still in hit set
  if (length(hit)>0) {
    for (j in 1:length(hit)) {     
      K[[j+nm0]] <- Read_GRM(scan(hit[j]))
      Z[[j+nm0]] <- Ztest 
    }   
  }
  if (!is.null(MDrop)) {
     emm    <- emmremlMultiKernel(y,X,Z,K)
     llik0  <- -2.0*emm$loglik
     temp   <- data.frame(Iter=iter,LLIK=llik0,Hit=MDrop,Type='Drop')
     write.table(temp,flog,row.names=F,col.names=F,quote=F,sep='\t',append=TRUE)
  }
   # if all matrices are in the model, it is not necessary to iterate further. 
   if (sum(mat2%in%hit)==length(mat2)) break    
   # otherwise, continue
   imat   <- mat2[!(mat2%in%hit)] # test matrices excluding hitted matrices
   # to read the matrices fixed in K and Z, n order to avoid this situation
   # sometimes, if one matrix was dropped, the lengthes of K and Z changed,
   # but the matrices 1:nm0 will not.   
   VCA <- function(mat_file){ 
     # test matrices   
     K[[length(hit)+nm0+1]] <- Read_GRM(scan(mat_file))    
     Z[[length(hit)+nm0+1]] <- Ztest
     emm   <- emmremlMultiKernel(y,X,Z,K)
     c(-2.0*emm$loglik, emm$weight)  
   } 
   llik_weight <- mclapply(X=imat,FUN=VCA,mc.cores=ncore,mc.preschedule=TRUE)
   llik_weight <- matrix(unlist(llik_weight),nrow=length(imat),byrow=T)
   llik        <- llik_weight[,1]
   weight      <- as.matrix(llik_weight[,-1])  
   test <- llik0-llik   # test statistic is -2*logLik, ~Chisq distribution.
   # the statistic in iteration 1 will be printed for plot afterwards.
   if (iter==1) write.table(test,flik,row.names=F,col.names=F,quote=F,sep='\t')
   if (!any(test>threshold)) break  # if no sig., then exit
   if ( any(test>threshold)) {
     lmax   <- which.max(test) #decide which matrix should be selected.
     MPause <- imat[lmax]      # factor entered into the model
     if ( MPause%in%MDrop) break
     if (!MPause%in%MDrop) MDrop <- NULL    
     hit   <- c(hit, MPause)
     llik0 <- llik[lmax] 
     temp  <- data.frame(Iter=iter,LLIK=llik0,Hit=MPause,Type='Enter')
     write.table(temp,flog,row.names=F,col.names=F,quote=F,sep='\t',append=TRUE)     
     wt    <- weight[lmax,][((1:length(hit))+nm0)]
     if (any(wt<varpct)) {
        if (MPause==hit[which.min(wt)]) MDrop <- NULL
        if (MPause!=hit[which.min(wt)]) MDrop <- hit[which.min(wt)]
        hit <- hit[!hit%in%MDrop]
     }
   }
}
# final step, reestimate the variance components
if (length(hit)==0) {
  write("No significant one was hitted.",flog)
}
if (length(hit)>0) {
  if (nm0>0){
    K <- K[1:nm0]
    Z <- Z[1:nm0]
  }
  rFactor <- c(rFactor, hit, 'VarE')
  for (j in 1:length(hit)) {     
     K[[j+nm0]] <- Read_GRM(scan(hit[j]))
     Z[[j+nm0]] <- Ztest 
     nlev       <- nrow(K[[j+nm0]]) 
     ran_code   <- rbind(ran_code, cbind(rep(j+nm0,nlev),1:nlev))            
  }
  emm <- emmremlMultiKernel(y,X,Z,K)
  ## write summary file
  write('hitted factors:',fsum)
  write(hit, fsum, append=T)
  write('\n', fsum, append=T)
  write(paste('-2*logLik:',-2*emm$loglik,sep='\t'), fsum, append=T)
  write('\n',fsum, append=T)
  var_comp <- c(emm$weights*emm$Vu, emm$Ve)
  write('variance components:', fsum, append=T)
  var_comp <- data.frame(rFactor, var_comp) 
  write.table(var_comp,fsum,append=T,row.names=F,col.names=F,quote=F,sep='\t') 
  emm_uhat <- cbind(ran_code, emm$uhat)
  write.table(emm_uhat,fsol,row.names=F,col.names=F,quote=F,sep='\t')
  emm_beta <- cbind(fix_code,emm$betahat)
  write.table(emm_beta, fbeta, row.names=F,col.names=F,quote=F,sep='\t')       
}
q('no')
