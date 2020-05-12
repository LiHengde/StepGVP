#!/usr/bin/env Rscript
# linux version
arg  <- commandArgs(TRUE)
Help_Function <- function(){
cat(
 "\t************************************************************************\n
  \t*                Genetic Relationship Matrix                           *\n
  \t*                  Author: Hengde Li, Ph.D.                            *\n
  \t*                  Email: hengde@aliyun.com                            *\n
  \t*             Chinese Academy of Fishery Sciences                      *\n
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
  \t      GDGRM.R -genotype [or -geno,-g] genotype_file\n
  \t              -format   [or -form,-f] gentype_format\n 
  \t              -pedigree [or -ped, -p] pedigree_file\n
  \t              -method   [or -met, -m] method\n  
  \t              -snp      [or -map, -s] map_information\n    
  \t              -output   [or -out, -o] output_file_prefix\n   
  \t              -level    [or       -l] chromosome_level\n 
  \t              -type     [or       -t] type_of_matrix\n 
  \t              -unit     [or       -u] unit_megabase \n   
  \t              -chr      [or       -c] chr or 'all']\n 

  \t The input files consists:\n
  \t    Genotype file: id, allele1.1 allele1.2 allele2.1 allele2.2 ...\n
  \t               or: id, geno1 geno2 geno3 geno4 ... \n

  \t    Pedigree file: id, father, mother, sex \n
  \t    Map file: MarkerName, Chromosome, Position\n

  \t The output files consists:\n
  \t genome level:\n
  \t    prefix_type\
  \t chromosome level:\n 
  \t    prefix_chr_type, prefix_chr_type_i\n 

  \t segment level:\n 
  \t    prefix_chr.segi\n 
  
  \t For the details of usage, please see manuals.  
  \tIf any bug was found, please contact the author.\n")
}


if (any(arg=='-h') | any(arg=='-help')) Help_Function()


fgen <- NULL        # genotype file
fped <- NULL        # pedigree file, id, father, mother
fmap <- NULL        # Marker_Name, Chr, Position, no header
fout <- NULL        # output file prefix
meth <- '1'         # method: 1, 2, 3
unit <- 2           # every unit megabase as a subset
form <- 'genotype'  # 'allele'
lev  <- 'chromosome' #'segment' or 'seg', genome or gen
chrs <- 'all'         # all chromosomes 
rang <- NULL

for (i in 1:length(arg)) {

   if (arg[i]=='-geno' | arg[i]=='-genotype'| arg[i]=='-g') fgen <- arg[i+1]
   if (arg[i]=='-form' | arg[i]=='-format'  | arg[i]=='-f') form <- arg[i+1]
   if (arg[i]=='-ped'  | arg[i]=='-pedigree'| arg[i]=='-p') fped <- arg[i+1]
   if (arg[i]=='-out'  | arg[i]=='-output'  | arg[i]=='-o') fout <- arg[i+1]
   if (arg[i]=='-met'  | arg[i]=='-method'  | arg[i]=='-m') meth <- arg[i+1]
   if (arg[i]=='-map'  | arg[i]=='-snp'     | arg[i]=='-s') fmap <- arg[i+1] 
     
   if (arg[i]=='-level'| arg[i]=='-l') lev  <- arg[i+1]   
   if (arg[i]=='-type' | arg[i]=='-t') type <- arg[i+1]     # A or D
   if (arg[i]=='-unit' | arg[i]=='-u') unit <- arg[i+1]
   
   if (arg[i]=='-chr'  | arg[i]=='-c') {
      if (arg[i+1]!='ALL' & arg[i+1]!='all') chrs <- as.numeric(arg[i+1])
   }
   if (arg[i]=='-range'  | arg[i]=='-r')  {
      rang <- rep(NA,2)
      for (j in 1:2) rang[j] <- as.numeric(arg[i+j])
   }
}

#  incomplete command error
if (is.null(fgen)|is.null(fped)|is.null(fout)|is.null(fmap)) {
   cat('\tPlease check the command!\n')
   q('no')
}


map <- read.table(fmap,stringsAsFactors=FALSE)
names(map) <- c('Marker','Chr','Pos')


# chromosome specification error or unit error
if (lev=='seg' | lev=='segment') {
   unit <- as.numeric(unit)
   if (is.na(unit)) {
     cat('\tPlease check the command!\n')
     q('no')
   }
}

if (lev=='poi' | lev=='point') {
  if (!is.null(rang) & chrs=='all') {
    cat('\tPlease check [chr] and [rang] in the command line!\n')
    q('no')
  }
} 
 
#if (!is.null(rang) & (lev!='seg' & lev!='segment')) {
#  cat('\tPlease check [range] or [level] in the command line!\n')
#  q('no')
#} 

if (!is.null(rang)) {
  if (length(rang)!=2) {
    cat('\tPlease check [range] in the command line!\n')
    q('no')  
  }
  rang1 <- rang[2]-rang[1]
  if ((lev=='seg' | lev=='segment') & rang1<unit) {
    cat('\tPlease check [unit] or [range] in the command line!\n')
    q('no')  
  }
  if ((lev=='poi' | lev=='point') & rang1<=0) {
    cat('\tPlease check [range] in the command line!\n')
    q('no')  
  }
}


ped <- as.matrix(read.table(fped)[,1:3])
gen <- as.matrix(read.table(fgen)[,-1])
fun1 <- function(i) rep(i,i)
fun2 <- function(i) 1:i
u    <- unlist(sapply(1:nrow(gen),fun2))
v    <- unlist(sapply(1:nrow(gen),fun1))
uv   <- cbind(u,v)
ord  <- order(u,v)
uv   <- uv[ord,]
vu   <- cbind(v,u)


Marker_Mean_Var <- function(ped, gg) {
  f <- which(ped[,2]==0 & ped[,3]==0)      # f: founder
  #p: the second allele frequency
  p <- sapply(1:ncol(gg),function(i)(sum(gg[f,i]==2)+0.5*sum(gg[f,i]==1))/length(f))
  list(mean=2*p, variance=2*p*(1-p))
}

'SNP_Grouping' <- function(chr_map, unit){     # unit: Megabase
   start_pos <- floor(with(chr_map, min(Pos))/(unit*1000000))
   end_pos <- ceiling(with(chr_map, max(Pos)+1)/(unit*1000000))
   break_pos <- 1000000*(start_pos:end_pos)
   out <- list()
   for (i in 1:(length(break_pos)-1)) {
      out[[i]] <- with(chr_map, which(Pos>=break_pos[i] & Pos<break_pos[i+1])) 
   }
   return(out)
}


'GRM_Mat' <- function(g, type='a', ped, method='1') { 
   mv <- Marker_Mean_Var(ped, g)
   if (method=='1'){
      if (type=='a'|type=='A'){
        for(i in 1:ncol(g)) g[,i] <- g[,i]-mv$mean[i]
        G     <- g%*%t(g)/sum(mv$variance)
      }
      if (type=='d'|type=='D'){
        g[g==2] <- 0
        v <- apply(g,2,var)
        dmv <- list()
        dmv$mean <- mv$variance[v>0]
        dmv$variance <- dmv$mean*(1-dmv$mean)
        g <- g[,v>0]
        for(i in 1:ncol(g)) g[,i] <- g[,i]-dmv$mean[i]
        G     <- g%*%t(g)/sum(dmv$variance)
      }

   }
   if (method=='2'){
      if (type=='a'|type=='A'){
        for(i in 1:ncol(g)) g[,i] <- g[,i]-mv$mean[i]
        G     <- g%*%diag(1/mv$variance)%*%t(g)/ncol(g)
      }
      if (type=='d'|type=='D'){
        g[g==2] <- 0
        v <- apply(g,2,var)
        dmv <- list()
        dmv$mean <- mv$variance[v1>0]
        dmv$variance <- dmv$mean*(1-dmv$mean)
        g <- g[,v>0]
        for(i in 1:ncol(g)) g[,i] <- g[,i]-dmv$mean[i]
        G     <- g%*%diag(1/dmv$variance)%*%t(g)/ncol(g)
      }
   }
   if (method=='3'){
      if (type=='a'|type=='A'){
        g <- scale(g)
        G <- g%*%t(g)/ncol(g)
      }
      if (type=='d'|type=='D'){
        g[g==2] <- 0
        v <- apply(g,2,var)
        g <- g[,v>0]
        g <- scale(g)
        G <- g%*%t(g)/ncol(g)
      }
   }    
   G[vu]
}



uchr <- unique(map$Chr)
if (chrs!='all') uchr <- chrs
nchr <- length(uchr)



if (lev=='seg' | lev=='segment') {
  for (j in 1:nchr){
    chr <- uchr[j]
    gen1 <- gen[,with(map,which(Chr==chr))] 
    map1 <- subset(map, Chr==chr) 
    Var1 <- apply(gen1, 2, var)  
    gen1 <- gen1[,Var1>0 ]
    map1 <- map1[ Var1>0,] 
    
    if (!is.null(rang)) {
      s1 <- floor(min(map1$Pos)/1000000)
      o1 <-with(map1,which(Pos/1000000>s1+rang[1] & Pos/1000000<s1+rang[2]))
      gen1 <- gen1[,o1 ]   
      map1 <- map1[ o1,]
    }
      
    snp_set <- SNP_Grouping(map1,unit)
    
    for (i in 1:length(snp_set)) {
      gg  <- gen1[,snp_set[[i]]]
      si  <- paste(chr,i,sep='.s')    
      fo  <- paste(fout,si,type,sep='_')
      out <- GRM_Mat(gg,type,ped,method=meth)
      write.table(out,fo,row.names=F,col.names=F,quote=F,sep='\t')
    }
  }
}


if (lev=='chr' | lev=='chromosome') {
  for (i in 1:nchr) {
     chr <- uchr[i]
     ggg <- gen[,with(map,which(Chr==chr))]  
     Var <- apply(ggg, 2, var) 
     ggg <- ggg[,Var>0 ]  
     fo  <- paste(fout,chr,type,sep='_')
     out <- GRM_Mat(ggg,type,ped,method=meth)
     write.table(out,fo,row.names=F,col.names=F,quote=F,sep='\t')
  }  
}

if (lev=='gen' | lev=='genome') {  
   Var <- apply(gen, 2, var) 
   ggg <- gen[,Var>0 ]  
   fo  <- paste(fout,type,sep='_')
   out <- GRM_Mat(ggg,type,ped,method=meth)
   write.table(out,fo,row.names=F,col.names=F,quote=F,sep='\t')  
}



if (lev=='poi' | lev=='point') {
  for (j in 1:nchr){
    chr <- uchr[j]
    gen1 <- gen[,with(map,which(Chr==chr))] 
    Map1 <- subset(map, Chr==chr) 
    Var1 <- apply(gen1, 2, var)  
    gen1 <- gen1[,Var1>0 ]
    map1 <- Map1[ Var1>0,] 
    
    if (!is.null(rang)) {
      s1 <- floor(min(map1$Pos)/1000000)
      o1 <-with(map1,which(Pos/1000000>s1+rang[1] & Pos/1000000<s1+rang[2]))
      gen1 <- gen1[,o1 ]   
      map1 <- map1[ o1,]
    }
          
    for (i in 1:nrow(map1)) {
      gg  <- gen1[,i]
      ii  <- which(Map1$Marker==map1$Marker[i])
      si  <- paste(chr,ii,sep='.p') 
      flag <- 0
      if (type=='D'|type=='d'){
        GG1 <- gg
        GG1[GG1==2] <- 0
        if (var(GG1)==0) flag <- 1 
      }
      if (flag==0) {   
        fo  <- paste(fout,si,type,sep='_')
        out <- GRM_Mat(gg,type,ped,method=meth)
        write.table(out,fo,row.names=F,col.names=F,quote=F,sep='\t')
      }
    }
  }
}

q('no')

