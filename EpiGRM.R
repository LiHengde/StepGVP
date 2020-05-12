#!/usr/bin/env Rscript
# linux version
arg  <- commandArgs(TRUE)
Help_Function <- function(){
cat(
 "\t************************************************************************\n
  \t*            Epistasis Genetic Relationship Matrix                     *\n
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
  \t     EpiGRM.R -genotype [or -geno,-g] genotype_file\n
  \t              -format   [or -form,-f] gentype_format\n 
  \t              -pedigree [or -ped, -p] pedigree_file\n
  \t              -method   [or -met, -m] method\n  
  \t              -snp      [or -map, -s] map_information\n    
  \t              -output   [or -out, -o] output_file_prefix\n   
  \t              -level    [or       -l] chromosome_level\n   
  \t              -type     [or       -t] epistasis_type\n   
  \t              -unit     [or       -u] unit_megabase \n   
  \t              -chr      [or       -c] chr1 chr2 [or 'all']\n 
  \t              -range    [or       -r] range of chr1 and chr2, unit Mb  \n
  
  \t The input files consists:\n
  \t    Genotype file: id, allele1.1 allele1.2 allele2.1 allele2.2 ...\n
  \t               or: id, geno1 geno2 geno3 geno4 ... \n

  \t    Pedigree file: id, father, mother, sex \n
  \t    Map file: MarkerName, Chromosome, Position\n

  \t The output files consists:\n
  \t chromosome level:\n 
  \t    prefix_chr1_chr2_type, prefix_chr1_chr2_type_i\n 

  \t segment level:\n 
  \t    prefix_chr1.segi_chr2.segj_type, prefix_chr1.segi_chr2.segj_type_i\n 
  
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
lev  <- 'chromosome' #'segment' or 'seg',
type <- c('aa','ad','da','dd')         # epi-type, 'aa','ad','da', 'dd' 
chrs <- NULL         # all chromosomes or specific two chromosomes
chr1 <- chr2 <- NA
rang <- NULL

for (i in 1:length(arg)) {

   if (arg[i]=='-geno' | arg[i]=='-genotype'| arg[i]=='-g') fgen <- arg[i+1]
   if (arg[i]=='-form' | arg[i]=='-format'  | arg[i]=='-f') form <- arg[i+1]
   if (arg[i]=='-ped'  | arg[i]=='-pedigree'| arg[i]=='-p') fped <- arg[i+1]
   if (arg[i]=='-out'  | arg[i]=='-output'  | arg[i]=='-o') fout <- arg[i+1]
   if (arg[i]=='-met'  | arg[i]=='-method'  | arg[i]=='-m') meth <- arg[i+1]
   if (arg[i]=='-map'  | arg[i]=='-snp'     | arg[i]=='-s') fmap <- arg[i+1] 
     
   if (arg[i]=='-level'| arg[i]=='-l') lev  <- arg[i+1]   
   if (arg[i]=='-type' | arg[i]=='-t') typ  <- arg[i+1]
   if (arg[i]=='-unit' | arg[i]=='-u') unit <- arg[i+1]
   
   if (arg[i]=='-chr'  | arg[i]=='-c') {
      if (arg[i+1]=='ALL' | arg[i+1]=='all') chrs <- 'all'
      if (arg[i+1]!='ALL' & arg[i+1]!='all') { 
         chr1 <- as.numeric(arg[i+1])
         chr2 <- as.numeric(arg[i+2])
      }
   }
   if (arg[i]=='-range'  | arg[i]=='-r')  {
      rang <- rep(NA,4)
      for (j in 1:4) rang[j] <- as.numeric(arg[i+j])
   }

}

#  incomplete command error
if (is.null(fgen)|is.null(fped)|is.null(fout)|is.null(fmap)) {
 cat('\tPlease check [genotype],[pedigree],[snp] or [output] in the command!\n')
 q('no')
}
if (typ!='all' & typ!='ALL') type <- typ
#  epistasis type error
if (!type%in%c('aa','ad','da','dd','AA','AD','DA','DD')) {
   cat('\tPlease check [type] in the command line!\n')
   q('no')
}

if ((lev=='seg'|lev=='segment') & !is.null(chrs)) {
   cat('\tPlease check [chr] or [level] in the command line!\n')
   q('no')
}

# chromosome  specification error 1 
if(lev!='gen' & lev!='genome') {
  if (is.null(chrs) & (is.na(chr1) | is.na(chr2))) {
     cat('\tPlease check [chr] in the the command line!\n')
     q('no')
  }
}

map <- read.table(fmap,stringsAsFactors=FALSE)
names(map) <- c('Marker','Chr','Pos')

# chromosome specification error 2  
if(lev!='gen' & lev!='genome') {
  if (!is.null(chrs) & !is.na(chr1)) {
     if (!any(map$Chr==chr1)) {
        cat('\tPlease check [chr] in the command line!\n')
        q('no')
     }
  }
}
# chromosome specification error 2 
if(lev!='gen' & lev!='genome') {
  if (!is.na(chr2)) {
     if (!any(map$Chr==chr2)) {
        cat('\tPlease check [chr] in the command line!\n')
        q('no')
     }
  }
}
if (lev=='poi' | lev=='point') {
  if (is.na(chr1) | is.na(chr2)) {
    cat('\tPlease check [chr] in the command line!\n')
    q('no')
  }
}
# chromosome specification error or unit error
if (lev=='seg' | lev=='segment') {
   unit <- as.numeric(unit)
   if (is.na(unit)) {
     cat('\tPlease check [unit] in the command line!\n')
     q('no')
   }
}

 
if (!is.null(rang) & (lev!='seg' & lev!='segment')) {
  cat('\tPlease check [range] or [level] in the command line!\n')
  q('no')
} 

if (!is.null(rang)) {
  if (length(rang)!=4) {
    cat('\tPlease check [range] in the command line!\n')
    q('no')  
  }
}
if (!is.null(rang) ) {
  rang1 <- rang[2]-rang[1]
  rang2 <- rang[4]-rang[3]
  if (lev=='seg' | lev=='segment'){
    if (rang1<unit | rang2<unit) {
      cat('\tPlease check [unit] or [range] in the command line!\n')
      q('no')  
    }
  }
  if (lev=='poi' | lev=='point') {
    if (rang1<=0 | rang2<=0) {
      cat('\tPlease check [range] in the command line!\n')
      q('no')  
    }
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


'Epi_Mat' <- function(g1, g2, type='aa', ped, method='1') { 
   mv1 <- Marker_Mean_Var(ped, g1)
   mv2 <- Marker_Mean_Var(ped, g2)
   if (method=='1'){
      if (substr(type,1,1)=='a' | substr(type,1,1)=='A'){
        for(i in 1:ncol(g1)) g1[,i] <- g1[,i]-mv1$mean[i]      
        G1     <- g1%*%t(g1)/sum(mv1$variance) 
      }
      if (substr(type,2,2)=='a' | substr(type,2,2)=='A'){
        for(i in 1:ncol(g2)) g2[,i] <- g2[,i]-mv2$mean[i]
        G2     <- g2%*%t(g2)/sum(mv2$variance) 
      }
      if (substr(type,1,1)%in%c('d','D')){  
        g1[g1==2] <- 0
        v1 <- apply(g1,2,var)        
        dmv1 <- list()
        dmv1$mean <- mv1$variance[v1>0]
        dmv1$variance <- dmv1$mean*(1-dmv1$mean)
        g1 <- g1[,v1>0]
        for(i in 1:ncol(g1)) g1[,i] <- g1[,i]-dmv1$mean[i]
        G1     <- g1%*%t(g1)/sum(dmv1$variance) 
      }
      if (substr(type,2,2)%in%c('d','D')){ 
        g2[g2==2] <- 0
        v2 <- apply(g2,2,var)
        g2 <- g2[,v2>0]
        dmv2 <- list()
        dmv2$mean <- mv2$variance[v2>0]
        dmv2$variance <- dmv2$mean*(1-dmv2$mean)
        for(i in 1:ncol(g2)) g2[,i] <- g2[,i]-dmv2$mean[i]
        G2     <- g2%*%t(g2)/sum(dmv2$variance) 
      }           
   }
      if (method=='2'){
      if (substr(type,1,1)%in%c('a','A')){
        for(i in 1:ncol(g1)) g1[,i] <- g1[,i]-mv1$mean[i]
        G1     <- g1%*%diag(1/mv1$variance)%*%t(g1)/ncol(g1) 
      }
      if (substr(type,2,2)%in%c('a','A')){
        for(i in 1:ncol(g2)) g2[,i] <- g2[,i]-mv2$mean[i]
        G2     <- g2%*%diag(1/mv2$variance)%*%t(g2)/ncol(g2) 
      }
      if (substr(type,1,1)%in%c('d','D')){  
        g1[g1==2] <- 0
        v1 <- apply(g1,2,var)        
        dmv1 <- list()
        dmv1$mean <- mv1$variance[v1>0]
        dmv1$variance <- dmv1$mean*(1-dmv1$mean)
        g1 <- g1[,v1>0]
        for(i in 1:ncol(g1)) g1[,i] <- g1[,i]-dmv1$mean[i]
        G1     <- g1%*%diag(1/dmv1$variance)%*%t(g1)/ncol(g1) 
      }
      if (substr(type,2,2)%in%c('d','D')){ 
        g2[g2==2] <- 0
        v2 <- apply(g2,2,var)
        g2 <- g2[,v2>0]
        dmv2 <- list()
        dmv2$mean <- mv2$variance[v2>0]
        dmv2$variance <- dmv2$mean*(1-dmv2$mean)
        for(i in 1:ncol(g2)) g2[,i] <- g2[,i]-dmv2$mean[i]
        G2     <- g2%*%diag(1/dmv2$variance)%*%t(g2)/ncol(g2) 
      }           
   }
      if (method=='3'){
      if (substr(type,1,1)%in%c('a','A')){
        g1 <- scale(g1)
        G1 <- g1%*%t(g1)/ncol(g1) 
      }
      if (substr(type,2,2)%in%c('a','A')){
        g2 <- scale(g2)
        G2 <- g2%*%t(g2)/ncol(g2) 
      }
      if (substr(type,1,1)%in%c('d','D')){  
        g1[g1==2] <- 0
        v1 <- apply(g1,2,var)        
        g1 <- g1[,v1>0]
        g1 <- scale(g1)
        G1 <- g1%*%t(g1)/ncol(g1) 
      }
      if (substr(type,2,2)%in%c('d','D')){ 
        g2[g2==2] <- 0
        v2 <- apply(g2,2,var)
        g2 <- g2[,v2>0]
        g2 <- scale(g2)
        G2 <- g2%*%t(g2)/ncol(g2) 
      }           
   }    
   (G1*G2)[vu]
}


if (lev=='gen' | lev=='genome') {
  gen1 <- gen 
  map1 <- map 
  var1 <- apply(gen1, 2, var)
  gen1 <- gen1[,var1>0 ]
  map1 <- map1[ var1>0,]
  if (sum(c('aa','ad','dd','da')%in%type)==4) type <- c('aa','ad','dd') 
  for (t in type) {   
     fo  <- paste(fout,t,sep='_')
     out <- Epi_Mat(gen1,gen1,t,ped,method=meth)
     write.table(out,fo,row.names=F,col.names=F,quote=F,sep='\t')
  }
}

if (lev!='gen' & lev!='genome'){
  uchr <- unique(map$Chr)
  if (is.null(chrs)) uchr <- c(chr1,chr2)
  nchr <- length(uchr)
}

if ((lev=='seg' | lev=='segment')) {
  gen1 <- gen[,with(map,which(Chr==chr1))]
  gen2 <- gen[,with(map,which(Chr==chr2))]    
  map1 <- subset(map, Chr==chr1) 
  map2 <- subset(map, Chr==chr2)     
  var1 <- apply(gen1, 2, var)  
  var2 <- apply(gen2, 2, var)
  gen1 <- gen1[,var1>0 ]
  gen2 <- gen2[,var2>0 ]
  map1 <- map1[ var1>0,]  
  map2 <- map2[ var2>0,]
  
  if (!is.null(rang)) {
    s1 <- floor(min(map1$Pos)/1000000)
    s2 <- floor(min(map2$Pos)/1000000)
    o1 <-with(map1,which(Pos/1000000>s1+rang[1] & Pos/1000000<s1+rang[2]))
    o2 <-with(map2,which(Pos/1000000>s2+rang[3] & Pos/1000000<s2+rang[4]))  
    gen1 <- gen1[,o1 ]
    gen2 <- gen2[,o2 ]
    map1 <- map1[ o1,]
    map2 <- map2[ o2,]
  }
  
  snp_set1 <- SNP_Grouping(map1,unit)    
  snp_set2 <- SNP_Grouping(map2,unit)
  
  for (i in 1:length(snp_set1)) {
    gg1 <- gen1[,snp_set1[[i]]]
    si  <- paste(chr1,i,sep='.s')
    
    for (j in 1:length(snp_set2)) {
      gg2 <- gen2[,snp_set2[[j]]]
      sj  <- paste(chr2,j,sep='.s')
  
      for (t in type) {
         flag <- 0
         if (chr1==chr2 & (t%in%c('aa','dd')) & i>j) flag <- 1
         if (flag==0) {   
           fo  <- paste(fout,si,sj,t,sep='_')
           out <- Epi_Mat(gg1,gg2,t,ped,method=meth)
           write.table(out,fo,row.names=F,col.names=F,quote=F,sep='\t')
         }
      }
    }
  }
}

if (lev=='chr' & !is.null(chrs)) {
  for (i in 1:nchr) {
     gen1 <- gen[,with(map,which(Chr==uchr[i]))]
     map1 <- subset(map, Chr==uchr[i])   
     var1 <- apply(gen1, 2, var) 
     gen1 <- gen1[,var1>0 ] 
     map1 <- map1[ var1>0,] 
     for (j in 1:nchr) {
        gen2 <- gen[,with(map,which(Chr==uchr[j]))]
        map2 <- subset(map, Chr==uchr[j])  
        var2 <- apply(gen2, 2, var)
        gen2 <- gen2[,var2>0 ]
        map2 <- map2[ var2>0,]
        
        for (t in type) { 
          flag <- 0  
          if (!(t%in%c('aa','dd') & i<=j)) flag <- 1
          if ((sum(c('ad','da')%in%type)==2 & i==j & t=='da')) flag <- 1 
          if (flag==0) {
            fo  <- paste(fout,i,j,t,sep='_')
            out <- Epi_Mat(gen1,gen2,t,ped,method=meth)
            write.table(out,fo,row.names=F,col.names=F,quote=F,sep='\t')
          }
        }
     }
  }  
}

if ((lev=='chr' | lev=='chromosome')&is.null(chrs))  {
  gen1 <- gen[,with(map,which(Chr==chr1))]
  gen2 <- gen[,with(map,which(Chr==chr2))]
  map1 <- subset(map, Chr==chr1)
  map2 <- subset(map, Chr==chr2)
  var1 <- apply(gen1, 2, var)
  var2 <- apply(gen2, 2, var)
  gen1 <- gen1[,var1>0 ]
  gen2 <- gen2[,var2>0 ]
  map1 <- map1[ var1>0,]
  map2 <- map2[ var2>0,]
  if (chr1==chr2 & (typ=='all' | typ=='ALL')) type <- c('aa','ad','dd')
  for (t in type) {   
    fo  <- paste(fout,chr1,chr2,t,sep='_')
    out <- Epi_Mat(gen1,gen2,t,ped,method=meth)
    write.table(out,fo,row.names=F,col.names=F,quote=F,sep='\t')
  }
}

if ((lev=='poi' | lev=='point')) {
  gen1 <- gen[,with(map,which(Chr==chr1))]
  gen2 <- gen[,with(map,which(Chr==chr2))]    
  Map1 <- subset(map, Chr==chr1) 
  Map2 <- subset(map, Chr==chr2)     
  var1 <- apply(gen1, 2, var)  
  var2 <- apply(gen2, 2, var)
  gen1 <- gen1[,var1>0 ]
  gen2 <- gen2[,var2>0 ]
  map1 <- Map1[ var1>0,]  
  map2 <- Map2[ var2>0,]
  
  if (!is.null(rang)) {
    s1 <- floor(min(map1$Pos)/1000000)
    s2 <- floor(min(map2$Pos)/1000000)
    o1 <-with(map1,which(Pos/1000000>s1+rang[1] & Pos/1000000<s1+rang[2]))
    o2 <-with(map2,which(Pos/1000000>s2+rang[3] & Pos/1000000<s2+rang[4]))  
    gen1 <- gen1[,o1 ]
    gen2 <- gen2[,o2 ]
    map1 <- map1[ o1,]
    map2 <- map2[ o2,]
  }

  for (i in 1:nrow(map1)) {
    gg1 <- gen1[,i]
    ii  <- which(Map1$Marker==map1$Marker[i])
    si  <- paste(chr1,ii,sep='.p')
  
    for (j in 1:nrow(map2)) {
      gg2 <- gen2[,j]
      jj  <- which(Map2$Marker==map1$Marker[j])   
      sj  <- paste(chr2,jj,sep='.p')
  
      for (t in type) {
        flag <- 0      
        if (chr1==chr2 & (t%in%c('aa','dd')) & i>j) flag <- 1
        if (flag==0) {
          if (t=='dd' | t=='da') {
            GG1 <- gg1
            GG1[GG1==2] <- 0
            if (var(GG1)==0) flag <- 1 
          }
        }
        if (flag==0) {
          if (t=='dd' | t=='ad') {
            GG2 <- gg2
            GG2[GG2==2] <- 0
            if (var(GG2)==0) flag <- 1        
          } 
        }
        if (flag==0) {   
          fo  <- paste(fout,si,sj,t,sep='_')
          out <- Epi_Mat(gg1,gg2,t,ped,method=meth)
          write.table(out,fo,row.names=F,col.names=F,quote=F,sep='\t')
        }
      }
    }
  }
}

q('no')

