###Here, the construction of LD-correction GRM (Gld), the application of Gld to 'rrBLUP' Package for genome prediction and genetic parameter estimation are presented

#Select for a specific number of SNPs that are evenly distributed across the genome  (here 10 K SNPs are selected for example)
library(data.table)
library("doParallel")
bimfile<-c(NULL,NULL,NULL,NULL,NULL,NULL)
equ_distance_extract_snp = foreach(j = 1:29,.packages = "data.table" ) %dopar% {
  bim = fread(paste0("chr",j,".bim"))
  bim<-as.data.frame(bim)
  bp<-bim[,4]
  i<-1
  m=NULL
  n=1
  while ( i<= nrow(bim)){
  m =min(which(bp-bp[i]  >= 183000))
  n =c(n,m)
  i= m
  }
  n<-n[1:(length(n)-1)]
  rowbim<-bim[n,]
  bimfile<-rbind(bimfile,rowbim)
  }
write.table(bimfile[,2],"position",row.names=F,col.names=F,quote=F)

mv position position_10k

#creates a PLINK 1 binary fileset with the selected 10 K SNPs  
plink --bfile /WORK/scau_ljq_1/rendy/grm/ld/smat/ld_matr/ld_pruned/pruneddata --chr-set 29 no-x no-y no-xy no-mt  --extract position_10k --make-bed --out 10k

#Calculate the LD (r squared) between two SNPs on each chromosome
for i in {1..29}
do
plink --bfile 10k --chr-set 29 no-x no-y no-xy no-mt --r2  --chr "${i}" --matrix
mv plink.ld plink.ld_chr"${i}"
done


#Construction of LD-correction matrix and its inverse matrix
library("data.table")
chr1<-as.data.frame(fread("plink.ld_chr1"))
chr2<-as.data.frame(fread("plink.ld_chr2"))
chr3<-as.data.frame(fread("plink.ld_chr3"))
chr4<-as.data.frame(fread("plink.ld_chr4"))
chr5<-as.data.frame(fread("plink.ld_chr5"))
chr6<-as.data.frame(fread("plink.ld_chr6"))
chr7<-as.data.frame(fread("plink.ld_chr7"))
chr8<-as.data.frame(fread("plink.ld_chr8"))
chr9<-as.data.frame(fread("plink.ld_chr9"))
chr10<-as.data.frame(fread("plink.ld_chr10"))
chr11<-as.data.frame(fread("plink.ld_chr11"))
chr12<-as.data.frame(fread("plink.ld_chr12"))
chr13<-as.data.frame(fread("plink.ld_chr13"))
chr14<-as.data.frame(fread("plink.ld_chr14"))
chr15<-as.data.frame(fread("plink.ld_chr15"))
chr16<-as.data.frame(fread("plink.ld_chr16"))
chr17<-as.data.frame(fread("plink.ld_chr17"))
chr18<-as.data.frame(fread("plink.ld_chr18"))
chr19<-as.data.frame(fread("plink.ld_chr19"))
chr20<-as.data.frame(fread("plink.ld_chr20"))
chr21<-as.data.frame(fread("plink.ld_chr21"))
chr22<-as.data.frame(fread("plink.ld_chr22"))
chr23<-as.data.frame(fread("plink.ld_chr23"))
chr24<-as.data.frame(fread("plink.ld_chr24"))
chr25<-as.data.frame(fread("plink.ld_chr25"))
chr26<-as.data.frame(fread("plink.ld_chr26"))
chr27<-as.data.frame(fread("plink.ld_chr27"))
chr28<-as.data.frame(fread("plink.ld_chr28"))
chr29<-as.data.frame(fread("plink.ld_chr29"))

M1<-cbind(chr1,matrix(0,nrow=nrow(chr1),ncol=(ncol(chr2)+ncol(chr3)+ncol(chr4)+ncol(chr5)+ncol(chr6)+ncol(chr7)+ncol(chr8)+ncol(chr9)+ncol(chr10)+ncol(chr11)+ncol(chr12)+ncol(chr13)+ncol(chr14)+ncol(chr15)+ncol(chr16)+ncol(chr17)+ncol(chr18)+ncol(chr19)+ncol(chr20)+ncol(chr21)+ncol(chr22)+ncol(chr23)+ncol(chr24)+ncol(chr25)+ncol(chr26)+ncol(chr27)+ncol(chr28)+ncol(chr29))))
M2<-cbind(matrix(0,nrow=nrow(chr2),ncol=ncol(chr1)),chr2,matrix(0,nrow=nrow(chr2),ncol=(ncol(chr3)+ncol(chr4)+ncol(chr5)+ncol(chr6)+ncol(chr7)+ncol(chr8)+ncol(chr9)+ncol(chr10)+ncol(chr11)+ncol(chr12)+ncol(chr13)+ncol(chr14)+ncol(chr15)+ncol(chr16)+ncol(chr17)+ncol(chr18)+ncol(chr19)+ncol(chr20)+ncol(chr21)+ncol(chr22)+ncol(chr23)+ncol(chr24)+ncol(chr25)+ncol(chr26)+ncol(chr27)+ncol(chr28)+ncol(chr29))))
M3<-cbind(matrix(0,nrow=nrow(chr3),ncol=(ncol(chr1)+ncol(chr2))),chr3,matrix(0,nrow=nrow(chr3),ncol=(ncol(chr4)+ncol(chr5)+ncol(chr6)+ncol(chr7)+ncol(chr8)+ncol(chr9)+ncol(chr10)+ncol(chr11)+ncol(chr12)+ncol(chr13)+ncol(chr14)+ncol(chr15)+ncol(chr16)+ncol(chr17)+ncol(chr18)+ncol(chr19)+ncol(chr20)+ncol(chr21)+ncol(chr22)+ncol(chr23)+ncol(chr24)+ncol(chr25)+ncol(chr26)+ncol(chr27)+ncol(chr28)+ncol(chr29))))
M4<-cbind(matrix(0,nrow=nrow(chr4),ncol=(ncol(chr1)+ncol(chr2)+ncol(chr3))),chr4,matrix(0,nrow=nrow(chr4),ncol=(ncol(chr5)+ncol(chr6)+ncol(chr7)+ncol(chr8)+ncol(chr9)+ncol(chr10)+ncol(chr11)+ncol(chr12)+ncol(chr13)+ncol(chr14)+ncol(chr15)+ncol(chr16)+ncol(chr17)+ncol(chr18)+ncol(chr19)+ncol(chr20)+ncol(chr21)+ncol(chr22)+ncol(chr23)+ncol(chr24)+ncol(chr25)+ncol(chr26)+ncol(chr27)+ncol(chr28)+ncol(chr29))))
M5<-cbind(matrix(0,nrow=nrow(chr5),ncol=(ncol(chr1)+ncol(chr2)+ncol(chr3)+ncol(chr4))),chr5,matrix(0,nrow=nrow(chr5),ncol=(ncol(chr6)+ncol(chr7)+ncol(chr8)+ncol(chr9)+ncol(chr10)+ncol(chr11)+ncol(chr12)+ncol(chr13)+ncol(chr14)+ncol(chr15)+ncol(chr16)+ncol(chr17)+ncol(chr18)+ncol(chr19)+ncol(chr20)+ncol(chr21)+ncol(chr22)+ncol(chr23)+ncol(chr24)+ncol(chr25)+ncol(chr26)+ncol(chr27)+ncol(chr28)+ncol(chr29))))
M6<-cbind(matrix(0,nrow=nrow(chr6),ncol=(ncol(chr1)+ncol(chr2)+ncol(chr3)+ncol(chr4)+ncol(chr5))),chr6,matrix(0,nrow=nrow(chr6),ncol=(ncol(chr7)+ncol(chr8)+ncol(chr9)+ncol(chr10)+ncol(chr11)+ncol(chr12)+ncol(chr13)+ncol(chr14)+ncol(chr15)+ncol(chr16)+ncol(chr17)+ncol(chr18)+ncol(chr19)+ncol(chr20)+ncol(chr21)+ncol(chr22)+ncol(chr23)+ncol(chr24)+ncol(chr25)+ncol(chr26)+ncol(chr27)+ncol(chr28)+ncol(chr29))))
M7<-cbind(matrix(0,nrow=nrow(chr7),ncol=(ncol(chr1)+ncol(chr2)+ncol(chr3)+ncol(chr4)+ncol(chr5)+ncol(chr6))),chr7,matrix(0,nrow=nrow(chr7),ncol=(ncol(chr8)+ncol(chr9)+ncol(chr10)+ncol(chr11)+ncol(chr12)+ncol(chr13)+ncol(chr14)+ncol(chr15)+ncol(chr16)+ncol(chr17)+ncol(chr18)+ncol(chr19)+ncol(chr20)+ncol(chr21)+ncol(chr22)+ncol(chr23)+ncol(chr24)+ncol(chr25)+ncol(chr26)+ncol(chr27)+ncol(chr28)+ncol(chr29))))
M8<-cbind(matrix(0,nrow=nrow(chr8),ncol=(ncol(chr1)+ncol(chr2)+ncol(chr3)+ncol(chr4)+ncol(chr5)+ncol(chr6)+ncol(chr7))),chr8,matrix(0,nrow=nrow(chr8),ncol=(ncol(chr9)+ncol(chr10)+ncol(chr11)+ncol(chr12)+ncol(chr13)+ncol(chr14)+ncol(chr15)+ncol(chr16)+ncol(chr17)+ncol(chr18)+ncol(chr19)+ncol(chr20)+ncol(chr21)+ncol(chr22)+ncol(chr23)+ncol(chr24)+ncol(chr25)+ncol(chr26)+ncol(chr27)+ncol(chr28)+ncol(chr29))))
M9<-cbind(matrix(0,nrow=nrow(chr9),ncol=(ncol(chr1)+ncol(chr2)+ncol(chr3)+ncol(chr4)+ncol(chr5)+ncol(chr6)+ncol(chr7)+ncol(chr8))),chr9,matrix(0,nrow=nrow(chr9),ncol=(ncol(chr10)+ncol(chr11)+ncol(chr12)+ncol(chr13)+ncol(chr14)+ncol(chr15)+ncol(chr16)+ncol(chr17)+ncol(chr18)+ncol(chr19)+ncol(chr20)+ncol(chr21)+ncol(chr22)+ncol(chr23)+ncol(chr24)+ncol(chr25)+ncol(chr26)+ncol(chr27)+ncol(chr28)+ncol(chr29))))
M10<-cbind(matrix(0,nrow=nrow(chr10),ncol=(ncol(chr1)+ncol(chr2)+ncol(chr3)+ncol(chr4)+ncol(chr5)+ncol(chr6)+ncol(chr7)+ncol(chr8)+ncol(chr9))),chr10,matrix(0,nrow=nrow(chr10),ncol=(ncol(chr11)+ncol(chr12)+ncol(chr13)+ncol(chr14)+ncol(chr15)+ncol(chr16)+ncol(chr17)+ncol(chr18)+ncol(chr19)+ncol(chr20)+ncol(chr21)+ncol(chr22)+ncol(chr23)+ncol(chr24)+ncol(chr25)+ncol(chr26)+ncol(chr27)+ncol(chr28)+ncol(chr29))))
M11<-cbind(matrix(0,nrow=nrow(chr11),ncol=(ncol(chr1)+ncol(chr2)+ncol(chr3)+ncol(chr4)+ncol(chr5)+ncol(chr6)+ncol(chr7)+ncol(chr8)+ncol(chr9)+ncol(chr10))),chr11,matrix(0,nrow=nrow(chr11),ncol=(ncol(chr12)+ncol(chr13)+ncol(chr14)+ncol(chr15)+ncol(chr16)+ncol(chr17)+ncol(chr18)+ncol(chr19)+ncol(chr20)+ncol(chr21)+ncol(chr22)+ncol(chr23)+ncol(chr24)+ncol(chr25)+ncol(chr26)+ncol(chr27)+ncol(chr28)+ncol(chr29))))
M12<-cbind(matrix(0,nrow=nrow(chr12),ncol=(ncol(chr1)+ncol(chr2)+ncol(chr3)+ncol(chr4)+ncol(chr5)+ncol(chr6)+ncol(chr7)+ncol(chr8)+ncol(chr9)+ncol(chr10)+ncol(chr11))),chr12,matrix(0,nrow=nrow(chr12),ncol=(ncol(chr13)+ncol(chr14)+ncol(chr15)+ncol(chr16)+ncol(chr17)+ncol(chr18)+ncol(chr19)+ncol(chr20)+ncol(chr21)+ncol(chr22)+ncol(chr23)+ncol(chr24)+ncol(chr25)+ncol(chr26)+ncol(chr27)+ncol(chr28)+ncol(chr29))))
M13<-cbind(matrix(0,nrow=nrow(chr13),ncol=(ncol(chr1)+ncol(chr2)+ncol(chr3)+ncol(chr4)+ncol(chr5)+ncol(chr6)+ncol(chr7)+ncol(chr8)+ncol(chr9)+ncol(chr10)+ncol(chr11)+ncol(chr12))),chr13,matrix(0,nrow=nrow(chr13),ncol=(ncol(chr14)+ncol(chr15)+ncol(chr16)+ncol(chr17)+ncol(chr18)+ncol(chr19)+ncol(chr20)+ncol(chr21)+ncol(chr22)+ncol(chr23)+ncol(chr24)+ncol(chr25)+ncol(chr26)+ncol(chr27)+ncol(chr28)+ncol(chr29))))
M14<-cbind(matrix(0,nrow=nrow(chr14),ncol=(ncol(chr1)+ncol(chr2)+ncol(chr3)+ncol(chr4)+ncol(chr5)+ncol(chr6)+ncol(chr7)+ncol(chr8)+ncol(chr9)+ncol(chr10)+ncol(chr11)+ncol(chr12)+ncol(chr13))),chr14,matrix(0,nrow=nrow(chr14),ncol=(ncol(chr15)+ncol(chr16)+ncol(chr17)+ncol(chr18)+ncol(chr19)+ncol(chr20)+ncol(chr21)+ncol(chr22)+ncol(chr23)+ncol(chr24)+ncol(chr25)+ncol(chr26)+ncol(chr27)+ncol(chr28)+ncol(chr29))))
M15<-cbind(matrix(0,nrow=nrow(chr15),ncol=(ncol(chr1)+ncol(chr2)+ncol(chr3)+ncol(chr4)+ncol(chr5)+ncol(chr6)+ncol(chr7)+ncol(chr8)+ncol(chr9)+ncol(chr10)+ncol(chr11)+ncol(chr12)+ncol(chr13)+ncol(chr14))),chr15,matrix(0,nrow=nrow(chr15),ncol=(ncol(chr16)+ncol(chr17)+ncol(chr18)+ncol(chr19)+ncol(chr20)+ncol(chr21)+ncol(chr22)+ncol(chr23)+ncol(chr24)+ncol(chr25)+ncol(chr26)+ncol(chr27)+ncol(chr28)+ncol(chr29))))
M16<-cbind(matrix(0,nrow=nrow(chr16),ncol=(ncol(chr1)+ncol(chr2)+ncol(chr3)+ncol(chr4)+ncol(chr5)+ncol(chr6)+ncol(chr7)+ncol(chr8)+ncol(chr9)+ncol(chr10)+ncol(chr11)+ncol(chr12)+ncol(chr13)+ncol(chr14)+ncol(chr15))),chr16,matrix(0,nrow=nrow(chr16),ncol=(ncol(chr17)+ncol(chr18)+ncol(chr19)+ncol(chr20)+ncol(chr21)+ncol(chr22)+ncol(chr23)+ncol(chr24)+ncol(chr25)+ncol(chr26)+ncol(chr27)+ncol(chr28)+ncol(chr29))))
M17<-cbind(matrix(0,nrow=nrow(chr17),ncol=(ncol(chr1)+ncol(chr2)+ncol(chr3)+ncol(chr4)+ncol(chr5)+ncol(chr6)+ncol(chr7)+ncol(chr8)+ncol(chr9)+ncol(chr10)+ncol(chr11)+ncol(chr12)+ncol(chr13)+ncol(chr14)+ncol(chr15)+ncol(chr16))),chr17,matrix(0,nrow=nrow(chr17),ncol=(ncol(chr18)+ncol(chr19)+ncol(chr20)+ncol(chr21)+ncol(chr22)+ncol(chr23)+ncol(chr24)+ncol(chr25)+ncol(chr26)+ncol(chr27)+ncol(chr28)+ncol(chr29))))
M18<-cbind(matrix(0,nrow=nrow(chr18),ncol=(ncol(chr1)+ncol(chr2)+ncol(chr3)+ncol(chr4)+ncol(chr5)+ncol(chr6)+ncol(chr7)+ncol(chr8)+ncol(chr9)+ncol(chr10)+ncol(chr11)+ncol(chr12)+ncol(chr13)+ncol(chr14)+ncol(chr15)+ncol(chr16)+ncol(chr17))),chr18,matrix(0,nrow=nrow(chr18),ncol=(ncol(chr19)+ncol(chr20)+ncol(chr21)+ncol(chr22)+ncol(chr23)+ncol(chr24)+ncol(chr25)+ncol(chr26)+ncol(chr27)+ncol(chr28)+ncol(chr29))))
M19<-cbind(matrix(0,nrow=nrow(chr19),ncol=(ncol(chr1)+ncol(chr2)+ncol(chr3)+ncol(chr4)+ncol(chr5)+ncol(chr6)+ncol(chr7)+ncol(chr8)+ncol(chr9)+ncol(chr10)+ncol(chr11)+ncol(chr12)+ncol(chr13)+ncol(chr14)+ncol(chr15)+ncol(chr16)+ncol(chr17)+ncol(chr18))),chr19,matrix(0,nrow=nrow(chr19),ncol=(ncol(chr20)+ncol(chr21)+ncol(chr22)+ncol(chr23)+ncol(chr24)+ncol(chr25)+ncol(chr26)+ncol(chr27)+ncol(chr28)+ncol(chr29))))
M20<-cbind(matrix(0,nrow=nrow(chr20),ncol=(ncol(chr1)+ncol(chr2)+ncol(chr3)+ncol(chr4)+ncol(chr5)+ncol(chr6)+ncol(chr7)+ncol(chr8)+ncol(chr9)+ncol(chr10)+ncol(chr11)+ncol(chr12)+ncol(chr13)+ncol(chr14)+ncol(chr15)+ncol(chr16)+ncol(chr17)+ncol(chr18)+ncol(chr19))),chr20,matrix(0,nrow=nrow(chr20),ncol=(ncol(chr21)+ncol(chr22)+ncol(chr23)+ncol(chr24)+ncol(chr25)+ncol(chr26)+ncol(chr27)+ncol(chr28)+ncol(chr29))))
M21<-cbind(matrix(0,nrow=nrow(chr21),ncol=(ncol(chr1)+ncol(chr2)+ncol(chr3)+ncol(chr4)+ncol(chr5)+ncol(chr6)+ncol(chr7)+ncol(chr8)+ncol(chr9)+ncol(chr10)+ncol(chr11)+ncol(chr12)+ncol(chr13)+ncol(chr14)+ncol(chr15)+ncol(chr16)+ncol(chr17)+ncol(chr18)+ncol(chr19)+ncol(chr20))),chr21,matrix(0,nrow=nrow(chr21),ncol=(ncol(chr22)+ncol(chr23)+ncol(chr24)+ncol(chr25)+ncol(chr26)+ncol(chr27)+ncol(chr28)+ncol(chr29))))
M22<-cbind(matrix(0,nrow=nrow(chr22),ncol=(ncol(chr1)+ncol(chr2)+ncol(chr3)+ncol(chr4)+ncol(chr5)+ncol(chr6)+ncol(chr7)+ncol(chr8)+ncol(chr9)+ncol(chr10)+ncol(chr11)+ncol(chr12)+ncol(chr13)+ncol(chr14)+ncol(chr15)+ncol(chr16)+ncol(chr17)+ncol(chr18)+ncol(chr19)+ncol(chr20)+ncol(chr21))),chr22,matrix(0,nrow=nrow(chr22),ncol=(ncol(chr23)+ncol(chr24)+ncol(chr25)+ncol(chr26)+ncol(chr27)+ncol(chr28)+ncol(chr29))))
M23<-cbind(matrix(0,nrow=nrow(chr23),ncol=(ncol(chr1)+ncol(chr2)+ncol(chr3)+ncol(chr4)+ncol(chr5)+ncol(chr6)+ncol(chr7)+ncol(chr8)+ncol(chr9)+ncol(chr10)+ncol(chr11)+ncol(chr12)+ncol(chr13)+ncol(chr14)+ncol(chr15)+ncol(chr16)+ncol(chr17)+ncol(chr18)+ncol(chr19)+ncol(chr20)+ncol(chr21)+ncol(chr22))),chr23,matrix(0,nrow=nrow(chr23),ncol=(ncol(chr24)+ncol(chr25)+ncol(chr26)+ncol(chr27)+ncol(chr28)+ncol(chr29))))
M24<-cbind(matrix(0,nrow=nrow(chr24),ncol=(ncol(chr1)+ncol(chr2)+ncol(chr3)+ncol(chr4)+ncol(chr5)+ncol(chr6)+ncol(chr7)+ncol(chr8)+ncol(chr9)+ncol(chr10)+ncol(chr11)+ncol(chr12)+ncol(chr13)+ncol(chr14)+ncol(chr15)+ncol(chr16)+ncol(chr17)+ncol(chr18)+ncol(chr19)+ncol(chr20)+ncol(chr21)+ncol(chr22)+ncol(chr23))),chr24,matrix(0,nrow=nrow(chr24),ncol=(ncol(chr25)+ncol(chr26)+ncol(chr27)+ncol(chr28)+ncol(chr29))))
M25<-cbind(matrix(0,nrow=nrow(chr25),ncol=(ncol(chr1)+ncol(chr2)+ncol(chr3)+ncol(chr4)+ncol(chr5)+ncol(chr6)+ncol(chr7)+ncol(chr8)+ncol(chr9)+ncol(chr10)+ncol(chr11)+ncol(chr12)+ncol(chr13)+ncol(chr14)+ncol(chr15)+ncol(chr16)+ncol(chr17)+ncol(chr18)+ncol(chr19)+ncol(chr20)+ncol(chr21)+ncol(chr22)+ncol(chr23)+ncol(chr24))),chr25,matrix(0,nrow=nrow(chr25),ncol=(ncol(chr26)+ncol(chr27)+ncol(chr28)+ncol(chr29))))
M26<-cbind(matrix(0,nrow=nrow(chr26),ncol=(ncol(chr1)+ncol(chr2)+ncol(chr3)+ncol(chr4)+ncol(chr5)+ncol(chr6)+ncol(chr7)+ncol(chr8)+ncol(chr9)+ncol(chr10)+ncol(chr11)+ncol(chr12)+ncol(chr13)+ncol(chr14)+ncol(chr15)+ncol(chr16)+ncol(chr17)+ncol(chr18)+ncol(chr19)+ncol(chr20)+ncol(chr21)+ncol(chr22)+ncol(chr23)+ncol(chr24)+ncol(chr25))),chr26,matrix(0,nrow=nrow(chr26),ncol=(ncol(chr27)+ncol(chr28)+ncol(chr29))))
M27<-cbind(matrix(0,nrow=nrow(chr27),ncol=(ncol(chr1)+ncol(chr2)+ncol(chr3)+ncol(chr4)+ncol(chr5)+ncol(chr6)+ncol(chr7)+ncol(chr8)+ncol(chr9)+ncol(chr10)+ncol(chr11)+ncol(chr12)+ncol(chr13)+ncol(chr14)+ncol(chr15)+ncol(chr16)+ncol(chr17)+ncol(chr18)+ncol(chr19)+ncol(chr20)+ncol(chr21)+ncol(chr22)+ncol(chr23)+ncol(chr24)+ncol(chr25)+ncol(chr26))),chr27,matrix(0,nrow=nrow(chr27),ncol=(ncol(chr28)+ncol(chr29))))
M28<-cbind(matrix(0,nrow=nrow(chr28),ncol=(ncol(chr1)+ncol(chr2)+ncol(chr3)+ncol(chr4)+ncol(chr5)+ncol(chr6)+ncol(chr7)+ncol(chr8)+ncol(chr9)+ncol(chr10)+ncol(chr11)+ncol(chr12)+ncol(chr13)+ncol(chr14)+ncol(chr15)+ncol(chr16)+ncol(chr17)+ncol(chr18)+ncol(chr19)+ncol(chr20)+ncol(chr21)+ncol(chr22)+ncol(chr23)+ncol(chr24)+ncol(chr25)+ncol(chr26)+ncol(chr27))),chr28,matrix(0,nrow=nrow(chr28),ncol=(ncol(chr29))))
M29<-cbind(matrix(0,nrow=nrow(chr29),ncol=(ncol(chr1)+ncol(chr2)+ncol(chr3)+ncol(chr4)+ncol(chr5)+ncol(chr6)+ncol(chr7)+ncol(chr8)+ncol(chr9)+ncol(chr10)+ncol(chr11)+ncol(chr12)+ncol(chr13)+ncol(chr14)+ncol(chr15)+ncol(chr16)+ncol(chr17)+ncol(chr18)+ncol(chr19)+ncol(chr20)+ncol(chr21)+ncol(chr22)+ncol(chr23)+ncol(chr24)+ncol(chr25)+ncol(chr26)+ncol(chr27)+ncol(chr28))),chr29)

cname<-paste0("snp",1:10113)
#cname<-noquote("cname")
write.table(cname,"cname",row.names=F,col.names=F,quote=F)
cname<-read.table("cname")[,1]
colnames(M1)<-cname
colnames(M2)<-cname
colnames(M3)<-cname
colnames(M4)<-cname
colnames(M5)<-cname
colnames(M6)<-cname
colnames(M7)<-cname
colnames(M8)<-cname
colnames(M9)<-cname
colnames(M10)<-cname
colnames(M11)<-cname
colnames(M12)<-cname
colnames(M13)<-cname
colnames(M14)<-cname
colnames(M15)<-cname
colnames(M16)<-cname
colnames(M17)<-cname
colnames(M18)<-cname
colnames(M19)<-cname
colnames(M20)<-cname
colnames(M21)<-cname
colnames(M22)<-cname
colnames(M23)<-cname
colnames(M24)<-cname
colnames(M25)<-cname
colnames(M26)<-cname
colnames(M27)<-cname
colnames(M28)<-cname
colnames(M29)<-cname

s<-rbind(M1,M2,M3,M4,M5,M6,M7,M8,M9,M10,M11,M12,M13,M14,M15,M16,M17,M18,M19,M20,M21,M22,M23,M24,M25,M26,M27,M28,M29)
si<-solve(s)

my_list <- rep("",29)
for(i in 1:length(my_list)) {assign(paste0("M", i), my_list[i])}
for(i in 1:length(my_list)) {assign(paste0("chr", i), my_list[i])}
gc()

fwrite(x = data.table::data.table(si),file = ("si_10k"),col.names = FALSE,row.names = FALSE,quote=F)

#Creat SNP genotype file, SNP genotypes are coded as a single Allele dosage number of additive components
plink --bfile 10k --chr-set 29 no-x no-y no-xy no-mt --recodeA --out A_10k


#Construction of LD-correction GRM
library("data.table")
s<-as.data.frame(fread("si_10k"))
A<-as.data.frame(fread("A_10k.raw",head=T))
A<-A[,-c(1:6)]
scaled.A <- scale(A,scale=FALSE)
a<-scaled.A
a<-as.matrix(a)
s<-as.matrix(s)
g0<-a%*%s%*%t(a)
p<-colMeans(A)/2
scale<-2*t(p)%*%s%*%(1-p)
g<-(1/scale[1,1])*g0
fwrite(x = data.table::data.table(g),file = ("g_10k"),col.names = FALSE,row.names = FALSE,quote=F)


#LD-correction GRM based GBLUP
library("data.table")
library(rrBLUP)
A <- fread("g_10k")
fam<-as.data.frame(fread("10k.fam"))
id<-fam[,2]
id<-as.character(id)
rownames(A)<-id
colnames(A)<-id
A<-as.matrix(A)
rownames(A)<-id
mean(diag(A))
dat<-read.table("sim.pheno")
colnames(dat)<-c("gid","id","y")
datt<-dat[1:1800,]
ans <- kin.blup(data=datt,geno="id",pheno="y",K=A)
cor(dat$y[1801:2000],ans$g[1801:2000])
cor(dat$y[1:1800],ans$g[1:1800])
ans$Vg/(ans$Vg+ans$Ve)


#VanRaden G based GBLUP
M<-as.data.frame(fread("A_10k.raw",head=T))
M<-M[,-c(1:6)]
rownames(M) <-id
m<-M-1
A <- A.mat(m)
mean(diag(A))
ans <- kin.blup(data=datt,geno="id",pheno="y",K=A)
cor(dat$y[1801:2000],ans$g[1801:2000])
cor(dat$y[1:1800],ans$g[1:1800])
ans$Vg/(ans$Vg+ans$Ve)
