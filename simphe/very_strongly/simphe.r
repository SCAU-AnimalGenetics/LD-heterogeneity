#Simulation of phenotypes controlled by very strongly tagged causal variants. Phenotypic heritability was 0.8 and QTL number was 100
library("data.table")
#Read the SNP genotype file “aseq.raw”, SNP genotypes are coded as a single allele dosage number of additive components, this SNP genotype file was generated using PLINK --recodeA
geno = as.data.frame(fread("/WORK/scau_ljq_1/rendy/density/cattle/data/qc/maf0.01/aseq.raw"))
#Remove columns 1, 3, 4, 5, and 6 of ASeq.raw, keep the second column (individual ID)
cleang<-geno[,-c(1,3,4,5,6)]
rownames(cleang) = cleang[,1]
#Name columns
colname<-c("ID",paste0("SNP", c(1:336977)))
colnames(cleang)<- colname
id<-read.table("/WORK/scau_ljq_1/rendy/grm/simulation/simphe/id")
clean<-cleang[match(id[,1], cleang[,1]),]
cleang<-clean
cleang1<-cleang[,-1]
#Read "weights" file, "weights" is calculated using the --calc-weights-all of LDAK software, where the column named "Eff_Tagging" is the tagging condition of each markers
weights<-as.data.frame(fread("/WORK/scau_ljq_1/rendy/density/cattle/data/qc/maf0.01/weights"))
#Sort the markers from strong to weak tagging
weights[,7]<-c(1:336977)
colnames(weights)[7] <- "SNP_Index"
weights<-weights[order(-weights$Eff_Tagging), ]
#The first 1/5 SNPs are labeled as very_strong tagging, the first 2/5 SNPs are labeled as strong tagging, the last 2/5 SNPs are labeled as weak tagging, and the last 1/5 SNPs are labeled as very_weak tagging
very_strong<-weights[1:67395,7]
strong<-weights[1:134790,7]
weak<-weights[202187:336977,7]
very_weak<-weights[269582:336977,7]
#calculate the variance of each SNP genotype (the true variance)
genvar <- apply(cleang1, 2, var)
#Calculate the frequency of the given alleles for each SNP
p<-apply(cleang1, 2, mean)/2
#calculate the variance of each SNP genotype (the expected variance, suppose the population is in hardy-Weinberg equilibrium)
mult <- 2*p*(1-p)
diff<-mult-genvar
data<-cbind(c(1:336977),diff)
colnames(data)<- c("position","Difference_of_variance")
chip_position<-read.table("/WORK/scau_ljq_1/rendy/density/cattle/data/qc/maf0.01/chip_position")
data<-data[which(data[,1]%in%chip_position[,1]),]
#Filter loci based on differences between True and Expected variance
data2<-data[which(abs(data[,2])<0.005),]
#Simulation of phenotypes controlled by very strongly tagged causal variants, 100 QTL were randomly selected on the genome
locinum<-sample(data2[data2[,1]%in%very_strong,1],100,replace=FALSE)
qtlpos<-locinum[order(locinum)]
#QTL genotypes matrix
QTLloci<-cleang1[,qtlpos]
#Calculate the frequency of the given alleles for each QTL
p<-apply(QTLloci, 2, mean)/2
#calculate the variance of each QTL genotype
mult <- 2*p*(1-p)
qtl<-as.matrix(QTLloci)
#Centralize each column
qtlcen<-sweep(qtl, 2, colMeans(qtl))
#Simulation of allele substitution effects of QTL
eff<-sqrt((0.8*(mult^(-1)))/(length(mult)))
#Half of the QTL allele substitution effects were set as negative and half as positive, so that the mean breeding value of the population is zero
ef<-cbind(1:length(eff),eff)
rows<-sample(1:length(eff),length(eff)/2,replace=FALSE)
minus<-ef[rows,]
minus[,2]<-minus[,2]*-1
effe<-rbind(ef[-rows,],minus)
effec<-effe[order(effe[,1]),2]
#Simulated the real breeding values of individuals
tbv<-qtlcen%*%effec
snpvar<-mult*(effec^2)
res<-rnorm(2000,mean=0,sd=sqrt(0.2))
#Simulate residual effects
re<-sqrt((0.2*var(tbv))/(0.8*var(res)))[1]*(res-rep(mean(res),2000))
#Simulate phenotype
phe<-tbv+re
indinfo<-cbind(cleang[,1],tbv,re,phe)
qtlinfo<-cbind(qtlpos,p,mult,effec,snpvar)
write.table(qtlinfo,"qtlinfo",row.names=FALSE,col.names=FALSE,quote=F)
write.table(indinfo,"indinfo",row.names=FALSE,col.names=FALSE,quote=F)
vid<-sample(cleang[,1],200,replace=FALSE)
v.pheno<-indinfo[which(indinfo[,1]%in%vid),c(1,1,4)]
t.pheno<-indinfo[-which(indinfo[,1]%in%vid),c(1,1,4)]
write.table(t.pheno,"t.pheno",row.names=FALSE,col.names=FALSE,quote=F)
write.table(v.pheno,"v.pheno",row.names=FALSE,col.names=FALSE,quote=F)
