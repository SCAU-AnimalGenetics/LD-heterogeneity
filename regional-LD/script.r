#GCTA software was used to calculate LD score, LD score is defined as the sum of LD r2 between a variant and all the variants in a region (10000 Kb was used in here).
gcta64 --bfile clean --autosome-num 29 --ld-score --ld-wind 10000 --ld-rsq-cutoff 0.01 --out ld_score
#Output file: ld_score.score.ld
SNP chr bp MAF mean_rsq snp_num max_rsq ldscore
1_100260 1 100260 0.177 0.137505 1046 1 144.83
1_120183 1 120183 0.15675 0.139788 1101 1 154.907
1_146011 1 146011 0.15675 0.139788 1101 1 154.907
1_147231 1 147231 0.15675 0.139788 1101 1 154.907
#we used the mean LD Score to represent the regional LD. The mean LD Score was fitted by segments with an average length of 100 Kb using a sliding window approach
library("data.table")
ldscore<-as.data.frame(fread("ld_score.score.ld"))
het<-NULL
for (j in 1:29)
{
chr<-ldscore[which(ldscore[,2]==j),]
mn<-NULL
for(i in chr[,3])
{
a<-chr[which(chr[,3]>=i-50000&chr[,3]<=i+50000),8]
b<-mean(a)
mn<-c(mn,b)
}
het<-c(het,mn)
}
ld_het<-cbind(ldscore,het)
colnames(ld_het)[9]<-"het"
write.table(ld_het,"ld_het",row.names=F,quote=F)
#Output
SNP chr bp MAF mean_rsq snp_num max_rsq ldscore het
1_100260 1 100260 0.177 0.137505 1046 1 144.83 152.8916
1_120183 1 120183 0.15675 0.139788 1101 1 154.907 153.990909090909
1_146011 1 146011 0.15675 0.139788 1101 1 154.907 153.990909090909
#The last column (het) is the mean LD Score

