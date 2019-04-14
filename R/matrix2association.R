# 0. Librairies needed (to install before if not)
library(grid)
library(ggplot2)
library(ggbeeswarm)
source("snp_distribution_plot.r")
source("multiplot.r")

# 1. Path : To ADAPT !!!!!!!!!!

setwd("C:/Users/l.guillier/Desktop/Projets/COMPARE")

# 2. Read tsv file of snp pairwise distance
distances_matrix_1<-read.table("cgMLST_.tsv")

# 3. Define group of related strains
group<-c("FI_2017_HUMAN18",
"FI_2017_HUMAN26",
"FR_2017_FOOD8",
"FI_2017_HUMAN19",
"FI_2017_HUMAN20",
"FI_2017_HUMAN21",
"FI_2017_HUMAN22",
"FI_2017_HUMAN23",
"FI_2017_HUMAN24",
"FI_2017_HUMAN25",
"AT_2016_FOOD16","AT_2017_FOOD17")


# 4. Matrix to list function (do not modify !!)
matrix_to_list<-function(distances_matrix)
{li<-1
distances_list<-c()
list_strains_row<-rownames(distances_matrix)
list_strains_col<-list_strains_row
strain1<-c()

for (col in 1:length(distances_matrix)-1)
{for (line in li:length(distances_matrix))              
{distances_list<-c(distances_list,distances_matrix[col,line])
#print(distances_matrix[col,line])
strain1<-c(strain1,list_strains_col[col])}
  li<-li+1}

strain2<-c()
li<-0
for (col in 1:c(length(distances_matrix)-1))
{li<-col+1
#print(li)
  while (li<length(distances_matrix)+1)              
{strain2<-c(strain2,list_strains_col[li])
li<-li+1}
}

list_pairs<-data.frame(s1=strain1,s2=strain2,d=distances_list)
return(list_pairs)}

# 5. Get lists for comparison (do not modify !!)
distances_list_1<-matrix_to_list(distances_matrix_1)

write.csv2(distances_list_1,"list_pairwise_snp.csv",row.names = FALSE)

# 6. Sort Strains of interest (do not modify !!)
strain_names<-rownames(distances_matrix_1)


pos<-matrix(NA,length(group),length(strain_names))
for (i in 1:c(length(group)))
{
  pos1<-c(grep(group[i],distances_list_1$s1))
pos2<-c(grep(group[i],distances_list_1$s2))
pos12<-c(pos1,pos2)
pos[i,1:length(pos12)]<-pos12
#pos1<-c()
#pos2<-c()
}

pairs<-c()
for (i in 1:c(length(group))) { j<-i+1;
while (j<length(group)+1){ pair<-intersect(pos[i,],pos[j,]) 
  pair<-pair[!is.na(intersect(pos[i,],pos[j,]))] 
  pairs<-c(pairs,pair)
  j<-j+1  }
}

pairs_control<-distances_list_1$d[pairs]

# 6. Test each candidate strains  (do not modify !!)
result1<-rep(NA,length(strain_names))
result2<-rep(NA,length(strain_names))
result3<-rep(NA,length(strain_names))
bin_with<-2  # To be adapted by the user

for (i in 1:length(strain_names))
{ pairs_c<-c() 
test<-c()
  if (length((which(strain_names[i]!=group)))==length(group))
    {pos1<-c(grep(strain_names[i],distances_list_1$s1))
pos2<-c(grep(strain_names[i],distances_list_1$s2))
pos12<-c(pos1,pos2)
pos_candidate<-pos12
j<-1;
while (j<length(group)+1){ pair_c<-intersect(pos_candidate,pos[j,]) 
pair_c<-pair_c[!is.na(intersect(pos_candidate,pos[j,]))] 
pairs_c<-c(pairs_c,pair_c)
j<-j+1  }
pairs_candidate<-distances_list_1$d[pairs_c]
limits<-max(pairs_control)   # to be adapted by the user
pa<-snp_distribution_plot1(pairs_control,bin_with,limits,max(pairs_candidate,pairs_control)+5)
pb<-snp_distribution_plot1(pairs_candidate+0.001,bin_with,limits,max(pairs_candidate,pairs_control)+5)
pdf(file=paste(strain_names[i], '.pdf', sep=''))
multiplot(pa, pb,cols=1)
dev.off()

test1<-wilcox.test(pairs_candidate,pairs_control)
test2<-ks.test(pairs_control,pairs_candidate)
test3<-kruskal.test(list(pairs_control,pairs_candidate))
print(strain_names[i])
print(pairs_candidate)
print(pb)
result1[i]<-test1$p.value
result2[i]<-test2$p.value
result3[i]<-test3$p.value
}

#else {result[i]<-c(NA)}
}

associated_strains<-strain_names[which(result1>0.01)]
pvalue_associated<-result1[which(result1>0.01)]
associated<-cbind.data.frame(associated_strains,pvalue_associated)
write.csv2(associated,"associated_w.csv")


nonassociated_strains<-strain_names[which(result1<0.01)]
pvalue_nonassociated<-result1[which(result1<0.01)]
nonassociated<-cbind.data.frame(nonassociated_strains,pvalue_nonassociated)
write.csv2(nonassociated,"nonassociated_w.csv")


associated_strains<-strain_names[which(result2>0.01)]
pvalue_associated<-result2[which(result2>0.01)]
associated<-cbind.data.frame(associated_strains,pvalue_associated)
write.csv2(associated,"associated_ks.csv")


nonassociated_strains<-strain_names[which(result2<0.01)]
pvalue_nonassociated<-result2[which(result2<0.01)]
nonassociated<-cbind.data.frame(nonassociated_strains,pvalue_nonassociated)
write.csv2(nonassociated,"nonassociated_ks.csv")


associated_strains<-strain_names[which(result3>0.01)]
pvalue_associated<-result3[which(result3>0.01)]
associated<-cbind.data.frame(associated_strains,pvalue_associated)
write.csv2(associated,"associated_kw.csv")


nonassociated_strains<-strain_names[which(result3<0.01)]
pvalue_nonassociated<-result3[which(result3<0.01)]
nonassociated<-cbind.data.frame(nonassociated_strains,pvalue_nonassociated)
write.csv2(nonassociated,"nonassociated_kw.csv")


# 7. Global representation  (do not modify !!)
alpha=0.01
result1<-rep(NA,length(strain_names))
result2<-rep(NA,length(strain_names))
result3<-rep(NA,length(strain_names))

pairs_candidate_conca<-c()
id_conca<-c()
p_values_conca1<-c()
p_values_conca2<-c()
p_values_conca3<-c()
asso_conca1<-c()
asso_conca2<-c()
asso_conca3<-c()

for (i in 1:length(strain_names))
{ pairs_c<-c() 
test<-c()
if (length((which(strain_names[i]!=group)))==length(group))
{pos1<-c(grep(strain_names[i],distances_list_1$s1))
pos2<-c(grep(strain_names[i],distances_list_1$s2))
pos12<-c(pos1,pos2)
pos_candidate<-pos12
j<-1;
while (j<length(group)+1){ pair_c<-intersect(pos_candidate,pos[j,]) 
pair_c<-pair_c[!is.na(intersect(pos_candidate,pos[j,]))] 
pairs_c<-c(pairs_c,pair_c)
j<-j+1  }
pairs_candidate<-distances_list_1$d[pairs_c]
pairs_candidate_conca<-c(pairs_candidate_conca,pairs_candidate)

id_conca<-c(id_conca,rep(strain_names[i],length(pairs_candidate)))

test1<-wilcox.test(pairs_candidate,pairs_control)
test2<-ks.test(pairs_control,pairs_candidate)
test3<-kruskal.test(list(pairs_control,pairs_candidate))
result1[i]<-test1$p.value
result2[i]<-test2$p.value
result3[i]<-test3$p.value

p_values_conca1<-c(p_values_conca1,rep(test1$p.value,length(pairs_candidate)))
p_values_conca2<-c(p_values_conca2,rep(test2$p.value,length(pairs_candidate)))
p_values_conca3<-c(p_values_conca3,rep(test3$p.value,length(pairs_candidate)))

if (result1[i]<alpha)
{asso_conca1=c(asso_conca1,rep("unlinked",length(pairs_candidate)))} else {asso_conca1=c(asso_conca1,rep("linked",length(pairs_candidate)))}

if (result2[i]<alpha)
{asso_conca2=c(asso_conca2,rep("unlinked",length(pairs_candidate)))} else {asso_conca2=c(asso_conca2,rep("linked",length(pairs_candidate)))}

if (result3[i]<alpha)
{asso_conca3=c(asso_conca3,rep("unlinked",length(pairs_candidate)))} else {asso_conca3=c(asso_conca3,rep("linked",length(pairs_candidate)))}


}
}

data1<-data.frame(id=id_conca,distance=pairs_candidate_conca,p=p_values_conca1,link=asso_conca1)
data2<-data.frame(id=id_conca,distance=pairs_candidate_conca,p=p_values_conca2,link=asso_conca2)
data3<-data.frame(id=id_conca,distance=pairs_candidate_conca,p=p_values_conca3,link=asso_conca3)

#write.csv("concatenate.csv",data2)
#data3<-read.csv("concatenate.csv",stringsAsFactors = FALSE)

#par( mar=c(5, 5, 2.5, 1)) 
#(distance ~ id, data=data, col=sample(colors(), 27),ylab='',xlab='distance',vertical=FALSE, method="swarm", cex=0.2,pch=19,cex.axis=0.5)
#beeswarm(distance ~ id, data=data, col="grey60",xlab='',ylab='distance',vertical=TRUE, method="swarm", cex=0.2,pch=19,cex.axis=0.5,las=2)

pdf(file="wilcoxon.pdf")
ggplot(data1,aes(distance, id,color=link)) + geom_quasirandom()
dev.off()

pdf(file="kolmogorov.pdf")
ggplot(data2,aes(distance, id,color=link)) + geom_quasirandom()
dev.off()

pdf(file="kruskal.pdf")
ggplot(data3,aes(distance, id,color=link)) + geom_quasirandom()
dev.off()  
  
