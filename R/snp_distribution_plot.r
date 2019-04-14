snp_distribution_plot1<-function(distances_list,bin_width,limits,max)
{data<-data.frame(distances=distances_list)
library(ggplot2)
taille_histo<-bin_width

g1<-ggplot(data,aes(distances))+geom_histogram(origin = 0,binwidth = 1) +
  geom_vline(xintercept=c(limits),size=0.5,color='red',linetype = 20)
g1<-g1+xlim(c(0, max))
return(g1)
}

# Within [0,max_range], that is, a zoom
snp_distribution_plot2<-function(distances_list,bin_width,limits,max_range)
{data<-data.frame(distances=distances_list)
library(ggplot2)
taille_histo<-bin_width

data50<-data.frame(distances=distances_list[which(distances_list<max_range)])
g2<-ggplot(data50,aes(snp_d))+geom_histogram(binwidth=taille_histo)+
  geom_vline(xintercept=c(limits[which(limits<max_range)]),size=0.5,color='red',linetype = 20)
return(g2)
}