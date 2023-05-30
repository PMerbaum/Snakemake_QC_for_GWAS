hets<-read.table(snakemake@input[[1]], header=TRUE) #calculates 3 SD's from the het mean and excludes all the rest 
stdev=sd(hets$F)
stdev3=sd(hets$F)*3
avg=mean(hets$F)

F_min<- avg-stdev3
F_plus<- avg+stdev3
hets_good<-subset(hets, hets$F>F_min & hets$F <F_plus)
hets_fail<-subset(hets, hets$F<F_min | hets$F >F_plus)

hets_fail2=hets_fail[,1:2]
write.table(hets_fail2, file=snakemake@output[[1]], row.names=FALSE, quote=FALSE, sep="\t")