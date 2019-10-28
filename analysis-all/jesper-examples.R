
#hitta på lite p-värden och fc och symboler

#tusen random värden enligt normalfördelningen
all <- matrix(NA, nrow=1000, ncol=100)
for (i in 1:1000){
  if (i%%100==0){
     all[i,1:50] <- rnorm(50, mean=0, sd=1)
     all[i,51:100] <- rnorm(50, mean=1, sd=1)
  }else{
     all[i,] <- rnorm(100, mean=0, sd=1)
  }
}
#gör två grupper
grp1 <- all[,1:50]
grp2 <- all[,51:100]

#gör t.test
p <- unlist(lapply(1:1000,function(i,x,y){t.test(x[i,],y[i,])$p.value}, x=grp1, y=grp2))
fc <- unlist(lapply(1:1000,function(i,x,y){mean(x[i,], na.rm=TRUE) - mean(y[i,], na.rm=TRUE)}, x=grp1, y=grp2))

symbols <- paste("gene",1:1000,sep="")

#select genes to label
symbolvec <- c("gene23","gene143","gene2","gene1","gene100", "gene200")

#try function
setwd("/Users/hanbjo/Desktop/2019-08-21-OLINK/repo/jesper")
source("functions-all/2019-09-05-funktionsfil.R")

#VOLCANO TIME
Volvano_plot_pro(p,fc,symbols,symbolvec, basecutoff=1)


