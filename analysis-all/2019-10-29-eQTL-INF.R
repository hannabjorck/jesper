#ställ dig i rätt mapp (användarspecifik, så en annan mapp på jespers system)
setwd("/Users/hanbjo/Desktop/2019-08-21-OLINK/repo/jesper")


#where is data stored (relative path)
datastorage <- "../../data/"
scriptsFunctions <- "functions-all"

#läs in SNP data
#header=TRUE tar bort V:na och gör så att första kolumnen blir heading

#df_short <- read.table("../../data/191025_CHB_1E04_list.txt", header=TRUE, stringsAsFactors=FALSE) KORTARE LISTA (50 SNP ca)
df <- read.table("../../data/191025_CHB_5E04_list.txt", header=TRUE, stringsAsFactors=FALSE)

#Ladda paketet GenomicRanges
library("GenomicRanges")

#ploppa in all information du behöver/som finns
#Ojekt gr som innehåller nu alla positioner för de varianter du vill ha. 
#Detta kan man använda som input till readVcf som då kan hämta exakt de raderna från VCF filen
  
gr <- GRanges(seqnames=df$CHR, IRanges(start=df$BP,end=df$BP), strand="*") 
#gr <- GRanges(seqnames=paste("chr",df$CHR,sep=""), IRanges(start=df$BP,end=df$BP), strand="*")

    
#SNPs åker in som namn    
names(gr) <- df$SNP

#gr2 <- gr[seqnames(gr)=="22"]

#Ladda paketet VariantAnnotation
library("VariantAnnotation")
#------------------------------------------------------------------------------------------

#En chr för att testa
#fl <- "../../data/vcfs/ASMADA_ASAP-MASAP_1KG_chr22.dose_maf-0.001_rsq-0.3-overlapping-variants.vcf.gz"

#LOOPA sen iställt in alla chr samtidigt


  lstVcf <- list()
      for (i in 1:22){
        cat("iteration number", i, "\n")
        gr2 <- gr[seqnames(gr)==i] 
        fl <- paste("../../data/vcfs/ASMADA_ASAP-MASAP_1KG_chr", i, ".dose_maf-0.001_rsq-0.3-overlapping-variants.vcf.gz", sep="")

#-------------------------------------------------------------------------------------------      

      #param <- ScanVcfParam(which=gr2)
      #vcf <- readVcf(fl, "hg19", param=gr2)
    
      vcf <- readVcf(fl, "hg19", param=ScanVcfParam(which=gr2))
      lstVcf[[i]] <- vcf   

      }

vcfAll <- do.call("rbind", lstVcf)
vcf <- vcfAll

#Läs in genexpressionsdata
heartexp <- read.table("../../data/2018-06-25-bed-expression-ASAP_H.bed", header=TRUE, comment.char = "&", stringsAsFactors=FALSE)

#Gör ett GRanges objekt och ploppa in informationen.
gr_exp <- GRanges(seqnames=heartexp$X.Chr, IRanges(start=heartexp$start,end=heartexp$end), strand="*") 

#Ladda paketet SummarizedExperiment
library("SummarizedExperiment")

#Ta ut en rad (dvs. SNP) med dess genotyper (dubbelbrackets [[]] om du vill få ut ngt frpn en lista)
gt <- geno(vcf)[["GT"]]

#för att visualisera på skärmen/se vad du har, plocka ut gt för SNPs (rad 1-4 och 6), och individer (col 1-5). 
#gt[c(1:4,6), 1:5]

#Subset heartexpdatan och ta bort kol du inte behöver 
heartexp_clean <- heartexp[ ,grep("A", colnames(heartexp))]

#omvandla datafram (heartexp_clean) till matris med kommandot
counts=as.matrix(heartexp_clean)

inplace <- DataFrame(row.names=colnames(counts))

#Skapa SE som sedan håller koll på alla vektorer
se <-SummarizedExperiment(assays=list(counts=as.matrix(heartexp_clean)), rowRanges=gr_exp, colData=inplace)


#om du i framtiden vill lägga in exp data från alla vävnade, använd
#assays=list(counts1=as.matrix(heartexp_clean), counts2=as.matrix(liverexp_clean))
#och döp counts1 till heart,osv..


#Ta ut en random SNP och förvandlavärdena till 1,2,3
ranSNP <- geno(vcf)[["GT"]][1,]


#FÄRGLÄGG prickarna enligt deras INF-score. Börja med att importera data
infscore <- read.table("../../data/2019-10-25-INF score.txt", header=TRUE, stringsAsFactors=FALSE)

#Gör om data/ID så att det matchar utseendemässigt
id_infscore <- as.character(infscore$ID)
id_infscore2 <- paste("A", infscore[ ,"ID"], sep="")

#finns alla x i y, och vice versa
#tf2 <- id_se %in% id_infscore2 & id_se %in% id_vcf2

id_se <- colnames(se)

tf <- id_infscore2 %in% id_se
tf2 <- id_se %in% id_infscore2

#plocka ut ALL ID-data som överlappar.
infscore[ ,"ID"] <- id_infscore2
se2 <- se[,tf2]

#matcha så att ID ligger i samma ordning
ID1 <- infscore$ID
m <- match(colnames(se2), infscore$ID)

tf <- infscore[m,]$ID==id_se
if(!all(tf)) stop("something is wrong")
infscore_m <-infscore[m,]

#ret <- cbind(infscore2_m,id_se)

#lägg in infscore under INF
colData(se2)[["INF"]] <- infscore_m[,2]

###########################
# Här matcha vcf-informationen
###########################

#fixa till dubbelnamnen i vcf ID
id_vcf <- colnames(vcf)
id_vcf2 <- sub(".*_","", id_vcf)
#byt ut de gamla id mot dessa nya i id_vcf2
colnames(vcf) <- id_vcf2

#finns alla x i y, och vice versa
id_se2 <- colnames(se2)
tf <- id_vcf2 %in% id_se2
tf2 <- id_se2 %in% id_vcf2

se3 <- se2[,tf2]

#matcha på så sätt att individerna är i samma ordning
vcf2 <- vcf[,tf]
m <- match(colnames(se3), colnames(vcf2))

tf <- colnames(vcf2)[m]==colnames(se3)
if(!all(tf)) stop("something is wrong") 
vcf2m <- vcf2[,m]

#vi har genexpression, vi har genotyp, alla individer är matchade, DÄRFÖR eQTL
#Ta ut en SNP från den matchade datan
snp1 <- geno(vcf2m)[["GT"]][1,]

#skapa en tom vektor där du stoppar in alla värden, snp1X
snp1X <- rep(NA, length(snp1))

#lägg ihop heterozygoterna i samma grupp, .
snp1X[snp1=="0|0"] <- 1
snp1X[snp1=="1|0" | snp1=="0|1" ] <- 2
snp1X[snp1=="1|1"] <- 3

#lm gillar inte characters, gör därför värdena till numeric innan du gör lm
snp1X <- as.numeric(snp1X)

#Ta ut en gen från den matchade datan
gene1 <- assays(se3)[["counts"]][1,]
#Ta ut INF från den matchade datan, vi ska använda den för att göra färgkaos
inf <- colData(se3)[["INF"]] 
colData(se3)$INF <- sub(",",".",colData(se3)$INF)
inf<- as.numeric(colData(se3)$INF)
as.numeric(colData(se3)$INF)

#gör färgkaoset. Först, ladda paketet och "set colours using RColorBrewer"
library(RColorBrewer)
cols = brewer.pal(4, "Blues")
# Define colour pallete
pal = colorRampPalette(c("white", "red"))

# Rank variable for colour assignment
order = findInterval(inf, sort(inf))
order = findInterval(as.numeric(colData(se3)$INF), sort(as.numeric(colData(se3)$INF)))
pal(length(order))[order]

#plotta. SNP måste göras om till factor, för att ta bort snuffarna. EDIT: snp1X måste göras om till as.integer för att få scatter plot. 
library(gplots)
better <- space(as.integer(snp1X), gene1, s=1/30, na.rm=TRUE, direction="x")
x <- better$x
y <- better$y

#plot(as.integer(snp1X), gene1, main="SNP1", xlab= "gt", ylab="expression", pch=19, col=pal(length(order))[order] )
plot(x, y, main="SNP1", xlab= "gt", ylab="expression", pch=19, col=pal(length(order))[order] )



#--------------------------------------------------------------------------------------
#analysera
res1 <- lm(gene1~snp1X)

#För att kolla alla gener kring en snp. HANNA, KOMMENTERA HÄR!!
#giantSNP <- granges(vcf2m[1,])+200000
hits <- findOverlaps(granges(vcf2m[1,])+200000, se2)

# --> Hits object with 2 hits and 0 metadata columns, betyder att i se2 så matchade element 14255 och 14256 (i detta fall två gener som "träffar" giantSNP)
#plocka ut dessa två hits och spara dem i ett nytt objekt (se4)
se4 <- se3[subjectHits(hits)]
#--------------------------------------------------------------------------------------

#for loop
#TOM VEKTOR
lst <- list()


  #Ta ut alla GT från den matchade datan
  for (i in 1:nrow(vcf2m)){

    snp_each <- geno(vcf2m)[["GT"]][i, ]
    res <- rep(NA, length(snp_each))
    
    #lägg ihop heterozygoterna i samma grupp
    res[snp_each=="0|0"] <- 1
    res[snp_each=="1|0" | snp1=="0|1" ] <- 2
    res[snp_each=="1|1"] <- 3
    
    #lm gilla rinte characters, gör därför värdena till numeric 
    # res <- as.integer(snp_each)
    
    #För att kolla alla gener kring en snp. 
    hits <- findOverlaps(granges(vcf2m[i,])+200000, se3)
    
    # plocka ut hits och spara dem i ett nytt objekt (se3)
    se4 <- se3[subjectHits(hits)]
    
    
        #en till for loop
        res_gene <- rep(NA, nrow(se4))
        
       #pdf("../../out/eQTL_test.pdf")
        
        filename <- paste("../../out/eQTL_test2_", i, ".pdf", sep="" )
        pdf(filename)
        
        for(j in 1:nrow(se4)){
    
        #Ta ut en gen från den matchade datan
        gene <- assays(se4[j,])[["counts"]]
        gene <- as.numeric(gene)
        
        #Ta ut en INF från den matchade datan
        colData(se4)$INF <- sub(",",".",colData(se4)$INF)
        inf <- colData(se4)[["INF"]]
        inf<- as.numeric(colData(se4)$INF)
        
       
        #gör färgkaoset. Först, ladda paketet och "set colours using RColorBrewer"
        library(RColorBrewer)
        cols = brewer.pal(4, "Blues")
        # Define colour pallete
        pal = colorRampPalette(c("beige", "red"))
        
        #Rank variable for colour assignment
        #order = findInterval(inf, sort(inf))
        order = findInterval(as.numeric(colData(se4)$INF), sort(as.numeric(colData(se4)$INF)))
        
        
        #plotta. SNP måste göras om till factor, för att ta bort snuffarna. EDIT: snp1X måste göras om till as.integer för att få scatter plot. 
        library(gplots)
        better <- space(as.integer(res), gene, s=1/30, na.rm=TRUE, direction="x")
        x <- better$x
        y <- better$y
        
       
        #analysera
        res2 <- lm(gene~res)
        res2
        
        res_gene[j] <- summary(res2)$coefficients[2,4]
        
        #plot(x, y, main="SNP", xlab= "Genotype", ylab="Expression", pch=19, col=pal(length(order))[order] )
        
        
        
        plot(x, y, main=paste(mcols(vcf)[["paramRangeID"]], names(vcf), sep=", ")[i], xlab= "Genotype", ylab="Expression", pch=19, col=pal(length(order))[order] )
        abline(res2)
        #mtext(summary(res2)$coefficients[2,4], summary(paste(res2)$coefficients[2,4], summary(res2)$adj.r.squared, side=3, padj=-1)
       
        #textAttSkrivaTillPlot <- paste(summary(res2)$coefficients[2,4], summary(res2)$adj.r.squared, sep=", ")
        textAttSkrivaTillPlot <- paste("Pvalue:", round(summary(res2)$coefficients[2,4], 3), "R2:", round(summary(res2)$adj.r.squared, 3))
        mtext(textAttSkrivaTillPlot, side=3)
        
        
        #mtext(summary(res2)$adj.r.squared, side=4, padj=1)
        #mtext(summary(res2)$coefficients[2,4],side=3, padj=0)
       

      }
    
        
        dev.off()
        
    lst[[i]] <- res_gene 
    
  }

#dev.off()

