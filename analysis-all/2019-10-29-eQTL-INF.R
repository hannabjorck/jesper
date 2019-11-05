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

gr2 <- gr[seqnames(gr)=="22"]
  

#Ladda paketet VariantAnnotation
library("VariantAnnotation")
  
fl <- "../../data/vcfs/ASMADA_ASAP-MASAP_1KG_chr22.dose_maf-0.001_rsq-0.3-overlapping-variants.vcf.gz"
#param <- ScanVcfParam(which=gr2)
#vcf <- readVcf(fl, "hg19", param=gr2)
    
vcf <- readVcf(fl, "hg19", param=ScanVcfParam(which=gr2))
vcf
    

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

#nu är det så att du måste matcha individ i vcf med individ i se, du kan smidigast göra det såhär: GÖR FUNKTION AV DETTA ""match_vcf_object_to_another_vcf_object"
#först GT

#fixa till dubbelnamnen i vcf ID
id_vcf <- colnames(vcf)
id_vcf2 <- sub(".*_","", id_vcf)
#byt ut de gamla id mot dessa nya i id_vcf2
colnames(vcf) <- id_vcf2

#finns alla x i y, och vice versa
id_se <- colnames(se)
tf <- id_vcf2 %in% id_se
tf2 <- id_se %in% id_vcf2

se2 <- se[,tf2]

#matcha på så sätt att individerna är i samma ordning
vcf2 <- vcf[,tf]
m <- match(colnames(se2), colnames(vcf2))


tf <- colnames(vcf2)[m]==colnames(se2)
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

#lm gilla rinte characters, gör därför värdena till numeric innan du gör lm
snp1X <- as.numeric(snp1X)

#Ta ut en gen från den matchade datan
gene1 <- assays(se2)[["counts"]][1,]

#plotta. SNP måste göras om till factor, för att ta bort snuffarna
plot(factor(snp1X), gene1)

#analysera
res1 <- lm(gene1~snp1X)

#För att kolla alla gener kring en snp. HANNA, KOMMENTERA HÄR!!
#giantSNP <- granges(vcf2m[1,])+200000
hits <- findOverlaps(granges(vcf2m[1,])+200000, se2)

# --> Hits object with 2 hits and 0 metadata columns, betyder att i se2 så matchade element 14255 och 14256 (i detta fall två gener som "träffar" giantSNP)
#plocka ut dessa två hits och spara dem i ett nytt objekt (se3)
se3 <- se2[subjectHits(hits)]


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
    hits <- findOverlaps(granges(vcf2m[i,])+200000, se2)
    
    # plocka ut hits och spara dem i ett nytt objekt (se3)
    se3 <- se2[subjectHits(hits)]
    
    
        #en till for loop
        res_gene <- rep(NA, nrow(se3))
        
        for(j in 1:nrow(se3)){
    
        #Ta ut en gen från den matchade datan
        gene <- assays(se3[j,])[["counts"]]
        gene <- as.numeric(gene)
        
        
        #analysera
        res2 <- lm(gene~res)
        res2
        
        res_gene[j] <- summary(res2)$coefficients[2,4]
      }
    
    lst[[i]] <- res_gene 
    
  }



