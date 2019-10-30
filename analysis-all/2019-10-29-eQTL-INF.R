#ställ dig i rätt mapp (användarspecifik, så en annan mapp på jespers system)
setwd("/Users/hanbjo/Desktop/2019-08-21-OLINK/repo/jesper")


#where is data stored (relatie path)
datastorage <- "../../data/"
scriptsFunctions <- "functions-all"

#läs in SNP data
#header=TRUE tar bort V:na och gör så att första kolumnen blir heading

#df <- read.table("../../data/191025_CHB_1E04_list.txt", header=TRUE) KORTARE LISTA (50 SNP ca)
df <- read.table("../../data/191025_CHB_5E04_list.txt", header=TRUE)

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
  gr2
  

#Ladda paketet VariantAnnotation
  library("VariantAnnotation")
  
    fl <- "../../data/vcfs/ASMADA_ASAP-MASAP_1KG_chr22.dose_maf-0.001_rsq-0.3-overlapping-variants.vcf.gz"
    #param <- ScanVcfParam(which=gr2)
    #vcf <- readVcf(fl, "hg19", param=gr2)
    
    vcf <- readVcf(fl, "hg19", param=ScanVcfParam(which=gr2))
    vcf
    

  
  