#ställ dig i rätt mapp
setwd("/Users/hanbjo/Desktop/2019-08-21-OLINK")

#laddar in fenotypdata
load("data/2019-06-16-phex.rdata")

#laddar in proteindata
load("data/raw.RData")

library("readxl")
prot2 <- read_excel("data/ASAP_Olink_clinical_May_2018.xlsx")
prot2 <- prot2[ ,c(1, 94:551)]


#läs in funktionerna i funktionsfilen
source("scripts/2019-09-05-funktionsfil.R")

#Plocka bort col 1-4 ur protfilen 
#prot <- raw.clean[,-c(1,2,3,4)]

#transponera col till rader. 
protTransposed <- t(prot)
nyaradnamn <- sub("patient.","A", rownames(protTransposed))

#Ersätt det gamla radnamnet patient. mot det nya radnamnet A
rownames(protTransposed) <- nyaradnamn

#gör protTransposed till en data fram
protTransposed <- as.data.frame(protTransposed)

colnames(protTransposed) <- paste("OLINK", colnames(protTransposed), sep="_")

#tryck in de nya radnamnen i en NY kolumn som du kallar ID. (skapas så > protTransposed$ID)
protTransposed$ID <- nyaradnamn

#använd funktionen från funktionsfilen

df <- subset_match_and_combine_on_id(df1=phex, df2=protTransposed)

#Fixa interestanta covariates så att rätt format
df$Age <- as.numeric(df$Age)
df$BMI <- as.numeric(df$BMI)

df_p_fc <- regression_cusp(df=df)
volcanoplot(df_p_fc)

df_p_fc <- regression_cusp(df=df, covs=c("Age","BMI"))
volcanoplot(df_p_fc)


