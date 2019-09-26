#ställ dig i rätt mapp
setwd("/Users/hanbjo/Desktop/2019-08-21-OLINK")

#laddar in fenotypdata
load("data/2019-06-16-phex.rdata")

#laddar in proteindata
load("data/raw.RData")

library("readxl")
prot2 <- read_excel("data/ASAP_Olink_clinical_May_2018.xlsx")
prot2 <- as.data.frame(prot2)
prot2 <- prot2[ ,c(1, 94:551)]


#läs in funktionerna i funktionsfilen
source("scripts/2019-09-05-funktionsfil.R")

#addera A på ID för att kunna matcha på fenotypdata 
nyaIDnamn <- paste("A", prot2[ ,"ID"], sep="")

#Ersätt det gamla radnamnet patient. mot det nya radnamnet A
#rownames(prot2) <- nyaradnamn
prot2[ ,"ID"] <- nyaIDnamn

colnames(prot2)[2:459] <- paste("OLINK", colnames(prot2)[2:459], sep="_")

#filter prot2 on proteins with too many NAs (reomve proteins with less than 10 individual measures)
tooManyNA <- apply(prot2, 2, function(x){sum(!is.na(x)) < 10})
prot2 <- prot2[, !tooManyNA]

#använd funktionen från funktionsfilen
df <- subset_match_and_combine_on_id(df1=phex, df2=prot2)

#Fixa interestanta covariates så att rätt format
df$Age <- as.numeric(df$Age)
df$BMI <- as.numeric(df$BMI)
#df$Sex <- as.integer(as.factor(df$Sex))-1
df$'Preop_TEE_Ao-asc-max,_Ao_asc_max_width' <- as.numeric(df$'Preop_TEE_Ao-asc-max,_Ao_asc_max_width')
df$aodia <-df$'Preop_TEE_Ao-asc-max,_Ao_asc_max_width'
df$Preop_Aortic_Stenosis <- as.integer(df$Preop_Aortic_Stenosis)
df$AS <- df$Preop_Aortic_Stenosis
df$Preop_Aortic_Regurgitation <- as.integer(df$Preop_Aortic_Regurgitation)
df$AR <- df$Preop_Aortic_Regurgitation

df_p_fc <- regression_cusp(df=df)
volcanoplot(df_p_fc)

df_p_fc <- regression_cusp(df=df, covs=c("Age","BMI", "Sex", "aodia"))
volcanoplot(df_p_fc)

hopp <- whichProteinsAreAboveCutoff(df_p_fc, fdrcutoffFromPvalues(df_p_fc[, "res"]))


