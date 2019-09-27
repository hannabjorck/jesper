#ställ dig i rätt mapp
setwd("/Users/hanbjo/Desktop/2019-08-21-OLINK")

#laddar in fenotypdata
load("data/2019-06-16-phex.rdata")

#laddar in proteindata
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
df$Perioperative_Data__Number_of_cusps <- as.integer(df$Perioperative_Data__Number_of_cusps)
df$cusp <- df$Perioperative_Data__Number_of_cusps


#analys bav mot tav, inga covariates
df_p_fc <- regression_cusp(df=df)
volcanoplot(df_p_fc)

#analys bav mot tav, inkl covariates
df_p_fc <- regression_cusp(df=df, covs=c("Age","BMI", "Sex", "aodia"))
volcanoplot(df_p_fc, main="BAV vs. TAV")

hopp <- whichProteinsAreAboveCutoff(df_p_fc, fdrcutoffFromPvalues(df_p_fc[, "res"]))

#VULCANO olika kombinationer

#kategorisera dilatation i nonDIL och DIL
tf_dil<- df$aodia>=45
tf_dil[is.na(tf_dil)] <- FALSE
tf_nonDIL <- df$aodia<=40
tf_nonDIL[is.na(tf_nonDIL)] <- FALSE

#specify which cusps to use in the analysis (only cusp==2 and cusp==3)
tf_BAV<- df$cusp==2
tf_BAV[is.na(tf_BAV)] <- FALSE

tf_TAV<- df$cusp==3
tf_TAV[is.na(tf_TAV)] <- FALSE

tf_BAV_DIL <- tf_BAV & tf_dil
tf_BAV_nonDIL <- tf_BAV & tf_nonDIL
tf_TAV_DIL <- tf_TAV & tf_dil
tf_TAV_nonDIL <- tf_TAV & tf_nonDIL

pdf("out/volcano.all.groups.pdf")

#utan covs
df_p_fc <- regression_for_volcano_two_groups(df=df, tf1=tf_BAV_DIL, tf2=tf_BAV_nonDIL , totest="aodia")
volcanoplot(df_p_fc, main="BD vs BND unadj")

df_p_fc <- regression_for_volcano_two_groups(df=df, tf1=tf_TAV_DIL, tf2=tf_TAV_nonDIL , totest="aodia")
volcanoplot(df_p_fc, main="TD vs TND unadj")

df_p_fc <- regression_for_volcano_two_groups(df=df, tf1=tf_BAV_nonDIL, tf2=tf_TAV_nonDIL , totest="cusp")
volcanoplot(df_p_fc, main="BND vs TND unadj")

df_p_fc <- regression_for_volcano_two_groups(df=df, tf1=tf_BAV_DIL, tf2=tf_TAV_DIL , totest="cusp")
volcanoplot(df_p_fc, main="BD vs TD unadj")

pdf("out/volcano.all.groups.adjusted_Age.pdf")

#med covs
df_p_fc <- regression_for_volcano_two_groups(df=df, tf1=tf_BAV_DIL, tf2=tf_BAV_nonDIL , totest="aodia", covs=c("Age"))
volcanoplot(df_p_fc, main="BD vs BND adj")

df_p_fc <- regression_for_volcano_two_groups(df=df, tf1=tf_TAV_DIL, tf2=tf_TAV_nonDIL , totest="aodia", covs=c("Age"))
volcanoplot(df_p_fc, main="TD vs TND adj")

df_p_fc <- regression_for_volcano_two_groups(df=df, tf1=tf_BAV_nonDIL, tf2=tf_TAV_nonDIL , totest="cusp", covs=c("Age", "AS"))
volcanoplot(df_p_fc, main="BND vs TND adj")

df_p_fc <- regression_for_volcano_two_groups(df=df, tf1=tf_BAV_DIL, tf2=tf_TAV_DIL , totest="cusp", covs=c("Age"))
volcanoplot(df_p_fc, main="BD vs TD adj")

dev.off()

hopp <- whichProteinsAreAboveCutoff(df_p_fc, fdrcutoffFromPvalues(df_p_fc[, "res"]))


#PEARSON CORRELATION PROT VS AO.DIA

#subset and grab only the columns with protein measures (OLINK)
prot <- df[ ,grep("OLINK", colnames(df))]

#correlationen med forloop
#create result vectors, to collect the results from the for-loop
res_P <- rep(NA, ncol(prot))
resR_P <- rep(NA, ncol(prot))

res_S <- rep(NA, ncol(prot))
resR_S <- rep(NA, ncol(prot))

for (i in 1:ncol(prot)){
  
  #pick out protein expression for i
  protexp <- as.numeric(prot[,i])
  
pearson <- cor.test(protexp, df$aodia, method = c("pearson"))
res_P[i] <- pearson$p.value
resR_P[i] <- pearson$estimate

spearman <- cor.test(protexp, df$aodia, method = c("spearman"))
res_S[i] <- spearman$p.value
resR_S[i] <- spearman$estimate

}


plot(resR_P, resR_S, col="purple")

