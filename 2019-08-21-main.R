#ställ dig i rätt mapp
setwd("/Users/hanbjo/Desktop/2019-08-21-OLINK")

#laddar in fenotypdata
load("data/2019-06-16-phex.rdata")

#laddar in proteindata
load("data/raw.RData")

#plotta ålder mot cusp
age <- as.integer(phex$Age)
cusp <- phex$Perioperative_Data__Number_of_cusps
BMI <- phex$BMI
plot(cusp, age)
boxplot(age~cusp)

#färglägga boxplot
boxplot(age~cusp, col="green")
boxplot(age~cusp, col="green", xlab="antal cuspar", ylab="ålder")
boxplot(age~cusp, col="green", xlab="antal cuspar", ylab="ålder", main="ålder vs. cusp")


# gör t-test BAV mot TAV, ett protein

# matcha ID-nummer penotypfil och proteinfil
# Titta vilka ID finns i phex OCH prot
idf <- (phex$ID)
idp <- colnames(raw.clean)[-c(1,2,3,4)]

# gör så att de ser lika ut, dvs byt ut patient. mot A. Det "nya ID prot" heter nu idp2
idp2 <- sub("patient.","A",idp)

# vilka idf finns i idp2. Matcha ihop prot och phex ID, dvs bara plocka ut de ID från phex som även finns i prot. 
tf <- idf %in% idp2
tf2 <- idp2 %in% idf

# titta hur vektorn ser ut
tf

# Titta hur mga som är TRUE. N du får ut är antalet sanna. Antalet falska får man genom att sätta utropstecken framför tf
sum(tf)
sum(!tf)

# välj ut de som matchar, dvs välj ut ett subsett >> dessa läggs under variablen idf2
#idf2 <- idf[tf] 
#idf2

# plocka ut ALL fenotypsdata från phex baserat på ID som överlappar prot och phex.
phex2 <- phex[tf, ]

# matcha proteinfilen och phex2 så att ID ligger i samma ordning
ID <- phex2$ID
m <- match(ID,idp2)

# säkerhetskontroll, matchar ID:s? Alla ska vara TRUE
idp2==ID
idp2[m]==ID
idp2==ID[m]

# slutgiltig matchning. Plocka bort col 1-4 > prot. 
prot <- raw.clean[,-c(1,2,3,4)]

# ta alla rader, men sortera om col baserat på m
wbASAPprot <- prot[ , m]

# transponera col till rader. 
wbASAPprotTransposed <- t(wbASAPprot)
nyaradnamn <- sub("patient.","A", rownames(wbASAPprotTransposed))

# Ersätt det gamla radnamnet patient. mot det nya radnamnet A
rownames(wbASAPprotTransposed) <- nyaradnamn


# ==säger att är det 2 (cusp) så är det SANT, dvs BAV. Det är den som väljs senare. 
cusp <- phex2$Perioperative_Data__Number_of_cusps==2
bav <- phex2[cusp, ]
bavProt <- wbASAPprotTransposed[cusp, ]

# Gör samma för TAV
cuspTAV <- phex2$Perioperative_Data__Number_of_cusps==3
tav <- phex2[cuspTAV, ]
tavProt <- wbASAPprotTransposed[cuspTAV, ]

# gör samma med protein filen (se ovan, rad 3) för BAV och TAV respectively. 

#t-test för ett protein, det som ligger i kol 1, BAV mot TAV. Det är tomt innan komma, 
#vilket innebär att vi anv alla rader. Patienter i det här fallet.
ttt <- t.test(tavProt[,1], bavProt[,1])

#str=structure, för att se strukturen i ttt och se vad du vill plocka ut. 
str(ttt)

#gör t-test för alla prot mha fourloop. Det är där NA-värdet ligger som p-values i det här fallet stoppas in. 
res <- rep(NA, ncol(tavProt))
resfc <- rep(NA, ncol(tavProt))

#själva fourloopen. i kommer att besöka varje värde i listan "(1:ncol(tavProt)" i tur och ordning. plocka ut p-value (första raden)
for (i in 1:ncol(tavProt)){
  res[i] <- t.test(tavProt[,i], bavProt[,i])$p.value
  resfc[i] <- mean(tavProt[,i],na.rm=TRUE) - mean(bavProt[,i],na.rm = TRUE)
}

#res > 454 p-values
#resfc > 454 fold change values

#volcanoplotten

#lägg till bonferronilinjen och färglägg allt under linjen med grått
bonf <- -log10(0.05/ncol(tavProt))
colv=rep("black", ncol(tavProt))
colv[-log10(res) < bonf] <- "grey" 
plot(resfc,-log10(res), col=colv)
abline(h=bonf)


#färglägg allt över linjen med blått el rött 
bonf <- -log10(0.05/ncol(tavProt))
colv <- rep("black", ncol(tavProt))
#colv[-log10(res) < bonf] <- "grey" 
colv[-log10(res) > bonf&resfc>0] <- "blue"
colv[-log10(res) > bonf&resfc<0] <- "red"
plot(resfc,-log10(res), col=colv, pch=1)
abline(h=bonf)

#linear regression för att korregera för kovariater. Simulera två variabler för att förstå koncept
#använd funtionen + (i detta fall c), så +c för att korrigera för variablen c. 
a <- c(1,2,3,4, 5)
b <- c(2,4,6,9,22)
c <- c(2,2,4,8,6)

plot(a, b)
res <- lm(formula=a~b+c)
summary2<- summary(res)

str(summary2)

mat<- (summary2)$coefficients 
pmat <- mat['b','Pr(>|t|)']

#linear regression PÅ VÅR DATA för att korregera för kovariater.
protdata <- wbASAPprotTransposed
phex2

exp <- protdata[,1]
cusp <- phex2$Perioperative_Data__Number_of_cusps
tf<- cusp==3
tf2<- cusp==2
cusp2 <- cusp[tf|tf2]  
age <- phex2$Age
sex <- phex2$Sex
bmi <- as.numeric(phex2$BMI)
bsa <- as.numeric(phex2$BSA_DuBois)
dm <- phex2$Medical_History_Diabetes

df<- data.frame(exp, cusp, age, sex)
model1<- lm(formula=exp~cusp2+age)
plot(cusp2, exp)
abline(model1)

#skapa en volcanoplot av p-values från linjär regression
res <- rep(NA, ncol(tavProt))
resfc <- rep(NA, ncol(tavProt))
  
for (i in 1:ncol(tavProt)){
  exp <- as.numeric(protdata[,i])
  cusp <- as.integer(phex2$Perioperative_Data__Number_of_cusps)
  tf<- cusp==3
  tf2<- cusp==2
  cusp2 <- cusp[tf|tf2]  
  exp2 <- exp[tf|tf2] 
  age2 <- age[tf|tf2]
  bmi2 <- bmi[tf|tf2] 
  bsa2 <- bsa[tf|tf2] 
  dm2 <- dm[tf|tf2] 
  #age <- phex2$Age
  #sex <- phex2$Sex
  #model1<- lm(formula=exp2~cusp2+age2)
  #model2 <-lm(formula=exp2~cusp2+age2+bmi2)
  #model3 <-lm(formula=exp2~cusp2+age2+bsa2)
  model4 <-lm(formula=exp2~cusp2+age2+bsa2+dm2)
  #model5 <-lm(formula=exp2~cusp2+age2+bmi2+dm2)
  
  summary2<- summary(model4)
  mat<- (summary2)$coefficients 
  pmat <- mat['cusp2','Pr(>|t|)']
  res[i]  <-pmat
  
  resfc[i] <- mean(tavProt[,i],na.rm=TRUE) - mean(bavProt[,i],na.rm = TRUE)
}

#bonf <- -log10(0.05/ncol(tavProt))
logpvals <- -log10(res)
fdr5 <- 

colv=rep("black", ncol(tavProt))
colv[-log10(res) < bonf] <- "grey" 
plot(resfc,-log10(res), col=colv)
abline(h=bonf)


#vilket/vilka protein klarar bonferroni?
tf <- -log10(res) > bonf
raw.clean[tf, ]

# sum(tf) --> hur mga klara kriteriet?
# which(tf) --> vilket protein/rad i detta fall






