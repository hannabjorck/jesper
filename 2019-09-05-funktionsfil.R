
exempelFunktion <- function(x, y){
  z <- x+y
  q <- z*2
  # sista värdet returneras från funktionen
  q/3
}


#funktioner
subset_match_and_combine_on_id <- function(df1, df2){

  #ploka ut ID  
  id1 <- df1$ID
  id2 <- df2$ID
# vilka idf finns i idp2. Matcha ihop ID, dvs bara plocka ut de ID från df2 som även finns i df2. 
  tf1 <- id2 %in% id1
  tf2 <- id1 %in% id2
# plocka ut ALL ID-data som överlappar.
  dfid2 <- df2[tf1, ]
  dfid1 <- df1[tf2, ] 

# matcha så att ID ligger i samma ordning
  ID1 <- dfid1$ID
  ID2 <- dfid2$ID
  m <- match(ID2,ID1)
 
  tf <- dfid1[m,]$ID==dfid2$ID
  if(!all(tf)) stop("something is wrong")
  dfid1m <-dfid1[m,]
  #returnera den ihopslagna matrisen
  ret <- cbind(dfid1m,dfid2)
  ret
}


#regressionsfunktion
regression_cusp <- function(df, covs=NULL){
  
  #specify which cusps to use in the analysis (only cusp==2 and cusp==3)
  cusp <- as.integer(df$Perioperative_Data__Number_of_cusps)
  tf<- cusp==3
  tf2<- cusp==2
  
  #add the cusp variable, which has been transformed to integers
  df$cusp <- cusp
  
  #subset df to only includef the cusps we are interested in
  df2 <- df[tf|tf2, ]
 
  #subset and grab only the columns with protein measures (OLINK)
  prot <- df2[ ,grep("OLINK", colnames(df2))]

  
  #create the correct formula for the linear regression lm function, based on covs input variable
  if(is.null(covs)){ 
    form <- formula("exp~cusp")
  }else{
    form <- formula(paste("exp~cusp", paste(covs, collapse = "+"), sep="+"))
  }
  
  #create result vectors, to collect the results from the for-loop
  res <- rep(NA, ncol(prot))
  resfc <- rep(NA, ncol(prot))
  
  for (i in 1:ncol(prot)){
    
    #pick out protein expression for i
    exp <- as.numeric(prot[,i])
    
    #add or replace exp in the df2 dataframe (in order to use the i:th protein)
    df2$exp <- exp
    
    #run model generation basedn on 
    model1 <- lm(formula=form, data=df2)
    
    #use the summary function to get a richer output
    summary2 <- summary(model1)
    #pick out the matrix of coefficients from the regression
    mat <- (summary2)$coefficients 
    #only pick out the p-vcalue element
    pmat <- mat['cusp','Pr(>|t|)']
    #store pvalues in vector
    res[i]  <-pmat
    #store fc:s in vector
    resfc[i] <- mean(exp[df2$cusp==3],na.rm=TRUE) - mean(exp[df2$cusp==2],na.rm = TRUE)
  }

  #add pvalues and fold change to the same data frame, which we can return as one object
df_p_fc <- data.frame(res, resfc)
rownames(df_p_fc) <- colnames(prot)

#return object
df_p_fc
    
}

#volcano plot-funktion
 
volcanoplot <- function(df_p_fc, main=""){

plot(df_p_fc[, "resfc"], -log10(df_p_fc[, "res"]), main=main, xlab="fold-change" , ylab="-log10 pvalue")
bonf <- 0.05/length(df_p_fc[, "res"])

abline(h=-log10(bonf))

fdrcutoff <- fdrcutoffFromPvalues(df_p_fc[, "res"])
abline(h=-log10(fdrcutoff), col="red")

}





fdrcutoffFromPvalues <- function(p){
  
  p.sort <- sort(p)
  fdr <- p.adjust(p.sort, "fdr")
  
  #vilka fdr-värden är under 0.05?
  tf <- fdr<0.05
  max <- max(which(tf))
  
  #din fdr cutoff (är alltså pmax)
  pmax <- p.sort[max]
  pmax
}

whichProteinsAreAboveCutoff <- function(regresult, cutoff){
 tf <- regresult[, "res"]<=cutoff
 hej <- regresult[tf, ]
  hej
}


  
#regressionsfunktion
regression_for_volcano_two_groups <- function(df, tf1, tf2, totest, covs=NULL){
  
  #subset df to only include the sets we are interested in
  df2 <- df[tf1|tf2, ]
  
  #subset and grab only the columns with protein measures (OLINK)
  prot <- df2[ ,grep("OLINK", colnames(df2))]
  
  
  #create the correct formula for the linear regression lm function, based on covs input variable
  if(is.null(covs)){ 
    form <- formula(paste("exp~", totest))
  }else{
    form <- formula(paste(paste("exp~",totest, sep=""), paste(covs, collapse = "+"), sep="+"))
  }
  
  #create result vectors, to collect the results from the for-loop
  res <- rep(NA, ncol(prot))
  resfc <- rep(NA, ncol(prot))
  
  for (i in 1:ncol(prot)){
    
    #pick out protein expression for i
    exp <- as.numeric(prot[,i])
    
    #add or replace exp in the df2 dataframe (in order to use the i:th protein)
    df2$exp <- exp
    
    #run model generation basedn on 
    model1 <- lm(formula=form, data=df2)
    
    #use the summary function to get a richer output
    summary2 <- summary(model1)
    #pick out the matrix of coefficients from the regression
    mat <- (summary2)$coefficients 
    #only pick out the p-vcalue element
    pmat <- mat[totest,'Pr(>|t|)']
    #store pvalues in vector
    res[i]  <-pmat
    #store fc:s in vector
    resfc[i] <- mean(exp[tf1],na.rm=TRUE) - mean(exp[tf2],na.rm = TRUE)
  }
  
  #add pvalues and fold change to the same data frame, which we can return as one object
  df_p_fc <- data.frame(res, resfc)
  rownames(df_p_fc) <- colnames(prot)
  
  #return object
  df_p_fc
  
}
