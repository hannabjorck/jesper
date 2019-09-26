biv(6, 8)

#Gör en matchningsfuntion. Matcha ID från två filer. 
#hitta på data
source("scripts/2019-09-05-funktionsfil.R")
df1 <- data.frame(ID=c(2,4,6,8,9), age=c(55, 77, 23, 90, 1))
df2 <- data.frame(ID=c(9,6,1,8,4), protein=c(66, 78, 13, 190, 11))

df <- subset_match_and_combine_on_id(df1, df2)


