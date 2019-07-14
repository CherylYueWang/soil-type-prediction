# load data
cover_data <- read.csv('cleaned_data.csv', header = TRUE)
cover_data = cover_data[-c(1,14)]

# look at the correlation matrix
#install.packages('corrplot')
#library(corrplot)
#corrplot(cover_data, method="circle",use = "na.or.complete")
#M = cor(cover_data)
#corrplot(M, method="number")

# look at the conditional independence of response (cover-type) and elevation given soil-type

#cover_data$Elevation_cat = cut(cover_data$Elevation,breaks=10, labels=c(1:10))
cover_data$Elevation_cat = cut(cover_data$Elevation,breaks=2)

cover_data$soil_cat = cut(cover_data$Soil_type, breaks =4)
cover_data$cover_cat = cut(cover_data$Cover_Type,breaks=2)
table1 = table(cover_data$cover_cat,cover_data$Elevation_cat, cover_data$soil_cat)
names(dimnames(table1)) <- c("cover_type", "elevation_level","soil_type")

# chi-squre test for independence
E=array(NA,dim(table1))
for (i in 1:dim(table1)[1]) {
  for (j in 1:dim(table1)[2]) {
    for (k in 1:dim(table1)[3]) {
      E[i,j,k]=(sum(table1[i,,k])*sum(table1[,j,k]))/margin.table(table1,3)[k]
    }}}
names(dimnames(E)) <- c("cover_type", "elevation_level","soil_type")

#X2=sum((temp-E)^2/E)
X2 = 0
for (i in 1:dim(table1)[1]) {
  for (j in 1:dim(table1)[2]) {
    for (k in 1:dim(table1)[3]) {
      X2=X2+(table1[i,j,k]-E[i,j,k])^2/E[i,j,k]
    }}}
df = 4*(2-1)*(2-1)
p_value=pchisq(X2,df,lower.tail = FALSE)

table1[,,1]
chisq.test(table1[,,1], correct=FALSE)
table1[,,2]
chisq.test(table1[,,2], correct=FALSE)
table1[,,3]
chisq.test(table1[,,3], correct=FALSE)
table1[,,4]
chisq.test(table1[,,4], correct=FALSE)
X2=sum(chisq.test(table1[,,1], correct=FALSE)$statistic+chisq.test(table1[,,2], correct=FALSE)$statistic+
         chisq.test(table1[,,3], correct=FALSE)$statistic+chisq.test(table1[,,4], correct=FALSE)$statistic)
pchisq(X2,df,lower.tail = FALSE)


# logit model for testing same conditional independence
cover_data$Elevation_cat = cut(cover_data$Elevation,breaks=2)
cover_data$soil_cat = cut(cover_data$Soil_type, breaks =4)
cover_data$cover_cat = cut(cover_data$Cover_Type,breaks=2)
table1 = table(cover_data$cover_cat,cover_data$Elevation_cat, cover_data$soil_cat)

cover=array(NA,c(8,2))
for (i in 1:dim(table1)[1]) {
  for (k in 1:dim(table1)[3]) {
    for (j in 1:dim(table1)[2]) {
      cover[j+2*(k-1),i]=table1[,,k][i,j]
    }}}

soiltype=c('1-10','1-10',
           '11-20','11-20',
           '21-30','21-30',
           '31-40','31-40')
elevation=c('1860-2860','1861-3850',
            '1860-2860','1861-3850',
            '1860-2860','1861-3850',
            '1860-2860','1861-3850')
cover1=cover[,1]
cover2=cover[,2]


temp1=data.frame(soiltype=soiltype,elevation=elevation,cover_1to4=cover1,cover_5to7=cover2)
cover.fit = glm(cbind(cover_1to4,cover_5to7)~ elevation + soiltype,family=binomial(),data=temp1)
summary(cover.fit)
cover.fit2 = update(cover.fit,formula = ~.-elevation) # fit reduce model without elevation 
anova(cover.fit2,cover.fit,test="Chisq") #  compute difference in deviances (LR ratio statistic)


# test for aspect and elevation



# test for cover_type and hillshade_9am independence
cover_data$Shade9_cat = cut(cover_data$Hillshade_9am,c(0,127,255))
cover_data$cover_cat = cut(cover_data$Cover_Type,breaks=2)

data_two = table(cover_data$Shade9_cat,cover_data$Cover_Type)
names(dimnames(data_two)) <- c("hillshade_9am", "cover_type")
result<-chisq.test(data_two)
result


#A function for computing Pearson correlation for IxJ tables & Mantel-Haenszel, M2
#pearson correlation for IxJ tables 
#table = IxJ table or a matrix 
#rscore=vector of row scores 
#cscore=vector of column scores 
pears.cor=function(table, rscore, cscore)
{ 
  dim=dim(table) 
  rbar=sum(margin.table(table,1)*rscore)/sum(table) 
  rdif=rscore-rbar 
  cbar=sum(margin.table(table,2)*cscore)/sum(table) 
  cdif=cscore-cbar 
  ssr=sum(margin.table(table,1)*(rdif^2)) 
  ssc=sum(margin.table(table,2)*(cdif^2)) 
  ssrc=sum(t(table*rdif)*cdif) 
  pcor=ssrc/(sqrt(ssr*ssc)) 
  pcor 
  M2=(sum(table)-1)*pcor^2
  M2
  result=c(pcor, M2)
  result
} 

### For the Pearson correlation coefficent 
### and Mantel-Haenszel, 
### for IxJ tables, you can also use 
### pears.cor() function. 
### Mak sure you run this function first!
### c(1,2) and c(1,2,3,4), are the vectors of score values 
test = pears.cor(data_two, c(1,2),c(1,2,3,4,5,6,7)) 
pvalue = pchisq(test[2],df=6,lower.tail = FALSE)


# can we remove some x variables?
# catogorize all continuous data
cover_data$Elevation_cat = cut(cover_data$Elevation,breaks=2)
cover_data$Aspect_cat = cut(cover_data$Aspect,breaks=2) # 360 degree
cover_data$Slope_cat = cut(cover_data$Aspect,breaks=2)
cover_data$HorToHydr_cat = cut(cover_data$Horizontal_Distance_To_Hydrology,breaks=2)
cover_data$VertToHydr_cat = cut(cover_data$Vertical_Distance_To_Hydrology,breaks=2)
cover_data$HorToRoad_cat = cut(cover_data$Horizontal_Distance_To_Roadways,breaks=2)
#cover_data$HorToHydr_cat = cut(cover_data$Horizontal_Distance_To_Hydrology,breaks=5)
cover_data$Shade9_cat = cut(cover_data$Hillshade_9am,breaks=2)
cover_data$ShadeNoon_cat = cut(cover_data$Hillshade_Noon,breaks=2)
cover_data$Shade3_cat = cut(cover_data$Hillshade_3pm,breaks=2)
cover_data$HorToFire_cat = cut(cover_data$Horizontal_Distance_To_Fire_Points,breaks=2)
cover_data$EucToHydr_cat = cut(cover_data$Euc_Dis_To_Hydrology,breaks=2)
cover_data$Wilderness = cut(cover_data$Wilderness,breaks=2)


# response
cover_data$cover_cat = cut(cover_data$Cover_Type,breaks=2)

#Aspect & Elevation
cover_data$Elevation_cat = cut(cover_data$Elevation,breaks=10)
cover_data$Aspect_cat = cut(cover_data$Aspect,breaks=4) # 360 degree
data_3 = table(cover_data$Elevation_cat,cover_data$Aspect_cat)
names(dimnames(data_3)) <- c("Elevation", "Aspect")
result<-chisq.test(data_3)
result
