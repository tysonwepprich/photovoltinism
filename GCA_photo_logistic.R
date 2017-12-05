# runs logistic regression on unsummarized data (zeros and ones)

# setwd('C:\\R\\Aphalara diapause')

data=read.csv('data/combined 13_15.csv')
str(data)
#data=subset(data1,Year==2014)

# what about latitude covariate?
# (could also depend on first frost, gdd remaining, etc.)
geo <- data.frame(Site = unique(data$Site), Lat = c(39.435, 43.39,
                                                    47.32, 48.75))
data <- merge(data, geo)

mod1 <- glm(diapause~Treatment * Lat,data=data,family=binomial)
mod2 <- glm(diapause~Treatment + Site,data=data,family=binomial)

# what about just North and Southern populations?

agdata=aggregate(data$diapause,by=list(data$Site,data$Treatment),FUN=mean)


#South
output=glm(diapause~Treatment, subset=Site %in% c('Palermo', 'Sutherlin'),data=data,family=binomial)
summary(output)
anova(output,test="Chisq")
coefs =coef(output)

d1=subset(agdata,Group.1 %in% c('Palermo', 'Sutherlin'))

plot(d1$x~d1$Group.2, ylab="Proportion reproductive",xlab="Photoperiod (hrs)", xlim=c(13,18),col="red", cex=1.5, lwd=3,pch=0)
curve(exp(coefs[1]+coefs[2]*x)/(1+exp(coefs[1]+coefs[2]*x)),col="red",add=TRUE)

#North
output=glm(diapause~Treatment, subset=Site %in% c('Bellingham', 'Ephrata'),data=data,family=binomial)
summary(output)
anova(output,test="Chisq")
coefs =coef(output)

d1=subset(agdata,Group.1 %in% c('Bellingham', 'Ephrata'))


curve(exp(coefs[1]+coefs[2]*x)/(1+exp(coefs[1]+coefs[2]*x)),col="blue", add=TRUE)
points(d1$x~d1$Group.2, ylab="Proportion reproductive",xlab="Photoperiod (hrs)", col="blue", cex=1.5, lwd=3, pch=1, add=TRUE)



# predicting diapause based on raster of photoperiods

source('CDL_funcs.R')
library(raster)
doy <- 170
GDD <- brick("meanGDD_2007to2013/meanGDD_07_13.grd")
template <- crop(GDD[[1]], extent(-124, -122.5, 44, 45))
template[!is.na(template)] <- 0
photo <- RasterPhoto(template, doy, 100)

prop_diap <- 1 - exp(coefs[1]+coefs[2]*photo)/(1+exp(coefs[1]+coefs[2]*photo))
plot(prop_diap)

#Palermo
output=glm(diapause~Treatment, subset=Site=='Palermo',data=data,family=binomial)
summary(output)
anova(output,test="Chisq")
coefs =coef(output)

d1=subset(agdata,Group.1=="Palermo")

plot(d1$x~d1$Group.2, ylab="Proportion reproductive",xlab="Photoperiod (hrs)", xlim=c(13,18),col="red", cex=1.5, lwd=3,pch=0)
curve(exp(coefs[1]+coefs[2]*x)/(1+exp(coefs[1]+coefs[2]*x)),col="red",add=TRUE)

#Bellingham

output=glm(diapause~Treatment, subset=Site=='Bellingham',data=data,family=binomial)
summary(output)
anova(output,test="Chisq")
coefs =coef(output)

d1=subset(agdata,Group.1=="Bellingham")

curve(exp(coefs[1]+coefs[2]*x)/(1+exp(coefs[1]+coefs[2]*x)),col="blue", add=TRUE)
points(d1$x~d1$Group.2, ylab="Proportion reproductive",xlab="Photoperiod (hrs)", col="blue", cex=1.5, lwd=3, pch=1, add=TRUE)


#Ephrata

output=glm(diapause~Treatment, subset=Site=='Ephrata',data=data,family=binomial)
summary(output)
anova(output,test="Chisq")
coefs =coef(output)

d1=subset(agdata,Group.1=="Ephrata")

points(d1$x~d1$Group.2, ylab="Proportion reproductive",xlab="Photoperiod (hrs)", col="green", cex=1.5, lwd=3, pch=2,add=TRUE)
curve(exp(coefs[1]+coefs[2]*x)/(1+exp(coefs[1]+coefs[2]*x)), col="green",add=TRUE)


#Sutherlin

output=glm(diapause~Treatment, subset=Site=='Sutherlin',data=data,family=binomial)
summary(output)
anova(output,test="Chisq")
coefs =coef(output)

d1=subset(agdata,Group.1=="Sutherlin")

points(d1$x~d1$Group.2, ylab="Proportion reproductive",xlab="Photoperiod (hrs)", col="orange", cex=1.5, lwd=3,pch=3,add=TRUE)
curve(exp(coefs[1]+coefs[2]*x)/(1+exp(coefs[1]+coefs[2]*x)),col="orange", add=TRUE)

