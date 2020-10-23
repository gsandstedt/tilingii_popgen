###install and load packages

install.packages("mice")
install.packages("VIM")
install.packages("caret")
install.packages("klaR")
install.packages("rattle.data")


library(mice)
library(VIM)
library(MASS)
library(ggplot2)
library(tidyverse)
library(caret)
library(e1071)
library(klaR)
library(ellipse)

# set directory and load data
setwd("")

df.mis <- read.csv("morphological_data.csv", header = TRUE)

#visualize variables with missing data

mice_plot <- aggr(df.mis, col=c('navyblue','yellow'),
                  numbers=TRUE, sortVars=TRUE,
                  labels=names(df.mis), cex.axis=.7,
                  gap=3, ylab=c("Missing data","Pattern"))

#impute missing data 50 reiterations

imputed_Data <- mice(df.mis, m=1, maxit = 50, method = 'pmm', seed = 500)
summary(imputed_Data)

#new data set with imputed numbers
com <- mice::complete(imputed_Data,1)

model <- lda(population~., data = com)
plot(model)


#Estimate Parameters for standardizing data
preproc.param <- com %>% 
  preProcess(method = c("center", "scale"))

# standardize the data using the estimated parameters
com.transformed <- preproc.param %>% predict(com)

#check
mean(com.transformed$leaf_width)


# Fit the model
model <- lda(population~., data = com.transformed)
model

# model predictions
predict(model)

# compute LDA

lda.data <- cbind(com.transformed, predict(model)$x)

# graph with ggplot

p <- ggplot(lda.data, aes(x=LD1, y=LD2, col=population) ) + geom_point( size = 3, aes(color = population))
q <- p + scale_color_manual(values=c("darkorange1", "darkorchid3", "dodgerblue"))
q + theme(axis.text.x = element_text(color = "black", 
                                     size = 15, angle = 0),
          axis.text.y = element_text(color = "black", 
                                     size = 15, angle = 0)) + scale_x_discrete(breaks=c(-6, -3, 0, 3, 6), limits=c(-6, -3, 0, 3, 6)) + expand_limits(x=c(-6,6)) + scale_y_discrete(breaks=c(-8, -4, 0, 4), limits=c(-8, -4, 0, 4)) + expand_limits(y=c(-8,4)) 


# LDA with 95% ellipses

df<-lda(population ~ .,  data = com.transformed) 

datPred<-data.frame(population=predict(df)$class,predict(df)$x) 

dat_ell <- data.frame() 
for(g in levels(datPred$population)){ 
  dat_ell <- rbind(dat_ell, cbind(as.data.frame(with(datPred[datPred$population==g,], ellipse(cor(LD1, LD2), 
                                                                                              scale=c(sd(LD1),sd(LD2)), 
                                                                                              centre=c(mean(LD1),mean(LD2))))),population=g)) 
} 


p <- ggplot(datPred, aes(x=LD1, y=LD2, col=population) ) + geom_point( size = 3, aes(color = population))+theme_bw()+geom_path(data=dat_ell,aes(x=x,y=y),size=1,linetype=2) 
q <- p + scale_color_manual(values=c("darkorange1", "darkorchid3", "dodgerblue"))
q + theme(axis.text.x = element_text(color = "black", 
                                     size = 15, angle = 0),
          axis.text.y = element_text(color = "black", 
                                     size = 15, angle = 0)) + scale_x_discrete(breaks=c(-6, -3, 0, 3, 6), limits=c(-6, -3, 0, 3, 6)) + expand_limits(x=c(-6,6)) + scale_y_discrete(breaks=c(-8, -4, 0, 4), limits=c(-8, -4, 0, 4)) + expand_limits(y=c(-8,4)) 


# discriminant factors

df.lda <- lda(population~ ., data= com.transformed)
df.lda.values <- predict(df.lda)
head(df.lda)

