####  This is the script to analyses Nanda's project module data
library(rethinking)
library(data.table)


rm(list=ls(all=TRUE))

# Functions
LOS <- function(x=NULL){
  
  out =   (length(which(x > 0)) / length(x))*100
  
  return(round(out, 3))
  
}


###########################################################################################################
# First get the data
getwd() # check your working directory

Data <- (read.csv("data/Data.csv"))


library("plyr")
vars = names(Data)[c(1:5,7,10)]

vars
df = ddply(Data, vars ,summarise, Size= mean(Size), Parasites= mean(Parasites))



# So, something we want to look at is if parasite exposure and infection alter life history traits such somathic growth and survival
# of copepods
# You have look at a set of copepods for four weeks, so we can compare the growth in those time intervals.
# we need to organize the data in a way that we can use, so I will split it into weeks for now.

w1 <- subset(df, Week == 1)
w2 <- subset(df, Week == 2)
w3 <- subset(df, Week == 3)
w4 <- subset(df, Week == 4)

# Something we need to do is to calculate how the size of an individual change from week 1 to 2, from week 2 to 3, and week 3 to 4,
# so, we need to add the size of the same individual from week 2 in front of week 1, and so oo.

# so, I can create a new data for each one

w2a <- w2
w3a <- w3
w4a <- w4

# now we change the "Size" variable to "Size2"
names(w2a)[8]  
names(w2a)[8] <- "Size2"
names(w3a)[8] <- "Size2"
names(w4a)[8] <- "Size2"

# Merge the week sub data sets ---------------------------------------------
names(Data)

w1$name = paste(w1$Plate, w1$Well, sep = "-")
w2a$name = paste(w2a$Plate, w2a$Well, sep = "-")
w1c = merge.data.table(w1, w2a[,c("name", "Size2")], by = "name" )
head(w1c)

w2$name = paste(w2$Plate, w2$Well, sep = "-")
w3a$name = paste(w3a$Plate, w3a$Well, sep = "-")
w2c = merge.data.table(w2, w3a[,c("name", "Size2")], by = "name" )
head(w2c)


w3$name = paste(w3$Plate, w3$Well, sep = "-")
w4a$name = paste(w4a$Plate, w4a$Well, sep = "-")
w3c = merge.data.table(w3, w4a[,c("name", "Size2")], by = "name" )
head(w3c)


## Now put the weeks c into a single dataframe

df_growth = data.frame(rbind(w1c, w2c, w3c))
head(df_growth)

# calculate growth

df_growth$Treatment = factor(df_growth$Treatment)
levels(df_growth$Treatment) # what does it mean "C/E in Treatment?

df_growth$Gz = log(df_growth$Size2/df_growth$Size)
df_growth$Infection = ifelse(df_growth$Parasites > 0, "I+", "I-")
df_growth$Infection[is.na(df_growth$Infection)] = "I-"


library(ggplot2)

ggplot(df_growth, aes(x = Size, y = Gz, colour = Infection)) + geom_point() +
  geom_smooth(method = "lm")


library(brms)

df_growth$z = log(df_growth$Size) - log(200)
df_growth$z2 = log(df_growth$Size^2) - log(300^2)
df_growth$W = df_growth$Week - 3
mod1 = brm(Gz ~ z * Treatment   +  W + (1|name) +  (1|Block), df_growth)
summary(mod1)
conditional_effects(mod1, effects = "z:Treatment")

