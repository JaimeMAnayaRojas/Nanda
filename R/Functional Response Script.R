####  This is the script to analyses Nanda's project module data

# https://xcelab.net/rm/statistical-rethinking/

# Install the rethinking package
# https://github.com/rmcelreath/rethinking




library(rethinking)
rm(list=ls(all=TRUE))

data <- read.csv("data/FRdata_2.csv")
data$Eaten <- data$R - data$Left

# remove dead
str(data)
data <- data[-which(is.na(data$Eaten)),]

data$Expose = factor(data$Infection.status)

levels(data$Expose)
levels(data$Expose) <- c("E", "C","E","E") 
data$Expose


# run the analyses in Batch 2

B4 <- subset(data, Batch == 4)



B4$Expose
B4$Exp = factor(B4$Expose)
levels(B4$Expose)
levels(B4$Exp) <- c(1,0) # recode the exposure to  exopse = 1 and non-expose or c = 0


B4$ID = factor(B4$Individual)

levels(B4$ID) <- 1:length(levels(B4$ID))


df <- B4[,c("Eaten", "R", "Exp", "ID")]

df$Exp = as.numeric(as.character(df$Exp))

# get formula for rethinking
#glimmer(Eaten ~ R + Exp + (1|ID), df)


dlist = list(
  N = length(df$Eaten),
  Nid = (length(levels(df$ID))),
  Eaten = df$Eaten,
  ID = as.integer(df$ID),
  Exp = df$Exp,
  R = df$R
  
)


mod2 <- stan("R/FR_mod.stan", data= dlist, iter = 1000,# warmup = 1000, 
             chains =4, cores = 4,
             control = list(adapt_delta = 0.9, max_treedepth = 13))


Females_Sum = precis(mod2, depth = 2, prob=.95)

write.csv(Females_Sum, "results/Females_Sum.csv")
# tracerplot(mod2)

post = extract.samples(mod2)

names(post)

length(post$a)

p_link <- function(post, Exp = 0, R){
  
  FR <- with(post, 
             
             ((10^(a + a_Exp*Exp)) * R) / (1 + (10^(a + a_Exp*Exp)) * (10^(h + h_Exp)) * R)
             
             )
  return(FR)
  
}

median(p_link(post, Exp= 1, 5))

mean(p_link(post, Exp= 1, 25))

xR = seq(from=0.1, to= 25, length= 100)


Exp0 <- sapply(1:length(xR), function(i) p_link(post, Exp = 0, xR[i]))
Exp1 <- sapply(1:length(xR), function(i) p_link(post, Exp = 1, xR[i]))


LOS <- function(x=NULL){
  
  out =   (length(which(x > 0)) / length(x))*100
  
  return(round(out, 3))
  
}


LOS(Exp1-Exp0)# 63.4% that the curve of exposure is higher than the control

apply(Exp1,2, mean)



sdf <- data.frame(R= c(xR,xR), 
                  Expose = c(rep("C",length(xR)),rep("E",length(xR))),
                  Eaten = c(apply(Exp0, 2, median),apply(Exp1, 2, median)),
                  lc = c(t(apply(Exp0, 2, HPDI, prob = .68))[,1], t(apply(Exp1, 2, HPDI, prob = .68))[,1]),
                  uc = c(t(apply(Exp0, 2, HPDI, prob = .68))[,2], t(apply(Exp1, 2, HPDI, prob = .68))[,2]))


C = ggplot(B4, aes(x = R, y = Eaten, group= Expose, fill = Expose)) + 
  geom_point(size = 1.5, shape = 21, position = position_jitterdodge(dodge.width=1.2), alpha = 0.2, show.legend = F) + 
  geom_line(data=sdf, aes(x=R, y=Eaten, colour = Expose), size = 1) + 
  geom_ribbon(data=sdf, aes(ymin = lc, ymax = uc, fill= Expose, group= Expose,), 
              alpha = 0.2, size = 0.01, show.legend = F) +
  scale_fill_manual(values=c("deepskyblue",
                             "black"))+
  scale_color_manual(values=c("deepskyblue",
                              "black")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ylab("Eaten prey \n (mean, 68% CI)") +
  xlab("Prey density") + ylim(c(0,13))


## Results females
LOS(10^(post$a + post$a_Exp) - 10^(post$a)) # 71.7% that expose copepods had a higher attack rate
LOS(10^(post$h + post$h_Exp) - 10^(post$h)) # 86.85% that exposed copepods have a higher handling time... 



#  then handling time of the expose copepods is a bit longer.

### Model for the males

B2 <- subset(data, Batch == 2)
B2$Expose
B2$Exp = factor(B2$Expose)
levels(B2$Expose)
levels(B2$Exp) <- c(1,0) # recode the exposure to  exopse = 1 and non-expose or c = 0


B2$ID = factor(B2$Individual)

levels(B2$ID) <- 1:length(levels(B2$ID))


df <- B2[,c("Eaten", "R", "Exp", "ID")]

df$Exp = as.numeric(as.character(df$Exp))

# get formula for rethinking
#glimmer(Eaten ~ R + Exp + (1|ID), df)


dlist = list(
  N = length(df$Eaten),
  Nid = (length(levels(df$ID))),
  Eaten = df$Eaten,
  ID = as.integer(df$ID),
  Exp = df$Exp,
  R = df$R
  
)


mod2 <- stan("R/FR_mod.stan", data= dlist, iter = 1000,# warmup = 1000, 
             chains =4, cores = 4,
             control = list(adapt_delta = 0.9, max_treedepth = 13))


malesSum = precis(mod2, depth = 2, prob = .95)
write.csv(malesSum, "results/Males_vals.csv")
# tracerplot(mod2)

post = extract.samples(mod2)

names(post)

length(post$a)

p_link <- function(post, Exp = 0, R){
  
  FR <- with(post, 
             
             ((10^(a + a_Exp*Exp)) * R) / (1 + (10^(a + a_Exp*Exp)) * (10^(h + h_Exp)) * R)
             
  )
  return(FR)
  
}

median(p_link(post, Exp= 1, 5))

mean(p_link(post, Exp= 1, 25))

xR = seq(from=0.1, to= 25, length= 100)


Exp0 <- sapply(1:length(xR), function(i) p_link(post, Exp = 0, xR[i]))
Exp1 <- sapply(1:length(xR), function(i) p_link(post, Exp = 1, xR[i]))


LOS <- function(x=NULL){
  
  out =   (length(which(x > 0)) / length(x))*100
  
  return(round(out, 3))
  
}


LOS(Exp1-Exp0)# 63.4% that the curve of exposure is higher than the control

apply(Exp1,2, mean)



sdf <- data.frame(R= c(xR,xR), 
                  Expose = c(rep("C",length(xR)),rep("E",length(xR))),
                  Eaten = c(apply(Exp0, 2, median),apply(Exp1, 2, median)),
                  lc = c(t(apply(Exp0, 2, HPDI, prob = .68))[,1], t(apply(Exp1, 2, HPDI, prob = .68))[,1]),
                  uc = c(t(apply(Exp0, 2, HPDI, prob = .68))[,2], t(apply(Exp1, 2, HPDI, prob = .68))[,2]))


MP = ggplot(B2, aes(x = R, y = Eaten, group= Expose, fill = Expose)) + 
  geom_point(size = 1.5, shape = 21, position = position_jitterdodge(dodge.width=1.2), alpha = 0.2, show.legend = F) + 
  geom_line(data=sdf, aes(x=R, y=Eaten, colour = Expose), size = 1) + 
  geom_ribbon(data=sdf, aes(ymin = lc, ymax = uc, fill= Expose, group= Expose,), 
              alpha = 0.2, size = 0.01, show.legend = F) +
  scale_fill_manual(values=c("deepskyblue",
                             "black"))+
  scale_color_manual(values=c("deepskyblue",
                              "black")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ylab("Eaten prey \n (mean, 68% CI)") +
  xlab("Prey density") + ylim(c(0,13))


### Results females
LOS(10^(post$a + post$a_Exp) - 10^(post$a)) # 66.95% that the expose have higher attack rate.  there are no major difference in the attack rate
LOS(10^(post$h + post$h_Exp) - 10^(post$h)) # 46.95% that the expose have a higher handing time. no differences
###

library("ggpubr")
# theme_set(
#   # theme_bw() +
#   #   theme(legend.position = "top")
# )



figure <- ggarrange(C, MP,
                    labels = c("A)", "B)"),
                    ncol = 2, nrow = 1)
figure


svg("plots/PlotB2 and B4.svg", width = 7, height = 5 )
figure 
graphics.off()