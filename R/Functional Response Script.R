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
glimmer(Eaten ~ R + Exp + (1|ID), df)

# 
mod1 <- ulam(

  alist(
    Eaten ~ dnorm( mu , sigma ),
    mu <- Intercept +
      b_R*R +
      b_Exp*Exp +
      v_Intercept[ID],
    Intercept ~ dnorm(0,10),
    b_R ~ dnorm(0,10),
    b_Exp ~ dnorm(0,10),
    v_Intercept[ID] ~ dnorm(0,sigma_ID),
    sigma_ID ~ dcauchy(0,2),
    sigma ~ dcauchy(0,2)
  ), data =df

)
# 
# summary(mod1)
# 
rethinking::stancode(mod1)

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


precis(mod2)

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


ggplot(B4, aes(x = R, y = Eaten, group= Expose, fill = Expose)) + 
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
  xlab("Prey density")


LOS(10^(post$a + post$a_Exp) - 10^(post$a)) # there are no major difference in the attack rate

LOS(10^(post$h + post$h_Exp) - 10^(post$h)) # but there are stronger differences in the handling time... 


6*46

#  then handling time of the expose copepods is a bit longer.



