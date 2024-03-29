#############################################################################
## This script was developed by Dr. Erin Sauer for Sauer et al. (submitted) #
# "Male pathology regardless of behaviour drives transmission in an         #
#  avian host-pathogen system"                                              #
# The script fits the SEIR model described in the manuscript                #
#############################################################################

library(deSolve)
library(tidyverse)

#beta = transmission coefficient OR average number of infectious contacts
Bmt <- 0.133928571
Bft <- 0.127513228
Bf <- 0.1375
Bm <- 0.335185185

#sigma = is the rate at which exposed individuals become infectious
sigma.male <- 1/3.66666667
sigma.female <- 1/5

#gamma = rate at which infectious individuals recover
gamma.male <- 1/21
gamma.female <- 1/11

#mu = mortality rate of infectious individuals
#zero in controls
mu.male = -log(1-3/7)/35 
mu.female = -log(1-1/6)/35

######### SEIR model ################

SEIR2 <- function(time, current_state, params){
  
with(as.list(c(current_state, params)),{ 
    N <- S.m+E.m+I.m+R.m+S.f+E.f+I.f+R.f
    dS.m <- -(beta.m*S.m*I.m+
              beta.ft*S.m*I.f)/N
    dE.m <- (beta.m*S.m*I.m+
               beta.ft*S.m*I.f)/N - sigma.m*E.m
    dI.m <- sigma.m*E.m - gamma.m*I.m - mu.m*I.m
    dR.m <- gamma.m*I.m
    dM.m <- mu.m*I.m    
   
    dS.f <- -(beta.f*S.f*I.f+
              beta.mt*S.f*I.m)/N
    dE.f <- (beta.f*S.f*I.f+
               beta.mt*S.f*I.m)/N - sigma.f*E.f
    dI.f <- sigma.f*E.f - gamma.f*I.f - mu.f*I.f
    dR.f <- gamma.f*I.f
    dM.f <- mu.f*I.f

    return(list(c(dS.m, dE.m, dI.m, dR.m, dM.m, dS.f, dE.f, dI.f, dR.f, dM.f)))
  })
}

times <- 0:200
params <- c(beta.m=Bmm, beta.f=Bff, beta.mt=Bme, beta.ft=Bfe,
            gamma.m=gamma.male, gamma.f=gamma.female, 
            sigma.m=sigma.male, sigma.f=sigma.female,
            mu.m=mu.male, mu.f=mu.female)

E.m=0; R.m=0; M.m=0; R.f=0; M.f=0; E.f=0

# varying sex ratios of flocks 
S.m.list <- c(99,74,49,24,0)
S.f.list <- c(0,25,50,75,99)
I.m.list <- c(1,1,1,1,0)
I.f.list <- c(0,0,0,0,1)

plotn <- list()
modeln <- list()
s.modeln <- list()

# run the SEIR2 function and make databases and plots

for(i in 1:5){
S.m = S.m.list[i]
S.f = S.f.list[i]
I.m = I.m.list[i]
I.f = I.f.list[i]
initial_state <- c(S.m=S.m, E.m=E.m, I.m=I.m, R.m=R.m, 
                     M.m=M.m, S.f=S.f, E.f=E.f, I.f=I.f, R.f=R.f, M.f=M.f)
model <- ode(initial_state, times, SEIR2, params)
summary(model)
model <- as.data.frame(model)
model$S <- model$S.m + model$S.f
model$E <- model$E.m + model$E.f
model$I <- model$I.m + model$I.f
model$R <- model$R.m + model$R.f
model$M <- model$M.m + model$M.f
s.modeln[[i]] <- model
model <- model[c(1,12:16)]
model <- model %>% 
  pivot_longer(!time, names_to = "Parameters", values_to = "Individuals")
modeln[[i]] <- model
p <- ggplot(model, aes(x = time, y = Individuals)) + 
  geom_line(aes(color = Parameters)) +
  theme_bw()+
  theme(axis.text=element_text(size=13, color="black"),
        panel.border = element_blank(),
        axis.title=element_text(size=15),
        legend.position="none",
        legend.text= element_blank(),
        legend.title= element_blank())+
  ylab("") + xlab("Days")

plotn[[i]] <- p

}

#get the plots and model outputs
plotn.100 <- plotn[[1]]
plotn.75 <- plotn[[2]]
plotn.50 <- plotn[[3]]
plotn.25 <- plotn[[4]]
plotn.0 <- plotn[[5]]

model.100 <- modeln[[1]]
model.75 <- modeln[[2]]
model.50 <- modeln[[3]]
model.25 <- modeln[[4]]
model.0 <- modeln[[5]]
