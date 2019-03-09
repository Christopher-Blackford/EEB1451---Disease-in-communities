rm(list = ls())

require(deSolve)   
require(rootSolve)
require(FME)
require(ggplot2)
library(plyr)

###################TABLE OF CONTENTS

###[1] Setting up the model
###[2] Running model

################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
### [1] Setting up the model

##### Model description  #####
SISmodel=function(t,y,parameters) 
{ 
  ## Variables
  S=y[1]; I=y[2];
  
  ## Parameters
  p = parameters[1];
  K = parameters[2];
  beta = parameters[3];
  mu = parameters[4];
  alpha = parameters[5];

  
  ## Ordinary differential equations
  logistic_growth <- (1 - (S + I)/K)
  dSdt <- p*(S+I)*logistic_growth - beta*S*I - S * mu
  dIdt <- beta*S*I - I*(alpha + mu)

  return(list(c(dSdt,
                dIdt))); 
}  

###Paramater values (expect p and beta because that's what we're testing)
K <- 5000
mu <- 0.02
alpha <- 0.2

## Initial state
variables0=c(S0=100, I0=1)

## Times at which estimates of the variables are returned
timevec=seq(0,20,1)

################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
### [2] Relationship between p and beta

#Looping to get range of beta values
Beta_values <- seq(0.001, 0.2, by = 0.001)
Suscep_num <- NULL
Infect_num <- NULL

time <- c(999)
S1 <- c(99)
Infect_dyn_null <- data.frame(time, S1)

#constant describing relationship to beta
p_constant <- 10

beta_exp <- 1


for (i in Beta_values){
       beta <- i
       p = p_constant*(beta^beta_exp)
       
       parameters=c(p,
               K,
               beta,
               mu,
               alpha
               )
       
       
       ## Numerical integration through time
       output=lsoda(y = variables0,    # intial values  
                    times = timevec,   # time vector
                    func = SISmodel,   # model
                    parms = parameters # constant parameters
       )
       
       ## Plotting
       colnames(output)=c("time","S1","I1")
       infectionDynamics=as.data.frame(output)
       
       infectionDynamics <- rename(infectionDynamics, c("S1" = paste0("Beta_", i)))
       infectionDynamics <- subset(infectionDynamics, select= -c(I1, time))
       
       
       Infect_dyn_null <- cbind(infectionDynamics, Infect_dyn_null)
     
       assign(paste0("Infectiondynamics", i), infectionDynamics)
       
       plot_title <- paste0("Beta =", beta, ", alpha =", alpha, ", p_constant = ", p_constant)
       plot_title
       
       Linear_plot <- ggplot(infectionDynamics, aes(time)) +
         geom_line(aes(y = S1, colour = "Susceptibles")) + 
         geom_line(aes(y = I1, colour = "Infecteds")) +
         labs(title = plot_title, x = "Time", y = "Population size")
       
       Linear_plot + theme(
         plot.title = element_text(size = 16), 
         axis.text = element_text(size = 16),
         axis.title = element_text(size = 16),
         axis.line = element_line("black"),
         legend.title = element_blank(),
         panel.background = element_blank()
       )
   
       #ggsave(paste0("./plots/SI_direct/Population_dynamics/", plot_title, ".png"), width = 10, height = 6)
       
       }

Infect_dyn_null <- subset(Infect_dyn_null, select= -c(S1, time))

Infect_dyn_null2 <- Infect_dyn_null
Infect_dyn_null <- Infect_dyn_null2

time1_2 <- (Infect_dyn_null[1,] + Infect_dyn_null[2,])/2
time1_3 <- (Infect_dyn_null[1,] + Infect_dyn_null[2,] + Infect_dyn_null[3,])/3
time1_4 <- (Infect_dyn_null[1,] + Infect_dyn_null[2,] + Infect_dyn_null[3,] + Infect_dyn_null[4,])/4


#Time step 1
time1_2 <- t(time1_2)
colnames(time1_2)[1] <- "Average number susceptibles"
time1_2 <- as.data.frame(time1_2)
#time1_2["sorted"]<- c(1:length(time1_2$`Average number susceptibles`))
#time1_2 <- time1_2[with(time1_2, order(-sorted)), ]
#time1_2 <- time1_2["Average number susceptibles"]

write.csv(time1_2, file = "./plots/SI_direct/timestep1_2.csv")

#Time step 2
time1_3 <- t(time1_3)
colnames(time1_3)[1] <- "Average number susceptibles"
time1_3 <- as.data.frame(time1_3)
#time1_3["sorted"]<- c(1:length(time1_2$`Average number susceptibles`))
#time1_3 <- time1_3[with(time1_3, order(-sorted)), ]
#time1_3 <- time1_3["Average number susceptibles"]

write.csv(time1_3, file = "./plots/SI_direct/timestep1_3.csv")

#Time step 3
time1_4 <- t(time1_4)
colnames(time1_4)[1] <- "Average number susceptibles"
time1_4 <- as.data.frame(time1_4)
#time1_4["sorted"]<- c(1:length(time1_2$`Average number susceptibles`))
#time1_4 <- time1_4[with(time1_4, order(-sorted)), ]
#time1_4 <- time1_4["Average number susceptibles"]

write.csv(time1_4, file = "./plots/SI_direct/timestep1_4.csv")


#END
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
