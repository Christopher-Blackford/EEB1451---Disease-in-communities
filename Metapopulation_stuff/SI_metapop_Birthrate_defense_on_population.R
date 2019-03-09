rm(list = ls())

require(deSolve)   
require(rootSolve)
require(FME)
require(ggplot2)
setwd("C:/Users/fortinlab/Dropbox/Courses/2016-2017/EEB1451 - Disease in Communities/Term_paper/R")
###################TABLE OF CONTENTS
getwd()
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
  S1=y[1]; I1=y[2]; S2=y[3]; I2=y[4];
  
  ## Parameters
  p = parameters[1];
  K = parameters[2];
  beta = parameters[3];
  mu = parameters[4];
  alpha = parameters[5];
  m12 = parameters[6];
  m21 = parameters[7];
  i12 = parameters[8];
  i21 = parameters[9]
  
  ## Ordinary differential equations
  logistic_growth1 <- (1 - (S1 + I1)/K)
  dS1dt <- p*(S1+I1)*logistic_growth1 - beta*S1*I1 - S1*mu - S1*m12 + S2*m21 #migration

  dI1dt <- beta*S1*I1 - I1*(alpha + mu) - I1*i12 + I2*i21 #migration
  
  logistic_growth2 <- (1 - (S2 + I2)/K)
  dS2dt <- p*(S2+I2)*logistic_growth2 - beta*S2*I2 - S2*mu - S2*m21 + S1*m12 #migration
  
  dI2dt <- beta*S2*I2 - I2*(alpha + mu) - I2*i21 + I1*i12 #migration
  
  
  return(list(c(dS1dt,
                dI1dt,
                dS2dt,
                dI2dt))); 
}  

## Initial state
variables0=c(S1=100, I1=1, S2=0, I2=0) #These values don't link up to model, they just represent first 2 equations

## Times at which estimates of the variables are returned
timevec=seq(0,30,1)


################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
### [1] Describing parameters

###Paramater values (expect p and beta because that's what we're testing)
K <- 5000
mu <- 0.02
alpha <- 0.2
m12 <- 0.01
m21 <- m12
i12 <- 0.005
i21 <- i12


#Looping to get range of beta values
Beta_values <- seq(0.005, 0.1, by = 0.005)

#constant describing relationship to beta
p_constant <- 10
beta_exp <- 1

for (i in Beta_values){
  beta <- i
  p <- p_constant*(beta^beta_exp)
  m12 <- beta
  m21 <- m12
  
  parameters=c(p,
             K,
             beta,
             mu,
             alpha,
             m12,
             m21,
             i12,
             i21
             )
  
  output=lsoda(y = variables0,    # intial values  
             times = timevec,   # time vector
             func = SISmodel,   # model
             parms = parameters # constant parameters
             )
  
  ## Plotting
  colnames(output)=c("time","S1","I1", "S2", "I2")
  infectionDynamics=as.data.frame(output)
  #par(mfrow=c(1,2),mar=c(5,5,2,2))
  
  
  
  plot_title <- paste0("Population size, alpha ", alpha, ", p = ", p_constant, "_", i, "_", beta_exp)
  plot_title
  
  Linear_plot <- ggplot(infectionDynamics, aes(time)) +
    geom_line(aes(y = S1, colour = "Susceptibles1")) + 
    geom_line(aes(y = I1, colour = "Infecteds1")) +
    geom_line(aes(y = S2, colour = "Susceptibles2")) +
    geom_line(aes(y = I2, colour = "Infecteds2")) +
    labs(title = plot_title, x = "Time", y = "Population size")
  Linear_plot
  
  Linear_plot + theme(
    plot.title = element_text(size = 16), 
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16),
    axis.line = element_line("black"),
    legend.title = element_blank(),
    panel.background = element_blank()
  )
 
  #Write out plots
  ggsave(paste0("./plots/SI_2metapop/Population_dynamics/2patch_alpha", alpha, "p =", p_constant, "_", i, "_", beta_exp, ".png"))
  
}