rm(list = ls())

require(deSolve)   
require(rootSolve)
require(FME)
require(ggplot2)
setwd("C:/Users/fortinlab/Dropbox/Courses/2016-2017/EEB1451 - Disease in Communities/Term_paper")
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
variables0=c(S1=1000, I1=1, S2=0, I2=0) #These values don't link up to model, they just represent first 2 equations

## Times at which estimates of the variables are returned
timevec=seq(0,20,1)


################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
### [1] Describing parameters

###Paramater values (expect p and beta because that's what we're testing)
p <- 2
K <- 5000
mu <- 0.05
alpha <- 0.1
m12 <- 0.01
m21 <- m12
i12 <- 0.005
i21 <- i12


#Looping to get range of beta values
Beta_values <- seq(0.005, 0.2, by = 0.005)
Suscep_num1 <- NULL
Infect_num1 <- NULL
Suscep_num2 <- NULL
Infect_num2 <- NULL


#constant describing relationship to beta
p_constant <- 1
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
  
  steadyState=runsteady(y = variables0, 
                        fun = SISmodel, 
                        parms = parameters)
  
  Suscep1 <- steadyState$y[1]
  Suscep_num1 <- append(Suscep_num1, Suscep1)
  
  Infect1 <- steadyState$y[2]
  Infect_num1 <- append(Infect_num1, Infect1)
  
  Suscep2 <- steadyState$y[3]
  Suscep_num2 <- append(Suscep_num2, Suscep2)
  
  Infect2 <- steadyState$y[4]
  Infect_num2 <- append(Infect_num2, Infect2)
  
  ## Numerical integration through time
  output=lsoda(y = variables0,    # intial values  
               times = timevec,   # time vector
               func = SISmodel,   # model
               parms = parameters # constant parameters
  )
  
  ## Plotting
  colnames(output)=c("time","S1","I1", "S2", "I2")
  infectionDynamics_meta2=as.data.frame(output)
  
  assign(paste0("Infection2meta_beta", i), infectionDynamics_meta2)
  
  plot_title <- paste0("Population size, alpha ", alpha, ", p = ", p_constant, "_", i, "_", beta_exp)
  plot_title
  
  Linear_plot <- ggplot(infectionDynamics_meta2, aes(time)) +
    geom_line(aes(y = S1, colour = "Susceptibles1")) + 
    geom_line(aes(y = I1, colour = "Infecteds1")) +
    geom_line(aes(y = S2, colour = "Susceptibles2")) +
    geom_line(aes(y = I2, colour = "Infecteds2")) +
    labs(title = plot_title, x = "Time", y = "Population size")
  
  Linear_plot + theme(
    plot.title = element_text(size = 16), 
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16),
    axis.line = element_line("black"),
    legend.title = element_blank(),
    panel.background = element_blank()
  )

  ggsave(paste0("./plots/SI_Birthrate_defense_on_population/mig_trade/", plot_title, ".png"))
  
}



#Getting data into dataframe for plotting - figure out how to do in one line
tempdf <- cbind(Beta_values, Suscep1, Infect1, Suscep2, Infect2)
Metapop_equil2 <- as.data.frame(tempdf, row.names = T)


View(Metapop_equil)






#Susceptible hosts1
plot(infectionDynamics$S1~infectionDynamics$time,
     main="Susceptibles1",ylab="Density",xlab="Time (days)",
     type="l",las=1,
     xlim=c(0,length(timevec)))

#Infected hosts1
plot(infectionDynamics$I1~infectionDynamics$time,
     main="Infecteds1",ylab="Density",xlab="Time (days)",
     type="l",las=1,
     xlim=c(0,length(timevec)))


#Susceptible hosts2
plot(infectionDynamics$S2~infectionDynamics$time,
     main="Susceptibles2",ylab="Density",xlab="Time (days)",
     type="l",las=1,
     xlim=c(0,length(timevec)))

#Infected hosts2
plot(infectionDynamics$I2~infectionDynamics$time,
     main="Infecteds2",ylab="Density",xlab="Time (days)",
     type="l",las=1,
     xlim=c(0,length(timevec)))




################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
### [2] Relationship between p and beta

#Looping to get range of beta values
Beta_values <- seq(0.001, 0.2, by = 0.002)
Suscep_num <- NULL
Infect_num <- NULL
Environ_num <- NULL

#constant describing relationship to beta
p_constant <- 1

beta_exp <- 1


for (i in Beta_values){
       beta <- i
       p = p_constant*(beta^beta_exp)
       
       parameters=c(p,
               K,
               beta,
               mu,
               alpha,
               lambda,
               tau)
       
       steadyState=runsteady(y = variables0, 
                      fun = SISmodel, 
                      parms = parameters)
       
       Suscep <- steadyState$y[1]
       Suscep_num <- append(Suscep_num, Suscep)
       
       Infect <- steadyState$y[2]
       Infect_num <- append(Infect_num, Infect)
       
       Environ <- steadyState$y[3]
       Environ_num <- append(Environ_num, Environ)
       
       }
  
#Getting data into dataframe for plotting - figure out how to do in one line
tempdf <- cbind(Beta_values, Suscep_num, Infect_num, Environ_num)
Linear_data <- as.data.frame(tempdf, row.names = T)




######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
