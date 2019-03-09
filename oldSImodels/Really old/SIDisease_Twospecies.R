rm(list = ls())

library(deSolve)   
library(rootSolve)
library(FME)


##### Model description  #####
SISmodel=function(t,y,parameters) 
{ 
  ## Variables
  S1=y[1]; I1=y[2]; S2=y[3]; I2=y[4]; E=y[5];
  
  ## Parameters
  birth_rate1 = parameters[1];
  birth_rate2 = parameters[2];
  K = parameters[3];
  beta1 = parameters[4];
  beta2 = parameters[5];
  natural_mortality = parameters[6];
  disease_morality = parameters[7];
  lambda = parameters[8];
  mu = parameters[9] 
    
  
  ## Ordinary differential equations
  dS1dt <- birth_rate1 * (S1+I1) * (K - (S1 + I1)/K) - beta1*S1*E - S1 * natural_mortality
  dI1dt <- beta1*S1*E - I1 * disease_morality - I1 * natural_mortality
  dS2dt <- birth_rate2 * (S2+I2) * (K - (S2 + I2)/K) - beta2*S2*E - S2 * natural_mortality
  dI2dt <- beta2*S2*E - I2 * disease_morality - I2 * natural_mortality
  dEdt <- lambda * (I1 + I2) - mu * E
  
  return(list(c(dS1dt,
                dI1dt,
                dS2dt,
                dI2dt,
                dEdt))); 
} 


#####  Intial condition  #####
## Parameters
parameters=c(birth_rate1 <- 1.05,
             birth_rate2 <- 1.01,
             K <- 5000,
             beta1 <- 0.01,
             beta2 <- 0.005,
             natural_mortality <- 0.001,
             disease_morality <- 0.05,
             lambda <- 0.001,
             mu <- 0.2)



## Initial state
variables0=c(S1_0=500,I1_0=1,S2_0=500,I2_0=1,E0=10)

## Times at which estimates of the variables are returned
timevec=seq(0,20,1)
###############################

#####  Running the model   #####
## Numerical integration through time
output=lsoda(y = variables0,    # intial values  
             times = timevec,   # time vector
             func = SISmodel,   # model
             parms = parameters # constant parameters
)

## Plotting
colnames(output)=c("time","S1","I1", "S2","I2", "E")
infectionDynamics=as.data.frame(output)
par(mfrow=c(1,2),mar=c(5,5,2,2))

#Susceptible species1
plot(infectionDynamics$S1~infectionDynamics$time,
     main="Susceptibles1",ylab="Density",xlab="Time (days)",
     type="l",las=1,
     xlim=c(0,length(timevec)))

#Infected species1
plot(infectionDynamics$I1~infectionDynamics$time,
     main="Infecteds1",ylab="Density",xlab="Time (days)",
     type="l",las=1,
     xlim=c(0,length(timevec)))

#Susceptible species2
plot(infectionDynamics$S2~infectionDynamics$time,
     main="Susceptibles2",ylab="Density",xlab="Time (days)",
     type="l",las=1,
     xlim=c(0,length(timevec)))

#Infected species2
plot(infectionDynamics$I2~infectionDynamics$time,
     main="Infecteds2",ylab="Density",xlab="Time (days)",
     type="l",las=1,
     xlim=c(0,length(timevec)))


#Environmental parasite
plot(infectionDynamics$E~infectionDynamics$time,
     main="Env Parasite",ylab="Density",xlab="Time (days)",
     type="l",las=1,
     xlim=c(0,length(timevec)))
################################

##### Finding steady state #####
steadyState=runsteady(y = variables0, 
                      fun = SISmodel, 
                      parms = parameters)
steadyState$y
################################

##### Assessing the stability of the steady state #####
## Stability
Jmatrix=jacobian.full(y = steadyState$y, 
                      func = SISmodel, 
                      parms = as.numeric(parameters))
Eigenvalues=Re(eigen(Jmatrix)$values)
Eigenvalues

