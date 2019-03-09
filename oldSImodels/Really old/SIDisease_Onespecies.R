rm(list = ls())

library(deSolve)   
library(rootSolve)
library(FME)


##### Model description  #####
SISmodel=function(t,y,parameters) 
{ 
  ## Variables
  S=y[1]; I=y[2]; E=y[3];
  
  ## Parameters
  p = parameters[1];
  K = parameters[2];
  beta = parameters[3];
  mu = parameters[4];
  alpha = parameters[5];
  lambda = parameters[6];
  tau = parameters[7]  
  
  ## Ordinary differential equations
  logistic_growth <- (1 - (S + I)/K)
  dSdt <- p * logistic_growth * S - beta*S*E - S * mu
  dIdt <- beta*S*E - I * alpha - I * mu
  dEdt <- lambda * I - tau * E
  
  return(list(c(dSdt,
                dIdt,
                dEdt))); 
} 

## Initial state
variables0=c(S0=1000,I0=1,E0=10)

## Times at which estimates of the variables are returned
timevec=seq(0,100,5)

######################################################################################

#####  Intial condition  #####
## Parameters
parameters=c(p <- 10,
K <- 5000,
beta <- 0.005,
mu <- 0.05,
alpha <- 0.05,
lambda <- 0.001,
tau <- 0.2
)

###############################

#####  Running the model   #####
## Numerical integration through time
output=lsoda(y = variables0,    # intial values  
             times = timevec,   # time vector
             func = SISmodel,   # model
             parms = parameters # constant parameters
)

## Plotting
colnames(output)=c("time","S","I","E")
infectionDynamics=as.data.frame(output)
par(mfrow=c(1,2),mar=c(5,5,2,2))

#Susceptible hosts
plot(infectionDynamics$S~infectionDynamics$time,
     main="Susceptibles",ylab="Density",xlab="Time (days)",
     type="l",las=1,
     xlim=c(0,length(timevec)))

#Infected hosts
plot(infectionDynamics$I~infectionDynamics$time,
     main="Infecteds",ylab="Density",xlab="Time (days)",
     type="l",las=1,
     xlim=c(0,length(timevec)))

##### Finding steady state #####
steadyState=runsteady(y = variables0, 
                      fun = SISmodel, 
                      parms = parameters)
steadyState$y


#Environmental parasite
plot(infectionDynamics$E~infectionDynamics$time,
     main="Env Parasite",ylab="Density",xlab="Time (days)",
     type="l",las=1,
     xlim=c(0,length(timevec)))
################################



################################

##### Assessing the stability of the steady state #####
## Stability
Jmatrix=jacobian.full(y = steadyState$y, 
                      func = SISmodel, 
                      parms = as.numeric(parameters))
Eigenvalues=Re(eigen(Jmatrix)$values)
Eigenvalues

