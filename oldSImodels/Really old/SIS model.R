rm(list = ls())

library(deSolve)   
library(rootSolve)
library(FME)

#####  Model description  #####
SISmodel=function(t,y,parameters) 
{ 
  ## Variables
  S=y[1]; I=y[2];
  
  ## Parameters
  beta = parameters[1]; gamma = parameters[2]
  
  ## Ordinary differential equations
  dSdt = - beta * S * I + gamma * I 
  dIdt = beta * S * I - gamma * I
  
  return(list(c(dSdt,dIdt))); 
} 
###############################

#####  Intial condition  #####
## Parameters
parameters=c(beta=0.001, gamma=1/10)

## Initial state
variables0=c(S0=999,I0=1)

## Times at which estimates of the variables are returned
timevec=seq(0,30,1)
###############################

#####  Running the model   #####
## Numerical integration through time
output=lsoda(y = variables0,    # intial values  
             times = timevec,   # time vector
             func = SISmodel,   # model
             parms = parameters # constant parameters
)

## Plotting
colnames(output)=c("time","S","I")
infectionDynamics=as.data.frame(output)

par(mfrow=c(1,2),mar=c(5,5,2,2))
plot(infectionDynamics$S~infectionDynamics$time,
     main="Susceptibles",ylab="Density",xlab="Time (days)",
     type="l",las=1,
     ylim=c(0,1000),xlim=c(0,length(timevec)))
plot(infectionDynamics$I~infectionDynamics$time,
     main="Infecteds",ylab="Density",xlab="Time (days)",
     type="l",las=1,
     ylim=c(0,1000),xlim=c(0,length(timevec)))
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

