rm(list = ls())

library(deSolve)   
library(rootSolve)
library(FME)


##### Model description  #####
SISmodel=function(t,y,parameters) 
{ 
  ## Variables
  S=y[1]; I=y[2]; R=y[3]; E=y[4];
  
  ## Parameters
  birth_rate = parameters[1];
  K = parameters[2];
  beta = parameters[3];
  natural_mortality = parameters[4];
  disease_morality = parameters[5];
  recovered = parameters[6];
  lambda = parameters[7];
  mu = parameters[8] 
    
  
  ## Ordinary differential equations
  dSdt <- birth_rate * (S+I+R) * (K - (S+I+R)/K) - beta*S*E - S * natural_mortality
  dIdt <- beta * S * E - I * disease_morality - I * natural_mortality
  dRdt <- recovered * I - R * natural_mortality
  dEdt <- lambda * I - mu * E
  
  return(list(c(dSdt,
                dIdt,
                dRdt,
                dEdt))); 
} 


#####  Intial condition  #####
## Parameters
parameters=c(birth_rate <- 1.5,
K <- 5000,
beta <- 0.005,
natural_mortality <- 0.001,
disease_morality <- 0.005,
recovered  <- 0.005,
lambda <- 0.01,
mu <- 0.2)



## Initial state
variables0=c(S0=999,I0=1,R0=0,E0=10)

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
colnames(output)=c("time","S","I","E")
infectionDynamics=as.data.frame(output)
par(mfrow=c(1,2),mar=c(5,5,2,2))

#Susceptible hosts
plot(infectionDynamics$S~infectionDynamics$time,
     main="Susceptibles",ylab="Density",xlab="Time (days)",
     type="l",las=1,
     ylim = c(0,1000), xlim=c(0,length(timevec)))

#Infected hosts
plot(infectionDynamics$I~infectionDynamics$time,
     main="Infecteds",ylab="Density",xlab="Time (days)",
     type="l",las=1,
     xlim=c(0,length(timevec)))

#Recovered hosts
plot(infectionDynamics$R~infectionDynamics$time,
     main="Recovered",ylab="Density",xlab="Time (days)",
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

