#extra code that may be useful

######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
#OLD STUFF

#Trying to loop across multiple p constants


for (z in c(1/100, 1/10, 1, 10, 100)){
  Beta_values <- seq(0.005, 0.2, by = 0.005)
  Suscep_num <- NULL
  Infect_num <- NULL
  Environ_num <- NULL
  
  for (i in Beta_values){
    beta <- i
    p = z*beta
    
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
  
  temp_name <- paste0("Suscep_num", z)
  assign(temp_name, Suscep_num) 
  
  temp_name <- paste0("Infect_num", z)
  assign(temp_name, Infect_num)
  
  temp_name <- paste0("Environ_num", z)
  assign(temp_name, Environ_num)
  
}




#####  Running the model   #####

## Times at which estimates of the variables are returned
timevec=seq(0,100,5)

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

