rm(list = ls())

require(deSolve)   
require(rootSolve)
require(FME)
require(ggplot2)
#setwd("C:/Users/fortinlab/Dropbox/Courses/2016-2017/EEB1451 - Disease in Communities/Term_paper")
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
  S=y[1]; I=y[2]; R=y[3];
  
  ## Parameters
  p = parameters[1];
  K = parameters[2];
  beta = parameters[3];
  mu = parameters[4];
  alpha = parameters[5];
  gamma = parameters[6]  
  
  ## Ordinary differential equations
  logistic_growth <- (1 - (S + I + R)/K)
  dSdt <- p*(S+I+R)*logistic_growth - beta*S*I - S*mu
  dIdt <- beta*S*I - I*(alpha + mu + gamma)
  dRdt <- I*gamma - R*mu
  
  return(list(c(dSdt,
                dIdt,
                dRdt))); 
} 

###Paramater values (expect p and beta because that's what we're testing)
K <- 5000
mu <- 0.02
alpha <- 0.4
gamma <- 0.01

## Initial state
variables0=c(S0=1000,I0=1,R=0)


################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
### [2] Linear relationship between p and beta

#Looping to get range of beta values
Beta_values <- seq(0.001, 0.2, by = 0.001)
Suscep_num <- NULL
Infect_num <- NULL
Recovered_num <- NULL

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
               gamma
               )
  
  steadyState=runsteady(y = variables0, 
                      fun = SISmodel, 
                      parms = parameters)
  
  Suscep <- steadyState$y[1]
  Suscep_num <- append(Suscep_num, Suscep)
  
  Infect <- steadyState$y[2]
  Infect_num <- append(Infect_num, Infect)
  
  Recover <- steadyState$y[3]
  Recovered_num <- append(Recovered_num, Recover)
  }

#Getting data into dataframe for plotting - figure out how to do in one line
tempdf <- cbind(Beta_values, Suscep_num, Infect_num, Recovered_num)
Linear_data <- as.data.frame(tempdf, row.names = T)

#Finding optimal beta value for susceptibles
Opt.row <- which(Linear_data$Suscep_num == max(Linear_data$Suscep_num), arr.ind=TRUE)
Opt.beta.value <- Linear_data$Beta_values[Opt.row]

Opt.Infect.beta <- which(Linear_data$Infect_num == max(Linear_data$Infect_num), arr.ind=TRUE)
Opt.Infect.beta <- Linear_data$Beta_values[Opt.Infect.beta]
View(Linear_data)


#Building plot with data
plot_title <- paste0("Optimal beta value =", Opt.beta.value, ", alpha =", alpha, ", p_constant = ", p_constant, "beta^", beta_exp)
plot_title

Linear_plot <- ggplot(Linear_data, aes(Beta_values)) +
  geom_line(aes(y = Suscep_num, colour = "Susceptibles")) + 
  geom_line(aes(y = Infect_num, colour = "Infecteds")) +
  geom_line(aes(y = Recovered_num, colour = "Recovered")) +
  labs(title = plot_title, x = "Beta values", y = "Population size")

#Setting theme for plot
Linear_plot + theme(
  plot.title = element_text(size = 16), 
  axis.text = element_text(size = 16),
  axis.title = element_text(size = 16),
  axis.line = element_line("black"),
  legend.title = element_blank(),
  panel.background = element_blank()
)
getwd()
ggsave(paste0("./plots/SIR_direct/", plot_title, ".png"), width = 10, height = 8)


##########
#########
########
#######
######
#####
###
##
#END