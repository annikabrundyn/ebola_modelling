###### Liberia: Model fitting, parameter estimation ------

# clear environment and set working directory
rm(list = ls())
wd = "~/Liberia code"
setwd(wd)

# install required libraries
library(deSolve)
library(gtools)    # for logit function
library(ggplot2)
library(gridExtra)

# read in data, create dates sequence
lib.data = readRDS("lib_data_23Sept.rds")
data.f = lib.data[,c(2:4)]
dates.seq = data.f$Day # save dates for which case counts are available 
death.dates.seq = data.f$Day[!is.na(data.f$Deaths)] # save dates for which death counts are available

# starting date: "2014-06-02" = day 1

# set initial values
InitPop = 4294000
E0 = 80
I0 = 14
R0 = 0
D0 = 3
B0 = 34
Inc0 = 51
S0 = InitPop - E0 - I0 - R0 - D0 - B0
start = c(S = S0, E = E0, I = I0, R = R0, D = D0, B = B0, Inc = Inc0)

# create model times vector
startday = 1
endday = 800   # data stops at day
times = seq(startday, endday, 1)                                        

#### Define functions -----

# SEIRDB function for estimation
seirdb.est = function(t, x, parms){
  with(as.list(c(parms,x)), {
    betaI = exp(logbetaI)           #effective contact rate with infectious people (alive)
    betaD = exp(logbetaD)           #effective contact rate with dead but infectious people
    alpha = exp(logalpha)           #1/latency period
    lambda1 = exp(loglambda1)       #1/period of infection to survival - still infectious
    lambda2 = exp(loglambda2)       #1/period of infection to death
    rho = exp(logrho)               #1/time to dispose of a body
    mu = inv.logit(logitmu)         #fatality rate
    eta = inv.logit(logiteta)       #factor to decrease betaI for t > tc
    tc = exp(logtc)                 #time of intervention/control measures implemented
    
    if (t >= tc) {
      etat = eta
    } else {
      etat = 1
    }
    
    N = S + E + I + R + D
    
    dS = - betaI*etat*(I/N)*S - betaD*etat*(D/N)*S
    
    dE = betaI*etat*(I/N)*S + betaD*etat*(D/N)*S - alpha*E
    
    dI = alpha*E - (1 - mu)*(lambda1)*I - mu*(lambda2)*I
    
    dR = (1 - mu)*lambda1*I
    
    dD = mu*lambda2*I - rho*D
    
    dB = rho*D
    
    dInc = betaI*etat*(I/N)*S + betaD*etat*(D/N)*S
    
    output = c(dS, dE, dI, dR, dD, dB, dInc)
    list(output)
  })
}

# SEIRD function - without transformed parameters for estimation
seirdb = function(t, x, parms){
  with(as.list(c(parms,x)), {
    
    if (t >= tc) {
      etat = eta
      #rhot = rho2
    } else {
      etat = 1
      #rhot = rho1
    }
    
    N = S + E + I + R + D
    
    dS = - betaI*etat*(I/N)*S - betaD*etat*(D/N)*S
    
    dE = betaI*etat*(I/N)*S + betaD*etat*(D/N)*S - alpha*E
    
    dI = alpha*E - (1 - mu)*(lambda1)*I - mu*(lambda2)*I
    
    dR = (1 - mu)*lambda1*I
    
    dD = mu*lambda2*I - rho*D
    
    dB = rho*D
    
    dInc = betaI*etat*(I/N)*S + betaD*etat*(D/N)*S
    
    output = c(dS, dE, dI, dR, dD, dB, dInc)
    list(output)
  })
}


# Function for calculating sum of squared errors from case and death data
seirdb.sse = function(varparms, fixparms, times, start, data) {
  seirdb.lse = ode(times = times, y = start, func = seirdb.est, parms = c(varparms, fixparms))
  
  error.cum.cases = (seirdb.lse[dates.seq, 8] - data$Cases)^2
  error.cum.cases[-c(1:48)] = 2*error.cum.cases[-c(1:48)]     #values after 180 days is observation 48
  
  error.cum.deaths = (seirdb.lse[death.dates.seq, 7] - data$Deaths[data$Day %in% death.dates.seq])^2
  
  sse.cases = sum(error.cum.cases)
  sse.deaths = sum(error.cum.deaths)
  sse = sse.cases + 1.5*sse.deaths
  
  return(sse)
}

# Base reproductive number
r0.fn2 = function(estms){
  betaI = estms['betaI']
  betaD = estms['betaD']
  lambda1 = estms['lambda1']
  lambda2 = estms['lambda2']
  mu = estms['mu']
  rho = estms['rho']
  r0 = ( (betaI / (lambda1 + mu*(lambda2 - lambda1))) + ((betaD * mu * lambda2)/(rho*(lambda1 + mu*(lambda2 - lambda1))) ) )
  return(r0)
}

#### Estimate parameter values from data -----

# only need to run this the first time to initialize values for total error and estimates
min.err = 1000000000000
min.start = min.estms = NULL

# set number of iterations to run with different starting values
nsim = 500

# provide values for fixed parameters
fixparms = c(logalpha = log(1/10),   
             loglambda1 = log(1/9.4), 
             loglambda2 = log(1/7.5))

#start_time <- Sys.time()  #to measure run time

# create for loop to generate random starting values for variable parameters,
# optimize variable parameters using L-BFGS-B method
# calculate total error with these model parameters, if the total error is less
# than the current saved minimum error, save the estimated parameters as the 
# best parameters ('min.estms')

for (i in 1:nsim) {
  
  varparms = c(logbetaI = log(runif(1, 0.15, 0.2)), 
               logbetaD = log(runif(1, 0.1, 0.4)),
               logrho = log(runif(1, 0.25, 1.5)), 
               logtc = log(runif(1, 40, 120)),
               logiteta = logit(runif(1, 0.01, 0.6)),
               logitmu = logit(runif(1, 0.25, 0.7)))
  
  sl.optim = optim(par = varparms, seirdb.sse, fixparms = fixparms, method = "L-BFGS-B",
                   times = times, start = start, data = data.f,
                   lower = c(-10, -10, log(0.2), -10, -10, logit(0.2)),
                   upper = c(log(2), log(2), log(1.5), log(300), logit(0.8), logit(0.8)))
  sl.sse = sl.optim$value
  
  if (sl.sse < min.err) {
    min.err = sl.sse
    min.start = varparms
    min.estms = sl.optim$par
  }
}

#end_time <- Sys.time()
#end_time - start_time


# back transform estimates from log/logit scale to original scale

all.estms = c( fixparms, min.estms)
all.estms = setNames(all.estms, c("alpha", "lambda1", "lambda2",
                                  "betaI", "betaD",  "rho", "tc", "eta", "mu"))
all.estms
estms.no.tr = all.estms

all.estms[c("alpha", "lambda1", "lambda2",
            "betaI", "betaD","rho", "tc")] = sapply(all.estms[c("alpha", "lambda1", "lambda2",
                                                                "betaI", "betaD","rho", "tc")], exp)
all.estms[c("eta", "mu") ] = sapply(all.estms[c("eta", "mu")], inv.logit)
all.estms

# save estimates in rds and csv format
saveRDS(all.estms, file = 'lib_final_24Sept.rds')
write.table(all.estms, file = "final_lib.csv")

model.estms = all.estms
model.estms

# calculate r0 value
r0_lib = r0.fn2(model.estms)
r0_lib

# fit model using estimated parameters
model.f = ode(times = times, y = start, func = seirdb, parms = model.estms)

#### plot predicted behaviour of model compartments -----

model.f2 = as.data.frame(model.f)

# cumulative cases
pInc = ggplot(data.f, aes(x = Day, y = Cases)) +
  geom_point(shape = 1, color="gray35") + 
  geom_line(data = model.f2, aes(x = time, y = Inc), col = "tomato", size = 0.7) +
  ggtitle("Liberia Cumulative Cases")
#theme_minimal()
pInc

#deaths
pB = ggplot(data.f, aes(x = Day, y = Deaths)) +
  geom_point(shape = 1, color="gray35") + 
  geom_line(data = model.f2, aes(x = time, y = B), col = "tomato", size = 0.7) +
  ggtitle("Liberia Cumulative Deaths (Class B)") 
#theme_minimal()
pB

#incidence
pI = ggplot(lib.data, aes(x = Day, y = inc)) +
  geom_point(shape = 1, color="gray35")  + ylab("Cases") +
  geom_line(data = model.f2, aes(x = time, y = I), col = "tomato", size = 0.7) +
  ggtitle("Liberia Case Incidence (Class I)")
#theme_minimal()
pI

#grid.arrange(pInc, pB, pI, ncol = 1)

#S
pS = ggplot(model.f2, aes(x = time, y = S)) + 
  geom_line(color = "tomato", size = 1) +
  xlab("Days") + ylab("S") +
  theme_minimal()
#theme(axis.text.x=element_blank())

#E
pE = ggplot(model.f2, aes(x = time, y = E)) + 
  geom_line(color = "tomato", size = 1) +
  xlab("Days") + ylab("E") +
  theme_minimal()

#R
pR = ggplot(model.f2, aes(x = time, y = R)) + 
  geom_line(color = "tomato", size = 1) +
  xlab("Days") + ylab("R") +
  theme_minimal()

#D
pD = ggplot(model.f2, aes(x = time, y = D)) + 
  geom_line(color = "tomato", size = 1) +
  xlab("Days") + ylab("D") +
  theme_minimal()

grid.arrange(pS, pE, pR, pD, nrow = 2, ncol = 2)




