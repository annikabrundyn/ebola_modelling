###### Sierra Leone: Vaccination model ------

# set working directory, read in saved parameter estimates
wd = "~/Desktop/Project/Drafts/mathematical-modelling-spread/code/sl"
setwd(wd)
model.estms = read.csv("final_sl.csv", sep = " ", skip = 1, header = FALSE)
model.estms = setNames( model.estms$V2, model.estms$V1)  #turn df into named vector
model.estms

# required libraries
library(deSolve)
library(ggplot2)

#SVEIRDV
sveirdb = function(t, x, parms){
  with(as.list(c(parms,x)), {
    
    N = S + E + I + R + D + Vac
    
    dS = - betaI*(I/N)*S - betaD*(D/N)*S - omega*S
    
    dE = betaI*(I/N)*S + betaD*(D/N)*S - alpha*E
    
    dI = alpha*E - (1 - mu)*(lambda1)*I - mu*(lambda2)*I
    
    dR = (1 - mu)*lambda1*I
    
    dD = mu*lambda2*I - rho*D
    
    dB = rho*D
    
    dInc = betaI*(I/N)*S + betaD*(D/N)*S
    
    dVac = omega*S
    
    output = c(dS, dE, dI, dR, dD, dB, dInc, dVac)
    list(output)
  })
}

InitPop = 6092000
E0 = 2
I0 = 3
R0 = 0
D0 = 0
B0 = 0
Inc0 = 5 
Vac0 = 0
S0 = InitPop - E0 - I0 - R0 - D0 - B0 - Vac0
start = c(S = S0, E = E0, I = I0, R = R0, D = D0, B = B0, Inc = Inc0, Vac = Vac0)

# create model times vector
startday = 1
endday = 800   #Data stops at day 682
times = seq(startday, endday, 1)    


# plot omega ranges
omega_range = c(0.0006, 0.0008, 0.001, 0.0015)
model.estms2 = model.estms
par(mfrow = c(2,2))

for ( i in 1:length(omega_range)){
  model.estms2['omega'] = omega_range[i]
  fit = data.frame(ode(times = times, y = start, func = sveirdb, parms = model.estms2))
  fit$propv = fit[,9]/InitPop
  plot(fit$propv, fit[,8], type = "l", ylim = c(0,4000), xlim = c(0,0.6), 
       main = bquote(omega == .(format(omega_range[i], scientific = FALSE))), 
       xlab = "Proportion of total population vaccinated", ylab = "Cumulative Cases", 
       col = "#00AFBB", lwd = 1.2, cex.main = 1.5, cex.lab = 1.25, cex.axis = 1.25)
  print(c(omega_range[i], fit[800,8]))
}


