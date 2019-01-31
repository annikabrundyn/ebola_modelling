###### Sensitivity analysis for Liberia ------

# set working directory, read in saved parameter estimates
wd = "~/Desktop/Project/Drafts/mathematical-modelling-spread/code/lib"
setwd(wd)
model.estms = read.csv("final_lib.csv", sep = " ", skip = 1, header = FALSE)
model.estms = setNames( model.estms$V2, model.estms$V1)  #turn df into named vector
model.estms

# install required libraries
library(deSolve)
library(lhs)
library(ppcor)
library(data.table)
library(ggplot2)

#### Defining functions -----

# base reproductive number
r0.calc = function(params){
  
  betaI = params['betaI']
  betaD = params['betaD']
  lambda1 = params['lambda1']
  lambda2 = params['lambda2']
  mu = params['mu']
  rho = params['rho']
  
  r0 = ( (betaI / (lambda1 + mu*(lambda2 - lambda1))) + ((betaD * mu * lambda2)/(rho*(lambda1 + mu*(lambda2 - lambda1))) ) )
  return(r0)
}

#### Univariate relationship between R0 and key parameters -----

## betaI
betaI_range = seq(0.12, 0.22, length.out = 500)
model.estms2 = model.estms
betaI_r0_range = NULL
for (i in 1:length(betaI_range)){
  model.estms2['betaI'] = betaI_range[i]
  betaI_r0_range = c(betaI_r0_range, r0.calc(model.estms2))
}
betaI_info = data.frame(betaI = betaI_range, r0 = betaI_r0_range)

r0_betaI = ggplot(betaI_info, aes(x = betaI, y = r0)) +
  geom_line(size = 1.2, colour = "tomato") + xlab(expression(beta[I])) + ylab(expression('R'[0])) +
  theme_minimal() +
  theme( axis.title.x = element_text(size=rel(1.2)),
         axis.title.y = element_text(size=rel(1.2)))
#theme(axis.title.y = element_text(angle=0)) 
#r0_betaI


## betaD
betaD_range = seq(0.24, 0.42, length.out = 500)
model.estms2 = model.estms
betaD_r0_range = NULL
for (i in 1:length(betaD_range)){
  model.estms2['betaD'] = betaD_range[i]
  betaD_r0_range = c(betaD_r0_range, r0.calc(model.estms2))
}
betaD_info = data.frame(betaD = betaD_range, r0 = betaD_r0_range)

r0_betaD = ggplot(betaD_info, aes(x = betaD, y = r0)) +
  geom_line(size = 1.2, colour = "tomato") + xlab(expression(beta[D])) + ylab(expression('R'[0])) +
  theme_minimal() +
  theme( axis.title.x = element_text(size=rel(1.2)),
         axis.title.y = element_text(size=rel(1.2)))
#theme(axis.title.y = element_text(angle=0)) 
#r0_betaD


## lambda1
lambda1_range = seq(0.06, 0.17, length.out = 500)
model.estms2 = model.estms
lambda1_r0_range = NULL
for (i in 1:length(lambda1_range)){
  model.estms2['lambda1'] = lambda1_range[i]
  lambda1_r0_range = c(lambda1_r0_range, r0.calc(model.estms2))
}
lambda1_info = data.frame(lambda1 = lambda1_range, r0 = lambda1_r0_range)

r0_lambda1 = ggplot(lambda1_info, aes(x = lambda1, y = r0)) +
  geom_line(size = 1.2, colour = "tomato") + xlab(expression(lambda[1])) + ylab(expression('R'[0])) +
  theme_minimal() +
  theme( axis.title.x = element_text(size=rel(1.2)),
         axis.title.y = element_text(size=rel(1.2)))
#theme(axis.title.y = element_text(angle=0)) 
#r0_lambda1


## lambda2
lambda2_range = seq(0.06, 0.17, length.out = 500)
model.estms2 = model.estms
lambda2_r0_range = NULL
for (i in 1:length(lambda2_range)){
  model.estms2['lambda2'] = lambda2_range[i]
  lambda2_r0_range = c(lambda2_r0_range, r0.calc(model.estms2))
}
lambda2_info = data.frame(lambda2 = lambda2_range, r0 = lambda2_r0_range)

r0_lambda2 = ggplot(lambda2_info, aes(x = lambda2, y = r0)) +
  geom_line(size = 1.2, colour = "tomato") + xlab(expression(lambda[2])) + ylab(expression('R'[0])) +
  theme_minimal() +
  theme( axis.title.x = element_text(size=rel(1.2)),
         axis.title.y = element_text(size=rel(1.2)))
#theme(axis.title.y = element_text(angle=0)) 
#r0_lambda2


## mu
mu_range = seq(0.2, 0.7, length.out = 500)
model.estms2 = model.estms
mu_r0_range = NULL
for (i in 1:length(mu_range)){
  model.estms2['mu'] = mu_range[i]
  mu_r0_range = c(mu_r0_range, r0.calc(model.estms2))
}
mu_info = data.frame(mu = mu_range, r0 = mu_r0_range)

r0_mu = ggplot(mu_info, aes(x = mu, y = r0)) +
  geom_line(size = 1.2, colour = "tomato") + xlab(expression(mu)) + ylab(expression('R'[0])) +
  theme_minimal() +
  theme( axis.title.x = element_text(size=rel(1.2)),
         axis.title.y = element_text(size=rel(1.2)))
#theme(axis.title.y = element_text(angle=0)) 
r0_mu


## rho
rho_range = seq(0.33, 1.5, length.out = 500)
model.estms2 = model.estms
rho_r0_range = NULL
for (i in 1:length(rho_range)){
  model.estms2['rho'] = rho_range[i]
  rho_r0_range = c(rho_r0_range, r0.calc(model.estms2))
}
rho_info = data.frame(rho = rho_range, r0 = rho_r0_range)

r0_rho = ggplot(rho_info, aes(x = rho, y = r0)) +
  geom_line(size = 1.2, colour = "tomato") + xlab(expression(rho)) + ylab(expression('R'[0])) +
  theme_minimal() +
  theme( axis.title.x = element_text(size=rel(1.2)),
         axis.title.y = element_text(size=rel(1.2)))
#theme(axis.title.y = element_text(angle=0)) 
#r0_rho


grid.arrange(r0_betaI, r0_betaD, r0_lambda1, r0_lambda2, r0_mu, r0_rho, ncol = 2, nrow = 3)

#### Assessing impact of earlier intervention (changes in tC) -----

#4 plots - 3 monthly changes
tc_range = c(13, 43, 73, 103)
model.estms2 = model.estms
par(mfrow = c(2,2))

for ( i in 1:length(tc_range)){
  model.estms2['tc'] = tc_range[i]
  fit = ode(times = times, y = start, func = seirdb, parms = model.estms2)
  diff_inc = c(NA, fit[2:nrow(fit),8] - fit[1:(nrow(fit)-1),8] )
  #peak = which.max(diff_inc) #which(diff_inc <= 0)[1]
  slow_down = which(diff_inc < 1)[1]
  plot(fit[,1], fit[,8], type = "l", ylim = c(0,12000), xlim = c(0,800), 
       main = paste("tc =", tc_range[i]), xlab = "Days", ylab = "Cumulative Cases", col = "tomato", lwd = 1.2)
  abline(v = tc_range[i], lty = 2, col = "darkblue")
  abline(v = slow_down, lty = 2, col = "darkred")
  print(c(tc_range[i], slow_down, fit[800,8]))
}

# now look at relationship numerically
tc_diff_range = seq(0, 90, 0.25)
cases_avoided = NULL
total_cases = model.f[800,8]

for (i in 1:length(tc_diff_range)){
  model.estms2['tc'] = model.estms['tc'] - tc_diff_range[i]
  fit = ode(times = times, y = start, func = seirdb, parms = model.estms2)
  cases_avoided = c(cases_avoided, (total_cases - fit[800,8]))
}
df.cases_avoided = data.frame("days" = tc_diff_range, "cases" = cases_avoided)

tc_early = ggplot(df.cases_avoided, aes(x = days, y = cases, group = 1)) + 
  geom_line(size = 1.2, color = "tomato") + xlab("Days sooner") + ylab("Cases avoided") +
  ggtitle("Liberia") + 
  theme_minimal()

tc_early

#### Relationship between 'disease threshold' and duration of the epidemic -----
tc_range2 = seq(30, 103, 1)
model.estms2 = model.estms
threshold_vals = c()
duration_vals = c()

for (i in 1:length(tc_range2)){
  model.estms2['tc'] = tc_range2[i]
  fit = ode(times = times, y = start, func = seirdb, parms = model.estms2)
  diff_inc = c(NA, fit[2:nrow(fit),8] - fit[1:(nrow(fit)-1),8] )
  #peak = which.max(diff_inc) #which(diff_inc <= 0)[1]
  slow_down = which(diff_inc < 1)[1]
  
  threshold_vals = c(threshold_vals, fit[tc_range2[i], 8])
  
  duration_vals = c(duration_vals, (slow_down))
}
thresh_data = data.frame(threshold = threshold_vals, duration = duration_vals)

p_thresh = ggplot(thresh_data, aes(x = threshold, y = duration)) + 
  geom_line(size = 1.2, color = "tomato") + xlab("Threshold (number of cases)") + ylab("Duration of epidemic") +
  ggtitle( expression( paste("Liberia (", eta, " = 0.502)"))) + 
  scale_x_continuous(breaks = seq(300, 6300, 1000)) + 
  scale_y_continuous(minor_breaks = seq(0, 800, 25), breaks = seq(200, 800, 50)) +
  theme_minimal() +
  theme(text = element_text(size = 14))
p_thresh


#### Multivariate sensitivity analysis on R0 -----

n = 10000# number to sample
k = 6 #number of variables

lhsp = randomLHS(n, k)

samples = data.frame( betaI = qunif( lhsp[,1], 0.12, 0.2),   #check these two
                      betaD = qunif( lhsp[,2], 0.24, 0.41),     #check
                      lambda1 = qnorm( lhsp[,3], 0.1, 0.0225),
                      lambda2 = qnorm( lhsp[,4], 0.1, 0.0225),
                      mu = qunif( lhsp[,5], 0.2, 0.7),
                      rho = qunif( lhsp[,6], 0.25, 1.5)
)

save.vals = NULL

for (i in 1:n) {
  #run = ode( times = times, y = start, func = seirdb, parms = samples[i,])
  save.vals[i] = r0.calc(samples[i,])
}

save.vals = unlist(save.vals)

save.vals2 = data.frame(values = save.vals)
ggplot(save.vals2, aes(values)) + 
  geom_histogram(binwidth = 0.2, fill = "tomato", color = "grey45", alpha = 0.5) +
  xlab(expression("R"[0])) + ylab("Frequency") +
  ggtitle(expression(paste("Uncertainty analysis of ", "R"[0], " for Liberia"))) +
  theme_minimal()

t.test(save.vals)
summary(save.vals)

#hist(save.vals, main = expression(paste("Uncertainty analysis of ", "R"[0], " for Sierra Leone")), xlab = expression("R"[0]))

pcor_vals = c()

for ( i in 1:k){
  vals = pcor.test(x = samples[,i], y = save.vals, z = samples[,-i], method = "spearman")
  pcor_vals = c(pcor_vals, vals$estimate)
}

pcor_vals
pcor_vals = setNames(pcor_vals, colnames(samples))
barplot(pcor_vals, ylim = c(-1, 1))



#### Univariate sensitivity analysis on cumulative cases -----

# Using LHS
n = 10000 # number to sample
k = 9 #number of variables

lhsp = randomLHS(n, k)

samples = data.frame( alpha = qnorm( lhsp[,1], 0.1, 0.0125),
                      lambda1 = qnorm( lhsp[,2], 0.1, 0.0225),
                      lambda2 = qnorm( lhsp[,3], 0.1, 0.0225),
                      betaI = qunif( lhsp[,4], 0.12, 0.22),
                      betaD = qunif( lhsp[,5], 0.24, 0.42),
                      rho = qunif( lhsp[,6], 0.25, 1.5),
                      tc = qunif( lhsp[,7], 73, 133),
                      eta = qunif( lhsp[,8], 0.3, 0.7),
                      mu = qunif( lhsp[,9], 0.2, 0.7)
                      #E0 = qunif( lhsp[,9], 30, 70),
                      #I0 = qunif( lhsp[,10], 10, 40),
                      #D0 = qunif( lhsp[,11], 10, 30)
)

save.vals50 = NULL
save.vals200 = NULL
save.vals400 = NULL
save.vals800 = NULL

for (i in 1:n) {
  run = ode( times = times, y = start, func = seirdb, parms = samples[i,])
  save.vals50[i] = run[50, 8]
  save.vals200[i] = run[200,8]
  save.vals400[i] = run[400,8]
  save.vals800[i] = run[800,8]
}

save.vals50 = unlist(save.vals50)
save.vals200 = unlist(save.vals200)
save.vals400 = unlist(save.vals400)
save.vals800 = unlist(save.vals800)
#hist(save.vals, main = expression(paste("Uncertainty analysis of ", "R"[0], " for Liberia")), xlab = expression("R"[0]))

pcor_vals50 = c()
pcor_vals200 = c()
pcor_vals400 = c()
pcor_vals800 = c()

for ( i in 1:k){
  vals50 = pcor.test(x = samples[,i], y = save.vals50, z = samples[,-i], method = "spearman") #removed method = "s"
  vals200 = pcor.test(x = samples[,i], y = save.vals200, z = samples[,-i], method = "spearman")
  vals400 = pcor.test(x = samples[,i], y = save.vals400, z = samples[,-i], method = "spearman")
  vals800 = pcor.test(x = samples[,i], y = save.vals800, z = samples[,-i], method = "spearman")
  pcor_vals50 = c(pcor_vals50, vals50$estimate)
  pcor_vals200 = c(pcor_vals200, vals200$estimate)
  pcor_vals400 = c(pcor_vals400, vals400$estimate)
  pcor_vals800 = c(pcor_vals800, vals800$estimate)
}

#pcor_vals
pcor_vals50 = setNames(pcor_vals50, colnames(samples))
pcor_vals200 = setNames(pcor_vals200, colnames(samples))
pcor_vals400 = setNames(pcor_vals400, colnames(samples))
pcor_vals800 = setNames(pcor_vals800, colnames(samples))

vals800 = data.frame(pcor_vals800)
setDT(vals800, keep.rownames = TRUE)[]
colnames(vals800) = c("variable", "value")

# try do with ggplot
#need to put into long format
sens_bar800 = ggplot(vals800, aes(x = variable, y = value)) +
  geom_bar(stat = "identity", fill = "tomato", col = "grey45", alpha = 0.5) + 
  ylab("PRCC") + xlab("Parameter") + ggtitle("Sensitivity analysis of cumulative cases at day 800 for Liberia") +
  theme_minimal()

sens_bar800

#other 4 plots
vals50 = data.frame(pcor_vals50)
setDT(vals50, keep.rownames = TRUE)[]
colnames(vals50) = c("variable", "value")

sens_bar50 = ggplot(vals50, aes(x = variable, y = value)) +
  geom_bar(stat = "identity", fill = "tomato", col = "grey45", alpha = 0.5) + 
  ylab("PRCC") + xlab("Parameter") + ggtitle("Day 50") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 7.5))

#sens_bar50

vals200 = data.frame(pcor_vals200)
setDT(vals200, keep.rownames = TRUE)[]
colnames(vals200) = c("variable", "value")

sens_bar200 = ggplot(vals200, aes(x = variable, y = value)) +
  geom_bar(stat = "identity", fill = "tomato", col = "grey45", alpha = 0.5) + 
  ylab("PRCC") + xlab("Parameter") + ggtitle("Day 200") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 7.5))

#sens_bar200

vals400 = data.frame(pcor_vals400)
setDT(vals400, keep.rownames = TRUE)[]
colnames(vals400) = c("variable", "value")

sens_bar400 = ggplot(vals400, aes(x = variable, y = value)) +
  geom_bar(stat = "identity", fill = "tomato", col = "grey45", alpha = 0.5) + 
  ylab("PRCC") + xlab("Parameter") + ggtitle("Day 400") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 7.5))


sens_bar800.v2 = ggplot(vals800, aes(x = variable, y = value)) +
  geom_bar(stat = "identity", fill = "tomato", col = "grey45", alpha = 0.5) + 
  ylab("PRCC") + xlab("Parameter") + ggtitle("Day 800") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 7.5))


grid.arrange(sens_bar50, sens_bar200, sens_bar400, sens_bar800.v2, nrow = 2, ncol = 2)

