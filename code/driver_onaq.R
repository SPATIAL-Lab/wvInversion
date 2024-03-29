library(R2jags)

#precip
p = read.csv("data/pcp_onaq_12hour_exdset3_05072020.csv")
##check for NAs
sum(is.na(p$bulk_precip))
##this is the vector of precip values, in mm
p = p$bulk_precip

#evapotranspiration
et = read.csv("data/le_nsae_onaq_12hour_exdset3_05072020.csv")
##complete timeseries, convert from w/m^2 to mm
etdf = data.frame("et" = et$mean_nsae_h2o * 2.257e-6 * 3600, 
                  "etsd" = rep(NA))
##made up sd for observed values
etdf$etsd[!is.na(etdf$et)] = 0.05
##fill missing values with means and sds of day and night obs
dmean = mean(etdf$et[seq(2, nrow(etdf), by = 2)], na.rm = TRUE)
dsd = sd(etdf$et[seq(2, nrow(etdf), by = 2)], na.rm = TRUE)
nmean = mean(etdf$et[seq(1, nrow(etdf), by = 2)], na.rm = TRUE)
nsd = sd(etdf$et[seq(1, nrow(etdf), by = 2)], na.rm = TRUE)
etmean = rep(c(nmean, dmean), nrow(etdf) / 2)
etsd = rep(c(nsd, dsd), nrow(etdf) / 2)
etdf$et[is.na(etdf$et)] = etmean[is.na(etdf$et)]
etdf$etsd[is.na(etdf$etsd)] = etsd[is.na(etdf$etsd)]
##make it a matrix
etp = as.matrix(etdf)

#Soil water content
swc = read.csv("data/swc_onaq_12hour_exdset3_05072020.csv")
##pull depths
ts = et$t
swc$ti = match(swc$t, ts)
swc$di = match(swc$v, c(501, 502, 503, 504, 505))
l = cbind(swc$ti, swc$di, swc$mean_vswc)

#Set up input
dat = list("delta_t" = 0.5, "ET" = etp*1e-3, "pre_pri" = p*1e-3,
           "phi.data" = l, "thick" = c(0.06, 0.1, 0.1, 0.2, 1),
           "alpha" = c(2, 2, 1, 1, 1), "nl" = 5, "nt" = length(p))

parms = c("phi", "adv", "diff", "et", "k_s", "phi_s", "psi_s", 
          "k_dt", "pre", "e_frac", "t_frac")

#Run inversion
pt = proc.time()
rmod = jags.parallel(model.file = "code/model_generic.R", 
                     parameters.to.save = parms, 
                     data = dat, inits = NULL, 
                     n.chains = 4, n.iter = 1000, n.burnin = 500, n.thin = 1)
(proc.time() - pt)[3]

#Shorthand for posterior parameters
sl = rmod$BUGSoutput$sims.list

#Some useful plots
##Compare soil phi timeseries
plot(apply(sl$phi[,,1], 2, mean), type = "l", ylim = c(0.05, 0.3))
lines(l[l[,2] == 1,1], l[l[,2] == 1,3], lty = 3)
for(i in 2:5){
  lines(apply(sl$phi[,,i], 2, mean), col = i)
  lines(l[l[,2] == i,1], l[l[,2] == i,3], lty = 3, col = i)
}

##Compare ET timeseries
plot(etp[,1]*1e-3, type = "l")
lines(2:length(p), apply(sl$et, 2, mean), col = "red")

##Compare PCP timeseries
plot(p*1e-3, type = "l")
lines(2:length(p), apply(sl$pre, 2, mean), col = "red")

##Soil phi w efrac overlay
par(mar= c(5,4,1,4))
plot(apply(sl$phi[,,1], 2, mean), type = "l", col = "red", ylim = c(0.05, 0.3), 
     ylab = "SWC")
lines(l[,1])
lines(apply(sl$phi[,,2], 2, mean), lty = 2, col = "red")
lines(l[,2], lty = 2)
par(new = TRUE)
plot(2:length(p), apply(sl$e_frac, 2, mean), type = "l", col = "blue", 
     lty = 3, xlab = "", ylab = "", axes = FALSE)
axis(4)
mtext("Evap fraction", side = 4, line = 2.2, col = "blue")
