library(R2jags)

#precip - prior
p = read.csv("data/pcp_onaq_12hour_exdset3_05072020.csv")
##check for NAs
sum(is.na(p$bulk_precip))
##timeseries for all other data
ts = p$t
##this is the vector of precip values, in mm
p = p$bulk_precip

#evapotranspiration - prior
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

#Soil temperature - prior
st = read.csv("data/st_onaq_12hour_05072020.csv")
st = st$mean_temp + 273

#Relative humidity - prior
rh = read.csv("data/rh_onaq_12Hour_05072020.csv")
rh = rh$mean_rh / 100

#Precipitation d18O - prior
p_o_pri = read.csv("data/d18Op_onaq_12hour_05072020.csv")
p_o_pri = p_o_pri[,3:4]

#BL water vapor isotopes - prior
wviso = read.csv("data/wvIso_onaq_12hour_05072020.csv")
wviso$t = as.POSIXct(wviso$t)
bliso = wviso[wviso$level == 1,]
bliso.ts = data.frame("t" = as.POSIXct(ts), "d18O_m" = rep(NA), "d18O_sd" = rep(NA))
bliso.ts[,2:3] = bliso[match(as.POSIXct(ts), bliso$t),4:5]
for(i in seq_along(bliso.ts$t)){
  if(is.na(bliso.ts$d18O_m[i])){
    j = k = i
    while(j >= 1 & is.na(bliso.ts$d18O_m[j])){
      j = j - 1
    }
    while(k <= nrow(bliso.ts) & is.na(bliso.ts$d18O_m[k])){
      k = k + 1
    }
    bliso.ts$d18O_m[i] = mean(bliso.ts$d18O_m[c(j,k)], na.rm = TRUE)
    #uncertainty - might need to revise
    bliso.ts$d18O_sd[i] = sqrt(sum(bliso.ts$d18O_sd[c(j,k)], na.rm = TRUE))
  }
}

#Profile water vapor isotopes & SH - data
##time index
wviso$ti = match(wviso$t, as.POSIXct(ts))
##for now, only keep times where the top level has a measurement
##and the lower leve has a measured sh higher than top
top = max(unique(wviso$level))
for(i in unique(wviso$ti)){
  if(!(top %in% wviso$level[wviso$ti == i])){
    wviso = wviso[wviso$ti != i,]
  }else{
    top.sh = wviso$sh_m[wviso$ti == i & wviso$level == top]
    wviso = wviso[wviso$ti != i | wviso$sh_m >= top.sh,]
  }
}
##remove ti = 1, model doesn't calculate fluxes at that step
wviso = wviso[wviso$ti != 1,]
wviso_top_pri = wviso[wviso$level == top, c("t", "ti", "sh_m",
                                            "sh_sd", "d18O_m",
                                            "d18O_sd")]
##reset the ti for the data to refer to the top_pri
wviso$ti = match(wviso$t, wviso_top_pri$t)
##remove the time from top_pri
wviso_top_pri = wviso_top_pri[, -1]
##finally the data
wviso_data = wviso[wviso$level != top, c("ti", "sh_m",
                                            "sh_sd", "d18O_m",
                                            "d18O_sd")]


#Soil water content - data
swc = read.csv("data/swc_onaq_12hour_exdset3_05072020.csv")
##create time and depth indices
swc$ti = match(swc$t, ts)
swc$di = match(swc$v, c(501, 502, 503, 504, 505))
l = cbind(swc$ti, swc$di, swc$mean_vswc, swc$meas_uncert)


#Soil water isotopes - data
swi = read.csv("data/ONAQsoil.csv")
swi$Date = as.POSIXct(swi$Date, format = "%m/%d/%Y %H:%M")
swi$Date = swi$Date + 3600 * 15
##time index
swi = swi[swi$Date > min(ts) & swi$Date < max(ts),]
swi$ti = match(swi$Date, as.POSIXct(ts))
##depth index
swi$di = rep(NA)
for(i in 1:nrow(swi)){
  swi$di[i] = switch(match(swi$Depth..cm.[i], c(5, 10, 15, 25, 40, 60, 100)),
         1, 2, 2, 3, 4, 5, 5)
}
d_o = data.frame(swi$ti, swi$di, swi$Calibrated.O, swi$Calibrated.O.SE)

#Set up input
dat = list("delta_t" = 0.5, "ET" = etp*1e-3, "pre_pri" = p*1e-3, "st" = st,
           "phi.data" = l, #"phi.calUC.pri" = swc$cal_uncert[1],
           "thick" = c(0.06, 0.1, 0.1, 0.2, 1), "rh" = rh,
           "r_bl_pri" = bliso.ts[,2:3], "wviso.top.pri" = wviso_top_pri,
           "wviso.data" = wviso_data,
           "alpha" = c(1, 2, 2, 2, 1), "p_o_pri" = p_o_pri, "d_o.data" = d_o,
           "nl" = 5, "nt" = length(p))

parms = c("phi", "et", "k_s", "phi_s", "psi_s", "d_o.et", "evap_r",
          "k_dt", "pre", "e_frac", "t_frac", "d_o", "p_o", "diff_o")

#Run inversion
pt = proc.time()
rmod = jags.parallel(model.file = "code/model_generic_oIso.R", 
                     parameters.to.save = parms, 
                     data = dat, inits = NULL, 
                     n.chains = 4, n.iter = 7500, n.burnin = 2500, n.thin = 5)
(proc.time() - pt)[3]

save(rmod, file = "~/rmod3.Rdata")

View(rmod$BUGSoutput$summary)
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

##Compare soil isotope timeseries
plot(apply(sl$d_o[,,1], 2, mean), type = "l", ylim = c(-13, 20))
points(d_o[d_o[,2] == 1,1], d_o[d_o[,2] == 1,3])
for(i in 2:5){
  lines(apply(sl$d_o[,,i], 2, mean), col = i)
  points(d_o[d_o[,2] == i,1], d_o[d_o[,2] == i,3], col = i)
}

