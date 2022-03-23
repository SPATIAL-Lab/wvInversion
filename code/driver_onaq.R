library(R2jags)

#precip
p = read.csv("data/pcp_onaq_12hour_exdset2_05072020.csv")
sum(is.na(p$bulk_precip))
ts = p$t
p = p$bulk_precip

et = read.csv("data/le_nsae_onaq_12hour_exdset2_05072020.csv")
#complete timeseries
etdf = data.frame(ts, "et" = rep(NA), "etsd" = rep(NA))
#plug in observed values
etdf$et[etdf$ts %in% et$t] = et$mean_nsae_h2o * 2.257e-6 * 3600
#made up observed sd
etdf$etsd[etdf$ts %in% et$t] = 0.05
#fill missing values with means and sds of day and night obs
dmean = mean(etdf$et[seq(2, nrow(etdf), by = 2)], na.rm = TRUE)
dsd = sd(etdf$et[seq(2, nrow(etdf), by = 2)], na.rm = TRUE)
nmean = mean(etdf$et[seq(1, nrow(etdf), by = 2)], na.rm = TRUE)
nsd = sd(etdf$et[seq(1, nrow(etdf), by = 2)], na.rm = TRUE)
etmean = rep(c(nmean, dmean), nrow(etdf) / 2)
etsd = rep(c(nsd, dsd), nrow(etdf) / 2)
etdf$et[is.na(etdf$et)] = etmean[is.na(etdf$et)]
etdf$etsd[is.na(etdf$etsd)] = etsd[is.na(etdf$etsd)]
#make it a matrix
etp = as.matrix(etdf[,2:3])


swc = read.csv("data/swc_onaq_12hour_exdset2_05072020.csv")
#pull depths
l1 = swc$mean_vswc[swc$verticalPosition == 501]
l2 = swc$mean_vswc[swc$verticalPosition == 502]
l3 = swc$mean_vswc[swc$verticalPosition == 503]
#something's missing from depth 3, its time step #66
!swc$t[swc$verticalPosition == 501] %in% swc$t[swc$verticalPosition == 503]
#just replicate ts 65
l3 = c(l3[1:65], l3[65], l3[66:175])

#set up input
dat = list("delta_t" = 0.5, "ET" = etp*1e-3, "pre_1" = p*1e-3,
           "phi_1.data" = l1, "phi_2.data" = l2, "phi_3.data" = l3)

parms = c("phi_1", "phi_2", "adv_1_2", "adv_2_3", "adv_3_4", 
          "diff_1_2", "diff_2_3", "diff_3_4", "et", 
          "e_frac", "t_frac_1", "k_s", "b")

rmod = jags.parallel(model.file = "code/model_onaq.R", 
                     parameters.to.save = parms, 
                     data = dat, inits = NULL, 
                     n.chains = 3, n.iter = 50, n.burnin = 10, n.thin = 1)

sl = rmod$BUGSoutput$sims.list

#Compare soil phi timeseries
plot(apply(sl$phi_1, 2, mean), type = "l", col = "red")
lines(phi_1.data)
lines(apply(sl$phi_2, 2, mean), lty = 3, col = "red")
lines(phi_2.data, lty = 3)

#Compare ET timeseries
plot(etp[,1], type = "l")
lines(apply(sl$et, 2, mean), col = "red")

