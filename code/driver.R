library(R2jags)

pcp = read.csv("data/pcp_onaq_12hour_exdset2_05072020.csv")
swc = read.csv("data/swc_onaq_12hour_exdset2_05072020.csv")
le = read.csv("data/le_nsae_onaq_12hour_exdset2_05072020.csv")

#all timesteps in model
ts = pcp$t

#remove any NAs from pcp
pcp.narm = pcp[!is.na(pcp$bulk_precip),]

#remove NAs from le
le.narm = le[!(is.na(le$mean_nsae_h2o)),]


#get individual soil depth series
swc1 = swc[swc$verticalPosition == 501,]
swc2 = swc[swc$verticalPosition == 502,]
swc3 = swc[swc$verticalPosition == 503,]

#remove NAs from soil series
swc1 = swc1[!(is.na(swc1$mean_vswc)),]
swc2 = swc2[!(is.na(swc2$mean_vswc)),]
swc3 = swc3[!(is.na(swc3$mean_vswc)),]




#time parameters
d = 60 * 24
delta_t = 1/24

#precipitation timeseries
pre_1 = rep(0, d)
pre_1[30:50] = 0.1 * delta_t
pre_1[800:815] = 0.1 * delta_t
pre_1[1000:1045] = 0.04 * delta_t

#ET timeseries
ETt = seq(0, 0.0003, length.out = d)
ETp = rep(mean(ETt), d)
#for(i in 2:d){
#  ET[i] = rnorm(1, ET[1] + (ET[i-1] - ET[1]) * 0.9, ET[1] / 40)  
#}
#plot(ET, type = "l")

#storage for variables
phi_1 = phi_2 = phi_3 = SWD_2 = SWD_1 = SWD_eff = pre_x = numeric(d)
runoff_1 = adv_1_2 = adv_2_3 = adv_3_4 = diff_1_2 = diff_2_3 = numeric(d)
e_frac = etx = et = transp_1 = transp_2 = evap_1 = numeric(d)
phi_1.delta = phi_2.delta = phi_3.delta = numeric(d)

#initial condition
phi_1[1] = phi_2[1] = phi_3[1] = 0.2

#fixed model parameters
t_frac_1 = 0.1
e_frac = rep(0.3, d)
thick_1 = 0.2
thick_2 = 0.8
thick_3 = 2

phi_s = 0.41
k_s = 7.1e-6
b = 8
k_dt = 3
psi_s = -0.62 

#W diffusion
for(i in 2:d){
  SWD_2[i] = (phi_s - phi_2[i-1]) * thick_2
  SWD_eff[i] = SWD_2[i] * (1 - exp(-k_dt * delta_t))
  SWD_1[i] = (phi_s - phi_1[i-1]) * thick_1
  #no negative pre_x
  pre_x[i] = max(pre_1[i] - SWD_1[i], 0)
  runoff_1[i] = pre_x[i] * (pre_x[i] / (pre_x[i] + SWD_eff[i]))
  
  #advection
  adv_1_2[i] = k_s * 8.64e4 * delta_t * (phi_1[i-1] / phi_s) ^ (2 * b + 3)
  adv_2_3[i] = k_s * 8.64e4 * delta_t * (phi_2[i-1] / phi_s) ^ (2 * b + 3)
  adv_3_4[i] = k_s * 8.64e4 * delta_t * (phi_3[i-1] / phi_s) ^ (2 * b + 3)
  
  #diffusion - swap phi terms
  diff_1_2[i] = adv_1_2[i] * (-2 * psi_s * phi_s ^ 2 / phi_1[i-1] ^ 3) * 
    (phi_1[i-1] - phi_2[i-1]) / ((thick_1 + thick_2) / 2)
  diff_2_3[i] = adv_2_3[i] * (-2 * psi_s * phi_s ^ 2 / phi_2[i-1] ^ 3) * 
    (phi_2[i-1] - phi_3[i-1]) / ((thick_2 + thick_3) / 2)
  
  et[i] = ETt[i]
  
  evap_1[i] = et[i] * e_frac[i]
  transp_1[i] = et[i] * (1 - e_frac[i]) * t_frac_1
  transp_2[i] = et[i] * (1 - e_frac[i]) * (1 - t_frac_1)
  
  phi_1.delta[i] = -diff_1_2[i] - adv_1_2[i] + 
    min(SWD_1[i], pre_1[i] - runoff_1[i]) - evap_1[i] - transp_1[i]
  phi_1[i] = phi_1[i-1] + phi_1.delta[i] / thick_1
  
  phi_2.delta[i] = diff_1_2[i] - diff_2_3[i] + adv_1_2[i] -
    adv_2_3[i] + max(0, pre_1[i] - runoff_1[i] - SWD_1[i]) - transp_2[i]
  phi_2[i] = phi_2[i-1] + phi_2.delta[i] / thick_2
  
  phi_3.delta[i] = diff_2_3[i] + adv_2_3[i] - adv_3_4[i]
  phi_3[i] = phi_3[i-1] + phi_3.delta[i] / thick_3
  
}
plot(phi_1, type = "l", ylim = c(0,0.4))
lines(phi_2, lty = 2)

ar = 12
inds = seq(1, d, by = ar)

#Now assign to data
phi_1.data = phi_1[inds]
phi_2.data = phi_2[inds]

dat = list("delta_t" = delta_t * ar, "ET" = ETp[inds] * ar, "pre_1" = pre_1[inds] * ar,
           "phi_1.data" = phi_1.data, "phi_2.data" = phi_2.data)

parms = c("phi_1", "phi_2", "adv_1_2", "adv_2_3", "adv_3_4", "et", 
          "e_frac", "t_frac_1", "k_s", "b")

rmod = jags.parallel(model.file = "code/model.R", 
                     parameters.to.save = parms, 
                     data = dat, inits = NULL, 
                     n.chains = 3, n.iter = 5000, n.burnin = 50, n.thin = 1)

sl = rmod$BUGSoutput$sims.list

#Compare soil phi timeseries
plot(apply(sl$phi_1, 2, mean), type = "l", col = "red")
lines(phi_1.data)
lines(apply(sl$phi_2, 2, mean), lty = 3, col = "red")
lines(phi_2.data, lty = 3)

#Compare ET timeseries
plot(ETt[inds], type = "l", ylim = c(0, 1e-2))
lines(apply(sl$et, 2, mean), col = "red")

