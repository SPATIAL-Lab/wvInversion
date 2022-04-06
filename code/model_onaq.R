model {
  
  #data model 
  for(i in 2:length(phi_1.data)){
    phi_1.data[i] ~ dnorm(phi_1[i], 1/(0.005^2))
    phi_2.data[i] ~ dnorm(phi_2[i], 1/(0.005^2))
    phi_3.data[i] ~ dnorm(phi_3[i], 1/(0.005^2))
  }
  
  #process model
  for(i in 2:length(phi_1.data)){
    phi_1[i] = phi_1[i-1] + phi_1.delta[i] / thick_1
    phi_1.delta[i] = -diff_1_2[i] - adv_1_2[i] + 
      min(SWD_1[i-1], pre[i] - runoff_1[i]) - evap_1[i] - transp_1[i]
    
    phi_2[i] = phi_2[i-1] + phi_2.delta[i] / thick_2
    phi_2.delta[i] = diff_1_2[i] - diff_2_3[i] + adv_1_2[i] -
      adv_2_3[i] + max(0, pre[i] - runoff_1[i] - SWD_1[i-1]) - transp_2[i]
    
    phi_3[i] = phi_3[i-1] + phi_3.delta[i] / thick_3
    phi_3.delta[i] = diff_2_3[i] + adv_2_3[i] - adv_3_4[i]
    
    #transpiration and evaporation
    evap_1[i] = et[i] * e_frac[i]
    transp_1[i] = et[i] * (1 - e_frac[i]) * t_frac_1
    transp_2[i] = et[i] * (1 - e_frac[i]) * (1 - t_frac_1)
    #total ET - the min statement keeps et from drying upper layer below phi = 0.005
    et[i] = min(etx[i], thick_1 * phi_1[i-1] / 
                  (e_frac[i] + (1 - e_frac[i]) * t_frac_1) - 0.005)
    e_frac[i] ~ dunif(0, 1)    
    etx[i] ~ dnorm(ET[i,1], 1/ET[i,2]^2)I(0,)
    
    #diffusion
    diff_1_2[i] = mean(c(adv_1_2[i], adv_2_3[i])) * 
      (-2 * psi_s * phi_s ^ 2 / phi_1[i-1] ^ 3) * 
      (phi_1[i-1] - phi_2[i-1]) / ((thick_1 + thick_2) / 2)
    diff_2_3[i] = mean(c(adv_2_3[i], adv_3_4[i])) * 
      (-2 * psi_s * phi_s ^ 2 / phi_2[i-1] ^ 3) * 
      (phi_2[i-1] - phi_3[i-1]) / ((thick_2 + thick_3) / 2)
    
    #advection
    adv_1_2[i] = (phi_1[i-1] - phi_drain_1[i]) * thick_1
    adv_2_3[i] = (phi_2[i-1] - phi_drain_2[i]) * thick_2
    adv_3_4[i] = (phi_3[i-1] - phi_drain_3[i]) * thick_3
    
    phi_drain_1[i] = ((a / thick_1 * phi_s ^ (-c) * delta_t + z_1[i]) *
                          (1 - c)) ^ (1 / (1 - c))
    phi_drain_2[i] = ((a / thick_2 * phi_s ^ (-c) * delta_t + z_2[i]) *
                        (1 - c)) ^ (1 / (1 - c))
    phi_drain_3[i] = ((a / thick_3 * phi_s ^ (-c) * delta_t + z_3[i]) *
                        (1 - c)) ^ (1 / (1 - c))
    
    z_1[i] = phi_1[i-1] ^ (1 - c) / (1 - c)
    z_2[i] = phi_2[i-1] ^ (1 - c) / (1 - c)
    z_3[i] = phi_3[i-1] ^ (1 - c) / (1 - c)

    #runoff
    runoff_1[i] = pre_x[i] * (pre_x[i] / (pre_x[i] + SWD_eff[i-1]))
    pre_x[i] = max(pre[i] - SWD_1[i-1], 0)
    SWD_1[i-1] = (phi_s - phi_1[i-1]) * thick_1
    SWD_eff[i-1] = SWD_2[i-1] * (1 - exp(-k_dt * delta_t))
    SWD_2[i-1] = (phi_s - phi_2[i-1]) * thick_2
    pre[i] ~ dnorm(pre_1[i], 1/1e-4^2)T(0,)
  }
  
  #initial condition
  phi_1[1] = phi_1.data[1]
  phi_2[1] = phi_2.data[1]
  phi_3[1] = phi_3.data[1]

  #process model priors
  #vertical distribution of transpiration
  #et.tau = 1 / 0.0002 ^ 2

  thick_1 = 0.06
  thick_2 = 0.1
  thick_3 = 2
  
  #advection
  c = 2 * b + 3
  a = -k_s * 8.64e4

#  t_frac_1 = 0.1
#  phi_s = 0.4
#  k_s = 7.1e-6
  b = 8
#  k_dt = 3
#  psi_s = -0.62 

  t_frac_1 ~ dbeta(2, 4)    
  phi_s ~ dunif(0.4, 0.45)
  k_s ~ dunif(1e-7, 1e-5)
#  b ~ dunif(4, 12)
  k_dt ~ dunif(2, 4)
  psi_s ~ dunif(-0.66, -0.58) 
}