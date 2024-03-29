model {
  
  # data model 
  for(i in 1:length(phi.data[,1])){
    phi.data[i, 3] ~ dnorm(phi[phi.data[i, 1], phi.data[i, 2]], 1/(0.005^2))
  }
  
  # process model
  for(i in 2:nt){
    for(j in 1:nl){
      phi[i, j] = phi[i-1, j] + phi.delta[i, j] / thick[j]
    }

    # PDEs
    phi.delta[i, 1] = -diff[i, 1] - adv[i, 1] + infil[i, 1] -
      evap[i] - transp[i, 1]
    
    for(j in 2:(nl-1)){
      phi.delta[i, j] = diff[i, j-1] - diff[i, j] + adv[i, j-1] - adv[i, j] +
        infil[i, j] - transp[i, j]
    }
    
    phi.delta[i, nl] = diff[i, nl-1] + adv[i, nl-1] - adv[i, nl] - transp[i, nl]
    
    # transpiration
    for(j in 1:nl){
      transp[i, j] = et[i] * (1 - e_frac[i]) * t_frac[j]
    }

    # evaporation
    evap[i] = et[i] * e_frac[i]
    
    # total ET - the min statement keeps et from drying upper layer below phi = 0.005
    et[i] = min(etx[i], thick[1] * phi[i-1, 1] / 
                  (e_frac[i] + (1 - e_frac[i]) * t_frac[1]) - 0.005)
    e_frac[i] ~ dunif(0, 1)    
    etx[i] ~ dnorm(ET[i,1], 1/ET[i,2]^2)I(0,)
    
    # diffusion
    for(j in 1:(nl-1)){
      diff[i, j] = mean(c(adv[i, j], adv[i, j+1])) * 
        (-2 * psi_s * phi_s ^ 2 / phi[i-1, j] ^ 3) * 
        (phi[i-1, j] - phi[i-1, j+1]) / ((thick[j] + thick[j+1]) / 2)
    }

    # advection    
    for(j in 1:nl){
      adv[i, j] = (phi[i-1, j] - phi_drain[i, j]) * thick[j]
      phi_drain[i, j] = ((a / thick[j] * phi_s ^ (-c) * delta_t + z[i, j]) *
                          (1 - c)) ^ (1 / (1 - c))
      z[i, j] = phi[i-1, j] ^ (1 - c) / (1 - c)
    }

    # runoff
    for(j in 3:(nl-1)){
      infil[i, j] = 0
    }
    infil[i, 2] = pre[i] - runoff[i] - infil[i, 1]
    infil[i, 1] = min(SWD_1[i-1], pre[i] - runoff[i])
    runoff[i] = pre_x[i] * (pre_x[i] / (pre_x[i] + SWD_eff[i-1]))
    pre_x[i] = max(pre[i] - SWD_1[i-1], 0)
    SWD_1[i-1] = (phi_s - phi[i-1, 1]) * thick[1]
    SWD_eff[i-1] = SWD_2[i-1] * (1 - exp(-k_dt * delta_t))
    SWD_2[i-1] = (phi_s - phi[i-1, 2]) * thick[2]
    pre[i] ~ dnorm(pre_pri[i], 1/1e-4^2)T(0,)
  }
  
  # initial condition
  for(j in 1:nl){
    phi[1, j] ~ dunif(0.15, 0.3)
  }

  # derived parms
  c = 2 * b + 3
  a = -k_s * 8.64e4

  # parameter priors
  t_frac ~ ddirch(alpha)  
  #t_frac_1 = 0.1

  phi_s ~ dunif(0.4, 0.45)
  #phi_s = 0.4

  k_s ~ dgamma(10, 1e5)
  #k_s = 7.1e-6
  
  #b ~ dunif(4, 12)
  b = 8
  
  k_dt ~ dunif(2, 4)
  #k_dt = 3

  psi_s ~ dunif(-0.66, -0.58) 
  #psi_s = -0.62 
}