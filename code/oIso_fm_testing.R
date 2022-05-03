thick = c(0.06, 0.1, 0.1, 0.2, 1)

pv = rmod$BUGSoutput$mean

d_o = r_o = phi_o = phi_o.delta = diff_o = matrix(nrow = 175, ncol = 5)
p_o.del = p_o = evap_a = h = r_a = evap_r = numeric(175)

d_o[1,] = rep(-11)
r_o[1,] = (d_o[1,] / 1000 + 1) * 0.0020052 
phi_o[1,] = pv$phi[1,] * r_o[1,]

for(i in 2:175){
  # precip
  p_o.del[i] = rnorm(1, p_o_pri[i], 1)  
  p_o[i] = (p_o.del[i] / 1000 + 1) * 0.0020052
  
  # evaporation
  #evap_a is a function of temperature, for now fixed
  evap_a[i] = 1 / 1.010
  h[i] = 0.20
  r_a[i] = 0.001955
  evap_r[i] = (evap_a[i] * r_o[i-1, 1] - h[i] * r_a[i]) / 
    ((1 - h[i]) * 0.9723)
  
  # diffusion
  for(j in 1:(nl-1)){
    diff_o[i, j] = mean(c(pv$adv[i, j], pv$adv[i, j+1])) * 
      (-2 * pv$psi_s * pv$phi_s ^ 2 / pv$phi[i-1, j] ^ 3) * 0.9723 *
      (phi_o[i-1, j] - phi_o[i-1, j+1]) / ((thick[j] + thick[j+1]) / 2)
  }
  
  phi_o.delta[i, 1] = pv$pre[i] * p_o[i] - pv$et[i] * pv$e_frac[i] * evap_r[i] - 
    (pv$adv[i, 1] + pv$et[i] * (1 - pv$e_frac[i]) * pv$t_frac[1]) * r_o[i-1, 1] - 
    diff_o[i, 1]
  
  for(j in 2:(nl-1)){
    phi_o.delta[i, j] = pv$adv[i, j-1] * r_o[i-1, j-1] -
      (pv$adv[i, j] + pv$et[i] * (1 - pv$e_frac[i]) * pv$t_frac[j]) * r_o[i-1, j] +
      diff_o[i, j-1] - diff_o[i, j] 
  }
  
  phi_o.delta[i, nl] = pv$adv[i, nl-1] * r_o[i-1, nl-1] - 
    (pv$adv[i, nl] + pv$et[i] * (1 - pv$e_frac[i]) * pv$t_frac[nl]) * r_o[i-1, nl] + 
    diff_o[i, nl-1] 

  for(j in 1:nl){
    phi_o[i, j] = phi_o[i-1, j] + phi_o.delta[i, j] / thick[j]
    r_o[i, j] = phi_o[i, j] / pv$phi[i, j]
    d_o[i, j] = (r_o[i, j] / 0.0020052 - 1) * 1000
  }
}
