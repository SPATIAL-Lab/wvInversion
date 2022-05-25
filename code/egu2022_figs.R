#First load all data from driver_onaq_Oiso.R

#Resources
library(R2jags)
source("code/helpers.R")

#Results
load("C:/Users/u0133977/Dropbox/Hypomirror/NEON_wvInversion/rmod.Rdata")
sl = rmod$BUGSoutput$sims.list

#timeseries dates as dates
tsw = as.POSIXct(ts)

#reload wviso in more useful way
wviso = read.csv("data/wvIso_onaq_12hour_05072020.csv")
wviso$t = as.POSIXct(wviso$t)
wviso$ti = match(wviso$t, as.POSIXct(ts))

#Data plots
##1
png("C:/Users/u0133977/Dropbox/Hypomirror/Talks/2022/EGU/data.png",
    units = "in", width = 7, height = 7, res = 300)
layout(matrix(c(1, 2), nrow = 2), heights = lcm(c(3.1, 3.9) * 2.54))
par("mai"= c(0.2,1,0.2,0.2)) 
plot(tsw[l[l[,2] == 1,1]], l[l[,2] == 1,3], type = "l", lty = 3, 
     lwd = 2, xlim = c(1.5885e9, 1.59653e9), ylim = c(0.05, 0.3), 
     xlab = "", ylab = "VSWC (unitless)", axes = FALSE)
axis(1, labels = FALSE, at = as.POSIXct(c("2020-05-01", "2020-06-01",
                                          "2020-07-01", "2020-08-01")))
axis(2)
box()
for(i in 2:5){
  lines(tsw[l[l[,2] == i,1]], l[l[,2] == i,3], 
        lty = 3, lwd = 2, col = i)
}
legend("topright",
       c("6 cm", "16 cm", "26 cm", "46 cm", "106 cm"),
       bty = "n", ncol = 2,
       lty = 3, col = c(1:5))

par("mai"= c(1,1,0.2,0.2)) 
plot(tsw[d_o[d_o[,2] == 1,1]], d_o[d_o[,2] == 1,3], col = 1, 
     xlim = c(1.5885e9, 1.59653e9), ylim = c(-30, 12),
     xlab = "", ylab = expression("Soil or vapor  "*delta^{18}*"O"))
for(i in 2:5){
  points(tsw[d_o[d_o[,2] == i,1]], d_o[d_o[,2] == i,3], 
         col = i)
}
lines(tsw[wviso$ti[wviso$level == 1]], wviso$d18O_m[wviso$level ==1], 
     type = "l", lwd = 2, lty = 3)
for(i in 2:4){
  lines(tsw[wviso$ti[wviso$level == i]], 
        wviso$d18O_m[wviso$level == i], col = i, lwd = 2, lty = 3)
}
dev.off()

##2
png("C:/Users/u0133977/Dropbox/Hypomirror/Talks/2022/EGU/vapZoom.png",
    units = "in", width = 7, height = 6, res = 300)
par("mai"= c(1,1,0.2,0.2)) 

plot(tsw[wviso$ti[wviso$level == 1]], wviso$d18O_m[wviso$level ==1], 
     type = "l", lwd = 2, lty = 3, xlab = "", ylab = "",
     xlim = as.POSIXct(c("2020-06-06", "2020-06-15")),
     ylim = c(-29, -16), axes = FALSE)
for(i in 2:4){
  lines(tsw[wviso$ti[wviso$level == i]], 
        wviso$d18O_m[wviso$level == i], col = i, lwd = 2, lty = 3)
}
axis(1, at = as.POSIXct(c("2020-06-06", "2020-06-09", 
                          "2020-06-12", "2020-06-15")),
     labels = c("June 6", "June 9", "June 12", "June 15"))
axis(2)
box()
legend("topright", c("top", "", "", "base"), col = c(1:4),
       lty = 3, bty = "n", inset = 0.05)
dev.off()

#Result plots
##1
png("C:/Users/u0133977/Dropbox/Hypomirror/Talks/2022/EGU/results1.png",
    units = "in", width = 7, height = 7, res = 300)
layout(matrix(c(1, 2), nrow = 2), heights = lcm(c(3.1, 3.9) * 2.54))
par("mai"= c(0.2,1,0.2,0.2)) 
plot(tsw[l[l[,2] == 1,1]], l[l[,2] == 1,3], type = "l", lty = 3, 
     xlim = c(1.5885e9, 1.59653e9), ylim = c(0.05, 0.3), 
     xlab = "", ylab = "VSWC (unitless)", axes = FALSE)
axis(1, labels = FALSE, at = as.POSIXct(c("2020-05-01", "2020-06-01",
                                          "2020-07-01", "2020-08-01")))
axis(2)
box()
ed = apply(sl$phi[,,1], 2, quantile, probs = c(0.5, 0.25, 0.5, 0.75, 0.95))
ed = t(ed)
tsdens(cbind(tsw, ed))

for(i in 2:5){
  lines(tsw[l[l[,2] == i,1]], l[l[,2] == i,3], lty = 3, col = 3)
  ed = apply(sl$phi[,,i], 2, quantile, probs = c(0.5, 0.25, 0.5, 0.75, 0.95))
  ed = t(ed)
  tsdens(cbind(tsw, ed), base = i)
}
for(i in 2:5){
  lines(tsw[l[l[,2] == i,1]], l[l[,2] == i,3], 
        lty = 3, lwd = 2, col = i)
}
legend("topright",
       c("6 cm", "16 cm", "26 cm", "46 cm", "106 cm"),
       bty = "n", ncol = 2,
       lty = 3, col = c(1:5))

par("mai"= c(1,1,0.2,0.2)) 
plot(tsw[d_o[d_o[,2] == 1,1]], d_o[d_o[,2] == 1,3], col = 1, 
     xlim = c(1.5885e9, 1.59653e9), ylim = c(-15, 21),
     xlab = "", ylab = expression("Soil or vapor  "*delta^{18}*"O"))
ed = apply(sl$d_o[,,1], 2, quantile, probs = c(0.5, 0.25, 0.5, 0.75, 0.95))
ed = t(ed)
tsdens(cbind(tsw, ed))
for(i in 2:5){
  points(tsw[d_o[d_o[,2] == i,1]], d_o[d_o[,2] == i,3], 
         col = i)
  ed = apply(sl$d_o[,,i], 2, quantile, probs = c(0.5, 0.25, 0.5, 0.75, 0.95))
  ed = t(ed)
  tsdens(cbind(tsw, ed), base = i)
}
dev.off()

##2
ed = apply(sl$e_frac, 2, quantile, probs = c(0.5, 0.25, 0.5, 0.75, 0.95))
ed = t(ed)

png("C:/Users/u0133977/Dropbox/Hypomirror/Talks/2022/EGU/results2.png",
    units = "in", width = 7, height = 5, res = 300)
layout(matrix(c(1, 2), nrow = 2), heights = lcm(c(1.1, 3.9) * 2.54))
par("mai"= c(0,1,0,0.2)) 
plot(tsw[l[l[,2] == 1,1]], l[l[,2] == 1,3], type = "l", lty = 3, 
     xlim = c(1.5885e9, 1.59653e9), ylim = c(0.05, 0.3), 
     xlab = "", ylab = "", axes = FALSE)

par("mai"= c(1,1,0.2,0.2))
plot(tsw[-1], ed[,1], ylim = range(ed), type = "n",
     xlim = c(1.5885e9, 1.59653e9),
     xlab = "Date", ylab = "E / (E + T)")
tsdens(cbind(tsw[-1], ed), "blue")
dev.off()


lines(tsw[wviso$ti[wviso$level == 1]], wviso$d18O_m[wviso$level ==1], 
      type = "l", lwd = 2, lty = 3)
for(i in 2:4){
  lines(tsw[wviso$ti[wviso$level == i]], 
        wviso$d18O_m[wviso$level == i], col = i, lwd = 2, lty = 3)
}
dev.off()

plot(apply(sl$d_o[,,1], 2, mean), type = "l", ylim = c(-13, 20))
points(d_o[d_o[,2] == 1,1], d_o[d_o[,2] == 1,3])
for(i in 2:5){
  lines(apply(sl$d_o[,,i], 2, mean), col = i)
  points(d_o[d_o[,2] == i,1], d_o[d_o[,2] == i,3], col = i)
}


plot(apply(sl$phi[,,1], 2, mean), type = "l", ylim = c(0.05, 0.3))
lines(l[l[,2] == 1,1], l[l[,2] == 1,3], lty = 3)
for(i in 2:5){
  lines(apply(sl$phi[,,i], 2, mean), col = i)
  lines(l[l[,2] == i,1], l[l[,2] == i,3], lty = 3, col = i)
}