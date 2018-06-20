library(imudata)
library(simts)
library(av)

# Real data
data("cont.imu1")
class(cont.imu1)
attributes(cont.imu1)$freq
test = avar(cont.imu1)
print(test)
summary(test)
plot(test)

length(test$clusters)
x = avlr(test, wn = 1:12, rw = 12:17)

# Simulated data
N = 100000
ts = gen_gts(N, WN(sigma2 = 2) + RW(gamma2 = 1))
av_mat_mo = avar(ts, type = "mo", freq = 100)
av_mat_tau = avar(ts, type = "to")
print(av_mat_mo)
summary(av_mat_mo)
plot(av_mat_mo)

fit = avlr(av, wn = 1:8, rw = 10:15)
plot(fit, decomp = TRUE)
fit
fit = avlr(av, wn = 1:8, rw = 10:15, ci = TRUE)
