library(imudata)
library(simts)
library(av)

data("cont.imu1")

class(cont.imu1)

attributes(cont.imu1)$freq

test = avar(cont.imu1)

length(test$clusters)
x = avlr(test, wn = 1:12, rw = 12:17)


