## Total area (in km^2)
A0<-219000
obset<-data.frame("Cell.area" = c(40000,10000,2500,400,100,25,4,1,0.25),
                 "Occ" = c(0.36529680371,0.18264840177,0.05707762554,0.02009132420,0.00821917808,0.00229452055,0.00047488584,0.00014155251,0.00004423516))
#obset[,"Cell.area"] <- obset[,"Cell.area"] ^ 2 / A0

model.list<-c("Nachman","PL","Logis","Poisson","NB","GNB","INB","FNB","Thomas")
observed<-obset[1:4,]
estimated<-obset[,1]

## Loop to run and plot all downscaling functions
par(mfrow=c(3,3),mar=c(5,4,2,1))
for(i in 1:9){
  model.run<-model.list[i]

  mod<-downscale(occupancy = observed[,2],
                 area = observed[,1],
                 model = model.run,
                 extent = A0,
                 tolerance = 1e-10)
  est<-predict.downscale(object = mod,
                         newdata = estimated,
                         extent = A0,
                         tolerance = 1e-10)
  plot.downscale(est)
  
  plot(obset[,"Occ"]~obset[,"Cell.area"],type="b",log="xy",col="grey",
       xlim=c(min(estimated),max(observed)),ylim=c(0.0001,1),
       xlab="Log cell area",ylab="Log occupancy",main=model.run)
  points(est$predicted[,"Occupancy"]~est$predicted[,"Cell.area"],type="l",col="red",lwd=2)
  points(observed[,"Occ"]~observed[,"Cell.area"],type="b",lwd=2)
  #print(i)
}

model.list<-c("Nachman","PL","Logis","Poisson","NB","GNB","INB","FNB","Thomas")

all<-ensemble.downscale(occupancy = observed[,2],
                   area = observed[,1],
                   newdata = estimated,
                   extent = A0,
                   tolerance = 1e-10,
                   models = "all")
par(mfrow=c(3,3))
for(i in 1:9){
  model.run<-model.list[i]

  plot(obset[,"Occ"]~obset[,"Cell.area"],type="b",log="x",col="grey",
       xlim=c(min(all[,1]),100),ylim=c(0.00001,0.01),
       xlab="Log cell area",ylab="Log occupancy",main=model.run)
  points(all[,"Means"]~all[,"Cell.area"],type="b",lwd=2)
  points(all[,i+1]~all[,"Cell.area"],type="l",col="red",lwd=2)
}
