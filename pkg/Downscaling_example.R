require(minpack.lm)
source("Downscaling.Functions.R")

## List of models to be tried
model.list<-c("Nachman","PL","Logis","Poisson","NB","GNB","INB")
names(all.pars)<-model.list

## Input data - column 1 = cell area; column 2 = % occupancy
## Some dummy data (cell area in km^2)
observed<-data.frame(Cell.area=c(10,40,160,640),Occ=c(0.6,0.85,0.95,1))

## Total area (in km^2)
A0<-640000 

## Convert cell areas to proportion of total area (should standardise any data to fit within starting parameter range)
observed[,"Cell.area"]<-observed[,"Cell.area"]/A0

## Convert occupancy = 1 to NAs
observed[observed[,"Occ"]==1,2]<-NA 

## Determine the number of observed scales and predicted scales 
## (here automatically generated but could be inputed as a vector of cell sizes)
obs.scales<-length(observed[,"Cell.area"])
est.scales<-c((1:9/A0),observed[,"Cell.area"])

## Dataframe to store predicted results
estimated<-matrix(NA,ncol=length(model.list)+1,nrow=length(est.scales))
colnames(estimated)<-c("Cell.area",model.list)
estimated[,"Cell.area"]<-est.scales

## Loop to run and plot all downscaling functions
par(mfrow=c(3,3),mar=c(5,4,2,1))
for(i in 1:7){
  model.run<-model.list[i]
  resid.fun<-getFunction(paste("resid",model.run,sep=""))
  pred.fun<-getFunction(model.run)  
  parStart<-get(paste("par",model.run,sep=""))

  pars<-nls.lm(par=parStart,fn=resid.fun,
               A=observed[!is.na(observed[,2]),"Cell.area"],
               observed=log(observed[!is.na(observed[,"Occ"]),"Occ"]),
               control=nls.lm.control(maxiter=1000))
  estimated[,i+1]<-exp(pred.fun(A=estimated[,"Cell.area"],as.list(coef(pars))))
  
  plot(observed[,"Occ"]~observed[,"Cell.area"],type="n",log="xy",
       xlim=c(min(est.scales),max(est.scales)),ylim=c(0.001,1),
       xlab="Log cell area",ylab="Log occupancy",main=model.run)
  points(estimated[,i+1]~estimated[,"Cell.area"],type="l",col="red",lwd=2)
  points(observed[,"Occ"]~observed[,"Cell.area"],type="b",lwd=2)
}
estimated

