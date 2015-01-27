#############################################################################
#               Downscaling.Functions.R 27/11/2014                          #
#                                                                           #
# Functions for downscaling methods                                         #
# residFunction = calculation of residuals during optimisation              #
# parFunction = default starting parameters                                 #
#############################################################################

### Nachman
Nachman<-function(par,A){log(1-exp(-par$C*A^par$z))}
residNachman<-function(par,observed,A){observed-Nachman(par,A)}
parNachman<-data.frame("C"=0.01,"z"=0.01)

### Power Law
PL<-function(par,A){log(par$C*A^par$z)}
residPL<-function(par,observed,A){(observed-PL(par,A))}
parPL<-data.frame("C"=0.01,"z"=0.01)

### Logistic
Logis<-function(par,A){log((par$C*(A^par$z))/(1+(par$C*(A^par$z))))}
residLogis<-function(par,observed,A){observed-Logis(par,A)}
parLogis<-data.frame("C"=0.01,"z"=0.01)

### Poisson
Poisson<-function(par,A){log(1-(exp(-par$lambda*A)))}
residPoisson<-function(par,observed,A){observed-Poisson(par,A)}
parPoisson<-data.frame("lambda"=0.001)

### Negative binomial
NB<-function(par,A){log(1-(1+(par$C*A)/par$k)^-par$k)}
residNB<-function(par,observed,A){observed-NB(par,A)}
parNB<-data.frame("C"=0.01,"k"=0.01)

### Generalised negative binomial
GNB<-function(par,A){log(1-(1+(par$C*A^par$z)/par$k)^-par$k)}
residGNB<-function(par,observed,A){observed-GNB(par,A)}
parGNB<-data.frame("C"=0.1,"z"=0.001,"k"=1)

### Improved negative binomial
INB<-function(par,A){log(1-((par$C*A^(par$b-1))^((par$r*A)/(1-par$C*A^(par$b-1)))))}
residINB<-function(par,observed,A){observed-INB(par,A)}
parINB<-data.frame("C"=0.01,"r"=1e-5,"b"=1)

### Finite negative binomial
FNB<- function(par,A,A0=219000){
  require(Rmpfr)
  as.numeric(log(1-(
    gamma(as(par$W+((A0*par$k)/A)-par$k,"mpfr")) *
      gamma(as((A0*par$k)/A,"mpfr")) /
      (gamma(as(par$W+((A0*par$k)/A),"mpfr")) *
         gamma(as(((A0*par$k)/A)-par$k,"mpfr")))
  )))
}
residFNB<-function(par,observed,A,A0){observed-FNB(par,A,A0)}     # W = 0.01, k = 0.01
parFNB<-data.frame("W"=10,"k"=10)
# need lower=c("W"=-Inf,"k"=0), upper=c("W"=Inf,"k"=min(obset[,1],na.rm=TRUE)*100)
