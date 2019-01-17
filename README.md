# hierTs
#Bayesian reconcialition of hierarchies of time series 

1. The function hierRec reconciles a hierarchy containing forecast for grouped time series.
2. The function temporalRec reconciles a hierarchy containing forecast for the same time series sampled at different frequencies.

Example:
```R
 hierRec(dset=infantgts, fmethod="ets", h=1)
```
