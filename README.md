# Bayesian reconcialition of grouped time series and temporal hierarchies 

The main functionalities are provided by two functions:
1. The function `hierRec` reconciles forecasts for grouped time series.
2. The function `temporalRec` reconciles a temporal hierarchy (i.e., forecast for the same time series sampled at different frequencies).

## Required packages (available from CRAN):
*`forecast` (produce base forecasts)
*`hts` (reconciliation algorithms)
*`thief` (temporal hierarchies)
*`huge` (glasso estimation of the covariance matrix)
*`mcomp` (time series of the M3 competition)
*`fpp2` (various time series utilities)


Example of reconciliation of the infantgts data set, using ets and arima:
```R
 hierRec(dset=infantgts, fmethod="ets", h=1)
 hierRec(dset=infantgts, fmethod="arima", h=1)
```

Example of reconciliation of a monthly temporal hierarchies:

```R
 temporalRec(tsObj, fmethod="ets", periodType="monthly")
```
where tsObj is a list organized as follows:
```R
tsObj$x #training data
tsObj$xx #test data
tsObj$sn #name
```   
The fmethod parameters can be set to either "ets" or "arima".
The periodType parameters can be set to either "monthly", "quarterly" or "weekly".

Experiments on the whole set of monthly time series from the Mcomp package can be performed:
```R
 batchM3 <- function(type="monthly",fmethod="ets"){
```
The Mcomp package has to be installed.
