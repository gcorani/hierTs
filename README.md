# Bayesian reconcialition of grouped time series and temporal hierarchies 

The main functionalities are provided by two functions:
1. The function `hierRec` reconciles forecasts for grouped time series.
2. The function `temporalRec` reconciles a temporal hierarchy (i.e., forecast for the same time series sampled at different frequencies).

## Required packages (available from CRAN):
* `forecast` (produce base forecasts)
* `hts` (reconciliation algorithms)
* `thief` (temporal hierarchies)
* `huge` (glasso estimation of the covariance matrix)
* `mcomp` (time series of the M3 competition)
* `fpp2` (various time series utilities)


## Example of reconciliation of grouped time series
Two data sets can be used: `infantgts` (available from `hts`) or `tourism` (raw data available from [https://robjhyndman.com/publications/mint/](https://robjhyndman.com/publications/mint/); the csv file is available in this repository. When  `tourism` is passed as a parameter, the raw data are parsed by a function.)

The base forecasts can be created using either `arima` or `ets`, both available from `forecast`.

Arguments:
* `dset` : can be either `infantgts` or `tourism`
* `fmethod` : can be either `arima` or `ets`
* `h`: forecast horizon. We use between 1 and 4 in the experiments of the paper.
* `iTest`: control of the split the data between train and test. Each different `iTest` yields a different split between training and test. Useful for parallelizing the experiments. Admissible values are between 1 and 45.


Some examples:
```R
 hierRec(dset="infantgts", fmethod="ets", h=1, iTest=1)
 hierRec(dset="infantgts", fmethod="arima", h=1, iTest=2)
 hierRec(dset="tourism", fmethod="ets", h=2, iTest=1)
 hierRec(dset="tourism", fmethod="arima", h=1)
```
The raw results are written in the file `results/mseHierReconc[dsetName].csv`.

Moreover 


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
