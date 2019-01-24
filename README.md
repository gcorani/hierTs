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
* `iTest`: control of the split the data between train and test. Each different `iTest` yields a different split between training and test. Useful for parallelizing the experiments. Admissible values are between 1 and 45. If unspecified, we resort to the default value iTest=1.


Some examples:
```R
 hierRec(dset="infantgts", fmethod="ets", h=1, iTest=1) #use the first training/test split
 hierRec(dset="infantgts", fmethod="arima", h=1, iTest=2) #use the second training/test split
 hierRec(dset="tourism", fmethod="ets", h=2, iTest=1)
 hierRec(dset="tourism", fmethod="arima", h=1) #use the first training/test split
```
The raw results are written in the file `results/mseHierReconc[dsetName].csv`.


## Reconciliation of temporal hierarchies


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

## Extensive experiments with temporal hierarchies
To reconcile  the whole set of monthly time series from the Mcomp package:
```R
 batchM3(type="monthly",fmethod="ets")
 batchM3(type="quarterly",fmethod="ets")
 batchM3(type="quarterly",fmethod="arima")
 batchM3(type="monthly",fmethod="arima")
```
The raw results are written in the file `results/mseHierReconc[dsetname].csv`.

The Mcomp package has to be installed.

## Analyzing the results
The previous functions save the raw results in a csv file within the `results/` directory (if missing, the directory is created.
If for instance you have reconciled the `infantgts` data set multiple times using different values of iTest, you will find the results within the file `mseHierReconcinfantgts.csv`, containing the mse of each method; each row of the file corresponds to a different iTest.

The results can be analyzed via
```R
parseHierResults_aggregatedH ("name of the data set")
```
For instance:

```R
parseHierResults_aggregatedH ("infantgts")
parseHierResults_aggregatedH ("tourism")
```

This creates another text file (`summarytourism.csv` or `summaryinfantgts.csv`) containing the mse statistics given in the paper (% of times better than other methods, median mse ratio). Moreover it creates the boxplot of the log (relative mse), saved as pdf files.

The same analysis, performed separately for each h can be done as follows:

```R
parseHierResults_eachH ("infantgts")
parseHierResults_eachH ("tourism")
```



