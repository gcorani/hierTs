# Bayesian reconcialition of grouped time series and temporal hierarchies 

This package implements all the functions necessary for Bayesian reconciliation of hierarchical forecasts.

In particular:

1. the function `hierRec` reconciles forecasts for hierachical / grouped time series.

2. The function `temporalRec` reconciles a temporal hierarchy (i.e., forecast for the same time series sampled at different frequencies).

## Required packages for reproducing the paper experiments

* `forecast` (produces base forecasts with either auto.arima or ets)
* `hts`   (algorithms for reconciling hierarchical time series, including minT)
* `thief` (algorithms for reconciling temporal hierarchies)
* `huge`  (glasso estimation of the covariance matrix)
* `mcomp` (time series of the M3 competition)
* `fpp2`  (various time series utilities)


## Running reconciliation of grouped time series
Two data sets can be used: `infantgts` (available from `hts`) or `tourism` (raw data available from [https://robjhyndman.com/publications/mint/](https://robjhyndman.com/publications/mint/); the csv file is available in this repository. When  `tourism` is passed as a parameter, another function of this directory imports the raw data into R.)

The base forecasts can be created using either `auto.arima` or `ets`, both available from `forecast`.

Arguments:

* `dset` : can be either `infantgts` or `tourism`

* `fmethod` : method for generating the base forecasts: it can be either `arima` or `ets`. Default: 'ets'.

* `h`: forecast horizon for which we reconcile the forecasts. We use between 1 and 4 in the experiments of the paper.
The reconciliation is performed independently for different values of h. Default: 1.

* `iTest`: control of the split the data between train and test. The train contains the data from the first observation up to the observation in position (length(timeSeries) - h - (iTest - 1)). This is  especially useful for parallelizing the experiments. Admissible values are between 1 and 45. Default: 1.

An additional set of parameters should be specified for running experiments with synthetic time series.

* `seed` : seed (default:0)

* `synthCorrel` : correlation of the noise affecting the bottom time series. Default: 0.5.

* `synth_n` : how many data point each generated time series should contain. Default: 100.

* `howManyBottom` : how many bottom time series in the hierarchy (either 2 or 4).


Examples with real data sets:

```R
 hierRec(dset="infantgts", fmethod="ets", h=1, iTest=1) 
 hierRec(dset="infantgts", fmethod="arima", h=1, iTest=2) 
 hierRec(dset="tourism", fmethod="ets", h=2, iTest=1)
 hierRec(dset="tourism", fmethod="arima", h=1) 
```

The raw results are written in the file `results/mseHierReconc[dsetName].csv`.

Examples with generated data sets:
```R
 hierRec(dset="synthetic", h=1, synth_n=100, synthCorrel=0.2, howManyBottom=2) 
 hierRec(dset="synthetic", h=3, synth_n=300, synthCorrel=0.8, howManyBottom=4)  
```



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



