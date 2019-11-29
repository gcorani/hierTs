# Reconciling hierarchical forecasts via Bayes' rule

This package implements  the Bayesian reconciliation of hierarchical forecasts.

In particular the function `hierRec` reconciles forecasts for hierachical / grouped time series.


## Required packages for reproducing the paper experiments

* `forecast` (produces base forecasts with either auto.arima or ets)
* `hts`   (algorithms for reconciling hierarchical time series, including minT)
* `fpp2`  (various time series utilities)


## Choice of data sets and base forecasts
Two real data sets can be used: `infantgts` (available from `hts`) or `tourism` (raw data available from [https://robjhyndman.com/publications/mint/](https://robjhyndman.com/publications/mint/). The csv file of tourism is available in this repository. When the   `tourism` data set is selected, function `loadTourism.R` reads the csv file and builds the hierarchy.

Alternatively, the data can be synthetically generated:
* from a hierarchy having two bottom time series and one top time series  (`synthetic`).
* from a hierarchy having four bottom time series, two intermediate and one top time series  (`syntheticLarge`).

The base forecasts can be created using either `auto.arima` or `ets`, both available from `forecast`.

## Examples

```R
 hierRec(dset="infantgts", fmethod="ets", h=1, iTest=1) 
 hierRec(dset="infantgts", fmethod="arima", h=1, iTest=2) 
 hierRec(dset="tourism", fmethod="ets", h=2, iTest=1)
 hierRec(dset="tourism", fmethod="arima", h=1) 
```


Examples with generated data sets:
```R
 hierRec(dset="synthetic", h=1, synth_n=100, synthCorrel=0.2, howManyBottom=2) 
 hierRec(dset="synthetic", h=3, synth_n=300, synthCorrel=0.8, howManyBottom=4)  
```
Arguments:

* `dset` : can be either `infantgts`, `tourism`, `synthetic` or `syntheticLarge` 

* `fmethod` : method for generating the base forecasts: it can be either `arima` or `ets`. Default: 'ets'.

* `h`: forecast horizon for which we reconcile the forecasts. We use between 1 and 4 in the experiments of the paper.
Default: 1.

* `iTest`: controls how the split the data between train and test. The training set contains the data from the first observation up to the observation in position (length(timeSeries) - h - (iTest - 1)). This is useful for parallelizing the experiments. Admissible values are between 1 and 45. Default: 1.

An additional set of parameters applies to set of experiments with synthetic time series.

* `seed` : seed (default:0)

* `synth_n` : length of the generated time series. Default: 100.

* `synthCorrel` : correlation of the noise affecting the bottom time series. Applies only to the `synthetic` case; default: 0.5.
For the `syntheticLarge` case, the covariance matrix of the noise is set as in the MinT paper (Wickramasuriya et al., 2019)



The reconciliation results  are written in the file `results/mseHierReconc[dsetName].csv`.
The file shows the mean absolute error (mae) of the base forecasts, minT, and the Bayesian reconciliation (with and without estimating correlation). 


## Analyzing the results
The previous functions save the raw results in a csv file within the `results/` directory (if missing, the directory is created.
If for instance you have reconciled the `infantgts` data set multiple times using different values of iTest, you will find the results within the file `mseHierReconcinfantgts.csv`, containing the mse of each method; each row of the file corresponds to a different iTest.

The results can be analyzed via
```R
parseHierResults("name of the data set")
```
For instance:

```R
parseHierResults("infantgts")
parseHierResults("tourism")
```

This creates another text file (`summarytourism.csv` or `summaryinfantgts.csv`) containing the mse statistics given in the paper. Moreover it creates the boxplot of the log (relative mse), saved as pdf files.




