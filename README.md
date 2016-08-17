# expmixdist
Matlab functions to fit parameters for a mixture of exponential distribution

## Introduction

Matlab has some good functions for dealing with exponential mixture distributions, e.g. `expfit`.  These functions extend this support to fit distributions that are mixtures of exponential components.

## Getting started

Run the [getting started script](getting_started.m) for an example of usage.

## Fitting an exponential mixture distribution

```Matlab
[ taushat, weightshat ] = emdistfit(wait_times, tau_guess, weight_guess);
```

## Plotting a histogram and the fitted pdf

```Matlab
emhist(wait_times, taushat, weightshat)
```
