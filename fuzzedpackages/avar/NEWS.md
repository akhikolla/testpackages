# avar 0.1.1

## Features

This package provides the tools necessary to compute the empirical Allan Variance (AV) and use it to estimate the parameters of (latent) time series models. The estimation of the AV is performed through the estimator proposed by Allan (1966) and, based on this quantity, the Allan Variance Linear Regression (AVLR) approach is often used by engineers to retrieve the parameters of time series models which are assumed to underlie the observed signals (see for example Guerrier, Molinari, and Stebler 2016). These estimators are implemented in this package along with the relevant plotting and summary functions.

Compared to the 0.1.0 version, we add the application of Allan variance on IMU data. Specifically, we add the following new features: 
- New function avar.imu() which allows the computation of Allan variance based on IMU data;
- New function plot.imu_avar() which allows to plot the Allan variance computed based on IMU data;
- New function avlr.imu_avar() which allows to compute the Allan Variance Linear Regression estimator based on IMU data;
- A few datasets of Allan variance based on real IMU data as examples to test the new functions.
