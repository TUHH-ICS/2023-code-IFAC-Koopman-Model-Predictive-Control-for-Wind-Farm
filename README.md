# 2023-code-IFAC-Koopman Model Predictive Control for Wind Farm Yield Optimization with Combined Thrust and Yaw Control

## General

This repository contains an implementation of the simulations described in 

> A. Dittmer, B. Sharan and H. Werner, "Koopman Model Predictive Control for Wind Farm Yield Optimization with Combined Thrust and Yaw Control"

submitted to IFAC World Congress, 2023.

It may be used to recreate the simulation results and figures from the paper. To do so, run the script `mainGeneratePlots.m`.

Running the simulation for the first time takes roughly 20 minutes using a laptop with an 11th Gen Intel(R) Core(TM) i7-11850H @ 2.50GHz processor

## Simulation Environment

The simulations WindFarmSimulator (WFSim) developed by Doekemeijer and Boersma [on GitHub](https://github.com/TUDelft-DataDrivenControl/WFSim) is utilized.
It used both for data generation in open loop and for closed loop testing.


## Evaluation 

The code in this repository was tested in the following environment:

* *Windows 10 Enterprise 21H2
* *Matlab* 2020a (Toolboxes used: Control System Toolbox and Optimization Toolbox)


