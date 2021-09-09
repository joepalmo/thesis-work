# thesis-work

### Table of Contents
**[Overview](#overview)**<br>
**[Set-Up and Build Dependencies](#set-up-and-build-dependencies)**<br>
**[Code Outline](#code-outline)**<br>
**[Tutorial](#tutorial)**<br>
**[To-Do](#to-do)**<br>
**[Miscellaneous](#miscellaneous)**<br>

## Overview

This repository contains the code related to my senior honors thesis in astronomy, titled: Simulated Scatter: Computational Modeling of (Sub)Stellar Accretion. My project consists of a Monte Carlo error propagation, and will generate a simulated distribution of accretion rates given uncertainty values for each leg of the process that goes into measuring and calculating the value. With our simulated distribution, we will be able to separate systematic uncertainties such as observational error, parameter estimation error, and scaling relationship uncertainty, from physical uncertainties. By isolating the physical scatter in the distribution, we learn more about deep questions surrounding the formation of planetary systems in our universe.

## Set-Up and Build Dependencies
Create and activate new python environment (Optional)

1. With Anaconda Navigator
```bash
$ conda create -n newenvname python=3.6
$ conda activate newenvname
```
2. Or
```bash
$ source create -n newenvname python=3.6
$ source activate newenvname
```

Then, install necessary dependencies:
```bash
(newenvname) $ pip install -r requirements.txt
```

Your environment should be set up and ready to run! Just remember to activate it using `conda` or `source` on the command line each time you are working with the code. And try to add any new dependencies to the requirements.txt file.

## Files Outline

1. `accretion.py` 
- Mdot Calculations - functions used to calculate mass accretion rate (Mdot) values for accreting objects
- Object Parameter Estimation - models and functions used to estimate various object parameters
- Uncertainty Distributions - code written to build Gaussian uncertainty distributions for a parameter given a mean and standard deviation
- Inverse Transform Sampling
- IMFs

2. `accretion_objects.py`
- class Accretion - single accreting star/BD/planet object
- class AccretionDistribution - distribution of accreting objects
- MC error propagation - these classes contain the functionality for running a Monte Carlo error propagation simulation
- plotting functions - these classes have built-in plotting functions for visualizing the results of a simulation

3. `accdb_updated.csv`
- old version of the CASPAR database. I know that there will be a newer version of this that you can use!

3. `StellarParams` folder
- contains the Baraffe models used within the simulation

4. `ageinterceptfunc.pickle`
- text encoded pickle file containing the original exponential fit model I derived between age and MvMdot intercept. This is now deprecated, because I eventually opted to hard-code the fit into `accretion.py` (see function age_intercept_exponential())

5. `requirements.txt` 
- dependencies

6. `Palmo_Thesis_Final.pdf`
- PDF document of my written up thesis


## Tutorial

For a very quick walkthrough tutorial on how to use the code, open up `tutorial.ipynb`.

For a bunch of (messy) notebooks that I used during my thesis for creating figures, tables, and developing certain sections of the modules, check out the notebooks folder.

## To-Do

Some of this is covered in the **Future Work** chapter of `Palmo_Thesis_Final.pdf`, or elsewhere.

1. 2D model interpolation
    - Currently, I use a 1-dimensional interpolation approach coupled with the discrete age Baraffe models.
    - Adapting this to a continuous age approach might make things more accurate
2. Uncertainty Interaction
    - Using a Monte Carlo approach for an error propagation makes it possible to build in interre- lationships between parameters and uncertainties. 
    - Currently, my simulation does not make use of any property covariances, but there are some key places where it may be implemented to make the results more accurate.
3. Add disk fragmentation
    - **Future Work** - **Model "Truth"**

4. Detection Limits 
    - **Future Work** - **Detection Limits**

5. Packaging Code
    - pip installable package

6. Keyword Arguments
    - while completing this thesis I learned a lot and made a TON of mistakes. There are places where my code is inefficient. One that comes to mind is making the inputs for the Accretion and AccretionDistribution classes follow a keyword arguments (**kwargs) structure.
    - There are many other things that might be made more efficient!

## Miscellaneous

- If you are interested in seeing the code for the Plotly-Dash app that I built as an accessible extension to this code, contact me!

- See **Appendix D: Scatter Metric Script** for the code that I wrote to generate the scatter statistics tables in Chapter 9.

