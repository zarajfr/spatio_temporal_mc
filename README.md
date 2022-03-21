# Markov chain analysis of anti-refugee incidents in Germany
***
## Authors

 * Zahra Jafari, <zahra.jafari.17@ucl.ac.uk>
   * UCL Jill Dando Institute of Security and Crime Science, University College London
 * Michael Frith,
   * UCL Jill Dando Institute of Security and Crime Science, University College London
 * Toby Davies,
   * UCL Jill Dando Institute of Security and Crime Science, University College London
 * Shane Johnson,
   * UCL Jill Dando Institute of Security and Crime Science, University College London

# Introduction

This is a code repository for the paper "Spatio-Temporal Analysis of Anti-Refugee Incidents in Germany".
It provides all the information, data and code that are required in order to replicate the analyses presented in the paper.
This document provides the instruction and the description of the content of the repository. This analysis can be applied to any other spatio-temporal data once put in the format used in the source code. The details of the analyses are documented within the source file as comments.

# Requirements

The project uses Python and the source code should run on any standard operating system (i.e. Linux/Unix, MacOS, Windows).

## Python Dependencies

 - Python3.6+
 - Other dependencies can be found in the `requirements.txt`.

 To install them run `pip3 install -r requirements.txt`.

# Content of the repository

## Scripts
  - `germany_mc.py`
    - General module with the core functions to simulate the three analysis presented in the paper.

## Directories

  - `incidents_data`
    - Here the data has been presented at various spatio-temporal granularity. The analysis used 4 files from the 16 files presented here.
### Input format

The data is geocoded at four spatial granularity of state level, district level, association level, and municipality level. The temporal intervals are daily, monthly, biweekly, weekly, and daily. The directory includes 16 files corresponding the possible combinations, named by the initials of their space-time level; for example _SM.csv_ is state-monthly data. The temporal units are across the columns and the space is across the rows.

## Results

Simulation results are printed. The result is plotted for spatial and LISA Markov chains.
