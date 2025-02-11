# About 
This repository contains code to take genus 2 curves, search for genus 1 curves which may be gluable to them, and perform the actual gluings.

Unfortunately it is not possible for the entire workflow to be done with only one of SageMath and Magma, so this code uses SageMath's Magma interface, which means that the console output for the Magma code is invisible during normal usage.

# Contents


## src

Contains the code that does the actual work. By file:
- **compatible_curves.py**. The button that contains a method used to find all compatible elliptic curves in a single call.
- **compute_traces.py** Functions related to trace of Frobenius.
- **fast_gluing.m** An implementation of faster gluing algorithm as described in Section 6.3 of our paper.
- **gluing_utils.m** A wrapper of fast_gluing.m that simplifies a call to gluing functions.
- **gluing_wrapper.py** A Sage wrapper of gluing_wrapper.m
- **helpers.m** Contains some extra helper functions in Magma, mostly used in the reducible case.
- **helpers.py** A Sage wrapper of helpers.m
- **search_db.py** Functions that search LMFDB for compatible elliptic curves by the Frobenius trace constraints.
- **symplectic_test.m** A Magma implementation of symplectic test described in Section 5 of the paper.
- **symplectic_test.py** A Sage wrapper of symplectic_test.m

### 

## test

Contains tests for the code in src. Tests for the SageMath code are test.ipynb, while tests for the Magma code are in test_magma.ipynb. This directory also contains timing tests for sample gluings.