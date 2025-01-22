# About 
This repository contains code to take genus 2 curves, search for genus 1 curves which may be gluable to them, and perform the actual gluings.

Unfortunately it is not possible for the entire workflow to be done with only one of SageMath and Magma, so this code uses SageMath's Magma interface, which means that the console output for the Magma code is invisible during normal usage.

# Contents

## large_scale_search

Contains code to perform searches for gluings as well as the output of said code.

## src

Contains the code that does the actual work

## test

Contains tests for the code in src. Tests for the actual 