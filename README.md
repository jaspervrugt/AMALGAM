# AMALGAM: Multi-criteria optimization using multimethod adaptive search in MATLAB and Python

## Description

This toolbox infers the Pareto trade-off surface corresponding to a vector of competing objective functions. In principle, each Pareto solution is a different weighting of the objectives used. Therefore, one could use differential weighting and approximte the Pareto front using many different trials with a single objective optimization algorithm. This approach, however, has shown to not be particularly efficient (exception is MOEA/D). This MATLAB and Python toolbox implements adaptive multimethod search to ensure a fast, reliable and computationally efficient solution to multiobjective optimization problems. This method is called a multi-algorithm, genetically adaptive multiobjective, or AMALGAM, method, to evoke the image of a procedure that blends the attributes of the best available individual optimization algorithms. AMALGAM finds a well-distributed set of Pareto solutions within a single optimization run, and achieves a good performance compared to commonly used methods such as SPEA2, NSGA-II and MOEA/D. The AMALGAM toolbox implements multi-core (thread) evaluation of the offspring to speed-up inference of CPU-intensive system models, and provides convergence diagnostics and graphical output. Built-in case studies involving (among others) multimodality, high-dimensionality, bounded parameter spaces, dynamic simulation models, and distributed multi-core computation illustrate the main capabilities and functionalities of the AMALGAM toolbox. 

## Getting Started

### Installing: MATLAB

* Download and unzip the zip file 'MATLAB_code_AMALGAM_V2.0.zip' in a directory 'AMALGAM'.
* Add the toolbox to your MATLAB search path by running the script 'install_AMALGAM.m' available in the root of 'AMALGAM'
* You are ready to run the examples.

### Executing program

* After intalling, you can simply direct to each example folder and execute the local 'example_X.m' file.
* Please Make sure you read carefully the instructions (i.e., green comments) in 'install_AMALGAM.m' and the manual !!!  

### Installing: Python

* Download and unzip the zip file 'Python_code_AMALGAM_V2.0.zip' to a directory called 'AMALGAM'.

### Executing program

* Go to Command Prompt and directory of example_X in the root of 'AMALGAM'
* Now you can execute this example by typing 'python example_X.py'.
* Instructions can be found in the file 'AMALGAM.py' and in the manual !!!

## Authors

* Vrugt, Jasper A. (jasper@uci.edu) 

## Literature

1. Vrugt, J.A. (2024), Distribution-Based Model Evaluation and Diagnostics: Elicitability, Propriety, and Scoring Rules for Hydrograph Functionals, _Water Resources Research_, 60, e2023WR036710, https://doi.org/10.1029/2023WR036710
2. Vrugt, J.A., Multi-criteria optimization using the AMALGAM software package: Theory, concepts, and MATLAB implementation, UCI, 2015
3. Vrugt, J.A., B.A. Robinson, and J.M. Hyman (2009), Self-adaptive multimethod search for global optimization in real-parameter spaces, _IEEE Transactions on Evolutionary Computation_, 13(2), pp. 243-259, https://doi.org/10.1109/TEVC.2008.924428
4. Vrugt, J.A., P.H. Stauffer, T. Wöhling, B.A. Robinson, and V.V. Vesselinov (2008), Inverse modeling of subsurface flow and transport properties: A review with new developments, _Vadose Zone Journal_, 7(2), pp. 843-864, https://doi.org/10.2136/vzj2007.0078
5. Wöhling, T., J.A. Vrugt, and G.F. Barkle (2008), Comparison of three multiobjective optimization algorithms for inverse modeling of vadose zone hydraulic properties, _Soil Science Society of America Journal_, 72, 305 - 319, https://doi.org/10.2136/sssaj2007.0176
6. Wöhling, T., and J.A. Vrugt (2008), Combining multi-objective optimization and Bayesian model averaging to calibrate forecast ensembles of soil hydraulic models, _Water Resources Research_, 44, W12432, https://doi.org/10.1029/2008WR007154
7. Vrugt, J.A., and B.A. Robinson (2007), Improved evolutionary optimization from genetically adaptive multimethod search, _Proceedings of the National Academy of Sciences of the United States of America_, 104, pp. 708-711, https://doi.org/10.1073/pnas.061047110407
8. Vrugt, J.A., and B.A. Robinson (2007), Treatment of uncertainty using ensemble methods: Comparison of sequential data assimilation and Bayesian model averaging, _Water Resources Research_, 43, W01411, https://doi.org/10.1029/2005WR004838
9. Schoups, G.H., J.W. Hopmans, C.A. Young, J.A. Vrugt, and W.W. Wallender (2006), Multi-objective optimization of a regional spatially-distributed subsurface water flow model, _Journal of Hydrology_, pp. 20-48, 311(1-4), https://doi.org/10.1016/j.jhydrol.2005.01.001
10. Vrugt, J.A., M.P. Clark, C.G.H. Diks, Q. Duan, and B. A. Robinson (2006), Multi-objective calibration of forecast ensembles using Bayesian model averaging, _Geophysical Research Letters_, 33, L19817, https://doi.org/10.1029/2006GL027126
11. Vrugt, J.A., H.V. Gupta, L.A. Bastidas, W. Bouten, and S. Sorooshian (2003), Effective and efficient algorithm for multi-objective optimization of hydrologic models, _Water Resources Research_, 39(8), art. No. 1214, https://doi.org/10.1029/2002WR001746

## Version History

* 1.0
    * Initial Release
* 2.0
    * Improved postprocessing capabilities
    * Python implementation

## Built-in Case Studies

1. Example 1: ZDT1: test function
2. Example 2: ZDT2: test function
3. Example 3: ZDT3: test function
4. Example 4: ZDT4: test function
5. Example 5: ZDT6: test function
6. Example 6: ZDT6: test function, discrete parameter space
7. Example 7: DTLZ1: test function, 3 objectives
8. Example 8: DTLZ2: test function, 3 objectives
9. Example 9: DTLZ3: test function, 3 objectives
11. Example 11: Real-world example using rainfall-discharge modeling
12. Example 12: Watershed modeling using driven & nondriven hydrograph
13. Example 13: Bayesian model averaging: RMSE, IS and CRPS
14. Example 14: Multi-criteria BMA training temperature ensemble
15. Example 15: Multi-criteria BMA training sea-level pressure ensemble

## Acknowledgments

