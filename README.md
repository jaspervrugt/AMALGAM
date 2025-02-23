# AMALGAM: Multi-criteria optimization using multimethod adaptive search in MATLAB and Python

## Description

The evolutionary algorithm AMALGAM implements the novel concept of adaptive multimethod search to ensure a fast, reliable and computationally efficient solution to multiobjective optimization problems. The method finds a well-distributed set of Pareto solutions within a single optimization run, and achieves an excellent performance compared to commonly used methods such as SPEA2, NSGA-II and MOEA/D. In this paper, I review the basic elements of AMALGAM, provide a pseudo-code of the algorithm, and introduce MATLAB and Python toolboxes which provide scientists and engineers with an arsenal of options and utilities to solve multiobjective optimization problems involving (among others) multimodality, high-dimensionality, bounded parameter spaces, dynamic simulation models, and distributed multi-core computation. The AMALGAM toolbox supports parallel computing to permit inference of CPU-intensive system models, and provides convergence diagnostics and graphical output. Different built-in case studies are used to illustrate the main capabilities and functionalities of the AMALGAM toolbox.

## Getting Started

### Installing: MATLAB

* Download and unzip the zip file 'MATLAB_code_AMALGAM_V2.0.zip' in a directory 'AMALGAM'.
* Add the toolbox to your MATLAB search path by running the script 'install_AMALGAM.m' available in the root of 'AMALGAM'
* You are ready to run the examples.

### Executing program

* After intalling, you can simply direct to each example folder and execute the local 'example_X.m' file.
* Please Make sure you read carefully the instructions (i.e., green comments) in 'install_AMALGAM.m' and the manual in PDF !!!  

### Installing: Python

* Download and unzip the zip file 'Python_code_AMALGAM_V2.0.zip' to a directory called "AMALGAM".

### Executing program

* Go to Command Prompt and directory of example_X in the root of 'AMALGAM'
* Now you can execute this example by typing "python example_X.py".
* Instructions can be found in the file "AMALGAM.py" and in the manual !!!

## Authors

* Vrugt, Jasper A. (jasper@uci.edu) 

## Version History

* 1.0
    * Initial Release
* 2.0
    * Improved postprocessing capabilities
    * Python implementation
  
## Acknowledgments

