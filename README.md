# Inverse-Nonlinear-Schrodinger


| Build Status                                                 | Version                                                      |
| ------------------------------------------------------------ | ------------------------------------------------------------ |
| [![Build Status](https://travis-ci.com/nsoedjak/Inverse-Nonlinear-Schrodinger.svg?branch=master)](https://travis-ci.com/nsoedjak/Inverse-Nonlinear-Schrodinger) | ![Version](https://img.shields.io/badge/version-0.1.0-blue) |



This repo is built for nonlinear Schrödinger inversion problems. The package implements the forward simulation of nonlinear Schrödinger equations and enables flexible inversions of parameters using **technique**. Key features of this package include:

- **Feature 1**: Simultaneous reconstruction of multiple parameters.
- **Feature 2**: Experiments on multiple data acquisition regimes.


## Installation



To set up the Inverse-Nonlinear-Schrodinger project in MATLAB, you will need to run the `startup.m` script. This script adds the necessary directories (`src`, `test`, and `examples`) to the MATLAB search path.

Follow these steps to run the `startup.m` script:

1. **Clone the Repository (Optional)**  
    If you haven't already, clone the repository to your local machine.
    ```
    git clone https://github.com/nsoedjak/Inverse-Nonlinear-Schrodinger
    ```
    
2. **Open MATLAB**  
    Launch MATLAB on your computer.

3. **Set the Current Directory**  
    In the MATLAB interface, navigate to the current directory toolbar and set it to the `Inverse-Nonlinear-Schrodinger` directory where the `startup.m` file is located. 

4. **Run the `startup.m` Script**  
    In the MATLAB Command Window, enter the following command and press Enter:
    ```matlab
    startup
    ```

That's it! The necessary paths have been added for this MATLAB session. Whenever you restart MATLAB, make sure to run the `startup.m` script again or set the `Inverse-Nonlinear-Schrodinger` directory as your default directory to run the script automatically at startup.