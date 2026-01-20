# TCST-Sequential-Quadratic-Sum-of-squares-Programming-
This repository supplements the paper "Sequential Quadratic Sum-of-squares Programming for Nonlinear Systems Control Design and Analysis" submitted to IEEE  Transactions on Control Systems Technology. It contains the case studies for different analysis and control problems.

### Requirements and Setup
The example package itself does not not need to be installed. Only a [stable version](https://github.com/ifr-acso/casos/releases/latest) of CaΣoS, [MOSEK](https://www.mosek.com/downloads/) v11.1 is needed and CasADi v3.6.7 are required. 

  1.  Download CasADi v3.6.x and add it to your Matlab path.
  2. Download and install MOSEK, and add the solver to the Matlab path.
  3.   Add the CaΣoS root folder (the one that contains the directory +casos) to your Matlab path.

If installed correctly, you can excute the corresponding files as outlined below for reproduction.

### Folder Structure


```text
00_ROA/		  # Folder containing the region-of-attraction cases
├── F18       # region of attraction for the 7D F18 aircraft
├── RobotArm  # Contains the N-link robot arm example. 
01_control_design # Folder containing 
```
