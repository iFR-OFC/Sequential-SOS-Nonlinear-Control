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
constraintViolationCheck_comparison.m   # Comparison from Section IV-A
00_ROA		                            # Folder containing the region-of-attraction cases
├── F18                                 # region of attraction for the 7D F18 aircraft
├── RobotArm                            # Contains the N-link robot arm example. 
01_Control_design                       # Folder containing for nonlinear systems
├── Aircraft_4D                         # Contains the control design for nonlinear longitudinal motion of an aircraft
├── Satellite_6D                        # Control desing for control-affine nonlinear dynamics of a satellite
02_Reachabilty                          # Folder contains two cases to compute an inner-approximation of the reachable set
├──Aircarft_4D                          # Inner-approximation of the reachable set for an aircraft lonitudinal motion
├──VanDerPol_2D                         # Inner-approximation of the reachable set for the Van-der-Pol Oscillator
03_CBF_CLF                              # Files for sequential SOS to compute a compatible pair of CBF-CLF
```
All folder contain the script for either the sequential SOS or coordinate-descent implementation. Except for the N-link robot arm we also provide screenshots showing the computational stats. Note that these might differ on other machines.
