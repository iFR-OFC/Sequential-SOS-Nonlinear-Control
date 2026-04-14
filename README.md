# Sequential-SOS-Nonlinear-Control
[![Identifier](https://img.shields.io/badge/doi-10.18419%2Fdarus--5677-d45815.svg)](https://doi.org/10.18419/darus-5677)

This repository supplements the paper On the Practical Implementation of a Filter Line
Search Sequential Quadratic Programming Method
for Nonconvex Sum-of-squares Problems submitted to Mathematical Programming Computation. The preprint can be found on [arxiv](https://arxiv.org/abs/2602.02394).
It contains the case studies for different analysis and control problems.

### Requirements and Setup

  1. Download or clone [CaΣoS v1.0.0](https://github.com/iFR-OFC/casos/releases/tag/v1.0.0 )
  2. Download CasADi v3.6.7 and add it to your Matlab path.
  3. Download and install MOSEK v11.1.8, and add the solver to the Matlab path.  An acaddemic license can be obtained [here](https://www.mosek.com/products/academic-licenses/) 
  4. Add the CaΣoS root folder (the one that contains the directory +casos) to your Matlab path.

If installed correctly, you can excute the corresponding files as outlined below for reproduction.

In case you encounter any issues or for questions, please contact the corresponding authors.

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

## Onboarding Tutorial
For a simple onboarding process, we provide a tutorial.m file. This file implements a simpl example and guides the user through the solver build process. Due to its simplicity, this file is recommended for debugging purposes of the proposes filter line search algorithm.

## Running the Benchmarks:
1. Navigate to the corresponding folder.
2. Run the MATLAB script either for the sequential approach or for coordinate descent.

An exception is the N-link robot arm example (00_ROA/RobotArm). Use run_Nlink_benchmark.m to either run the benchmark from scratch or to just inspect the results, provided in the corresponding .mat file. (2026-04-07) 




## Citation
Please cite the following paper, if you use the sequential algorithm.

> J. Olucak and T. Cunis , ‘On the Practical Implementation of a Sequential
Quadratic Programming Algorithm for Nonconvex
Sum-of-squares Problems’,  submitted to Mathematical Programming Computation, 2026. Pre-print:(https://arxiv.org/abs/2602.02394)


<details>

<summary>Bibtex entry</summary>

```bibtex
@misc{olucak2026practicalimplementationsequentialquadratic,
      title={On the Practical Implementation of a Sequential Quadratic Programming Algorithm for Nonconvex Sum-of-squares Problems}, 
      author={Jan Olucak and Torbjørn Cunis},
      year={2026},
      eprint={2602.02394},
      archivePrefix={arXiv},
      primaryClass={math.OC},
      url={https://arxiv.org/abs/2602.02394}, 
}
```



