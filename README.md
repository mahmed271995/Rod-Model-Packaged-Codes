# Rod Model Code (Packaged)

This repository contains MATLAB codes for a fully computational nonlinear rod model. For a detailed description of the rod model algorithm, please refer to the following [link](https://www.overleaf.com/read/qwysjxhrmtrw#b7f7cd)  
The code is structured for ease of use, allowing users to modify boundary conditions, numerical parameters, and material properties. However, certain critical files are protected to prevent accidental changes that could disrupt the functionality of the model.

Organization is as follows:
1. Linear: Has the rod model packaged codes considering linear constitutive relationship (`$$\vec{q} = [\mathbf{B}] \vec{\kappa}$$`)

## Instructions for running code
[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=mahmed271995/Rod-Model-Packaged-Codes) - Open the git repository on MATLAB online or Download the zip file  
Alternatively clone this repo using 
```
git clone https://github.com/mahmed271995/Rod-Model-Packaged-Codes
```
1. Linear:
  * Navigate to the "Linear" folder
  * Open the file main.mlx
  * Click "Run" to execute the simulation. This will generate the plots, save the variable data, and create video file
**Important:**
  * Do **not** change the numerical parameters (`at`, `bt`, `gt`, `as`, `bs`, `gs`). These are set to ensure the model runs in an unconditionally stable condition. Changing them may cause the simulation to diverge.
  * The user can change the time spatial (Length, ds, N) and temporal (dt, nT) simulation parameters
  * The user can adjust the material properties: poison ratio (poiss), moment of Inertial (I1, I2 and I3, currently defined for a circular cross section), mass per unit lenght (m) and Young's modulus (E)








