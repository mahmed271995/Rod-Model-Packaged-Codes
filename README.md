# Rod Model Code (Packaged)

This repository contains MATLAB codes for a fully computational nonlinear rod model. For a detailed description of the rod model algorithm, please refer to the following [link](https://www.overleaf.com/read/qwysjxhrmtrw#b7f7cd).  
The code is structured for ease of use, allowing users to modify boundary conditions, numerical parameters, and material properties. However, certain critical files are protected to prevent accidental changes that could disrupt the functionality of the model.

Organization is as follows:
1. Forward Code: Has the packaged codes that solve the rod model equations in the forward approach i.e.: given a known constitutive law and known boundary conditions solve for the unknown parameters
    1. Linear: Has the rod model packaged codes considering linear constitutive relationship ($$\vec{q} = [\mathbf{B}] \vec{\kappa}$$)
    2. Nonlinear: Has the rod model packaged codes considering nonlinear constitutive relationship ($$\vec{q} = \sum_{i=0}^{5} \mathbf{B_i} \vec{\kappa}^i$$) more specifically,

    $$  \left\lbrace \begin{array}{c} q_1 \newline q_2 \newline q_3 \end{array} \right\rbrace = \left\lbrace \begin{array}{c} E I_1 \kappa_1 \newline E I_2 \kappa_2 + a \kappa_2^3 + b \kappa_2^5 \newline G J \kappa_3  \end{array}\right\rbrace $$
    
    3. Nonlinear and Non-homogeneous: Has the rod model packaged codes considering nonlinear and non-homogeneous constitutive relationship ($$\vec{q} = \sum_{i=0}^{5} \mathbf{B(s)_i} \vec{\kappa}^i$$) more specifically,  
    
    $$  \left\lbrace \begin{array}{c} q_1 \newline q_2 \newline q_3 \end{array} \right\rbrace = \left\lbrace \begin{array}{c} E I_1 \kappa_1 \newline a s^3 - b s^2 + c s + d)*(e \kappa^3 + f \kappa) \newline G J \kappa_3  \end{array}\right\rbrace $$   
2. Inverse Code: Has the code that solve the rod model equations in the inverse approach i.e.: given the deformation data law and known boundary conditions solve for the arguments of the constitutive law

## Instructions for running code
[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=mahmed271995/Rod-Model-Packaged-Codes) - Open the git repository on MATLAB online or Download the zip file  
Alternatively clone this repo using 
```
git clone https://github.com/mahmed271995/Rod-Model-Packaged-Codes
```
* Navigate to the `Linear` or `Nonlinear` folder
* Open `main.mlx` and click **Run** to execute the simulation
  - This generates plots, saves variable data, and creates a video file
* **Do not modify** the numerical parameters (`at`, `bt`, `gt`, `as`, `bs`, `gs`), as these ensure the model runs under an unconditionally stable condition. Changing them may cause the simulation to diverge
* The user can adjust the following parameters:
  - Spatial and temporal parameters: rod length (`Length`), spatial step size (`ds`), number of nodes (`N`), time step (`dt`), and total time steps (`nT`).
  - Material properties: Poisson’s ratio (`poiss`), moments of inertia (`I1`, `I2`, `I3`, defined for a circular cross-section), mass per unit length (`m`), and Young’s modulus (`E`).   
* Boundary conditions can be modified in the following files:
  - `leftbound.m` — defines the first boundary condition (currently fixed):
  ```
  v1 = 0, v2 = 0, v3 = 0, w1 = 0, w2 = 0, w3 = 0
  ```
  - `YG_new_direct.m` — defines the second boundary condition (currently sinusoidal bending load):  
   ```
   BV7  = 0, BV8  = 4 * sin(2 * pi * d * dt), BV9  = 0, BV10 = 0, BV11 = 0, BV12 = 0
   ```
    Here, `BV7–BV9` correspond to curvature components ($\kappa_1$, $\kappa_2$, $\kappa_3$) that map to bending moments **q** through the constitutive law, while `BV10–BV12` represent force components ($f_1$, $f_2$, $f_3$).








































































