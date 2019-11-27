# Optimal Power Flow with Convex Restriction

This code provides a Julia+JuMP based implementation for Convex Restriction of AC Optimal Power Flow problem.


TL;DR

[1] Convex Restriction of Power Flow Feasibility Sets: Convex Restriction for OPF problem provides a sufficient convex condition for solving AC optimal power flow while satisfying operational constraints. 

[2] Feasible Path Identification in OPF with Sequential Convex Restriction: Using the Convex Restriction conditions, the AC OPF problem can be solved by a sequence of convex optimization problems. The output of the algorithm provides a guaranteed feasible path that satisfy the operational constraints.

[3] Robust Optimal Power Flow with Convex Restriction: 


ConVeX ReStriction (CVXRS) is a general technique that can analyse and optimize nonlinear systems using convex optimization.
You can find implementations of convex restrictions for other problems such as MPC, Neural Network and Robotics in the [project page](https://dclee131.github.io/research/2019/10/07/CVXRS.html).

### Installation Requirements

The script has been tested in Julia v1.1.
The primary script is in "CVXRS_OPF.jl", and it can be run without installation. 
To run the code, the follwing lists the necessary packages.

```julia
using JuMP, PowerModels, Ipopt, MosekTools, SparseArrays, LinearAlgebra, Plots
```
The project uses [JuMP](https://github.com/JuliaOpt/JuMP.jl) for modeling the QCQP problem and uses MOSEK as the solver.
Power flow data and equations are stored and solved based on [PowerModels](https://github.com/lanl-ansi/PowerModels.jl). 
Figures in Julia are created with [Plots](https://github.com/JuliaPlots/Plots.jl). Color-coded figures were drawn in MATLAB, and you can find relavant code in `plots/draw_boundary.m`.

## Quick Start

You can find example codes in the folder `example`.

```julia
## Import CVXRS_OPF functions
include("src/CVXRS_OPF.jl")

## Read network data using PowerModels.jl
network_data = PowerModels.parse_file("../../pglib-opf-master/pglib_opf_case118_ieee.m");

## Initiailize the network data by solving OPF problem
network_data=opf_initialization(network_data)

```

## Citing CVXRS_OPF

If you find this content useful for your research, please consider citing: 

[1] Convex Restriction of Power Flow Feasibility Set

    @article{lee2019convex,
      author={D. {Lee} and H. D. {Nguyen} and K. {Dvijotham} and K. {Turitsyn}},
      journal={IEEE Transactions on Control of Network Systems},
      title={Convex Restriction of Power Flow Feasibility Sets},
      year={2019}, volume={6}, number={3}, pages={1235-1245}, month={Sep.}
    }

[2] Feasible Path Identification in Optimal Power Flow with Sequential Convex Restriction

    @article{lee2019feasible,
      title={Feasible Path Identification in Optimal Power Flow with Sequential Convex Restriction},
      author={Lee, Dongchan and Turitsyn, Konstantin and Molzahn, Daniel K and Roald, Line A},
      journal={arXiv preprint arXiv:1906.09483},
      year={2019}
    }


