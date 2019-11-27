# Optimal Power Flow with Convex Restriction

This code provides a Julia+JuMP based implementation for Convex Restriction of AC Power Flow Equations.

[1] Convex Restriction of Power Flow Feasibility Sets: Convex Restriction for OPF problem provides a sufficient convex condition for solving AC optimal power flow while satisfying operational constraints. 

[2] Feasible Path Identification: The method has been exteded to solve AC Optimal Power Flow problems with Sequential Convex Restriction [2]. The algorithm solves a sequence of convex optimization problems, in particular QCQP, to compute the optimal operating point.

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
The power flow data and equations are stored and solved based on [PowerModels](https://github.com/lanl-ansi/PowerModels.jl).

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

    @ARTICLE{lee2019convex,
      author={D. {Lee} and H. D. {Nguyen} and K. {Dvijotham} and K. {Turitsyn}},
      journal={IEEE Transactions on Control of Network Systems},
      title={Convex Restriction of Power Flow Feasibility Sets},
      year={2019}, volume={6}, number={3}, pages={1235-1245}, month={Sep.}
    }

[2] Feasible Path Identification of Optimal Power Flow

    @ARTICLE{lee2019convex,
      author={D. {Lee} and H. D. {Nguyen} and K. {Dvijotham} and K. {Turitsyn}},
      journal={IEEE Transactions on Control of Network Systems},
      title={Convex Restriction of Power Flow Feasibility Sets},
      year={2019}, volume={6}, number={3},
      pages={1235-1245}, doi={10.1109/TCNS.2019.2930896}, month={Sep.}
    }


