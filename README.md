# Optimal Power Flow with Convex Restriction

This code provides implementation code for Convex Restriction of AC Optimal Power Flow problem. Please see citation for more details on methodology.

You can also find other applications of convex restriction in my [project page](https://dclee131.github.io/research/2019/10/07/CVXRS.html).

### Installation Requirements

The script has written in Julia v1.1, and can be run without installation. The code uses [JuMP](https://github.com/JuliaOpt/JuMP.jl) for modeling the QCQP problem and uses MOSEK as the solver.
Power flow data and equations are stored and solved based on [PowerModels.jl](https://github.com/lanl-ansi/PowerModels.jl).
Follwing packages are summary of necessary packages to run the code.

```julia
using JuMP, PowerModels, Ipopt, MosekTools, Gurobi, SparseArrays, LinearAlgebra, Plots
```

## Quick Start


```julia
## Import CVXRS_OPF functions
include("../src/main.jl")

## Read network data using PowerModels.jl
network_data = PowerModels.parse_file("../../pglib-opf-master/pglib_opf_case118_ieee.m");

## Initiailize the network data by solving OPF problem
network_data=opf_initialization(network_data)

## Run Sequential Convex Restriction with 5 iterations
network_data, result_cvxr = scrs(network_data, 5);
```

You can check the folder `tutorials` to find more examples.

## Citation

If you find this content useful for your research, please consider citing: 

[1] Convex Restriction of Power Flow Feasibility Set: proposed convex restriction of OPF constraints.

    @article{lee2019convex,
      author={Lee, Dongchan and Nguyen, Hung D. and Dvijotham, K. and Turitsyn, Konstantin},
      journal={IEEE Transactions on Control of Network Systems},
      title={Convex Restriction of Power Flow Feasibility Sets},
      year={2019}, volume={6}, number={3}, pages={1235-1245}, month={Sep.}
    }

[2] Feasible Path Identification in Optimal Power Flow with Sequential Convex Restriction: uses convex restriction sequentially to identify a feasible path.

    @article{lee2019feasible,
      title={Feasible Path Identification in Optimal Power Flow with Sequential Convex Restriction},
      author={Lee, Dongchan and Turitsyn, Konstantin and Molzahn, Daniel K and Roald, Line A},
      journal={arXiv preprint arXiv:1906.09483},
      year={2019}
    }

[3] Robust Optimal Power Flow with Convex Restriction: convex restriction provides a tractable way to solve robust OPF problem.

