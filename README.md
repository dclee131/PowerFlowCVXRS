# Optimal Power Flow with Convex Restriction

This code provides a Julia+JuMP based implementation for solving the AC Optimal Power Flow problem with Sequential Convex Restriction. The algorithm solves a sequence of convex optimization, in particular QCQP, to compute the optimal operating point.

### Installation Requirements

The script has been tested in Julia v1.1.
The primary script is in "CVXRS_OPF.jl", and it can be run without installation. 
To run the code, the follwing lists the necessary packages. The following packages are used currently:

```julia
JuMP
MosekTools
SparseArrays
LinearAlgebra
Plots
PowerModels
Ipopt
```

## Quick Start

Go to the root folder of `CVXRS_OPF`.

```julia
include("src/PowerCVXRS_polar.jl")

network_data = PowerModels.parse_file("../../pglib-opf-master/pglib_opf_case118_ieee.m");
network_data=opf_initialization(network_data, 0.0)

```

## Citing CVXRS_OPF

If you find this content useful for your research, please consider citing: 

Convex Restriction of Power Flow Feasibility Set

    @ARTICLE{lee2019convex,
      author={D. {Lee} and H. D. {Nguyen} and K. {Dvijotham} and K. {Turitsyn}},
      journal={IEEE Transactions on Control of Network Systems},
      title={Convex Restriction of Power Flow Feasibility Sets},
      year={2019}, volume={6}, number={3}, pages={1235-1245}, month={Sep.}
    }

Feasible Path Identification of Optimal Power Flow

    @ARTICLE{lee2019convex,
      author={D. {Lee} and H. D. {Nguyen} and K. {Dvijotham} and K. {Turitsyn}},
      journal={IEEE Transactions on Control of Network Systems},
      title={Convex Restriction of Power Flow Feasibility Sets},
      year={2019}, volume={6}, number={3},
      pages={1235-1245}, doi={10.1109/TCNS.2019.2930896}, month={Sep.}
    }


