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

```python
# An example that shows how to use the UCB1 learning policy
# to make decisions between two arms based on their expected rewards.

# Import MABWiser Library
from mabwiser.mab import MAB, LearningPolicy, NeighborhoodPolicy

# Data
arms = ['Arm1', 'Arm2']
decisions = ['Arm1', 'Arm1', 'Arm2', 'Arm1']
rewards = [20, 17, 25, 9]

# Model 
mab = MAB(arms, LearningPolicy.UCB1(alpha=1.25))

# Train
mab.fit(decisions, rewards)

# Test
mab.predict()
```

```julia
include("src/PowerCVXRS_polar.jl")

#network_data = PowerModels.parse_file("cases/case9.m");
#network_data = PowerModels.parse_file("../../pglib-opf-master/pglib_opf_case14_ieee.m");
network_data = PowerModels.parse_file("../../pglib-opf-master/pglib_opf_case118_ieee.m");
#network_data = PowerModels.parse_file("../../pglib-opf-master/pglib_opf_case240_pserc.m");

result_opf=run_opf(network_data, ACPPowerModel, with_optimizer(Ipopt.Optimizer,print_level=0))
println("Local search OPF objective solution: ",round(result_opf["objective"],digits=2))

network_data=opf_initialization(network_data, 0.0)
#network_data=opf_initialization(network_data, 0.01)
#network_data=opf_initialization(network_data, "uniform")

[network_data["gen"][i]["ptc_factor"]=Int(network_data["bus"][string(gen["gen_bus"])]["bus_type"]==3) for (i,gen) in network_data["gen"]]
#[network_data["gen"][i]["ptc_factor"]=1/length(network_data["gen"])  for (i,gen) in network_data["gen"]]
#pg0=zeros(Float64, num_gen); [pg0[gen["index"]]=gen["gen_status"]*gen["pg"] for (i,gen) in network_data["gen"]]
#[network_data["gen"][i]["ptc_factor"]=gen["gen_status"]*gen["pg"]/sum(pg0) for (i,gen) in network_data["gen"]]

network_data=runpf(network_data)
test_runpf(network_data)
println("Initial OPF objective solution: ",round(network_data["cost"],digits=2))
violation_status, margin=check_violation(network_data)



```
    Local search OPF objective solution: 97213.61
    Initial OPF objective solution: 97213.61





    (Dict("line"=>[],"qg"=>[],"pg"=>[],"angle"=>[],"vmag"=>[]), 0.1)



```julia
max_iter_SCRS=10

(Σ0,γ0)=network_data["uncertainty"]
network_data["uncertainty"]=(Σ0,0.1)

num_gen=length(network_data["gen"])
result_cvxr=Dict("pg"=>zeros(num_gen,max_iter_SCRS+1),"vg"=>zeros(num_gen,max_iter_SCRS+1),"alpha"=>zeros(num_gen,max_iter_SCRS+1),
    "obj"=> zeros(max_iter_SCRS+1),"solver_time"=>zeros(max_iter_SCRS+1),"issue"=>0,"max_iter"=>max_iter_SCRS)

result_cvxr["obj"][1]=network_data["cost"]
[result_cvxr["pg"][gen["index"],1]=gen["gen_status"].*gen["pg"] for (i,gen) in network_data["gen"]]
[result_cvxr["vg"][gen["index"],1]=gen["vg"] for (i,gen) in network_data["gen"]]
[result_cvxr["alpha"][gen["index"],1]=gen["ptc_factor"] for (i,gen) in network_data["gen"]]

for iter=2:result_cvxr["max_iter"]
    network_data,sanity_check=cvxrs(network_data)
    network_data=runpf(network_data)
    test_runpf(network_data)
    test_cvxrs(network_data,sanity_check)
    violation_status, margin=check_violation(network_data)
    
    result_cvxr["obj"][iter]=network_data["cost"]
    [result_cvxr["pg"][gen["index"],iter]=gen["gen_status"].*gen["pg"] for (i,gen) in network_data["gen"]]
    [result_cvxr["vg"][gen["index"],iter]=gen["vg"] for (i,gen) in network_data["gen"]]
    [result_cvxr["alpha"][gen["index"],iter]=gen["ptc_factor"] for (i,gen) in network_data["gen"]]
    result_cvxr["solver_time"][iter]=sanity_check["solve_time"]
    
    dpg_step=norm(result_cvxr["pg"][:,iter]-result_cvxr["pg"][:,iter-1])
    dvg_step=norm(result_cvxr["vg"][:,iter]-result_cvxr["vg"][:,iter-1])
    println("Iteration ",iter ,": ",round(result_cvxr["obj"][iter],digits=2),"  ",round(sanity_check["obj"],digits=3),"   ", round(dpg_step,digits=4),"   ", round(dvg_step,digits=4))

    if dpg_step<1e-4 && dvg_step<1e-4; result_cvxr["max_iter"]=iter; break; end
end
```


## Citing CVXRS_OPF

If you found this content useful for your research, please consider citing: 

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


```julia

```


