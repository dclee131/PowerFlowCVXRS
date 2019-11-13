# Optimal Power Flow with Convex Restriction

## Code

Go to the root folder of `CVXRS_OPF`.

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

    ┌ Info: Precompiling Ipopt [b6b21f68-93f8-5de0-b562-5493be1d77c9]
    └ @ Base loading.jl:1186


    [32m[info | PowerModels]: removing 3 cost terms from generator 32: Float64[][39m
    [32m[info | PowerModels]: removing 1 cost terms from generator 29: [3266.88, 0.0][39m
    [32m[info | PowerModels]: removing 3 cost terms from generator 1: Float64[][39m
    [32m[info | PowerModels]: removing 3 cost terms from generator 54: Float64[][39m
    [32m[info | PowerModels]: removing 3 cost terms from generator 2: Float64[][39m
    [32m[info | PowerModels]: removing 3 cost terms from generator 41: Float64[][39m
    [32m[info | PowerModels]: removing 1 cost terms from generator 51: [3504.34, 0.0][39m
    [32m[info | PowerModels]: removing 3 cost terms from generator 53: Float64[][39m
    [32m[info | PowerModels]: removing 3 cost terms from generator 27: Float64[][39m
    [32m[info | PowerModels]: removing 3 cost terms from generator 42: Float64[][39m
    [32m[info | PowerModels]: removing 3 cost terms from generator 33: Float64[][39m
    [32m[info | PowerModels]: removing 1 cost terms from generator 28: [3478.18, 0.0][39m
    [32m[info | PowerModels]: removing 3 cost terms from generator 50: Float64[][39m
    [32m[info | PowerModels]: removing 3 cost terms from generator 52: Float64[][39m
    [32m[info | PowerModels]: removing 1 cost terms from generator 26: [1605.6, 0.0][39m
    [32m[info | PowerModels]: removing 3 cost terms from generator 10: Float64[][39m
    [32m[info | PowerModels]: removing 3 cost terms from generator 24: Float64[][39m
    [32m[info | PowerModels]: removing 1 cost terms from generator 25: [2486.19, 0.0][39m
    [32m[info | PowerModels]: removing 3 cost terms from generator 23: Float64[][39m
    [32m[info | PowerModels]: removing 3 cost terms from generator 49: Float64[][39m
    [32m[info | PowerModels]: removing 1 cost terms from generator 5: [2498.34, 0.0][39m
    [32m[info | PowerModels]: removing 3 cost terms from generator 31: Float64[][39m
    [32m[info | PowerModels]: removing 3 cost terms from generator 43: Float64[][39m
    [32m[info | PowerModels]: removing 1 cost terms from generator 39: [3407.26, 0.0][39m
    [32m[info | PowerModels]: removing 3 cost terms from generator 34: Float64[][39m
    [32m[info | PowerModels]: removing 3 cost terms from generator 44: Float64[][39m
    [32m[info | PowerModels]: removing 3 cost terms from generator 17: Float64[][39m
    [32m[info | PowerModels]: removing 1 cost terms from generator 37: [2460.08, 0.0][39m
    [32m[info | PowerModels]: removing 3 cost terms from generator 47: Float64[][39m
    [32m[info | PowerModels]: removing 3 cost terms from generator 9: Float64[][39m
    [32m[info | PowerModels]: removing 1 cost terms from generator 12: [2222.1, 0.0][39m
    [32m[info | PowerModels]: removing 1 cost terms from generator 20: [2420.23, 0.0][39m
    [32m[info | PowerModels]: removing 1 cost terms from generator 6: [12458.2, 0.0][39m
    [32m[info | PowerModels]: removing 1 cost terms from generator 14: [2599.4, 0.0][39m
    [32m[info | PowerModels]: removing 3 cost terms from generator 7: Float64[][39m
    [32m[info | PowerModels]: removing 3 cost terms from generator 8: Float64[][39m
    [32m[info | PowerModels]: removing 3 cost terms from generator 19: Float64[][39m
    [32m[info | PowerModels]: removing 3 cost terms from generator 4: Float64[][39m
    [32m[info | PowerModels]: removing 1 cost terms from generator 22: [2727.73, 0.0][39m
    [32m[info | PowerModels]: removing 1 cost terms from generator 11: [2894.83, 0.0][39m
    [32m[info | PowerModels]: removing 3 cost terms from generator 35: Float64[][39m
    [32m[info | PowerModels]: removing 3 cost terms from generator 13: Float64[][39m
    [32m[info | PowerModels]: removing 3 cost terms from generator 15: Float64[][39m
    [32m[info | PowerModels]: removing 3 cost terms from generator 16: Float64[][39m
    [32m[info | PowerModels]: removing 1 cost terms from generator 40: [2460.51, 0.0][39m
    [32m[info | PowerModels]: removing 1 cost terms from generator 21: [1667.39, 0.0][39m
    [32m[info | PowerModels]: removing 3 cost terms from generator 38: Float64[][39m
    [32m[info | PowerModels]: removing 1 cost terms from generator 46: [2864.95, 0.0][39m
    [32m[info | PowerModels]: removing 1 cost terms from generator 45: [1261.22, 0.0][39m
    [32m[info | PowerModels]: removing 3 cost terms from generator 36: Float64[][39m
    [32m[info | PowerModels]: removing 3 cost terms from generator 48: Float64[][39m
    [32m[info | PowerModels]: removing 3 cost terms from generator 18: Float64[][39m
    [32m[info | PowerModels]: removing 1 cost terms from generator 30: [2575.84, 0.0][39m
    [32m[info | PowerModels]: removing 3 cost terms from generator 3: Float64[][39m
    
    ******************************************************************************
    This program contains Ipopt, a library for large-scale nonlinear optimization.
     Ipopt is released as open source code under the Eclipse Public License (EPL).
             For more information visit http://projects.coin-or.org/Ipopt
    ******************************************************************************
    
    Local search OPF objective solution: 97213.61
    Initial OPF objective solution: 97213.61





    (Dict("line"=>[],"qg"=>[],"pg"=>[],"angle"=>[],"vmag"=>[]), 0.1)




```julia

```


```julia

```


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

    MOSEK warning 705: #1 (nearly) zero elements are specified in sparse row ''(3196) of matrix 'A'.
    MOSEK warning 705: #1 (nearly) zero elements are specified in sparse row ''(3198) of matrix 'A'.
    MOSEK warning 705: #1 (nearly) zero elements are specified in sparse row ''(3365) of matrix 'A'.
    MOSEK warning 705: #1 (nearly) zero elements are specified in sparse row ''(3409) of matrix 'A'.
    MOSEK warning 705: #1 (nearly) zero elements are specified in sparse row ''(3447) of matrix 'A'.
    MOSEK warning 705: #1 (nearly) zero elements are specified in sparse row ''(3449) of matrix 'A'.
    MOSEK warning 705: #1 (nearly) zero elements are specified in sparse row ''(3616) of matrix 'A'.
    MOSEK warning 705: #1 (nearly) zero elements are specified in sparse row ''(3660) of matrix 'A'.


    ┌ Warning: CVXRS was not optimal. (SLOW_PROGRESS)
    └ @ Main /Users/dc/Documents/Github/CVXRS_master/PowerCVXRS/PowerCVXRS.jl:650


    Iteration 2: 104879.07  113487.182   8.4549   0.0523
    MOSEK warning 705: #1 (nearly) zero elements are specified in sparse row ''(3196) of matrix 'A'.
    MOSEK warning 705: #1 (nearly) zero elements are specified in sparse row ''(3198) of matrix 'A'.
    MOSEK warning 705: #1 (nearly) zero elements are specified in sparse row ''(3365) of matrix 'A'.
    MOSEK warning 705: #1 (nearly) zero elements are specified in sparse row ''(3409) of matrix 'A'.
    MOSEK warning 705: #1 (nearly) zero elements are specified in sparse row ''(3447) of matrix 'A'.
    MOSEK warning 705: #1 (nearly) zero elements are specified in sparse row ''(3449) of matrix 'A'.
    MOSEK warning 705: #1 (nearly) zero elements are specified in sparse row ''(3616) of matrix 'A'.
    MOSEK warning 705: #1 (nearly) zero elements are specified in sparse row ''(3660) of matrix 'A'.


    ┌ Warning: CVXRS was not optimal. (SLOW_PROGRESS)
    └ @ Main /Users/dc/Documents/Github/CVXRS_master/PowerCVXRS/PowerCVXRS.jl:650


    Iteration 3: 105861.82  113886.578   0.3324   0.0147
    MOSEK warning 705: #1 (nearly) zero elements are specified in sparse row ''(3196) of matrix 'A'.
    MOSEK warning 705: #1 (nearly) zero elements are specified in sparse row ''(3198) of matrix 'A'.
    MOSEK warning 705: #1 (nearly) zero elements are specified in sparse row ''(3365) of matrix 'A'.
    MOSEK warning 705: #1 (nearly) zero elements are specified in sparse row ''(3409) of matrix 'A'.
    MOSEK warning 705: #1 (nearly) zero elements are specified in sparse row ''(3447) of matrix 'A'.
    MOSEK warning 705: #1 (nearly) zero elements are specified in sparse row ''(3449) of matrix 'A'.
    MOSEK warning 705: #1 (nearly) zero elements are specified in sparse row ''(3616) of matrix 'A'.
    MOSEK warning 705: #1 (nearly) zero elements are specified in sparse row ''(3660) of matrix 'A'.


    ┌ Warning: CVXRS was not optimal. (SLOW_PROGRESS)
    └ @ Main /Users/dc/Documents/Github/CVXRS_master/PowerCVXRS/PowerCVXRS.jl:650


    Iteration 4: 106809.24  115314.091   0.6049   0.0084
    MOSEK warning 705: #1 (nearly) zero elements are specified in sparse row ''(3196) of matrix 'A'.
    MOSEK warning 705: #1 (nearly) zero elements are specified in sparse row ''(3198) of matrix 'A'.
    MOSEK warning 705: #1 (nearly) zero elements are specified in sparse row ''(3365) of matrix 'A'.
    MOSEK warning 705: #1 (nearly) zero elements are specified in sparse row ''(3447) of matrix 'A'.
    MOSEK warning 705: #1 (nearly) zero elements are specified in sparse row ''(3449) of matrix 'A'.
    MOSEK warning 705: #1 (nearly) zero elements are specified in sparse row ''(3616) of matrix 'A'.


    ┌ Warning: CVXRS was not optimal. (SLOW_PROGRESS)
    └ @ Main /Users/dc/Documents/Github/CVXRS_master/PowerCVXRS/PowerCVXRS.jl:650


    Iteration 5: 106557.38  115422.539   10.1719   0.0074
    MOSEK warning 705: #1 (nearly) zero elements are specified in sparse row ''(3196) of matrix 'A'.
    MOSEK warning 705: #1 (nearly) zero elements are specified in sparse row ''(3198) of matrix 'A'.
    MOSEK warning 705: #1 (nearly) zero elements are specified in sparse row ''(3365) of matrix 'A'.
    MOSEK warning 705: #1 (nearly) zero elements are specified in sparse row ''(3409) of matrix 'A'.
    MOSEK warning 705: #1 (nearly) zero elements are specified in sparse row ''(3447) of matrix 'A'.
    MOSEK warning 705: #1 (nearly) zero elements are specified in sparse row ''(3449) of matrix 'A'.
    MOSEK warning 705: #1 (nearly) zero elements are specified in sparse row ''(3616) of matrix 'A'.
    MOSEK warning 705: #1 (nearly) zero elements are specified in sparse row ''(3660) of matrix 'A'.


    ┌ Warning: CVXRS was not optimal. (SLOW_PROGRESS)
    └ @ Main /Users/dc/Documents/Github/CVXRS_master/PowerCVXRS/PowerCVXRS.jl:650


    Iteration 6: 106051.56  111768.095   10.1712   0.0059
    MOSEK warning 705: #1 (nearly) zero elements are specified in sparse row ''(3196) of matrix 'A'.
    MOSEK warning 705: #1 (nearly) zero elements are specified in sparse row ''(3198) of matrix 'A'.
    MOSEK warning 705: #1 (nearly) zero elements are specified in sparse row ''(3365) of matrix 'A'.
    MOSEK warning 705: #1 (nearly) zero elements are specified in sparse row ''(3447) of matrix 'A'.
    MOSEK warning 705: #1 (nearly) zero elements are specified in sparse row ''(3449) of matrix 'A'.
    MOSEK warning 705: #1 (nearly) zero elements are specified in sparse row ''(3616) of matrix 'A'.


    ┌ Warning: CVXRS was not optimal. (SLOW_PROGRESS)
    └ @ Main /Users/dc/Documents/Github/CVXRS_master/PowerCVXRS/PowerCVXRS.jl:650


    Iteration 7: 106735.71  115709.251   0.1555   0.0051
    MOSEK warning 705: #1 (nearly) zero elements are specified in sparse row ''(3196) of matrix 'A'.
    MOSEK warning 705: #1 (nearly) zero elements are specified in sparse row ''(3198) of matrix 'A'.
    MOSEK warning 705: #1 (nearly) zero elements are specified in sparse row ''(3365) of matrix 'A'.
    MOSEK warning 705: #1 (nearly) zero elements are specified in sparse row ''(3409) of matrix 'A'.
    MOSEK warning 705: #1 (nearly) zero elements are specified in sparse row ''(3447) of matrix 'A'.
    MOSEK warning 705: #1 (nearly) zero elements are specified in sparse row ''(3449) of matrix 'A'.
    MOSEK warning 705: #1 (nearly) zero elements are specified in sparse row ''(3616) of matrix 'A'.
    MOSEK warning 705: #1 (nearly) zero elements are specified in sparse row ''(3660) of matrix 'A'.


    ┌ Warning: CVXRS was not optimal. (SLOW_PROGRESS)
    └ @ Main /Users/dc/Documents/Github/CVXRS_master/PowerCVXRS/PowerCVXRS.jl:650


    Iteration 8: 105282.15  110637.618   0.2554   0.0047
    MOSEK warning 705: #1 (nearly) zero elements are specified in sparse row ''(3196) of matrix 'A'.
    MOSEK warning 705: #1 (nearly) zero elements are specified in sparse row ''(3198) of matrix 'A'.
    MOSEK warning 705: #1 (nearly) zero elements are specified in sparse row ''(3365) of matrix 'A'.
    MOSEK warning 705: #1 (nearly) zero elements are specified in sparse row ''(3409) of matrix 'A'.
    MOSEK warning 705: #1 (nearly) zero elements are specified in sparse row ''(3447) of matrix 'A'.
    MOSEK warning 705: #1 (nearly) zero elements are specified in sparse row ''(3449) of matrix 'A'.
    MOSEK warning 705: #1 (nearly) zero elements are specified in sparse row ''(3616) of matrix 'A'.
    MOSEK warning 705: #1 (nearly) zero elements are specified in sparse row ''(3660) of matrix 'A'.


    ┌ Warning: CVXRS was not optimal. (SLOW_PROGRESS)
    └ @ Main /Users/dc/Documents/Github/CVXRS_master/PowerCVXRS/PowerCVXRS.jl:650


    Iteration 9: 106177.66  115154.939   0.1746   0.0043
    MOSEK warning 705: #1 (nearly) zero elements are specified in sparse row ''(3196) of matrix 'A'.
    MOSEK warning 705: #1 (nearly) zero elements are specified in sparse row ''(3198) of matrix 'A'.
    MOSEK warning 705: #1 (nearly) zero elements are specified in sparse row ''(3365) of matrix 'A'.
    MOSEK warning 705: #1 (nearly) zero elements are specified in sparse row ''(3409) of matrix 'A'.
    MOSEK warning 705: #1 (nearly) zero elements are specified in sparse row ''(3447) of matrix 'A'.
    MOSEK warning 705: #1 (nearly) zero elements are specified in sparse row ''(3449) of matrix 'A'.
    MOSEK warning 705: #1 (nearly) zero elements are specified in sparse row ''(3616) of matrix 'A'.
    MOSEK warning 705: #1 (nearly) zero elements are specified in sparse row ''(3660) of matrix 'A'.


    ┌ Warning: CVXRS was not optimal. (SLOW_PROGRESS)
    └ @ Main /Users/dc/Documents/Github/CVXRS_master/PowerCVXRS/PowerCVXRS.jl:650


    Iteration 10: 105822.47  111112.751   0.1096   0.004



```julia
result_cvxr["alpha"]
```




    3×11 Array{Float64,2}:
     0.282114  0.282114  0.282114  0.282114  …  0.282114  0.282114  0.282114  0.0
     0.421985  0.421985  0.421985  0.421985     0.421985  0.421985  0.421985  0.0
     0.295901  0.295901  0.295901  0.295901     0.295901  0.295901  0.295901  0.0




```julia
### Procedure for plots
resolution=20
plot_rng=[-0.5 3 0 3]
plot_bus=["2","3"]
U1_plot,U2_plot,exact_plot,cvxrs_plot=plot2D(network_data,plot_bus,plot_rng,resolution)

u_plot0=[network_data["gen"][plot_bus[1]]["pg"], network_data["gen"][plot_bus[2]]["pg"]]

pyplot()
default(show = true)
#contour(U1_plot,U2_plot,exact_plot)
contour(U1_plot,U2_plot,exact_plot,levels=0,color=:blues)
contour!(U1_plot,U2_plot,cvxrs_plot,levels=0,color=:greens)
scatter!([u_plot0[1]],[u_plot0[2]],markersize=6, c=:red)
#plot!(γ_plot*cos.(0:pi/50:2*pi).+u_plot0[1],γ_plot*sin.(0:pi/50:2*pi).+u_plot0[2])


```

    progress: 1/20
    progress: 2/20
    progress: 3/20
    progress: 4/20
    progress: 5/20


    ┌ Warning: CVXRS was not optimal. (SLOW_PROGRESS)
    └ @ Main /Users/dc/Documents/Github/CVXRS_master/PowerCVXRS/PowerCVXRS.jl:638


    progress: 6/20


    ┌ Warning: CVXRS was not optimal. (SLOW_PROGRESS)
    └ @ Main /Users/dc/Documents/Github/CVXRS_master/PowerCVXRS/PowerCVXRS.jl:638


    progress: 7/20
    progress: 8/20


    ┌ Warning: CVXRS was not optimal. (SLOW_PROGRESS)
    └ @ Main /Users/dc/Documents/Github/CVXRS_master/PowerCVXRS/PowerCVXRS.jl:638


    progress: 9/20


    ┌ Warning: CVXRS was not optimal. (SLOW_PROGRESS)
    └ @ Main /Users/dc/Documents/Github/CVXRS_master/PowerCVXRS/PowerCVXRS.jl:638


    progress: 10/20
    progress: 11/20
    progress: 12/20
    progress: 13/20
    progress: 14/20
    progress: 15/20
    progress: 16/20
    progress: 17/20
    progress: 18/20
    progress: 19/20



![png](output_6_9.png)


    progress: 20/20





![png](output_6_11.png)




```julia

```
