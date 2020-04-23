using JuMP, MosekTools
using SparseArrays, LinearAlgebra, Plots
using PowerModels, Ipopt

#OPT_version="Gurobi"
OPT_version="Mosek"

const GRB_ENV = Gurobi.Env()

function makeYbus(network_data,phase_shift=false)
    num_bus=length(network_data["bus"])
    num_line=length(network_data["branch"])
    num_sh=length(network_data["shunt"])

    idx_fr=zeros(Int64, num_line); [idx_fr[branch["index"]]=network_data["bus"][string(branch["f_bus"])]["idx"] for (i,branch) in network_data["branch"]]
    idx_to=zeros(Int64, num_line); [idx_to[branch["index"]]=network_data["bus"][string(branch["t_bus"])]["idx"] for (i,branch) in network_data["branch"]]
    idx_sh=zeros(Int64, num_sh); [idx_sh[shunt["index"]]=network_data["bus"][string(shunt["shunt_bus"])]["idx"] for (i,shunt) in network_data["shunt"]]

    line_status=zeros(Int64, num_line); [line_status[branch["index"]]=branch["br_status"] for (i,branch) in network_data["branch"]]
    y_line=zeros(Complex{Float64},num_line); [y_line[branch["index"]]=branch["br_status"]/(branch["br_r"]+1im*branch["br_x"]) for (i,branch) in network_data["branch"]]
    b_c=zeros(Float64,num_line); [b_c[branch["index"]]=branch["br_status"]*branch["b_to"]+branch["b_fr"] for (i,branch) in network_data["branch"]]
    #tap=zeros(Complex{Float64},num_line); [tap[branch["index"]]=branch["tap"]*exp(1im*pi/180*branch["shift"]) for (i,branch) in network_data["branch"]]
    tap=zeros(Complex{Float64},num_line); [tap[branch["index"]]=branch["tap"]*exp(1im*branch["shift"]) for (i,branch) in network_data["branch"]]
    y_sh=zeros(Complex{Float64},num_sh); [y_sh[shunt["index"]]=shunt["gs"]+1im*shunt["bs"] for (i,shunt) in network_data["shunt"]]#/mpc_input.baseMVA

    E_fr=sparse(idx_fr,1:num_line,ones(num_line),num_bus,num_line)
    E_to=sparse(idx_to,1:num_line,ones(num_line),num_bus,num_line)
    E=E_fr-E_to
    Csh= sparse(idx_sh, (1:num_sh), ones(num_sh), num_bus, num_sh)
    y_tt=y_line+1im*b_c/2
    y_ff=y_tt./(tap.*conj(tap))
    y_ft=-y_line./conj(tap)
    y_tf=-y_line./tap
    Yf=sparse([1:num_line; 1:num_line], [idx_fr; idx_to], [y_ff; y_ft], num_line, num_bus)
    Yt=sparse([1:num_line; 1:num_line], [idx_fr; idx_to], [y_tf; y_tt], num_line, num_bus)
    Y=E_to*Diagonal(y_tt)*E_to'+E_to*Diagonal(y_tf)*E_fr'+E_fr*Diagonal(y_ft)*E_to'+E_fr*Diagonal(y_ff)*E_fr'+Csh*Diagonal(y_sh)*Csh'
    #Y=E_fr*Yf+E_to*Yt+Csh*Diagonal(y_sh)*Csh'


    if phase_shift==false
        Y_cos=sparse([idx_fr; idx_to],[1:num_line; 1:num_line],[y_ft;y_tf],num_bus,num_line) #E_fr*Diagonal(y_ft)+E_to*Diagonal(y_tf)
        Y_sin=sparse([idx_fr; idx_to],[1:num_line; 1:num_line],[y_ft;-y_tf],num_bus,num_line) #E_fr*Diagonal(y_ft)-E_to*Diagonal(y_tf)
        Y_diag=E_to*Diagonal(y_tt)*E_to'+E_fr*Diagonal(y_ff)*E_fr'+Csh*Diagonal(y_sh)*Csh'

        M=[Matrix{Float64}(I,num_bus,num_bus) zeros(num_bus,num_bus) -real(Y_cos) -imag(Y_sin) -real(Y_diag);
           zeros(num_bus,num_bus) Matrix{Float64}(I,num_bus,num_bus)  imag(Y_cos) -real(Y_sin)  imag(Y_diag)]

        M_line=[zeros(num_line,2*num_bus)  real(Diagonal(y_ft)) imag(Diagonal(y_ft))  real(Diagonal(y_ff)*E_fr');
                zeros(num_line,2*num_bus)  real(Diagonal(y_tf)) -imag(Diagonal(y_tf))   real(Diagonal(y_tt)*E_to');
                zeros(num_line,2*num_bus) -imag(Diagonal(y_ft)) real(Diagonal(y_ft)) -imag(Diagonal(y_ff)*E_fr');
                zeros(num_line,2*num_bus) -imag(Diagonal(y_tf)) -real(Diagonal(y_tf))  -imag(Diagonal(y_tt)*E_to')]
    elseif phase_shift==true
        θ0=zeros(Float64, num_bus); [θ0[bus["idx"]]=bus["va"] for (i,bus) in network_data["bus"]]#*pi/180
        φ0=E'*θ0

        Y_plus=E_fr*Diagonal(y_ft.*exp.(-1im*φ0))+E_to*Diagonal(y_tf.*exp.(1im*φ0))
        Y_minus=E_fr*Diagonal(y_ft.*exp.(-1im*φ0))-E_to*Diagonal(y_tf.*exp.(1im*φ0))
        Y_diag=E_to*Diagonal(y_tt)*E_to'+E_fr*Diagonal(y_ff)*E_fr'+Csh*Diagonal(y_sh)*Csh'

        M=[Matrix{Float64}(I,num_bus,num_bus) zeros(num_bus,num_bus) -real(Y_plus) -imag(Y_minus) -real(Y_diag);
           zeros(num_bus,num_bus) Matrix{Float64}(I,num_bus,num_bus)  imag(Y_plus) -real(Y_minus)  imag(Y_diag)]

       M_line=[zeros(num_line,2*num_bus)  real(Diagonal(y_ft.*exp.(-1im*φ0))) imag(Diagonal(y_ft.*exp.(-1im*φ0)))  real(Diagonal(y_ff)*E_fr');
               zeros(num_line,2*num_bus)  real(Diagonal(y_tf.*exp.(1im*φ0))) -imag(Diagonal(y_tf.*exp.(1im*φ0)))   real(Diagonal(y_tt)*E_to');
               zeros(num_line,2*num_bus) -imag(Diagonal(y_ft.*exp.(-1im*φ0))) real(Diagonal(y_ft.*exp.(-1im*φ0))) -imag(Diagonal(y_ff)*E_fr');
               zeros(num_line,2*num_bus) -imag(Diagonal(y_tf.*exp.(1im*φ0))) -real(Diagonal(y_tf.*exp.(1im*φ0)))  -imag(Diagonal(y_tt)*E_to')]
    end

    #M_eq=M[[1:num_bus; num_bus.+idx_pq],:]
    return Y,Yf,Yt,M,M_line
end

function runpf(network_data)
    num_bus=length(network_data["bus"])
    num_gen=length(network_data["gen"])
    num_load=length(network_data["load"])

    bus_type=zeros(Int64, num_bus); [bus_type[bus["idx"]]=bus["bus_type"] for (i,bus) in network_data["bus"]]
    slack_bus=findall(bus_type.==3)[1]
    idx_nslack=setdiff(1:num_bus,slack_bus)
    idx_gen=zeros(Int64, num_gen); [idx_gen[gen["index"]]=network_data["bus"][string(gen["gen_bus"])]["idx"] for (i,gen) in network_data["gen"]]
    idx_load=zeros(Int64, num_load); [idx_load[load["index"]]=network_data["bus"][string(load["load_bus"])]["idx"] for (i,load) in network_data["load"]]
    idx_pq= findall(bus_type.==1)

    gencost=zeros(Float64,num_gen,3); [gencost[gen["index"],4-size(gen["cost"],1):end]=gen["gen_status"]*gen["cost"] for (i,gen) in network_data["gen"]]

    Cg=sparse(idx_gen,(1:num_gen),ones(num_gen),num_bus,num_gen)
    Cl=sparse(idx_load,(1:num_load),ones(num_load),num_bus,num_load)

    gen_status=zeros(Int64, num_gen); [gen_status[gen["index"]]=gen["gen_status"] for (i,gen) in network_data["gen"]]
    load_status=zeros(Int64, num_load); [load_status[load["index"]]=load["status"] for (i,load) in network_data["load"]]
    pg0=zeros(Float64, num_gen); [pg0[gen["index"]]=gen["gen_status"].*gen["pg"] for (i,gen) in network_data["gen"]]
    #qg0=zeros(Float64, num_gen); [qg0[gen["index"]]=gen["gen_status"].*gen["qg"] for (i,gen) in network_data["gen"]]
    pl0=zeros(Float64, num_load); [pl0[load["index"]]=load["status"].*load["pd"] for (i,load) in network_data["load"]]
    ql0=zeros(Float64, num_load); [ql0[load["index"]]=load["status"].*load["qd"] for (i,load) in network_data["load"]]
    pinj0=Cg*pg0-Cl*pl0; qinj0=-Cl*ql0
    Sinj=pinj0+1im*qinj0

    Y,Yf,Yt=makeYbus(network_data)

    # Initialization
    vm=zeros(Float64, num_bus); [vm[bus["idx"]]=bus["vm"] for (i,bus) in network_data["bus"]]
    [vm[network_data["bus"][string(gen["gen_bus"])]["idx"]]=gen["vg"] for (i,gen) in network_data["gen"]]
    va=zeros(Float64, num_bus); [va[bus["idx"]]=bus["va"] for (i,bus) in network_data["bus"]]
    va[slack_bus]=0
    delta=0
    #if ~haskey(network_data["gen"][network_data["gen"].keys[1]],"ptc_factor")
    #    alpha=Cg[slack_bus,:]
    #    [network_data["gen"][i]["ptc_factor"]=alpha[gen["index"]]  for (i,gen) in network_data["gen"]]
    #end
    alpha=zeros(Float64, num_gen); [alpha[gen["index"]]=gen["ptc_factor"] for (i,gen) in network_data["gen"]]

    max_iter=30
    nf=zeros(max_iter); ndx=zeros(max_iter)

    for iter=1:max_iter
        v_cpx=vm.*cos.(va)+1im*vm.*sin.(va)
        S_bal=v_cpx.*conj(Y*v_cpx)-Sinj-Cg*alpha*delta
        f=[real(S_bal); imag(S_bal[idx_pq])]
        J1=real(Diagonal(v_cpx)*conj(Y.*(ones(num_bus)*transpose(1im*v_cpx)))+Diagonal(1im*v_cpx)*Diagonal(conj(Y*v_cpx)))
        J2=real(Diagonal(v_cpx)*conj(Y.*(ones(num_bus)*transpose(v_cpx./vm)))+Diagonal(v_cpx./vm)*Diagonal(conj(Y*v_cpx)))
        J3=imag(Diagonal(v_cpx)*conj(Y.*(ones(num_bus)*transpose(1im*v_cpx)))+Diagonal(1im*v_cpx)*Diagonal(conj(Y*v_cpx)))
        J4=imag(Diagonal(v_cpx)*conj(Y.*(ones(num_bus)*transpose(v_cpx./vm)))+Diagonal(v_cpx./vm)*Diagonal(conj(Y*v_cpx)))
        J =[J1[:,idx_nslack] J2[:,idx_pq] -Cg*alpha; J3[idx_pq,idx_nslack] J4[idx_pq,idx_pq] zeros(size(idx_pq))]
        dx=-Matrix(J)\Vector(f)

        va[idx_nslack]=va[idx_nslack]+dx[1:num_bus-1]
        vm[idx_pq]=vm[idx_pq]+dx[num_bus:end-1]
        delta=delta+dx[end]

        nf[iter]=norm(f); ndx[iter]=norm(dx);
        if isnan(nf[iter]) || ((nf[iter]<1e-8)&&(ndx[iter]<1e-8)); break; end;
    end

    [network_data["bus"][i]["vm"]=vm[bus["idx"]]  for (i,bus) in network_data["bus"]]
    [network_data["bus"][i]["va"]=va[bus["idx"]]  for (i,bus) in network_data["bus"]]
    network_data["delta"]=delta

    v_cpx=vm.*cos.(va)+1im*vm.*sin.(va)
    S_inj=v_cpx.*conj(Y*v_cpx)+1im*Cl*ql0
    #qg_inj=Cg'*Diagonal(1 ./(sum(Cg,dims=2)[:]).^2)*imag(S_inj)
    qg_inj=Cg'*Diagonal(1 ./sum(Cg,dims=2)[:])*imag(S_inj)
    [network_data["gen"][i]["qg"]=qg_inj[gen["index"]]  for (i,gen) in network_data["gen"]]

    network_data["cost"]=gencost[:,1]'*(pg0+alpha.*delta).^2+gencost[:,2]'*(pg0+alpha.*delta)+sum(gencost[:,3])

    if nf[end]>1e-8; network_data["delta"]=NaN; end #@warn "Newton-Raphson did not converge!";
    return network_data
end

function test_runpf(network_data)
    ## This is the power flow solving test
    num_bus=length(network_data["bus"])
    vm=zeros(Float64, num_bus); [vm[bus["idx"]]=bus["vm"] for (i,bus) in network_data["bus"]]
    va=zeros(Float64, num_bus); [va[bus["idx"]]=bus["va"] for (i,bus) in network_data["bus"]]

    network_data_test=deepcopy(network_data)
    [gen["pg"]=gen["pg"]+gen["ptc_factor"]*network_data_test["delta"] for (i,gen) in network_data_test["gen"]]
    result_test=run_pf(network_data_test, ACPPowerModel, with_optimizer(Ipopt.Optimizer,print_level=0))
    vm_test=zeros(Float64, num_bus); [vm_test[bus["idx"]]=result_test["solution"]["bus"][i]["vm"] for (i,bus) in network_data_test["bus"]]
    va_test=zeros(Float64, num_bus); [va_test[bus["idx"]]=result_test["solution"]["bus"][i]["va"] for (i,bus) in network_data_test["bus"]]

    if max(maximum(abs.(vm-vm_test)),maximum(abs.(va-va_test)))>1e-5; @error "PowerFlow Failed" end;
end

function check_violation(network_data,suppress_warning=0)
    num_bus=length(network_data["bus"])
    num_line=length(network_data["branch"])
    num_gen=length(network_data["gen"])
    num_load=length(network_data["load"])

    idx_fr=zeros(Int64, num_line); [idx_fr[branch["index"]]=network_data["bus"][string(branch["f_bus"])]["idx"] for (i,branch) in network_data["branch"]]
    idx_to=zeros(Int64, num_line); [idx_to[branch["index"]]=network_data["bus"][string(branch["t_bus"])]["idx"] for (i,branch) in network_data["branch"]]
    idx_gen=zeros(Int64, num_gen); [idx_gen[gen["index"]]=network_data["bus"][string(gen["gen_bus"])]["idx"] for (i,gen) in network_data["gen"]]
    idx_load=zeros(Int64, num_load); [idx_load[load["index"]]=network_data["bus"][string(load["load_bus"])]["idx"] for (i,load) in network_data["load"]]
    Cg=sparse(idx_gen,(1:num_gen),ones(num_gen),num_bus,num_gen)
    Cl=sparse(idx_load,(1:num_load),ones(num_load),num_bus,num_load)

    v0=zeros(Float64, num_bus); [v0[bus["idx"]]=bus["vm"] for (i,bus) in network_data["bus"]]
    θ0=zeros(Float64, num_bus); [θ0[bus["idx"]]=bus["va"] for (i,bus) in network_data["bus"]]
    v_cplx0=v0.*exp.(1im*θ0)
    Y,Yf,Yt=makeYbus(network_data)

    ql0=zeros(Float64, num_load); [ql0[load["index"]]=load["status"].*load["qd"] for (i,load) in network_data["load"]]
    qg_inj0=imag(v_cplx0.*conj(Y*v_cplx0))+Cl*ql0
    s_line_fr0=v_cplx0[idx_fr].*conj(Yf*v_cplx0)
    s_line_to0=v_cplx0[idx_to].*conj(Yt*v_cplx0)

    qg_max=zeros(Float64, num_gen); [qg_max[gen["index"]]=gen["gen_status"]*gen["qmax"] for (i,gen) in network_data["gen"]]
    qg_min=zeros(Float64, num_gen); [qg_min[gen["index"]]=gen["gen_status"]*gen["qmin"] for (i,gen) in network_data["gen"]]

    qinj_max=Cg*qg_max
    qinj_min=Cg*qg_min

    limit_tolerence=1e-5#1e-3
    margin=0.1
    digits_spec=2
    violation_status=Dict("vmag"=>[],"pg"=>[],"qg"=>[],"line"=>[],"angle"=>[])

    if isnan(network_data["delta"]) && suppress_warning==0; @error "Base point did not converge!"; end

    for (i,bus) in network_data["bus"]
        if (bus["vm"]>bus["vmax"]+limit_tolerence) || (bus["vm"]<bus["vmin"]-limit_tolerence)
            append!(violation_status["vmag"],i)
            if suppress_warning==0; @warn string("Vm violation at bus ",i,": ",round(bus["vmin"],digits=digits_spec)," ≤ ",round(bus["vm"],digits=digits_spec)," ≤ ",round(bus["vmax"],digits=digits_spec)); end
            margin=min(margin,bus["vmax"]-bus["vm"],bus["vm"]-bus["vmin"])
        end
        if (qg_inj0[bus["idx"]]>qinj_max[bus["idx"]]+limit_tolerence) || (qg_inj0[bus["idx"]]<qinj_min[bus["idx"]]-limit_tolerence)
            append!(violation_status["qg"],i)
            if suppress_warning==0; @warn string("Qg violation at bus ",i,": ",round(qinj_min[bus["idx"]],digits=digits_spec)," ≤ ",round(qg_inj0[bus["idx"]],digits=digits_spec)," ≤ ",round(qinj_max[bus["idx"]],digits=digits_spec)); end
            margin=min(margin,qinj_max[bus["idx"]]-qg_inj0[bus["idx"]],qg_inj0[bus["idx"]]-qinj_min[bus["idx"]])
        end
    end

    for (i,gen) in network_data["gen"]
        p_dispatch=gen["pg"]+gen["ptc_factor"]*network_data["delta"]
        if (p_dispatch>gen["gen_status"]*gen["pmax"]+limit_tolerence) || (p_dispatch<gen["gen_status"]*gen["pmin"]-limit_tolerence)
            append!(violation_status["pg"],i)
            if suppress_warning==0; @warn string("Pg violation at gen ",i,": ",round(gen["pmin"],digits=digits_spec)," ≤ ",round(p_dispatch,digits=digits_spec)," ≤ ",round(gen["pmax"],digits=digits_spec)); end
            margin=min(margin,gen["pmax"]-p_dispatch,p_dispatch-gen["pmin"])
        end
    end

    for (i,branch) in network_data["branch"]
        angle_diff=network_data["bus"][string(branch["f_bus"])]["va"]-network_data["bus"][string(branch["t_bus"])]["va"]
        if (angle_diff>branch["angmax"]+limit_tolerence) || (angle_diff<branch["angmin"]-limit_tolerence)
            append!(violation_status["angle"],i)
            if suppress_warning==0; @warn string("Angle violation at line ",i,": ",round(branch["angmin"],digits=digits_spec)," ≤ ",round(angle_diff,digits=digits_spec)," ≤ ",round(branch["angmax"],digits=digits_spec)); end
            margin=min(margin,branch["angmax"]-angle_diff,angle_diff-branch["angmin"])
        end
        if (abs(s_line_fr0[branch["index"]])>branch["rate_a"]+limit_tolerence) || (abs(s_line_to0[branch["index"]])>branch["rate_a"]+limit_tolerence)
            #println(abs(s_line_fr0[branch["index"]]),"    ",branch["rate_a"])
            append!(violation_status["line"],i)
            if suppress_warning==0; @warn string("Line Flow violation at line ",i,": ",round(abs(s_line_fr0[branch["index"]]),digits=digits_spec),", ",round(abs(s_line_to0[branch["index"]]),digits=digits_spec)," ≤ ",round(branch["rate_a"],digits=digits_spec)); end
            margin=min(margin,branch["rate_a"]-abs(s_line_fr0[branch["index"]]),branch["rate_a"]-abs(s_line_to0[branch["index"]]))
        end
    end
    return violation_status, margin
end

function opf_initialization(network_data, restriction_level)
    network_data_opf=deepcopy(network_data)

    if restriction_level=="uniform"
        [gen["cost"]=[1.,0.] for (i,gen) in network_data_opf["gen"]]
        [gen["ncose"]=2 for (i,gen) in network_data_opf["gen"]]
        restriction_level=0
    end

    for (i,bus) in network_data_opf["bus"]
        bus["vmax"]=bus["vmax"]-restriction_level*(bus["vmax"]-bus["vmin"])
        bus["vmin"]=bus["vmin"]+restriction_level*(bus["vmax"]-bus["vmin"])
    end

    for (i,gen) in network_data_opf["gen"]
        gen["pmax"]=gen["gen_status"]*gen["pmax"]-restriction_level*(gen["gen_status"]*gen["pmax"]-gen["pmin"])
        gen["pmin"]=gen["gen_status"]*gen["pmin"]+restriction_level*(gen["gen_status"]*gen["pmax"]-gen["pmin"])
        gen["qmax"]=gen["gen_status"]*gen["qmax"]-restriction_level*(gen["gen_status"]*gen["qmax"]-gen["qmin"])
        gen["qmin"]=gen["gen_status"]*gen["qmin"]+restriction_level*(gen["gen_status"]*gen["qmax"]-gen["qmin"])
    end

    for (i,branch) in network_data_opf["branch"]
        branch["angmax"]=branch["angmax"]-restriction_level*(branch["angmax"]-branch["angmin"])
        branch["angmin"]=branch["angmin"]+restriction_level*(branch["angmax"]-branch["angmin"])
        branch["rate_a"]=branch["rate_a"]*(1-restriction_level)
        branch["rate_b"]=branch["rate_b"]*(1-restriction_level)
        branch["rate_c"]=branch["rate_c"]*(1-restriction_level)
    end

    if restriction_level==-1
        [gen["cost"]=[1.,0.] for (i,gen) in network_data_opf["gen"]]
        [gen["ncose"]=2 for (i,gen) in network_data_opf["gen"]]
    end

    result_opf=run_opf(network_data_opf, ACPPowerModel, with_optimizer(Ipopt.Optimizer,print_level=0))
    if (result_opf["termination_status"]!=JuMP.MathOptInterface.OPTIMAL) && (result_opf["termination_status"]!=JuMP.MathOptInterface.LOCALLY_SOLVED)
        @error string("OPF Initialization failed. (",result_opf["termination_status"],")")
    end
    PowerModels.update_data!(network_data, result_opf["solution"])
    for (i,gen) in network_data["gen"] # set gen to 0 for buses with gen_status off. (otherwise it's set to NaN)
        if isnan(gen["pg"])
            network_data["gen"][i]["pg"]=0
        end
    end
    [network_data["gen"][i]["vg"]=network_data["bus"][string(gen["gen_bus"])]["vm"] for (i,gen) in network_data["gen"]]

    index_num=1
    for (i,bus) in network_data["bus"]
        network_data["bus"][i]["idx"]=index_num
        #network_data["bus"][i]["idx"]=network_data["bus"][i]["index"]
        index_num=index_num+1
    end

    network_data["delta"]=0
    #num_gen=length(network_data["gen"])
    #[network_data["gen"][i]["ptc_factor"]=1/num_gen for (i,gen) in network_data["gen"]]
    [network_data["gen"][i]["ptc_factor"]=Int(network_data["bus"][string(gen["gen_bus"])]["bus_type"]==3) for (i,gen) in network_data["gen"]]

    num_load=length(network_data["load"])
    pl0=zeros(Float64, num_load); [pl0[load["index"]]=load["status"].*load["pd"] for (i,load) in network_data["load"]]
    Σ0=Matrix(Diagonal(pl0.^2))
    γ0=0.
    network_data["uncertainty"]=(Σ0,γ0)
    return network_data
end

function cvxrs(network_data,option,target_network_data=nothing,phase_shift=false)
    num_bus=length(network_data["bus"])
    num_gen=length(network_data["gen"])
    num_load=length(network_data["load"])
    num_line=length(network_data["branch"])

    bus_type=zeros(Int64, num_bus); [bus_type[bus["idx"]]=bus["bus_type"] for (i,bus) in network_data["bus"]]
    slack_bus=findall(bus_type.==3)[1]
    idx_nslack=setdiff(1:num_bus,slack_bus)

    idx_fr=zeros(Int64, num_line); [idx_fr[branch["index"]]=network_data["bus"][string(branch["f_bus"])]["idx"] for (i,branch) in network_data["branch"]]
    idx_to=zeros(Int64, num_line); [idx_to[branch["index"]]=network_data["bus"][string(branch["t_bus"])]["idx"] for (i,branch) in network_data["branch"]]
    idx_gen=zeros(Int64, num_gen); [idx_gen[gen["index"]]=network_data["bus"][string(gen["gen_bus"])]["idx"] for (i,gen) in network_data["gen"]]
    idx_load=zeros(Int64, num_load); [idx_load[load["index"]]=network_data["bus"][string(load["load_bus"])]["idx"] for (i,load) in network_data["load"]]
    idx_pq=findall(bus_type.==1)
    onoff_pq=zeros(Float64,num_bus); onoff_pq[idx_pq].=1.
    num_pq=sum(bus_type.==1)

    gencost=zeros(Float64,num_gen,3); [gencost[gen["index"],4-size(gen["cost"],1):end]=gen["gen_status"]*gen["cost"] for (i,gen) in network_data["gen"]]

    E_fr=sparse(idx_fr,1:num_line,ones(num_line),num_bus,num_line)
    E_to=sparse(idx_to,1:num_line,ones(num_line),num_bus,num_line)
    E=E_fr-E_to
    Cg=sparse(idx_gen,(1:num_gen),ones(num_gen),num_bus,num_gen)
    Cl=sparse(idx_load,(1:num_load),ones(num_load),num_bus,num_load)

    v0=zeros(Float64, num_bus); [v0[bus["idx"]]=bus["vm"] for (i,bus) in network_data["bus"]]
    θ0=zeros(Float64, num_bus); [θ0[bus["idx"]]=bus["va"] for (i,bus) in network_data["bus"]]#*pi/180
    φ0=E'*θ0
    alpha0=zeros(Float64, num_gen); [alpha0[gen["index"]]=gen["ptc_factor"] for (i,gen) in network_data["gen"]]
    delta0=network_data["delta"]

    gen_status=zeros(Int64, num_gen); [gen_status[gen["index"]]=gen["gen_status"] for (i,gen) in network_data["gen"]]
    load_status=zeros(Int64, num_load); [load_status[load["index"]]=load["status"] for (i,load) in network_data["load"]]
    pg0=zeros(Float64, num_gen); [pg0[gen["index"]]=gen["gen_status"].*gen["pg"] for (i,gen) in network_data["gen"]]
    qg0=zeros(Float64, num_gen); [qg0[gen["index"]]=gen["gen_status"].*gen["qg"] for (i,gen) in network_data["gen"]]
    pl0=zeros(Float64, num_load); [pl0[load["index"]]=load["status"].*load["pd"] for (i,load) in network_data["load"]]
    ql0=zeros(Float64, num_load); [ql0[load["index"]]=load["status"].*load["qd"] for (i,load) in network_data["load"]]
    p_inj0=Cg*(pg0+alpha0*delta0)-Cl*pl0
    q_inj0=Cg*qg0-Cl*ql0
    power_factor=ql0./(pl0.+1e-3)

    vmax=zeros(Float64, num_bus); [vmax[bus["idx"]]=bus["vmax"] for (i,bus) in network_data["bus"]]
    vmin=zeros(Float64, num_bus); [vmin[bus["idx"]]=bus["vmin"] for (i,bus) in network_data["bus"]]
    φmax=zeros(Float64, num_line); [φmax[branch["index"]]=branch["angmax"] for (i,branch) in network_data["branch"]]
    φmin=zeros(Float64, num_line); [φmin[branch["index"]]=branch["angmin"] for (i,branch) in network_data["branch"]]
    pg_max=zeros(Float64, num_gen); [pg_max[gen["index"]]=gen["gen_status"]*gen["pmax"] for (i,gen) in network_data["gen"]]
    pg_min=zeros(Float64, num_gen); [pg_min[gen["index"]]=gen["gen_status"]*gen["pmin"] for (i,gen) in network_data["gen"]]
    qg_max=zeros(Float64, num_gen); [qg_max[gen["index"]]=gen["gen_status"]*gen["qmax"] for (i,gen) in network_data["gen"]]
    qg_min=zeros(Float64, num_gen); [qg_min[gen["index"]]=gen["gen_status"]*gen["qmin"] for (i,gen) in network_data["gen"]]
    s_line_max=zeros(Float64, num_line); [s_line_max[branch["index"]]=branch["rate_a"] for (i,branch) in network_data["branch"]]

    (Σ0,γ0)=network_data["uncertainty"]

    J1=-Diagonal(v0[idx_fr].*v0[idx_to].*sin.(E'*θ0))*E'
    J2=Diagonal(v0[idx_to].*cos.(E'*θ0))*E_fr'+Diagonal(v0[idx_fr].*cos.(E'*θ0))*E_to'
    J3=Diagonal(v0[idx_fr].*v0[idx_to].*cos.(E'*θ0))*E'
    J4=Diagonal(v0[idx_to].*sin.(E'*θ0))*E_fr'+Diagonal(v0[idx_fr].*sin.(E'*θ0))*E_to'

    J_psi0=[zeros(2*num_bus,size(idx_nslack,1)+num_pq) [Cg*alpha0; zeros(num_bus)];
        J1[:,idx_nslack] J2[:,idx_pq] zeros(num_line);
        J3[:,idx_nslack] J4[:,idx_pq] zeros(num_line);
        zeros(num_bus,size(idx_nslack,1)) 2*Diagonal(v0)[:,idx_pq] zeros(num_bus)]

    g_vvcos0=v0[idx_fr].*v0[idx_to].*cos.(φ0)-onoff_pq[idx_fr].*v0[idx_fr].*v0[idx_to].*cos.(φ0)-v0[idx_fr].*onoff_pq[idx_to].*v0[idx_to].*cos.(φ0)+v0[idx_fr].*v0[idx_to].*sin.(φ0).*φ0
    g_vvsin0=v0[idx_fr].*v0[idx_to].*sin.(φ0)-onoff_pq[idx_fr].*v0[idx_fr].*v0[idx_to].*sin.(φ0)-v0[idx_fr].*onoff_pq[idx_to].*v0[idx_to].*sin.(φ0)-v0[idx_fr].*v0[idx_to].*cos.(φ0).*φ0
    g_vv0=v0.^2-2*onoff_pq.*v0.*v0
    ψ_vvcos0=v0[idx_fr].*v0[idx_to].*cos.(φ0)
    ψ_vvsin0=v0[idx_fr].*v0[idx_to].*sin.(φ0)
    ψ_vv0=v0.^2

    function Psi(v,θ,p_inj,q_inj)
        return [p_inj; q_inj; v[idx_fr].*v[idx_to].*cos.(E'*θ); v[idx_fr].*v[idx_to].*sin.(E'*θ); v.^2]
    end

    function gg(v,θ,p_inj,q_inj)
        φ=E'*θ
        g_vvcos=v[idx_fr].*v[idx_to].*cos.(φ)-onoff_pq[idx_fr].*v[idx_fr].*v0[idx_to].*cos.(φ0)-v0[idx_fr].*onoff_pq[idx_to].*v[idx_to].*cos.(φ0)+v0[idx_fr].*v0[idx_to].*sin.(φ0).*φ
        g_vvsin=v[idx_fr].*v[idx_to].*sin.(φ)-onoff_pq[idx_fr].*v[idx_fr].*v0[idx_to].*sin.(φ0)-v0[idx_fr].*onoff_pq[idx_to].*v[idx_to].*sin.(φ0)-v0[idx_fr].*v0[idx_to].*cos.(φ0).*φ
        g_vv=v.^2-2*onoff_pq.*v0.*v
        return [p_inj-Cg*alpha0*delta0; q_inj; g_vvcos; g_vvsin; g_vv]
    end

    Y,Yf,Yt,M,M_line=makeYbus(network_data,phase_shift)
    M_eq=M[[1:num_bus; num_bus.+idx_pq],:]
    idx_pvs=setdiff(1:num_bus,idx_pq)
    M_ineq=M[num_bus.+idx_pvs,:]

    J_inv_M=-Matrix(M_eq*J_psi0)\Matrix(M_eq)
    C=[E[idx_nslack,:]' zeros(num_line,num_pq+1);
        zeros(num_pq,size(idx_nslack,1)) Matrix{Float64}(I,num_pq,num_pq) zeros(num_pq);
        zeros(1,size(idx_nslack,1)+num_pq) 1]
    A=[Matrix{Float64}(I,num_line+num_pq+1,num_line+num_pq+1); -Matrix{Float64}(I,num_line+num_pq+1,num_line+num_pq+1)]*C
    K=A*J_inv_M

    θ_test=copy(θ0); v_test=copy(v0); delta_test=copy(delta0)
    θ_test[idx_nslack]=rand(size(idx_nslack,1)); v_test[idx_pq]=rand(num_pq);
    v_cpx=v_test.*cos.(θ_test)+1im*v_test.*sin.(θ_test)
    p_inj_test=real(v_cpx.*conj(Y*v_cpx))
    q_inj_test=imag(v_cpx.*conj(Y*v_cpx))

    if maximum(abs.(M*Psi(v0,θ0,p_inj0,q_inj0)))>1e-4; @error "Power flow equation is incorrect (nominal point)"; end
    if maximum(abs.(M*Psi(v_test,θ_test,p_inj_test,q_inj_test)))>1e-5; @error "Power flow equation is Incorrect (test)"; end
    if maximum(abs.(J_inv_M*gg(v0,θ0,p_inj0,q_inj0)-[θ0[idx_nslack];v0[idx_pq];delta0]))>1e-5; @error "Function g is Incorrect"; end
    if maximum(abs.(J_inv_M*gg(v0,θ0,p_inj0,-Cl*ql0)-[θ0[idx_nslack];v0[idx_pq];delta0]))>1e-5; @error "Function g without q_pv is Incorrect"; end
    if maximum(abs.(J_inv_M*gg(v_test,θ_test,p_inj_test,q_inj_test)-[θ_test[idx_nslack];v_test[idx_pq];delta_test]))>1e-5; @error "Function g is Incorrect"; end

    tol1=0
    K_plus=max.(K,zeros(size(K))).*(K.>tol1)
    K_minus=min.(K,zeros(size(K))).*(K.<-tol1)
    M_ineq_plus=max.(M_ineq,zeros(size(M_ineq))).*(M_ineq.>tol1)
    M_ineq_minus=min.(M_ineq,zeros(size(M_ineq))).*(M_ineq.<-tol1)
    M_line_plus=max.(M_line,zeros(size(M_line))).*(M_line.>tol1)
    M_line_minus=min.(M_line,zeros(size(M_line))).*(M_line.<-tol1);

    if OPT_version=="Gurobi"
        m_cvx = Model(with_optimizer(Gurobi.Optimizer,GRB_ENV,OutputFlag=0))
    elseif OPT_version=="Mosek"
        m_cvx = Model(with_optimizer(Mosek.Optimizer,LOG=0))
    end

    @variable(m_cvx, vmin[i]<=vᵘ[i=1:num_bus]<=vmax[i], start=v0[i])
    @variable(m_cvx, vmin[i]<=vˡ[i=1:num_bus]<=vmax[i], start=v0[i])
    @variable(m_cvx, vmin[i]-v0[i]<=Δvᵘ[i=1:num_bus]<=vmax[i]-v0[i], start=0)
    @variable(m_cvx, vmin[i]-v0[i]<=Δvˡ[i=1:num_bus]<=vmax[i]-v0[i], start=0)
    @variable(m_cvx, φmin[i]<=φᵘ[i=1:num_line]<=φmax[i], start=φ0[i])
    @variable(m_cvx, φmin[i]<=φˡ[i=1:num_line]<=φmax[i], start=φ0[i])
    @variable(m_cvx, φmin[i]-φ0[i]<=Δφᵘ[i=1:num_line]<=φmax[i]-φ0[i], start=0)
    @variable(m_cvx, φmin[i]-φ0[i]<=Δφˡ[i=1:num_line]<=φmax[i]-φ0[i], start=0)
    @variable(m_cvx, deltaᵘ, start=delta0)
    @variable(m_cvx, deltaˡ, start=delta0)
    @variable(m_cvx, Δdeltaᵘ, start=0)
    @variable(m_cvx, Δdeltaˡ, start=0)
    @variable(m_cvx, gᵘ_pinj[1:num_bus], start=0)
    @variable(m_cvx, gˡ_pinj[1:num_bus], start=0)
    @variable(m_cvx, gᵘ_qinj[1:num_bus], start=0)
    @variable(m_cvx, gˡ_qinj[1:num_bus], start=0)
    @variable(m_cvx, gᵘ_vvsin[1:num_line], start=0)
    @variable(m_cvx, gˡ_vvsin[1:num_line], start=0)
    @variable(m_cvx, gᵘ_vvcos[1:num_line], start=0)
    @variable(m_cvx, gˡ_vvcos[1:num_line], start=0)
    @variable(m_cvx, gᵘ_vv[1:num_bus], start=0)
    @variable(m_cvx, gˡ_vv[1:num_bus], start=0)

    @variable(m_cvx, ψᵘ_vvcos[1:num_line])#,start=ψ_vvcos0
    @variable(m_cvx, ψˡ_vvcos[1:num_line])#,start=ψ_vvcos0
    @variable(m_cvx, ψᵘ_vvsin[1:num_line])#,start=ψ_vvsin0
    @variable(m_cvx, ψˡ_vvsin[1:num_line])#,start=ψ_vvsin0
    @variable(m_cvx, ψᵘ_vv[i=1:num_bus])#,start=ψ_vv0
    @variable(m_cvx, ψˡ_vv[i=1:num_bus])#,start=ψ_vv0

    @variable(m_cvx, Δvvᵘ_uu[i=1:num_line])
    @variable(m_cvx, Δvvᵘ_ll[i=1:num_line])
    @variable(m_cvx, Δvvᵘ_ul[i=1:num_line])
    @variable(m_cvx, Δvvᵘ_lu[i=1:num_line])
    @variable(m_cvx, Δvvˡ_uu[i=1:num_line])
    @variable(m_cvx, Δvvˡ_ll[i=1:num_line])
    @variable(m_cvx, Δvvˡ_ul[i=1:num_line])
    @variable(m_cvx, Δvvˡ_lu[i=1:num_line])

    @variable(m_cvx, Δcosᵘ_u[i=1:num_line])
    @variable(m_cvx, Δcosᵘ_l[i=1:num_line])
    @variable(m_cvx, Δcosˡ_u[i=1:num_line])
    @variable(m_cvx, Δcosˡ_l[i=1:num_line])
    @variable(m_cvx, Δsinᵘ_u[i=1:num_line])
    @variable(m_cvx, Δsinᵘ_l[i=1:num_line])
    @variable(m_cvx, Δsinˡ_u[i=1:num_line])
    @variable(m_cvx, Δsinˡ_l[i=1:num_line])

    @variable(m_cvx, p_line_frᵘ[1:num_line])
    @variable(m_cvx, q_line_frᵘ[1:num_line])
    @variable(m_cvx, p_line_toᵘ[1:num_line])
    @variable(m_cvx, q_line_toᵘ[1:num_line])

    @variable(m_cvx, pg_opt[i=1:num_gen], start=pg0[i])
    @variable(m_cvx, pl_opt[i=1:num_load], start=pl0[i])
    @variable(m_cvx, ql_opt[i=1:num_load], start=ql0[i])
    @variable(m_cvx, pg_optᵘ[i=1:num_gen], start=pg0[i])
    @variable(m_cvx, alpha_opt[i=1:num_gen], start=alpha0[i])
    @variable(m_cvx, Δalpha_opt[i=1:num_gen], start=0)
    @variable(m_cvx, cost_opt, start=network_data["cost"])
    @variable(m_cvx, γ_opt ,start=0)

    @variable(m_cvx, slack_pg, start=0)
    @variable(m_cvx, slack_qg, start=0)
    @variable(m_cvx, slack_line, start=0)
    @variable(m_cvx, margin_opt, start=0)

    @constraint(m_cvx, vᵘ[idx_gen].==vˡ[idx_gen])

    @constraint(m_cvx, vᵘ.==v0+Δvᵘ)
    @constraint(m_cvx, vˡ.==v0+Δvˡ)
    @constraint(m_cvx, φᵘ.==φ0+Δφᵘ)
    @constraint(m_cvx, φˡ.==φ0+Δφˡ)
    #@constraint(m_cvx, vᵘ.>=vˡ)
    #@constraint(m_cvx, φᵘ.>=φˡ)
    @constraint(m_cvx, deltaᵘ.==delta0+Δdeltaᵘ)
    @constraint(m_cvx, deltaˡ.==delta0+Δdeltaˡ)
    @constraint(m_cvx, alpha_opt.==alpha0+Δalpha_opt)

    @constraint(m_cvx, Δvvᵘ_uu.>=Δvᵘ[idx_fr].*v0[idx_to]+v0[idx_fr].*Δvᵘ[idx_to]+1/4*(Δvᵘ[idx_fr]+Δvᵘ[idx_to]).^2)
    @constraint(m_cvx, Δvvᵘ_ll.>=Δvˡ[idx_fr].*v0[idx_to]+v0[idx_fr].*Δvˡ[idx_to]+1/4*(Δvˡ[idx_fr]+Δvˡ[idx_to]).^2)
    @constraint(m_cvx, Δvvᵘ_ul.>=Δvᵘ[idx_fr].*v0[idx_to]+v0[idx_fr].*Δvˡ[idx_to]+1/4*(Δvᵘ[idx_fr]+Δvˡ[idx_to]).^2)
    @constraint(m_cvx, Δvvᵘ_lu.>=Δvˡ[idx_fr].*v0[idx_to]+v0[idx_fr].*Δvᵘ[idx_to]+1/4*(Δvˡ[idx_fr]+Δvᵘ[idx_to]).^2)

    @constraint(m_cvx, Δvvˡ_uu.<=Δvᵘ[idx_fr].*v0[idx_to]+v0[idx_fr].*Δvᵘ[idx_to]-1/4*(Δvᵘ[idx_fr]-Δvᵘ[idx_to]).^2)
    @constraint(m_cvx, Δvvˡ_ll.<=Δvˡ[idx_fr].*v0[idx_to]+v0[idx_fr].*Δvˡ[idx_to]-1/4*(Δvˡ[idx_fr]-Δvˡ[idx_to]).^2)
    @constraint(m_cvx, Δvvˡ_ul.<=Δvᵘ[idx_fr].*v0[idx_to]+v0[idx_fr].*Δvˡ[idx_to]-1/4*(Δvᵘ[idx_fr]-Δvˡ[idx_to]).^2)
    @constraint(m_cvx, Δvvˡ_lu.<=Δvˡ[idx_fr].*v0[idx_to]+v0[idx_fr].*Δvᵘ[idx_to]-1/4*(Δvˡ[idx_fr]-Δvᵘ[idx_to]).^2)

    @constraint(m_cvx, Δcosᵘ_u.>=-sin.(φ0).*Δφᵘ+1/2*(Δφᵘ).^2)
    @constraint(m_cvx, Δcosᵘ_l.>=-sin.(φ0).*Δφˡ+1/2*(Δφˡ).^2)
    @constraint(m_cvx, Δcosˡ_u.<=-sin.(φ0).*Δφᵘ-1/2*(Δφᵘ).^2)
    @constraint(m_cvx, Δcosˡ_l.<=-sin.(φ0).*Δφˡ-1/2*(Δφˡ).^2)

    @constraint(m_cvx, Δsinᵘ_u.>=cos.(φ0).*Δφᵘ+1/2*(Δφᵘ).^2)
    @constraint(m_cvx, Δsinᵘ_l.>=cos.(φ0).*Δφˡ+1/2*(Δφˡ).^2)
    @constraint(m_cvx, Δsinˡ_u.<=cos.(φ0).*Δφᵘ-1/2*(Δφᵘ).^2)
    @constraint(m_cvx, Δsinˡ_l.<=cos.(φ0).*Δφˡ-1/2*(Δφˡ).^2)

    @constraint(m_cvx, gᵘ_vvcos.>=ψ_vvcos0+cos.(φ0).*Δvvᵘ_uu+v0[idx_fr].*v0[idx_to].*Δcosᵘ_u+1/4*(Δvvᵘ_uu+Δcosᵘ_u).^2-onoff_pq[idx_fr].*vᵘ[idx_fr].*v0[idx_to].*cos.(φ0)-v0[idx_fr].*onoff_pq[idx_to].*vᵘ[idx_to].*cos.(φ0)+v0[idx_fr].*v0[idx_to].*sin.(φ0).*φᵘ)
    @constraint(m_cvx, gᵘ_vvcos.>=ψ_vvcos0+cos.(φ0).*Δvvᵘ_uu+v0[idx_fr].*v0[idx_to].*Δcosᵘ_l+1/4*(Δvvᵘ_uu+Δcosᵘ_l).^2-onoff_pq[idx_fr].*vᵘ[idx_fr].*v0[idx_to].*cos.(φ0)-v0[idx_fr].*onoff_pq[idx_to].*vᵘ[idx_to].*cos.(φ0)+v0[idx_fr].*v0[idx_to].*sin.(φ0).*φˡ)
    @constraint(m_cvx, gᵘ_vvcos.>=ψ_vvcos0+cos.(φ0).*Δvvᵘ_ul+v0[idx_fr].*v0[idx_to].*Δcosᵘ_u+1/4*(Δvvᵘ_ul+Δcosᵘ_u).^2-onoff_pq[idx_fr].*vᵘ[idx_fr].*v0[idx_to].*cos.(φ0)-v0[idx_fr].*onoff_pq[idx_to].*vˡ[idx_to].*cos.(φ0)+v0[idx_fr].*v0[idx_to].*sin.(φ0).*φᵘ)
    @constraint(m_cvx, gᵘ_vvcos.>=ψ_vvcos0+cos.(φ0).*Δvvᵘ_ul+v0[idx_fr].*v0[idx_to].*Δcosᵘ_l+1/4*(Δvvᵘ_ul+Δcosᵘ_l).^2-onoff_pq[idx_fr].*vᵘ[idx_fr].*v0[idx_to].*cos.(φ0)-v0[idx_fr].*onoff_pq[idx_to].*vˡ[idx_to].*cos.(φ0)+v0[idx_fr].*v0[idx_to].*sin.(φ0).*φˡ)
    @constraint(m_cvx, gᵘ_vvcos.>=ψ_vvcos0+cos.(φ0).*Δvvᵘ_lu+v0[idx_fr].*v0[idx_to].*Δcosᵘ_u+1/4*(Δvvᵘ_lu+Δcosᵘ_u).^2-onoff_pq[idx_fr].*vˡ[idx_fr].*v0[idx_to].*cos.(φ0)-v0[idx_fr].*onoff_pq[idx_to].*vᵘ[idx_to].*cos.(φ0)+v0[idx_fr].*v0[idx_to].*sin.(φ0).*φᵘ)
    @constraint(m_cvx, gᵘ_vvcos.>=ψ_vvcos0+cos.(φ0).*Δvvᵘ_lu+v0[idx_fr].*v0[idx_to].*Δcosᵘ_l+1/4*(Δvvᵘ_lu+Δcosᵘ_l).^2-onoff_pq[idx_fr].*vˡ[idx_fr].*v0[idx_to].*cos.(φ0)-v0[idx_fr].*onoff_pq[idx_to].*vᵘ[idx_to].*cos.(φ0)+v0[idx_fr].*v0[idx_to].*sin.(φ0).*φˡ)
    @constraint(m_cvx, gᵘ_vvcos.>=ψ_vvcos0+cos.(φ0).*Δvvᵘ_ll+v0[idx_fr].*v0[idx_to].*Δcosᵘ_u+1/4*(Δvvᵘ_ll+Δcosᵘ_u).^2-onoff_pq[idx_fr].*vˡ[idx_fr].*v0[idx_to].*cos.(φ0)-v0[idx_fr].*onoff_pq[idx_to].*vˡ[idx_to].*cos.(φ0)+v0[idx_fr].*v0[idx_to].*sin.(φ0).*φᵘ)
    @constraint(m_cvx, gᵘ_vvcos.>=ψ_vvcos0+cos.(φ0).*Δvvᵘ_ll+v0[idx_fr].*v0[idx_to].*Δcosᵘ_l+1/4*(Δvvᵘ_ll+Δcosᵘ_l).^2-onoff_pq[idx_fr].*vˡ[idx_fr].*v0[idx_to].*cos.(φ0)-v0[idx_fr].*onoff_pq[idx_to].*vˡ[idx_to].*cos.(φ0)+v0[idx_fr].*v0[idx_to].*sin.(φ0).*φˡ)

    @constraint(m_cvx, gˡ_vvcos.<=ψ_vvcos0+cos.(φ0).*Δvvˡ_uu+v0[idx_fr].*v0[idx_to].*Δcosˡ_u-1/4*(Δvvˡ_uu-Δcosˡ_u).^2-onoff_pq[idx_fr].*vᵘ[idx_fr].*v0[idx_to].*cos.(φ0)-v0[idx_fr].*onoff_pq[idx_to].*vᵘ[idx_to].*cos.(φ0)+v0[idx_fr].*v0[idx_to].*sin.(φ0).*φᵘ)
    @constraint(m_cvx, gˡ_vvcos.<=ψ_vvcos0+cos.(φ0).*Δvvˡ_uu+v0[idx_fr].*v0[idx_to].*Δcosˡ_l-1/4*(Δvvˡ_uu-Δcosˡ_l).^2-onoff_pq[idx_fr].*vᵘ[idx_fr].*v0[idx_to].*cos.(φ0)-v0[idx_fr].*onoff_pq[idx_to].*vᵘ[idx_to].*cos.(φ0)+v0[idx_fr].*v0[idx_to].*sin.(φ0).*φˡ)
    @constraint(m_cvx, gˡ_vvcos.<=ψ_vvcos0+cos.(φ0).*Δvvˡ_ul+v0[idx_fr].*v0[idx_to].*Δcosˡ_u-1/4*(Δvvˡ_ul-Δcosˡ_u).^2-onoff_pq[idx_fr].*vᵘ[idx_fr].*v0[idx_to].*cos.(φ0)-v0[idx_fr].*onoff_pq[idx_to].*vˡ[idx_to].*cos.(φ0)+v0[idx_fr].*v0[idx_to].*sin.(φ0).*φᵘ)
    @constraint(m_cvx, gˡ_vvcos.<=ψ_vvcos0+cos.(φ0).*Δvvˡ_ul+v0[idx_fr].*v0[idx_to].*Δcosˡ_l-1/4*(Δvvˡ_ul-Δcosˡ_l).^2-onoff_pq[idx_fr].*vᵘ[idx_fr].*v0[idx_to].*cos.(φ0)-v0[idx_fr].*onoff_pq[idx_to].*vˡ[idx_to].*cos.(φ0)+v0[idx_fr].*v0[idx_to].*sin.(φ0).*φˡ)
    @constraint(m_cvx, gˡ_vvcos.<=ψ_vvcos0+cos.(φ0).*Δvvˡ_lu+v0[idx_fr].*v0[idx_to].*Δcosˡ_u-1/4*(Δvvˡ_lu-Δcosˡ_u).^2-onoff_pq[idx_fr].*vˡ[idx_fr].*v0[idx_to].*cos.(φ0)-v0[idx_fr].*onoff_pq[idx_to].*vᵘ[idx_to].*cos.(φ0)+v0[idx_fr].*v0[idx_to].*sin.(φ0).*φᵘ)
    @constraint(m_cvx, gˡ_vvcos.<=ψ_vvcos0+cos.(φ0).*Δvvˡ_lu+v0[idx_fr].*v0[idx_to].*Δcosˡ_l-1/4*(Δvvˡ_lu-Δcosˡ_l).^2-onoff_pq[idx_fr].*vˡ[idx_fr].*v0[idx_to].*cos.(φ0)-v0[idx_fr].*onoff_pq[idx_to].*vᵘ[idx_to].*cos.(φ0)+v0[idx_fr].*v0[idx_to].*sin.(φ0).*φˡ)
    @constraint(m_cvx, gˡ_vvcos.<=ψ_vvcos0+cos.(φ0).*Δvvˡ_ll+v0[idx_fr].*v0[idx_to].*Δcosˡ_u-1/4*(Δvvˡ_ll-Δcosˡ_u).^2-onoff_pq[idx_fr].*vˡ[idx_fr].*v0[idx_to].*cos.(φ0)-v0[idx_fr].*onoff_pq[idx_to].*vˡ[idx_to].*cos.(φ0)+v0[idx_fr].*v0[idx_to].*sin.(φ0).*φᵘ)
    @constraint(m_cvx, gˡ_vvcos.<=ψ_vvcos0+cos.(φ0).*Δvvˡ_ll+v0[idx_fr].*v0[idx_to].*Δcosˡ_l-1/4*(Δvvˡ_ll-Δcosˡ_l).^2-onoff_pq[idx_fr].*vˡ[idx_fr].*v0[idx_to].*cos.(φ0)-v0[idx_fr].*onoff_pq[idx_to].*vˡ[idx_to].*cos.(φ0)+v0[idx_fr].*v0[idx_to].*sin.(φ0).*φˡ)


    @constraint(m_cvx, gᵘ_vvsin.>=ψ_vvsin0+sin.(φ0).*Δvvᵘ_uu+v0[idx_fr].*v0[idx_to].*Δsinᵘ_u+1/4*(Δvvᵘ_uu+Δsinᵘ_u).^2-onoff_pq[idx_fr].*vᵘ[idx_fr].*v0[idx_to].*sin.(φ0)-v0[idx_fr].*onoff_pq[idx_to].*vᵘ[idx_to].*sin.(φ0)-v0[idx_fr].*v0[idx_to].*cos.(φ0).*φᵘ)
    @constraint(m_cvx, gᵘ_vvsin.>=ψ_vvsin0+sin.(φ0).*Δvvᵘ_uu+v0[idx_fr].*v0[idx_to].*Δsinᵘ_l+1/4*(Δvvᵘ_uu+Δsinᵘ_l).^2-onoff_pq[idx_fr].*vᵘ[idx_fr].*v0[idx_to].*sin.(φ0)-v0[idx_fr].*onoff_pq[idx_to].*vᵘ[idx_to].*sin.(φ0)-v0[idx_fr].*v0[idx_to].*cos.(φ0).*φˡ)
    @constraint(m_cvx, gᵘ_vvsin.>=ψ_vvsin0+sin.(φ0).*Δvvᵘ_ul+v0[idx_fr].*v0[idx_to].*Δsinᵘ_u+1/4*(Δvvᵘ_ul+Δsinᵘ_u).^2-onoff_pq[idx_fr].*vᵘ[idx_fr].*v0[idx_to].*sin.(φ0)-v0[idx_fr].*onoff_pq[idx_to].*vˡ[idx_to].*sin.(φ0)-v0[idx_fr].*v0[idx_to].*cos.(φ0).*φᵘ)
    @constraint(m_cvx, gᵘ_vvsin.>=ψ_vvsin0+sin.(φ0).*Δvvᵘ_ul+v0[idx_fr].*v0[idx_to].*Δsinᵘ_l+1/4*(Δvvᵘ_ul+Δsinᵘ_l).^2-onoff_pq[idx_fr].*vᵘ[idx_fr].*v0[idx_to].*sin.(φ0)-v0[idx_fr].*onoff_pq[idx_to].*vˡ[idx_to].*sin.(φ0)-v0[idx_fr].*v0[idx_to].*cos.(φ0).*φˡ)
    @constraint(m_cvx, gᵘ_vvsin.>=ψ_vvsin0+sin.(φ0).*Δvvᵘ_lu+v0[idx_fr].*v0[idx_to].*Δsinᵘ_u+1/4*(Δvvᵘ_lu+Δsinᵘ_u).^2-onoff_pq[idx_fr].*vˡ[idx_fr].*v0[idx_to].*sin.(φ0)-v0[idx_fr].*onoff_pq[idx_to].*vᵘ[idx_to].*sin.(φ0)-v0[idx_fr].*v0[idx_to].*cos.(φ0).*φᵘ)
    @constraint(m_cvx, gᵘ_vvsin.>=ψ_vvsin0+sin.(φ0).*Δvvᵘ_lu+v0[idx_fr].*v0[idx_to].*Δsinᵘ_l+1/4*(Δvvᵘ_lu+Δsinᵘ_l).^2-onoff_pq[idx_fr].*vˡ[idx_fr].*v0[idx_to].*sin.(φ0)-v0[idx_fr].*onoff_pq[idx_to].*vᵘ[idx_to].*sin.(φ0)-v0[idx_fr].*v0[idx_to].*cos.(φ0).*φˡ)
    @constraint(m_cvx, gᵘ_vvsin.>=ψ_vvsin0+sin.(φ0).*Δvvᵘ_ll+v0[idx_fr].*v0[idx_to].*Δsinᵘ_u+1/4*(Δvvᵘ_ll+Δsinᵘ_u).^2-onoff_pq[idx_fr].*vˡ[idx_fr].*v0[idx_to].*sin.(φ0)-v0[idx_fr].*onoff_pq[idx_to].*vˡ[idx_to].*sin.(φ0)-v0[idx_fr].*v0[idx_to].*cos.(φ0).*φᵘ)
    @constraint(m_cvx, gᵘ_vvsin.>=ψ_vvsin0+sin.(φ0).*Δvvᵘ_ll+v0[idx_fr].*v0[idx_to].*Δsinᵘ_l+1/4*(Δvvᵘ_ll+Δsinᵘ_l).^2-onoff_pq[idx_fr].*vˡ[idx_fr].*v0[idx_to].*sin.(φ0)-v0[idx_fr].*onoff_pq[idx_to].*vˡ[idx_to].*sin.(φ0)-v0[idx_fr].*v0[idx_to].*cos.(φ0).*φˡ)

    @constraint(m_cvx, gᵘ_vvsin.>=ψ_vvsin0+sin.(φ0).*Δvvˡ_uu+v0[idx_fr].*v0[idx_to].*Δsinᵘ_u+1/4*(Δvvˡ_uu+Δsinᵘ_u).^2-onoff_pq[idx_fr].*vᵘ[idx_fr].*v0[idx_to].*sin.(φ0)-v0[idx_fr].*onoff_pq[idx_to].*vᵘ[idx_to].*sin.(φ0)-v0[idx_fr].*v0[idx_to].*cos.(φ0).*φᵘ)
    @constraint(m_cvx, gᵘ_vvsin.>=ψ_vvsin0+sin.(φ0).*Δvvˡ_uu+v0[idx_fr].*v0[idx_to].*Δsinᵘ_l+1/4*(Δvvˡ_uu+Δsinᵘ_l).^2-onoff_pq[idx_fr].*vᵘ[idx_fr].*v0[idx_to].*sin.(φ0)-v0[idx_fr].*onoff_pq[idx_to].*vᵘ[idx_to].*sin.(φ0)-v0[idx_fr].*v0[idx_to].*cos.(φ0).*φˡ)
    @constraint(m_cvx, gᵘ_vvsin.>=ψ_vvsin0+sin.(φ0).*Δvvˡ_ul+v0[idx_fr].*v0[idx_to].*Δsinᵘ_u+1/4*(Δvvˡ_ul+Δsinᵘ_u).^2-onoff_pq[idx_fr].*vᵘ[idx_fr].*v0[idx_to].*sin.(φ0)-v0[idx_fr].*onoff_pq[idx_to].*vˡ[idx_to].*sin.(φ0)-v0[idx_fr].*v0[idx_to].*cos.(φ0).*φᵘ)
    @constraint(m_cvx, gᵘ_vvsin.>=ψ_vvsin0+sin.(φ0).*Δvvˡ_ul+v0[idx_fr].*v0[idx_to].*Δsinᵘ_l+1/4*(Δvvˡ_ul+Δsinᵘ_l).^2-onoff_pq[idx_fr].*vᵘ[idx_fr].*v0[idx_to].*sin.(φ0)-v0[idx_fr].*onoff_pq[idx_to].*vˡ[idx_to].*sin.(φ0)-v0[idx_fr].*v0[idx_to].*cos.(φ0).*φˡ)
    @constraint(m_cvx, gᵘ_vvsin.>=ψ_vvsin0+sin.(φ0).*Δvvˡ_lu+v0[idx_fr].*v0[idx_to].*Δsinᵘ_u+1/4*(Δvvˡ_lu+Δsinᵘ_u).^2-onoff_pq[idx_fr].*vˡ[idx_fr].*v0[idx_to].*sin.(φ0)-v0[idx_fr].*onoff_pq[idx_to].*vᵘ[idx_to].*sin.(φ0)-v0[idx_fr].*v0[idx_to].*cos.(φ0).*φᵘ)
    @constraint(m_cvx, gᵘ_vvsin.>=ψ_vvsin0+sin.(φ0).*Δvvˡ_lu+v0[idx_fr].*v0[idx_to].*Δsinᵘ_l+1/4*(Δvvˡ_lu+Δsinᵘ_l).^2-onoff_pq[idx_fr].*vˡ[idx_fr].*v0[idx_to].*sin.(φ0)-v0[idx_fr].*onoff_pq[idx_to].*vᵘ[idx_to].*sin.(φ0)-v0[idx_fr].*v0[idx_to].*cos.(φ0).*φˡ)
    @constraint(m_cvx, gᵘ_vvsin.>=ψ_vvsin0+sin.(φ0).*Δvvˡ_ll+v0[idx_fr].*v0[idx_to].*Δsinᵘ_u+1/4*(Δvvˡ_ll+Δsinᵘ_u).^2-onoff_pq[idx_fr].*vˡ[idx_fr].*v0[idx_to].*sin.(φ0)-v0[idx_fr].*onoff_pq[idx_to].*vˡ[idx_to].*sin.(φ0)-v0[idx_fr].*v0[idx_to].*cos.(φ0).*φᵘ)
    @constraint(m_cvx, gᵘ_vvsin.>=ψ_vvsin0+sin.(φ0).*Δvvˡ_ll+v0[idx_fr].*v0[idx_to].*Δsinᵘ_l+1/4*(Δvvˡ_ll+Δsinᵘ_l).^2-onoff_pq[idx_fr].*vˡ[idx_fr].*v0[idx_to].*sin.(φ0)-v0[idx_fr].*onoff_pq[idx_to].*vˡ[idx_to].*sin.(φ0)-v0[idx_fr].*v0[idx_to].*cos.(φ0).*φˡ)

    @constraint(m_cvx, gˡ_vvsin.<=ψ_vvsin0+sin.(φ0).*Δvvˡ_uu+v0[idx_fr].*v0[idx_to].*Δsinˡ_u-1/4*(Δvvˡ_uu-Δsinˡ_u).^2-onoff_pq[idx_fr].*vᵘ[idx_fr].*v0[idx_to].*sin.(φ0)-v0[idx_fr].*onoff_pq[idx_to].*vᵘ[idx_to].*sin.(φ0)-v0[idx_fr].*v0[idx_to].*cos.(φ0).*φᵘ)
    @constraint(m_cvx, gˡ_vvsin.<=ψ_vvsin0+sin.(φ0).*Δvvˡ_uu+v0[idx_fr].*v0[idx_to].*Δsinˡ_l-1/4*(Δvvˡ_uu-Δsinˡ_l).^2-onoff_pq[idx_fr].*vᵘ[idx_fr].*v0[idx_to].*sin.(φ0)-v0[idx_fr].*onoff_pq[idx_to].*vᵘ[idx_to].*sin.(φ0)-v0[idx_fr].*v0[idx_to].*cos.(φ0).*φˡ)
    @constraint(m_cvx, gˡ_vvsin.<=ψ_vvsin0+sin.(φ0).*Δvvˡ_ul+v0[idx_fr].*v0[idx_to].*Δsinˡ_u-1/4*(Δvvˡ_ul-Δsinˡ_u).^2-onoff_pq[idx_fr].*vᵘ[idx_fr].*v0[idx_to].*sin.(φ0)-v0[idx_fr].*onoff_pq[idx_to].*vˡ[idx_to].*sin.(φ0)-v0[idx_fr].*v0[idx_to].*cos.(φ0).*φᵘ)
    @constraint(m_cvx, gˡ_vvsin.<=ψ_vvsin0+sin.(φ0).*Δvvˡ_ul+v0[idx_fr].*v0[idx_to].*Δsinˡ_l-1/4*(Δvvˡ_ul-Δsinˡ_l).^2-onoff_pq[idx_fr].*vᵘ[idx_fr].*v0[idx_to].*sin.(φ0)-v0[idx_fr].*onoff_pq[idx_to].*vˡ[idx_to].*sin.(φ0)-v0[idx_fr].*v0[idx_to].*cos.(φ0).*φˡ)
    @constraint(m_cvx, gˡ_vvsin.<=ψ_vvsin0+sin.(φ0).*Δvvˡ_lu+v0[idx_fr].*v0[idx_to].*Δsinˡ_u-1/4*(Δvvˡ_lu-Δsinˡ_u).^2-onoff_pq[idx_fr].*vˡ[idx_fr].*v0[idx_to].*sin.(φ0)-v0[idx_fr].*onoff_pq[idx_to].*vᵘ[idx_to].*sin.(φ0)-v0[idx_fr].*v0[idx_to].*cos.(φ0).*φᵘ)
    @constraint(m_cvx, gˡ_vvsin.<=ψ_vvsin0+sin.(φ0).*Δvvˡ_lu+v0[idx_fr].*v0[idx_to].*Δsinˡ_l-1/4*(Δvvˡ_lu-Δsinˡ_l).^2-onoff_pq[idx_fr].*vˡ[idx_fr].*v0[idx_to].*sin.(φ0)-v0[idx_fr].*onoff_pq[idx_to].*vᵘ[idx_to].*sin.(φ0)-v0[idx_fr].*v0[idx_to].*cos.(φ0).*φˡ)
    @constraint(m_cvx, gˡ_vvsin.<=ψ_vvsin0+sin.(φ0).*Δvvˡ_ll+v0[idx_fr].*v0[idx_to].*Δsinˡ_u-1/4*(Δvvˡ_ll-Δsinˡ_u).^2-onoff_pq[idx_fr].*vˡ[idx_fr].*v0[idx_to].*sin.(φ0)-v0[idx_fr].*onoff_pq[idx_to].*vˡ[idx_to].*sin.(φ0)-v0[idx_fr].*v0[idx_to].*cos.(φ0).*φᵘ)
    @constraint(m_cvx, gˡ_vvsin.<=ψ_vvsin0+sin.(φ0).*Δvvˡ_ll+v0[idx_fr].*v0[idx_to].*Δsinˡ_l-1/4*(Δvvˡ_ll-Δsinˡ_l).^2-onoff_pq[idx_fr].*vˡ[idx_fr].*v0[idx_to].*sin.(φ0)-v0[idx_fr].*onoff_pq[idx_to].*vˡ[idx_to].*sin.(φ0)-v0[idx_fr].*v0[idx_to].*cos.(φ0).*φˡ)

    @constraint(m_cvx, gˡ_vvsin.<=ψ_vvsin0+sin.(φ0).*Δvvᵘ_uu+v0[idx_fr].*v0[idx_to].*Δsinˡ_u-1/4*(Δvvᵘ_uu-Δsinˡ_u).^2-onoff_pq[idx_fr].*vᵘ[idx_fr].*v0[idx_to].*sin.(φ0)-v0[idx_fr].*onoff_pq[idx_to].*vᵘ[idx_to].*sin.(φ0)-v0[idx_fr].*v0[idx_to].*cos.(φ0).*φᵘ)
    @constraint(m_cvx, gˡ_vvsin.<=ψ_vvsin0+sin.(φ0).*Δvvᵘ_uu+v0[idx_fr].*v0[idx_to].*Δsinˡ_l-1/4*(Δvvᵘ_uu-Δsinˡ_l).^2-onoff_pq[idx_fr].*vᵘ[idx_fr].*v0[idx_to].*sin.(φ0)-v0[idx_fr].*onoff_pq[idx_to].*vᵘ[idx_to].*sin.(φ0)-v0[idx_fr].*v0[idx_to].*cos.(φ0).*φˡ)
    @constraint(m_cvx, gˡ_vvsin.<=ψ_vvsin0+sin.(φ0).*Δvvᵘ_ul+v0[idx_fr].*v0[idx_to].*Δsinˡ_u-1/4*(Δvvᵘ_ul-Δsinˡ_u).^2-onoff_pq[idx_fr].*vᵘ[idx_fr].*v0[idx_to].*sin.(φ0)-v0[idx_fr].*onoff_pq[idx_to].*vˡ[idx_to].*sin.(φ0)-v0[idx_fr].*v0[idx_to].*cos.(φ0).*φᵘ)
    @constraint(m_cvx, gˡ_vvsin.<=ψ_vvsin0+sin.(φ0).*Δvvᵘ_ul+v0[idx_fr].*v0[idx_to].*Δsinˡ_l-1/4*(Δvvᵘ_ul-Δsinˡ_l).^2-onoff_pq[idx_fr].*vᵘ[idx_fr].*v0[idx_to].*sin.(φ0)-v0[idx_fr].*onoff_pq[idx_to].*vˡ[idx_to].*sin.(φ0)-v0[idx_fr].*v0[idx_to].*cos.(φ0).*φˡ)
    @constraint(m_cvx, gˡ_vvsin.<=ψ_vvsin0+sin.(φ0).*Δvvᵘ_lu+v0[idx_fr].*v0[idx_to].*Δsinˡ_u-1/4*(Δvvᵘ_lu-Δsinˡ_u).^2-onoff_pq[idx_fr].*vˡ[idx_fr].*v0[idx_to].*sin.(φ0)-v0[idx_fr].*onoff_pq[idx_to].*vᵘ[idx_to].*sin.(φ0)-v0[idx_fr].*v0[idx_to].*cos.(φ0).*φᵘ)
    @constraint(m_cvx, gˡ_vvsin.<=ψ_vvsin0+sin.(φ0).*Δvvᵘ_lu+v0[idx_fr].*v0[idx_to].*Δsinˡ_l-1/4*(Δvvᵘ_lu-Δsinˡ_l).^2-onoff_pq[idx_fr].*vˡ[idx_fr].*v0[idx_to].*sin.(φ0)-v0[idx_fr].*onoff_pq[idx_to].*vᵘ[idx_to].*sin.(φ0)-v0[idx_fr].*v0[idx_to].*cos.(φ0).*φˡ)
    @constraint(m_cvx, gˡ_vvsin.<=ψ_vvsin0+sin.(φ0).*Δvvᵘ_ll+v0[idx_fr].*v0[idx_to].*Δsinˡ_u-1/4*(Δvvᵘ_ll-Δsinˡ_u).^2-onoff_pq[idx_fr].*vˡ[idx_fr].*v0[idx_to].*sin.(φ0)-v0[idx_fr].*onoff_pq[idx_to].*vˡ[idx_to].*sin.(φ0)-v0[idx_fr].*v0[idx_to].*cos.(φ0).*φᵘ)
    @constraint(m_cvx, gˡ_vvsin.<=ψ_vvsin0+sin.(φ0).*Δvvᵘ_ll+v0[idx_fr].*v0[idx_to].*Δsinˡ_l-1/4*(Δvvᵘ_ll-Δsinˡ_l).^2-onoff_pq[idx_fr].*vˡ[idx_fr].*v0[idx_to].*sin.(φ0)-v0[idx_fr].*onoff_pq[idx_to].*vˡ[idx_to].*sin.(φ0)-v0[idx_fr].*v0[idx_to].*cos.(φ0).*φˡ)

    @constraint(m_cvx, gᵘ_vv.>=vᵘ.^2-2*v0.*onoff_pq.*vᵘ)
    @constraint(m_cvx, gᵘ_vv.>=vˡ.^2-2*v0.*onoff_pq.*vˡ)
    #@constraint(m_cvx, gˡ_vv.<=ψ_vv0+2*v0.*Δvᵘ-2*v0.*onoff_pq.*vᵘ)
    #@constraint(m_cvx, gˡ_vv.<=ψ_vv0+2*v0.*Δvˡ-2*v0.*onoff_pq.*vˡ)
    @constraint(m_cvx, gˡ_vv.<=ψ_vv0+2*v0.*(vᵘ-v0)-2*v0.*onoff_pq.*vᵘ)
    @constraint(m_cvx, gˡ_vv.<=ψ_vv0+2*v0.*(vˡ-v0)-2*v0.*onoff_pq.*vˡ)

    @constraint(m_cvx, ψᵘ_vvcos.>=ψ_vvcos0+cos.(φ0).*Δvvᵘ_uu+v0[idx_fr].*v0[idx_to].*Δcosᵘ_u+1/4*(Δvvᵘ_uu+Δcosᵘ_u).^2)
    @constraint(m_cvx, ψᵘ_vvcos.>=ψ_vvcos0+cos.(φ0).*Δvvᵘ_uu+v0[idx_fr].*v0[idx_to].*Δcosᵘ_l+1/4*(Δvvᵘ_uu+Δcosᵘ_l).^2)
    @constraint(m_cvx, ψˡ_vvcos.<=ψ_vvcos0+cos.(φ0).*Δvvˡ_ll+v0[idx_fr].*v0[idx_to].*Δcosˡ_u-1/4*(Δvvˡ_ll-Δcosˡ_u).^2)
    @constraint(m_cvx, ψˡ_vvcos.<=ψ_vvcos0+cos.(φ0).*Δvvˡ_ll+v0[idx_fr].*v0[idx_to].*Δcosˡ_l-1/4*(Δvvˡ_ll-Δcosˡ_l).^2)

    @constraint(m_cvx, ψᵘ_vvsin.>=ψ_vvsin0+sin.(φ0).*Δvvᵘ_uu+v0[idx_fr].*v0[idx_to].*Δsinᵘ_u+1/4*(Δvvᵘ_uu+Δsinᵘ_u).^2)
    @constraint(m_cvx, ψᵘ_vvsin.>=ψ_vvsin0+sin.(φ0).*Δvvˡ_uu+v0[idx_fr].*v0[idx_to].*Δsinᵘ_u+1/4*(Δvvˡ_uu+Δsinᵘ_u).^2)
    @constraint(m_cvx, ψˡ_vvsin.<=ψ_vvsin0+sin.(φ0).*Δvvˡ_ll+v0[idx_fr].*v0[idx_to].*Δsinˡ_l-1/4*(Δvvˡ_ll-Δsinˡ_l).^2)
    @constraint(m_cvx, ψˡ_vvsin.<=ψ_vvsin0+sin.(φ0).*Δvvᵘ_uu+v0[idx_fr].*v0[idx_to].*Δsinˡ_l-1/4*(Δvvᵘ_uu-Δsinˡ_l).^2)

    @constraint(m_cvx, ψᵘ_vv.>=vᵘ.^2)
    @constraint(m_cvx, ψᵘ_vv.>=vˡ.^2)
    @constraint(m_cvx, ψˡ_vv.<=ψ_vv0+2*v0.*Δvᵘ)
    @constraint(m_cvx, ψˡ_vv.<=ψ_vv0+2*v0.*Δvˡ)

    @constraint(m_cvx, gᵘ_pinj.>=Cg*pg_opt-Cl*pl_opt)
    @constraint(m_cvx, gˡ_pinj.<=Cg*pg_opt-Cl*pl_opt)
    @constraint(m_cvx, alpha_opt.==alpha0)

    @constraint(m_cvx,pg_opt+alpha0.*deltaᵘ+slack_pg.<=pg_max)
    @constraint(m_cvx,pg_opt+alpha0.*deltaˡ-slack_pg.>=pg_min)

    zeta=sqrt.(sum((Cl[idx_pvs,:]*Diagonal(power_factor)*Matrix(cholesky(Σ0).U)).^2,dims=2))
    @constraint(m_cvx,M_ineq_minus*[zeros(num_bus); Cg*qg_max-Cl*ql_opt; ψᵘ_vvcos; ψᵘ_vvsin; ψᵘ_vv]+M_ineq_plus*[zeros(num_bus); Cg*qg_max-Cl*ql0; ψˡ_vvcos; ψˡ_vvsin; ψˡ_vv]-slack_qg-zeta.*γ_opt.>=0)
    @constraint(m_cvx,M_ineq_plus*[zeros(num_bus); Cg*qg_min-Cl*ql_opt; ψᵘ_vvcos; ψᵘ_vvsin; ψᵘ_vv]+M_ineq_minus*[zeros(num_bus); Cg*qg_min-Cl*ql0; ψˡ_vvcos; ψˡ_vvsin; ψˡ_vv]+slack_qg+zeta.*γ_opt.<=0)

    @constraint(m_cvx,M_line_plus*[zeros(2*num_bus); ψᵘ_vvcos; ψᵘ_vvsin; ψᵘ_vv]+M_line_minus*[zeros(2*num_bus); ψˡ_vvcos; ψˡ_vvsin; ψˡ_vv].<=[p_line_frᵘ; p_line_toᵘ; q_line_frᵘ; q_line_toᵘ])
    @constraint(m_cvx,-M_line_minus*[zeros(2*num_bus); ψᵘ_vvcos; ψᵘ_vvsin; ψᵘ_vv]-M_line_plus*[zeros(2*num_bus); ψˡ_vvcos; ψˡ_vvsin; ψˡ_vv].<=[p_line_frᵘ; p_line_toᵘ; q_line_frᵘ; q_line_toᵘ])
    @constraint(m_cvx,p_line_frᵘ.^2+q_line_frᵘ.^2+slack_line.<=s_line_max.^2);
    @constraint(m_cvx,p_line_toᵘ.^2+q_line_toᵘ.^2+slack_line.<=s_line_max.^2);

    xi=sqrt.(sum((K[:,1:2*num_bus]*[Cl; Cl*Diagonal(power_factor)]*cholesky(Σ0).U).^2,dims=2))
    #print(power_factor)
    @constraint(m_cvx,K_plus*[gᵘ_pinj; -Cl*ql_opt; gᵘ_vvcos; gᵘ_vvsin; gᵘ_vv]+K_minus*[gˡ_pinj; -Cl*ql_opt; gˡ_vvcos; gˡ_vvsin; gˡ_vv]+xi.*γ_opt.<=[φᵘ; vᵘ[idx_pq]; deltaᵘ; -φˡ; -vˡ[idx_pq]; -deltaˡ])

    @constraint(m_cvx,slack_pg>=margin_opt)
    @constraint(m_cvx,slack_qg>=margin_opt)
    @constraint(m_cvx,slack_line>=margin_opt)
    @constraint(m_cvx,γ_opt>=margin_opt)

    @constraint(m_cvx,pg_opt+alpha0.*deltaᵘ.<=pg_optᵘ)
    #@constraint(m_cvx,pg_opt+alpha0.*delta0+Δalpha_opt.*delta0+alpha0.*Δdeltaᵘ+1/4*(Δdeltaᵘ.+Δalpha_opt).^2 .<=pg_optᵘ)
    @constraint(m_cvx,gencost[:,1]'*(pg_optᵘ).^2+gencost[:,2]'*pg_optᵘ+sum(gencost[:,3])<=cost_opt)

    if target_network_data!=nothing
        pg_target=zeros(Float64, num_gen); [pg_target[gen["index"]]=gen["gen_status"].*gen["pg"] for (i,gen) in target_network_data["gen"]]
        v_target=zeros(Float64, num_bus); [v_target[bus["idx"]]=bus["vm"] for (i,bus) in target_network_data["bus"]]
        pl_target=zeros(Float64, num_load); [pl_target[load["index"]]=load["status"].*load["pd"] for (i,load) in target_network_data["load"]]
        ql_target=zeros(Float64, num_load); [ql_target[load["index"]]=load["status"].*load["qd"] for (i,load) in target_network_data["load"]]
    end

    if option=="margin" && target_network_data==nothing
        @constraint(m_cvx,pg_opt.==pg0)
        @constraint(m_cvx,pl_opt.==pl0)
        @constraint(m_cvx,ql_opt.==ql0)
        @constraint(m_cvx,vᵘ[idx_gen].==v0[idx_gen])
        @constraint(m_cvx,margin_opt>=0)
        @objective(m_cvx, Max, γ_opt)
        optimize!(m_cvx)
        γ0=JuMP.value(γ_opt)
        network_data["uncertainty"]=(Σ0,γ0)
    elseif option=="margin" && target_network_data!=nothing
        @constraint(m_cvx,pg_opt.==pg_target)
        @constraint(m_cvx,pl_opt.==pl_target)
        @constraint(m_cvx,ql_opt.==ql_target)
        @constraint(m_cvx,vᵘ[idx_gen].==v_target[idx_gen])
        @objective(m_cvx, Max, margin_opt)
        optimize!(m_cvx)
        γ0=JuMP.value(margin_opt)
        network_data["uncertainty"]=(Σ0,γ0)
    elseif option=="obj"
        @constraint(m_cvx,pl_opt.==pl0)
        @constraint(m_cvx,ql_opt.==ql0)
        @constraint(m_cvx,margin_opt>=0)
        @constraint(m_cvx,γ_opt>=γ0)
        @objective(m_cvx, Min, cost_opt)
        optimize!(m_cvx)
        [network_data["gen"][i]["pg"]=JuMP.value(pg_opt[gen["index"]])  for (i,gen) in network_data["gen"]]
        [network_data["gen"][i]["ptc_factor"]=JuMP.value(alpha_opt[gen["index"]])  for (i,gen) in network_data["gen"]]
        [network_data["gen"][i]["vg"]=JuMP.value(vᵘ[network_data["bus"][string(gen["gen_bus"])]["idx"]])  for (i,gen) in network_data["gen"]]
    end

    sanity_check=Dict()
    sanity_check["v0"]=v0
    sanity_check["φ0"]=φ0
    sanity_check["delta0"]=delta0
    sanity_check["alpha0"]=alpha0
    sanity_check["deltaᵘ"]=JuMP.value(deltaᵘ)
    sanity_check["deltaˡ"]=JuMP.value(deltaˡ)
    sanity_check["vᵘ"]=JuMP.value.(vᵘ)
    sanity_check["vˡ"]=JuMP.value.(vˡ)
    sanity_check["φᵘ"]=JuMP.value.(φᵘ)
    sanity_check["φˡ"]=JuMP.value.(φˡ)
    sanity_check["ψᵘ_vvcos"]=JuMP.value.(ψᵘ_vvcos)
    sanity_check["ψˡ_vvcos"]=JuMP.value.(ψˡ_vvcos)
    sanity_check["gᵘ_vvcos"]=JuMP.value.(gᵘ_vvcos)
    sanity_check["gˡ_vvcos"]=JuMP.value.(gˡ_vvcos)
    sanity_check["ψᵘ_vvsin"]=JuMP.value.(ψᵘ_vvsin)
    sanity_check["ψˡ_vvsin"]=JuMP.value.(ψˡ_vvsin)
    sanity_check["gᵘ_vvsin"]=JuMP.value.(gᵘ_vvsin)
    sanity_check["gˡ_vvsin"]=JuMP.value.(gˡ_vvsin)
    sanity_check["ψᵘ_vv"]=JuMP.value.(ψᵘ_vv)
    sanity_check["ψˡ_vv"]=JuMP.value.(ψˡ_vv)
    sanity_check["gᵘ_vv"]=JuMP.value.(gᵘ_vv)
    sanity_check["gˡ_vv"]=JuMP.value.(gˡ_vv)
    sanity_check["obj"]=JuMP.value(cost_opt)
    sanity_check["solve_time"]=solve_time(m_cvx)
    sanity_check["status"]=termination_status(m_cvx)

    if termination_status(m_cvx)!=JuMP.MathOptInterface.OPTIMAL; @warn string("CVXRS was not optimal. (", termination_status(m_cvx),")"); end
    return network_data, sanity_check
end

function test_cvxrs(network_data,sanity_check)
    # Sanity check and debugging

    limit_tolerence=1e-5
    digit_spec=6

    if (network_data["delta"]>sanity_check["deltaᵘ"]+limit_tolerence) || (network_data["delta"]<sanity_check["deltaˡ"]-limit_tolerence)
        @warn string("delta out of range: ",round(sanity_check["deltaˡ"],digits=digit_spec)," ≤ ",round(network_data["delta"],digits=digit_spec)," ≤ ",round(sanity_check["deltaᵘ"],digits=digit_spec))
    end

    for (i,bus) in network_data["bus"]
        if (bus["vm"]>sanity_check["vᵘ"][bus["idx"]]+limit_tolerence) || (bus["vm"]<sanity_check["vˡ"][bus["idx"]]-limit_tolerence)
            @warn string("Vm out of range at bus ",bus["bus_i"],"/",bus["idx"],": ",round(sanity_check["vˡ"][bus["idx"]],digits=digit_spec)," ≤ ",round(bus["vm"],digits=digit_spec)," ≤ ",round(sanity_check["vᵘ"][bus["idx"]],digits=digit_spec))
        end

        pq_i=Int(network_data["bus"][i]["bus_type"]==1)
        ψ_vv_i=bus["vm"].^2
        g_vv_i=bus["vm"].^2-2*pq_i.*sanity_check["v0"][network_data["bus"][i]["idx"]].*bus["vm"]
        if (ψ_vv_i>sanity_check["ψᵘ_vv"][bus["idx"]]+limit_tolerence) || (ψ_vv_i<sanity_check["ψˡ_vv"][bus["idx"]]-limit_tolerence)
            @warn string("vv out of range at bus ",bus["bus_i"],": ",round(sanity_check["ψˡ_vv"][bus["idx"]],digits=digit_spec)," ≤ ",round(ψ_vv_i,digits=digit_spec)," ≤ ",round(sanity_check["ψᵘ_vv"][bus["idx"]],digits=digit_spec))
        end
        if (g_vv_i>sanity_check["gᵘ_vv"][bus["idx"]]+limit_tolerence) || (g_vv_i<sanity_check["gˡ_vv"][bus["idx"]]-limit_tolerence)
            @warn string("g_vv out of range at bus ",bus["bus_i"],": ",round(sanity_check["gˡ_vv"][bus["idx"]],digits=digit_spec)," ≤ ",round(g_vv_i,digits=digit_spec)," ≤ ",round(sanity_check["gᵘ_vv"][bus["idx"]],digits=digit_spec))
        end
    end

    for (i,branch) in network_data["branch"]
        angle_diff=network_data["bus"][string(branch["f_bus"])]["va"]-network_data["bus"][string(branch["t_bus"])]["va"]
        if (angle_diff>sanity_check["φᵘ"][branch["index"]]+limit_tolerence) || (angle_diff<sanity_check["φˡ"][branch["index"]]-limit_tolerence)
            @warn string("Angle out of range at line ",i,": ",round(sanity_check["φˡ"][branch["index"]],digits=digit_spec)," ≤ ",round(angle_diff,digits=digit_spec)," ≤ ",round(sanity_check["φᵘ"][branch["index"]],digits=digit_spec))
        end

        v_fr=network_data["bus"][string(branch["f_bus"])]["vm"]
        v_to=network_data["bus"][string(branch["t_bus"])]["vm"]
        pq_fr=Int(network_data["bus"][string(branch["f_bus"])]["bus_type"]==1)
        pq_to=Int(network_data["bus"][string(branch["t_bus"])]["bus_type"]==1)
        v0_fr=sanity_check["v0"][network_data["bus"][string(branch["f_bus"])]["idx"]]
        v0_to=sanity_check["v0"][network_data["bus"][string(branch["t_bus"])]["idx"]]

        #angle_diff=JuMP.value(φᵘ[branch["index"]])
        #v_fr=JuMP.value(vᵘ[branch["f_bus"]])
        #v_to=JuMP.value(vᵘ[branch["t_bus"]])
        ψ_vvcos_i=v_fr.*v_to.*cos.(angle_diff)
        ψ_vvsin_i=v_fr.*v_to.*sin.(angle_diff)
        g_vvcos_i=v_fr*v_to*cos(angle_diff)-pq_fr*v_fr*v0_to*cos(sanity_check["φ0"][branch["index"]])-v0_fr*pq_to*v_to*cos(sanity_check["φ0"][branch["index"]])+v0_fr*v0_to*sin(sanity_check["φ0"][branch["index"]])*angle_diff
        g_vvsin_i=v_fr*v_to*sin(angle_diff)-pq_fr*v_fr*v0_to*sin(sanity_check["φ0"][branch["index"]])-v0_fr*pq_to*v_to*sin(sanity_check["φ0"][branch["index"]])-v0_fr*v0_to*cos(sanity_check["φ0"][branch["index"]])*angle_diff

        if (ψ_vvcos_i>sanity_check["ψᵘ_vvcos"][branch["index"]]+limit_tolerence) || (ψ_vvcos_i<sanity_check["ψˡ_vvcos"][branch["index"]]-limit_tolerence)
            @warn string("vvcos out of range at line ",i,": ",round(sanity_check["ψˡ_vvcos"][branch["index"]],digits=digit_spec)," ≤ ",round(ψ_vvcos_i,digits=digit_spec)," ≤ ",round(sanity_check["ψᵘ_vvcos"][branch["index"]],digits=digit_spec))
        end

        if (g_vvcos_i>sanity_check["gᵘ_vvcos"][branch["index"]]+limit_tolerence) || (g_vvcos_i<sanity_check["gˡ_vvcos"][branch["index"]]-limit_tolerence)
            @warn string("g_vvcos out of range at line ",i,": ",round(sanity_check["gˡ_vvcos"][branch["index"]],digits=digit_spec)," ≤ ",round(g_vvcos_i,digits=digit_spec)," ≤ ",round(sanity_check["gᵘ_vvcos"][branch["index"]],digits=digit_spec))
        end

        if (ψ_vvsin_i>sanity_check["ψᵘ_vvsin"][branch["index"]]+limit_tolerence) || (ψ_vvsin_i<sanity_check["ψˡ_vvsin"][branch["index"]]-limit_tolerence)
            @warn string("vvsin out of range at line ",i,": ",round(sanity_check["ψˡ_vvsin"][branch["index"]],digits=digit_spec)," ≤ ",round(ψ_vvsin_i,digits=digit_spec)," ≤ ",round(sanity_check["ψᵘ_vvsin"][branch["index"]],digits=digit_spec))
        end

        if (g_vvsin_i>sanity_check["gᵘ_vvsin"][branch["index"]]+limit_tolerence) || (g_vvsin_i<sanity_check["gˡ_vvsin"][branch["index"]]-limit_tolerence)
            @warn string("g_vvsin out of range at line ",i,": ",round(sanity_check["gˡ_vvsin"][branch["index"]],digits=digit_spec)," ≤ ",round(g_vvsin_i,digits=digit_spec)," ≤ ",round(sanity_check["gᵘ_vvsin"][branch["index"]],digits=digit_spec))
        end
    end
end

function plot2D(network_data,plot_bus,plot_rng,option,resolution=20)
    if option=="pg"
        u_plot0=[network_data["gen"][plot_bus[1]]["pg"], network_data["gen"][plot_bus[2]]["pg"]]
    elseif option=="pd"
        u_plot0=[network_data["load"][plot_bus[1]]["pd"], network_data["load"][plot_bus[2]]["pd"]]
    end

    u1_plot = range(plot_rng[1],stop=plot_rng[2],length=resolution)
    u2_plot = range(plot_rng[3],stop=plot_rng[4],length=resolution)
    U1_plot = repeat(u1_plot',length(u2_plot),1)
    U2_plot = repeat(u2_plot,1,length(u1_plot))
    exact_plot=zeros(size(U1_plot))
    cvxrs_plot=zeros(size(U1_plot))
    (Σ0,γ_plot)=network_data["uncertainty"]

    num_gen=length(network_data["gen"])
    pg_plot=zeros(Float64, num_gen); [pg_plot[gen["index"]]=gen["gen_status"].*gen["pg"] for (i,gen) in network_data["gen"]]
    for i=1:size(U1_plot,1)
        println("progress: ", i,"/",size(U1_plot,1))
        for j=1:size(U1_plot,2)
            network_data_plot=deepcopy(network_data)
            if option=="pg"
                network_data_plot["gen"][plot_bus[1]]["pg"]=U1_plot[i,j]
                network_data_plot["gen"][plot_bus[2]]["pg"]=U2_plot[i,j]
            elseif option=="pd"
                network_data_plot["load"][plot_bus[1]]["pd"]=U1_plot[i,j]
                network_data_plot["load"][plot_bus[2]]["pd"]=U2_plot[i,j]
                pf1=network_data["load"][plot_bus[1]]["qd"]/network_data["load"][plot_bus[1]]["pd"]
                pf2=network_data["load"][plot_bus[2]]["qd"]/network_data["load"][plot_bus[2]]["pd"]
                network_data_plot["load"][plot_bus[1]]["qd"]=pf1*U1_plot[i,j]
                network_data_plot["load"][plot_bus[2]]["qd"]=pf2*U2_plot[i,j]
            end
            network_data_plot=runpf(network_data_plot)
            #test_runpf(network_data_plot)
            violation_status, margin_plot=check_violation(network_data_plot,1)
            exact_plot[i,j]=margin_plot
            network_data_plot,sanity_check=cvxrs(network_data,"margin",network_data_plot)
            (Σ0,γ_plot)=network_data_plot["uncertainty"]
            cvxrs_plot[i,j]=γ_plot
        end
    end

    return U1_plot,U2_plot,exact_plot,cvxrs_plot,u_plot0
end


function scrs(network_data, max_iter_SCRS,phase_shift=false)
    num_gen=length(network_data["gen"])
    result_cvxr=Dict("pg"=>zeros(num_gen,max_iter_SCRS+1),"vg"=>zeros(num_gen,max_iter_SCRS+1),"alpha"=>zeros(num_gen,max_iter_SCRS+1),
        "obj"=> zeros(max_iter_SCRS+1),"worst_cost"=> zeros(max_iter_SCRS+1),"solver_time"=>zeros(max_iter_SCRS+1),"issue"=>0,"max_iter"=>max_iter_SCRS)

    result_cvxr["obj"][1]=network_data["cost"]
    [result_cvxr["pg"][gen["index"],1]=gen["gen_status"].*gen["pg"] for (i,gen) in network_data["gen"]]
    [result_cvxr["vg"][gen["index"],1]=gen["vg"] for (i,gen) in network_data["gen"]]
    [result_cvxr["alpha"][gen["index"],1]=gen["ptc_factor"] for (i,gen) in network_data["gen"]]

    for iter=2:result_cvxr["max_iter"]+1
        network_data,sanity_check=cvxrs(network_data,"obj",nothing,phase_shift)
        network_data=runpf(network_data)
        test_runpf(network_data)
        test_cvxrs(network_data,sanity_check)
        violation_status, margin=check_violation(network_data)

        result_cvxr["obj"][iter]=network_data["cost"]
        result_cvxr["worst_cost"][iter]=sanity_check["obj"]
        [result_cvxr["pg"][gen["index"],iter]=gen["gen_status"].*gen["pg"] for (i,gen) in network_data["gen"]]
        [result_cvxr["vg"][gen["index"],iter]=gen["vg"] for (i,gen) in network_data["gen"]]
        [result_cvxr["alpha"][gen["index"],iter]=gen["ptc_factor"] for (i,gen) in network_data["gen"]]
        result_cvxr["solver_time"][iter]=sanity_check["solve_time"]

        dpg_step=norm(result_cvxr["pg"][:,iter]-result_cvxr["pg"][:,iter-1])
        dvg_step=norm(result_cvxr["vg"][:,iter]-result_cvxr["vg"][:,iter-1])
        println("Iteration ",iter-1 ,": ",round(result_cvxr["obj"][iter],digits=2),"  ",round(sanity_check["obj"],digits=3),"   ", sanity_check["status"], "    ", round(dpg_step,digits=4),"   ", round(dvg_step,digits=4))

        #if dpg_step<1e-4 && dvg_step<1e-4; result_cvxr["max_iter"]=iter; break; end
    end

    return network_data, result_cvxr
end


function sample_pf(network_data,num_data,bounded=0)
    (Σ0,γ0)=network_data["uncertainty"]
    num_load=size(Σ0,1)

    pl0=zeros(Float64, num_load); [pl0[load["index"]]=load["status"].*load["pd"] for (i,load) in network_data["load"]]
    ql0=zeros(Float64, num_load); [ql0[load["index"]]=load["status"].*load["qd"] for (i,load) in network_data["load"]]
    power_factor=ql0./(pl0.+1e-3)

    u=randn(num_load, num_data)

    if bounded==1
        #http://extremelearning.com.au/how-to-generate-uniformly-random-points-on-n-spheres-and-n-balls/
        norm_u=sum(u.^2,dims=1).^(0.5)
        pl_sampled=pl0*ones(1,num_data)+γ0*cholesky(Σ0).L*((rand(1,num_data).^(1/num_load)).*u./norm_u)
        #pl_sampled=pl0*ones(1,num_data)+γ0*cholesky(Σ0).L*(u./norm_u)
    elseif bounded==0
        pl_sampled=pl0*ones(1,num_data)+γ0*cholesky(Σ0).L*u
    end
    ql_sampled=Diagonal(power_factor)*pl_sampled
    data_gencost=zeros(num_data)
    margin_all=zeros(num_data)
    sample_violation_status=Dict("vmag"=>zeros(num_data),"pg"=>zeros(num_data),
        "qg"=>zeros(num_data),"line"=>zeros(num_data),"angle"=>zeros(num_data))

    for idx_data=1:num_data
        network_data_run=deepcopy(network_data)

        [network_data_run["load"][i]["pd"]=pl_sampled[load["index"],idx_data] for (i,load) in network_data_run["load"]]
        [network_data_run["load"][i]["qd"]=ql_sampled[load["index"],idx_data] for (i,load) in network_data_run["load"]]

        network_data_run=runpf(network_data_run)
        #network_data_run=run_pf(network_data_run, ACPPowerModel, with_optimizer(Ipopt.Optimizer,print_level=0))

        #test_runpf(network_data_run)
        violation_status, margin_all[idx_data]=check_violation(network_data_run,1)
        data_gencost[idx_data]=network_data_run["cost"]
        sample_violation_status["vmag"][idx_data]=length(violation_status["vmag"])
        sample_violation_status["pg"][idx_data]=length(violation_status["pg"])
        sample_violation_status["qg"][idx_data]=length(violation_status["qg"])
        sample_violation_status["line"][idx_data]=length(violation_status["line"])
        sample_violation_status["angle"][idx_data]=length(violation_status["angle"])
    end
    violation_count=sum(margin_all.<=0)
    return pl_sampled,data_gencost,violation_count,sample_violation_status
end
