
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

    gencost=zeros(Float64,num_gen,3); [gencost[gen["index"],4-size(gen["cost"],1):end]=gen["cost"] for (i,gen) in network_data["gen"]]

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
    power_factor=ql0./pl0

    vmax=zeros(Float64, num_bus); [vmax[bus["idx"]]=bus["vmax"] for (i,bus) in network_data["bus"]]
    vmin=zeros(Float64, num_bus); [vmin[bus["idx"]]=bus["vmin"] for (i,bus) in network_data["bus"]]
    φmax=zeros(Float64, num_line); [φmax[branch["index"]]=branch["angmax"] for (i,branch) in network_data["branch"]]
    φmin=zeros(Float64, num_line); [φmin[branch["index"]]=branch["angmin"] for (i,branch) in network_data["branch"]]
    pg_max=zeros(Float64, num_gen); [pg_max[gen["index"]]=gen["pmax"] for (i,gen) in network_data["gen"]]
    pg_min=zeros(Float64, num_gen); [pg_min[gen["index"]]=gen["pmin"] for (i,gen) in network_data["gen"]]
    qg_max=zeros(Float64, num_gen); [qg_max[gen["index"]]=gen["qmax"] for (i,gen) in network_data["gen"]]
    qg_min=zeros(Float64, num_gen); [qg_min[gen["index"]]=gen["qmin"] for (i,gen) in network_data["gen"]]
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

    @variable(m_cvx, vmin[i]<=vᵘ[i=1:num_bus]<=vmax[i])#,start=v0[i])
    @variable(m_cvx, vmin[i]<=vˡ[i=1:num_bus]<=vmax[i])#,start=v0[i]
    @variable(m_cvx, vmin[i]-v0[i]<=Δvᵘ[i=1:num_bus]<=vmax[i]-v0[i])#,start=0
    @variable(m_cvx, vmin[i]-v0[i]<=Δvˡ[i=1:num_bus]<=vmax[i]-v0[i])#,start=0
    @variable(m_cvx, φmin[i]<=φᵘ[i=1:num_line]<=φmax[i])#,start=φ0[i]
    @variable(m_cvx, φmin[i]<=φˡ[i=1:num_line]<=φmax[i])#,start=φ0[i]
    @variable(m_cvx, φmin[i]-φ0[i]<=Δφᵘ[i=1:num_line]<=φmax[i]-φ0[i])#,start=0
    @variable(m_cvx, φmin[i]-φ0[i]<=Δφˡ[i=1:num_line]<=φmax[i]-φ0[i])#,start=0
    @variable(m_cvx, deltaᵘ)#,start=delta0
    @variable(m_cvx, deltaˡ)#,start=delta0
    @variable(m_cvx, Δdeltaᵘ)#,start=0
    @variable(m_cvx, Δdeltaˡ)#,start=0
    @variable(m_cvx, gᵘ_pinj[1:num_bus])#,start=0
    @variable(m_cvx, gˡ_pinj[1:num_bus])#,start=0
    @variable(m_cvx, gᵘ_qinj[1:num_bus])#,start=0
    @variable(m_cvx, gˡ_qinj[1:num_bus])#,start=0
    @variable(m_cvx, gᵘ_vvsin[1:num_line])#,start=0
    @variable(m_cvx, gˡ_vvsin[1:num_line])#,start=0
    @variable(m_cvx, gᵘ_vvcos[1:num_line])#,start=0
    @variable(m_cvx, gˡ_vvcos[1:num_line])#,start=0
    @variable(m_cvx, gᵘ_vv[1:num_bus])#,start=0
    @variable(m_cvx, gˡ_vv[1:num_bus])#,start=0

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

    @variable(m_cvx, pg_opt[i=1:num_gen])#,start=pg0[i]
    @variable(m_cvx, pl_opt[i=1:num_load])#,start=pl0[i]
    @variable(m_cvx, ql_opt[i=1:num_load])#,start=ql0[i]
    @variable(m_cvx, pg_optᵘ[i=1:num_gen])#,start=pg0[i]
    @variable(m_cvx, alpha_opt[i=1:num_gen])#,start=alpha0[i]
    @variable(m_cvx, Δalpha_opt[i=1:num_gen])#,start=0
    @variable(m_cvx, cost_opt)#,start=network_data["cost"]
    @variable(m_cvx, γ_opt)#,start=0

    @variable(m_cvx, slack_pg)#,start=0
    @variable(m_cvx, slack_qg)#,start=0
    @variable(m_cvx, slack_line)#,start=0
    @variable(m_cvx, margin_opt)#,start=0

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
