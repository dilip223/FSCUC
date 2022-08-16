using JuMP
using Ipopt
using UnitCommitment
using Random
Random.seed!(1234)
rand(376)

M=10
N=3

H_final  = zeros((M,1))
R_final = zeros((M,1))
F_final = zeros((M,1))
f_max = zeros((M,1))

for i in 1:M
   instance = UnitCommitment.read_benchmark("matpower/case118/2017-02-01")
   H_base = 3 .+ 6 .* rand(Float64,(length(instance.units)))
   R_base = 0.03 .+ 0.08 .* rand(Float64,(length(instance.units)))
   K_base =  0.8 .+ 0.4 .* rand(Float64,(length(instance.units)))
   F_base = 0.1 .+ 0.25 .* rand(Float64,(length(instance.units)))
   D_base = 0.6 .* ones(length(instance.units))
   T_r = 8.0
   pu = [instance.units[g].max_power[1]/sum(instance.units[g].max_power[1] for g in 1:length(instance.units))
            for g in 1:length(instance.units)]

   R_T= sum((K_base[g]/R_base[g])*pu[g] for g in 1:length(instance.units))
   F_T= sum((F_base[g]*K_base[g]/R_base[g])*pu[g] for g in 1:length(instance.units))
   H_T= sum(H_base[g]*pu[g] for g in 1:length(instance.units))
   D_T= sum(D_base[g]*pu[g] for g in 1:length(instance.units))

   omega_n = sqrt((D_T+R_T)/(2*H_T*T_r))
   zeta = (2*H_T+T_r*(D_T+F_T))/(2*sqrt(2*T_r*H_T*(D_T+R_T)))
   omega_d = omega_n*sqrt(1-zeta^2)
   t_max = (1/omega_d) *atan(omega_d/(zeta*omega_n-1/T_r))
   f_max[i] = 0.1/(D_T+R_T)*(1+sqrt(T_r*(R_T-F_T)/(2*H_T))*
              exp(-zeta*omega_n*t_max))
   H_final[i]=H_T
   F_final[i]=F_T
   R_final[i]=R_T
end
model = JuMP.Model(Ipopt.Optimizer)

@variable(model, c[1:M,1:N])
@variable(model, b[1:M])

@NLobjective(model, Min, sum((max(c[i,1]*H_final[i]+c[i,2]*R_final[i]+c[i,3]*F_final[i]+b[i])-f_max[i])^2 for i in 1:M))
JuMP.optimize!(model)
