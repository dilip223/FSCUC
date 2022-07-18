using Cbc
using JuMP
using UnitCommitment
using Plots
using CPLEX
using Combinatorics
using Random

Random.seed!(1234)
rand(376)

struct DynamicsData
   R::Vector{Float64}
   K::Vector{Float64}
   F::Vector{Float64}
   D::Vector{Float64}
   H::Vector{Float64}
   pu::Vector{Float64}
   T_r::Float64
end

include("freq_response.jl")
include("rocof_constraints.jl")
# Load benchmark instance
instance = UnitCommitment.read_benchmark("matpower/case118/2017-02-01")
#instance = UnitCommitment.read("C:\\Dilip\\ANL\\UnitCommitment.jl-dev\\UnitCommitment.jl-dev\\case14-flex_contingencies.json")
# Build JuMP model
model = UnitCommitment.build_model(
    instance = instance,
    optimizer = CPLEX.Optimizer,
)

#Find max
# @variable(model, b[1:length(instance.units)], Bin)
# @variable(model, z)
# @expression(model, pow_gen[g = 1:length(instance.units)], instance.units[g].min_power[1]
#                             +model[:prod_above][instance.units[g].name,1])
#
# for g in 1:length(instance.units)
#     @constraint(model, z >= pow_gen[g])
#     @constraint(model, z-pow_gen[g] <= (800-0)*(1-b[g]))
# end
#@constraint(model, z <=200)
DData = freq_initialize(instance)


n_cont = 1 #Number of generator contingencies to consider
if n_cont >= 1
   constraint_mode = "Dynamic" #Specify constraint mode as "Constant", "Static", "Dynamic"
   rocof_constraint(instance, model, DData, n_cont, constraint_mode)
end


# Solve the model
UnitCommitment.optimize!(model)
solution = UnitCommitment.solution(model)
UnitCommitment.write("C:/Dilip/ANL/FSCUC/output.json", solution)
UnitCommitment.validate(instance, solution)


function solve_time(model::Model)
    return MOI.get(model, MOI.SolveTime())
end
solvetime = solve_time(model)

f_max = zeros((instance.time))
rocof = zeros((instance.time))

if n_cont>=1
   freq_test(instance,model,f_max,rocof,DData,n_cont)
end
plot(rocof.*60, legend=false, xlabel= "Time (hours)", ylabel= "RoCoF (Hz/s)")
tot_units =zeros((instance.time,1))
for i in 1:instance.time
       tot_units[i] = sum(value(model[:is_on][instance.units[g].name,i]) for g in 1:length(instance.units))
end
plot(tot_units,legend=false, xlabel= "Time (hours)", ylabel= "No. of units online")
# function freq_est(instance,model,f_max,rocof,R,K,F,D,H,n_cont)
#     # R_T=zeros((T, 1))
#     # F_T=zeros((T, 1))
#     # H_T=zeros((T, 1))
#     # D_T=zeros((T, 1))
#     for t in 1:T
#         R_T=0
#         F_T=0
#         H_T=0
#         D_T=0
#
#         pow_mat = [(instance.units[g].min_power[t]+value(model[:prod_above][instance.units[g].name,t])).*value(model[:is_on][instance.units[g].name,t]) for g in 1:length(instance.units)]
#         pow_max = sort(pow_mat,rev=true)[1:n_cont]
#         ind_max = sortperm(pow_mat,rev=true)[1:n_cont]
#
#         pow_max = sum(pow_max)./sum(instance.units[g].max_power[t]
#                                 for g in 1:length(instance.units))
#         temp = [value(model[:is_on][instance.units[g].name,t])
#                             for g in 1:length(instance.units)]
#         # print(temp)
#         #temp = ones((length(instance.units),1))
#         temp[ind_max].=0
#         for g in 1:length(instance.units)
#             #on_status = value(model[:is_on][instance.units[g].name,t])
#             on_status = temp[g]
#             const_par = instance.units[g].max_power[T]/sum(instance.units[g].max_power[1] for g in 1:length(instance.units))
#             #tot_pow[t] = instance.units[g].max_power[T].*value(model[:is_on][instance.units[g].name,t])
#              R_T = R_T+(K[g]/R[g])*on_status*const_par
#              F_T = F_T+(F[g]*K[g]/R[g])*on_status*const_par
#              H_T = H_T+ H[g]*on_status*const_par
#              D_T = D_T+ D[g]*on_status*const_par
#         end
#         omega_n = sqrt((D_T+R_T)/(2*H_T*T_r))
#         zeta = (2*H_T+T_r*(D_T+F_T))/(2*sqrt(2*T_r*H_T*(D_T+R_T)))
#         omega_d = omega_n*sqrt(1-zeta^2)
#         t_max = (1/omega_d) *atan(omega_d/(zeta*omega_n-1/T_r))
#         f_max[t] = pow_max/(D_T+R_T)*(1+sqrt(T_r*(R_T-F_T)/(2*H_T))*
#         exp(-zeta*omega_n*t_max))
#         rocof[t] = pow_max/(2*H_T)
#      end
#
# end
# f_max = zeros((instance.time,1))
# rocof = zeros((instance.time,1))
# freq_test(instance,model,f_max,rocof,R_base,K_base,F_base,D_base,H_base,n_cont)
#
