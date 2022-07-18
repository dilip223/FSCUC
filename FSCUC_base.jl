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
