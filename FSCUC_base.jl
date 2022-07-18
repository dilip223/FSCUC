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

#Initialize the generator dynamics
DData = freq_initialize(instance)
#Invoke the rocof constraints considering a selected number of generator contingencies
n_cont = 0 #Number of generator contingencies to consider
if n_cont >= 1
   constraint_mode = "Constant" #Specify constraint mode as "Constant", "Static", "Dynamic"
   rocof_constraint(instance, model, DData, n_cont, constraint_mode)
end

# Solve the model
UnitCommitment.optimize!(model)
solution = UnitCommitment.solution(model)
UnitCommitment.write("C:/Dilip/ANL/FSCUC/output.json", solution)
UnitCommitment.validate(instance, solution)

#Get the solver time from JuMP
function solve_time(model::Model)
    return MOI.get(model, MOI.SolveTime())
end
solvetime = solve_time(model)

#Get the frequency nadir and maximum RoCoF value
f_max = zeros((instance.time))
rocof = zeros((instance.time))
freq_test(instance,model,f_max,rocof,DData,n_cont)
plot_figures(instance, model, f_max, rocof)
