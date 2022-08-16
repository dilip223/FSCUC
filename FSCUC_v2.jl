using Cbc
using JuMP
using UnitCommitment
using Plots
using CPLEX
using Combinatorics
using Random
using Ipopt
Random.seed!(1234)
rand(376)

include("initialize/initialize_data.jl")
include("test_response/freq_response.jl")
include("linear_constraints/rocof_constraints_lazy.jl")
include("linear_constraints/linear_constraints_lazy.jl")
include("linear_constraints/rocof_constraints.jl")
include("ESS_operation_constraints/ess_operation.jl")
include("nadir_constraints/lazy_freq_nadir_dynamic.jl")
include("initialize/warm_start_nadir.jl")
# Load benchmark instance
instance = UnitCommitment.read_benchmark("matpower/case118/2017-02-01")
#instance = UnitCommitment.read("C:\\Dilip\\ANL\\UnitCommitment.jl-dev\\UnitCommitment.jl-dev\\case14-flex_contingencies.json")
# Build JuMP model
instance.time = 24
model = UnitCommitment.build_model(
    instance = instance,
    optimizer = CPLEX.Optimizer,
)

#set_optimizer_attribute(model, "CPXPARAM_MIP_Tolerances_MIPGap", 0.001)
#Initialize the generator dynamics
DData = freq_initialize(instance)
#Invoke the rocof constraints considering a selected number of generator contingencies
n_cont = 1 #Number of generator contingencies to consider
ESS_ON = 0
# Total_percent = 0.05
# ESSdata = ESS_initialize(instance,Total_percent)
# add_ess(instance,model,ESSdata)

if n_cont >= 1
   constraint_mode = "Dynamic" #Specify constraint mode as "Constant", "Static", "Dynamic"

   if ESS_ON == 1
       #Add a certain number of ESS to the system. It's in percentage of total number of buses
       Total_percent = 0.04
       ESSdata = ESS_initialize(instance,Total_percent)
       add_ess(instance,model,ESSdata)
       linear_constraints_lazy(instance, model, DData, n_cont, constraint_mode, ESSdata)
   else
       linear_constraints_lazy(instance, model, DData, n_cont, constraint_mode)
   end
end
# Solve the model
MOI.set(model, MOI.NumberOfThreads(), 1)
start = time()
UnitCommitment.optimize!(model)
solution = UnitCommitment.solution(model)
if n_cont >=1
    model = UnitCommitment.build_model(
        instance = instance,
        optimizer = CPLEX.Optimizer,)
    warm_start_constraints(model,DData,solution)
    lazy_freq_nadir_v2(instance, model, DData, n_cont, constraint_mode)
    UnitCommitment.optimize!(model)
end
finish = time()
tottime = finish - start


#UnitCommitment.write("C:/Dilip/ANL/FSCUC/output.json", solution)
UnitCommitment.validate(instance, solution)

# #Get the solver time from JuMP
# function solve_time(model::Model)
#     return MOI.get(model, MOI.SolveTime())
# end
# solvetime = solve_time(model)

#Get the frequency nadir and maximum RoCoF value
f_max = zeros((instance.time))
rocof = zeros((instance.time))
steady_freq = zeros((instance.time))
freq_test(instance, model, f_max, rocof, steady_freq, DData, n_cont, ESS_ON)

tot_units =zeros((instance.time))

for i in 1:instance.time
       tot_units[i] = sum(value(model[:is_on][instance.units[g].name,i]) for g in 1:length(instance.units))
end

if ESS_ON & n_cont == 1
    soc=zeros((length(ESSdata.ess_bus), instance.time))

    for s in 1:length(ESSdata.ess_bus)
        for i in 1:instance.time
            soc[s,i]=0.5-sum(value(model[:P_dis][s,j])/(ESSdata.eta_d*ESSdata.EScap)+
            value(model[:P_ch][s,j])*ESSdata.eta_c/ESSdata.EScap for j in 1:i)
        end
    end
end
