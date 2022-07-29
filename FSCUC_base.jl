using Cbc
using JuMP
using UnitCommitment
using Plots
using CPLEX
using Combinatorics
using Random
Random.seed!(1234)
rand(376)

include("initialize_data.jl")
include("freq_response.jl")
include("rocof_constraints.jl")
include("linear_constraints_lazy.jl")
include("ess_operation.jl")

# Load benchmark instance
instance = UnitCommitment.read_benchmark("matpower/case118/2017-02-01")
#instance = UnitCommitment.read("C:\\Dilip\\ANL\\UnitCommitment.jl-dev\\UnitCommitment.jl-dev\\case14-flex_contingencies.json")
# Build JuMP model
instance.time = 12
model = UnitCommitment.build_model(
    instance = instance,
    optimizer = CPLEX.Optimizer,
)


#Add a certain number of ESS to the system. It's in percentage of total number of buses
eta_c = eta_d = 0.85
EScap = 200.0
power_p = 50.0
power_n = -50.0
Total_percent = 0.05
ESSdata = ESS_initialize(eta_c,eta_d,EScap,power_p,power_n,Total_percent)
S = add_ess(instance,model,ESSdata)

#Initialize the generator dynamics
DData = freq_initialize(instance)
#Invoke the rocof constraints considering a selected number of generator contingencies
n_cont = 1#Number of generator contingencies to consider

if n_cont >= 1
   constraint_mode = "Dynamic" #Specify constraint mode as "Constant", "Static", "Dynamic"
   linear_constraints_lazy(instance, model, DData, n_cont, constraint_mode)
end

# Solve the model
MOI.set(model, MOI.NumberOfThreads(), 1)
start = time()
UnitCommitment.optimize!(model)
finish = time()
tottime = finish - start

solution = UnitCommitment.solution(model)
UnitCommitment.write("C:/Dilip/ANL/FSCUC/output.json", solution)
UnitCommitment.validate(instance, solution)

# #Get the solver time from JuMP
# function solve_time(model::Model)
#     return MOI.get(model, MOI.SolveTime())
# end
# solvetime = solve_time(model)

#Get the frequency nadir and maximum RoCoF value
if n_cont==0
   n_cont = n_cont+1
end
f_max = zeros((instance.time))
rocof = zeros((instance.time))
steady_freq = zeros((instance.time))
freq_test(instance, model, f_max, rocof, steady_freq, DData, n_cont)

tot_units =zeros((instance.time,1))
for i in 1:instance.time
       tot_units[i] = sum(value(model[:is_on][instance.units[g].name,i]) for g in 1:length(instance.units))
end


soc=zeros((length(ESSdata.ess_bus), instance.time))

for s in 1:length(ESSdata.ess_bus)
    for i in 1:instance.time
        soc[s,i]=0.5-sum(value(model[:P_dis][s,j])/(ESSdata.eta_d*ESSdata.EScap)+
        value(model[:P_ch][s,j])*ESSdata.eta_c/ESSdata.EScap for j in 1:i)
    end
end
