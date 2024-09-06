struct DynamicsData
   R::Vector{Float64}
   K::Vector{Float64}
   F::Vector{Float64}
   D::Vector{Float64}
   H::Vector{Float64}
   pu::Vector{Float64}
   T_r::Float64
end
struct ESSData
   eta_c::Float64
   eta_d::Float64
   EScap::Float64
   power_p::Float64
   power_n::Float64
   ess_bus::Vector{String}
end

function freq_initialize(instance::UnitCommitmentInstance)
   thermal_units = instance.scenarios[1].thermal_units
   H_base = 3 .+ 6 .* rand(Float64, (length(thermal_units)))
   R_base = 0.03 .+ 0.08 .* rand(Float64, (length(thermal_units)))
   K_base = 0.8 .+ 0.4 .* rand(Float64, (length(thermal_units)))
   F_base = 0.1 .+ 0.25 .* rand(Float64, (length(thermal_units)))
   D_base = 0.6 .* ones(length(thermal_units))
   T_r = 8.0
   P_tot = sum(thermal_units[g].max_power[1] for g in eachindex(thermal_units))
   pu = [thermal_units[g].max_power[1] / P_tot for g in eachindex(thermal_units)]

   return DynamicsData(R_base, K_base, F_base, D_base, H_base, pu, T_r)

end

function ESS_initialize(
   instance::UnitCommitmentInstance,
   ess_percent::Float64
)
   eta_c = eta_d = 0.85
   EScap = 200.0
   power_p = 50.0
   power_n = -50.0
   S = round(Int64, ess_percent * length(instance.buses))
   rand_position = rand(1:length(instance.buses), S)
   ess_bus = [instance.buses[rand_position[i]].name for i in 1:S]

   return ESSData(eta_c, eta_d, EScap, power_p, power_n, ess_bus)
end
