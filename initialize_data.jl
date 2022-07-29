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

   H_base = 3 .+ 6 .* rand(Float64,(length(instance.units)))
   R_base = 0.03 .+ 0.08 .* rand(Float64,(length(instance.units)))
   K_base =  0.8 .+ 0.4 .* rand(Float64,(length(instance.units)))
   F_base = 0.1 .+ 0.25 .* rand(Float64,(length(instance.units)))
   D_base = 0.6 .* ones(length(instance.units))
   T_r = 8.0
   pu = [instance.units[g].max_power[1]/sum(instance.units[g].max_power[1] for g in 1:length(instance.units))
            for g in 1:length(instance.units)]

   DData = DynamicsData(R_base,K_base,F_base,D_base,H_base,pu,T_r)
   return DData

end

function ESS_initialize(eta_c::Float64,
                        eta_d::Float64,
                        capacity::Float64,
                        discharge_lim::Float64,
                        charge_lim::Float64,
                        ess_percent::Float64)


   S = round(Int64,ess_percent*length(instance.buses))
   rand_position = rand(1:length(instance.buses),S)
   ess_bus = [instance.buses[rand_position[i]].name for i in 1:S]

   ESS_details = ESSData(eta_c,eta_d,capacity,discharge_lim,charge_lim,ess_bus)
   return ESS_details

end
