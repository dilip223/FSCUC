function warm_start_constraints(
    model::JuMP.Model,
    DData::DynamicsData,
    solution::AbstractDict
)

    instance, T = model[:instance], model[:instance].time
    thermal_units = instance.scenarios[1].thermal_units
    is_on = model[:is_on]
    switch_on = model[:switch_on]
    switch_off = model[:switch_off]
    for g in thermal_units
        for t in 1:T
            JuMP.set_start_value(is_on[g.name, t], solution["Is on"][g.name][t])
            JuMP.set_start_value(
                switch_on[g.name, t],
                solution["Switch on"][g.name][t],
            )
            JuMP.set_start_value(
                switch_off[g.name, t],
                solution["Switch off"][g.name][t],
            )
        end
    end
    for g in thermal_units
        for t in 1:T
            JuMP.set_start_value(
                model[:prod_above][instance.scenarios[1].name, g.name, t],
                solution["Thermal production (MW)"][g.name][t] - g.min_power[1] * solution["Is on"][g.name][t])
        end
    end
    # instance = model[:instance]
    # H_T = DData.H .* DData.pu
    # D_T = [(DData.D[g])*DData.pu[g] for g in 1:length(thermal_units)]
    # R_T = [(DData.K[g]/DData.R[g])*DData.pu[g] for g in 1:length(thermal_units)]
    # F_T = [(DData.F[g]*DData.K[g]/DData.R[g])*DData.pu[g] for g in 1:length(thermal_units)]
    # P_tot = sum(thermal_units[g].max_power[1] for g in 1:length(thermal_units))
    #
    #
    # for t in 1:instance.time
    #     H_current = sum(solution["Is on"][thermal_units[g].name][t]*H_T[g] for g in 1:length(thermal_units))
    #     R_current = sum(solution["Is on"][thermal_units[g].name][t]*R_T[g] for g in 1:length(thermal_units))
    #     F_current = sum(solution["Is on"][thermal_units[g].name][t]*F_T[g] for g in 1:length(thermal_units))
    #     D_current = sum(solution["Is on"][thermal_units[g].name][t]*D_T[g] for g in 1:length(thermal_units))
    #     max_pow,ind_max = findmax([solution["Production (MW)"][thermal_units[g].name][1] for g in 1:length(thermal_units)])
    #
    #     # @constraint(model, sum(H_T[g]*model[:is_on][thermal_units[g].name,t] for g in 1:length(thermal_units)) >= H_current)
    #     # @constraint(model, sum(R_T[g]*model[:is_on][thermal_units[g].name,t]  for g in 1:length(thermal_units)) >= R_current)
    #     # @constraint(model, sum(F_T[g]*model[:is_on][thermal_units[g].name,t]  for g in 1:length(thermal_units)) >= F_current)
    #     # @constraint(model, sum(D_T[g]*model[:is_on][thermal_units[g].name,t]  for g in 1:length(thermal_units)) >= D_current)
    #
    # end

end
