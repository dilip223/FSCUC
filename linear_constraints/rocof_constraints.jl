function rocof_constraint(
    instance::UnitCommitmentInstance,
    model::JuMP.Model,
    DData::DynamicsData,
    n_cont::Int,
    constraint_mode::String
)

    thermal_units = instance.scenarios[1].thermal_units
    H = DData.H .* DData.pu
    comb_units = collect(combinations(1:length(thermal_units), n_cont))

    if constraint_mode == "Constant"
        #Constant model with k generator failures
        max_deltaP = sort([thermal_units[g].max_power[1] for g in eachindex(thermal_units)] ./
                          sum(thermal_units[g].max_power[1] for g in eachindex(thermal_units)), rev=true)[1:n_cont]
        ind_Hloss = sortperm([thermal_units[g].max_power[1] for g in eachindex(thermal_units)] ./
                             sum(thermal_units[g].max_power[1] for g in eachindex(thermal_units)), rev=true)[1:n_cont]
        H_loss = sum(H[ind_Hloss])
        for t in 1:instance.time
            @constraint(model, sum(H[i] * model[:is_on][thermal_units[i].name, t]
                                   for i in eachindex(thermal_units)) - H_loss - sum(max_deltaP) * 60 >= 0)
        end
    end



    if constraint_mode == "Static"
        #Static Model (rocof constraint) with n_cont generator failures
        for t in 1:instance.time
            for g in eachindex(comb_units)
                @constraint(model, sum(H[i] * model[:is_on][thermal_units[i].name, t] for i in eachindex(thermal_units))
                                   -
                                   sum(H[comb_units[g][i]] * model[:is_on][thermal_units[comb_units[g][i]].name, t] for i in 1:n_cont)
                                   -
                                   60 * sum(model[:is_on][thermal_units[comb_units[g][i]].name, t] * DData.pu[comb_units[g][i]] for i in 1:n_cont) >= 0)
            end
        end

    end


    if constraint_mode == "Dynamic"
        #Dynamic Model (rocof constraint) with n_cont generator failures

        @variable(model, temp_prod[1:instance.time, eachindex(thermal_units)])

        for t in 1:instance.time
            for g in eachindex(comb_units)
                for k in 1:n_cont
                    @constraint(model, temp_prod[t, comb_units[g][k]] <= model[:is_on][thermal_units[comb_units[g][k]].name, t] * thermal_units[comb_units[g][k]].max_power[t])
                    @constraint(model, temp_prod[t, comb_units[g][k]] >= 0)
                    @constraint(model, temp_prod[t, comb_units[g][k]] <= thermal_units[comb_units[g][k]].min_power[t]
                                                                         +
                                                                         model[:prod_above][instance.scenarios[1].name, thermal_units[comb_units[g][k]].name, t])
                    @constraint(model, temp_prod[t, comb_units[g][k]] >= thermal_units[comb_units[g][k]].min_power[t]
                                                                         +
                                                                         model[:prod_above][instance.scenarios[1].name, thermal_units[comb_units[g][k]].name, t]
                                                                         -
                                                                         (1 - model[:is_on][thermal_units[comb_units[g][k]].name, t]) * thermal_units[comb_units[g][k]].max_power[t])
                end

                @constraint(model, sum(H[i] * model[:is_on][thermal_units[i].name, t] for i in eachindex(thermal_units))
                                   -
                                   sum(H[comb_units[g][i]] * model[:is_on][thermal_units[comb_units[g][i]].name, t] for i in 1:n_cont)
                                   -
                                   60 * sum(temp_prod[t, comb_units[g][i]] / sum(thermal_units[i].max_power[1] for i in eachindex(thermal_units)) for i in 1:n_cont) >= 0)
            end
        end

    end
end
#Dynamic Model (rocof constraint)
# @variable(model, temp_prod[1:instance.time,1:length(comb_units),1:n_cont])
#
# for t in 1:instance.time
#     for g in 1:length(comb_units)
#             for k in 1:n_cont
#                 @constraint(model, temp_prod[t,g,k] <= model[:is_on][thermal_units[comb_units[g][k]].name,t]*thermal_units[comb_units[g][k]].max_power[t])
#                 @constraint(model, temp_prod[t,g,k] >= 0)
#                 @constraint(model, temp_prod[t,g,k] <= thermal_units[comb_units[g][k]].min_power[t]
#                                         + model[:prod_above][instance.scenarios[1].name, thermal_units[comb_units[g][k]].name,t])
#                 @constraint(model, temp_prod[t,g,k] >= thermal_units[comb_units[g][k]].min_power[t]
#                                         + model[:prod_above][instance.scenarios[1].name, thermal_units[comb_units[g][k]].name,t]
#                                         - (1-model[:is_on][thermal_units[comb_units[g][k]].name,t])*thermal_units[comb_units[g][k]].max_power[t])
#             end
#             #
#             # @constraint(model,  sum(H[i]*model[:is_on][thermal_units[i].name,t] for i in eachindex(thermal_units))
#             #        - sum(H[comb_units[g][i]]*model[:is_on][thermal_units[comb_units[g][i]].name,t] for i in 1:n_cont)
#             #             - 60*sum(temp_prod[t,g,i]/sum(thermal_units[i].max_power[1] for i in eachindex(thermal_units)) for i in 1:n_cont) >= 0)
#      end
# end
