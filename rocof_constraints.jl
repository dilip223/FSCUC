function rocof_constraint(instance::UnitCommitmentInstance,
                  model::JuMP.Model,
                  DData::DynamicsData,
                  n_cont::Int,
                  constraint_mode::String)

    H = DData.H .* DData.pu
    comb_units = collect(combinations(1:length(instance.units),n_cont))

    if constraint_mode == "Constant"
        #Constant model with k generator failures
        max_deltaP = sort([instance.units[g].max_power[1] for g in 1:length(instance.units)]./
                    sum(instance.units[g].max_power[1] for g in 1:length(instance.units)),rev=true)[1:n_cont]
        ind_Hloss = sortperm([instance.units[g].max_power[1] for g in 1:length(instance.units)]./
                    sum(instance.units[g].max_power[1] for g in 1:length(instance.units)),rev=true)[1:n_cont]
        H_loss = sum(H[ind_Hloss])
        for t in 1:instance.time
            @constraint(model,  sum(H[i]*model[:is_on][instance.units[i].name,t]
             for i in 1:length(instance.units)) - H_loss >= sum(max_deltaP)*60)
        end
    end



    if constraint_mode == "Static"
        #Static Model (rocof constraint) with k generator failures
        for t in 1:instance.time
            for g in 1:length(comb_units)
                @constraint(model,  sum(H[i]*model[:is_on][instance.units[i].name,t] for i in 1:length(instance.units))
                - sum(H[comb_units[g][i]]*model[:is_on][instance.units[comb_units[g][i]].name,t] for i in 1:n_cont)
                     - 60*sum(model[:is_on][instance.units[comb_units[g][i]].name,t]*DData.pu[comb_units[g][i]] for i in 1:n_cont) >= 0)
            end
        end
    end


    if constraint_mode == "Dynamic"
        #Dynamic Model (rocof constraint)
        @variable(model, temp_prod[1:instance.time,1:length(comb_units),1:n_cont])

        for t in 1:instance.time
            for g in 1:length(comb_units)
                    for k in 1:n_cont
                        @constraint(model, temp_prod[t,g,k] <= model[:is_on][instance.units[comb_units[g][k]].name,t]*instance.units[comb_units[g][k]].max_power[t])
                        @constraint(model, temp_prod[t,g,k] >= 0)
                        @constraint(model, temp_prod[t,g,k] <= instance.units[comb_units[g][k]].min_power[t]
                                                + model[:prod_above][instance.units[comb_units[g][k]].name,t])
                        @constraint(model, temp_prod[t,g,k] >= instance.units[comb_units[g][k]].min_power[t]
                                                + model[:prod_above][instance.units[comb_units[g][k]].name,t]
                                                - (1-model[:is_on][instance.units[comb_units[g][k]].name,t])*instance.units[comb_units[g][k]].max_power[t])
                    end

                    @constraint(model,  sum(H[i]*model[:is_on][instance.units[i].name,t] for i in 1:length(instance.units))
                            - sum(H[comb_units[g][i]]*model[:is_on][instance.units[comb_units[g][i]].name,t] for i in 1:n_cont)
                                 - 60*sum(temp_prod[t,g,i]/sum(instance.units[i].max_power[1] for i in 1:length(instance.units)) for i in 1:n_cont) >= 0)
             end
        end
    end
end
##Dynamic Model (rocof constraint)
# @variable(model, temp_prod[1:T,1:length(instance.units)])
# #
# for t in 1:T
#     for g in 1:length(instance.units)
#
#           @constraint(model, temp_prod[t,g] <= model[:is_on][instance.units[g].name,t]*50000)
#           @constraint(model, temp_prod[t,g] >= 0)
#           @constraint(model, temp_prod[t,g] <= instance.units[g].min_power[t] + model[:prod_above][instance.units[g].name,t])
#           @constraint(model, temp_prod[t,g] >= instance.units[g].min_power[t] + model[:prod_above][instance.units[g].name,t]
#                         - (1-model[:is_on][instance.units[g].name,t])*50000)
#           @constraint(model,  sum(H[i]*model[:is_on][instance.units[i].name,t]
#              for i in 1:length(instance.units)) - H[g]*model[:is_on][instance.units[g].name,t]
#              >= 60*temp_prod[t,g]/sum(instance.units[g].max_power[1] for g in 1:length(instance.units)))
#      end
# end
