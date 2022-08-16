function add_ess(instance::UnitCommitmentInstance,
                  model::JuMP.Model,
                  ESSData::ESSData)

    T = instance.time
    S = length(ESSData.ess_bus)

    @variable(model, P_ch[i=1:S,j=1:T])
    @variable(model, P_dis[i=1:S, j=1:T])
    @variable(model, is_discharge[i=1:S, j=1:T], Bin)


    for s in 1:S
        for t in 1:T

        # Add production costs to the objective function.
            set_objective_coefficient(model, P_dis[s,t], 15.0)

            # Attach the new component to ESS buses, by modifying the
            # constraint `eq_net_injection`.
            set_normalized_coefficient(
                model[:eq_net_injection][ESSData.ess_bus[s], t],
                P_dis[s,t],
                1.0,
            )
            set_normalized_coefficient(
                model[:eq_net_injection][ESSData.ess_bus[s], t],
                P_ch[s,t],
                1.0,
            )

            if t==1
                @constraint(model, P_ch[s,t]==0)
                @constraint(model, P_dis[s,t]==0)
            else
                @constraint(model,0.1 <= 0.5-sum(P_dis[s,i]/(ESSData.eta_d*ESSData.EScap)+
                P_ch[s,i]*ESSData.eta_c/ESSData.EScap for i in 1:t) <= 0.9)
            end

             @constraint(model,P_dis[s,t]>=0)
             @constraint(model,P_dis[s,t]<=is_discharge[s,t]*ESSData.power_p)

             @constraint(model,P_ch[s,t]<=0)
             @constraint(model,P_ch[s,t]>=(1-is_discharge[s,t])*ESSData.power_n)

             #SOS! implementation (slower)
            # @constraint(model,P_dis[s,t]<=ESSData.power_p)
            # @constraint(model,P_ch[s,t]>=ESSData.power_n)
            # @constraint(model,[P_dis[s,t],P_ch[s,t]] in SOS1())

        end

        @constraint(model,sum(P_dis[s,i]/(ESSData.eta_d*ESSData.EScap)+
        P_ch[s,i]*ESSData.eta_c/ESSData.EScap for i in 1:T)==0)

    end

end
