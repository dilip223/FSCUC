###Addition of the ESS
S = 2
ess_bus = ["b8","b9"]
eta_c = eta_d = 0.85
ESS_cap = 200
ESS_power_p = 50
ESS_power_n = -50

#create new continuous variables 0 ≤ x[t] ≤ 10.
@variable(model, P_ch[i=1:S,j=1:T])
@variable(model, P_dis[i=1:S, j=1:T])
@variable(model, is_discharge[i=1:S, j=1:T], Bin)

for s in 1:S
    for t in 1:T

    # Add production costs to the objective function.
        set_objective_coefficient(model, P_dis[s,t], 10.0)
        #set_objective_coefficient(model, P_ch[s,t], 15.0)
        #add_to_expression!(
            #     model[:obj],
            #     P_dis[s,t],
            #     5.0,
            # )
        # Attach the new component to bus b1, by modifying the
        # constraint `eq_net_injection`.
        set_normalized_coefficient(
        #add_to_expression!(
            model[:eq_net_injection][ess_bus[s], t],
            P_dis[s,t],
            1.0,
        )
        set_normalized_coefficient(
        #add_to_expression!(
            model[:eq_net_injection][ess_bus[s], t],
            P_ch[s,t],
            1.0,
        )

        if t==1
        #    @constraint(model,0.1 <= 0.5-P_ch[t]/(0.85*100)-P_dis[t]*0.85/100 <= 0.95)
        #    @constraint(model,0.5-x[t]/(0.85*100)-y[t]*0.85/100>=0.1)
            @constraint(model, P_ch[s,t]==0)
            @constraint(model, P_dis[s,t]==0)
        else
        #     @constraint(model,x[t-1]/(0.85*100)+y[t-1]*0.85/100
        #              -x[t]/(0.85*100)+y[t]*0.85/100<=0.9)
            @constraint(model,0.1 <= 0.5-sum(P_dis[s,i]/(eta_d*ESS_cap)+
            P_ch[s,i]*eta_c/ESS_cap for i in 1:t) <= 0.9)
            # @constraint(model,0.5-sum(P_dis[s,i]/(eta_d*ESS_cap)+
            # P_ch[s,i]*eta_c/ESS_cap for i in 1:t)>=0.1)
        end
        @constraint(model,P_dis[s,t]>=0)
        @constraint(model,P_dis[s,t]<=is_discharge[s,t]*ESS_power_p)
        @constraint(model,P_ch[s,t]<=0)
        @constraint(model,P_ch[s,t]>=(1-is_discharge[s,t])*ESS_power_n)
        #@constraint(model,is_charge[s,t]+is_discharge[s,t]<=1)
    end

    @constraint(model,sum(P_dis[s,i]/(eta_d*ESS_cap)+
    P_ch[s,i]*eta_c/ESS_cap for i in 1:T)==0)

end

JuMP.fix(
    model[:P_dis][1,6],
    50.0,
    force=true,
)
JuMP.fix(
    model[:P_dis][2,20],
    50.0,
    force=true,
)

# #Plot the SoC
# soc=zeros((S, T))
#
# for s in 1:S
#     for i in 1:T
#         soc[s,i]=0.5-sum(value(model[:P_dis][s,j])/(eta_d*ESS_cap)+
#         value(model[:P_ch][s,j])*eta_c/ESS_cap for j in 1:i)
#     end
# end
