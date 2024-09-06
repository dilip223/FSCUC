function rocof_constraint_lazy(
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

        function constant_callback_function(cb_data)
            status = callback_node_status(cb_data, model)

            if status == MOI.CALLBACK_NODE_STATUS_FRACTIONAL
                # `callback_value(cb_data, x)` is not integer (to some tolerance).
                # If, for example, your lazy constraint generator requires an
                # integer-feasible primal solution, you can add a `return` here.
                return
            elseif status == MOI.CALLBACK_NODE_STATUS_INTEGER
                # `callback_value(cb_data, x)` is integer (to some tolerance).
                #return
            else
                @assert status == MOI.CALLBACK_NODE_STATUS_UNKNOWN
                # `callback_value(cb_data, x)` might be fractional or integer.
            end

            for t in 1:instance.time
                # x_val = zeros((length(thermal_units)))
                # for g in length(x_val)
                #     x_val[g] = callback_value(cb_data, model[:is_on][thermal_units[g].name,t])
                # end
                viol = -(sum(H[g] * callback_value(cb_data, model[:is_on][thermal_units[g].name, t]) for g in eachindex(thermal_units))
                         -
                         H_loss - sum(max_deltaP) * 60)

                if viol > 1e-5
                    con = @build_constraint(sum(H[g] * model[:is_on][thermal_units[g].name, t] for g in eachindex(thermal_units))
                                            -
                                            H_loss - sum(max_deltaP) * 60 >= 0)
                    println("violation = $viol in period $t")
                    MOI.submit(model, MOI.LazyConstraint(cb_data), con)

                end
            end
        end

        MOI.set(model, MOI.LazyConstraintCallback(), constant_callback_function)
    end




    if constraint_mode == "Static"

        function static_callback_function(cb_data)
            status = callback_node_status(cb_data, model)
            if status == MOI.CALLBACK_NODE_STATUS_FRACTIONAL
                return

            elseif status == MOI.CALLBACK_NODE_STATUS_INTEGER

            else
                @assert status == MOI.CALLBACK_NODE_STATUS_UNKNOWN
            end

            # for t in 1:instance.time
            #     for g in 1:length(comb_units)
            #
            #         viol = sum(H[i]*callback_value(cb_data,model[:is_on][thermal_units[i].name,t]) for i in eachindex(thermal_units))-
            #                sum(H[comb_units[g][i]]*callback_value(cb_data,model[:is_on][thermal_units[comb_units[g][i]].name,t]) for i in 1:n_cont) - 60*sum(callback_value(cb_data,model[:is_on][thermal_units[comb_units[g][i]].name,t])*DData.pu[comb_units[g][i]] for i in 1:n_cont)
            #            status = callback_value(cb_data,model[:is_on][thermal_units[comb_units[g][1]].name,t])
            #            gen_value = comb_units[g][1]
            #            # term1 = sum(H[i]*callback_value(cb_data,model[:is_on][thermal_units[i].name,t]) for i in eachindex(thermal_units))
            #            # term2 = sum(H[comb_units[g][i]]*callback_value(cb_data,model[:is_on][thermal_units[comb_units[g][i]].name,t]) for i in 1:n_cont)
            #            # term3 = 60*sum(callback_value(cb_data,model[:is_on][thermal_units[comb_units[g][i]].name,t])*DData.pu[comb_units[g][i]] for i in 1:n_cont)
            #            if (viol <1e-4)
            #                #println("Violation $viol at time $t")
            #
            #               # println("Violation $viol at time $t")
            #               # println("Generator number $gen_value")
            #                con = @build_constraint(sum(H[i]*model[:is_on][thermal_units[i].name,t] for i in eachindex(thermal_units))-
            #                                        sum(H[comb_units[g][i]]*model[:is_on][thermal_units[comb_units[g][i]].name,t] for i in 1:n_cont)- 60*sum(model[:is_on][thermal_units[comb_units[g][i]].name,t]*DData.pu[comb_units[g][i]] for i in 1:n_cont) >= 0)
            #
            #                MOI.submit(model, MOI.LazyConstraint(cb_data), con)
            #           end
            #     end
            # end

            for t in 1:instance.time

                pow_mat = [DData.pu[g] * callback_value(cb_data, model[:is_on][thermal_units[g].name, t]) for g in eachindex(thermal_units)]
                pow_max = sort(pow_mat, rev=true)[1:n_cont]
                ind_max = sortperm(pow_mat, rev=true)[1:n_cont]
                viol = sum(H[i] * callback_value(cb_data, model[:is_on][thermal_units[i].name, t]) for i in eachindex(thermal_units)) - sum(H[ind_max[i]] for i in 1:n_cont) - 60 * sum(pow_max)

                if viol < 1e-6
                    println("Violation $viol at $t. Max pow = $pow_max, index= $ind_max")
                    con = @build_constraint(sum(H[g] * model[:is_on][thermal_units[g].name, t] for g in eachindex(thermal_units)) - sum(H[ind_max[i]] * model[:is_on][thermal_units[ind_max[i]].name, t] for i in 1:n_cont) - 60 * sum(DData.pu[ind_max[i]] * model[:is_on][thermal_units[ind_max[i]].name, t] for i in 1:n_cont) >= 0)
                    MOI.submit(model, MOI.LazyConstraint(cb_data), con)
                end


            end
        end

        MOI.set(model, MOI.LazyConstraintCallback(), static_callback_function)
    end

    # for g in 1:length(comb_units)
    #
    #     if (sum(H[i]*callback_value(cb_data,model[:is_on][thermal_units[i].name,t]) for i in eachindex(thermal_units))
    #     - sum(H[comb_units[g][i]]*callback_value(cb_data,model[:is_on][thermal_units[comb_units[g][i]].name,t]) for i in 1:n_cont)
    #          - 60*sum(callback_value(cb_data,model[:is_on][thermal_units[comb_units[g][i]].name,t])*DData.pu[comb_units[g][i]] for i in 1:n_cont) < 0)
    #
    #             con = @build_constraint(sum(H[i]*model[:is_on][thermal_units[i].name,t] for i in eachindex(thermal_units))
    #             - sum(H[comb_units[g][i]]*model[:is_on][thermal_units[comb_units[g][i]].name,t] for i in 1:n_cont)
    #                  - 60*sum(model[:is_on][thermal_units[comb_units[g][i]].name,t]*DData.pu[comb_units[g][i]] for i in 1:n_cont) >= 0)
    #
    #             MOI.submit(model, MOI.LazyConstraint(cb_data), con)
    #     end
    # end

    if constraint_mode == "Dynamic"
        #Dynamic Model (rocof constraint)
        @variable(model, temp_prod[1:instance.time, length(thermal_units)])
        for t in 1:instance.time
            for g in eachindex(thermal_units)
                @constraint(model, temp_prod[t, g] <= model[:is_on][thermal_units[g].name, t] * thermal_units[g].max_power[1])
                @constraint(model, temp_prod[t, g] >= 0)
                @constraint(model, temp_prod[t, g] <= thermal_units[g].min_power[1]
                                                      +
                                                      model[:prod_above][instance.scenarios[1].name, thermal_units[g].name, t])
                @constraint(model, temp_prod[t, g] >= thermal_units[g].min_power[1]
                                                      +
                                                      model[:prod_above][instance.scenarios[1].name, thermal_units[g].name, t]
                                                      -
                                                      (1 - model[:is_on][thermal_units[g].name, t]) * thermal_units[g].max_power[1])
            end
        end

        function dynamic_callback_function(cb_data)
            status = callback_node_status(cb_data, model)
            if status == MOI.CALLBACK_NODE_STATUS_FRACTIONAL

                return

            elseif status == MOI.CALLBACK_NODE_STATUS_INTEGER

            else
                @assert status == MOI.CALLBACK_NODE_STATUS_UNKNOWN

            end

            for t in 1:instance.time

                pow_mat = [(thermal_units[g].min_power[1] +
                            callback_value(cb_data, model[:prod_above][instance.scenarios[1].name, thermal_units[g].name, t])) *
                           callback_value(cb_data, model[:is_on][thermal_units[g].name, t]) for g in eachindex(thermal_units)] ./
                          sum(thermal_units[i].max_power[1] for i in eachindex(thermal_units))
                pow_max = sort(pow_mat, rev=true)[1:n_cont]
                ind_max = sortperm(pow_mat, rev=true)[1:n_cont]
                viol = -(sum(H[i] * callback_value(cb_data, model[:is_on][thermal_units[i].name, t]) for i in eachindex(thermal_units))
                         -
                         sum(H[ind_max[i]] for i in 1:n_cont)
                         -
                         60 * sum(pow_max))


                if viol > 1e-4
                    println("Violation $viol at $t. Max pow = $pow_max, index= $ind_max")
                    con = @build_constraint(sum(H[g] * model[:is_on][thermal_units[g].name, t] for g in eachindex(thermal_units))
                                            -
                                            sum(H[ind_max[i]] * model[:is_on][thermal_units[ind_max[i]].name, t] for i in 1:n_cont)
                                            -
                                            60 * sum(temp_prod[t, ind_max[i]] for i in 1:n_cont) ./ sum(thermal_units[i].max_power[1] for i in eachindex(thermal_units)) >= 0)

                    MOI.submit(model, MOI.LazyConstraint(cb_data), con)

                end
            end
        end
        MOI.set(model, MOI.LazyConstraintCallback(), dynamic_callback_function)
    end
end
