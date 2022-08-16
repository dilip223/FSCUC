function linear_constraints_lazy(instance::UnitCommitmentInstance,
                  model::JuMP.Model,
                  DData::DynamicsData,
                  n_cont::Int,
                  constraint_mode::String)


    H = DData.H .* DData.pu
    RD = [(DData.D[g]+DData.K[g]/DData.R[g])*DData.pu[g] for g in 1:length(instance.units)]
    P_tot = sum(instance.units[g].max_power[1] for g in 1:length(instance.units))

    if constraint_mode == "Constant"
        #Constant model with k generator failures
        max_deltaP = sort([DData.pu[g] for g in 1:length(instance.units)] , rev=true)[1:n_cont]
        ind_loss = sortperm([DData.pu[g] for g in 1:length(instance.units)] , rev=true)[1:n_cont]
        RD_loss = sum(RD[ind_loss[i]] for i in 1:n_cont)
        H_loss = sum(H[ind_loss[i]] for i in 1:n_cont)

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

                viol_ss = -(sum(callback_value(cb_data,model[:is_on][instance.units[g].name,t])*RD[g] for g in 1:length(instance.units))
                - RD_loss - sum(max_deltaP)*60/0.3)

                viol_rocof = -(sum(H[g]*callback_value(cb_data,model[:is_on][instance.units[g].name,t]) for g in 1:length(instance.units))-H_loss - sum(max_deltaP)*60)

                if viol_ss > 1e-4
                    con = @build_constraint(sum(model[:is_on][instance.units[g].name,t]*RD[g] for g in 1:length(instance.units))
                                                - RD_loss - sum(max_deltaP)*60/0.3 >= 0)
                    #println("SS violation = $viol_ss in period $t")
                    MOI.submit(model, MOI.LazyConstraint(cb_data), con)

                end

                if viol_rocof > 1e-4
                    con = @build_constraint(sum(H[g]*model[:is_on][instance.units[g].name,t] for g in 1:length(instance.units))
                                                - H_loss - sum(max_deltaP)*60 >= 0)
                    #println("RoCoF violation = $viol_rocof in period $t")
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


            for t in 1:instance.time

                pow_mat = [DData.pu[g] *callback_value(cb_data,model[:is_on][instance.units[g].name,t]) for g in 1:length(instance.units)]
                pow_max = sort(pow_mat,rev=true)[1:n_cont]
                ind_max = sortperm(pow_mat,rev=true)[1:n_cont]

                viol_ss = -(sum(callback_value(cb_data,model[:is_on][instance.units[g].name,t])*RD[g] for g in 1:length(instance.units))
                          - sum(RD[ind_max[i]] for i in 1:n_cont)
                          - sum(pow_max)*60/0.3)

                viol_rocof = -(sum(H[i]*callback_value(cb_data,model[:is_on][instance.units[i].name,t]) for i in 1:length(instance.units))
                             -sum(H[ind_max[i]] for i in 1:n_cont)
                             -60*sum(pow_max))


                if viol_ss > 1e-4
                    println("SS Violation $viol_ss at $t. Max pow = $pow_max, index= $ind_max")
                    con = @build_constraint(sum(RD[g]*model[:is_on][instance.units[g].name,t] for g in 1:length(instance.units))
                                                - sum(RD[ind_max[i]]*model[:is_on][instance.units[ind_max[i]].name,t] for i in 1:n_cont)
                                                - 60*sum(DData.pu[ind_max[i]]*model[:is_on][instance.units[ind_max[i]].name,t] for i in 1:n_cont)/0.3 >= 0)

                    MOI.submit(model, MOI.LazyConstraint(cb_data), con)

                end

                if viol_rocof > 1e-4
                    println("RoCoF Violation $viol_rocof at $t. Max pow = $pow_max, index= $ind_max")
                    con = @build_constraint(sum(H[g]*model[:is_on][instance.units[g].name,t] for g in 1:length(instance.units))
                                                - sum(H[ind_max[i]]*model[:is_on][instance.units[ind_max[i]].name,t] for i in 1:n_cont)
                                                - 60*sum(DData.pu[ind_max[i]]*model[:is_on][instance.units[ind_max[i]].name,t] for i in 1:n_cont) >= 0)

                    MOI.submit(model, MOI.LazyConstraint(cb_data), con)

                end

            end
        end

        MOI.set(model, MOI.LazyConstraintCallback(), static_callback_function)
    end

    if constraint_mode == "Dynamic"
        #Dynamic Model (rocof constraint)
        @variable(model, temp_prod[1:instance.time,1:length(instance.units)])
        for t in 1:instance.time
            for g in 1:length(instance.units)
                        @constraint(model, temp_prod[t,g] <= model[:is_on][instance.units[g].name,t]*instance.units[g].max_power[1])
                        @constraint(model, temp_prod[t,g] >= 0)
                        @constraint(model, temp_prod[t,g] <= instance.units[g].min_power[1]
                                                + model[:prod_above][instance.units[g].name,t])
                        @constraint(model, temp_prod[t,g] >= instance.units[g].min_power[1]
                                                + model[:prod_above][instance.units[g].name,t]
                                                - (1-model[:is_on][instance.units[g].name,t])*instance.units[g].max_power[1])
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

                pow_mat = [(instance.units[g].min_power[1]+
                            callback_value(cb_data,model[:prod_above][instance.units[g].name,t]))*
                            callback_value(cb_data,model[:is_on][instance.units[g].name,t]) for g in 1:length(instance.units)] ./
                            P_tot
                pow_max = sort(pow_mat,rev=true)[1:n_cont]
                ind_max = sortperm(pow_mat,rev=true)[1:n_cont]

                viol_ss = -(sum(callback_value(cb_data,model[:is_on][instance.units[g].name,t])*RD[g] for g in 1:length(instance.units))
                - sum(RD[ind_max[i]] for i in 1:n_cont)
                 - sum(pow_max)*60/0.3)



                viol_rocof = -(sum(H[i]*callback_value(cb_data,model[:is_on][instance.units[i].name,t]) for i in 1:length(instance.units))
                            -sum(H[ind_max[i]] for i in 1:n_cont)
                            -60*sum(pow_max))

                if viol_rocof > 1e-4
                    println("RoCoF violation $viol_rocof at $t. Max pow = $pow_max, index= $ind_max")
                    con = @build_constraint(sum(H[g]*model[:is_on][instance.units[g].name,t] for g in 1:length(instance.units))
                                                - sum(H[ind_max[i]]*model[:is_on][instance.units[ind_max[i]].name,t] for i in 1:n_cont)
                                                - 60*sum(temp_prod[t,ind_max[i]] for i in 1:n_cont)./P_tot >= 0)

                    MOI.submit(model, MOI.LazyConstraint(cb_data), con)
                end
                if viol_ss > 1e-4
                    println("SS violation $viol_ss at $t. Max pow = $pow_max, index= $ind_max")
                    con = @build_constraint(sum(RD[g]*model[:is_on][instance.units[g].name,t] for g in 1:length(instance.units))
                                                - sum(RD[ind_max[i]]*model[:is_on][instance.units[ind_max[i]].name,t] for i in 1:n_cont)
                                                - 60/0.3*sum(temp_prod[t,ind_max[i]] for i in 1:n_cont)./P_tot >= 0)

                    MOI.submit(model, MOI.LazyConstraint(cb_data), con)

                end
            end
        end

        MOI.set(model, MOI.LazyConstraintCallback(), dynamic_callback_function)
    end
end



function linear_constraints_lazy(instance::UnitCommitmentInstance,
                  model::JuMP.Model,
                  DData::DynamicsData,
                  n_cont::Int,
                  constraint_mode::String,
                  ESSData::ESSData)


    H = DData.H .* DData.pu
    RD = [(DData.D[g]+DData.K[g]/DData.R[g])*DData.pu[g] for g in 1:length(instance.units)]
    P_tot = sum(instance.units[g].max_power[1] for g in 1:length(instance.units))

    if constraint_mode == "Constant"
        #Constant model with k generator failures
        max_deltaP = sort([instance.units[g].max_power[1] for g in 1:length(instance.units)]./
                    P_tot ,rev=true)[1:n_cont]
        ind_loss = sortperm([instance.units[g].max_power[1] for g in 1:length(instance.units)]./
                    P_tot ,rev=true)[1:n_cont]
        RD_loss = sum(RD[ind_loss[i]] for i in 1:n_cont)
        H_loss = sum(H[ind_loss[i]] for i in 1:n_cont)

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

                viol_ss = -(sum(callback_value(cb_data,model[:is_on][instance.units[g].name,t])*RD[g] for g in 1:length(instance.units))
                          - RD_loss
                          - sum(max_deltaP)*60/0.3
                           + 60/0.3 * sum(callback_value(cb_data,model[:P_dis][s,t]) for s in 1:length(ESSData.ess_bus)) ./ P_tot)

                viol_rocof = -(sum(H[g]*callback_value(cb_data,model[:is_on][instance.units[g].name,t]) for g in 1:length(instance.units))
                             - H_loss
                             - sum(max_deltaP)*60
                             + 60 * sum(callback_value(cb_data,model[:P_dis][s,t]) for s in 1:length(ESSData.ess_bus)) ./ P_tot)

                if viol_ss > 1e-4
                    con = @build_constraint(sum(model[:is_on][instance.units[g].name,t]*RD[g] for g in 1:length(instance.units))
                                                - RD_loss
                                                - sum(max_deltaP)*60/0.3
                                                + 60/0.3 * sum(model[:P_dis][s,t] for s in 1:length(ESSData.ess_bus)) ./ P_tot >= 0)
                    println("SS violation = $viol_ss in period $t")
                    MOI.submit(model, MOI.LazyConstraint(cb_data), con)

                end

                if viol_rocof > 1e-4
                    con = @build_constraint(sum(H[g]*model[:is_on][instance.units[g].name,t] for g in 1:length(instance.units))
                                                - H_loss
                                                - sum(max_deltaP)*60
                                                + 60 * sum(model[:P_dis][s,t] for s in 1:length(ESSData.ess_bus)) ./ P_tot >= 0)
                    println("RoCoF violation = $viol_rocof in period $t")
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


            for t in 1:instance.time

                pow_mat = [DData.pu[g] *callback_value(cb_data,model[:is_on][instance.units[g].name,t]) for g in 1:length(instance.units)]
                pow_max = sort(pow_mat,rev=true)[1:n_cont]
                ind_max = sortperm(pow_mat,rev=true)[1:n_cont]

                viol_ss = -(sum(callback_value(cb_data,model[:is_on][instance.units[g].name,t])*RD[g] for g in 1:length(instance.units))
                          - sum(RD[ind_max[i]] for i in 1:n_cont)
                          - sum(pow_max)*60/0.3
                          + 60/0.3 * sum(callback_value(cb_data,model[:P_dis][s,t]) for s in 1:length(ESSData.ess_bus)) ./ P_tot)

                viol_rocof = -(sum(H[i]*callback_value(cb_data,model[:is_on][instance.units[i].name,t]) for i in 1:length(instance.units))
                             - sum(H[ind_max[i]] for i in 1:n_cont)
                             - 60 * sum(pow_max)
                             + 60 * sum(callback_value(cb_data,model[:P_dis][s,t]) for s in 1:length(ESSData.ess_bus)) ./ P_tot)


                if viol_ss > 1e-4
                    println("SS Violation $viol_ss at $t. Max pow = $pow_max, index= $ind_max")
                    con = @build_constraint(sum(RD[g]*model[:is_on][instance.units[g].name,t] for g in 1:length(instance.units))
                                                - sum(RD[ind_max[i]]*model[:is_on][instance.units[ind_max[i]].name,t] for i in 1:n_cont)
                                                - 60/0.3 * sum(DData.pu[ind_max[i]]*model[:is_on][instance.units[ind_max[i]].name,t] for i in 1:n_cont)
                                                + 60/0.3 * sum(model[:P_dis][s,t] for s in 1:length(ESSData.ess_bus)) ./ P_tot >= 0)

                    MOI.submit(model, MOI.LazyConstraint(cb_data), con)

                end

                if viol_rocof > 1e-4
                    println("RoCoF Violation $viol_rocof at $t. Max pow = $pow_max, index= $ind_max")
                    con = @build_constraint(sum(H[g]*model[:is_on][instance.units[g].name,t] for g in 1:length(instance.units))
                                                - sum(H[ind_max[i]]*model[:is_on][instance.units[ind_max[i]].name,t] for i in 1:n_cont)
                                                - 60 * sum(DData.pu[ind_max[i]]*model[:is_on][instance.units[ind_max[i]].name,t] for i in 1:n_cont)
                                                + 60 * sum(model[:P_dis][s,t] for s in 1:length(ESSData.ess_bus)) ./ P_tot >= 0)

                    MOI.submit(model, MOI.LazyConstraint(cb_data), con)

                end

            end
        end

        MOI.set(model, MOI.LazyConstraintCallback(), static_callback_function)
    end

    if constraint_mode == "Dynamic"
        #Dynamic Model (rocof constraint)
        @variable(model, temp_prod[1:instance.time,1:length(instance.units)])
        for t in 1:instance.time
            for g in 1:length(instance.units)
                        @constraint(model, temp_prod[t,g] <= model[:is_on][instance.units[g].name,t]*instance.units[g].max_power[1])
                        @constraint(model, temp_prod[t,g] >= 0)
                        @constraint(model, temp_prod[t,g] <= instance.units[g].min_power[1]
                                                + model[:prod_above][instance.units[g].name,t])
                        @constraint(model, temp_prod[t,g] >= instance.units[g].min_power[1]
                                                + model[:prod_above][instance.units[g].name,t]
                                                - (1-model[:is_on][instance.units[g].name,t])*instance.units[g].max_power[1])
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

                pow_mat = [(instance.units[g].min_power[1]+
                            callback_value(cb_data,model[:prod_above][instance.units[g].name,t]))*
                            callback_value(cb_data,model[:is_on][instance.units[g].name,t]) for g in 1:length(instance.units)] ./
                            sum(instance.units[i].max_power[1] for i in 1:length(instance.units))
                pow_max = sort(pow_mat,rev=true)[1:n_cont]
                ind_max = sortperm(pow_mat,rev=true)[1:n_cont]

                viol_ss = -(sum(callback_value(cb_data,model[:is_on][instance.units[g].name,t])*RD[g] for g in 1:length(instance.units))
                            - sum(RD[ind_max[i]] for i in 1:n_cont)
                            - sum(pow_max)*60/0.3
                            + 60/0.3 * sum(callback_value(cb_data,model[:P_dis][s,t]) for s in 1:length(ESSData.ess_bus)) ./ P_tot)

                if viol_ss > 1e-4
                    println("SS violation $viol_ss at $t. Max pow = $pow_max, index= $ind_max")
                    con = @build_constraint(sum(RD[g]*model[:is_on][instance.units[g].name,t] for g in 1:length(instance.units))
                                                - sum(RD[ind_max[i]]*model[:is_on][instance.units[ind_max[i]].name,t] for i in 1:n_cont)
                                                - 60/0.3 * sum(temp_prod[t,ind_max[i]] for i in 1:n_cont)./sum(instance.units[i].max_power[1] for i in 1:length(instance.units))
                                                + 60/0.3 * sum(model[:P_dis][s,t] for s in 1:length(ESSData.ess_bus)) ./ P_tot >= 0)

                    MOI.submit(model, MOI.LazyConstraint(cb_data), con)

                end

                viol_rocof = -(sum(H[i]*callback_value(cb_data,model[:is_on][instance.units[i].name,t]) for i in 1:length(instance.units))
                                - sum(H[ind_max[i]] for i in 1:n_cont)
                                - 60 * sum(pow_max)
                                + 60 * sum(callback_value(cb_data,model[:P_dis][s,t]) for s in 1:length(ESSData.ess_bus)) ./ P_tot)

                if viol_rocof > 1e-4
                    println("RoCoF violation $viol_rocof at $t. Max pow = $pow_max, index= $ind_max")
                    con = @build_constraint(sum(H[g]*model[:is_on][instance.units[g].name,t] for g in 1:length(instance.units))
                                                - sum(H[ind_max[i]]*model[:is_on][instance.units[ind_max[i]].name,t] for i in 1:n_cont)
                                                - 60*sum(temp_prod[t,ind_max[i]] for i in 1:n_cont)./sum(instance.units[i].max_power[1] for i in 1:length(instance.units))
                                                + 60/0.3 * sum(model[:P_dis][s,t] for s in 1:length(ESSData.ess_bus)) ./ P_tot >= 0)

                    MOI.submit(model, MOI.LazyConstraint(cb_data), con)
                end
            end
        end

        MOI.set(model, MOI.LazyConstraintCallback(), dynamic_callback_function)
    end
end
