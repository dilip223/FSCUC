function lazy_freq_nadir_v2(instance::UnitCommitmentInstance,
                  model::JuMP.Model,
                  DData::DynamicsData,
                  n_cont::Int,
                  constraint_mode::String)


  H_T = DData.H .* DData.pu
  D_T = [(DData.D[g])*DData.pu[g] for g in 1:length(instance.units)]
  R_T = [(DData.K[g]/DData.R[g])*DData.pu[g] for g in 1:length(instance.units)]
  F_T = [(DData.F[g]*DData.K[g]/DData.R[g])*DData.pu[g] for g in 1:length(instance.units)]
  P_tot = sum(instance.units[g].max_power[1] for g in 1:length(instance.units))

  if constraint_mode == "Dynamic"
      ##Dynamic Model (rocof constraint)
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

      function dynamic1_callback_function(cb_data)
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

                #const_par = instance.units[g].max_power[T]/sum(instance.units[g].max_power[1] for g in 1:length(instance.units))
              temp_RT = sum(R_T[g]*callback_value(cb_data,model[:is_on][instance.units[g].name,t]) for g in 1:length(instance.units))-sum(R_T[ind_max[i]] for i in 1:n_cont)
              temp_FT = sum(F_T[g]*callback_value(cb_data,model[:is_on][instance.units[g].name,t]) for g in 1:length(instance.units))-sum(F_T[ind_max[i]] for i in 1:n_cont)
              temp_HT = sum(H_T[g]*callback_value(cb_data,model[:is_on][instance.units[g].name,t]) for g in 1:length(instance.units))-sum(H_T[ind_max[i]] for i in 1:n_cont)
              temp_DT = sum(D_T[g]*callback_value(cb_data,model[:is_on][instance.units[g].name,t]) for g in 1:length(instance.units))-sum(D_T[ind_max[i]] for i in 1:n_cont)

              omega_n = sqrt((temp_DT+temp_RT)/(2*temp_HT*DData.T_r))
              zeta = (2*temp_HT+DData.T_r*(temp_DT+temp_FT))/(2*sqrt(2*DData.T_r*temp_HT*(temp_DT+temp_RT)))
              omega_d = omega_n*sqrt(1-zeta^2)
              t_max = (1/omega_d) *atan(omega_d/(zeta*omega_n-1/DData.T_r))

              viol_nadir =  60*sum(pow_max)/(temp_DT+temp_RT)*(1+sqrt(DData.T_r*(temp_RT-temp_FT)/(2*temp_HT))*exp(-zeta*omega_n*t_max)) - 0.6
              #println("Nadir violation $viol_nadir, pow_max $pow_max")
              RD = [(DData.D[g]+DData.K[g]/DData.R[g])*DData.pu[g] for g in 1:length(instance.units)]
              RD_loss = sum(RD[ind_max[i]] for i in 1:n_cont)

              viol_rocof = -(temp_HT-60*sum(pow_max))

              if viol_rocof > 1e-4
                  println("RoCoF violation $viol_rocof at $t. Max pow = $pow_max, index= $ind_max")
                  con = @build_constraint(sum(H_T[g]*model[:is_on][instance.units[g].name,t] for g in 1:length(instance.units))
                                              - sum(H_T[ind_max[i]]*model[:is_on][instance.units[ind_max[i]].name,t] for i in 1:n_cont)
                                              - 60*sum(temp_prod[t,ind_max[i]] for i in 1:n_cont)./P_tot >= 0)

                  MOI.submit(model, MOI.LazyConstraint(cb_data), con)
              end

              viol_ss = -(sum(callback_value(cb_data,model[:is_on][instance.units[g].name,t])*RD[g] for g in 1:length(instance.units))
              - sum(RD[ind_max[i]] for i in 1:n_cont)
               - sum(pow_max)*60/0.3)

              if viol_ss > 1e-4
                  println("SS violation $viol_ss at $t. Max pow = $pow_max, index= $ind_max")
                  con = @build_constraint(sum(RD[g]*model[:is_on][instance.units[g].name,t] for g in 1:length(instance.units))
                                              - sum(RD[ind_max[i]]*model[:is_on][instance.units[ind_max[i]].name,t] for i in 1:n_cont)
                                              - 60/0.3*sum(temp_prod[t,ind_max[i]] for i in 1:n_cont)./P_tot >= 0)

                  MOI.submit(model, MOI.LazyConstraint(cb_data), con)

              end

              if viol_nadir > 1e-4
                 #Dumb way
                 # temp_on,temp_off = zeros((length(instance.units)))
                 #  for g in 1:length(instance.units)
                 #      if callback_value(cb_data,model[:is_on][instance.units[g].name,t]) == 0
                 #          temp_off[g]=1
                 #      else
                 #          temp_on[g]=1
                 #      end
                 #  end


                  #Another way
                  model1 = JuMP.Model(Ipopt.Optimizer)

                  @variable(model1, c[i=1:4])
                  @variable(model1, b)
                  const_term = (temp_DT+temp_RT)/(1+sqrt(DData.T_r*(temp_RT-temp_FT)/(2*temp_HT))*exp(-zeta*omega_n*t_max))
                  @NLobjective(model1, Min, (c[1]*temp_HT+c[2]*temp_RT+c[3]*temp_FT+c[4]*temp_DT+b-const_term)^2)
                  JuMP.optimize!(model1)
                  scale_H = value(model1[:c][1])
                  scale_R = value(model1[:c][2])
                  scale_F = value(model1[:c][3])
                  scale_D = value(model1[:c][4])
                  const_val = value(model1[:b])

                  #println("max_power = $pow_max, Constant term $const_term, Nadir violation $viol_nadir, Reconstructed $reconstruct")
                  con = @build_constraint((scale_H*(sum(H_T[g]*model[:is_on][instance.units[g].name,t] for g in 1:length(instance.units))-
                                   sum(H_T[ind_max[i]]*model[:is_on][instance.units[ind_max[i]].name,t] for i in 1:n_cont))+
                                   scale_R*(sum(R_T[g]*model[:is_on][instance.units[g].name,t] for g in 1:length(instance.units))-
                                    sum(R_T[ind_max[i]]*model[:is_on][instance.units[ind_max[i]].name,t] for i in 1:n_cont))+
                                    scale_F*(sum(F_T[g]*model[:is_on][instance.units[g].name,t] for g in 1:length(instance.units))-
                                     sum(F_T[ind_max[i]]*model[:is_on][instance.units[ind_max[i]].name,t] for i in 1:n_cont))+
                                     scale_D*(sum(D_T[g]*model[:is_on][instance.units[g].name,t] for g in 1:length(instance.units))-
                                      sum(D_T[ind_max[i]]*model[:is_on][instance.units[ind_max[i]].name,t] for i in 1:n_cont))+const_val)*0.6/60
                                      - sum(model[:temp_prod][t,ind_max[i]] for i in 1:n_cont)./P_tot >= 0)
                 MOI.submit(model, MOI.LazyConstraint(cb_data), con)

              end
          end
      end

      MOI.set(model, MOI.LazyConstraintCallback(), dynamic1_callback_function)
  end


end
