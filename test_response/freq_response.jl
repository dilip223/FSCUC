function freq_test(instance::UnitCommitmentInstance,
                  model::JuMP.Model,
                  f_max::Vector{Float64},
                  rocof::Vector{Float64},
                  steady_freq::Vector{Float64},
                  DData::DynamicsData,
                  n_cont::Int,
                  ESS_ON::Int)

      T = instance.time
      P_tot = sum(instance.units[g].max_power[1] for g in 1:length(instance.units))
      if n_cont == 0
         n_cont = n_cont+1
      end
      for t in 1:T
         R_T=0
         F_T=0
         H_T=0
         D_T=0

         pow_mat = [(instance.units[g].min_power[1]+value(model[:prod_above][instance.units[g].name,t])).*value(model[:is_on][instance.units[g].name,t]) for g in 1:length(instance.units)]
         pow_max = sort(pow_mat,rev=true)[1:n_cont]
         x=sum(pow_mat)
         ind_max = sortperm(pow_mat,rev=true)[1:n_cont]
         gen_lost_name = [instance.units[ind_max[g]].name for g in 1:n_cont]
         println("Lost generators = $gen_lost_name Lost Power = $pow_max, Total Power =$x, Time = $t")
         pow_max = sum(pow_max) ./ P_tot
         temp = [value(model[:is_on][instance.units[g].name,t])
                              for g in 1:length(instance.units)]
         temp[ind_max].= 0

         for g in 1:length(instance.units)
              on_status = temp[g]
              #const_par = instance.units[g].max_power[T]/sum(instance.units[g].max_power[1] for g in 1:length(instance.units))
              R_T = R_T+(DData.K[g]/DData.R[g])*on_status*DData.pu[g]
              F_T = F_T+(DData.F[g]*DData.K[g]/DData.R[g])*on_status*DData.pu[g]
              H_T = H_T+ DData.H[g]*on_status*DData.pu[g]
              D_T = D_T+ DData.D[g]*on_status*DData.pu[g]
         end
         omega_n = sqrt((D_T+R_T)/(2*H_T*DData.T_r))
         zeta = (2*H_T+DData.T_r*(D_T+F_T))/(2*sqrt(2*DData.T_r*H_T*(D_T+R_T)))

         omega_d = omega_n*sqrt(1-zeta^2)
         t_max = (1/omega_d) *atan(omega_d/(zeta*omega_n-1/DData.T_r))
         f_max[t] = pow_max/(D_T+R_T)*(1+sqrt(DData.T_r*(R_T-F_T)/(2*H_T))*
                     exp(-zeta*omega_n*t_max))
         rocof[t] = pow_max/(2*H_T)
         steady_freq[t] = pow_max/(D_T+R_T)
         
         ##time series plot
         # phi = atan(zeta/(sqrt(1 - zeta^2)));
         # temp_l = range(0.1, 30, 2000)
         # delf_dev = zeros((length(temp_l),1))
         # rocof_dev = zeros((length(temp_l),1))
         # for i in 1:length(temp_l)
         #     delf_dev[i] = pow_max/(2*H_T*DData.T_r*omega_n^2)*(1-(1/sqrt(1-zeta^2))*
         #     exp(-zeta*omega_n*temp_l[i])*cos(omega_n*sqrt(1-zeta^2)*temp_l[i]-phi))+
         #     pow_max/(2*H_T*omega_d)*exp(-zeta*omega_n*temp_l[i])*sin(omega_n*sqrt(1-zeta^2)*temp_l[i])
         # end
         # for i=1:length(temp_l)
         #    rocof_dev[i] = pow_max/(2*H_T*DData.T_r*omega_d)*exp(-zeta*omega_n*temp_l[i])*sin(omega_d*temp_l[i])+
         #       pow_max/(2*H_T*sqrt(1-zeta^2))*exp(-zeta*omega_n*temp_l[i])*cos(omega_d*temp_l[i]+phi);
         # end
         # return rocof_dev,temp_l
       end


end
