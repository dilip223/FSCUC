function freq_test(instance::UnitCommitmentInstance,
                  model::JuMP.Model,
                  f_max::Vector{Float64},
                  rocof::Vector{Float64},
                  steady_freq::Vector{Float64},
                  DData::DynamicsData,
                  n_cont::Int)

      T = instance.time
      for t in 1:T
         R_T=0
         F_T=0
         H_T=0
         D_T=0

         pow_mat = [(instance.units[g].min_power[1]+value(model[:prod_above][instance.units[g].name,t])).*value(model[:is_on][instance.units[g].name,t]) for g in 1:length(instance.units)]
         pow_max = sort(pow_mat,rev=true)[1:n_cont]

         ind_max = sortperm(pow_mat,rev=true)[1:n_cont]

         pow_max = sum(pow_max)./sum(instance.units[g].max_power[1]
                                  for g in 1:length(instance.units))
         temp = [value(model[:is_on][instance.units[g].name,t])
                              for g in 1:length(instance.units)]
         temp[ind_max].= 0
         println("pow_max= $pow_max at index= $ind_max")
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

       end
end
