using ForwardDiff
using LinearAlgebra
using QuasiMonteCarlo
using UnitCommitment
using Random
Random.seed!(1234)
rand(376)

function minus_g(HT, RT, FT, DT; T_r)
    omega_n = sqrt((DT + RT) / (2 * HT * T_r))
    zeta = (2 * HT + (T_r * (DT + FT))) / (2 * sqrt(2 * T_r * HT * (DT + RT)))
    omega_d = omega_n * sqrt(1 - zeta^2)
    t_max = (1 / omega_d) * atan(omega_d / ((zeta * omega_n) - (1 / T_r)))

    return -(DT + RT) / (1 + (sqrt(T_r * (RT - FT) / (2 * HT)) * exp(-zeta * omega_n * t_max)))
end

function check_convexity_qmc(instance, DData)
    thermal_units = instance.scenarios[1].thermal_units
    H_T = DData.H .* DData.pu
    D_T = [(DData.D[g]) * DData.pu[g] for g in 1:length(thermal_units)]
    R_T = [(DData.K[g] / DData.R[g]) * DData.pu[g] for g in 1:length(thermal_units)]
    F_T = [(DData.F[g] * DData.K[g] / DData.R[g]) * DData.pu[g] for g in 1:length(thermal_units)]

    f(x) = minus_g(x[1], x[2], x[3], x[4]; T_r=DData.T_r)

    for g in eachindex(thermal_units)
        x = [sum(H_T[1:g]), sum(R_T[1:g]), sum(F_T[1:g]), sum(D_T[1:g])]
        H = Symmetric(ForwardDiff.hessian(f, x))
        @show(eigvals(H))
    end


    # lb = [0.02, 0.07, 0.009, 0.005]
    # ub = [6.01, 16.10, 3.26, 0.6]
    # n = 1000
    # samples = QuasiMonteCarlo.sample(n, lb, ub, SobolSample())
    # for j in 1:2
    #     x = samples[:, j]
    #     H = Symmetric(ForwardDiff.hessian(f, x))
    #     @show(H)
    #     eigvals(H)
    # end
end

instance = UnitCommitment.read_benchmark("matpower/case118/2017-02-01")
DData = freq_initialize(instance)
check_convexity_qmc(instance, DData)
