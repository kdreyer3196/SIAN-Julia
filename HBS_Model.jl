using SIAN
using Logging
using Nemo
using OrderedCollections
using Dates
println(Dates.format(now(), "HH:MM"))
@info "Setting up the problem"

ode = @ODEmodel(
    x0'(t) = k_txn2 - (27/10)*x0,
    x1'(t) = x0 - (7/20)*x1 - k_dHAF*x1 - (k_bHS/(38/5))*x1*x3,
    x2'(t) = 1 - (27/10)*x2,
    x3'(t) = x2 - (7/20)*x3 - (k_bHS/(38/5))*x1*x3,
    x4'(t) = (k_bHS/(38/5))*x1*x3 - (7/20)*x4 - k_bHH*x9*x4,
    x5'(t) = k_txnH*(x7 + x10) - (27/10)*x5,
    x6'(t) = 1 - (27/10)*x6  - k_dH1R*x5*x6,
    x7'(t) = x6 - (7/20)*x7 - k_dHP*(38/5)*x7 - k_dH1P*x7*(x1 + x4),
    x8'(t) = 1 + k_txnBH*(x7 + x10) - (27/10)*x8,
    x9'(t) = x8 - (7/20)*x9 - k_dHP*(38/5)*x9 - k_bHH*x9*x4,
    x10'(t) = k_bHH*x9*x4 - (7/20)*x10,
    x11'(t) = k_txnBH*(x7 + x10) - (27/10)*x11,
    x12'(t) = x11 - (29/1000)*x12,
    y1(t) = x12(t)
)

res = identifiability_ode(ode, get_parameters(ode); p=0.99, p_mod=0, local_only=true)

println(res, Dates.format(now(), "HH:MM"))

