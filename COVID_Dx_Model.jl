using SIAN
using SIAN.Logging
using SIAN.Nemo
using SIAN.OrderedCollections
@info "Setting up the problem"

# k_bds = k_cas13 #nM-1 min-1
# k_SSS = k_FSS #min-1
# k_degRrep = k_degv  #nM-1 min-1

# k_RTon = (3/125) #0.024 #nM-1 min-1
# k_RToff = (24/10) #2.4 #min-1
# k_T7on = (84/25) #3.36 #nM-1 min-1
# k_T7off = 12 #min-1
# k_RNaseon = (3/125) #0.024 #nM-1 min-1
# k_RNaseoff = (24/10) #2.4 #min-1


ode = @ODEmodel(
    x_v'(t) = k_degv*x_v*x_Cas13 - k_cas13*x_v*x_u - k_cas13*x_v*x_p1,
    x_p1'(t) = - k_cas13*x_v*x_p1- k_cas13*x_p1*x_cDNA2,
    x_p2'(t) = - k_cas13*x_u*x_p2 - k_cas13*x_p2*x_cDNA1,
    x_p1v'(t) = k_cas13*x_v*x_p1 - k_degv*x_p1v*x_Cas13 - (3/125)*x_p1v*x_RT + (24/10)*x_RTp1v,
    x_p2u'(t) = k_cas13*x_u*x_p2 - k_degv*x_p2u*x_Cas13 - (3/125)*x_p2u*x_RT + (24/10)*x_RTp2u,
    x_p1cv'(t) = k_degv*x_p1v*x_Cas13 - (3/125)*x_p1cv*x_RT + (24/10)*x_RTp1cv,
    x_p2cu'(t) = k_degv*x_p2u*x_Cas13 - (3/125)*x_p2cu*x_RT + (24/10)*x_RTp2cu,
    x_RT'(t) = (24/10)*x_RTp1v + (24/10)*x_RTp1cv + (24/10)*x_RTp2cDNA1 + (24/10)*x_RTp2u + 
        (24/10)*x_RTp2cu + (24/10)*x_RTp1cDNA2 - (3/125)*x_RT*x_p1v - (3/125)*x_RT*x_p1cv - 
        (3/125)*x_RT*x_p2cDNA1 - (3/125)*x_RT*x_p2u - (3/125)*x_RT*x_p2cu - (3/125)*x_RT*x_p1cDNA2 +
        k_FSS*x_RTp1v + k_FSS*x_RTp2u + k_FSS*x_RTp2cDNA1 + k_FSS*x_RTp1cDNA2,
    x_RNase'(t) = (24/10)*x_RNasecDNA1v + (24/10)*x_RNasecDNA2u- (3/125)*x_RNase*x_cDNA1v -
        (3/125)*x_RNase*x_cDNA2u + k_RHA*x_RNasecDNA1v + k_RHA*x_RNasecDNA2u,
    x_RTp1v'(t) = - (24/10)*x_RTp1v + (3/125)*x_RT*x_p1v - k_degv*x_RTp1v*x_Cas13 - k_FSS*x_RTp1v,
    x_RTp2u'(t) = - (24/10)*x_RTp2u + (3/125)*x_RT*x_p2u - k_degv*x_RTp2u*x_Cas13 - k_FSS*x_RTp2u,
    x_RTp1cv'(t) = - (24/10)*x_RTp1cv + (3/125)*x_RT*x_p1cv + k_degv*x_RTp1v*x_Cas13,
    x_RTp2cu'(t) = - (24/10)*x_RTp2cu + (3/125)*x_RT*x_p2cu + k_degv*x_RTp2u*x_Cas13,
    x_cDNA1v'(t) = k_FSS*x_RTp1v - (3/125)*x_cDNA1v*x_RNase + (24/10)*x_RNasecDNA1v,
    x_cDNA2u'(t) = k_FSS*x_RTp2u - (3/125)*x_cDNA2u*x_RNase + (24/10) *x_RNasecDNA2u,
    x_RNasecDNA1v'(t) = - k_RHA*x_RNasecDNA1v - (24/10)*x_RNasecDNA1v + 
        (3/125)*x_RNase*x_cDNA1v,
    x_RNasecDNA2u'(t) = - k_RHA*x_RNasecDNA2u - (24/10)*x_RNasecDNA2u + 
        (3/125)*x_RNase*x_cDNA2u,
    x_cDNA1'(t) = k_RHA*x_RNasecDNA1v - k_cas13*x_cDNA1*x_p2,
    x_cDNA2'(t) = k_RHA*x_RNasecDNA2u - k_cas13*x_cDNA2*x_p1,
    x_p2cDNA1'(t) = k_cas13*x_cDNA1*x_p2 + (24/10)*x_RTp2cDNA1 - (3/125)*x_RT*x_p2cDNA1,
    x_p1cDNA2'(t) = k_cas13*x_cDNA2*x_p1 + (24/10)*x_RTp1cDNA2 - (3/125)*x_RT*x_p1cDNA2,
    x_RTp2cDNA1'(t) = (3/125)*x_RT*x_p2cDNA1 - (24/10)*x_RTp2cDNA1 - k_FSS*x_RTp2cDNA1,
    x_RTp1cDNA2'(t) = (3/125)*x_RT*x_p1cDNA2 - (24/10)*x_RTp1cDNA2 - k_FSS*x_RTp1cDNA2,
    x_T7'(t) = + 12*x_T7pro - (84/25)*x_T7*x_pro + k_txn*x_T7pro,
    x_pro'(t) = k_FSS*x_RTp2cDNA1 + k_FSS*x_RTp1cDNA2 - (84/25)*x_T7*x_pro + 
        12*x_T7pro + k_txn*x_T7pro,
    x_T7pro'(t) = - 12*x_T7pro + (84/25)*x_T7*x_pro - k_txn*x_T7pro,
    x_u'(t) = k_txn*x_T7pro - k_cas13*x_u*x_v - k_degv*x_u*x_Cas13 - k_cas13*x_u*x_iCas13 -
        k_cas13*x_u*x_p2,
    x_iCas13'(t) = - k_cas13*x_u*x_iCas13,
    x_Cas13'(t) = k_cas13*x_u*x_iCas13,
    x_uv'(t) = k_cas13*x_u*x_v,
    x_qRf'(t) = - k_degv*x_Cas13*x_qRf,
    x_q'(t) = k_degv*x_Cas13*x_qRf,
    x_f'(t) = k_degv*x_Cas13*x_qRf,
    y1(t) = x_f(t)
)

@time res = identifiability_ode(ode, get_parameters(ode); p=0.99, p_mod=0, local_only=true)

for (k, v) in res
    println("$k = $v")
end