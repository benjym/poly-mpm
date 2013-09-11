from numpy import linspace, sin, cos, pi, zeros, outer, array, dot
from numpy import trunc, arctan, eye, trace, nan_to_num, tensordot
from numpy import sqrt, abs, ones, minimum, maximum
from numpy.linalg import norm

def rigid(MP,P,G,p):
    pass

def elastic(MP,P,G,p):
#    MP.p = -trace(MP.stress)
#    s_ij = -MP.stress - eye(3)*MP.p
#    MP.q = sqrt(3.*sum(sum(s_ij*s_ij))/2.)
    de_kk = trace(MP.dstrain)
    de_ij = MP.dstrain - de_kk*eye(3)/3.
    MP.dstress = P.S[p].K*de_kk*eye(3) + 2.*P.S[p].G*de_ij
    MP.stress += MP.dstress
    MP.gammadot = sqrt(sum(sum(de_ij**2)))
    MP.pressure = trace(MP.stress)/3.
    MP.sigmav = MP.stress[1,1]
    MP.sigmah = MP.stress[0,0]
    for r in xrange(4):
        n = G.nearby_nodes(MP.n_star,r,P)
        G.pressure[n] += MP.N[r]*MP.pressure
        G.gammadot[n] += MP.N[r]*MP.gammadot/P.dt
        G.sigmav[n] += MP.N[r]*MP.sigmav
        G.sigmah[n] += MP.N[r]*MP.sigmah

def von_mises(MP,P,G,p):
    de_kk = trace(MP.dstrain) # scalar
    de_ij = MP.dstrain - de_kk*eye(3)/3. # matrix
    dsigma_kk = 3.*P.S[p].K*de_kk # scalar
#            dev_work_norm = norm(tensordot(MP.dev_stress,de_ij)) # scalar
    dev_work_norm = abs(sum(sum(MP.dev_stress*de_ij))) # scalar
    dev_stress_norm = sqrt(sum(sum(MP.dev_stress**2))) # scalar
    MP.dev_dstress = 2.*P.S[p].G*(de_ij - dev_work_norm/(2.*(P.S[p].k**2))*
                  ((dev_stress_norm/(sqrt(2.)*P.S[p].k))**(P.S[p].s-2.))*MP.dev_stress) # matrix
    MP.sigma_kk += dsigma_kk # scalar
    MP.dev_stress += MP.dev_dstress # matrix
    MP.dstress = (MP.dev_stress + MP.sigma_kk*eye(3)/3.) - MP.stress # matrix
    MP.stress = MP.dev_stress + MP.sigma_kk*eye(3)/3. # matrix
    MP.yieldfunction = dev_stress_norm/(sqrt(2.)*P.S[p].k) - 1. # scalar
    MP.gammadot = sqrt(sum(sum(de_ij**2)))
    MP.pressure = trace(MP.stress)/3.
    MP.sigmav = MP.stress[1,1]
    MP.sigmah = MP.stress[0,0]

    for r in xrange(4):
        n = G.nearby_nodes(MP.n_star,r,P)
        G.pressure[n] += MP.N[r]*MP.pressure*MP.m
        G.yieldfunction[n] += MP.N[r]*MP.yieldfunction*MP.m
        G.gammadot[n] += MP.N[r]*MP.gammadot/P.dt*MP.m
        G.sigmav[n] += MP.N[r]*MP.sigmav*MP.m
        G.sigmah[n] += MP.N[r]*MP.sigmah*MP.m
        G.dev_stress[n] += MP.N[r]*norm(MP.dev_stress)*MP.m
        G.dev_stress_dot[n] += MP.N[r]*norm(MP.dev_dstress)/P.dt*MP.m

def dp(MP,P,G,p):
    dstrain = -MP.dstrain
    strain = -MP.strain
    stress = -MP.stress
    
    de_kk = trace(dstrain)
    de_ij = dstrain - de_kk*eye(3)/3.
    MP.p = trace(stress)/3.
    s_ij = stress - eye(3)*MP.p
    MP.q = sqrt(3.*sum(sum(s_ij*s_ij))/2.)
    
    MP.y = MP.q/(P.S[p].beta*MP.p + (P.S[p].mu-P.S[p].beta)*MP.p) - 1.
#    if MP.y > 0.:
#        print 'WARNING: y is: ' + str(MP.y)
    lambda_2 = ((3.*P.S[p].G*MP.p*sum(sum(s_ij*de_ij)) - P.S[p].K*MP.q**2.*de_kk)/
                (3.*P.S[p].G*P.S[p].mu*MP.p**2 + P.S[p].K*P.S[p].beta*MP.q**2))
    Gamma_2 = lambda_2*(lambda_2>0)
    dstress = (2.*P.S[p].G*(de_ij - 3./2.*s_ij/MP.q*Gamma_2*(MP.q/(P.S[p].mu*MP.p))**(P.S[p].s-1.)) +
               P.S[p].K*eye(3)*(de_kk + P.S[p].beta*Gamma_2*(MP.q/(P.S[p].mu*MP.p))**P.S[p].s))
               
    MP.dstress = -dstress
    MP.stress += MP.dstress
#     MP.yieldfunction = 
    MP.gammadot = sqrt(sum(sum(de_ij**2)))
    MP.pressure = trace(MP.stress)/3.
    MP.dev_stress = MP.stress - eye(3)*MP.pressure
    MP.dev_stress_dot = MP.dstress - eye(3)*trace(MP.dstress)/3.
    MP.sigmav = MP.stress[1,1]
    MP.sigmah = MP.stress[0,0]
    
    for r in xrange(4):
        n = G.nearby_nodes(MP.n_star,r,P)
        G.pressure[n] += MP.N[r]*MP.pressure*MP.m
#         G.yieldfunction[n] += MP.N[r]*MP.yieldfunction*MP.m
        G.gammadot[n] += MP.N[r]*MP.gammadot/P.dt*MP.m
        G.sigmav[n] += MP.N[r]*MP.sigmav*MP.m
        G.sigmah[n] += MP.N[r]*MP.sigmah*MP.m
        G.dev_stress[n] += MP.N[r]*norm(MP.dev_stress)*MP.m
        G.dev_stress_dot[n] += MP.N[r]*norm(MP.dev_stress_dot)/P.dt*MP.m

def viscous(MP,P,G,p): # VALIDATED!
    MP.de_kk = trace(MP.dstrain)/3.
    MP.de_ij = MP.dstrain - MP.de_kk*eye(3)
    MP.dp = P.S[p].K*MP.de_kk

    MP.pressure += MP.dp
    MP.dev_stress = 2.*P.S[p].mu_s*MP.de_ij/P.dt

    MP.stress = MP.pressure*eye(3) + MP.dev_stress - P.S[p].mu_v*MP.de_kk*eye(3)
    MP.gammadot = sqrt(sum(sum(MP.de_ij**2)))/P.dt
    
    for r in xrange(4):
        n = G.nearby_nodes(MP.n_star,r,P)
        G.pressure[n] += MP.N[r]*MP.pressure*MP.m
        G.gammadot[n] += MP.N[r]*MP.gammadot/P.dt*MP.m
        G.dev_stress[n] += MP.N[r]*norm(MP.dev_stress)*MP.m

def bingham(MP,P,G,p): # UNDER VALIDATION
    MP.de_kk = trace(MP.dstrain)/3.
    MP.de_ij = MP.dstrain - MP.de_kk*eye(3)
    MP.gammadot = sqrt(sum(sum(MP.de_ij**2)))/P.dt
    MP.dp = P.S[p].K*MP.de_kk

    MP.pressure += MP.dp
    if MP.gammadot <= P.S[p].gamma_c:
        MP.dev_stress = 2.*P.S[p].mu_0*MP.de_ij/P.dt
        isyielded = -1.
    else:
        MP.dev_stress = 2.*(P.S[p].mu_s + P.S[p].tau_0/MP.gammadot)*MP.de_ij/P.dt
        isyielded = 0.

    MP.stress = MP.pressure*eye(3) + MP.dev_stress

    for r in xrange(4):
        n = G.nearby_nodes(MP.n_star,r,P)
        G.pressure[n] += MP.N[r]*MP.pressure*MP.m
        G.gammadot[n] += MP.N[r]*MP.gammadot/P.dt*MP.m
        G.dev_stress[n] += MP.N[r]*norm(MP.dev_stress)*MP.m
        G.yieldfunction[n] += MP.N[r]*isyielded*MP.m

def pouliquen(MP,P,G,p): # UNVALIDATED
    MP.de_kk = trace(MP.dstrain)/3.
    MP.de_ij = MP.dstrain - MP.de_kk*eye(3)
    
    MP.pressure = trace(MP.stress)/3.
    MP.pressure += P.S[p].K*MP.de_kk

    MP.gammadot_ij = 2.*MP.de_ij/P.dt
    MP.gammadot = sqrt(0.5*sum(sum(MP.gammadot_ij**2))) # second invariant of strain rate tensor
    MP.gammadot = maximum(MP.gammadot,1e-20) # no zero
    d = 1. # MEAN DIAMETER
    I = sqrt(2)*d*MP.gammadot/sqrt(abs(MP.pressure)/MP.rho)
    mu = nan_to_num(P.S[p].mu_s + (P.S[p].mu_2 - P.S[p].mu_s)/nan_to_num(P.S[p].I_0/I + 1.))
    eta = maximum(mu*MP.pressure/(sqrt(2)*MP.gammadot),0)
#    eta_max = 250.*MP.rho # missing a g*H^3 !!
#    eta = minimum(eta,,eta_max)
    MP.dev_stress = eta*MP.gammadot_ij
    MP.stress = MP.pressure*eye(3) + MP.dev_stress

    for r in xrange(4):
        n = G.nearby_nodes(MP.n_star,r,P)
        G.pressure[n] += MP.N[r]*MP.pressure*MP.m
        G.gammadot[n] += MP.N[r]*MP.gammadot/P.dt*MP.m
        G.dev_stress[n] += MP.N[r]*norm(MP.dev_stress)*MP.m
