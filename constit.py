from numpy import linspace, sin, cos, pi, zeros, outer, array, dot
from numpy import trunc, arctan, eye, trace, nan_to_num, tensordot
from numpy import sqrt, abs, ones
from numpy.linalg import norm

def elastic(MP,P,G):
#    MP.p = -trace(MP.stress)
#    s_ij = -MP.stress - eye(3)*MP.p
#    MP.q = sqrt(3.*sum(sum(s_ij*s_ij))/2.)
    de_kk = trace(MP.dstrain)
    de_ij = MP.dstrain - de_kk*eye(3)/3.
    MP.dstress = P.S.K*de_kk*eye(3) + 2.*P.S.G*de_ij
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

def von_mises(MP,P,G):
    de_kk = trace(MP.dstrain) # scalar
    de_ij = MP.dstrain - de_kk*eye(3)/3. # matrix
    dsigma_kk = 3.*P.S.K*de_kk # scalar
#            dev_work_norm = norm(tensordot(MP.dev_stress,de_ij)) # scalar
    dev_work_norm = abs(sum(sum(MP.dev_stress*de_ij))) # scalar
    dev_stress_norm = sqrt(sum(sum(MP.dev_stress**2))) # scalar
    MP.dev_dstress = 2.*P.S.G*(de_ij - dev_work_norm/(2.*(P.S.k**2))*
                  ((dev_stress_norm/(sqrt(2.)*P.S.k))**(P.S.s-2.))*MP.dev_stress) # matrix
    MP.sigma_kk += dsigma_kk # scalar
    MP.dev_stress += MP.dev_dstress # matrix
    MP.dstress = (MP.dev_stress + MP.sigma_kk*eye(3)/3.) - MP.stress # matrix
    MP.stress = MP.dev_stress + MP.sigma_kk*eye(3)/3. # matrix
    MP.yieldfunction = dev_stress_norm/(sqrt(2.)*P.S.k) - 1. # scalar
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

def dp(MP,P,G):
#    if P.t < .1: # isotropic loading
#        dstrain = eye(3)*1e-4*P.dt
#    elif P.t < .2: # shear
#        dstrain = 3e-8*ones((3,3))
#        dstrain = dstrain - eye(3)*trace(dstrain)/3.
#    else: # unload shear
#        dstrain = -2e-8*ones((3,3))
#        dstrain = dstrain - eye(3)*trace(dstrain)/3.
    dstrain = -MP.dstrain
    strain = -MP.strain
    stress = -MP.stress
    
    de_kk = trace(dstrain)
    de_ij = dstrain - de_kk*eye(3)/3.
    MP.p = trace(stress)/3.
    s_ij = stress - eye(3)*MP.p
    MP.q = sqrt(3.*sum(sum(s_ij*s_ij))/2.)
    
    MP.y = MP.q/(P.S.beta*MP.p + (P.S.mu-P.S.beta)*MP.p) - 1.
#    if MP.y > 0.:
#        print 'WARNING: y is: ' + str(MP.y)
    lambda_2 = ((3.*P.S.G*MP.p*sum(sum(s_ij*de_ij)) - P.S.K*MP.q**2.*de_kk)/
                (3.*P.S.G*P.S.mu*MP.p**2 + P.S.K*P.S.beta*MP.q**2))
    Gamma_2 = lambda_2*(lambda_2>0)
    dstress = (2.*P.S.G*(de_ij - 3./2.*s_ij/MP.q*Gamma_2*(MP.q/(P.S.mu*MP.p))**(P.S.s-1.)) +
               P.S.K*eye(3)*(de_kk + P.S.beta*Gamma_2*(MP.q/(P.S.mu*MP.p))**P.S.s))
               
    MP.dstress = -dstress
    MP.stress += MP.dstress
        
def inviscid(MP,P,G):
    MP.dstress = 0
    MP.stress = -P.F.P_0*eye(3)

def viscous(MP,P,G):
    MP.stress = 2.*(P.F.mu*MP.dstrain/P.dt
                   - P.F.mu*trace(MP.dstrain/P.dt)*eye(3)/3.
                   - P.F.P*eye(3)/2.)

def compressible(MP,P,G):
    MP.rho = MP.rho/(1 + trace(MP.dstrain))
    MP.V = MP.m/MP.rho
    P_hat = MP.eq_of_state(P) # pressure
    MP.strain_rate = MP.dstrain/P.dt
    MP.stress = 2.*(P.F.mu*MP.strain_rate
                   - P.F.mu*trace(MP.strain_rate)*eye(3)/3.
                   - P_hat*eye(3)/2.)
    MP.int_energy += trace(MP.stress.T*MP.dstrain)/MP.rho

