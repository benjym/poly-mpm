"""
A set of methods for describing all of the constitutive models.

Each model operates on an individual material point, updating its' properties. The relevant properties are also mapped to the grid for visualisation.

The models which currently work reliably are:
    - Rigid
    - Elastic
    - Von mises
    - Newtonian viscosity

"""
import sys
from numpy import linspace, sin, cos, pi, zeros, outer, array, dot
from numpy import trunc, arctan, eye, trace, nan_to_num, tensordot
from numpy import sqrt, abs, ones, minimum, maximum, exp, isfinite
from numpy.linalg import norm

def rigid(MP,P,G,p):
    """Perfectly rigid particles. No stresses are calculated.

    This material model has been validated and is functional.

    :param MP: A particle.Particle instance.
    :param P: A param.Param instance.
    :param G: A grid.Grid instance
    :param p: The particle number.
    :type p: int

    """
    pass

def elastic(MP,P,G,p):
    """Linear elasticity.

    This material model has been validated and is functional.

    :param MP: A particle.Particle instance.
    :param P: A param.Param instance.
    :param G: A grid.Grid instance
    :param p: The particle number.
    :type p: int

    """
    de_kk = trace(MP.dstrain)
    de_ij = MP.dstrain - de_kk*eye(3)/3.
    MP.dstress = P.S[p].K*de_kk*eye(3) + 2.*P.S[p].G*de_ij

    # for visualisation
    MP.gammadot = sqrt(sum(sum(de_ij**2)))
    MP.pressure = trace(MP.stress)/3.
    MP.sigmav = MP.stress[1,1]
    MP.sigmah = MP.stress[0,0]
    for r in range(4):
        n = G.nearby_nodes(MP.n_star,r,P)
        G.pressure[n] += MP.N[r]*MP.pressure*MP.m
        G.gammadot[n] += MP.N[r]*MP.gammadot/P.dt*MP.m
        G.sigmav[n] += MP.N[r]*MP.sigmav*MP.m
        G.sigmah[n] += MP.N[r]*MP.sigmah*MP.m

def elastic_time(MP,P,G,p):
    """Linear elasticity that is time dependent.

    This material model has been validated and is functional.

    :param MP: A particle.Particle instance.
    :param P: A param.Param instance.
    :param G: A grid.Grid instance
    :param p: The particle number.
    :type p: int

    """

    if (p == 0) and (P.t > P.S[p].t_0):
        if MP.m > P.S[p].rho_f*MP.V:
            MP.m *= 0.95

    K_curr = P.S[p].E_0/(3*(1-2*P.S[p].nu))
    G_curr = P.S[p].E_0/(2*(1+P.S[p].nu))

#     if p == 0:
#         if P.t < P.S[p].t_0: E = P.S[p].E_0
#         else: E = P.S[p].E_f + (P.S[p].E_0-P.S[p].E_f)*exp(-3.*(P.t-P.S[p].t_0)**2/P.S[p].t_c**2)
#     else: E = P.S[p].E_0
#     K_curr = E/(3*(1-2*P.S[p].nu))
#     G_curr = E/(2*(1+P.S[p].nu))


    de_kk = trace(MP.dstrain)
    de_ij = MP.dstrain - de_kk*eye(3)/3.
    MP.dstress = K_curr*de_kk*eye(3) + 2.*G_curr*de_ij

    # for visualisation
    MP.gammadot = sqrt(sum(sum(de_ij**2)))
    MP.pressure = -trace(MP.stress)/3.
    MP.sigmav = -MP.stress[1,1]
    MP.sigmah = -MP.stress[0,0]
    for r in range(4):
        n = G.nearby_nodes(MP.n_star,r,P)
        G.pressure[n] += MP.N[r]*MP.pressure*MP.m
        G.gammadot[n] += MP.N[r]*MP.gammadot/P.dt*MP.m
        G.sigmav[n] += MP.N[r]*MP.sigmav*MP.m
        G.sigmah[n] += MP.N[r]*MP.sigmah*MP.m


def von_mises(MP,P,G,p):
    """Von Mises :math:`H^2` yield criterion --- see http://dx.doi.org/10.1016/j.ijsolstr.2012.02.003 Eqs 7.11 - 7.15

    This material model has been validated and is functional.

    :param MP: A particle.Particle instance.
    :param P: A param.Param instance.
    :param G: A grid.Grid instance
    :param p: The particle number.
    :type p: int

    """
    de_kk = trace(MP.dstrain) # scalar
    de_ij = MP.dstrain - de_kk*eye(3)/3. # matrix
    dsigma_kk = 3.*P.S[p].K*de_kk # scalar
    dev_work_norm = abs(sum(sum(MP.dev_stress*de_ij))) # scalar
    dev_stress_norm = sqrt(sum(sum(MP.dev_stress**2))) # scalar
    MP.dev_dstress = 2.*P.S[p].G*(de_ij - dev_work_norm/(2.*(P.S[p].k**2))*
                  ((dev_stress_norm/(sqrt(2.)*P.S[p].k))**(P.S[p].s-2.))*MP.dev_stress) # matrix
    MP.sigma_kk += dsigma_kk # scalar
    MP.dev_stress += MP.dev_dstress # matrix
    MP.dstress = (MP.dev_stress + MP.sigma_kk*eye(3)/3.) - MP.stress # matrix

    # for visualisation
    MP.yieldfunction = dev_stress_norm/(sqrt(2.)*P.S[p].k) - 1. # scalar
    MP.gammadot = sqrt(sum(sum(de_ij**2)))
    MP.pressure = trace(MP.stress)/3.
    MP.sigmav = MP.stress[1,1]
    MP.sigmah = MP.stress[0,0]

    for r in range(4):
        n = G.nearby_nodes(MP.n_star,r,P)
        G.pressure[n] += MP.N[r]*MP.pressure*MP.m
        G.yieldfunction[n] += MP.N[r]*MP.yieldfunction*MP.m
        G.gammadot[n] += MP.N[r]*MP.gammadot/P.dt*MP.m
        G.sigmav[n] += MP.N[r]*MP.sigmav*MP.m
        G.sigmah[n] += MP.N[r]*MP.sigmah*MP.m
        G.dev_stress[n] += MP.N[r]*norm(MP.dev_stress)*MP.m
        G.dev_stress_dot[n] += MP.N[r]*norm(MP.dev_dstress)/P.dt*MP.m

def dp(MP,P,G,p): # UNVALIDATED
    """Drucker-Prager :math:`H^2` yield criterion.

    This material model is ONLY PARTIALLY WORKING - NOT YET VALIDATED!

    :param MP: A particle.Particle instance.
    :param P: A param.Param instance.
    :param G: A grid.Grid instance
    :param p: The particle number.
    :type p: int

    """
    def dp_guts(dstrain2d,stress2d):
        de_kk = trace(dstrain2d)
        de_ij = dstrain2d - de_kk*eye(2)/2. # shear strain
        MP.p = trace(stress2d)/2. # pressure, positive for compression
        s_ij = stress2d - eye(2)*MP.p # shear stress
        MP.q = sqrt(3.*sum(sum(s_ij*s_ij))/2.)

        if MP.p > 0: K = P.S[p].K
        else: K = 0

        lambda_2 = ((3.*P.S[p].G*MP.p*sum(sum(s_ij*de_ij)) - K*MP.q**2.*de_kk)/
                    (3.*P.S[p].G*P.S[p].mu*MP.p**2 + K*P.S[p].beta*MP.q**2))
        Gamma_2 = lambda_2*(lambda_2>0) # Macauley bracket
        dstress = (2.*P.S[p].G*(de_ij - 3./2.*s_ij/MP.q*Gamma_2*(MP.q/(P.S[p].mu*MP.p))**(P.S[p].s-1.)) +
                   K*eye(2)*(de_kk + P.S[p].beta*Gamma_2*(MP.q/(P.S[p].mu*MP.p))**P.S[p].s))
        return dstress

    dstrain = -MP.dstrain[:2,:2] # convert from fluid mechanics to soil mechanics convention, just 2d
    stress = -MP.stress[:2,:2]   # convert from fluid mechanics to soil mechanics convention, just 2d

    dstress1 = dp_guts(dstrain,stress)
    dstress2 = dp_guts(dstrain,stress + 0.5*dstress1)
    dstress3 = dp_guts(dstrain,stress + 0.5*dstress2)
    dstress4 = dp_guts(dstrain,stress + dstress3)
    dstress = (dstress1 + 2.*dstress2 + 2.*dstress3 + dstress4)/6.

#     dstress = dp_guts(dstrain,stress)

    MP.dstress[:2,:2] = -dstress

    # For visualisation
    MP.gammadot = sqrt(sum(sum((dstrain - trace(dstrain)*eye(2)/2.)**2)))
    MP.pressure = trace(MP.stress)/2. # 3
#     MP.dev_stress = MP.stress - eye(3)*MP.pressure
#     MP.dev_stress_dot = MP.dstress - eye(3)*trace(MP.dstress)/2. # 3
    MP.sigmav = MP.stress[1,1]
    MP.sigmah = MP.stress[0,0]
#     MP.yieldfunction = MP.q/(P.S[p].beta*MP.p - (P.S[p].mu-P.S[p].beta)*K*trace(strain)) - 1
    MP.yieldfunction = MP.q/(P.S[p].mu*MP.p) - 1. # Associated flow only!

#     MP.p = -(MP.stress[0,0]+MP.stress[1,1])/2. # pressure, positive for compression
#     s_ij = -MP.stress[:2,:2] - eye(2)*MP.p # shear stress
#     MP.q = sqrt(3.*sum(sum(s_ij*s_ij))/2.)

    for r in range(4):
        n = G.nearby_nodes(MP.n_star,r,P)
        G.pressure[n] += MP.N[r]*MP.m*MP.p
#         G.yieldfunction[n] += MP.N[r]*MP.yieldfunction*MP.m
        G.gammadot[n] += MP.N[r]*MP.gammadot/P.dt*MP.m
#         G.sigmav[n] += MP.N[r]*MP.sigmav*MP.m
#         G.sigmah[n] += MP.N[r]*MP.sigmah*MP.m
        G.dev_stress[n] += MP.N[r]*MP.q*MP.m
#         G.dev_stress_dot[n] += MP.N[r]*norm(MP.dev_stress_dot)/P.dt*MP.m

def dp_rate(MP,P,G,p): # UNVALIDATED
    """Drucker-Prager yield criterion with rate dependent behaviour.

    This material model is PROBALBLY WORKING - NOT YET VALIDATED!

    :param MP: A particle.Particle instance.
    :param P: A param.Param instance.
    :param G: A grid.Grid instance
    :param p: The particle number.
    :type p: int

    """
    def dp_guts(dstrain,stress):
        # calculate strain components
        de_kk = trace(dstrain)
        de_ij = dstrain - de_kk*eye(2)/2.

        # calculate stress components
        MP.p = trace(stress)/2.
        s_ij = stress - eye(2)*MP.p
        MP.q = sqrt(3.*sum(sum(s_ij*s_ij))/2.)

        # find plasticity multiplier
        lambda_2 = MP.q - P.S[0].mu*MP.p

        # take macaulay
        Gamma_2 = lambda_2*(lambda_2>0)

        dstress = (2.*P.S[0].G*(de_ij - 3.*s_ij/(2.*MP.q*MP.p*P.S[0].t_star)*Gamma_2) +
                   P.S[0].K*eye(2)*(de_kk + P.S[0].beta*MP.q/(P.S[0].mu*MP.p*MP.p*P.S[0].t_star)*Gamma_2))
        return dstress

    dstrain = -MP.dstrain[:2,:2] # convert from fluid mechanics to soil mechanics convention, just 2d
    stress = -MP.stress[:2,:2]   # convert from fluid mechanics to soil mechanics convention, just 2d

#     dstress1 = dp_guts(dstrain,stress)
#     dstress2 = dp_guts(dstrain,stress + 0.5*dstress1)
#     dstress3 = dp_guts(dstrain,stress + 0.5*dstress2)
#     dstress4 = dp_guts(dstrain,stress + dstress3)
#     dstress = (dstress1 + 2.*dstress2 + 2.*dstress3 + dstress4)/6.

    dstress = dp_guts(dstrain,stress)

    MP.dstress[:2,:2] = -dstress

#     print(MP.stress)
#     print(MP.p,MP.q)
    # For visualisation
    MP.gammadot = sqrt(sum(sum((dstrain - trace(dstrain)*eye(2)/2.)**2)))
    MP.pressure = trace(MP.stress)/2. # 3
#     MP.dev_stress = MP.stress - eye(3)*MP.pressure
#     MP.dev_stress_dot = MP.dstress - eye(3)*trace(MP.dstress)/2. # 3
    MP.sigmav = MP.stress[1,1]
    MP.sigmah = MP.stress[0,0]
#     MP.yieldfunction = MP.q/(P.S[p].beta*MP.p - (P.S[p].mu-P.S[p].beta)*K*trace(strain)) - 1
    MP.yieldfunction = MP.q/(P.S[p].mu*MP.p) - 1. # Associated flow only!

    MP.p = -(MP.stress[0,0]+MP.stress[1,1])/2. # pressure, positive for compression
    s_ij = -MP.stress[:2,:2] - eye(2)*MP.p # shear stress
    MP.q = sqrt(3.*sum(sum(s_ij*s_ij))/2.)

    for r in range(4):
        n = G.nearby_nodes(MP.n_star,r,P)
        G.pressure[n] += MP.N[r]*MP.m*MP.p
#         G.yieldfunction[n] += MP.N[r]*MP.yieldfunction*MP.m
        G.gammadot[n] += MP.N[r]*MP.gammadot/P.dt*MP.m
#         G.sigmav[n] += MP.N[r]*MP.sigmav*MP.m
#         G.sigmah[n] += MP.N[r]*MP.sigmah*MP.m
        G.dev_stress[n] += MP.N[r]*MP.q*MP.m
#         G.dev_stress_dot[n] += MP.N[r]*norm(MP.dev_stress_dot)/P.dt*MP.m

def viscous(MP,P,G,p):
    """Linear viscosity.

    This material model has been validated and is functional.

    :param MP: A particle.Particle instance.
    :param P: A param.Param instance.
    :param G: A grid.Grid instance
    :param p: The particle number.
    :type p: int

    """
    MP.de_kk = trace(MP.dstrain)/3.
    MP.de_ij = MP.dstrain - MP.de_kk*eye(3)
    MP.dp = P.S[p].K*MP.de_kk

    MP.pressure += MP.dp
    MP.dev_stress = 2.*P.S[p].mu_s*MP.de_ij/P.dt

    viscous_volumetric = P.S[p].mu_v*MP.de_kk*eye(3)/P.dt
    # MP.stress = MP.pressure*eye(3) + MP.dev_stress - viscous_volumetric
    MP.dstress = MP.pressure*eye(3) + MP.dev_stress - viscous_volumetric - MP.stress
    MP.gammadot = sqrt(sum(sum(MP.de_ij**2)))/P.dt

    for r in range(4):
        n = G.nearby_nodes(MP.n_star,r,P)
        G.pressure[n] += MP.N[r]*MP.pressure*MP.m
#         G.gammadot[n] += MP.N[r]*MP.gammadot*MP.m
        G.dev_stress[n] += MP.N[r]*norm(MP.dev_stress)*MP.m

def bingham(MP,P,G,p):
    """Bingham fluid.

    This material model has not been validated but may be functional. YMMV.

    :param MP: A particle.Particle instance.
    :param P: A param.Param instance.
    :param G: A grid.Grid instance
    :param p: The particle number.
    :type p: int

    """
    MP.de_kk = trace(MP.dstrain)/3.
    MP.de_ij = MP.dstrain - MP.de_kk*eye(3)
    MP.dp = P.S[p].K*MP.de_kk
    MP.gammadot = sqrt(sum(sum(MP.de_ij**2)))/P.dt
    MP.pressure += MP.dp
    if MP.gammadot <= P.S[p].gamma_c: MP.dev_stress = 2.*P.S[p].mu_0*MP.de_ij/P.dt
    else: MP.dev_stress = 2.*(P.S[p].mu_s + P.S[p].tau_0/MP.gammadot)*MP.de_ij/P.dt
    viscous_volumetric = P.S[p].mu_v*MP.de_kk*eye(3)/P.dt
    MP.stress = MP.pressure*eye(3) + MP.dev_stress - viscous_volumetric

    for r in range(4):
        n = G.nearby_nodes(MP.n_star,r,P)
        G.pressure[n] += MP.N[r]*MP.pressure*MP.m
#         G.gammadot[n] += MP.N[r]*MP.gammadot*MP.m
        G.dev_stress[n] += MP.N[r]*norm(MP.dev_stress)*MP.m

def ken_kamrin(MP,P,G,p):
    """:math:`\mu(I)` rheology implemented from Ken Kamrin's MPM paper. Definitely _NOT_ functional.

    :param MP: A particle.Particle instance.
    :param P: A param.Param instance.
    :param G: A grid.Grid instance
    :param p: The particle number.
    :type p: int

    """
    MP.de_kk = trace(MP.dstrain)/3.
    MP.de_ij = MP.dstrain - MP.de_kk*eye(3)

#     MP.pressure = trace(MP.stress)/3.
    MP.pressure += P.S[p].K*MP.de_kk

    MP.gammadot_ij = 2.*MP.de_ij/P.dt
    MP.gammadot = sqrt(0.5*sum(sum(MP.gammadot_ij**2))) # second invariant of strain rate tensor
    MP.gammadot = maximum(MP.gammadot,1e-5) # no zero
    d = 1. # MEAN DIAMETER
    I = sqrt(2)*d*MP.gammadot/sqrt(abs(MP.pressure)/MP.rho)
    mu = nan_to_num(P.S[p].mu_s + (P.S[p].mu_2 - P.S[p].mu_s)/(nan_to_num(P.S[p].I_0/I) + 1.))
    eta = nan_to_num(mu*MP.pressure/(sqrt(2)*MP.gammadot)) # does this make sense for negative pressures???
    eta = maximum(eta,0)
    eta_max = 250.*MP.rho*P.max_g*(P.G.y_M-P.G.y_m)**3 # missing a g*H^3 !!
    eta = minimum(eta,eta_max)
#    print MP.gammadot, mu, eta, eta_max

    MP.dev_stress = eta*MP.gammadot_ij
    viscous_volumetric = P.S[p].mu_v*MP.de_kk*eye(3)/P.dt
    # MP.stress = MP.pressure*eye(3) + MP.dev_stress - viscous_volumetric
    MP.dstress = MP.pressure*eye(3) + MP.dev_stress - viscous_volumetric - MP.stress

    for r in range(4):
        n = G.nearby_nodes(MP.n_star,r,P)
        G.pressure[n] += MP.N[r]*MP.pressure*MP.m
        G.gammadot[n] += MP.N[r]*MP.gammadot/P.dt*MP.m
        G.dev_stress[n] += MP.N[r]*norm(MP.dev_stress)*MP.m

def pouliquen(MP,P,G,p):
    """:math:`\mu(I)` rheology implemented only for flowing regime.

    This material model has been validated and is PROBABLY functional. (needs more mileage)

    :param MP: A particle.Particle instance.
    :param P: A param.Param instance.
    :param G: A grid.Grid instance
    :param p: The particle number.
    :type p: int

    """

    # Issues with the mu(I) model in general:
    # 1. gammadot = 0
    # 2. pressure = 0
    # 3. large eta
    # 4. no feedback to density

    s_bar = 0.
    for i in range(P.G.ns): s_bar += MP.phi[i]*P.G.s[i]

    MP.de_kk = trace(MP.dstrain)/3. # tension positive
    MP.de_ij = MP.dstrain - MP.de_kk*eye(3) # shear strain increment

    MP.gammadot = sqrt(sum(sum((2.*MP.de_ij/P.dt)**2))) # norm of shear strain rate

    MP.I = MP.gammadot*s_bar*sqrt(P.S[p].rho_s/abs(MP.pressure))
    MP.I = nan_to_num(MP.I)
    MP.mu = P.S[p].mu_0 + P.S[p].delta_mu/(P.S[p].I_0/MP.I + 1.)
    MP.eta = 2.*sqrt(2)*MP.mu*abs(MP.pressure)/MP.gammadot # HACK: 2*SQRT(2) FIXES ISSUES WITH DEFINITION OF STRAIN
    MP.eta_limited = minimum(nan_to_num(MP.eta),P.S[p].eta_max) # COPYING FROM HERE: http://www.lmm.jussieu.fr/~lagree/TEXTES/PDF/JFMcollapsePYLLSSP11.pdf
    MP.dev_stress = MP.eta_limited*MP.de_ij/P.dt

    MP.dp = P.S[p].K*MP.de_kk # tension positive # FIXME do I need to multiply this by 3??
    MP.pressure += MP.dp

    if (MP.pressure > 0) or (MP.rho < 2200): # can't go into tension - this is really important!!
        MP.dp -= MP.pressure # set increment back to zero
        MP.pressure = 0.

    # if MP.rho < 2200:
    #     MP.dstress = - MP.stress # ADDED BY BENJY - CANNOT SUPPORT LOAD IF DENSITY LESS THAN CUTOFF
    #     MP.pressure = 0.
    # else: MP.dstress = MP.pressure*eye(3) + MP.dev_stress - MP.stress

    MP.dstress = MP.pressure*eye(3) + MP.dev_stress - MP.stress


    if not isfinite(MP.dstress).all():
        print('THIS IS GOING TO BE A PROBLEM! FOUND SOMETHING NON-FINITE IN CONSTITUTIVE MODEL')
        print(MP.de_kk)
        print(MP.de_ij)
        print(s_bar)
        print(MP.gammadot)
        print(MP.I)
        print(MP.mu)
        print(MP.eta)
        print(MP.dev_stress)
        print(MP.dp)
        print(MP.pressure)
        print(MP.dstress)
        sys.exit()

    for r in range(4):
        n = G.nearby_nodes(MP.n_star,r,P)
        G.pressure[n] += MP.N[r]*MP.pressure*MP.m
        G.dev_stress[n] += MP.N[r]*norm(MP.dev_stress)/sqrt(2.)*MP.m
        G.mu[n] += MP.N[r]*MP.mu*MP.m
        G.I[n] += MP.N[r]*MP.I*MP.m
        G.eta[n] += MP.N[r]*MP.eta*MP.m

def linear_mu(MP,P,G,p):
    """:math:`\mu(I)` rheology implemented only for flowing regime.

    This material model has been validated and is PROBABLY functional. (needs more mileage)

    :param MP: A particle.Particle instance.
    :param P: A param.Param instance.
    :param G: A grid.Grid instance
    :param p: The particle number.
    :type p: int

    """
    s_bar = 0.
    for i in range(P.G.ns): s_bar += MP.phi[i]*P.G.s[i]

    MP.de_kk = trace(MP.dstrain)/3. # tension positive
    MP.de_ij = MP.dstrain - MP.de_kk*eye(3) # shear strain increment

    MP.gammadot = sqrt(sum(sum((2.*MP.de_ij/P.dt)**2))) # norm of shear strain rate

    MP.I = MP.gammadot*s_bar*sqrt(P.S[p].rho_s/abs(MP.pressure))
    # MP.I = minimum(MP.I,10.) # HACK - IF I DO THIS DO I NEED TO FIDDLE WITH ETA LATER?
    # MP.mu = P.S[p].mu_0 + P.S[p].delta_mu/(P.S[p].I_0/MP.I + 1.)
    MP.mu = P.S[p].mu_0 + P.S[p].b*MP.I

    MP.eta = 2.*sqrt(2)*MP.mu*abs(MP.pressure)/MP.gammadot # HACK: 2*SQRT(2) FIXES ISSUES WITH DEFINITION OF STRAIN
    MP.eta = minimum(MP.eta,P.S[p].eta_max) # COPYING FROM HERE: http://www.lmm.jussieu.fr/~lagree/TEXTES/PDF/JFMcollapsePYLLSSP11.pdf
    MP.dev_stress = MP.eta*MP.de_ij/P.dt

    MP.dp = P.S[p].K*MP.de_kk # tension positive
    MP.pressure += MP.dp
    if MP.pressure > 0.: MP.pressure = 0. # can't go into tension - this is really important!!

    # MP.stress = MP.pressure*eye(3) + MP.dev_stress
    MP.dstress = MP.pressure*eye(3) + MP.dev_stress - MP.stress

    for r in range(4):
        n = G.nearby_nodes(MP.n_star,r,P)
        G.pressure[n] += MP.N[r]*MP.pressure*MP.m
        G.dev_stress[n] += MP.N[r]*norm(MP.dev_stress)/sqrt(2.)*MP.m
        # G.dev_stress[n] += MP.N[r]*MP.dev_stress[0,1]*MP.m
        G.mu[n] += MP.N[r]*MP.mu*MP.m
        G.I[n] += MP.N[r]*MP.I*MP.m

def ken_simple(MP,P,G,p):
    """:math: constant `\mu` rheology with stress = 0 when rho < rho_c.

    This material model is UNVALIDATED and is DEFINITELY NOT WORKING.

    :param MP: A particle.Particle instance.
    :param P: A param.Param instance.
    :param G: A grid.Grid instance
    :param p: The particle number.
    :type p: int

    """

    s_bar = 0.
    for i in range(P.G.ns): s_bar += MP.phi[i]*P.G.s[i]

    MP.de_kk = trace(MP.dstrain)/3. # tension positive
    MP.de_ij = MP.dstrain - MP.de_kk*eye(3) # shear strain increment

    MP.eta = 2.*sqrt(2)*P.S[p].mu_0*abs(MP.pressure)/MP.gammadot # HACK: 2*SQRT(2) FIXES ISSUES WITH DEFINITION OF STRAIN
    MP.eta = minimum(nan_to_num(MP.eta),P.S[p].eta_max) # COPYING FROM HERE: http://www.lmm.jussieu.fr/~lagree/TEXTES/PDF/JFMcollapsePYLLSSP11.pdf


    MP.dev_stress = MP.eta*MP.de_ij/P.dt

    MP.dp = P.S[p].K*MP.de_kk # tension positive # FIXME do I need to multiply this by 3??
    MP.pressure += MP.dp
    if MP.pressure > 0: # can't go into tension - this is really important!!
        MP.dp -= MP.pressure # set increment back to zero
        MP.pressure = 0.

    if MP.rho > 2500:
        MP.dstress = MP.pressure*eye(3) + MP.dev_stress - MP.stress
    else:
        MP.dstress = -MP.stress

    if not isfinite(MP.dstress).all():
        print('THIS IS GOING TO BE A PROBLEM! FOUND SOMETHING NON-FINITE IN CONSTITUTIVE MODEL')
        print(MP.de_kk)
        print(MP.de_ij)
        print(s_bar)
        print(MP.eta)
        print(MP.dev_stress)
        print(MP.dp)
        print(MP.pressure)
        print(MP.dstress)
        sys.exit()

    for r in range(4):
        n = G.nearby_nodes(MP.n_star,r,P)
        G.pressure[n] += MP.N[r]*MP.pressure*MP.m
        G.dev_stress[n] += MP.N[r]*norm(MP.dev_stress)/sqrt(2.)*MP.m
        G.eta[n] += MP.N[r]*MP.eta*MP.m
