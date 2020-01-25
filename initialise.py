import os
import numpy as np
from grid import Grid
from particles import Particles


def get_parameters(params):
    input_file = 'from inputs.' + params[0] +' import Params'
    exec(input_file, globals())
    P = Params(params) # Initialise parameters from input file

    if not hasattr(P, 't'): P.t = 0. # physical time (s)
    if not hasattr(P, 'tstep'): P.tstep = 0
    if not hasattr(P, 'grid_save'): P.grid_save = 0 # save counter
    if not hasattr(P, 'mp_save'): P.mp_save = 0 # save counter
    if not hasattr(P, 'M_tol'): P.M_tol = 1e-10 # very small mass (kg)
    if not hasattr(P, 'nt'): P.nt = int(P.t_f/P.dt) # number of timesteps
    if not hasattr(P, 'max_g'): P.max_g = 0.
    if not hasattr(P, 'max_q'): P.max_q = 0.
    if not hasattr(P, 'g'): P.g = P.max_g#array([0,0,0])
    if not hasattr(P, 'q_v'): P.q_v = P.max_q#array([0,0,0])
    if not hasattr(P, 'q_h'): P.q_h = P.max_q#array([0,0,0])
    if not hasattr(P, 'has_yielded'): P.has_yielded = False
    if not hasattr(P, 'damping'): P.damping = False
    if not hasattr(P, 'theta'): P.theta = 0. # vertical
    if not hasattr(P, 'pressure'): P.pressure = None # lithostatic
    if not hasattr(P, 'segregate'): P.segregate = False
    if not hasattr(P, 'initial_flow'): P.initial_flow = False # give initial velocities which correspond to 1D steady state
    if not hasattr(P, 'smooth_gamma_dot'): P.smooth_gamma_dot = False # use artifical smoothing on calculation of the shear strain rate
    if not hasattr(P, 'smooth_grad2'): P.smooth_grad2 = False # use artifical smoothing on calculation of the gradient of the shear strain rate
    if not hasattr(P, 'normalise_phi'): P.normalise_phi = False # normalise phi every time step HACK: THIS IS CHEATING!!! or is it?
    if not hasattr(P, 'time_stepping'): P.time_stepping = 'static' # do not update time stepping automatically - NOTE: ASSUMES A FUNCTION update_timestep() EXISTS FOR P.S
    if not hasattr(P, 'CFL'): P.CFL = 0.1 # default value of CFL condition is quite conservative

    if not hasattr(P.B, 'wall'): P.B.wall = False
    if not hasattr(P.B, 'two_walls'): P.B.two_walls = False
    if not hasattr(P.B, 'roughness'): P.B.roughness = False
    if not hasattr(P.B, 'wall_mu'): P.B.wall_mu = False
    if not hasattr(P.B, 'no_slip_bottom'): P.B.no_slip_bottom = False
    if not hasattr(P.B, 'cyclic_lr'): P.B.cyclic_lr = False
    if not hasattr(P.B, 'has_top'): P.B.has_top = False
    if not hasattr(P.B, 'has_bottom'): P.B.has_bottom = False
    if not hasattr(P.B, 'has_right'): P.B.has_right = False
    if not hasattr(P.B, 'has_left'): P.B.has_left = False
    if not hasattr(P.B, 'box'): P.B.box = False
    if not hasattr(P.B, 'outlet_left'): P.B.outlet_left = False
    if not hasattr(P.B, 'outlet_bottom'): P.B.outlet_bottom = False
    if not hasattr(P.B, 'inlet_right'): P.B.inlet_right = False
    if not hasattr(P.B, 'inlet_top'): P.B.inlet_top = False
    if not hasattr(P.B, 'force_boundaries'): P.B.force_boundaries = False
    if not hasattr(P.B, 'vertical_force'): P.B.vertical_force = False
    if not hasattr(P.B, 'horizontal_force'): P.B.horizontal_force = False
    if not hasattr(P.B, 'conveyor'): P.B.conveyor = False
    if not hasattr(P.B, 'silo_left'): P.B.silo_left = False
    if not hasattr(P.B, 'silo_bottom'): P.B.silo_bottom = False

    if not hasattr(P.G, 'thickness'): P.G.thickness = 1. # (m) into page
    if not hasattr(P.G, 'ns'): P.G.ns = 1 # number of grain sizes
    if not hasattr(P.G, 's'): P.G.s = [0.001] # 1mm grain size by default (m)
    if not hasattr(P.G, 's_m'): P.G.s_m = np.min(P.G.s)
    if not hasattr(P.G, 's_M'): P.G.s_M = np.max(P.G.s)

#         root_dir = '~/Documents/poly-mpm/'
    root_dir = ''
    # root_dir = os.path.expanduser(root_dir)
    if hasattr(P, 'supername'): P.save_dir = root_dir + P.supername + '/'
    else: P.save_dir = root_dir + 'im/' + params[0] + '/' + P.S.law + '/'

    if not hasattr(P.O, 'check_positions'): P.O.check_positions = False
    if not hasattr(P.O, 'measure_stiffness'): P.O.measure_stiffness = False
    if not hasattr(P.O, 'measure_energy'): P.O.measure_energy = False
    if not hasattr(P.O, 'plot_gsd_mp'): P.O.plot_gsd_mp = False
    if not hasattr(P.O, 'plot_gsd_grid'): P.O.plot_gsd_grid = False
    if not hasattr(P.O, 'plot_gsd_debug'): P.O.plot_gsd_debug = False
    if not hasattr(P.O, 'plot_material_points'): P.O.plot_material_points = False
    if not hasattr(P.O, 'plot_continuum'): P.O.plot_continuum = False
    if not hasattr(P.O, 'save_s_bar'): P.O.save_s_bar = False
    if not hasattr(P.O, 'save_density'): P.O.save_density = False
    if not hasattr(P.O, 'save_u'): P.O.save_u = False
    if not hasattr(P.O, 'save_phi_MP'): P.O.save_phi_MP = False
    if not hasattr(P.O, 'at_start'): P.O.at_start = null_func
    if not hasattr(P.O, 'at_end'): P.O.at_end = null_func
    if not hasattr(P.O, 'after_ever_timestep'): P.O.after_ever_timestep = null_func
    if not hasattr(P.O, 'after_every_savetime'): P.O.after_every_savetime = null_func

    P.mode = params[0]
    P.update_forces()
    G = Grid(P) # Initialise grid
    L = Particles(P,G)

    if not hasattr(P, 'update_timestep'):
        def update_timestep(P,G):
            distance = np.minimum(P.G.dx,P.G.dy)
            t_c = [0.1*distance/np.amax(np.abs([np.nan_to_num(G.q[:,0]/G.m),np.nan_to_num(G.q[:,1]/G.m)]))] # cell crossing condition
            if P.segregate:
                t_seg = distance/(P.c*(P.G.s_M/P.G.s_m - 1.)*np.amax(np.abs(np.nan_to_num(G.grad_pk)))) # NOTE: CHECK THIS
                t_c.append(t_seg)
            if hasattr(P.S, 'critical_time'): t_c.append(P.S.critical_time(P))
            P.dt = P.CFL*np.min(t_c)
            # if t_seg == min(t_c): print('\nTimestep limited by segregation', end='\n')
        P.update_timestep = update_timestep

    # Save input params
    # with open(P.save_dir + 'ini.json', 'w') as f:
    #     out = jsonpickle.encode(P)
    #     print(out)
    # print(P.G.__dict__)
    # print(P.B.__dict__)
    # print(P.O.__dict__)
    # print(P.S.__dict__)

    return P,G,L

def null_func(*argv):
    pass

if __name__ == '__main__':
    P,G,L = get_parameters(params)
