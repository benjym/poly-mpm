#!/usr/bin/python
# -*- coding: utf-8 -*-
# import colorcet as cc
import matplotlib.pyplot as plt
from matplotlib.cm import viridis, bwr
from matplotlib.colors import LinearSegmentedColormap, Normalize, LogNorm
import matplotlib.patches as patches
from numpy import *
from numpy.linalg import norm
import os
plt.viridis()

cdict = {'red': ((0.0, 1.0, 1.0),
                 (0.25, 1.0, 1.0),
                 (0.5, 1.0, 1.0),
                 (0.75, 0.902, 0.902),
                 (1.0, 0.0, 0.0)),
         'green': ((0.0, 0.708, 0.708),
                   (0.25, 0.302, 0.302),
                   (0.5, 0.2392, 0.2392),
                   (0.75, 0.1412, 0.1412),
                   (1.0, 0.0, 0.0)),
         'blue': ((0.0, 0.4, 0.4),
                  (0.25, 0.3569, 0.3569),
                  (0.5, 0.6078, 0.6078),
                  (0.75, 1., 1.),
                  (1.0, 1., 1.))}
orange_blue = LinearSegmentedColormap('orange_blue',cdict,256)
orange_blue.set_bad('w',1.0)
# print(orange_blue)
# viridis.set_under('w')
# viridis.set_over('w')
# bwr.set_under('w')
# bwr.set_over('w')
bwr.set_bad('k')

figs = {'backend': 'agg',
          'axes.labelsize': 10,
          'xtick.labelsize':10,
          'ytick.labelsize':10,
#           'text.usetex': True,
#           'font.serif' : 'serif',
#           'font.sans-serif' : 'cm',
#           'font.family' : 'sans-serif',
          }
plt.rcParams.update(figs)

class Plotting:
    """ Class containing all of the plotting functions.
    """

    def savefig(self,P,name,dpi=150):
        """ A method for saving figures in the right place """
#         root_dir = '~/Documents/poly-mpm/'
        root_dir = ''
        root_dir = os.path.expanduser(root_dir)
        if hasattr(P, 'supername'): P.save_dir = root_dir + P.supername + '/'
        else: P.save_dir = root_dir + 'im/' + P.mode + '/' + P.S[0].law + '/'
        if P.O.plot_material_points:
            for p in range(P.phases):
                save_dir_p = P.save_dir + 'MP_' + str(p) + '/'
                if not os.path.isdir(save_dir_p):
                    os.makedirs(save_dir_p)
        if P.O.plot_continuum:
            if not os.path.isdir(P.save_dir + 'Continuum'):
                os.makedirs(P.save_dir + 'Continuum/')
        if P.O.plot_gsd_mp or P.O.plot_gsd_grid:
            if not os.path.isdir(P.save_dir + 'GSD/'): os.makedirs(P.save_dir + 'GSD/')
            # if not os.path.isdir(P.save_dir + 'GSD/s_bar/'): os.makedirs(P.save_dir + 'GSD/s_bar/')
            # for i in range(P.G.ns):
                # if not os.path.exists(P.save_dir + 'GSD/phi_' + str(i) + '/'): os.makedirs(P.save_dir + 'GSD/phi_' + str(i))

        plt.subplots_adjust(hspace=0.4,wspace=0.4)
        plt.savefig(P.save_dir + name + '.png', dpi=dpi)
        plt.close()
#         print('Saved "' + name + '.png"                     ',end='\r')

    def draw_grid(self,G):
        """ Draw an overlay of the grid. Green for both directions, red crosses for horizontal boundaries and red plusses for vertical boundaries."""
        plt.plot(G.X*(1-G.boundary_tot),G.Y*(1-G.boundary_tot),'g+')
#         plt.plot(G.X*G.boundary_tot,G.Y*G.boundary_tot,'r+')
        plt.plot(G.X*G.boundary_h,G.Y*G.boundary_h,'rx')
        plt.plot(G.X*G.boundary_v,G.Y*G.boundary_v,'r+')

    def draw_pretty_grid(self,G):
        """ Draw an overlay of the grid that looks vaguely nice."""
        plt.plot(G.X*(1-G.boundary_tot),G.Y*(1-G.boundary_tot),'ko',alpha=0.5,markersize=1)
        plt.plot(G.X*G.boundary_h,G.Y*G.boundary_h,'ko',markersize=1)
        plt.plot(G.X*G.boundary_v,G.Y*G.boundary_v,'ko',markersize=1)

    def draw_gamma_dot(self,L,P,G):
        for p in range(P.phases):
            if P.S[p].law is not 'rigid':
                x = zeros((P.S[p].n,3))
                v = zeros((P.S[p].n,3))
                gammadot = zeros((P.S[p].n))
                for i in range(P.S[p].n):
                    x[i] = L.S[p][i].x
                    v[i] = L.S[p][i].v
                    gammadot[i] = L.S[p][i].gammadot/P.dt
                plt.clf()
                plt.figure(figsize=(12,10))
                self.draw_grid(G)
                plt.xlabel(r'$\dot\gamma$',rotation='horizontal')
                plt.scatter(x[:,0],x[:,1],c=gammadot,marker='s',edgecolor='None')
                plt.colorbar()
                self.savefig(P,str(P.save).zfill(5))
                P.save += 1

    def draw_continuum(self,G,P):
#         for p in range(P.phases):
        for p in [0,]: # all phases have same continuum properties!
            if P.S[p].law is not 'rigid':
                plt.clf()
                if hasattr(P.O, 'continuum_fig_size'):
                    figsize = P.O.continuum_fig_size
                else:
                    scale = 3
                    figsize=[2.*scale*(P.G.x_M-P.G.x_m),scale*(P.G.y_M-P.G.y_m)]
                plt.figure(figsize=figsize)

                if P.S[p].law == 'viscous' or P.S[p].law == 'viscous_size' or P.S[p].law == 'shear_maxwell' or P.S[p].law == 'marks2012':
                    plt.close()
                    fig, ax = plt.subplots(2,4,figsize=figsize)
                    ax = array(ax).flatten()

                    ax[0].plot((G.q[:,0]/G.m)[P.G.nx//2::P.G.nx],G.y,'k')
                    ax[0].plot(G.gammadot[P.G.nx//2::P.G.nx],G.y,'b')
                    ax[0].set_ylim(P.G.y_m,P.G.y_M)
                    ax[0].set_xlabel(r'$u$, $\dot\gamma$',rotation='horizontal')
                    ax[0].set_ylabel(r'$y$',rotation='horizontal')
                    ax[0].set_title(r'$|u_{max}| = ' + str(amax(abs(nan_to_num(G.q[:,0]/G.m)))) + '$',
                              rotation='horizontal')

                    titles = [r'$u$',r'$v$',r'$\bar s$',r'$\rho$',r'$|\dot\gamma|$',r'$\nabla(|\dot\gamma|)_{x}$',r'$\nabla(|\dot\gamma|)_{y}$']
                    z = [G.q[:,0],
                         G.q[:,1],
                         G.s_bar*G.m,
                         G.m*G.m/G.V,
                         abs(G.gammadot)*G.m,
                         G.grad_gammadot[:,0]*G.m,
                         G.grad_gammadot[:,1]*G.m,
                         ]

                    for i in range(len(titles)):
                        plt.sca(ax[i+1])
                        cax = plt.pcolormesh(G.x_plot,G.y_plot,
                                                   ma.masked_invalid(z[i]/G.m).reshape(P.G.ny,P.G.nx),
                                                   )
                        plt.title(titles[i],rotation='horizontal')
                        plt.xlim(P.G.x_m,P.G.x_M)
                        plt.ylim(P.G.y_m,P.G.y_M)
                        plt.colorbar()

                elif P.S[p].law == 'bingham' or P.S[p].law == 'pouliquen2D' or (P.S[p].law == 'pouliquen' or (P.S[p].law == 'HB' or P.S[p].law == 'linear_mu')):
                    plt.close()
                    fig, ax = plt.subplots(4,4,figsize=figsize)
                    ax = array(ax).flatten()

                    ax[0].plot((G.q[:,0]/G.m)[P.G.nx//2::P.G.nx],G.y,'k')
                    ax[1].plot(G.gammadot[P.G.nx//2::P.G.nx],G.y,'b')
                    ax[0].set_ylim(P.G.y_m,P.G.y_M)
                    ax[0].set_xlabel(r'$u$',rotation='horizontal')
                    ax[0].set_ylabel(r'$y$',rotation='horizontal')
                    if P.mode == 'bi_seg_test':
                        ax[0].set_title(r'$u_{max}^{pred} = ' + str(P.S[p].v_max) + '$\n$|u_{max}| = ' +
                                         str(amax(abs(nan_to_num(G.q[:,0]/G.m)))) + '$',
                                         rotation='horizontal')
                        ax[0].plot(sqrt(abs(P.max_g)*P.G.s_bar_0)*(2./3.)*P.S[p].I_0*
                                   (tan(abs(P.theta))-P.S[p].mu_0)/(P.S[p].mu_1-tan(abs(P.theta)))*
                                   sqrt(P.S[p].packing*cos(abs(P.theta)))*
                                   (P.G.y_M**1.5-(P.G.y_M-P.G.y)**1.5)/P.G.s_bar_0**1.5,
                                   G.y,
                                   'k--')
                    else:
                        ax[0].set_title(r'$|u_{max}| = ' + str(amax(abs(nan_to_num(G.q[:,0]/G.m)))) + '$',
                              rotation='horizontal')

                    # ax[1].plot(G.gammadot[P.G.nx//2::P.G.nx],G.y,'b')
                    # ax[1].set_ylim(P.G.y_m,P.G.y_M)
                    # ax[1].set_xlabel(r'$\dot\gamma$',rotation='horizontal')
                    # ax[1].set_ylabel(r'$y$',rotation='horizontal')
                    # ax[1].set_title(r'$|\dot\gamma_{max}| = ' + str(amax(abs(nan_to_num(G.gammadot)))) + '$',
                    #           rotation='horizontal')

                    titles = [r'$\rho$',r'$\bar s$',
                              r'$u$',r'$v$',#r'$\nabla(|\dot\gamma|)_{x}$',r'$\nabla(|\dot\gamma|)_{y}$',
                              r'$|\dot\gamma|$',r'$P$',r'$\sigma_{xy}$',r'$|\sigma_{xy}/P|$',
                              # r'$\mu$',r'$\log_{10}I$',r'$\dot\phi^{m}$',r'$\dot\phi^{M}$']
                              r'$\mu$',r'$\eta/\eta_{max}$',
                              r'$p_k$',r'$||\nabla p_k||$',r'$\dot\phi^{M}$',r'$\phi^{M}$']
                    norm = [Normalize(),Normalize(),
                            Normalize(),Normalize(),Normalize(),Normalize(),
                            Normalize(),Normalize(),Normalize(),LogNorm(),
                            Normalize(),Normalize(),Normalize(),Normalize()
                            ]
                    cmap = [viridis,orange_blue,#cc.cm.bmy_r,
                            viridis,viridis,#viridis,viridis,
                            viridis,viridis,viridis,viridis,
                            viridis,'bwr',
                            viridis,viridis,'bwr',viridis
                            ]
                    z = [G.m*G.m/G.V,
                         G.s_bar*G.m,
                         G.q[:,0],
                         G.q[:,1],
                         abs(G.gammadot)*G.m,
                         -G.pressure,
                         G.dev_stress,
                         abs(G.dev_stress/G.pressure)*G.m,
                         G.mu,
                         # log10(G.I/G.m)*G.m,
                         ma.masked_less_equal(G.eta/P.S[0].eta_max,0.), # just for log scaling nicely
                         # G.dphi[:,0]/P.dt*G.m,
                         G.pk*G.m,
                         sqrt(G.grad_pk[:,0]**2 + G.grad_pk[:,1]**2)*G.m,
                         # G.grad_pk[:,0]*G.m,
                         # G.grad_pk[:,1]*G.m,
                         G.dphi[:,-1]/P.dt*G.m,
                         G.phi[:,-1]*G.m
                         ]
                    for i in range(len(titles)):
                        plt.sca(ax[i+2])
                        cax = plt.pcolormesh(G.x_plot,
                                             G.y_plot,
                                             ma.masked_invalid(z[i]/G.m).reshape(P.G.ny,P.G.nx),
                                             norm=norm[i],
                                             cmap=cmap[i]
                                             )
                        plt.title(titles[i],rotation='horizontal')
                        plt.xlim(P.G.x_m,P.G.x_M)
                        plt.ylim(P.G.y_m,P.G.y_M)
                        if P.segregate_grid:
                            if titles[i] == r'$\bar s$':
                                plt.clim(P.G.s_m,P.G.s_M)
                                # plt.set_cmap('bwr')
                                # plt.set_cmap(cc.cm.bmy)
                            elif titles[i] == r'$\eta/\eta_{max}$':
                                # plt.clim(0,1)
                                plt.clim(1e-1,1e1)
                            elif titles[i] == r'$\dot\phi^{M}$':
                                plt.clim(-amax(abs(G.dphi[:,-1]))/P.dt,amax(abs(G.dphi[:,-1]))/P.dt)
                            elif titles[i] == r'$\phi^{M}$':
                                plt.clim(0,1)
                                # plt.set_cmap('bwr')
                            # else:
                                # plt.set_cmap('viridis')
                        plt.colorbar()

                elif P.S[p].law == 'dp' or P.S[p].law == 'dp_rate':
                    plt.close()
                    fig, ax = plt.subplots(2,4,figsize=figsize)
                    ax = array(ax).flatten()

                    titles = [r'$\rho$',r'$\bar s$',r'$u$',r'$v$',r'$p$',r'$q$',r'$|\dot\gamma|$',r'$q/p$']
                    z = [G.m*G.m/G.V,G.s_bar*G.m,G.q[:,0],G.q[:,1],G.pressure,G.dev_stress,abs(G.gammadot),G.m*G.dev_stress/G.pressure]
                    for i in range(8):
                        plt.sca(ax[i])
                        plt.pcolormesh(G.x_plot,G.y_plot,
                                       ma.masked_invalid(z[i]/G.m).reshape(P.G.ny,P.G.nx),
                                       )
                        plt.title(titles[i],rotation='horizontal')
                        plt.colorbar()
                else:
                    plt.subplot(331)
                    plt.xlabel(r'$\rho$',rotation='horizontal')
                    plt.contourf(G.x,G.y,(G.m/G.V).reshape(P.G.ny,P.G.nx))
                    plt.colorbar()
                    plt.subplot(363)
                    plt.xlabel(r'$u$',rotation='horizontal')
                    plt.contourf(G.x,G.y,
                        (G.q[:,0]/G.m).reshape(P.G.ny,P.G.nx))
                    plt.colorbar()
                    plt.subplot(364)
                    plt.xlabel(r'$v$',rotation='horizontal')
                    plt.contourf(G.x,G.y,
                        (G.q[:,1]/G.m).reshape(P.G.ny,P.G.nx))
                    plt.colorbar()
                    plt.subplot(333)
                    plt.xlabel(r'$y_\Phi$',rotation='horizontal')
                    plt.contourf(G.x,G.y,
                        (G.yieldfunction/G.m).reshape(P.G.ny,P.G.nx))
                    plt.colorbar()
                    plt.subplot(334)
                    plt.xlabel(r'$\sigma_v$',rotation='horizontal')
                    plt.contourf(G.x,G.y,
                        (G.sigmav/G.m).reshape(P.G.ny,P.G.nx))
                    plt.colorbar()
                    plt.subplot(335)
                    plt.xlabel(r'$\sigma_h$',rotation='horizontal')
                    plt.contourf(G.x,G.y,
                        (G.sigmah/G.m).reshape(P.G.ny,P.G.nx))
                    plt.colorbar()
                    plt.subplot(336)
                    plt.xlabel(r'$|\dot\varepsilon_{ij}|$',rotation='horizontal')
                    plt.contourf(G.x,G.y,
                        (G.gammadot/G.m).reshape(P.G.ny,P.G.nx))
                    plt.colorbar()
                    plt.subplot(337)
                    plt.xlabel(r'$P$',rotation='horizontal')
                    plt.contourf(G.x,G.y,
                        (-G.pressure/G.m).reshape(P.G.ny,P.G.nx))
                    plt.colorbar()
                    plt.subplot(338)
                    plt.xlabel(r'$|s_{ij}|$',rotation='horizontal')
                    plt.contourf(G.x,G.y,
                        (G.dev_stress/G.m).reshape(P.G.ny,P.G.nx))
                    plt.colorbar()
                    plt.subplot(339)
                    plt.xlabel(r'$|\dot s_{ij}|$',rotation='horizontal')
                    plt.contourf(G.x,G.y,
                        (G.dev_stress_dot/G.m).reshape(P.G.ny,P.G.nx))
                    plt.colorbar()
        #        for i in range(6):
        #            plt.subplot(321 + i)
        #            plt.xticks([])
        #            plt.yticks([])
                    #plt.plot(G.X*G.boundary_tot,G.Y*G.boundary_tot,'r+')
                self.savefig(P,'Continuum/'+str(P.grid_save).zfill(5))
        P.grid_save += 1

    def draw_material_points(self,L,P,G,name=False):
        plt.clf()
        if hasattr(P.O, 'mp_fig_size'):
            figsize = P.O.mp_fig_size
        else:
            scale = 3
            figsize=[2.*scale*(P.G.x_M-P.G.x_m),scale*(P.G.y_M-P.G.y_m)]
        plt.figure(figsize=figsize)

        for p in range(P.phases):
            if P.S[p].law == 'elastic':
                x = zeros((P.S[p].n,3))
                v = zeros((P.S[p].n,3))
                gammadot = zeros((P.S[p].n))
                pressure = zeros((P.S[p].n))
                sigmah = zeros((P.S[p].n))
                sigmav = zeros((P.S[p].n))
                for i in range(P.S[p].n):
                    x[i] = L.S[p][i].x
                    v[i] = L.S[p][i].v
                    gammadot[i] = L.S[p][i].gammadot/P.dt
                    sigmav[i] = L.S[p][i].sigmav
                    sigmah[i] = L.S[p][i].sigmah
                    pressure[i] = L.S[p][i].pressure
                size=25.

                plt.subplot(321)
                self.draw_grid(G)
                plt.xlabel(r"$MP's$")
                plt.scatter(x[:,0],x[:,1],s=size,marker='s',
                            edgecolor='None',facecolor='b')
                plt.subplot(322)
                plt.xlabel(r'${\bf u}$',rotation='horizontal')
                plt.quiver(x[:,0],x[:,1],v[:,0],v[:,1])#,scale=1.)
                plt.subplot(323)
                plt.xlabel(r'$\sigma_v$',rotation='horizontal')
                plt.scatter(x[:,0],x[:,1],s=size,c=sigmav,marker='s',edgecolor='None')
                plt.colorbar()
                plt.subplot(324)
                plt.xlabel(r'$\sigma_h$',rotation='horizontal')
                plt.scatter(x[:,0],x[:,1],s=size,c=sigmah,marker='s',edgecolor='None')
                plt.colorbar()
                plt.subplot(325)
                plt.xlabel(r'$P$',rotation='horizontal')
                plt.scatter(x[:,0],x[:,1],s=size,c=pressure,marker='s',edgecolor='None')
                plt.colorbar()
                plt.subplot(326)
                plt.xlabel(r'$\dot\gamma$',rotation='horizontal')
                plt.scatter(x[:,0],x[:,1],s=size,c=gammadot,marker='s',edgecolor='None')
                plt.colorbar()
            elif P.S[p].law == 'von_mises':
                x = zeros((P.S[p].n,3))
                v = zeros((P.S[p].n,3))
                gammadot = zeros((P.S[p].n))
                yieldfunction = zeros((P.S[p].n))
                pressure = zeros((P.S[p].n))
                dev_stress = zeros((P.S[p].n))
                dev_stress_dot = zeros((P.S[p].n))
                sigmah = zeros((P.S[p].n))
                sigmav = zeros((P.S[p].n))
                for i in range(P.S[p].n):
                    x[i] = L.S[p][i].x
                    v[i] = L.S[p][i].v
                    gammadot[i] = L.S[p][i].gammadot/P.dt
                    sigmav[i] = L.S[p][i].sigmav
                    sigmah[i] = L.S[p][i].sigmah
                    yieldfunction[i] = L.S[p][i].yieldfunction
                    pressure[i] = L.S[p][i].pressure
                    dev_stress[i] = norm(L.S[p][i].dev_stress)
                    dev_stress_dot[i] = norm(L.S[p][i].dev_dstress)/P.dt
                size=25.
                try:
                    scale = P.G.scale
                except AttributeError:
                    scale = 10
                plt.figure(figsize=[scale*(P.G.x_M-P.G.x_m),scale*(P.G.y_M-P.G.y_m)])

                plt.subplot(331)
                self.draw_grid(G)
                plt.xlabel(r"$MP's$")
                plt.scatter(x[:,0],x[:,1],s=size,marker='s',
                            edgecolor='None',facecolor='b')
                plt.subplot(332)
                plt.xlabel(r'${\bf u}$',rotation='horizontal')
                plt.quiver(x[:,0],x[:,1],v[:,0],v[:,1])#,scale=1.)
                plt.subplot(333)
                plt.xlabel(r'$y_\Phi$',rotation='horizontal')
                plt.scatter(x[:,0],x[:,1],s=size,c=yieldfunction,marker='s',edgecolor='None')
                plt.clim(-1,0)
                plt.colorbar()
                plt.subplot(334)
                plt.xlabel(r'$\sigma_v$',rotation='horizontal')
                plt.scatter(x[:,0],x[:,1],s=size,c=sigmav,marker='s',edgecolor='None')
                plt.colorbar()
                plt.subplot(335)
                plt.xlabel(r'$\sigma_h$',rotation='horizontal')
                plt.scatter(x[:,0],x[:,1],s=size,c=sigmah,marker='s',edgecolor='None')
                plt.colorbar()
                plt.subplot(336)
                plt.xlabel(r'$|\dot\varepsilon_{ij}|$',rotation='horizontal')
                plt.scatter(x[:,0],x[:,1],s=size,c=gammadot,marker='s',edgecolor='None')
                plt.colorbar()
                plt.subplot(337)
                plt.xlabel(r'$P$',rotation='horizontal')
                plt.scatter(x[:,0],x[:,1],s=size,c=pressure,marker='s',edgecolor='None')
                plt.colorbar()
                plt.subplot(338)
                plt.xlabel(r'$|\sigma_{ij}|$',rotation='horizontal')
                plt.scatter(x[:,0],x[:,1],s=size,c=dev_stress,marker='s',edgecolor='None')
                plt.colorbar()
                plt.subplot(339)
                plt.xlabel(r'$|\dot\sigma_{ij}|$',rotation='horizontal')
                plt.scatter(x[:,0],x[:,1],s=size,c=dev_stress_dot,
                            marker='s',edgecolor='None')
                plt.colorbar()
            elif P.S[p].law == 'dp':
                x = zeros((P.S[p].n,3))
                v = zeros((P.S[p].n,3))
                pressure = zeros((P.S[p].n))
                q = zeros((P.S[p].n))
                for i in range(P.S[p].n):
                    x[i] = L.S[p][i].x
                    v[i] = L.S[p][i].v
                    pressure[i] = L.S[p][i].pressure
                    q[i] = L.S[p][i].q
                size=25.
                try:
                    figsize = P.O.mp_fig_size
                except AttributeError:
                    figsize=[6,6]
#                 plt.add_subplot(figsize=figsize)
                plt.subplot(321)
                self.draw_grid(G)
                plt.xlabel(r"$MP's$")
                plt.scatter(x[:,0],x[:,1],s=size,marker='s',
                            edgecolor='None',facecolor='b')
                plt.subplot(322)
                plt.xlabel(r'${\bf u}$',rotation='horizontal')
                plt.quiver(x[:,0],x[:,1],v[:,0],v[:,1])#,scale=1.)
                plt.subplot(323)
                plt.xlabel(r'$p$',rotation='horizontal')
                plt.scatter(x[:,0],x[:,1],s=size,c=pressure,marker='s',edgecolor='None')
                plt.colorbar()
                plt.subplot(324)
                plt.xlabel(r'$q$',rotation='horizontal')
                plt.scatter(x[:,0],x[:,1],s=size,c=q,marker='s',edgecolor='None')
                plt.colorbar()
            elif P.S[p].law == 'rigid':
                # pass
                x = zeros((P.S[p].n,3))
                v = zeros((P.S[p].n,3))
                for i in range(P.S[p].n):
                    x[i] = L.S[p][i].x
                    v[i] = L.S[p][i].v
                size=15.

                plt.subplot(211)
                self.draw_grid(G)
                plt.xlabel(r"$MP's$")
                plt.scatter(x[:,0],x[:,1],s=size,marker='.',edgecolor='None',facecolor='b')
                plt.xlim(P.G.x_m-G.dx/2,P.G.x_M+G.dx/2)
                plt.ylim(P.G.y_m-G.dy/2,P.G.y_M+G.dy/2)
                plt.subplot(212)
                plt.xlabel(r'${\bf u}$',rotation='horizontal')
                plt.quiver(x[:,0],x[:,1],v[:,0],v[:,1])#,scale=1.)
            else:
                x = zeros((P.S[p].n,3))
                v = zeros((P.S[p].n,3))
                for i in range(P.S[p].n):
                    x[i] = L.S[p][i].x
                    v[i] = L.S[p][i].v
                size=15.

                plt.subplot(211)
                self.draw_grid(G)
                plt.xlabel(r"$MP's$")
                plt.scatter(x[:,0],x[:,1],s=size,marker='.',edgecolor='None',facecolor='b')
                plt.xlim(P.G.x_m-G.dx/2,P.G.x_M+G.dx/2)
                plt.ylim(P.G.y_m-G.dy/2,P.G.y_M+G.dy/2)
                plt.subplot(212)
                plt.xlabel(r'${\bf u}$',rotation='horizontal')
                plt.quiver(x[:,0],x[:,1],v[:,0],v[:,1])#,scale=1.)
            # if P.S[p].law is not 'rigid':
            if name:
                self.savefig(P,'MP_'+str(p)+'/'+str(name).zfill(5))
            else:
                self.savefig(P,'MP_'+str(p)+'/'+str(P.mp_save).zfill(5))
        P.mp_save += 1

    def draw_energy(self,P):
        plt.clf()
        t = P.O.energy[:,0]
        plt.plot(t,P.O.energy[:,1],'b.',label='KE') # kinetic energy
        plt.plot(t,P.O.energy[:,2],'r.',label='GPE') # gravitational potential energy
        plt.plot(t,P.O.energy[:,3],'g.',label='strain energy') # strain energy
        plt.plot(t,sum(P.O.energy[:,1:],1),'k.',label='total energy')
        plt.legend(loc=0)
        plt.xlabel('Time')
        plt.ylabel('Energy')
        self.savefig(P,'ns_' + str(P.S[0].n))
        plt.close()

    def draw_gsd_mp(self,L,P,G):
        w = 0.1/P.G.nx
        for p in range(P.phases):
            plt.clf()
            aspect_ratio = P.G.ny/P.G.nx
            fig = plt.figure(figsize=[6,6*aspect_ratio])
            ax = plt.subplot(111)
            self.draw_pretty_grid(G)
            for i in range(P.S[p].n):
                phi_prev = 0.
                for s in range(P.G.ns):
                # for s in [1]:
                    # rectangle
                    # rect = patches.Rectangle(L.S[p][i].x + array([w*(phi_prev-0.5),-0.5*w,0.]),
                    #                          w*L.S[p][i].phi[s],w,
                    #                          edgecolor='none',
                    #                          facecolor=cc.cm.bmy_r(P.G.s[s]/P.G.s[-1])
                    #                          )
                    # circular wedge
                    rect = patches.Wedge(L.S[p][i].x[:2], # center
                                         w, # radius
                                         phi_prev*360., # theta1
                                         (phi_prev + L.S[p][i].phi[s])*360., # theta2
                                         edgecolor='none',
                                         # facecolor=cc.cm.bmy_r(P.G.s[s]/P.G.s[-1])
                                         facecolor=orange_blue(P.G.s[s]/P.G.s[-1])

                                         )
                    ax.add_patch(rect)
                    phi_prev += L.S[p][i].phi[s]
            plt.axis('equal')
            ax.set_axis_off()
            plt.subplots_adjust(left=0.,right=1.,bottom=0.,top=1.)
            self.savefig(P,'GSD/MP_'+str(p)+'_'+str(P.mp_save).zfill(5),dpi=500)
        P.mp_save += 1

    def save_s_bar(self,L,P,G):
        if not os.path.isdir(P.save_dir + 'data/'): os.makedirs(P.save_dir + 'data/')
        save(P.save_dir + 'data/s_bar_' + str(P.grid_save).zfill(5)+'.npy',G.s_bar.reshape(P.G.ny,P.G.nx))

    def save_density(self,L,P,G):
        if not os.path.isdir(P.save_dir + 'data/'): os.makedirs(P.save_dir + 'data/')
        save(P.save_dir + 'data/density_' + str(P.grid_save).zfill(5)+'.npy',(G.m/G.V).reshape(P.G.ny,P.G.nx))

    def save_u(self,L,P,G):
        if not os.path.isdir(P.save_dir + 'data/'): os.makedirs(P.save_dir + 'data/')
        save(P.save_dir + 'data/u_' + str(P.grid_save).zfill(5)+'.npy',(G.q[:,0]/G.m).reshape(P.G.ny,P.G.nx))
        save(P.save_dir + 'data/v_' + str(P.grid_save).zfill(5)+'.npy',(G.q[:,1]/G.m).reshape(P.G.ny,P.G.nx))

    def save_phi_MP(self,L,P,G):
        if not os.path.isdir(P.save_dir + 'data/'): os.makedirs(P.save_dir + 'data/')
        phi = []
        for p in range(P.phases):
            for i in range(P.S[p].n):
                phi.append(hstack([L.S[p][i].x[0],L.S[p][i].x[1],L.S[p][i].phi]))
        save(P.save_dir + 'data/MP_phi_' + str(P.mp_save).zfill(5)+'.npy',array(phi))

    # def save_field(self,L,P,G,field,fieldname):
    #     if not os.path.isdir(P.save_dir + 'data/'): os.makedirs(P.save_dir + 'data/')
    #     save(P.save_dir + 'data/' + fieldname + '_' + str(P.grid_save).zfill(5)+'.npy',field.reshape(P.G.ny,P.G.nx))

    def draw_gsd_grid(self,L,P,G):
        plt.clf()
        plt.pcolormesh(G.x_plot,
                       G.y_plot,
                       ma.masked_where(G.m<P.M_tol,G.s_bar).reshape(P.G.ny,P.G.nx),
                       vmin=P.G.s[0],
                       vmax=P.G.s[-1],
                       cmap=bwr,
                       )
        plt.colorbar()
        plt.xlim(P.G.x_m,P.G.x_M)
        plt.ylim(P.G.y_m,P.G.y_M)
        plt.title(r'$\bar s=$' + str(mean(G.s_bar)))
        self.savefig(P,'GSD/s_bar/'+str(P.grid_save).zfill(5))

        for i in range(P.G.ns):
            plt.clf()
            plt.pcolormesh(G.x_plot,
                           G.y_plot,
                           ma.masked_where(G.m<P.M_tol,G.phi[:,i]).reshape(P.G.ny,P.G.nx),
                           vmin=0.0,
                           vmax=1.0,
                           )
            plt.colorbar()
            plt.title(mean(G.phi[:,i]))
            self.savefig(P,'GSD/phi_' + str(i) + '/' +str(P.grid_save).zfill(5))
#             if P.t > 0:
#                 plt.clf()
#                 U = sqrt(G.u_hat[:,i]**2 + G.v_hat[:,i]**2).reshape(P.G.ny,P.G.nx)
#
#                 plt.streamplot(G.x,
#                                G.y,
#                                G.u_hat[:,i].reshape(P.G.ny,P.G.nx),
#                                G.v_hat[:,i].reshape(P.G.ny,P.G.nx),
#                                linewidth = U
#                                )
#                 self.savefig(P,'GSD/phi_' + str(i) + '/u_' +str(P.grid_save).zfill(5))

#         plt.clf()
#         plt.pcolormesh(G.x,
#                        G.y,
#                        sum(G.phi,1).reshape(P.G.ny,P.G.nx),
# #                            vmin=0.0,
# #                            vmax=1.0,
#                        )
#         plt.colorbar()
#         self.savefig(P,'GSD/phi_sum_' +str(P.grid_save).zfill(5))

    def draw_voronoi(self,P,G):
#        plt.clf()
#        plt.plot(G.xpt,G.ypt,'b.')
#        plt.plot(G.x_list,G.y_list,'k-')
#        plt.fill(G.x_list,G.y_list,'g',alpha=0.25,edgecolor='none')
#        plt.axis([P.G.x_m-P.G.dx,P.G.x_M+P.G.dx,P.G.y_m-P.G.dy,P.G.y_M+P.G.dy])
#        plt.axis([P.G.x_m,P.G.x_M,P.G.y_m,P.G.y_M])
#        self.savefig(P,str(30000+P.save).zfill(5))
        import subprocess
        subprocess.Popen('povray +W800 +H600 +A0.01 -GA +O' + str(30000+P.save).zfill(5) + '.png cpp/import.pov',shell=True)
        if not P.O.plot_material_points:
            P.save += 1
