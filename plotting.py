import matplotlib.pyplot as plt
from numpy import cumsum, zeros, sign, linspace, sum, abs, maximum, amax, nan_to_num
from numpy.linalg import norm
import os

figs = {'backend': 'agg',
          'axes.labelsize': 10,
          'xtick.labelsize': 6,
          'ytick.labelsize': 6,
          'text.usetex': True,
          'text.usefont': 'serif'}
plt.rcParams.update(figs)

class Plotting:
    def savefig(self,P,name):
        root_dir = '~/Documents/Python/mpm/'
        root_dir = os.path.expanduser(root_dir)
        try:
            P.supername
            save_dir = root_dir + P.supername + '/'
        except:
            save_dir = root_dir + 'im/' + P.mode + '/' + P.S[0].law + '/'
        for p in range(P.phases):
            save_dir_p = save_dir + str(p) + '/'
            if not os.path.isdir(save_dir_p):
                os.makedirs(save_dir_p)
        plt.savefig(save_dir + name + '.png', dpi=100)
        plt.close()
        print 'Saved "' + name + '.png"'
        
    def draw_grid(self,G):
        plt.plot(G.X*(1-G.boundary_tot),G.Y*(1-G.boundary_tot),'g+')
        plt.plot(G.X*G.boundary_tot,G.Y*G.boundary_tot,'r+')
    
    def draw_gamma_dot(self,L,P,G):
        for p in range(P.phases):
            if P.S[p].law is not 'rigid':
                x = zeros((P.S[p].n,3))
                v = zeros((P.S[p].n,3))
                gammadot = zeros((P.S[p].n))
                for i in xrange(P.S[p].n):
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
        for p in range(P.phases):
            if P.S[p].law is not 'rigid':
                plt.clf()
                try:
                    figsize = P.O.continuum_fig_size
                except AttributeError:
                    scale = 6
                    figsize=[2.*scale*(P.G.x_M-P.G.x_m),scale*(P.G.y_M-P.G.y_m)]
                plt.figure(figsize=figsize) 
                
                if P.S[p].law == 'viscous':
                    plt.subplots_adjust(hspace=0.4,wspace=0.4)
                    plt.subplot(231)
                    plt.plot(P.G.y,(G.q[:,0]/G.m).reshape(P.G.ny,P.G.nx))[P.G.nx/2::P.G.nx]
                    plt.ylabel(r'$u$',rotation='horizontal')
                    plt.xlabel(r'$y$',rotation='horizontal')
                    plt.title(r'$|u_{max}| = ' + str(amax(abs(nan_to_num(G.q[:,0]/G.m)))) + '$',
                              rotation='horizontal')
                    plt.subplot(232)
                    plt.title(r'$u$',rotation='horizontal')
                    plt.contourf(G.x,G.y,
                        (G.q[:,0]/G.m).reshape(P.G.ny,P.G.nx))
                    plt.colorbar()
                    plt.subplot(233)
                    plt.title(r'$v$',rotation='horizontal')
                    plt.contourf(G.x,G.y,
                        (G.q[:,1]/G.m).reshape(P.G.ny,P.G.nx))
                    plt.colorbar()
                    plt.subplot(245)
                    plt.title(r'$\rho$',rotation='horizontal')
                    plt.contourf(G.x,G.y,(G.m/G.V).reshape(P.G.ny,P.G.nx))#,
        #                          levels=[P.S.rho-1,P.S.rho-0.5,P.S.rho+0.5,P.S.rho+1])
                    plt.colorbar()
                    plt.subplot(246)
                    plt.title(r'$P$',rotation='horizontal')
                    plt.contourf(G.x,G.y,
                        (G.pressure/G.m).reshape(P.G.ny,P.G.nx))
                    plt.colorbar()
                    plt.subplot(247)
                    plt.title(r'$|\dot\gamma|$',rotation='horizontal')
                    plt.contourf(G.x,G.y,
                        (G.gammadot/G.m).reshape(P.G.ny,P.G.nx))
                    plt.colorbar()
                    plt.subplot(248)
                    plt.title(r'$||s_{ij}||$',rotation='horizontal')
                    plt.contourf(G.x,G.y,
                        (G.dev_stress/G.m).reshape(P.G.ny,P.G.nx))
                    plt.colorbar()
                elif P.S[p].law == 'bingham' or (P.S[p].law == 'pouliquen' or P.S[p].law == 'HB'):
                    plt.subplots_adjust(hspace=0.4,wspace=0.4)
                    plt.subplot(331)
                    plt.plot(P.G.y,(G.q[:,0]/G.m).reshape(P.G.ny,P.G.nx))[P.G.nx/2::P.G.nx]
                    plt.ylabel(r'$u$',rotation='horizontal')
                    plt.xlabel(r'$y$',rotation='horizontal')
                    plt.title(r'$|u_{max}| = ' + str(amax(abs(G.q[:,0]/G.m))) + '$',
                              rotation='horizontal')
                    plt.subplot(332)
                    plt.title(r'$u$',rotation='horizontal')
                    plt.contourf(G.x,G.y,
                        (G.q[:,0]/G.m).reshape(P.G.ny,P.G.nx))
                    plt.colorbar()
                    plt.subplot(333)
                    plt.title(r'$v$',rotation='horizontal')
                    plt.contourf(G.x,G.y,
                        (G.q[:,1]/G.m).reshape(P.G.ny,P.G.nx))
                    plt.colorbar()
                    plt.subplot(334)
                    plt.title(r'$\rho$',rotation='horizontal')
                    plt.contourf(G.x,G.y,(G.m/G.V).reshape(P.G.ny,P.G.nx))#,
        #                          levels=[P.S.rho-1,P.S.rho-0.5,P.S.rho+0.5,P.S.rho+1])
                    plt.colorbar()
                    plt.subplot(335)
                    plt.title(r'$P$',rotation='horizontal')
                    plt.contourf(G.x,G.y,
                        (G.pressure/G.m).reshape(P.G.ny,P.G.nx))
                    plt.colorbar()
                    plt.subplot(336)
                    plt.title(r'$|\dot\gamma|$',rotation='horizontal')
                    plt.contourf(G.x,G.y,
                        (G.gammadot/G.m).reshape(P.G.ny,P.G.nx))
                    plt.colorbar()
                    plt.subplot(337)
                    plt.title(r'$||s_{ij}||$',rotation='horizontal')
                    plt.contourf(G.x,G.y,
                        (G.dev_stress/G.m).reshape(P.G.ny,P.G.nx))
                    plt.colorbar()
                    plt.subplot(338)
                    plt.title(r'$y$',rotation='horizontal')
                    plt.contourf(G.x,G.y,
                        (G.yieldfunction/G.m).reshape(P.G.ny,P.G.nx))
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
                        (G.pressure/G.m).reshape(P.G.ny,P.G.nx))
                    plt.colorbar()
                    plt.subplot(338)
                    plt.xlabel(r'$|\sigma_{ij}|$',rotation='horizontal')
                    plt.contourf(G.x,G.y,
                        (G.dev_stress/G.m).reshape(P.G.ny,P.G.nx))
                    plt.colorbar()
                    plt.subplot(339)
                    plt.xlabel(r'$|\dot\sigma_{ij}|$',rotation='horizontal')
                    plt.contourf(G.x,G.y,
                        (G.dev_stress_dot/G.m).reshape(P.G.ny,P.G.nx))
                    plt.colorbar()
        #        for i in xrange(6):
        #            plt.subplot(321 + i)
        #            plt.xticks([])
        #            plt.yticks([])
                    #plt.plot(G.X*G.boundary_tot,G.Y*G.boundary_tot,'r+')
                self.savefig(P,str(p)+'/'+str(10000+P.grid_save).zfill(5))
        P.grid_save += 1

    def draw_material_points(self,L,P,G,name=False):
        plt.clf()
        for p in range(P.phases):
            if P.S[p].law == 'elastic':
                x = zeros((P.S[p].n,3))
                v = zeros((P.S[p].n,3))
                gammadot = zeros((P.S[p].n))
                pressure = zeros((P.S[p].n))
                sigmah = zeros((P.S[p].n))
                sigmav = zeros((P.S[p].n))
                for i in xrange(P.S[p].n):
                    x[i] = L.S[p][i].x
                    v[i] = L.S[p][i].v
                    gammadot[i] = L.S[p][i].gammadot/P.dt
                    sigmav[i] = L.S[p][i].sigmav
                    sigmah[i] = L.S[p][i].sigmah
                    pressure[i] = L.S[p][i].pressure
                size=25.
                try:
                    scale = P.G.scale
                except AttributeError:
                    scale = 3
                plt.figure(figsize=[scale*(P.G.x_M-P.G.x_m),scale*(P.G.y_M-P.G.y_m)])    
                plt.subplot(321)
                self.draw_grid(G)
                plt.xlabel(r"$MP's$")
                plt.scatter(x[:,0],x[:,1],s=size,marker='s',
                            edgecolor='None',facecolor='b')
                plt.subplot(322)
                plt.xlabel(r'{\bf $u$}',rotation='horizontal')
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
                for i in xrange(P.S[p].n):
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
                plt.xlabel(r'{\bf $u$}',rotation='horizontal')
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
                p = zeros((P.S[p].n))
                q = zeros((P.S[p].n))
                y = zeros((P.S[p].n))
                for i in xrange(P.S[p].n):
                    x[i] = L.S[p][i].x
                    v[i] = L.S[p][i].v
                    p[i] = L.S[p][i].p
                    q[i] = L.S[p][i].q
                    y[i] = L.S[p][i].y
                size=25.
                try:
                    figsize = P.O.mp_fig_size
                except AttributeError:
                    figsize=[6,6]
                plt.add_subplot(figsize=figsize)
                plt.subplot(321)
                self.draw_grid(G)
                plt.xlabel(r"$MP's$")
                plt.scatter(x[:,0],x[:,1],s=size,marker='s',
                            edgecolor='None',facecolor='b')
                plt.subplot(322)
                plt.xlabel(r'{\bf $u$}',rotation='horizontal')
                plt.quiver(x[:,0],x[:,1],v[:,0],v[:,1])#,scale=1.)
                plt.subplot(323)
                plt.xlabel(r'$p$',rotation='horizontal')
                plt.scatter(x[:,0],x[:,1],s=size,c=p,marker='s',edgecolor='None')
                plt.colorbar()
                plt.subplot(324)
                plt.xlabel(r'$q$',rotation='horizontal')
                plt.scatter(x[:,0],x[:,1],s=size,c=q,marker='s',edgecolor='None')
                plt.colorbar()
                plt.subplot(325)
                plt.xlabel(r'$y$',rotation='horizontal')
                plt.scatter(x[:,0],x[:,1],s=size,c=y,marker='s',edgecolor='None')
                plt.colorbar()
            elif P.S[p].law == 'rigid':
                pass
            else:
                x = zeros((P.S[p].n,3))
                v = zeros((P.S[p].n,3))
                for i in xrange(P.S[p].n):
                    x[i] = L.S[p][i].x
                    v[i] = L.S[p][i].v
                size=15.
                try:
                    figsize = P.O.mp_fig_size
                except AttributeError:
                    figsize=[6,6]
                plt.figure(figsize=figsize)    
                plt.subplot(121)
                self.draw_grid(G)
                plt.xlabel(r"$MP's$")
                plt.scatter(x[:,0],x[:,1],s=size,marker='s',edgecolor='None',facecolor='b')
                plt.subplot(122)
                plt.xlabel(r'{\bf $u$}',rotation='horizontal')
                plt.quiver(x[:,0],x[:,1],v[:,0],v[:,1])#,scale=1.)
            if P.S[p].law is not 'rigid':
                if name:
                    self.savefig(P,str(p)+'/'+str(name).zfill(5))
                else:
                    self.savefig(P,str(p)+'/'+str(P.mp_save).zfill(5))
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
        
    def draw_fluid(self,L,P,G):
        x = zeros((P.F.n,3))
        v = zeros((P.F.n,3))
        for i in xrange(P.F.n):
            x[i] = L.F[i].x
            v[i] = L.F[i].v
        size=25.
        plt.clf()
        plt.figure()
        plt.subplot(121)
        self.draw_grid(G)
        plt.xlabel(r"$MP's$")
        plt.scatter(x[:,0],x[:,1],s=size,marker='s',
                    edgecolor='None',facecolor='b')
        plt.subplot(122)
        plt.xlabel(r'{\bf $u$}',rotation='horizontal')
        plt.quiver(x[:,0],x[:,1],v[:,0],v[:,1])#,scale=1.)
        self.savefig(P,str(P.save).zfill(5))
        P.save += 1
        
    def draw_gsd(self,L,P,G):
        x = zeros((P.S[0].n,3))
        v = zeros((P.S[0].n,3))
        s = zeros((P.S[0].n))
        for i in xrange(P.S[0].n):
            x[i] = L.S[0][i].x
            v[i] = L.S[0][i].v
            s[i] = L.S[0][i].s_bar
        size=25.
        plt.clf()
        plt.figure()
        self.draw_grid(G)
        plt.xlabel(r"$MP's$")
        plt.scatter(x[:,0],x[:,1],s=size,c=s,marker='s',
                    edgecolor='None')
        plt.colorbar()
        self.savefig(P,str(20000+P.save).zfill(5))
    
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
