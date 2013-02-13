import matplotlib.pyplot as plt
from numpy import cumsum, zeros, sign, linspace, sum, abs
from numpy.linalg import norm
import os

class Plotting:
    def savefig(self,P,name):
        try:
            P.G.scale
            save_dir = './im/' + P.mode + '/' + P.S.law + '/' + str(P.G.scale) + '/'
        except AttributeError:
            save_dir = './im/' + P.mode + '/' + P.S.law + '/'
        if not os.path.isdir(save_dir):
            os.makedirs(save_dir)
        plt.savefig(save_dir + name + '.png', dpi=100)
        plt.close()
        print 'Saved "' + name + '.png"'
    def draw_grid(self,G):
        plt.plot(G.X*(1-G.boundary),G.Y*(1-G.boundary),'g+')
        plt.plot(G.X*G.boundary,G.Y*G.boundary,'r+')
        
    def draw_gamma_dot(self,L,P,G):
        x = zeros((P.S.n,3))
        v = zeros((P.S.n,3))
        gammadot = zeros((P.S.n))
        for i in xrange(P.S.n):
            x[i] = L.S[i].x
            v[i] = L.S[i].v
            gammadot[i] = L.S[i].gammadot/P.dt
        plt.clf()
        plt.figure(figsize=(12,10))
        self.draw_grid(G)
        plt.xlabel(r'$\dot\gamma$',rotation='horizontal')
        plt.scatter(x[:,0],x[:,1],c=gammadot,marker='s',edgecolor='None')
        plt.colorbar()
        self.savefig(P,str(P.save).zfill(5))
        P.save += 1
        
    def draw_continuum(self,G,P):
        plt.clf()
#        try:
#            scale = P.G.scale
#        except AttributeError:
#            scale = 5
        scale = 2
        plt.figure(figsize=[scale*(P.G.x_M-P.G.x_m),scale*(P.G.y_M-P.G.y_m)]) 
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
            #plt.plot(G.X*G.boundary,G.Y*G.boundary,'r+')
        self.savefig(P,str(10000+P.save).zfill(5))
#        P.save += 1

    def draw_material_points(self,L,P,G,name=False):
        if P.S.law == 'elastic':
            x = zeros((P.S.n,3))
            v = zeros((P.S.n,3))
            gammadot = zeros((P.S.n))
            pressure = zeros((P.S.n))
            sigmah = zeros((P.S.n))
            sigmav = zeros((P.S.n))
            for i in xrange(P.S.n):
                x[i] = L.S[i].x
                v[i] = L.S[i].v
                gammadot[i] = L.S[i].gammadot/P.dt
                sigmav[i] = L.S[i].sigmav
                sigmah[i] = L.S[i].sigmah
                pressure[i] = L.S[i].pressure
            size=25.
            plt.clf()
            try:
                scale = P.G.scale
            except AttributeError:
                scale = 10
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
        elif P.S.law == 'von_mises':
            x = zeros((P.S.n,3))
            v = zeros((P.S.n,3))
            gammadot = zeros((P.S.n))
            yieldfunction = zeros((P.S.n))
            pressure = zeros((P.S.n))
            dev_stress = zeros((P.S.n))
            dev_stress_dot = zeros((P.S.n))
            sigmah = zeros((P.S.n))
            sigmav = zeros((P.S.n))
            Rsigmah = zeros((P.R.n))
            Rsigmav = zeros((P.R.n))
            for i in xrange(P.S.n):
                x[i] = L.S[i].x
                v[i] = L.S[i].v
                gammadot[i] = L.S[i].gammadot/P.dt
                sigmav[i] = L.S[i].sigmav
                sigmah[i] = L.S[i].sigmah
                yieldfunction[i] = L.S[i].yieldfunction
                pressure[i] = L.S[i].pressure
                dev_stress[i] = norm(L.S[i].dev_stress)
                dev_stress_dot[i] = norm(L.S[i].dev_dstress)/P.dt
            if P.R.n > 0:
                Rx = zeros((P.R.n,3))
                for i in xrange(P.R.n):
                    Rx[i] = L.R[i].x
                    Rsigmav[i] = L.R[i].stress[1,1]
                    Rsigmah[i] = L.R[i].stress[0,0]
            size=25.
            plt.clf()
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
            if P.R.n > 0:
                plt.scatter(Rx[:,0],Rx[:,1],s=size/2.,marker='s',
                        edgecolor='None',facecolor='r')
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
            if P.R.n == 0:
                plt.scatter(x[:,0],x[:,1],s=size,c=sigmav,marker='s',edgecolor='None')
            else:
                plt.scatter(x[:,0],x[:,1],s=size,c=sigmav,marker='s',
                        edgecolor='None', vmin=min(sigmav.min(),Rsigmav.min()),
                        vmax=max(sigmav.max(),Rsigmav.max()))
                plt.scatter(Rx[:,0],Rx[:,1],s=size/2.,c=Rsigmav,marker='o',
                            edgecolor='None', vmin=min(sigmav.min(),Rsigmav.min()),
                        vmax=max(sigmav.max(),Rsigmav.max()))
            plt.colorbar()
            plt.subplot(335)
            plt.xlabel(r'$\sigma_h$',rotation='horizontal')
            if  P.R.n == 0:
                plt.scatter(x[:,0],x[:,1],s=size,c=sigmah,marker='s',edgecolor='None')
            else:
                plt.scatter(x[:,0],x[:,1],s=size,c=sigmah,marker='s',
                            edgecolor='None', vmin=min(sigmah.min(),Rsigmah.min()),
                            vmax=max(sigmah.max(),Rsigmah.max()))
                plt.scatter(Rx[:,0],Rx[:,1],s=size/2.,c=Rsigmah,marker='o',
                            edgecolor='None', vmin=min(sigmah.min(),Rsigmah.min()),
                            vmax=max(sigmah.max(),Rsigmah.max()))
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
        elif P.S.law == 'dp':
            x = zeros((P.S.n,3))
            v = zeros((P.S.n,3))
            p = zeros((P.S.n))
            q = zeros((P.S.n))
            y = zeros((P.S.n))
            for i in xrange(P.S.n):
                x[i] = L.S[i].x
                v[i] = L.S[i].v
                p[i] = L.S[i].p
                q[i] = L.S[i].q
                y[i] = L.S[i].y
            size=25.
            plt.clf()
            try:
                scale = 2.*P.G.scale
            except AttributeError:
                scale = 10
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

        if name:
            self.savefig(P,str(name).zfill(5))
        else:
            self.savefig(P,str(P.save).zfill(5))
        P.save += 1
        
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
        self.savefig(P,'ns_' + str(P.S.n))
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
        x = zeros((P.S.n,3))
        v = zeros((P.S.n,3))
        s = zeros((P.S.n))
        for i in xrange(P.S.n):
            x[i] = L.S[i].x
            v[i] = L.S[i].v
            s[i] = L.S[i].s_bar
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
