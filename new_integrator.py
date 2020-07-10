import os, glob, sys
import numpy as np
import matplotlib.pyplot as plt

def limiter(a,b):
    return minmod(a,b) # more diffusive
    # return superbee(a,b) # less diffusive
    # return vanLeer(a,b)
    # return vanAlbada1(a,b)

def superbee(a,b): return maxmod(minmod(a,2.*b),minmod(2.*a,b))
def maxmod(a,b): return 0.5*(np.sign(a) + np.sign(b)) * np.maximum(np.abs(a),np.abs(b))
def minmod(a,b): return 0.5*(np.sign(a) + np.sign(b)) * np.minimum(np.abs(a),np.abs(b))
def vanLeer(a,b):
    r = div0(a,b)
    return (r + np.abs(r))/(1 + np.abs(r))
def vanAlbada1(a,b):
    r = div0(a,b)
    return (r**2 + r)/(r**2 +1)
def minmod2(a,b,c,theta=2):
    # theta=1 - no increase of total variation
    # theta=2 - least dissipative
    retval = np.zeros_like(a)
    positive_values = (a>0)*(b>0)*(c>0)
    negative_values = (a<0)*(b<0)*(c<0)
    retval[positive_values] = np.minimum(theta*a,b,theta*c)[positive_values]
    retval[negative_values] = np.maximum(theta*a,b,theta*c)[negative_values]
    return retval

def div0( a, b ):
    """This function replaces nans with zeros when dividing by zero.

    :param a: array - Numerator
    :param b: array - Demoniator
    :returns: array - a/b, with infs and nans replaced by 0

    """
    with np.errstate(divide='ignore', invalid='ignore'):
        c = np.true_divide( a, b )
        c[ ~ np.isfinite( c )] = 0  # -inf inf NaN
    return c


def roll(c,step,ax):
    '''
    A step of +1 gives c_{j-1}. A step of -1 gives c_{j+1}
    '''
    return np.roll(c,step,ax)
#     if ax == 0:
#         if step == 1:
#             return np.vstack([c[0],c[:-1]]) ## <BC TAG>
#         elif step == -1:
#             return np.vstack([c[1:],c[-1]]) ## <BC TAG>
#     # elif ax == 1: # THIS DIRECTION UNTESTED
#         # if step == 1:
#             # return np.hstack([c[:,0],c[:,:-1]])
#         # elif step == -1:
#             # return np.hstack([c[:,1:],c[:,-1]])
#     else: print('WARNING: THIS AXIS NOT IMPLEMENTED.')

def flux(c,v):
    flux = c*velocity(c,v)
    return flux

def velocity(c,v):
    # vel = v*(1-c) # simple segregation model
    R = 2.5 # size ratio
    vel = v*(1-1/(c + (1-c)*R))
    return vel

def KT(c,v,dx,dt,ax):
    cx = limiter((c - roll(c,1,ax))/dx, (roll(c,-1,ax) - c)/dx)
    # cx = minmod2((roll(c,-1,ax) - c)/dx,
    #              (roll(c,-1,ax) - roll(c,1,ax))/(2*dx),
    #              (c - roll(c,1,ax))/dx
    #             )
    cplusleft   = c - dx/2*cx
    cminusright = c + dx/2*cx
    cplusright  = roll(cplusleft,-1,ax)
    cminusleft  = roll(cminusright,1,ax)

    vleft  = roll(v,1,ax)
    vright = roll(v,-1,ax)

    aright = np.maximum(np.abs(velocity(cminusright,vright)),np.abs(velocity(cplusright,vright)))
    aleft  = np.maximum(np.abs(velocity( cminusleft, vleft)),np.abs(velocity( cplusleft, vleft)))

    RHS = -(( flux( cplusright,vright) +
              flux(cminusright,vright) -
              flux(  cplusleft, vleft) -
              flux( cminusleft, vleft) ) -
              ( aright*(cplusright - cminusright) -
                 aleft*(cplusleft  -  cminusleft) )
            )/(2*dx)
    return RHS

def pad(c,v):
    # constant default value is zero
    C = np.pad(c,1,mode='edge') # this doesnt appear to matter at all
    V = np.pad(v,1,mode='constant') # zero velocity outside of domain - this sets up the no flux BC!
    # print(V[-1])
    return C,V

def BC(c,v): # apply no flux boundaries  ## <BC TAG>
    # KIND OF WORKS
    # c[:padwidth] = 0
    # c[-padwidth:] = 1

    # v[:padwidth] = 0
    # v[-padwidth:] = 0
    # c[1] = div0(c[2]*v[2],v[1])
    # c[0] = div0(c[1]*v[1],v[0])
    # c[-2] = div0(c[-3]*v[-3],v[-2])
    # c[-1] = div0(c[-2]*v[-2],v[-1])
    # v[1] = -div0(c[2]*v[2],c[1])
    # v[0] = 0#-div0(c[1]*v[1],c[0])
    # v[-2] = -div0(c[-3]*v[-3],c[-2])
    # v[-1] = 0#-div0(c[-2]*v[-2],c[-1])
    # print(c[-1,0])
    # c[0] = c[2]
    # c[-1] = c[-3]
    # v[0] = 0
    # v[1] = 0
    # v[2] = 0
    # v[3] = 0
    # v[-2] = 0
    # v[-1] = 0
    # v[1] = -div0(c[2]*v[2],c[1])
    # v[0] = -div0(c[1]*v[1],c[0])
    # v[-2] = -div0(c[-3]*v[-3],c[-2])
    # v[-1] = -div0(c[-2]*v[-2],c[-1])
    # c[-1] = 1
    # c[0] = c[1] = 0
    # c[-1] = c[-2] = 1
    return c, v

def RK4(C,V,dx,dt,ax):
    c,v = pad(C,V)
    k1 = KT(c,        v,dx,dt,ax)
    k2 = KT(c+dt/2*k1,v,dx,dt,ax)
    k3 = KT(c+dt/2*k2,v,dx,dt,ax)
    k4 = KT(c+dt*k3,  v,dx,dt,ax)
    dc = dt/6*(k1+2*k2+2*k3+k4)
    return dc[1:-1,1:-1]

def RK3(C,V,dx,dt,ax):
    c,v = pad(C,V)
    c1 =         c +               dt*KT( c,v,dx,dt,ax)
    c2 =    0.75*c +    0.25*(c1 + dt*KT(c1,v,dx,dt,ax))
    c3 = 0.33333*c + 0.66667*(c2 + dt*KT(c2,v,dx,dt,ax))
    return c3[1:-1,1:-1]

def Euler(c,v,dx,dt,ax):
    c,v = pad(c,v)
    dc = dt*KT(c,v,dx,dt,ax)
    return dc[1:-1,1:-1]

def diffusion(c,D,dx,dt,ax):
    dDc_dy = np.gradient(D*c,dx,axis=ax)
    # dc_dy[0] = 0
    # dc_dy[-1] = 0
    d2Dc_dy2 = np.gradient(dDc_dy,dx,axis=ax)
    return c + dt*d2Dc_dy2

def main():
    nx = 1
    ny = 201
    L = 1.0
    # c = 0.5*np.ones([ny,nx])

    c = np.zeros([ny,nx])
    c[ny//4:3*ny//4] = 1.0

    v = -np.ones_like(c)
    CFL = 0.025
    dx = dy = L/(np.maximum(nx,ny)-1)
    dt = CFL*dx*4/(np.max(np.abs(v)))
    t_max = 2.0
    nt = np.int(t_max/dt)
    D = 5e-3
    # nt = 10
    # c_old = c.copy()
    for i in range(nt):
        u = np.zeros_like(v)
        # c_old = c.copy()
        # c += Euler(c,v,dy,dt,ax=0)
        # c += RK4(c_pad.T,v_pad.T,dx,dt,ax=0).T
        c = RK3(c,v,dy,dt,ax=0)
        c = diffusion(c,D,dy,dt,ax=0)

        if i%(nt//10)==0:
            plt.plot(c[:,0])
        print('    t = ' + str(i*dt) + '                ', end='\r')
    plt.show()

if __name__=='__main__':
    main()
    print('\nAll done.')
