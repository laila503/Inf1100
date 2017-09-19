# -*- coding: utf-8 -*-
"""
Created on Sat Sep 16 14:13:10 2017

@author: laila
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import ArtistAnimation
from mpl_toolkits.mplot3d import Axes3D

"""2.3 Implementation"""

#Scalar (pointvise) solution:
def solver(I, V, q, b, f, c, Lx, Ly, Nx, Ny, dt, T, animate=False,frames=10):
    
    if animate:
        frame = []
        fig = plt.figure()
                   
    dx = Lx/float(Nx)
    dy = Ly/float(Ny)
    
    C2x = (dt/dx)**2 
    C2y = (dt/dy)**2

    Nt = int(T/dt)
    
    const = 1/(1+b*dt/2.)
    
    x = np.linspace(0,Lx,Nx)
    y = np.linspace(0,Ly,Ny)
    time = np.linspace(0,T,Nt)
    
    #u = u_n+1, up = u_n, upp = u_n-1
    u = np.zeros((int(Nx)+2,int(Ny)+2))
    up = np.zeros((int(Nx)+2,int(Ny)+2))
    upp = np.zeros((int(Nx)+2,int(Ny)+2))
    
    Ix = np.arange(0, u.shape[0])
    Iy = np.arange(0, u.shape[1])
    It = np.arange(0, time.shape[0])
    
    def fill_ghost(u_):
        u_[0,:] = np.copy(u_[2,:])
        u_[-1,:] = np.copy(u_[-3,:])
        u_[:,0] = np.copy(u_[:,2])
        u_[:,-1] = np.copy(u_[:,-3])
        
    #Initialconditions
    for i in Ix[1:-1]:
        for j in Iy[1:-1]:
            upp[i,j] = I(i*dx,j*dy)
       
    fill_ghost(upp)
    
    #Boundary conditions
    for i in Ix[1:-1]:
        for j in Iy[1:-1]:
            coeff_x = C2x*(q((i+0.5)*dx,j*dy)*(upp[i+1,j] - I(i*dx,j*dy)) - q((i-0.5)*dx,j*dy)*(I(i*dx,j*dy) - upp[i-1,j]))
            coeff_y = C2y*(q((i)*dx,(j+0.5)*dy)*(upp[i,j+1] - I(i*dx,j*dy)) - q((i)*dx,(j-0.5)*dy)*(I(i*dx,j*dy) - upp[i,j-1]))
            up[i,j] = 0.5*(2*I(i*dx,j*dy) + 2*dt*V(i*dx,j*dy)*(1-b*dt/2.) +coeff_x + coeff_y + dt**2*f(i*dx,j*dy,0))
    
    fill_ghost(up)   
    
    #Main scheme
    for n in It:
        for i in Ix[1:-1]:
            for j in Iy[1:-1]:
               coeff_x = C2x*(q((i+0.5)*dx,j*dy)*(up[i+1,j] - up[i,j]) - q((i-0.5)*dx,j*dy)*(up[i,j] - up[i-1,j]))
               coeff_y = C2y*(q((i)*dx,(j+0.5)*dy)*(up[i,j+1] - up[i,j]) - q((i)*dx,(j-0.5)*dy)*(up[i,j] - up[i,j-1])) 
               u[i,j] = const*(2*up[i,j] - upp[i,j]*(1-b*dt/2.) + coeff_x + coeff_y + dt**2*f(i*dx,j*dy,n*dt))
        fill_ghost(u)
        upp = np.copy(up)
        up = np.copy(u)

        if animate and n%frames == 0:
            frame.append(plt.plot(upp,"r-"))
    
    if animate:
        an = ArtistAnimation(fig,frame,interval=1/dt,blit=True)
        plt.show() 
        
    else:
        plt.plot(x,u[1:-1,1], label="With constant y=1")
        plt.xlabel("x")
        plt.ylabel("u")
        plt.hold("on")
        plt.plot(x,u[1,1:-1], label="With constant x=1")
        plt.xlabel("y")
        plt.ylabel("u")
        plt.title("Plot of solution")
        plt.show()
      
    return u[1:-1,1:-1],x,y,time

#Vectorized solution:
def solver_vec(I, V, q, b, f, c, Lx, Ly, Nx, Ny, dt, T, animate=False,frames=20):    

    if animate:
        frame = []
        fig = plt.figure()
    
    def make_f_array(f,x,y,t):
        if t.shape[0] != 0:
            a = np.zeros((x.shape[0],y.shape[0],t.shape[0]))
            for n in range(a.shape[2]):
                for i in range(a.shape[0]):
                    for j in range(a.shape[1]):
                        a[i,j,n] = f(x[i],y[j],t[n])
        return a
                        
    def make_func_to_array(f,x,y):
        a = np.zeros((x.shape[0],y.shape[0]))
        for i in range(a.shape[0]):
                for j in range(a.shape[1]):
                    a[i,j] = f(x[i],y[j])
        return a
        
    def make_q_array(q,x,y,dx,dy):
        #Makes q double the size, to be able to get half steps
        a = np.zeros((2*x.shape[0]+2,2*y.shape[0]+2))
        for i in range(a.shape[0]):
                for j in range(a.shape[1]):
                    a[i,j] = q((i-1)*dx/2.,(j-1)*dy/2.)
                    
        return a
                           
    dx = Lx/float(Nx)
    dy = Ly/float(Ny)
    
    C2x = (dt/dx)**2 
    C2y = (dt/dy)**2

    Nt = int(T/dt)
    
    const = 1/(1+b*dt/2.)
    
    x = np.linspace(0,Lx,Nx)
    y = np.linspace(0,Ly,Ny)
    time = np.linspace(0,T,Nt)
    
    #u = u_n+1, up = u_n, upp = u_n-1
    u = np.zeros((int(Nx)+2,int(Ny)+2))
    up = np.zeros((int(Nx)+2,int(Ny)+2))
    upp = np.zeros((int(Nx)+2,int(Ny)+2))
    
    Ix = np.arange(0, u.shape[0])
    Iy = np.arange(0, u.shape[1])
    It = np.arange(0, time.shape[0])
    
    f_array = make_f_array(f,x,y,time)
    I_array = make_func_to_array(I,x,y)
    V_array = make_func_to_array(V,x,y) 
    q_array = make_q_array(q,x,y,dx,dy) 
    
    def fill_ghost(u_):
        u_[0,:] = np.copy(u_[2,:])
        u_[-1,:] = np.copy(u_[-3,:])
        u_[:,0] = np.copy(u_[:,2])
        u_[:,-1] = np.copy(u_[:,-3])
        
    #Initialconditions vectorized:
    upp[1:-1,1:-1] = I_array
   
    fill_ghost(upp)

    #Boundary conditions
    coeff_x = C2x*(q_array[2*Ix[1:-1] +1 ,Iy[1:-1]]*(upp[2:,1:-1] - I_array)  - q_array[2*Ix[1:-1] -1 ,Iy[1:-1]]*(I_array - upp[:-2,1:-1]))
    coeff_y = C2y*(q_array[2*Iy[1:-1],Iy[1:-1]+1]*(upp[1:-1,2:] - I_array) - q_array[2*Iy[1:-1],Iy[1:-1]-1]*(I_array - upp[1:-1,:-2]))
    up[1:-1,1:-1] = 0.5*(2*I_array + 2*dt*V_array*(1-b*dt/2.) +coeff_x + coeff_y + dt**2*f_array[:,:,0])
        
    fill_ghost(up)   
    
    #Main scheme:
    for n in It:       
        coeff_x = C2x*(q_array[2*Ix[1:-1] +1 ,Iy[1:-1]]*(up[2:,1:-1] - up[1:-1,1:-1])  -  q_array[2*Ix[1:-1] -1 ,Iy[1:-1]]*(up[1:-1,1:-1] - up[:-2,1:-1]))
        coeff_y = C2y*(q_array[2*Iy[1:-1],Iy[1:-1]+1]*(up[1:-1,2:] - up[1:-1,1:-1]) - q_array[2*Iy[1:-1],Iy[1:-1]-1]*(up[1:-1,1:-1] - up[1:-1,:-2]))
        u[1:-1,1:-1] = const*(2*up[1:-1,1:-1] - upp[1:-1,1:-1]*(1-b*dt/2.) +coeff_x + coeff_y + dt**2*f_array[:,:,n])
               
        fill_ghost(u)
        upp = np.copy(up)
        up = np.copy(u)

        if animate and n%frames == 0:
            frame.append(plt.plot(up[:,4],"r-"))

    if animate:
        an = ArtistAnimation(fig,frame,interval=(dt),blit=True)
        plt.show()
        
    else:
        plt.plot(x,u[1:-1,1], label="With constant y=1")
        plt.xlabel("x")
        plt.ylabel("u")
        plt.hold("on")
        plt.plot(x,u[1,1:-1], label="With constant x=1")
        plt.xlabel("y")
        plt.ylabel("u")
        plt.title("Plot of solution")
        plt.show()
  
    return u[1:-1,1:-1],x,y,time

"""3.1 Constant solution """
def test_case():
    c = 5  
    u,x,y,time = solver_vec(I=lambda x,y:c,V=lambda x,y:0, q=lambda x,y:0, b=10, f=lambda x,y,t:0, c=c, Lx=1, Ly=1, Nx=1e2, Ny=1e2, dt=0.001, T=0.5,animate=False,frames=20)
    #u_const = lambda x,y,time:c*2
    u_const = lambda x,y,time:c
    tol = 1E-12
    diff = np.max(np.abs(u-u_const(x,y,time)))
    assert diff<tol, "Error in contant solution"   

"""3.3 Exact 1D plug-wave solution in 2D"""
def plug_wave_test():
    Lx=1.0; Ly=1.0
    Nx=9; Ny=9; T=1.0

    x = np.linspace(0,Lx,Nx+1)
    y = np.linspace(0,Ly,Ny+1)
 
    dx = x[1] - x[0]
    dt = dx  
    
    c =1.0; b=0.0

    #Defining the functions:
    def I(x, y):
        if y <= 0.6 and  y >= 0.4:
            I = 1
        else:
            I = 0
        return I

    V = lambda x, y: 0.0    
    f = lambda x, y, t: 0.0   
    q = lambda x, y: c #np.ones((Nx+1,Ny+1))*c
   
    u_num,x,y,t = solver_vec(I=I,V=V, q=q, b=0, f=f, c=c, Lx=Lx, Ly=Ly, Nx=Nx, Ny=Ny, dt=dt, T=T, animate=False, frames=20)
    u_exact = np.array([[I(x, y_) for y_ in y] for x_ in x])
    print np.size(u_num), np.size(u_exact)
    print u_num, u_exact
    tol = 1E-12
    step_num = t/dt
    diff = np.abs(u_exact - u_num).max()
    for i in t:
        msg = 'diff=%g at time=%g' % (diff, t[i])
        assert diff < tol, msg

    print 'Exact at initial step:\n', u_exact
    print 'Computed at time %s step:\n'% t, u_num

"""3.4 Standing, undamped waves"""
def standing_undamped():
    mx = 3.0; my = 3.0
    Lx = 1.0; Ly = 1.0
    dt = 0.001
    Nx = 100 ; Ny = 100; 
    
    # wavenumbers
    ky = (my*np.pi)/Ly
    kx = (mx*np.pi)/Lx

    c = 2.0; b = 0.0
    q_const = 0.4; A = 1.0
    T = 0.5
    #Nt = int(T/dt); t = np.linspace(0,T,Nt)
    
    w = c*np.sqrt(kx**2 + ky**2)# angular frequency

    def u_e(x,y,t):
        return A*np.cos(kx*x)*np.cos(ky*y)*np.cos(w*t)
    
    def q(x,y):
        return q_const

    def f(x,y,t):
       return np.exp(-b*t)*np.cos(kx*x)*np.cos(ky*y)*(np.cos(w*t)*(q_const*kx**2 + q_const*ky**2 - w**2) + b*w*np.sin(w*t))

    def I(x,y):
        return np.cos(mx*x*np.pi/Lx)*np.cos(my*y*np.pi/Ly)

    def V(x,y):
        return -b*np.cos(mx*x*np.pi/Lx)*np.cos(my*y*np.pi/Ly)
   
    
    u_num,x,y,t = solver_vec(I=I,V=V, q=q, b=b, f=f, c=c, Lx=Lx, Ly=Ly, Nx=Nx, Ny=Ny, dt=dt, T=T,animate=True,frames=5)

    
    h_values = [0.1*2**(-i) for i in range(7)]
    E_values = []

    for h in h_values:
        dt = h/2.
        Nx = 1.0/h; Ny = 1.0/h
        u_ex = u_e(x, y, T)
        #e = np.sqrt(np.sum((u_e - u[1:-1,1:-1])**2))
        e = np.abs(u_ex - u_num[:,:]).max()
        #E = np.sqrt(dx*dy*e)
        E_values.append(e)


    def compute_rates(h_values, E_values):
        """compute convergence rates"""
        m = len(h_values)
        r = [np.log(E_values[i-1]/E_values[i])/
             np.log(h_values[i-1]/h_values[i])
             for i in range(0, m, 1)]
        # Round to two decimals
        r = [round(r_, 2) for r_ in r]
        return r
        
    print 'E: %s' % E_values
    r = compute_rates(h_values, E_values)
    print 'r: %s' % r
    print 'h', h_values
    
"""4 Investigate a physical problem"""    
def physical_problem(hill, movie):
    ## Defining variables:
    T = 0.5
    dt = 0.001
    Lx = 1.0; Ly = 1.0
    Nx = 100; Ny = 100
    g = 9.81     # Gravity
    b = 0.0      # Initally no damping
    f = lambda x, y, t: 0.0
    V = lambda x, y: 0.0


    """ Different bottoms """
    ## Gaussian hill:
    if hill == 'gaussian':
        I0 = 1.1     # Ocean depth
        Ia = 0.4
        Im = Lx/4    # Loaction of peak, Start in origo
        Is = 0.1     # The width of function
        B0 = 0.0     # Ocean floor
        Ba = 0.7     # Height of top
        Bmx = Lx/2.0 # Where top is located at x axis
        Bmy = Ly/2.0 # Where top is located at y axis
        Bs = 0.2     # How steep it is
        b = 1.0      # With damping
    
        I = lambda x, y: I0 + Ia*np.exp(-((x - Im)/Is)**2)
        B = lambda x, y: B0 + Ba*np.exp(-((x - Bmx)/Bs)**2 - ((y - Bmy)/(b*Bs))**2)

    ## Cosine hat:
    elif hill == 'cosine_hat': # Steep hill
        I0 = 1.1; Ia = 0.4 ; Im = 0.0 ; Is = 0.1
        B0 = 0.0; Ba = 0.7; Bmx = Lx/2.0; Bmy = Ly/2.0 ; Bs = 0.1
        I = lambda x, y: I0 + Ia*np.exp(-((x - Im)/Is)**2)
        B = np.vectorize(lambda x, y: B0 + Ba*np.cos(np.pi*(x - Bmx)/(2*Bs))*np.cos(np.pi*(y - Bmy)/(2*Bs)) \
                         if 0 <= np.sqrt((x - Bmx)**2+(y - Bmy)**2) <= Bs else B0)

    ## Box:
    elif hill == 'box':
        I0 = 1.1; Ia = 0.4 ; Im = 0.0; Is = 0.1
        B0 = 0.0; Ba = 0.7; Bmx = Lx/2.0; Bmy = Ly/2.0 ; Bs = 0.1; b = 1.0
        I = lambda x, y: I0 + Ia*np.exp(-((x - Im)/Is)**2)
        B = np.vectorize(lambda x, y: B0 + Ba \
                         if Bmx - Bs <= x <= Bmx + Bs and Bmy-b*Bs <=y<= Bmy+b*Bs else B0)

    q = lambda x, y: np.sqrt(g*(I(x,y)-B(x,y)))

    ## Calling solver function:
    if movie == True:
        u_num,x,y,t = solver_vec(I=I,V=V, q=q, b=b, f=f, c=1, Lx=Lx, Ly=Ly, Nx=Nx, Ny=Ny, dt=dt, T=T,animate=True,frames=2)

    else:
        u_num,x,y,t = solver_vec(I=I,V=V, q=q, b=b, f=f, c=1, Lx=Lx, Ly=Ly, Nx=Nx, Ny=Ny, dt=dt, T=T,animate=False,frames=20)

   
if __name__=="__main__":    
#    solver(I=lambda x,y:x,V=lambda x,y:0, q=lambda x,y:1 , b=10, f=lambda x,y,t:0, c=0.5, Lx=1, Ly=1, Nx=100, Ny=100, dt=0.001, T=2,animate=True,frames=20)
#    solver_vec(I=lambda x,y:x,V=lambda x,y:0, q=lambda x,y:1, b=10, f=lambda x,y,t:0, c=0.5, Lx=1, Ly=1, Nx=10, Ny=10, dt=0.001, T=3,animate=False,frames=20)
    test_case()
#    plug_wave_test()
#    standing_undamped()
#    physical_problem('gaussian',True)
#    physical_problem('cosine_hat',True)
#    physical_problem('box',True)
#    pass

    



    