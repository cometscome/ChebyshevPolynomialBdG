from scipy.sparse import lil_matrix
from numpy.linalg import norm
import numpy as np
from multiprocessing import Pool
from multiprocessing import Process
from scipy import linalg
import time



def xy2i(ix,iy,Nx):
    ii = (iy)*Nx+ix
    return ii

def init_delta(Nx,Ny,delta):
    vec_delta = lil_matrix((Nx*Ny,Nx*Ny))
    for ii in range(Nx*Ny):
        vec_delta[ii,ii] = delta
    return vec_delta


def calc_A(Nx,Ny,mu,vec_delta,aa):
    Ln = Nx*Ny*2
    A = lil_matrix((Ln,Ln))
#    A.setdiag(-mu)
    for ix in range(Nx):        
        for iy in range(Ny):
            ii = xy2i(ix,iy,Nx)
            jx = ix
            jy = iy
            jj = xy2i(jx,jy,Nx)
            A[ii,jj] = -mu
#            print(ii)

# +1 in x direction
            jx = ix +1
            if jx == Nx:
                jx = 0            
            jy = iy
            jj = xy2i(jx,jy,Nx)
            A[ii,jj] = -1.0
            
# -1 in x direction
            jx = ix -1
            if jx == -1:
                jx = Nx-1           
            jy = iy
            jj = xy2i(jx,jy,Nx)
            A[ii,jj] = -1.0

# + 1 in y direction
            jx = ix 
            jy = iy +1
            if jy == Ny:
                jy = 0         
            jj = xy2i(jx,jy,Nx)
            A[ii,jj] = -1.0   

# -1 in y direction
            jx = ix 
            jy = iy -1
            if jy == -1:
                jy = Ny-1
            jj = xy2i(jx,jy,Nx)
            A[ii,jj] = -1.0
    for ii in range(Nx*Ny):
#        print(ii)
        for jj in range(Nx*Ny):
#            print("jj",jj)
            A[ii+Nx*Ny,jj+Nx*Ny] = -A[ii,jj]
            A[ii,jj+Nx*Ny] = vec_delta[ii,jj]
            A[ii+Nx*Ny,jj] = vec_delta[jj,ii]
    
    A = A/aa

#    print(A)
    A = A.tocsr()
#    print(A)
    return A

def calc_A2(A2,Nx,Ny,mu,vec_delta,aa):
    Ln = Nx*Ny*2
    A = A2.tolil()
    A = A*aa
    for ii in range(Nx*Ny):
#        print(ii)
        for jj in range(Nx*Ny):
#            print("jj",jj)
            A[ii,jj+Nx*Ny] = vec_delta[ii,jj]
            A[ii+Nx*Ny,jj] = vec_delta[jj,ii]
    
    A = A/aa

#    print(A)
    A = A.tocsr()
#    print(A)
    return A


def iteration(nc,Nx,Ny,aa,bb,omegac,U,full):
    Ln = Nx*Ny*2
    vec_delta = init_delta(Nx,Ny,0.1)
    vec_delta_old = vec_delta

    A = calc_A(Nx,Ny,1e-12,vec_delta,10.0)
    itemax = 100

    for ite in range(itemax):
        if full == False:
            vec_delta = calc_meanfields(nc,A,Nx,Ny,Ln,aa,bb,omegac)    
        else:
            vec_delta = calc_meanfields_full(nc,A,Nx,Ny,Ln,aa,bb,omegac)    

        vec_delta = vec_delta*U

#        A = calc_A(Nx,Ny,1e-12,vec_delta,10.0)
        A = calc_A2(A,Nx,Ny,1e-12,vec_delta,10.0)
        eps = 0.0
        nor = 0.0

        for i in range(Nx*Ny):
            eps += abs(vec_delta[i,i] -vec_delta_old[i,i])**2
            nor += abs(vec_delta_old[i,i])**2
        eps = eps/nor
        print "ite = ",ite,eps
        if eps <= 1e-6:
            print "End",vec_delta[Nx/2,Ny/2]
            break

        vec_delta_old = vec_delta

    return vec_delta
#        print(vec_delta)



        


def calc_meanfields(nc,A,Nx,Ny,Ln,aa,bb,omegac):    
    vec_delta = lil_matrix((Nx*Ny,Nx*Ny))
    for ix in range(Nx):
        for iy in range(Ny):
            ii = xy2i(ix,iy,Nx)
            jj = ii + Nx*Ny
            right_j = jj
            left_i = ii
            vec_ai = calc_polynomials(nc,left_i,right_j,A,Ln)
            density = calc_meanfield(vec_ai,aa,bb,omegac,nc)
            vec_delta[ii,ii] = density
    return vec_delta


def calc_meanfields_full(nc,A,Nx,Ny,Ln,aa,bb,omegac):    
    a = A.todense()*aa
    w, v = linalg.eigh(a, lower=True, eigvals_only=False, overwrite_a=False, eigvals=None, check_finite=True)    
    vec_delta = lil_matrix((Nx*Ny,Nx*Ny))
    for ix in range(Nx):
        for iy in range(Ny):
            ii = xy2i(ix,iy,Nx)
            jj = ii + Nx*Ny
            delta = 0.0

            for i in range(Nx*Ny*2):
                if w[i] <= 0.0:
                    if abs(w[i]) <= omegac:
                        delta += v[ii,i]*v[jj,i]
            vec_delta[ii,ii] = delta

    return vec_delta


def calc_meanfield(vec_ai,aa,bb,omegac,nc):    
    ba = np.arccos(-bb/aa)
    omeb = np.arccos(-(omegac+bb)/aa)
    pi = np.arctan(1.0)*4
    density = 0.0
    for j in range(nc-1):
        i = j + 1
        density +=  vec_ai[i]*(np.sin(i*omeb)-np.sin(i*ba))/i
    density += vec_ai[0]*(omeb-ba)/2.0

    density = density*2/pi

    return density



def calc_polynomials(nc,left_i,right_j,A,Ln):
    vec_jnmm = np.zeros(Ln)
    vec_jnm = np.zeros(Ln)
    vec_jn = np.zeros(Ln)
    vec_jn[right_j] = 1.0
    vec_ai = np.zeros(nc)

    for nn in range(nc):
        if nn == 0:
            vec_jnmm = 0.0
            vec_jnm = 0.0
            vec_jn[right_j] = 1.0
        elif nn == 1:
            vec_jn = A.dot(vec_jn)
        else:
            vec_jn = 2.0*A.dot(vec_jnm) - vec_jnmm
        vec_ai[nn] = vec_jn[left_i]
        vec_jnmm = vec_jnm
        vec_jnm = vec_jn
    return vec_ai



def main():
#    for i in range(1,6):
#        print("test")

    nc = 1000
    nx = 20
    ny = 20
    vec_delta = init_delta(nx,ny,0.1)
    A = calc_A(nx,ny,1e-12,vec_delta,10.0)

    a = A.todense()
    w, v = linalg.eigh(a, lower=True, eigvals_only=False, overwrite_a=False, eigvals=None, check_finite=True)
#    print(w 

    vec_ai = calc_polynomials(nc,1,1+nx*ny,A,nx*ny*2)
#    print(vec_ai)
    density = calc_meanfield(vec_ai,10.0,0.0,10.0,nc)
#    print(density)

    aa = 10.0
    bb = 0.0
    omegac = 10.0
    Ln = nx*ny*2
    
 #   vec_delta = calc_meanfields(nc,A,nx,ny,Ln,aa,bb,omegac)    
#    print(vec_delta)
    U = -2.0

    start = time.time()
    iteration(nc,nx,ny,aa,bb,omegac,U,True)
    elapsed_time = time.time()-start
    print "elapsed_time",elapsed_time 


    start = time.time()
    iteration(nc,nx,ny,aa,bb,omegac,U,False)
    elapsed_time = time.time()-start
    print "elapsed_time",elapsed_time


#    x = np.zeros((nx*ny*2))
#    for i in range(nx*ny*2):
#        x[i] = i + 1
#    print(x)
#    x = matrix([[1],[1],[1],[2],[2],[3]])

#    y = A.dot(x)
#    print(y)

    

if __name__ == "__main__":
    main()
