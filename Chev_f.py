from scipy.sparse import lil_matrix
from numpy.linalg import norm
from ctypes import *
import numpy as np
from multiprocessing import Pool
from multiprocessing import Process
from scipy import linalg

import time

def calc_meanfields_f(nc,A,Nx,Ny,Ln,aa,bb,omegac):    
    vec_delta = lil_matrix((Nx*Ny,Nx*Ny))
    subchev = np.ctypeslib.load_library("sub_chev.so",".")
#calcsgap(i,val,row,col,valsize,N,aa,bb,nc,omegac,sgap)
#  integer(8),intent(in)::i,N,valsize,nc
#  real(8),intent(in)::aa,bb,omegac
#  real(8),dimension(0:valsize-1),intent(in)::val
#  integer(8),dimension(0:valsize-1),intent(in)::col
#  integer(8),dimension(0:N),intent(in)::row
#  real(8),intent(out)::sgap
    subchev.calcsgap_.argtypes = [
        POINTER(c_int32),
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.int32),
        np.ctypeslib.ndpointer(dtype=np.int32),
        POINTER(c_int32),
        POINTER(c_int32),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_int32),
        POINTER(c_double),
        POINTER(c_double),
        ]
    subchev.calcsgap_.restype = c_void_p

    val = A.data
    col =  A.indices
    row =  A.indptr
    valsize = val.size
    N = Ln
#calcsgap(i,val,row,col,valsize,N,aa,bb,nc,omegac,sgap)

    valsizec = c_int32(valsize)
#    valsize = byref(c_int32(valsize))
    N_c = c_int32(N)
#    N = byref(c_int32(N))
    aac = c_double(aa)
#    aa = byref(c_double(aa))
    bbc = c_double(bb)
#    bb = byref(c_double(bb))
    omegacc = c_double(omegac)
#    omegac = byref(c_double(omegac))
    sgap = 0.0
#    sgap = byref(c_double(sgap))
    sgapc =c_double(sgap)
    ncc = c_int32(nc)
#    nc = byref(c_int32(nc))

    for ix in range(Nx):
        for iy in range(Ny):
            ii = xy2i(ix,iy,Nx)
            i = ii
            ic = c_int32(i)

#            subchev.calcsgap_(i,val,row,col,valsize,N,aa,bb,nc,omegac,sgap)
            subchev.calcsgap_(byref(ic),val,row,col,byref(valsizec),byref(N_c),byref(aac),byref(bbc),byref(ncc),byref(omegacc),byref(sgapc))
#            print "sgap ",sgapc.value

            vec_delta[ii,ii] = sgapc.value

            
    return vec_delta



 
    print "yes!!"
    return vec_delta




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

def calc_A3(A2,Nx,Ny,mu,vec_delta,aa):
    Ln = Nx*Ny*2
#    vec_delta2 = vec_delta.tocsr()

    for i in range(Nx*Ny):
        j = A2.indptr[i]
        j2 = A2.indptr[i+1]
        for j3 in range(j2-j):
            j4 = j + j3
#            print "j4",j4
            jj = A2.indices[j4]
#            print jj

            if jj>Nx*Ny-1:
                A2.data[j4] = vec_delta[i,jj-Nx*Ny]/aa

#    print "(^_-)"

    for i in range(Nx*Ny):
        j = A2.indptr[i+Nx*Ny]
        j2 = A2.indptr[i+1+Nx*Ny]
#        print j,j2
        for j3 in range(j2-j):
            j4 = j + j3
#            print "j4",j4
            jj = A2.indices[j4]
#            print jj
            if jj <Nx*Ny:
                A2.data[j4] = vec_delta[jj,i]/aa                
#            print jj





    return A2



def iteration(nc,Nx,Ny,aa,bb,omegac,U,full,fortran):
    Ln = Nx*Ny*2
    vec_delta = init_delta(Nx,Ny,0.1)
    vec_delta_old = vec_delta

    A = calc_A(Nx,Ny,1e-12,vec_delta,10.0)
    itemax = 100

    for ite in range(itemax):
        if full == False:
            if fortran == False:
                vec_delta = calc_meanfields(nc,A,Nx,Ny,Ln,aa,bb,omegac)    
            else: 
                vec_delta = calc_meanfields_f(nc,A,Nx,Ny,Ln,aa,bb,omegac)    
        else:
            vec_delta = calc_meanfields_full(nc,A,Nx,Ny,Ln,aa,bb,omegac)    

        vec_delta = vec_delta*U

#        A = calc_A(Nx,Ny,1e-12,vec_delta,10.0)
#        A = calc_A2(A,Nx,Ny,1e-12,vec_delta,10.0)
        A = calc_A3(A,Nx,Ny,1e-12,vec_delta,10.0)
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
    nx = 40
    ny = 40
    vec_delta = init_delta(nx,ny,0.1)


    print "yes!"

    A = calc_A(nx,ny,1e-12,vec_delta,10.0)
#    print A.data
#    print A.indices
#    print A.indptr

#    a = A.todense()
#    w, v = linalg.eigh(a, lower=True, eigvals_only=False, overwrite_a=False, eigvals=None, check_finite=True)
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
    iteration(nc,nx,ny,aa,bb,omegac,U,False,False)
    elapsed_time = time.time()-start
    print "Chebyshev: elapsed_time",elapsed_time

    start = time.time()
    iteration(nc,nx,ny,aa,bb,omegac,U,False,True)
    elapsed_time = time.time()-start
    print "Chebyshev with Fortran:elapsed_time",elapsed_time


    start = time.time()
    iteration(nc,nx,ny,aa,bb,omegac,U,True,False)
    elapsed_time = time.time()-start
    print "Full diagonalization:elapsed_time",elapsed_time 




#    x = np.zeros((nx*ny*2))
#    for i in range(nx*ny*2):
#        x[i] = i + 1
#    print(x)
#    x = matrix([[1],[1],[1],[2],[2],[3]])

#    y = A.dot(x)
#    print(y)

    

if __name__ == "__main__":
    main()
