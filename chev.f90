!This code is not finished
!I have not complied this code yet.


 subroutine init_delta(nx,ny,mat_delta,delta)
   implicit none
   integer,intent(in)::nx,ny
   real(8),intent(in)::delta
   real(8),intent(out)::mat_delta(nx*ny)

   mat_delta = 0d0
   mat_delta = delta

   return

 end subroutine init_delta

 function xy2i(ix,iy,Nx)
   implicit none
   integer::xy2i
   integer::ix,iy,Nx

   xy2i = (iy-1)*Nx + ix

 end function xy2i

 subroutine calc_val(n,mat_h,val,row,col)
   implicit none
   integer,intent(in)::n
   real(8),intent(in)::mat_h(:,:)
   real(8),intent(out),allocatable::val(:)
   integer,intent(out),allocatable::col(:),row(:)
   integer,allocatable::vec_i(:)
   integer::i,j,ip,im

   allocate(row(n+1))
   allocate(vec_i(n))


   row(1) = 1
   vec_i = 0

   im = 0
   do i = 1,n
      ip = 0
      do j = 1,n
         if(abs(mat_h(i,j)) .ne. 0d0) then
!            write(*,*) i,j,mat_h(i,j)
            ip = ip + 1
         end if
      end do
      vec_i(i) = ip
!      write(*,*) ip
      im = im + ip
      row(i+1) = im
   end do

   im = sum(vec_i)
   allocate(val(im),col(im))
   write(*,*) im

   do i = 1,n
      ip = 0
      do j = 1,n
         if(abs(mat_h(i,j)) .ne. 0d0) then
            ip = ip + 1
            val(ip) = mat_h(i,j)
            col(ip) = j
         end if
      end do
   end do

   deallocate(vec_i)
      

   return
 end subroutine calc_val

 subroutine calc_h(nx,ny,mat_h,Ln,mu,mat_delta)
   implicit none
   integer,intent(in)::nx,ny,Ln
   real(8),intent(in)::mu,mat_delta(nx*ny)
   real(8),intent(out)::mat_h(1:Ln,1:Ln)
   integer::jx,jy,ii,jj,ix,iy
   integer::xy2i
   external xy2i

   mat_h = 0d0

   do ix = 1,Nx
      do iy = 1,Ny
         ii = xy2i(ix,iy,nx)
         jx = ix
         jy = iy
         jj = xy2i(jx,jy,nx)

         mat_h(ii,jj) = -mu

         jx = ix +1 
         if(jx > Nx) then
            jx = 1
         end if
         jy = iy
         jj = xy2i(jx,jy,nx)
         mat_h(ii,jj) = -1d0
         
         jx = ix - 1
         if(jx < 1) then
            jx = Nx
         end if
         jy = iy
         jj = xy2i(jx,jy,nx)
         mat_h(ii,jj) = -1d0

         jx = ix
         jy = iy + 1
         if(jy > Ny) then
            jy = 1
         end if
         jj = xy2i(jx,jy,nx)
         mat_h(ii,jj) = -1d0

         jx = ix
         jy = iy -1
         if(jy < 1) then
            jy = Nx
         end if
         jj = xy2i(jx,jy,nx)
         mat_h(ii,jj) = -1d0
         
      end do
   end do

   do ii = 1,Nx*Ny
      do jj = 1,Nx*Ny
         mat_h(ii+Nx*Ny,jj+Nx*Ny) = - mat_h(ii,jj)
      end do
      mat_h(ii,ii+Nx*Ny) = mat_delta(ii)
      mat_h(ii+Nx*Ny,ii) = mat_delta(ii)
 
   end do
!   write(*,*) mat_h


   return
 end subroutine calc_h


 program main

   implicit none
   interface
      subroutine calc_val(Ln,mat_h,val,row,col)
        integer,intent(in)::Ln
        real(8),intent(in)::mat_h(:,:)
        real(8),allocatable,intent(out)::val(:)
        integer,allocatable,intent(out)::row(:),col(:)
      end subroutine calc_val
   end interface

   integer::nx,ny,Ln
   real(8),allocatable::mat_delta(:)
   real(8),allocatable::mat_h(:,:)
   real(8)::mu
   real(8),allocatable::val(:)
   integer,allocatable::col(:),row(:)

   nx = 20
   ny = 20
   Ln = nx*ny*2
   allocate(mat_delta(nx*ny))

   call init_delta(nx,ny,mat_delta,0.1d0)   
   allocate(mat_h(Ln,Ln))
   mu = 1d-12

   call calc_h(nx,ny,mat_h,Ln,mu,mat_delta)
   call calc_val(Ln,mat_h,val,row,col)

!   write(*,*) val

 end program main

