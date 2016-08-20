!sub_chev.f90
subroutine calcsgap(i,val,row,col,valsize,N,aa,bb,nc,omegac,sgap)
  implicit none
  integer,intent(in)::i,N,valsize,nc
  real(8),intent(in)::aa,bb,omegac
  real(8),dimension(0:valsize-1),intent(in)::val
  integer,dimension(0:valsize-1),intent(in)::col
  integer,dimension(0:N),intent(in)::row
  real(8),intent(inout)::sgap

  integer::left_i,right_j
  real(8)::vec_ai(0:nc)
!  write(*,*) "fortran"

  left_i = i
  right_j = i + N/2

  sgap = 0d0
  call calc_polynomials_f(nc,N,valsize,left_i,right_j,val,row,col,vec_ai)
  call calc_meanfield(vec_ai,aa,bb,omegac,nc,sgap)

  return

  contains

    subroutine calc_meanfield(vec_ai,aa,bb,omegac,nc,sgap)
      implicit none
      integer,intent(in)::nc
      real(8),intent(in)::vec_ai(0:nc)
      real(8),intent(in)::aa,bb,omegac
      real(8),intent(out)::sgap
      real(8)::ba,omeb,pi     
      integer::j,i
      pi = atan(1d0)*4d0

      ba = acos(-bb/aa)
      omeb = acos(-(omegac+bb)/aa)
      sgap = 0d0
      do j = 1,nc
         sgap = sgap + vec_ai(j)*(sin(dble(j)*omeb)-sin(dble(j)*ba))/dble(j)
      end do
      sgap = sgap + vec_ai(0)*(omeb-ba)/2d0
      sgap = sgap*2d0/pi

      return
    end subroutine calc_meanfield
    
    subroutine calc_polynomials_f(nc,N,valsize,left_i,right_j,val,row,col,vec_ai)
      implicit none
      integer,intent(in)::nc,left_i,right_j,N,valsize
      real(8),intent(in)::val(0:valsize-1)
      integer,intent(in)::row(0:N)
      integer,intent(in)::col(0:valsize-1)
      real(8),intent(out)::vec_ai(0:nc)
      real(8),dimension(0:N-1)::vec_jn,vec_jnm,vec_jnmm,vec_temp
      integer::nn
      vec_jnmm = 0d0
      vec_jnm = 0d0
      vec_jn = 0d0
      vec_jn(right_j) = 1d0
      vec_ai = 0d0

!      write(*,*) "polynomial"

      do nn = 0,nc
         if(nn == 0) then
            vec_jnmm = 0d0
            vec_jnm = 0d0
            vec_jn(right_j) = 1d0
         else if(nn == 1) then
            call matvec(vec_jn,val,row,col,N,valsize,vec_temp)
            vec_jn = vec_temp 
         else 
            call matvec(vec_jnm,val,row,col,N,valsize,vec_temp)
            vec_jn = 2d0*vec_temp - vec_jnmm
         end if
         vec_ai(nn) = vec_jn(left_i)
         vec_jnmm = vec_jnm
         vec_jnm = vec_jn
      end do
      
      return
    end subroutine calc_polynomials_f
   

    subroutine matvec(vec_a,val,row,col,N,valsize,vec_temp)
      implicit none
      integer,intent(in)::N,valsize
      real(8),intent(in)::val(0:valsize-1)
      integer,intent(in)::col(0:valsize-1)
      integer,intent(in)::row(0:N)
      real(8),intent(in)::vec_a(0:N-1)
      real(8),intent(out)::vec_temp(0:N-1)
      integer::i,jj,j

      vec_temp = 0d0
!      write(*,*) "matvec",valsize
!      do i = 0,N-1
!         write(*,*) row(i),vec_a(i)
!      end do
!      write(*,*) "r",row(N)
!      write(*,*) "v",val
!      write(*,*) "r",row
!      write(*,*) "c",col

      do i = 0,N-1
!         write(*,*) "i",i
         do jj = row(i),row(i+1)-1
!            write(*,*) jj
            j = col(jj)
            vec_temp(i) = vec_temp(i) + val(jj)*vec_a(j)
!            write(*,*) i,vec_temp(i)
         end do
      end do
!      write(*,*) "end"
      

      
      return
    end subroutine matvec
 


  end subroutine calcsgap
