Program Real_Eigenvec
  implicit none

  integer, parameter:: n=3
  integer:: i,j,k, flag,m,ni,nk
  real*8::anali,gammab,r,dx,R_rel,pla
  real*8:: A(n,n),X(n,n) 
  real*8:: lambda(n)
  real*8, parameter:: p=0.5d0/dsqrt(2.d0), q=1.d0-1.d0/dsqrt(2.d0)


  Print*,'Enter flag: 0 for DSYEV, 1 for DSYEVD'
  Read*, flag
  !print*,'digite o valor de m'
  !read*,m
  dx=10.0/1000.0
  do m=0,-3,-1

    gammab=sqrt((2.0*abs(m)+1)*0.5)!variação de gamma
    do i=1,n
      do j=1,n
        A(i,j)=gammab*anali(i-1,j-1,m)
      enddo
    enddo


    !calculando os autovetores auto valores
    call Eigen(A,lambda,X,n,flag)
    
    do nk=0,2


      do k=0,1000
        R_rel=0.0d0
        r=k*dx

        !print*,'Dimension of the matrix, n=',int(sqrt(float(size(A))))
        do ni=1,n
          R_rel=R_rel+X(nk+1,ni)*sqrt(gamma((ni)*1.0d0)/gamma((ni)*1.0d0+abs(m)*1.0d0))*pla((ni-1),abs(m),r*r*0.5d0)
        enddo

        write(10+nk*4+abs(m),*)r,R_rel
    
      enddo

    enddo

  enddo
100 FORMAT (1X,10(:2X,F10.2))
200 FORMAT (1X, 'Eigenvalue', 16X, 'Eigenvector^T')
201 FORMAT (1X, F10.2,4X,10(:f10.2))

  End Program Real_Eigenvec



!!!      SUBROUTINES -----------------------------------------


Subroutine Eigen(A,lambda,X,n,flag)
  implicit none

  integer:: i,n,flag
  real*8:: A(n,n),Ap(n,n),X(n,n)
  real*8:: lambda(n)
  real*8, allocatable  :: work(:)
      integer, allocatable :: iwork(:)
  integer:: lwork,liwork,info

  !print*,'n in Eiegen routine=',n

  lwork=3*n-1      ! DSYEV for flag=0
  if (flag==1) then ! DSYEVD for flag=1
    lwork=1+6*n+2*n**2
  end if 

  liwork=3+5*n


  allocate(work(lwork))
  allocate(iwork(liwork))

  Ap=A

  if (flag==0) then
    CALL DSYEV ('v', 'l', n, Ap, n, lambda, work, lwork, info)
  else
    CALL  DSYEVD ('V', 'U', n, Ap, n, lambda, work, &
              &   lwork, iwork, liwork, info)
    ! For doumentation visit: http://www.netlib.org/lapack/explore-html/d1/da2/dsyevd_8f.html
  end if

  X=Ap

  !print*,'info=',info

  deallocate(work)
  deallocate(iwork)

End Subroutine Eigen

!resultado analitico, bate com resultado do Python
!engraçado como eu falo comigo mesmo.
function anali(n,ni,m)
  integer , intent (in):: n,m,ni
  real*8  :: anali
  anali=0.0
  do i=0,ni!resultado analitco
    anali = anali + (-1)**i*(gamma(i+abs(m)+0.5)*gamma(n-i+0.5))/(gamma(i+1.0)*gamma(ni-i+1.0)*&
    gamma(i+abs(m)+1.0)*gamma(0.5-i)) 
  enddo
  anali=anali*sqrt(gamma(ni+1.0)*gamma(ni+abs(m)+1.0)/(gamma(n+1.0)*gamma(n+abs(m)+1.0)))
end function

!função de laguerre, está funcionando testado
function pla(n,m,x)
  real*8, intent (in) :: x
  integer , intent (in):: n,m
  real*8  :: pla
  pla=0.0
  do i=0,n
      pla = pla+(-x)**i*gamma(n+abs(m)+1.0)/(gamma(n-i+1.0)*gamma(i+abs(m)+1.0)*gamma(i+1.0)) 
  enddo
end function