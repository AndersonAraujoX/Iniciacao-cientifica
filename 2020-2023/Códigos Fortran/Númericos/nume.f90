Program Real_Eigenvec
  implicit none

  integer, parameter:: n=25
  integer:: l,i,j,k, flag,m,ni,nk
  integer :: temp, bubble, lsup
  real*8::anali,gammab,r,dx,R_rel,pla
  real*8:: A(n,n),X(n,n) 
  real*8:: lambda(n),vec(n)
  real*8, parameter:: p=0.5d0/dsqrt(2.d0), q=1.d0-1.d0/dsqrt(2.d0)


  Print*,'Enter flag: 0 for DSYEV, 1 for DSYEVD'
  Read*, flag
  m=0
  !print*,'digite o valor de m'
  !read*,m
  dx=10.0/1000.0
  do i=0,50
  nk=0
  enddo
  !do l=0,0
    gammab=0.5d0!sqrt((4.0*abs(m)+3.0)*0.5)!variação de gamma dx*l*10.0!
    !construindo a matriz
    do i=1,n
      do j=1,n
        if(i==j) then
          A(i,j)=(j-0.5d0)+gammab*anali(i-1,j-1,m)
        else
          A(i,j)=gammab*anali(i-1,j-1,m)
        end if
      enddo
    enddo
    !print
    !do i=0,ni
    !  write(*,*)anali(0,i,0),i,gamma(0.5)
    !enddo
    !write(*,*)m
    !DO i  = 1, n
    !  PRINT 201, (A(i,j), j=1,n)
    !END DO

    !calculando os autovetores auto valores
    call Eigen(A,lambda,X,n,flag)
    !PRINT 200
    !print *, lambda(1),gammab
    !DO i  = 1, n
    PRINT *, (X(j,1), j=1,n)
    !END DO
    do k=0,1000
      R_rel=0.0d0
      r=k*dx
      !print*,'Dimension of the matrix, n=',int(sqrt(float(size(A))))
      do ni=1,n!calculando a função de onda
        R_rel=R_rel+X(1,1)/abs(X(1,1))*X(ni,1)*sqrt(gamma((ni)*1.0d0)*&
        gamma(nk+abs(m)+1.0d0)/(gamma(ni+abs(m)*1.0d0)*gamma(nk+1.0d0)))*&
        pla((ni-1),abs(m),r*r*0.5d0)
      enddo

      write(10+abs(m),*)r,R_rel

    enddo
    !enddo

  !enddo
100 FORMAT (1X,10(:2X,F10.2))
200 FORMAT (1X, 'Eigenvalue', 16X, 'Eigenvector^T')
201 FORMAT (1X, F10.7,4X,10(:f10.7))

  End Program Real_Eigenvec



!!!      SUBROUTINES -----------------------------------------
function DET(aa)
!real*8 function DET(aa)
real*8 aa(:,:)
real*8 tmp,c(size(aa,dim=1),size(aa,dim=2))
real*8 max
	integer i,j,k,l,m,num(size(aa,dim=1))
	n=size(aa,dim=1)
	det=1.	
	do k=1,n
		max=aa(k,k);num(k)=k;
		do i=k+1,n 
			if(abs(max)<abs(aa(i,k))) then
				max=aa(i,k)
				num(k)=i
			endif
		enddo
		if (num(k)/=k) then
			do l=k,n 
				tmp=aa(k,l)
				aa(k,l)=aa(num(k),l)
				aa(num(k),l)=tmp
			enddo
			det=-1.*det
		endif
		do m=k+1,n
			c(m,k)=aa(m,k)/aa(k,k)
			do l=k,n 
				aa(m,l)=aa(m,l)-c(m,k)*aa(k,l)
			enddo
		enddo !There we made matrix triangular!	
	enddo

	do i=1,n
	det=det*aa(i,i)
	enddo
	return
end function
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

    anali = anali + (-1.0)**i*(dgamma(i+abs(m)+0.5d0)*dgamma(n-i+0.5d0))/(dgamma(i+1.0d0)*dgamma(ni-i+1.0d0)*&
    dgamma(i+abs(m)+1.0d0)*dgamma(0.5d0-i)) 
  enddo
  anali=anali*sqrt(dgamma(ni+1.0d0)*dgamma(ni+abs(m)+1.0d0)/(dgamma(n+1.0d0)*2*gamma(n+abs(m)+1.0d0)))
end function

!função de laguerre, está funcionando testado


!function pla(n,m,x)
!  real*8, intent (in) :: x
!  integer , intent (in):: n,m
!  real*8  :: pla
!  pla=0.0
!  do i=0,n
!      pla = pla+(-x)**i*gamma(n+abs(m)+1.0)/(gamma(n-i+1.0)*gamma(i+abs(m)+1.0)*gamma(i+1.0)) 
!  enddo
!end function

RECURSIVE function pla(n,m,s) RESULT (lag)
  real*8:: s,lag
  integer:: n,m
  if(n>=2)then

    lag=((2*(n-1)+1+m-s)*pla(n-1,m,s)-(n-1+m)*pla(n-2,m,s))/n
  elseif(n==1)then

    lag=1+m-s
  else

    lag=1
  endif
end function
