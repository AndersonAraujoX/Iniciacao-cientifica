Program Real_Eigenvec
  implicit none

  integer, parameter:: n=5
  integer::i,j,k,flag,m,ni,nk,maxn,p
  integer :: temp, bubble, lsup
  real*8::anali,gammab,r,dx,R_rel,pla,DET
  REAL*8:: mu,tol,norm
  REAL*8,DIMENSION(n,n) :: A
  REAL*8,DIMENSION(n) :: x,xn,y
  real*8:: lambda(n),vec(n)
  !real*8, parameter:: p=0.5d0/dsqrt(2.d0), q=1.d0-1.d0/dsqrt(2.d0)


  Print*,'Enter flag: 0 for DSYEV, 1 for DSYEVD'
  Read*, flag
  !print*,'digite o valor de m'
  !read*,m
  dx=10.0/1000.0
  do m=0,-3,-1

    gammab=sqrt((4.0*abs(m)+3.0)*0.5)!variação de gamma
    
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
    !encontrar os autovalores

    WRITE(*,7)"A=",((A(i,j),j=1,n),i=1,n)

    maxn=100
    tol=0.0001

    WRITE(*,*)'step          x1              x2           x3            mu'
    WRITE(*,8)0,(x(i),i=1,n)

    DO p=1,n
        IF(ABS(x(p))==MAXVAL(ABS(x))) EXIT
    END DO

    DO k=1,maxn
        y=MATMUL(a,x)
        mu=y(p)

        DO p=1,n
            IF(ABS(y(p))==MAXVAL(ABS(y))) EXIT
        END DO

        IF(y(p)==0) THEN
            WRITE(*,*)'eigenvector',x
            WRITE(*,*)'A has the eigenvalue 0, select a new vector x and restart.'
            STOP
        END IF

        xn=y/y(p)
        norm=MAXVAL(ABS(x-xn))
        x=xn

        WRITE(*,8)k,(x(i),i=1,n),mu

        IF(norm<tol)THEN
            WRITE(*,9)"Dominant eigenvalue=",mu
            WRITE(*,9)"eigenvector=",(x(i),i=1,n)
            STOP
        END IF

    END DO
    WRITE(2,*)'Max number of iteration exceeded.'

    7 FORMAT(a,/,5(5(f12.7,3x),/))
    8 FORMAT(i4,6(3x,f12.7))
    9 FORMAT(a,/,5(f12.7,/))

    
    !calculando os autovetores auto valores
    !call Eigen(A,lambda,X,n,flag)
    !PRINT 200
    !DO i  = 1, n
    !  PRINT 201, lambda(i), (X(i,j), j=1,n)
    !END DO
    !organizando os autovalores 
    !print *, vec

    do nk=0,2


      do k=0,1000
        R_rel=0.0d0
        r=k*dx

        !print*,'Dimension of the matrix, n=',int(sqrt(float(size(A))))
        !do ni=1,n
        !  R_rel=R_rel+X(nk+1,ni)*sqrt(gamma((ni)*1.0d0)/gamma((ni)*1.0d0+abs(m)*1.0d0))*pla((ni-1),abs(m),r*r*0.5d0)
        !enddo

        write(10+nk*4+abs(m),*)r,R_rel
    
      enddo

    enddo

  enddo
!100 FORMAT (1X,10(:2X,F10.2))
!200 FORMAT (1X, 'Eigenvalue', 16X, 'Eigenvector^T')
!201 FORMAT (1X, F10.7,4X,10(:f10.2))

End Program Real_Eigenvec



!!!      SUBROUTINES -----------------------------------------

function DET(aa)
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
function pla(n,m,x)
  real*8, intent (in) :: x
  integer , intent (in):: n,m
  real*8  :: pla
  pla=0.0
  do i=0,n
      pla = pla+(-x)**i*gamma(n+abs(m)+1.0)/(gamma(n-i+1.0)*gamma(i+abs(m)+1.0)*gamma(i+1.0)) 
  enddo
end function