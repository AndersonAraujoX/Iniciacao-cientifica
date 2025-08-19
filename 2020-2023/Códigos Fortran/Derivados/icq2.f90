program IC
    implicit real*8(a-h,o-z)
    do m=0,5
        dx=10.0/1000
        gab=sqrt((4*abs(m)+3.0))
        snorm=0
        a0=sqrt(2**(abs(m)+0.5)*gamma(abs(m)*1.0+0.5))
        a1=2*gab*a0/(2*abs(m)+1)
        a2=(2*gab*a1+(2.0-3.0-2.0)*a0)/(2*abs(m)+2.0)*2.0

        do i=0,1000
            r=i*dx
            p=a2*r*r+a1*r+a0
            snorm=snorm+(r**abs(m)*exp(-r*r*0.25)*p)**2*dx
        enddo

        do i=0,1000
            r=i*dx
            p=a2*r*r+a1*r+a0
            write(10+m,*)r,(r**abs(m)*exp(-r*r*0.25)*p)/snorm
        enddo
    
    enddo
    end program
    