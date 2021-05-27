      program wilson_king_dynamical
      implicit none
      integer N,kount,i,Nmax,m,opt
      parameter (Nmax=10000)
      real*8 xp(Nmax),yp(Nmax,2),dens_vec(Nmax),den,deno,W0,pi
      real*8 dens(Nmax),W_vec(Nmax),r_vec(Nmax),rking,dens_surf(Nmax)
      real c_ind,den_rt,rt,step
      character*5 W0_char
 
      pi=4*atan(1.)

      opt=0
100   W0=2.
      step=0.1
      
      do i=1,131!11
           call generate_king(N,W0,rking,xp,yp,dens_vec,kount,den,c_ind,rt,den_rt,opt)
           call densty(W0,deno,opt)
           dens(1)=deno
           W_vec(1)=W0
           r_vec(1)=0
           dens(2:kount+1)=dens_vec(1:kount)
           W_vec(2:kount+1)=xp(1:kount)
           r_vec(2:kount+1)=yp(1:kount,1)
           call surf_dens(kount+1,r_vec(1:kount+1),dens(1:kount+1),dens_surf)
           if (W0.ge.10) write(W0_char,'(f4.1)') W0
           if (W0.lt.10) write(W0_char,'(f3.1)') W0
           if (opt.eq.0) open(unit=7,file='nprofit_library/king_W0_'//trim(W0_char)//'.dat')
           if (opt.eq.1) open(unit=7,file='nprofit_library/wilson_W0_'//trim(W0_char)//'.dat')
           do m=1,kount+1
                  write(7,*) sqrt(r_vec(m)),W_vec(m),dens(m),dens_surf(m)
           enddo
           W0=W0+step
      enddo
      if (opt.eq.1) goto 101
      if (opt.eq.0) opt=1
      goto 100

101   end

      subroutine densty(W,den,opt)
      implicit none
      integer opt
      real*8 W,den,pi
      pi=4*atan(1.0)
      if (opt.eq.0) den=-sqrt(W)*(W+1.5)+0.75*sqrt(pi)*exp(W)*erf(sqrt(W))
      if (opt.eq.1) den=-(2*W/3.+1.)*sqrt(4*W/pi)+exp(W)*erf(sqrt(W))-4*sqrt(4*W/pi)*W**2/15.
      return
      end
 
      subroutine surf(t1,t2,f0,f1,f2,fsurf)
      implicit none
      real t1,t2,f0,f1,f2,fsurf
      fsurf=4./15*sqrt(t1)*((5-t1/t2)*f0+(3*t1-5*t2)/(2*(t1-t2))*f1+(t1**2*f2)/(t2*(t1-t2)))
      return 
      end


      subroutine surf_dens(N,x,dens_vec,dens_surf)
      implicit none
      integer N,i_tmp,i
      real*8 x(N),dens_vec(N),x1,y1,y2(N),x2,dens_surf(N)
      real*8 y_int(N),rho_int(N),int_trapz,t_int(N)
      real t1,t2,f0,f1,f2,fsurf
     
      dens_surf(:)=0 
      do i=1,N-1
      x1=(x(i)+x(i+1))/2.
      x2=(x(i)+x1)/2.
      y_int(1)=x1
      i_tmp=N-i+1
      y_int(2:i_tmp)=x(i+1:N)
      t_int=(y_int(1:i_tmp)-x(i))/(1+y_int(1:i_tmp))
      rho_int(1)=(dens_vec(i+1)-dens_vec(i))/(x(i+1)-x(i))*(x1-x(i))+dens_vec(i)
      rho_int(2:i_tmp)=dens_vec(i+1:N)
      call trapz(i_tmp,y_int(1:i_tmp),rho_int(1:i_tmp),t_int(1:i_tmp),int_trapz,t_int(1),t_int(i_tmp))
      t2=(x1-x(i))/(1+x1)
      t1=(x2-x(i))/(1+x2)
      f0=(1+x(i))**1.5*dens_vec(i)
      f1=(1+x2)**1.5*((dens_vec(i+1)-dens_vec(i))/(x(i+1)-x(i))*(x2-x(i))+dens_vec(i))
      f2=(1+x1)**1.5*((dens_vec(i+1)-dens_vec(i))/(x(i+1)-x(i))*(x1-x(i))+dens_vec(i))
      call surf(t1,t2,f0,f1,f2,fsurf)
      dens_surf(i)=(int_trapz+fsurf)/(1+x(i))
      enddo
      return 
      end

      subroutine trapz(N,y,rho,t,int_trapz,a,b)
      implicit none
      integer N
      real*8 y(N),rho(N),t(N),int_trapz,a,b,delta
      real*8 f1(N-1),f2(N-1)

      f1=(1+y(1:N-1)**1.5)*rho(1:N-1)*t(1:N-1)**(-0.5)
      f2=(1+y(2:N)**1.5)*rho(2:N)*t(2:N)**(-0.5)
      delta=(b-a)/(N*1.)
      int_trapz=sum((f1+f2)/2.*delta)
      return
      end
      

      subroutine generate_king(N,W0,rking,xp,yp,dens_vec,kount,den,c_ind,rt,den_rt,opt)
      implicit none
      integer N,symmetry,M,Nmax,i,j,k,kount,kdens,opt
      parameter (M=10001,Nmax=10000)
      real*8 W0,star,rvir,rh,rking,D,pi
      real*8 h,den,xstart,ystart0,ystart1,x1,x2
      real*8 xp(Nmax),x(M),yking(M,2),mass(M),rmin,yp(Nmax,2)    
      real*8 pot,totmas,zh,ve,cg1,fmas,r2,costh,sinth,phi,r
      real*8 xstar,ystar,zstar,w1,w,wj1,wj,vmax,speed,fstar
      real*8 ustar,vstar,wstar,mstar,ri,coord(N,6)
      real*8 dens_vec(Nmax)
      real*8 sigma2,rc,rho_rt
      real c_ind,rt,den_rt
      character*4 conc

      rmin = 0.2
      pot=0.0
      pi=4*atan(1.0)
      kount=1

      if (W0.gt.15.1) then
              write(*,*) "W0 too large"
              goto 13
      endif
      if (W0.lt.2) then
              write(*,*) "W0 too small"
              goto 13
      endif

      ! * INTERPOLATE KING VALUES *

      h = (W0-rmin)**2/(1.0*M-1.0) 
      call densty(W0,den,opt)                

      xstart = 0.000001
      ystart0 = 2.0*xstart/3.0
      ystart1 = -2.0/3.0

      x1 = W0 - xstart
      x2 = 0.0

      call odeint(ystart0, ystart1, x1, x2, den, kount, xp, yp, M, Nmax,opt)

      kdens=1
19    if (xp(kdens).gt.0.0) kdens=kdens+1
      if (xp(kdens).gt.0.0) goto 19
      
      !//interpolate yking 
      
      do k=1,M
           x(k) = W0-rmin-sqrt(h*k)

           if (x(k).gt.xp(1)) then
                   yking(k,1) = (W0-x(k))*yp(1,1)/(W0 - xp(1))
                   yking(k,2) = yp(1,2)
           else 
                   i = 1         

                   do 14 while (i.le.kount)

                           if ((x(k)-xp(i))*(x(k)-xp(i+1)).le. 0.0) then 
                                   yking(k,1) = yp(i,1) + (yp(i+1,1)-yp(i,1))*(x(k)-xp(i))/(xp(i+1)-xp(i))
                                   yking(k,2) = yp(i,2) + (yp(i+1,2)-yp(i,2))*(x(k)-xp(i))/(xp(i+1)-xp(i))
                                   goto 15
                           else 
                                   i=i+1
                           endif
14                 continue
           endif

15         if (i.gt. kount) then
              yking(k,1) = yp(kount,1)
              yking(k,2) = yp(kount,2)
           endif
      enddo

      rt=yking(M-2,1)
      call densty(x(M-2),rho_rt,opt)
      den_rt=rho_rt
      
      rking = sqrt(yking(M-2,1))

      k=1
24    call densty(xp(k),dens_vec(k),opt)
      k=k+1

      if (k.le.kount) goto 24

13    return
      end

      subroutine rk4(x,y,dydx,h,yout,den,opt)
      implicit none
      integer i,opt
      real*8 x,y(2),dydx(2),h,yout(2),den
      real*8 hh,h6,xh
      real*8 yt(2),dyt(2),dym(2)

      hh=h*0.5
      h6=h/6.0
      xh=x+hh

      do i=1,2 
           yt(i) = y(i) + hh*dydx(i)
      enddo

      call derivs(xh,yt,dyt,den,opt)

      do i=1,2
           yt(i) = y(i) + hh*dyt(i)
      enddo

      call derivs(xh,yt,dym,den,opt)

      do i=1,2
           yt(i) = y(i) + h*dym(i)
           dym(i) = dym(i)+ dyt(i)
      enddo

      call derivs(x+h,yt,dyt,den,opt)
   
      do i=1,2
           yout(i) = y(i) + h6*(dydx(i)+dyt(i)+2.0*dym(i))
      enddo
      return
      end

      subroutine derivs(x,y,dydx,den,opt)
      integer opt
      real*8 x,y(2),dydx(2),den,rhox,pi

      pi=4*atan(1.0)

      rhox=0.0
      if ((x.ge.0.0).and.(opt.eq.0)) rhox =-sqrt(x)*(x+1.5)+0.75*sqrt(pi)*exp(x)*erf(sqrt(x))
      if ((x.ge.0.0).and.(opt.eq.1)) rhox =-(2*x/3.+1.)*sqrt(4*x/pi)+exp(x)*erf(sqrt(x))-4*sqrt(4*x/pi)*x**2/15.
      !else rhox = 0.0

      dydx(1)= y(2)
      dydx(2) = 0.25*y(2)**2*(6.0+9.0*y(2)*rhox/den)/y(1)

      return
      end

      subroutine rkqc(y,dydx,x,h,den,yscal,hdid,hnext,TOL,opt)
      implicit none 
      integer i,opt
      real*8 y(2),dydx(2),x,h,den,yscal(2),hdid,hnext,TOL
      real*8 safety,fcor,errcon,pgrow,pshrnk,xsav,ysav(2)
      real*8 dysav(2),ytemp(2),errmax,hh

      safety = 0.9
      fcor = 0.0666666667
      errcon = 6.0E-4
      pgrow = -0.20
      pshrnk = -0.25
      xsav = x

      do i=1,2
           ysav(i) = y(i)
           dysav(i) = dydx(i)
      enddo

10    hh = 0.5*h
      call rk4(xsav, ysav, dysav, hh, ytemp, den,opt)

      x = xsav + hh
      call derivs(x,ytemp,dydx,den,opt)
      call rk4(x, ytemp, dydx, hh, y, den,opt)

      x = xsav + h
      if (x.eq.xsav) then
              write(*,*) "ERROR: Stepsize not significant in RKQC."
              goto 11
      endif
      call rk4(xsav, ysav, dysav, h, ytemp, den,opt)

      errmax = 0.0

      do i=1,2
              ytemp(i) = y(i) - ytemp(i)
              errmax = max(errmax,sqrt((ytemp(i)/yscal(i))**2))
      enddo
      errmax = errmax/TOL

      if (errmax.gt.1.0) h = safety*h*(errmax**pshrnk) 
      if (errmax.gt.1.0) goto 10

      hdid = h
      hnext = 4.0*h
      if (errmax.gt.errcon) hnext = safety*h*(errmax**pgrow)

      do i=1,2
           y(i) = y(i)+ytemp(i)*fcor
      enddo

11    return
      end  

      subroutine odeint(ystart0,ystart1,x1,x2,den,kount,xp,yp,M,Nmax,opt)
      integer kount,Nmax,M,MAXSTP,i,j,opt
      parameter (MAXSTP=10000)
      real*8 ystart0,ystart1,x1,x2,den,xp(MAXSTP),yp(MAXSTP,2),y(2),dydx(2)
      real*8 HMIN,H1,tin,DXSAV,TOL,x,h,hdid,hnext,dydx2,y2,yscal2,yscal(2)

      HMIN = 0.0    
      H1 = 0.0001   
      tin = 1e-30
      DXSAV = 0.0001
      TOL = 1e-12     

      x = x1 
      h = - sqrt(H1**2)
      if ((x2-x1).ge.0.0) h = sqrt(H1**2)

      y(1) = ystart0  
      y(2) = ystart1  

      xsav = x-DXSAV*2.0

      do i=1,MAXSTP

           call derivs(x,y,dydx,den,opt) 

           do j=1,2
                 yscal(j)=sqrt(y(j)**2)+sqrt((h*dydx(j)))+tin
           enddo

           if (sqrt((x-xsav)**2).gt.sqrt(DXSAV**2)) then
                 if (kount.lt.(Nmax-1)) then
                     xp(kount) = x
                     do j=1,2
                          yp(kount,j) = y(j)
                     enddo
                     kount = kount + 1
                     xsav = x
                  endif
           endif
           !}  //store x, y1 and y2 if the difference in x is smaller as the desired output step size DXSAV

           if (((x+h-x2)*(x+h-x1)).gt.0.0) h = x2-x

           call rkqc(y,dydx,x,h,den,yscal,hdid,hnext,TOL,opt)

           if ((x-x2)*(x2-x1).ge.0.0) then 
                   ystart0 = y(1)
                   ystart1 = y(2)

                   xp(kount) = x
                    
                   do j=1,2
                        yp(kount,j) = y(j)
                   enddo
                   goto 12
                   kount = kount +1
           endif

           if (sqrt(hnext**2).lt.HMIN) then 
                   write(*,*) "Stepsize smaller than minimum."
                   goto 12
           endif

           h = hnext
      enddo
12    return
      end
