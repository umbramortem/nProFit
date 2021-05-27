      program nprofit
      implicit none
      integer N,i,N_tmp,j,N_prof,N_arr,N_mirr,N_psf,im_size,i_tmp
      integer two_pow,op,N_mirr_psf,m,opt,N_rfit,i_down,i_up
      parameter (N=500001)
      character*5 W0_char,W0_char2
      character file_name*30,psf_name*30,obj_name*30,lib_path*100
      character key_model*8
      real W0,r_c,r_tmp,b(N),norm_psf,tmp,chi_sq_vec(N)
      real r_vec(N),W_vec(N),rho_vec(N),dens_surf_vec(N)
      real r_vec_down(N),W_vec_down(N),rho_vec_down(N),dens_surf_vec_down(N)
      real r_vec_up(N),W_vec_up(N),rho_vec_up(N),dens_surf_vec_up(N)
      real y2(N),r_prof(N),intens(N),int_err(N),n_pix(N)
      real r_arr(N),surf_vec_arr(N),surf_vec_mirr(N)
      real r_psf(N),psf(N),psf_zp(N),psf_mirr(N),r_psf2(N)
      real r_psf_mirr(N),r_mirr(N),conv(N),psf_step,y_conv(N)
      real tmp_mirr(N),psf_tmp_mirr(N),bg_value,rfit,texp
      real y_err(N),chi_sq,chi_best,rc_best,W0_best,model(N)
      real model_best(N),rc_step,W0_step,W0_down,W0_up,W0_tmp
      real rc_down,rc_up,W0_step2,W0_int,W0_int0
      integer N_arr_down,N_arr_up,N_tmp_down,N_tmp_up,N_arr_int
      real r_arr_down(N),surf_vec_down_arr(N),r_arr_up(N)
      real surf_vec_up_arr(N),surf_vec_int_arr(N),rc_tmp,rc_step2
      real rc_vec(N),W0_vec(N),chi_vec(N),chi_best2,rd,beta
      real beta_step,beta_step2,rd_step,rd_step2,rd_best,beta_best
      real beta_down,beta_up,rd_down,rd_up,beta_tmp,rd_tmp
      real rd_vec(N),beta_vec(N)
      character*8 formatW,formatW2
      complex ans(2*N)

      !For convolution using the FFT data should be stored in an array of length n
      !with n a power of two, the PSF in an array of length m<=n with m and odd number
      !and the output will be given in an array of length at least 2*n
      !The response function (psf) should be wrapped around, with the negative part
      !Copied in the positive end of the response function
 
      read(5,*) im_size
      read(5,*) opt
      read(5,*) psf_name
      read(5,*) file_name
      read(5,*) bg_value
      read(5,*) rfit
      read(5,*) obj_name
      read(5,'(a)') lib_path

      if (opt.ne.0) then
          W0_step=0.1
          W0_step2=0.01
          rc_step=0.1
          rc_step2=0.01
      else
          beta_step=0.1
          beta_step2=0.01
          rd_step=0.1
          rd_step2=0.01
      end if

      !Reading SBP profile and PSF data
      open(unit=8,file=file_name,status='old')
      j=1
11    read(8,*,end=98) r_prof(j),intens(j),int_err(j),n_pix(j)
      if (j.eq.1) int_err(j)=1
      y_err(j)=sqrt(int_err(j)+(bg_value**2)/n_pix(j))
      j=j+1
      goto 11      
98    continue      
      N_prof=j-1
      close(8) 

      N_rfit=0
      i=1
13    if ((r_prof(i)-rfit).le.0) then
           N_rfit=N_rfit+1
           i=i+1
           goto 13
      endif   

      open(unit=9,file=psf_name,status='old')
      j=1
12    read(9,*,end=97) r_psf(j),psf(j)
      j=j+1
      goto 12
97    continue
      N_psf=j-1
      close(9)

      !Computing the closest power of two to the data dimesion to create the data vector
      two_pow=1
40    if (im_size-2**two_pow.gt.0) then
          two_pow=two_pow+1
          goto 40
      endif

      !Creating psf radius vector with the dimensions of the required data vector
      r_psf2(1:N_psf)=r_psf(1:N_psf)
      psf_step=r_psf2(2)-r_psf2(1)
      i_tmp=1
41    r_psf2(N_psf+i_tmp)=r_psf(N_psf)+psf_step*i_tmp
      if ((i_tmp+N_psf).lt.(2**two_pow+1)) then
          i_tmp=i_tmp+1     
          goto 41
      endif

      psf_zp(:)=0.
      psf_zp(1:N_psf)=psf(1:N_psf)

      !Reading theoretical profiles up to the required size (power of two+1)
      if (opt.ne.0) then
          W0=2.
          r_c=0.5
      else
          rd=1.0
          beta=2.0
          call moff_x_vec(N_psf+i_tmp,r_psf2,N_mirr,r_mirr)
      endif
      chi_best=1e10
      chi_vec(:)=chi_best
      j=1
14    if (opt.eq.0) then
           call moffat(N_mirr,r_mirr,rd,beta,surf_vec_mirr)
           op=-1
           call mirr(N_psf+i_tmp,r_psf2,psf_zp,N_mirr_psf,r_psf_mirr,psf_mirr,op)
      end if

      if (opt.ne.0) then
      if (W0.ge.10) write(W0_char,'(f4.1)') W0
      if (W0.lt.10) write(W0_char,'(f3.1)') W0
          if (opt.eq.1) then
              open(unit=7,file=trim(adjustl(lib_path))//'king_W0_'//trim(adjustl(W0_char))//'.dat',status='old')
          endif
          if (opt.eq.2) then
              open(unit=7,file=trim(adjustl(lib_path))//'wilson_W0_'//trim(adjustl(W0_char))//'.dat',status='old')
          endif
      i=1
10    read(7,*,end=99) r_vec(i),W_vec(i),rho_vec(i),dens_surf_vec(i)
      if ((r_vec(i)*r_c).ge.(r_psf2(i_tmp+N_psf)+psf_step)) goto 99
      i=i+1
      goto 10
99    continue      
      close(7)

      N_tmp=i-1 
      call arrange(N_psf+i_tmp,r_psf2,N_tmp,r_vec*r_c,dens_surf_vec,N_arr,r_arr,surf_vec_arr)

      if (N_arr-1.ne.2**two_pow) then
          N_arr=2**two_pow+1
      endif
       
      op=1
      call mirr(N_arr,r_arr,surf_vec_arr,N_mirr,r_mirr,surf_vec_mirr,op)

      op=-1
      call mirr(N_arr,r_psf2,psf_zp,N_mirr_psf,r_psf_mirr,psf_mirr,op)
      end if

      call convlv(surf_vec_mirr/maxval(surf_vec_mirr(1:N_mirr)),N_mirr,psf_mirr/maxval(psf_mirr(1:N_mirr_psf)),N_mirr_psf,+1,ans) 

      m=1
42    y_conv(m)=sqrt(real(ans(m))**2+aimag(ans(m))**2)
      m=m+1
      if (m.le.(N_mirr/2)) goto 42
 
43    norm_psf=intens(1)/maxval(y_conv(1:N_mirr/2))

      do m=1,N_rfit
           tmp=y_conv(N_mirr/4+m-1)*norm_psf
           if (IsNaN(tmp).eqv..true.) tmp=0.
           chi_sq_vec(m)=(intens(m)-tmp)**2/((y_err(m)))**2
           model(m)=tmp
      enddo

      chi_sq=sum(chi_sq_vec(1:N_rfit))
      chi_vec(j)=chi_sq
      if (opt.eq.0) then
          rd_vec(j)=rd
          beta_vec(j)=beta
      else
          rc_vec(j)=r_c
          W0_vec(j)=W0
      endif
     
      if (chi_sq.lt.chi_best) then
           if (opt.eq.0) then
               rd_best=rd
               beta_best=beta
           else
               rc_best=r_c
               W0_best=W0
           endif
           chi_best=chi_sq
           model_best(1:N_mirr/4)=y_conv(N_mirr/4:N_mirr/2)*norm_psf
      endif
      j=j+1

      if (opt.eq.0) then
          rd=rd+rd_step
          if (rd.le.5.0) goto 14
          if ((rd.gt.5.0).and.(beta.lt.10)) then
              rd=0.5
              beta=beta+beta_step
              goto 14
          endif
          !write(*,*) chi_best,rd_best,beta_best
      else
          r_c=r_c+rc_step
          if (r_c.le.5.0) goto 14
          if ((r_c.gt.5.0).and.(W0.lt.15)) then
              r_c=0.5
              W0=W0+W0_step
              goto 14
          endif
      endif
 
      chi_best2=chi_best
      if (opt.ne.0) goto 56
      beta_down=beta_best-5*beta_step
      if (beta_down.lt.2.0) beta_down=2.0
      beta_up=beta_best+5*beta_step
      if (beta_up.lt.2.0) beta_up=2.1
      rd_down=rd_best-5*rd_step
      if (rd_down.lt.0.5) rd_down=0.5
      rd_up=rd_best+5*rd_step
      beta_tmp=beta_down 
      rd_tmp=rd_down
55    call moffat(N_mirr,r_mirr,rd_tmp,beta_tmp,surf_vec_mirr)
      go to 57
56    W0_down=W0_best-5*W0_step
      if (W0_down.lt.2) W0_down=2
      W0_up=W0_best+5*W0_step
      if (W0_up.gt.15) W0_up=15.
      rc_down=rc_best-5*rc_step
      if (rc_down.lt.0.5) rc_down=0.5
      rc_up=rc_best+5*rc_step
      W0_tmp=W0_down 
      rc_tmp=rc_down
     
17    formatW='(f4.1)'
      write(W0_char,formatW) W0_tmp
      formatW2='(f4.2)'
      read(W0_char,formatW2) W0_int
     
      !if (W0_int.le.W0_tmp) write(W0_char2,formatW) W0_tmp+0.1
      if ((W0_tmp+0.1).gt.15) then  
           if (W0_int.le.W0_tmp) write(W0_char2,formatW) 15.0
      else
           if (W0_int.le.W0_tmp) write(W0_char2,formatW) W0_tmp+0.1
      endif

      if ((W0_tmp-0.1).lt.2) then  
           if (W0_int.gt.W0_tmp) then
                W0_char2=W0_char
                write(W0_char,formatW) 2.0
           endif
      else
           if (W0_int.gt.W0_tmp) then
                W0_char2=W0_char
                write(W0_char,formatW) W0_tmp-0.1
           endif
      endif
 
      if (opt.eq.1) then
          open(unit=7,file=trim(adjustl(lib_path))//'king_W0_'//trim(adjustl(W0_char))//'.dat',status='old')
      end if
      if (opt.eq.2) then
          open(unit=7,file=trim(adjustl(lib_path))//'wilson_W0_'//trim(adjustl(W0_char))//'.dat',status='old')
      end if
      i_down=1
15    read(7,*,end=96) r_vec_down(i_down),W_vec_down(i_down),rho_vec_down(i_down),dens_surf_vec_down(i_down)
      if ((r_vec_down(i_down)*rc_tmp).ge.(r_psf2(i_tmp+N_psf)+psf_step)) goto 96
      i_down=i_down+1
      goto 15
      continue
96    close(7) 
      N_tmp_down=i_down-1
      if (opt.eq.1) then
          open(unit=8,file=trim(adjustl(lib_path))//'king_W0_'//trim(adjustl(W0_char2))//'.dat',status='old')
      endif
      if (opt.eq.2) then
          open(unit=8,file=trim(adjustl(lib_path))//'wilson_W0_'//trim(adjustl(W0_char2))//'.dat',status='old')
      endif
      i_up=1
16    read(8,*,end=95) r_vec_up(i_up),W_vec_up(i_up),rho_vec_up(i_up),dens_surf_vec_up(i_up)
      if ((r_vec_up(i_up)*rc_tmp).ge.(r_psf2(i_tmp+N_psf)+psf_step)) goto 95
      i_up=i_up+1
      goto 16
      continue
95    close(8) 
      N_tmp_up=i_up-1
      read(W0_char2,formatW) W0_int
      read(W0_char,formatW) W0_int0
18    call arrange(N_psf+i_tmp,r_psf2,N_tmp_down,r_vec_down*rc_tmp,dens_surf_vec_down,N_arr_down,r_arr_down,surf_vec_down_arr)
      call arrange(N_psf+i_tmp,r_psf2,N_tmp_up,r_vec_up*rc_tmp,dens_surf_vec_up,N_arr_up,r_arr_up,surf_vec_up_arr)
      if (N_arr_up.eq.N_arr_down) then
          surf_vec_int_arr(1:N_arr_up)=surf_vec_down_arr(1:N_arr_up)+(W0_tmp-W0_int0)*(surf_vec_up_arr(1:N_arr_up)-surf_vec_down_arr(1:N_arr_up))/(W0_int-W0_int0)
      endif
      if (N_arr_up.ne.N_arr_down) then
          N_arr_int=N_arr_down
          surf_vec_int_arr(1:N_arr_int)=surf_vec_down_arr(1:N_arr_int)+(W0_tmp-W0_int0)*(surf_vec_up_arr(1:N_arr_int)-surf_vec_down_arr(1:N_arr_int))/(W0_int-W0_int0)
          surf_vec_int_arr(N_arr_int+1:N_arr_up)=(W0_tmp-W0_int0)*(surf_vec_up_arr(N_arr_int+1:N_arr_up))/(W0_int-W0_int0)
      endif

      if (N_arr_up-1.ne.2**two_pow) then
          N_arr=2**two_pow+1
      endif

      op=1
      call mirr(N_arr,r_arr_up,surf_vec_int_arr,N_mirr,r_mirr,surf_vec_mirr,op)

57    call convlv(surf_vec_mirr/maxval(surf_vec_mirr(1:N_mirr)),N_mirr,psf_mirr/maxval(psf_mirr(1:N_mirr_psf)),N_mirr_psf,+1,ans) 

      m=1
44    y_conv(m)=sqrt(real(ans(m))**2+aimag(ans(m))**2)
      m=m+1
      if (m.le.(N_mirr/2)) goto 44
 
45    norm_psf=intens(1)/maxval(y_conv(1:N_mirr/2))

      do m=1,N_rfit
           tmp=y_conv(N_mirr/4+m-1)*norm_psf
           if (IsNaN(tmp).eqv..true.) tmp=0.
           chi_sq_vec(m)=(intens(m)-tmp)**2/((y_err(m)))**2
           model(m)=tmp
      enddo

      chi_sq=sum(chi_sq_vec(1:N_rfit))
      chi_vec(j)=chi_sq
      if (opt.eq.0) then
           rd_vec(j)=rd_tmp
           beta_vec(j)=beta_tmp
       else
           rc_vec(j)=rc_tmp
           W0_vec(j)=W0_tmp
       endif

      if (chi_sq.lt.chi_best) then
           if (opt.eq.0) then
               rd_best=rd_tmp
               beta_best=beta_tmp
           else
               rc_best=rc_tmp
               W0_best=W0_tmp
           endif
           chi_best=chi_sq
           model_best(1:N_mirr/4)=y_conv(N_mirr/4:N_mirr/2)*norm_psf
      endif


      if (opt.ne.0) goto 58

      rd_tmp=rd_tmp+rd_step2
      j=j+1
      if (rd_tmp.le.rd_up) goto 55
      if ((rd_tmp.gt.rd_up).and.(beta_tmp.le.beta_up)) then
          rd_tmp=rd_down
          beta_tmp=beta_tmp+beta_step2
          goto 55
      endif
      if (beta_tmp.gt.beta_up) goto 81

58    rc_tmp=rc_tmp+rc_step2
      j=j+1
      if (rc_tmp.le.rc_up) goto 18
      if ((rc_tmp.gt.rc_up).and.(W0_tmp.le.W0_up)) then
          rc_tmp=rc_down
          W0_tmp=W0_tmp+W0_step2
          goto 17
      endif

      if (W0_tmp.gt.W0_up) goto 81
      if (W0_tmp.gt.W0_int) goto 17 
     
81    continue
 
      if (opt.eq.0) key_model='moffat'
      if (opt.eq.1) key_model='king'
      if (opt.eq.2) key_model='wilson'
      open(unit=9,file=trim(adjustl(key_model))//'_model_object_'//trim(adjustl(obj_name))//'.dat')
      i=1
82    write(9,*) model_best(i)
      i=i+1 
      if (i.le.N_mirr/4) goto 82
      close(9)

      open(unit=8,file=trim(adjustl(key_model))//'_pars_object_'//trim(adjustl(obj_name))//'.dat')

      if (opt.eq.0) then
          write(8,*)'#chi rfit rd gamma rd_down rd_up gam_down gam_up'
          rd_down=minval(pack(rd_vec,chi_vec.le.(chi_best2+2.30)))
          rd_up=maxval(pack(rd_vec,chi_vec.le.(chi_best2+2.30)))
          beta_down=minval(pack(beta_vec,chi_vec.le.(chi_best2+2.30)))
          beta_up=maxval(pack(beta_vec,chi_vec.le.(chi_best+2.30)))
          write(8,*) chi_best2,rfit,rd_best,beta_best,rd_down,rd_up,beta_down,beta_up
      else
          write(8,*)'#chi rfit r_c W0 r_c_down r_c_up W0_down W0_up'
          rc_down=minval(pack(rc_vec,chi_vec.le.(chi_best2+2.30)))
          rc_up=maxval(pack(rc_vec,chi_vec.le.(chi_best2+2.30)))
          W0_down=minval(pack(W0_vec,chi_vec.le.(chi_best2+2.30)))
          W0_up=maxval(pack(W0_vec,chi_vec.le.(chi_best+2.30)))
          write(8,*) chi_best2,rfit,rc_best,W0_best,rc_down,rc_up,W0_down,W0_up
      endif

      close(8)

      end

      subroutine arrange(N_prof,r_prof,N_vec,r_vec,dens_surf_vec,N_arr,r_arr,surf_vec_arr)
      implicit none
      integer N_prof,N_vec,N_arr,s,i,j,k,k2,i_tmp,i_tmp2
      real r_prof(N_prof),r_vec(N_vec),dens_surf_vec(N_vec)
      real r_arr(N_vec),surf_vec_arr(N_vec),step,r_tmp,step2,dens_tmp
      real r_tmp2,r_it,dens_tmp2

      step=r_prof(3)-r_prof(2)
      j=1
      i=2
      k=1
      k2=1
      r_arr(:)=0
      surf_vec_arr(:)=0
      r_arr(1)=0
      surf_vec_arr(1)=dens_surf_vec(1)

      r_tmp=r_vec(1)
20    r_it=step*j
      if (r_it.gt.r_prof(N_prof)) goto 25
      if (i.ge.N_vec) then
         r_tmp=r_vec(N_vec)
         k=1
      end if
      if (r_tmp.le.r_it) then
          if (i.gt.N_vec) then 
             r_tmp=r_vec(N_vec)
          else 
             r_tmp=r_vec(i+1)
             i=i+1
             goto 20
          end if
      else
          r_arr(j+1)=r_it
          if (r_arr(j+1).gt.r_vec(N_vec)) goto 25
22        i_tmp=i-k
          r_tmp=r_vec(i_tmp)
          if (r_tmp.gt.r_it) then
               k=k+2
               goto 22  
          endif

23        i_tmp2=i_tmp+k2
          if (i_tmp2.gt.N_vec) i_tmp2=N_vec
          r_tmp2=r_vec(i_tmp2)
          if (r_tmp2.le.r_it) then
               k2=k2+1
               goto 23  
          endif
          i_tmp2=i_tmp2+1
          if (i_tmp2.gt.N_vec) i_tmp2=N_vec
          
       endif
       if (r_it.gt.r_vec(i_tmp2)) goto 25
       dens_tmp=dens_surf_vec(i_tmp)
       dens_tmp2=dens_surf_vec(i_tmp2)
       surf_vec_arr(j+1)=(dens_tmp*(r_vec(i_tmp2)-r_arr(j+1))+dens_tmp2*(r_arr(j+1)-r_vec(i_tmp)))/(r_vec(i_tmp2)-r_vec(i_tmp))
       j=j+1
       if (r_it.lt.r_vec(N_vec)) goto 20
25     N_arr=j
       return
       end

       subroutine mirr(N_arr,r_arr,surf_vec_arr,N_mirr,r_mirr,surf_vec_mirr,op)
       implicit none
       integer N_arr,N_mirr,i,op
       real r_arr(N_arr),surf_vec_arr(N_arr),r_mirr(N_mirr)
       real surf_vec_mirr(N_mirr)

       if (op.eq.1) N_mirr=2*N_arr-2
       if (op.eq.-1) N_mirr=2*N_arr-3

       r_mirr(:)=0
       surf_vec_mirr(:)=0
       i=1
       if (op.eq.1) then
27          r_mirr(i)=-r_arr(N_arr-i)
            surf_vec_mirr(i)=surf_vec_arr(N_arr-i)
            i=i+1
            if (i.lt.(N_arr)) goto 27 
            r_mirr(N_arr:N_mirr)=r_arr(2:N_arr)
            surf_vec_mirr(N_arr:N_mirr)=surf_vec_arr(2:N_arr)
       end if
       if (op.eq.-1) then
28          r_mirr(i)=-r_arr(N_arr-i)
            surf_vec_mirr(N_mirr-i+1)=surf_vec_arr(i)
            i=i+1
            if (i.lt.(N_arr-1)) goto 28
            r_mirr(N_arr-1:N_mirr)=r_arr(1:N_arr)
            surf_vec_mirr(1:N_arr)=surf_vec_arr(1:N_arr)
       endif
       end

       subroutine moff_x_vec(N_arr,r_prof,N_mirr,r_mirr)
       implicit none
       integer N_mirr,N_arr,i
       real step,r_prof(N_arr),r_mirr(N_mirr)

       step=r_prof(3)-r_prof(2)
       
       N_mirr=2*N_arr-2
35     r_mirr(i)=-r_prof(N_arr-i)
       i=i+1
       if (i.lt.(N_arr)) goto 35
       r_mirr(N_arr:N_mirr)=r_prof(2:N_arr)
       end

       subroutine moffat(N_prof,r_prof,rd,beta,prof)
       implicit none
       integer N_prof
       real r_prof(N_prof),rd,beta,prof(N_prof)

       prof(:)=0.
       prof(1:N_prof)=(1+(r_prof(1:N_prof)/rd)**2)**(-beta/2.)
       end
       
       subroutine convol(N_mirr,r_vec,surf_vec_mirr,r_psf_mirr,psf_mirr,conv)
       implicit none
       integer N_mirr,N_conv,i,j,i_tmp,j_tmp,g_tmp
       integer n_vec(N_mirr),m_vec(N_mirr),mid
       real r_vec(N_mirr),surf_vec_mirr(N_mirr),conv(N_mirr)
       real r_psf_mirr(N_mirr),psf_mirr(N_mirr),aux
       
       i=1
       n_vec(1)=-N_mirr/2
28     n_vec(i+1)=n_vec(i)+1
       i=i+1
       if (i.lt.N_mirr) goto 28
       m_vec=n_vec

       conv(:)=0
       i=N_mirr/2+1
       mid=i
29     aux=0.
       j=1      
30     i_tmp=n_vec(i)+i!+1
       j_tmp=m_vec(j)+i!+1 
       if ((j_tmp.lt.1).or.(j_tmp.gt.mid)) g_tmp=0
       if ((j_tmp.ge.1).or.(j_tmp.le.mid)) g_tmp=surf_vec_mirr(j_tmp)
       aux=aux+psf_mirr(i_tmp)*g_tmp
       j=j+1
       if (j.le.mid) goto 30
       if (j.gt.mid) then
           conv(i-(mid-1))=aux
           i=i+1
           if (i.gt.N_mirr) goto 31
           goto 29
       endif

31     conv=conv*(r_vec(N_mirr)-r_vec(1))/N_mirr
       end
