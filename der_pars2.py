import numpy as np
import scipy
import pylab as plt
import os
import scipy.special
import scipy.integrate
import scipy.interpolate

def ebin_func(r,gamma,rd,rho_o):
	G=0.0043
	tmp1=-(1/((-1+gamma)*(-1+gamma**2)))*(1+r**2/rd**2)**(1/2.*(-1-gamma))
	tmp_f=scipy.special.hyp2f1((1+gamma)/2.,(1+gamma)/2.,(3+gamma)/2.,-(rd**2/r**2))
	tmp2=tmp1*(-(1+gamma)*(r**2+rd**2*(1-(1+r**2/rd**2)**((1+gamma)/2)))-rd**2*(-1+gamma)*(1+rd**2/r**2)**((1+gamma)/2)*tmp_f) 
	rho_f=(1+(r/rd)**2)**(-(gamma+1)/2.)                                                
	return -tmp2*rho_f*r**2*4*np.pi**2*G*2*rho_o      

def Ir(r,a,b):
	tmp=(1./3)*r**3*scipy.special.hyp2f1(1.5,0.5*(1+b),2.5,-r**2/a**2)
	return tmp
	
def sig_rt(a,b,Io,rg,rho_o):
	inv= lambda x : Ir(x/a,a,b)/((x/a)**2*(1+(x/a)**2)**((b+1)/2.))
	tid_f=2*(100)**2/(rg*3.086e13)**2
	
	sigs=np.array([])
	sigs_cuad=np.array([])
	x_sigs_cuad=np.array([])
	flag=0
	step=0.1
	x=0
	while flag==0:
	    y, err = scipy.integrate.quad(inv,x*a, np.inf)
	    sig_cuad=4*np.pi*0.0043*a**2/a*(1+x**2/a**2)**((b+1)/2.)*np.abs(y)*rho_o-(a*3.086e13)**2*tid_f/(b-1)*(1+x**2/a**2)
	    sigs_cuad=np.append(sigs_cuad,sig_cuad)
	    x_sigs_cuad=np.append(x_sigs_cuad,x)
	    x=x+step
	    if sigs_cuad[-1]>=0:
	        flag=0
	    else:
	        flag=1
	int_rt=scipy.interpolate.interp1d(sigs_cuad,x_sigs_cuad)
	rt_mof=int_rt(0)
	return rt_mof,np.sqrt(sigs_cuad[0])

class derived_pars:

	def __init__(self,filters,data_dir,models,m_l_dic,Msun_dic):
		G=0.0043
		eb_ergs=1.989*1e43
		rg=1000.

		for k in filters:
			m_l=m_l_dic[k]
			for m in models:
				dir_pars=np.genfromtxt(os.path.join(data_dir,m+'_pars_'+k+'.dat'),comments='#',dtype='S')
				if np.size(dir_pars)<=10:
					ids=np.array([dir_pars[0]])
					chi=np.array([dir_pars[1]])
					npts=np.array([dir_pars[2]])
					mu=np.array([np.float(dir_pars[9])])
				else:
					ids=dir_pars[:,0]
					chi=dir_pars[:,1]
					npts=dir_pars[:,2]
					mu=dir_pars[:,9].astype(float)
				Io_vec=10**(0.4*(-mu+21.57+Msun_dic[k]))
				dIo_vec=np.sqrt(Io_vec)
				
				ofile=open('derived_parameters_'+m+'_'+k+'.dat','w')
				
				if 'moffat' in m:
					ofile.write('#ID\t c\t rc\t rc_up\t rc_down\t rh\t rh_up\t rh_down\t log_jo\t jo_up\t jo_down\t log_rho_o\t rho_o_up\t rho_o_down\t logLtot\t Ltot_up\t Ltot_down\t logMtot\t Mtot_up\t Mtot_down\t rt\t rt_up\t rt_down\t sigma\t sigma_up\t sigma_down\t logEb\t Eb_up\t Eb_down\n')
					for i in np.arange(0,len(ids),1):
						if len(ids)<=1:
							rd=float(dir_pars[3])
							rd_down=float(dir_pars[4])
							rd_up=float(dir_pars[5])
							gamma=float(dir_pars[6])
							gamma_down=float(dir_pars[7])
							gamma_up=float(dir_pars[8])
						else:
							rd=float(dir_pars[i,3])
							rd_down=float(dir_pars[i,4])
							rd_up=float(dir_pars[i,5])
							gamma=float(dir_pars[i,6])
							gamma_down=float(dir_pars[i,7])
							gamma_up=float(dir_pars[i,8])
						c=np.inf
						rc=rd*np.sqrt(2**(2./gamma)-1)
						drc_up=np.abs(rd_up*np.sqrt(2**(2./gamma)-1)-rd*0.5/np.sqrt(2**(2./gamma)-1)*2**(1+(2./gamma))*np.log(2)/gamma**2*gamma_up)
						drc_down=np.abs(rd_down*np.sqrt(2**(2./gamma)-1)+rd*0.5/np.sqrt(2**(2./gamma)-1)*2**(2./gamma)*np.log(2)-2./gamma**2*gamma_down)
						Io=Io_vec[i]
						dIo=dIo_vec[i]
						if gamma==2.:
							reff=np.inf
							dreff_up=np.inf
							dreff_down=np.inf
							Mtot=np.inf
							dMtot_up=np.inf
							dMtot_down=np.inf
							Ltot=np.inf
							dLtot_up=np.inf
							dLtot_down=np.inf
						else:
							reff=rd*np.sqrt(np.abs(0.5**(1./(1-gamma/2.))-1))
							dreff_up=rd_up*np.sqrt(np.abs(0.5**(1./(1-gamma/2.))-1))-0.25/np.sqrt(np.abs(0.5**(1./(1-gamma/2.))-1))*0.5**(1./(1-gamma/2.))*np.log(0.5)/(1./(1-gamma/2.)**2)*gamma_up*rd
							dreff_down=rd_down*np.sqrt(np.abs(0.5**(1./(1-gamma/2.))-1))-0.25/np.sqrt(np.abs(0.5**(1./(1-gamma/2.))-1))*0.5**(1./(1-gamma/2.))*np.log(0.5)/(1./(1-gamma/2.)**2)*gamma_down*rd
							Ltot=np.abs(2*np.pi*Io*rd**2/(gamma-2))
							dLtot_up=2*np.pi*((dIo*rd**2+Io*2*rd*rd_up)/(gamma-2)-Io*rd**2/(gamma-2)**2*gamma_up)
							dLtot_down=2*np.pi*((dIo*rd**2+Io*2*rd*rd_down)/(gamma-2)-Io*rd**2/(gamma-2)**2*gamma_down)
							Mtot=Ltot*m_l
							dMtot_up=dLtot_up*m_l
							dMtot_down=dLtot_down*m_l
						jo=Io*scipy.special.gamma(gamma/2.+0.5)/(np.sqrt(np.pi)*rd*scipy.special.gamma(gamma/2.))
						djo_up=np.abs((scipy.special.gamma(gamma/2.+0.5)/(np.sqrt(np.pi)*scipy.special.gamma(gamma/2.)))*(dIo/rd-Io/rd**2*rd_up))
						djo_down=np.abs((scipy.special.gamma(gamma/2.+0.5)/(np.sqrt(np.pi)*scipy.special.gamma(gamma/2.)))*(dIo/rd-Io/rd**2*rd_down))
					        rho_o=jo*m_l	
						drho_o_up=djo_up*m_l
						drho_o_down=djo_down*m_l
						x=0
						rt_mof,sigma=sig_rt(rd,gamma,Io,rg,rho_o)	
						drt_mof_up=np.sqrt(rd_up**2+gamma_up**2)
						drt_mof_down=np.sqrt(rd_down**2+gamma_down**2)
						dsigma_up=np.sqrt(drt_mof_up)
						dsigma_down=np.sqrt(drt_mof_down)
						inv= lambda x: ebin_func(x,rd,gamma,rho_o)
						eb_dic=scipy.integrate.quad(inv,0, rt_mof,full_output=1)
						eb=eb_dic[0]
						eb_err_up=np.sqrt((rd_up**2+gamma_up**2+drt_mof_up**2))*eb_ergs
						eb_err_down=np.sqrt((rd_down**2+gamma_down**2+drt_mof_down**2))*eb_ergs
						ebin=np.log10(np.abs(eb_ergs*eb))
						ofile.write(ids[i]+'\t'+str('%7.4f\t' %(c)))
						ofile.write(str('%7.4f\t' %(rc))+str('%7.4f\t' %(drc_up))+str('%7.4f\t' %(drc_down)))
						ofile.write(str('%10.4f\t' %(reff))+str('%10.4f\t' %(dreff_up))+str('%10.4f\t' %(dreff_down)))
						ofile.write(str('%10.4f\t' %(np.log10(jo)))+str('%7.4f\t' %(djo_up/jo))+str('%7.4f\t' %(djo_down/jo)))
						ofile.write(str('%10.4f\t' %(np.log10(rho_o)))+str('%7.4f\t' %(drho_o_up/rho_o))+str('%7.4f\t' %(drho_o_down/rho_o)))
						ofile.write(str('%10.4f\t' %(np.log10(Ltot)))+str('%7.4f\t' %(dLtot_up/Ltot))+str('%7.4f\t' %(dLtot_down/Ltot)))
						ofile.write(str('%10.4f\t' %(np.log10(Mtot)))+str('%7.4f\t' %(dMtot_up/Mtot))+str('%7.4f\t' %(dMtot_down/Mtot)))
						ofile.write(str('%10.4f\t' %(rt_mof))+str('%10.4f\t' %(drt_mof_up))+str('%10.4f\t' %(drt_mof_down)))
						ofile.write(str('%10.4f\t' %(sigma))+str('%10.4f\t' %(dsigma_up))+str('%10.4f\t' %(dsigma_down)))
						
						ofile.write(str('%7.4f\t' %(ebin))+str('%7.4f\t' %(np.abs(eb_err_up/(eb_ergs))))+str('%7.4f\t' %(np.abs(eb_err_down/(eb_ergs)))))
						ofile.write('\n')
						
						
				if ('king' in m) | ('wilson' in m):	
					ofile.write('#ID\t c\t c_up\t c_down\t rc\t rc_up\t rc_down\t rt\t rt_up\t rt_down\t rh\t rh_up\t rh_down\t log_jo\t jo_up\t jo_down\t log_rho_o\t rho_o_up\t rho_o_down\t logMtot\t Mtot_up\t Mtot_down\t logLtot\t Ltot_up\t Ltot_down\t sigma_o\t sigmao_up\t sigmao_down\t sigma\t sigma_up\t sigma_down\t logEb\t Eb_up\t Eb_down\n') 
					for i in np.arange(0,len(ids),1):
						if len(ids)<=1:
							ro=float(dir_pars[3])
							ro_down=float(dir_pars[4])
							ro_up=float(dir_pars[5])
							Wo=float(dir_pars[6])
							Wo_down=float(dir_pars[7])
							Wo_up=float(dir_pars[8])
						else:
							ro=float(dir_pars[i,3])
							ro_down=float(dir_pars[i,4])
							ro_up=float(dir_pars[i,5])
							Wo=float(dir_pars[i,6])
							Wo_down=float(dir_pars[i,7])
							Wo_up=float(dir_pars[i,8])
						Io=Io_vec[i]
						dIo=dIo_vec[i]
		
						if 'king' in m:
							c=3.622e-05*Wo**5-0.001362*Wo**4+0.01731*Wo**3-0.08129*Wo**2+0.3034*Wo+0.1208
							dc=5*3.622e-05*Wo**4-4*0.001362*Wo**3+3*0.01731*Wo**2-2*0.08129*Wo+0.3034
							log_reff=0.07105*c**5-0.7433*c**4+2.814*c**3-4.552*c**2+3.754*c-1.226
							dlog_reff=5*0.07105*c**4-4*0.7433*c**3+3*2.814*c**2-2*4.552*c+3.754
							rc_ro=-0.007654*c**6+0.104916*c**5-0.588633*c**4+1.732471*c**3-2.831729*c**2+2.457260*c-0.140352
							drc_ro=-6*0.007654*c**5+5*0.104916*c**4-4*0.588633*c**3+3*1.732471*c**2-2*2.831729*c+2.457260
							Io_joro=0.03026049*c**5-0.35282952*c**4+1.6248955*c**3-3.72541653*c**2+4.32147491*c-0.08966539
							dIo_joro=5*0.03026049*c**4-4*0.35282952*c**3+3*1.6248955*c**2-2*3.72541653*c+4.32147491
							int_mu_mass=0.01473853*c**5-0.28028575*c**4+1.63422783*c**3-3.91483438*c**2+4.74478698*c-2.24501706
							dint_mu_mass=5*0.01473853*c**4-4*0.28028575*c**3+3*1.63422783*c**2-2*3.91483438*c+4.74478698
		
							sigma_p2=6./5*(-np.sqrt(Wo)*(Wo**2+2.5*Wo+15./8)+15*np.sqrt(np.pi)/4.*scipy.special.erf(np.sqrt(Wo))*np.exp(Wo))/(-np.sqrt(Wo)*(Wo+1.5)+3./4*np.sqrt(np.pi)*np.exp(Wo)*scipy.special.erf(np.sqrt(Wo)))
							dsigma_p2_aux=-0.5/np.sqrt(Wo)*(Wo**2+2.5*Wo+15./8)-np.sqrt(Wo)*(2*Wo+2.5)+15/2.*np.exp(-Wo**2+Wo)+15*np.sqrt(np.pi)/4.*scipy.special.erf(np.sqrt(Wo))*np.exp(Wo)
							dsigma_p2_aux2=-0.5/np.sqrt(Wo)*(Wo+1.5)-np.sqrt(Wo)+3./4*np.sqrt(np.pi)*np.exp(Wo)*scipy.special.erf(np.sqrt(Wo))+1.5*np.exp(-Wo**2+Wo)
							dsigma_p2=6./5*dsigma_p2_aux/dsigma_p2_aux2
							eb_tmp=0.0033*c**5-0.2038*c**4+1.6194*c**3-4.7538*c**2+6.6155*c-2.9114
							deb_tmp=5*0.0033*c**4-4*0.2038*c**3+3*1.6194*c**2-2*4.7538*c+6.6155
						if 'wilson' in m:
							if Wo<8:
								c=-0.0006509*Wo**5+0.01645*Wo**4-0.1481*Wo**3+0.6292*Wo**2-1.086*Wo+1.387
								dc=-5*0.0006509*Wo**4+4*0.01645*Wo**3-3*0.1481*Wo**2+2*0.6292*Wo-1.086
							else:
								c=-0.001507*Wo**4+0.07088*Wo**3-1.211*Wo**2+9.014*Wo-21.51
								dc=-4*0.001507*Wo**3+3*0.07088*Wo**2-2*1.211*Wo+9.014
							if c<3.25:
								log_reff=0.05389*c**5-0.4788*c**4+1.708*c**3-3.007*c**2+2.79*c-1.253
								dlog_reff=5*0.05389*c**4-4*0.4788*c**3+3*1.708*c**2-2*3.007*c+2.79
							else:
								log_reff=5.03330062e+01*c**5-9.53958025e+02*c**4+7.22317915e+03*c**3-2.73123861e+04*c**2+5.15738755e+04*c-3.89063596e+04
								dlog_reff=5*5.03330062e+01*c**4-4*9.53958025e+02*c**3+3*7.22317915e+03*c**2-2*2.73123861e+04*c+5.15738755e+04
		
							if c<=3.3:
								rc_ro=0.014852*c**5-0.180891*c**4+0.878488*c**3-2.138649*c**2+2.641158*c-0.604565
								drc_ro=5*0.014852*c**4-4*0.180891*c**3+3*0.878488*c**2-2*2.138649*c+2.641158
		
							if c>3.3:
								rc_ro=0.488588*c**5-9.250956*c**4+69.973897*c**3-264.301886*c**2+498.517124*c-374.878421
								drc_ro=5*0.488588*c**4-4*9.250956*c**3+3*69.973897*c**2-2*264.301886*c+498.517124
		
							if c<=3.29:
								Io_joro=0.02551*c**5-0.3203*c**4+1.625*c**3-4.194*c**2+5.598*c-1.253
								dIo_joro=5*0.02551*c**4-4*0.3203*c**3+3*1.625*c**2-2*4.194*c+5.598
								int_mu_mass=0.06008356*c**5-0.64546529*c**4+2.80346846*c**3-6.21127072*c**2+7.40618805*c-3.81087141
								dint_mu_mass=5*0.06008356*c**4-4*0.64546529*c**3+3*2.80346846*c**2-2*6.21127072*c+7.40618805
							if c>3.29:
								Io_joro=-1.15235287e+01*c**6+2.61421371e+02*c**5-2.46839389e+03*c**4+1.24169313e+04*c**3-3.50964804e+04*c**2+5.28491926e+04*c-3.31209713e+04
								dIo_joro=-6*1.15235287e+01*c**5+5*2.61421371e+02*c**4-4*2.46839389e+03*c**3+3*1.24169313e+04*c**2-2*3.50964804e+04*c+5.28491926e+04
								int_mu_mass=-1.47993983e+02*c**6+3.36156870e+03*c**5-3.17813269e+04*c**4+1.60082454e+05*c**3-4.53086601e+05*c**2+6.83222747e+05*c-4.28823737e+05
								dint_mu_mass=-6*1.47993983e+02*c**5+5*3.36156870e+03*c**4-4*3.17813269e+04*c**3+3*1.60082454e+05*c**2-2*4.53086601e+05*c+6.83222747e+05
							sigma_p2=6./7*(-np.sqrt(Wo)*(Wo**3+3.5*Wo**2-35/4.*Wo+105./8)+105*np.sqrt(np.pi)/16.*np.exp(Wo)*scipy.special.erf(np.sqrt(Wo)))/(-np.sqrt(Wo)*(Wo**2+2.5*Wo+15./4)+15.*np.sqrt(np.pi)/8.*np.exp(Wo)*scipy.special.erf(np.sqrt(Wo)))
							dsigma_p2_aux=-0.5/np.sqrt(Wo)*(Wo**3+3.5*Wo**2-35/4.*Wo+105./8)-np.sqrt(Wo)*(3*Wo**2+7*Wo-35/4.)+105*np.sqrt(np.pi)/16.*np.exp(Wo)*scipy.special.erf(np.sqrt(Wo))+105./8*np.exp(-Wo**2+Wo)
							dsigma_p2_aux2=-0.5/np.sqrt(Wo)*(Wo**2+2.5*Wo+15./4)-np.sqrt(Wo)*(2*Wo+2.5)+15.*np.sqrt(np.pi)/8.*np.exp(Wo)*scipy.special.erf(np.sqrt(Wo))+15./4*np.exp(-Wo**2+Wo)
							dsigma_p2=6./7*dsigma_p2_aux/dsigma_p2_aux2
							
							if c<=3.3:
								eb_tmp=0.0784*c**5-0.8882*c**4+4.0641*c**3-9.4989*c**2+11.7447*c-5.6377
								deb_tmp=5*0.0784*c**4-4*0.8882*c**3+3*4.0641*c**2-2*9.4989*c+11.7447
							if c>3.3:
								eb_tmp=44.4913*c**5-843.3386*c**4+6386.5987*c**3-24154.2746*c**2+45623.7760*c-34431.0629
								deb_tmp=5*44.4913*c**4-4*843.3386*c**3+3*6386.5987*c**2-2*24154.2746*c+45623.7760
				
						dc_down=dc*Wo_down
						dc_up=dc*Wo_up
						if (dc_down==0) & (dc_up!=0):
							dc_down=dc_up
						if (dc_up==0) & (dc_down!=0):
							dc_up=dc_down
						rt=10**c*ro	
						drt_up=10**c*(np.log(10)*ro*dc_up+ro_up)/rt
						drt_down=10**c*(np.log(10)*ro*dc_down+ro_down)/rt
						reff=10**log_reff*ro
						dlog_reff_down=dlog_reff*dc_down
						dlog_reff_up=dlog_reff*dc_up
						dreff_up=(10**log_reff*np.log(10)*dlog_reff_up*ro+10**log_reff*ro_up)/reff
						dreff_down=(10**log_reff*np.log(10)*dlog_reff_down*ro+10**log_reff*ro_down)/reff
						rc=rc_ro*ro
						drc_up=drc_ro*ro*dc_up+rc_ro*ro_up
						drc_down=drc_ro*ro*dc_down+rc_ro*ro_down
						jo=Io/Io_joro/ro
						djo_up=(dIo*ro+Io*ro_up)/Io_joro-Io*ro/Io_joro**2*dIo_joro*dc_up
						djo_down=(dIo*ro+Io*ro_down)/Io_joro-Io*ro/Io_joro**2*dIo_joro*dc_down
						rho_o=jo*m_l
						drho_o_up=djo_up*m_l
						drho_o_down=djo_down*m_l
						Mtot=4*np.pi*rho_o*ro**3*10**int_mu_mass
						dint_mu_mass_up=dint_mu_mass*dc_up
						dint_mu_mass_down=dint_mu_mass*dc_down
						dMtot_up=4*np.pi*(drho_o_up*ro**3*10**int_mu_mass+rho_o*ro**2*ro_up*10**int_mu_mass+rho_o*ro**3*10**int_mu_mass*np.log(10)*dint_mu_mass_up)
						dMtot_down=4*np.pi*(drho_o_down*ro**3*10**int_mu_mass+rho_o*ro**2*ro_down*10**int_mu_mass+rho_o*ro**3*10**int_mu_mass*np.log(10)*dint_mu_mass_down)
						Ltot=Mtot/m_l
						dLtot_up=dMtot_up/m_l
						dLtot_down=dMtot_down/m_l
						sigma_o=np.sqrt(4*np.pi*G*rho_o*ro**2)/3.
						dsigma_o_up=2./3*np.pi*G/np.sqrt(4*np.pi*G*rho_o*ro**2)*(drho_o_up*ro**2+rho_o*2*ro*ro_up)
						dsigma_o_down=2./3*np.pi*G/np.sqrt(4*np.pi*G*rho_o*ro**2)*(drho_o_down*ro**2+rho_o*2*ro*ro_down)
						sigma=np.sqrt(sigma_p2)*sigma_o	
						dsigma_up=0.5/np.sqrt(sigma_p2)*dsigma_p2*Wo_up*sigma_o+np.sqrt(sigma_p2)*dsigma_o_up
						dsigma_down=0.5/np.sqrt(sigma_p2)*dsigma_p2*Wo_down*sigma_o+np.sqrt(sigma_p2)*dsigma_o_down
						ebin=2*np.pi*G*rho_o*ro**2*(Mtot/rt*int_mu_mass+sigma_o**2*eb_tmp)*eb_ergs
						debin_up=2*np.pi*G*eb_ergs*((drho_o_up*ro**2+rho_o*2*ro*ro_up)*(Mtot/rt*int_mu_mass+sigma_o**2*eb_tmp)+(rho_o*ro**2)*((dMtot_up/rt-Mtot/rt**2*drt_up)*int_mu_mass+Mtot/rt*dint_mu_mass_up+sigma_o*dsigma_o_up*eb_tmp+sigma_o**2*deb_tmp*dc_up))
						debin_down=2*np.pi*G*eb_ergs*((drho_o_down*ro**2+rho_o*2*ro*ro_down)*(Mtot/rt*int_mu_mass+sigma_o**2*eb_tmp)+(rho_o*ro**2)*((dMtot_down/rt-Mtot/rt**2*drt_down)*int_mu_mass+Mtot/rt*dint_mu_mass_down+sigma_o*dsigma_o_down*eb_tmp+sigma_o**2*deb_tmp*dc_down))
						ofile.write(ids[i]+'\t'+str('%7.4f\t' %(c))+str('%7.4f\t' %(dc_up))+str('%7.4f\t' %(dc_down)))
						ofile.write(str('%7.4f\t' %(rc))+str('%7.4f\t' %(drc_up))+str('%7.4f\t' %(drc_down)))
						ofile.write(str('%12.4f\t' %(rt))+str('%7.4f\t' %(drt_up))+str('%7.4f\t' %(drt_down)))
						ofile.write(str('%10.4f\t' %(reff))+str('%7.4f\t' %(dreff_up))+str('%7.4f\t' %(dreff_down)))
						ofile.write(str('%7.4f\t' %(np.log10(jo)))+str('%7.4f\t' %(np.abs(djo_up)/jo))+str('%7.4f\t' %(np.abs(djo_down)/jo)))
						ofile.write(str('%7.4f\t' %(np.log10(rho_o)))+str('%7.4f\t' %(np.abs(drho_o_up)/rho_o))+str('%7.4f\t' %(np.abs(drho_o_down)/rho_o)))
						ofile.write(str('%7.4f\t' %(np.log10(Mtot)))+str('%7.4f\t' %(np.abs(dMtot_up)/Mtot))+str('%7.4f\t' %(np.abs(dMtot_down)/Mtot)))
						ofile.write(str('%7.4f\t' %(np.log10(Ltot)))+str('%7.4f\t' %(np.abs(dLtot_up)/Ltot))+str('%7.4f\t' %(np.abs(dLtot_down)/Ltot)))
						ofile.write(str('%7.4f\t' %(sigma_o))+str('%7.4f\t' %(dsigma_o_up))+str('%7.4f\t' %(dsigma_o_down)))
						ofile.write(str('%7.4f\t' %(sigma))+str('%7.4f\t' %(dsigma_up))+str('%7.4f\t' %(dsigma_down)))
						ofile.write(str('%7.4f\t' %(np.log10(np.abs(ebin))))+str('%7.4f\t' %(np.abs(debin_up/ebin)))+str('%7.4f\t' %(np.abs(debin_down/ebin))))
						ofile.write('\n')
