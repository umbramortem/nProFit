import numpy as np
import scipy
import scipy.stats
import scipy.interpolate
import pylab as plt
import os
import logging
from pyraf import iraf
import warnings

class pixel_maps:

	def __init__(self,filters,mask_file,names,fit_box,box_opt,data_dir,suff_mask):
		self.id_mask=scipy.genfromtxt(mask_file,comments='#',usecols=0,dtype='S')
		self.n_mask=scipy.genfromtxt(mask_file,comments='#',usecols=2).astype(int)
		self.mask_file=mask_file
		self.filters=filters
		self.names=names
		self.box_opt=box_opt
                self.data_dir=data_dir
                if self.box_opt==1:
                        self.box_tmp=str(int(fit_box)+1)
                if self.box_opt==2:
                        self.box_tmp_dic=fit_box
		self.suff_mask=suff_mask

	def mask(self):

		for k in self.filters:
			for i in np.arange(0,len(self.id_mask),1):
				iraf.fields.lines=str(i+2)

                		if self.n_mask[i]>0:

                		        temp=self.n_mask[i]
                		        j=0

                		        while temp!=0:
						if self.box_opt==1:
                                        		self.fit_box=self.box_tmp
                                		if self.box_opt==2:
                                        		self.fit_box=str(int(self.box_tmp_dic[k][:,1][self.box_tmp_dic[k][:,0]==self.id_mask[i]][0])+1)
                		                x_mask=iraf.fields(self.mask_file,str(int(4+3*j)),Stdout=1)
                		                y_mask=iraf.fields(self.mask_file,str(int(5+3*j)),Stdout=1)
                		                z_mask=iraf.fields(self.mask_file,str(int(6+3*j)),Stdout=1)
                		                r_mask=float(z_mask[0])/2.0
                		                expr='(sqrt((I-'+str(x_mask[0])+')**2+(J-'+str(y_mask[0])+')**2)<='+str(r_mask)+')'

						if self.n_mask[i]==1:
                		                        iraf.imexpr(expr,os.path.join(self.data_dir,self.names[k]+'_'+str(self.id_mask[i])+self.suff_mask+'.pl'),dims=self.fit_box+','+self.fit_box)

                		                else:
                		                        iraf.imexpr(expr,'temp'+str(self.id_mask[i])+'.'+str(j+1)+'.pl',dims=self.fit_box+','+self.fit_box)

                		                        if temp==self.n_mask[i]:
                		                                iraf.imcopy.verbose='no'
                		                                iraf.imcopy('temp'+str(self.id_mask[i])+'.'+str(j+1)+'.pl','salida_temp_1.pl')
                		                                iraf.imdelete('temp'+str(self.id_mask[i])+'.'+str(j+1)+'.pl')

                		                        if temp<self.n_mask[i]:
                		                                s1='temp'+str(self.id_mask[i])+'.'+str(j+1)+'.pl'
                		                                s2="salida_temp_1.pl"
                		                                iraf.imexpr.dims=self.fit_box+','+self.fit_box
								iraf.imexpr('a+b','salida_temp_2.pl','temp'+str(self.id_mask[i])+'.'+str(j+1)+'.pl','salida_temp_1.pl')
                		                                iraf.imdelete(s1)
                		                                iraf.imdelete(s2)
                		                                iraf.imcopy('salida_temp_2.pl',s2)
                		                                iraf.imdelete('salida_temp_2.pl')

                		                                if temp-1==0:
                		                                        iraf.imcopy(s2,os.path.join(self.data_dir,self.names[k]+'_'+str(self.id_mask[i])+self.suff_mask+'.pl'))

                		                                        iraf.imdelete(s2)

                		                j=j+1
                		                temp=temp-1

class ellipse:

	def __init__(self,filters,opt_ellip,names,id_list,mag_zero,fit_box,ellip_f,instat_dic,coord_file,suff,data_dir):
		self.filters=filters
		self.opt_ellip=opt_ellip
		self.names=names
                self.suff=suff
                self.suff2=suff[:-1]
		self.id_list=id_list
		self.mag_zero=mag_zero
		self.fit_box=fit_box
                self.instat_dic=instat_dic
                self.data_dir=data_dir
		iraf.imcopy.verbose='no'
		
		if (self.opt_ellip==1):
                	self.id_coor=coord_file[:,0]
                	self.x0=coord_file[:,1].astype(int)
                	self.y0=coord_file[:,2].astype(int)
			self.ellip_f=ellip_f

		if (self.opt_ellip==3) or (self.opt_ellip==4):
			self.ellipse_data=coord_file
                        self.list_file=open(os.path.join(self.data_dir,'intensity_prof_info.dat'),'w')
                        self.list_file.write('#Filter\t list_file\n')

		if (self.opt_ellip==5):
			self.ellipse_data={}
           		self.ellipse_data=ellip_f 
                        self.list_file=open(os.path.join(self.data_dir,'psf_profiles_info.dat'),'w')
                        self.list_file.write('#Filter\t file\n')

	def task_ellip(self):
		ellip_dic={}
		rfit_dic={}
                error_log=open('nprofit.log','a')
                logging.basicConfig(filename = 'nprofit.log', level = logging.DEBUG)
		logger=logging.getLogger(__name__)
		for k in self.filters:

			if self.opt_ellip==1:
				id_ellip=np.array([])
				ellip_array=np.array([])
                                pa_array=np.array([])
				ellip_file=open('ellip_'+k+'.dat','w')	
                                ellip_file.write('#ID\t ELLIP\t P.A.\n')

			if ((self.opt_ellip==3) or (self.opt_ellip==4)) or (self.opt_ellip==5):
				sma_array=np.array([])
				intens_array=np.array([])
				int_err_array=np.array([])
				npix_array=np.array([])

			if (self.opt_ellip==3) or (self.opt_ellip==4):
				list_intens_prof=open(os.path.join(self.data_dir,'list_intens_profiles_'+k+'.dat'),'w')
                                self.list_file.write(k+'\t list_intens_profiles_'+k+'.dat\n')


			for i in np.arange(0,len(self.id_list),1): 
                		iraf.stsdas(_doprint=0)
                		iraf.analysis(_doprint=0)
                		iraf.isophote(_doprint=0)
		
				if (self.opt_ellip==1):
					if np.size(self.instat_dic[k])<=3:
						bg_value=self.instat_dic[k][1][0]
                                        	rms_value=self.instat_dic[k][2][0]
					else:
                                        	bg_value=self.instat_dic[k][:,1][self.instat_dic[k][:,0]==self.id_list[i]]
                                        	rms_value=self.instat_dic[k][:,2][self.instat_dic[k][:,0]==self.id_list[i]]

					try:
						iraf.imcopy(os.path.join(self.data_dir,self.names[k]+'_'+str(self.id_list[i])+self.suff),os.path.join(self.data_dir,self.names[k]+'_'+str(self.id_list[i])+self.suff2))
						h_flag=0
					except iraf.IrafError,iraf_err:
						iraf.imcopy(os.path.join(self.data_dir,self.names[k]+'_'+str(self.id_list[i])+self.suff+'[1]'),os.path.join(self.data_dir,self.names[k]+'_'+str(self.id_list[i])+self.suff2))
						h_flag=1
                			e_0=0.05
                			pa_0=5.0
                			hcenter='yes'
                			hellip='no'
                			hpa='no'
                                        recente='no'
                                        minsma=1.
					suff_ellip='.fit'

				if (self.opt_ellip==3) or (self.opt_ellip==4):
					if os.path.exists(os.path.join(self.data_dir,self.names[k]+'_'+str(self.id_list[i])+self.suff2))==False:
						try:
							h_flag=0
							iraf.imcopy(os.path.join(self.data_dir,self.names[k]+'_'+str(self.id_list[i])+self.suff),os.path.join(self.data_dir,self.names[k]+'_'+str(self.id_list[i])+self.suff2))
						except iraf.IrafError,iraf_err:
							iraf.imcopy(os.path.join(self.data_dir,self.names[k]+'_'+str(self.id_list[i])+self.suff+'[1]'),os.path.join(self.data_dir,self.names[k]+'_'+str(self.id_list[i])+self.suff2))
							h_flag=1
                			hcenter='yes'
                			hellip='yes'
                			hpa='yes'
                                        minsma=0.
				        linear='yes'
					suff_ellip=self.suff2

				if (self.opt_ellip==5):
					h_flag=0
					if os.path.exists(os.path.join(self.data_dir,self.names[k]+self.suff))==False:
						try:
							iraf.imcopy(self.names[k]+self.suff,os.path.join(self.data_dir,self.names[k]+self.suff))
							h_flag=0
						except iraf.IrafError,iraf_err:
							iraf.imcopy(self.names[k]+self.suff+'[1]',os.path.join(self.data_dir,self.names[k]+self.suff))
							h_flag=1
					if (h_flag==1):
                                		hselect_out=iraf.hselect(os.path.join(self.data_dir,self.names[k]+self.suff+'[1]'),"NAXIS1 NAXIS2",'yes',Stdout=2)
					else:
                                		hselect_out=iraf.hselect(os.path.join(self.data_dir,self.names[k]+self.suff),"NAXIS1 NAXIS2",'yes',Stdout=2)

				if (self.opt_ellip==3) or (self.opt_ellip==4):
					h_flag=0
					mask_eo=self.ellipse_data[k][:,0]==self.id_list[i]
					e_0=float(self.ellipse_data[k][:,3][mask_eo][0])
					pa_0=float(self.ellipse_data[k][:,4][mask_eo][0])
			                if pa_0<-90.:
                                                pa_0=-90.
                                        if pa_0>90.:
                                                pa_0=90.

				if (self.opt_ellip==1) or ((self.opt_ellip==3) or (self.opt_ellip==4)):
					if (h_flag==1):
                                		hselect_out=iraf.hselect(os.path.join(self.data_dir,self.names[k]+'_'+str(self.id_list[i])+self.suff+'[1]'),"NAXIS1 NAXIS2 EXPTIME",'yes',Stdout=2)
					else:
                                		hselect_out=iraf.hselect(os.path.join(self.data_dir,self.names[k]+'_'+str(self.id_list[i])+self.suff),"NAXIS1 NAXIS2 EXPTIME",'yes',Stdout=2)
						if hselect_out==[]:
                                			hselect_out=iraf.hselect(os.path.join(self.data_dir,self.names[k]+'_'+str(self.id_list[i])+self.suff+'[0]'),"NAXIS1 NAXIS2 EXPTIME",'yes',Stdout=2)
				
				try:		
                                	self.texp=str(hselect_out[0]).split()[2]
				except:
					self.texp='1'
	

                                x0=str(int(str(hselect_out[0]).split()[0])/2+1)
                                y0=str(int(str(hselect_out[0]).split()[1])/2+1)

                                if (self.opt_ellip==1):
					maxsma=x0
					linear='no'
					step=0.1

				if ((self.opt_ellip==3) or (self.opt_ellip==4)):
					maxsma=x0
					linear='yes'
					step=1.0
					minsma=0.

                		mag0=self.mag_zero[k][0]
				olthresh=0
				sma0=5.
                		verbose='no'
				flag_err=False
				len_arr=-99

				try:
					if (self.opt_ellip==1) or ((self.opt_ellip==3) or (self.opt_ellip==4)):
						mask_name=os.path.join(self.data_dir,self.names[k]+'_'+str(self.id_list[i])+self.suff2).split('.fit')[0]+'.pl'
						if (os.path.exists(os.path.join(self.data_dir,self.names[k]+'_'+str(self.id_list[i])+'.pl'))==True) and (os.path.exists(mask_name)==False):
							iraf.imcopy(os.path.join(self.data_dir,self.names[k]+'_'+str(self.id_list[i])+'.pl'),mask_name)
                				iraf.ellipse(os.path.join(self.data_dir,self.names[k]+'_'+str(self.id_list[i])+self.suff2),os.path.join(self.data_dir,self.names[k]+'_'+str(self.id_list[i])+self.suff[:-4]+'tab'),x0=x0,y0=y0,ellip0=e_0,pa0=pa_0,sma0=sma0,minsma=minsma,maxsma=maxsma,step=step,linear=linear,hcenter=hcenter,hellip=hellip,hpa=hpa,olthresh=olthresh,mag0=mag0,verbose=verbose)
					if self.opt_ellip==5:
						e_0=0.05
						pa_0=20.
						sma0=5.
						minsma=0.
						maxsma=x0
						step=0.5
						linear='yes'
						hcenter='yes'
						hellip='yes'
						hpa='yes'
						mag0=self.mag_zero[k][0]
						verbose='no'
						iraf.ellipse(os.path.join(self.data_dir,self.names[k]+self.suff),os.path.join(self.data_dir,self.names[k]+'.tab'),x0=x0,y0=y0,ellip0=e_0,pa0=pa_0,sma0=sma0,minsma=minsma,maxsma=maxsma,step=step,linear=linear,hcenter=hcenter,hellip=hellip,hpa=hpa,olthresh=olthresh,mag0=mag0,verbose=verbose)
				except iraf.IrafError, iraf_err:

                                        error_str=str(self.id_list[i])+' filter '+k
                                        logger.error('error was caught in object %s',error_str)
                                        logger.error(iraf_err)
					flag_err=True

				if flag_err==False:
					if (self.opt_ellip==1) or ((self.opt_ellip==3) or (self.opt_ellip==4)):
						len_arr=int(iraf.tinfo(os.path.join(self.data_dir,self.names[k]+'_'+str(self.id_list[i])+self.suff[:-4]+'tab'),Stdout=1)[1].split('rows written to table')[0])
					
					if self.opt_ellip==5:
						len_arr=int(iraf.tinfo(os.path.join(self.data_dir,self.names[k]+'.tab'),Stdout=1)[1].split('rows written to table')[0])

				if (flag_err==True) | ((len_arr<=1) & (len_arr!=-99)):
					if self.opt_ellip==1:
                                        	iraf.delete(os.path.join(self.data_dir,self.names[k]+'_'+str(self.id_list[i])+self.suff[:-4]+'tab'))
                                        	iraf.delete(os.path.join(self.data_dir,self.names[k]+'_'+str(self.id_list[i])+self.suff2))
						suff_ellip=self.suff
						try:
                					iraf.ellipse(os.path.join(self.data_dir,self.names[k]+'_'+str(self.id_list[i])+suff_ellip),os.path.join(self.data_dir,self.names[k]+'_'+str(self.id_list[i])+self.suff[:-4]+'tab'),x0=x0,y0=y0,ellip0=e_0,pa0=pa_0,sma0=sma0,minsma=minsma,maxsma=maxsma,step=step,linear=linear,hcenter=hcenter,hellip=hellip,hpa=hpa,olthresh=olthresh,mag0=mag0)
						except iraf.IrafError, iraf_err:
                					iraf.ellipse(os.path.join(self.data_dir,self.names[k]+'_'+str(self.id_list[i])+suff_ellip+'[1]'),os.path.join(self.data_dir,self.names[k]+'_'+str(self.id_list[i])+self.suff[:-4]+'tab'))
						len_arr=int(iraf.tinfo(os.path.join(self.data_dir,self.names[k]+'_'+str(self.id_list[i])+self.suff[:-4]+'tab'),Stdout=1)[1].split('rows written to table')[0])
                        
				if (self.opt_ellip==1):
					tab_name=os.path.join(self.data_dir,self.names[k]+'_'+str(self.id_list[i])+self.suff[:-4]+'tab')
					if (os.path.exists(tab_name)==True):
                                        	iraf.tprint.rows='1-'+str(len_arr)
						sky_thres=float(rms_value[0])
                                        	iraf.tprint.showro='no'
                                        	iraf.tprint.showhdr='no'
                                        	iraf.tprint.columns="SMA"
                                        	sma_tab=iraf.tprint(tab_name,Stdout=1)
                                        	sma_interp=np.array(sma_tab).astype(float)
                                        	iraf.tprint.columns="INTENS"
                                        	intens_tab=iraf.tprint(tab_name,Stdout=1)
                                        	intens_interp=np.array(intens_tab).astype(float)
                                        	id_ellip=np.append(id_ellip,self.id_list[i])

				if self.opt_ellip==1:
					if (os.path.exists(tab_name)==True):
                                        	iraf.tprint.columns="ELLIP"
						ellip_tab=iraf.tprint(tab_name,Stdout=1)
                                        	ellip_interp=np.array(ellip_tab).astype(float)
                                        	iraf.tprint.columns="PA"
						pa_tab=iraf.tprint(tab_name,Stdout=1)
                                        	pa_interp=np.array(pa_tab).astype(float)
						grad_vec=np.ones(len_arr)*99
						for l in np.arange(0,len_arr-5,1):
							grad_vec[l]=abs(scipy.stats.linregress(sma_interp[l:l+5],ellip_interp[l:l+5])[0])
						grad_vec[grad_vec==0.]=99
                                        	mask_grad=grad_vec==np.min(grad_vec)
						ellip_out=ellip_interp[mask_grad][0]
						pa_out=pa_interp[mask_grad][0]
						if ellip_out<0.05:
                                        	       ellip_out=0.05
                                        	if pa_out<-90.:
                                        	       pa_out=-90.
                                        	if pa_out>90.:
                                        	       pa_out=90.
						ellip_array=np.append(ellip_array,ellip_out)
						pa_array=np.append(pa_array,pa_out)
                                        	os.system('mv '+tab_name+' '+os.path.join(self.data_dir,self.names[k]+'_'+str(self.id_list[i])+self.suff[:-5]+'_calc_ellip.tab')) 

					
				if ((self.opt_ellip==3) or (self.opt_ellip==4)) or (self.opt_ellip==5):
					if (self.opt_ellip==3) or (self.opt_ellip==4):
						image_name=os.path.join(self.data_dir,self.names[k]+'_'+str(self.id_list[i])+self.suff[:-4]+'tab')
                                                data_file_name=os.path.join(self.data_dir,self.names[k]+'_'+str(self.id_list[i])+self.suff[:-4]+'dat')
					if self.opt_ellip==5:
						image_name=os.path.join(self.data_dir,self.names[k]+'.tab')
						data_file_name=os.path.join(self.data_dir,self.names[k]+'.dat')
					if os.path.exists(image_name)==True:
						iraf.tprint.rows='1-'+str(len_arr)
						iraf.tprint.showro='no'
                                        	iraf.tprint.showhdr='no'
                                        	iraf.tprint.columns="SMA"
						sma_tmp=iraf.tprint(image_name,Stdout=1)
                                        	sma_array=np.array(sma_tmp).astype(float)
						iraf.tprint.columns="INTENS"
						intens_tmp=iraf.tprint(image_name,Stdout=1)
						iraf.tprint.columns="INT_ERR"
						int_err_tmp=iraf.tprint(image_name,Stdout=1)
						iraf.tprint.columns="NPIX_E"
						npix_tmp=iraf.tprint(image_name,Stdout=1)
						index=[m for m, e in enumerate(intens_tmp) if 'INDEF' in e]
						index2=[m for m, e in enumerate(int_err_tmp) if 'INDEF' in e]
						index3=[m for m, e in enumerate(npix_tmp) if 'INDEF' in e]
						len_indef=len(index)
						len_indef2=len(index2)
						len_indef3=len(index3)
						
						if (len_indef>=1):
							for m in np.arange(0,len_indef,1):
								intens_tmp[m]=0
						if (len_indef2>=1):
							for m in np.arange(0,len_indef2,1):
								int_err_tmp[m]=0.
						if (len_indef3>=1):
							for m in np.arange(0,len_indef3,1):
								npix_tmp[m]=1.
                                        	intens_array=np.array(intens_tmp).astype(float)
                                        	int_err_array=np.array(int_err_tmp).astype(float)
                                        	npix_array=np.array(npix_tmp).astype(float)
                                        	ellipse_out=open(data_file_name,'w')
						ellipse_out.write('#SMA\t INTENS\t INT_ERR\t NPIX\n')
						np.savetxt(ellipse_out,np.array([sma_array,intens_array,int_err_array,npix_array]).T,fmt='%14.5f %14.5f %14.5f %10i')
						ellipse_out.close()
						if (self.opt_ellip==3) or (self.opt_ellip==4):
							list_intens_prof.write(str(self.id_list[i])+'\t'+self.names[k]+'_'+str(self.id_list[i])+self.suff[:-4]+'dat\n')

			if self.opt_ellip==5:
				self.list_file.write(k+'\t'+self.names[k]+self.suff[:-4]+'dat\n')

			if self.opt_ellip==2:
				rfit_dic[k]=np.array([id_ellip,rfit_array,reff_array]).T
				fit_opt=np.ones(len(id_ellip)).astype(int).astype(str)
				np.savetxt(analysis_file,np.array(zip(id_ellip,x0_array,y0_array,rfit_array,reff_array,pa_array,ellip_array,fit_opt),dtype=[('id_ellip','S16'),('x0_array',int),('y0_array',int),('rfit_array',float),('reff_array',float),('pa_array',float),('ellip_array',float),('fit_opt','S16')]),fmt="%10s %10d %10i %9.5f %9.5f %9.5f %9.5f %4s")
			if self.opt_ellip==1:
				ellip_dic[k]=np.array([id_ellip,ellip_array,pa_array]).T
				np.savetxt(ellip_file,np.array(zip(id_ellip,ellip_array,pa_array),dtype=[('id_ellip','S16'),('ellip_array',float),('pa_array',float)]),fmt="%10s %9.5f %9.5f")

		if self.opt_ellip==1:
			return ellip_dic
		if (self.opt_ellip==3) or (self.opt_ellip==4):
                	self.list_file.close()
		if self.opt_ellip==5:
			self.list_file.close()

class fitting_radius:

	def __init__(self,names,data_dir,sky_dic,scale,opt_moff,opt_king,opt_wilson,lib_path):
		self.names=names
		self.data_dir=data_dir
		self.sky_dic=sky_dic
		self.scale=scale
		self.list=np.genfromtxt(os.path.join(self.data_dir,'intensity_prof_info.dat'),comments='#',dtype=str)
		self.list_psf=np.genfromtxt(os.path.join(self.data_dir,'psf_profiles_info.dat'),comments='#',dtype=str)
		self.lib_path=lib_path
		os.system('mkdir '+data_dir+'/rfit_plots')
		
		
		if np.size(self.list)==2:
			self.list=np.array([self.list])
			self.list_psf=np.array([self.list_psf])

		for i in np.arange(0,len(self.list[:,0]),1):
			rfit_file=open(data_dir+'/rfit_plots/rfit_data_'+self.list[i,0],'w')
			rfit_file.write('#ID\t Rip\t R3s\n')
			self.prof_list=np.genfromtxt(os.path.join(self.data_dir,self.list[i,1]),comments='#',dtype='S')
			if self.prof_list.size==2:
				self.prof_list=np.array([self.prof_list])
			ids=self.prof_list[:,0]
			
			files=self.prof_list[:,1]
			psf_file=self.list_psf[:,1][self.list_psf[:,0]==self.list[i,0]][0]
			psf_data=np.genfromtxt(os.path.join(self.data_dir,psf_file),comments='#')
			sma_psf=psf_data[:,0]*self.scale[self.list[i,0]]
			intens_psf=psf_data[:,1]
			psf_name=psf_file.split('.dat')[0]+'_file.dat'
			psf_path=os.path.join(self.data_dir,psf_name)
			np.savetxt(psf_path,np.array([sma_psf,intens_psf]).T,fmt='%9.3f %9.2e')
			for j in np.arange(0,len(ids),1):
				data=scipy.genfromtxt(os.path.join(self.data_dir,files[j]),comments='#')
				if np.size(sky_dic[self.list[i,0]])<=3:
					sky_rms=np.float(sky_dic[self.list[i,0]][2][0])
				else:
					mask_sky=sky_dic[self.list[i,0]][:,0]==ids[j]
					sky_rms=np.float(sky_dic[self.list[i,0]][mask_sky,2][0])
				resta_rfit=data[:,1]-np.float(sky_rms)	
				if len(resta_rfit[resta_rfit<0])>0:
                                	mask_rfit=resta_rfit[resta_rfit<0][0]
                                	rfit_index=np.where(resta_rfit==mask_rfit)[0][0]-1
				else:
					rfit_index=len(data[:,1])-1
				r3s=data[rfit_index,0]*self.scale[self.list[i,0]]
						
				t3=0
				m3=1
				m_ant=0.
				while t3==0:
					if m3<len(data[:,1])-2:

                                                pend=(data[m3+2,1]-data[m3,1])/(data[m3+2,0]-data[m3,0])

                                                if (data[m3,0]>5.):
                                                        if (m_ant<=pend) and (pend<0.):
                                                                m3=m3+1
                                                                m_ant=pend
                                                                t3=0
                                                        else:
                                                                m_ant=pend
                                                                t3=1
                                                else:
                                                        m_ant=pend
                                                        m3=m3+1
                                                        t3=0
                                        if m3==len(data[:,1])-2:
                                                m3=len(data[:,1])-1
                                                t3=1
				rip=data[:,0][m3]*self.scale[self.list[i,0]]
				sma=data[:,0]*self.scale[self.list[i,0]]
				im_size=np.int(np.round(sma[-1]*2+1))
				intens=data[:,1]
				int_err=data[:,2]
				npix=data[:,3]
				rfit=np.min([r3s,rip])
				plt.figure(1,(8,8))
				plt.plot(sma[1:],intens[1:],'ko')
				plt.axvline(r3s,label='r3s',color='yellow')
				plt.axvline(rip,label='rip',color='red')
				plt.gca().set_xscale('log')
				plt.legend(loc=0)
				tmp_name=self.list[i,0]+'_'+ids[j]
				file_name='object_'+tmp_name
				rfit_file.write(ids[j]+'\t'+str('%7.4f\t' %(rip))+str('%7.4f\n' %(r3s)))
				plt.savefig(os.path.join(data_dir,'rfit_plots/'+file_name+'.png'))
				plt.close()

				np.savetxt(os.path.join(data_dir,file_name+'.dat'),np.array([sma,intens,int_err,npix]).T,fmt='%9.2f %9.3f %9.4f %10d')

				if opt_moff==1:
					file_moffat='moffat_input_'+file_name+'.dat'
					moffat_file=open(os.path.join(self.data_dir,file_moffat),'w')
					moffat_file.write(str('%d\n' %(im_size)))
					moffat_file.write('0\n')
					moffat_file.write(psf_name+'\n')
					moffat_file.write(file_name+'.dat\n')
					moffat_file.write(str('%f\n' %(sky_rms)))
					moffat_file.write(str('%f\n' %(rfit)))
					moffat_file.write(tmp_name+'\n')
					moffat_file.write(self.lib_path+'\n')
					moffat_file.close()

				if opt_king==1:
					file_king='king_input_'+file_name+'.dat'
					king_file=open(os.path.join(self.data_dir,file_king),'w')
					king_file.write(str('%d\n' %(im_size)))
					king_file.write('1\n')
					king_file.write(psf_name+'\n')
					king_file.write(file_name+'.dat\n')
					king_file.write(str('%f\n' %(sky_rms)))
					king_file.write(str('%f\n' %(rfit)))
					king_file.write(tmp_name+'\n')
					king_file.write(self.lib_path+'\n')
					king_file.close()

				if opt_wilson==1:
					file_wilson='wilson_input_'+file_name+'.dat'
					wilson_file=open(os.path.join(self.data_dir,file_wilson),'w')
					wilson_file.write(str('%d\n' %(im_size)))
					wilson_file.write('2\n')
					wilson_file.write(psf_name+'\n')
					wilson_file.write(file_name+'.dat\n')
					wilson_file.write(str('%f\n' %(sky_rms)))
					wilson_file.write(str('%f\n' %(rfit)))
					wilson_file.write(tmp_name+'\n')
					wilson_file.write(self.lib_path+'\n')
					wilson_file.close()
			rfit_file.close()

class psf_data:

	def __init__(self,filters,input_file,data_dir):
		self.filters=filters
                self.input_file=input_file
		self.data_dir=data_dir
		self.psf_dic={}
                self.psf_info=scipy.genfromtxt(self.input_file[0],dtype='S',comments='#',usecols=0)
                self.psf_list=scipy.genfromtxt(self.input_file[0],dtype='S',comments='#',usecols=6)

	def psf_arrange(self):
                for k in self.filters:
                        self.psf_dic[k]=self.psf_list[self.psf_info==k][0]
		
		return self.psf_dic
