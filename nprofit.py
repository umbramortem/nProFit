import numpy as np
import scipy
import pylab as plt
import sys
import os
import glob
import sky as skyclass
import ellipse as ellipclass
import results as resu
import der_pars2 as der
#import plots as plotclass
#import profiles as profitclass

#--------------------------------Checking input file-----------------------------

try:
        nombreinput=sys.argv[1]
except IndexError:
        sys.exit('Syntax: %s inputfile\n See README for instructions' % (sys.argv[0]))

iso_path=sys.argv[0].split('nprofit.py')[0]

print '#----------------------Reading input files---------------------------'

#scale=0.88 #me falta agregar este en el input
input_file=scipy.genfromtxt(sys.argv[1],comments='#',dtype='S')
filters=scipy.genfromtxt(input_file[0],comments='#',usecols=0,dtype='S')
mag_zero={}
names={}
scale={}
arc_pix={}
m_l={}
M_sun={}
magzero_arr=scipy.genfromtxt(input_file[0],comments='#',usecols=1,dtype='S')
image_name=scipy.genfromtxt(input_file[0],comments='#',usecols=2,dtype='S')
ds9_name=scipy.genfromtxt(input_file[0],comments='#',usecols=3,dtype='S')
scale_file=scipy.genfromtxt(input_file[0],comments='#',usecols=7)
arc_pix_file=scipy.genfromtxt(input_file[0],comments='#',usecols=8)
m_l_file=scipy.genfromtxt(input_file[0],comments='#',usecols=9)
M_sun_file=scipy.genfromtxt(input_file[0],comments='#',usecols=10)

if (filters.size==1):
	filters=np.array([filters])
	magzero_arr=np.array([magzero_arr])
	image_name=np.array([image_name])
	ds9_name=np.array([ds9_name])
	scale_file=np.array([scale_file])
	arc_pix_file=np.array([arc_pix_file])
	m_l_file=np.array([m_l_file])
	M_sun_file=np.array([M_sun_file])

for k in filters:
	mag_zero[k]=magzero_arr[filters==k]
	names[k]=image_name[filters==k][0].split('.fits')[0]
	scale[k]=scale_file[filters==k]
	arc_pix[k]=arc_pix_file[filters==k]
	m_l[k]=m_l_file[filters==k]
	M_sun[k]=M_sun_file[filters==k]

cur_dir=os.getcwd()
data_dir=os.path.join(os.getcwd(),'nprofit_data')

if os.path.exists(data_dir)==False:
	os.system('mkdir nprofit_data')

#-----------------------------Fitting box size-----------------------------------

if '1' in input_file[1]:
	fit_box=input_file[2]	
	box_opt=1

if '2' in input_file[1]:
	fit_box={}
	box_info=scipy.genfromtxt(input_file[2],comments='#',dtype='S')
	box_opt=2

	for k in box_info[:,0]:
		fit_box[k]=scipy.genfromtxt(box_info[:,1][box_info[:,0]==k][0],comments='#',dtype='S')


#--------------------------------Coordinates transformation---------------------

if '1' in input_file[4]:

	print '#--------------------------Reading coordinates file--------------------'
	coord_file=scipy.genfromtxt(input_file[3],comments='#',dtype='S')

if '2' in input_file[4]:

	print '#----------------------------Transforming to WCS-------------------------'
	wcsxy=skyclass.coords_trans(input_file[3],filters,names)
	wcs_out=wcsxy.wcsxymatch()

	coord_file=scipy.genfromtxt(wcs_out,comments='#',dtype='S')

if coord_file.size<=3:
	coord_file=np.array([coord_file])
	
#--------------------------------Sky substraction--------------------------------

sky=skyclass.cut_image(coord_file,filters,fit_box,names,box_opt,cur_dir,data_dir)

if 'yes' in input_file[5]:
	print '#--------------------Cutting image in sub-images---------------------'
	sky.boxes()

if (('no' in input_file[6]) & ('yes' in input_file[7])):
        if '1' in input_file[8]:
                print '#--------------Measuring sky rms and background values with instat--------------' 
        	sky_dic=sky.instat()
        if '2' in input_file[8]:
                print '#-----------------------Measuring sky rms and background values with median method------------------' 
		sky_dic=sky.median_sky()

if 'yes' in input_file[6]:

	if 'yes' in input_file[7]:
                if '1' in input_file[8]:
                        print '#--------------Measuring sky rms and background values with instat--------------' 
                	sky_dic=sky.instat()
                if '2' in input_file[8]:
                        print '#-----------------------Measuring sky rms and background values with median method------------------' 
			sky_dic=sky.median_sky()

	#print sky_dic[k]
	if 'yes' in input_file[9]:
	        print '#----------------------Reading sky file information-----------------'
		#sky_info=scipy.genfromtxt(input_file[10],comments='#',dtype='S')
                sky_info=scipy.genfromtxt(input_file[0],comments='#',usecols=4,dtype='S')
		if sky_info.size==1:
			sky_info=np.array([sky_info])
			#print sky_info,'bla'
		sky_dic={}
		for k in filters:
			#print scipy.genfromtxt(sky_info[filters==k][0],comments='#',dtype='S')
			sky_dic[k]=scipy.genfromtxt(sky_info[filters==k][0],comments='#',dtype='S')
	
	print '#------------------------Substracting sky------------------------'
	substract=skyclass.substract_sky(sky_dic,fit_box,names,box_opt,data_dir)	

	substract.read_sky()			

if 'no' in input_file[6]:
	#if ('yes' in input_file[7]) | ('yes' in input_file[12]) | ('yes' in input_file[10]) | ('yes' in input_file[13]) | ('yes' in input_file[14]) | ('yes' in input_file[15]) | ('yes' in input_file[15]) | ('yes' in input_file[25]) | ('yes' in input_file[26]):
	if ('yes' in input_file[7]) | ('yes' in input_file[12]) | ('yes' in input_file[10]) | ('yes' in input_file[13]):

	        print '#----------------------Reading sky file information-----------------'
                sky_info=scipy.genfromtxt(input_file[0],comments='#',usecols=4,dtype='S')
		if sky_info.size==1:
			sky_info=np.array([sky_info])
		sky_dic={}
		for k in filters:
			sky_dic[k]=scipy.genfromtxt(sky_info[filters==k][0],comments='#',dtype='S')
                suff='_sky_sub.fits'

if 'yes' in input_file[6]:
        suff='_sky_sub.fits'

#--------------------------------Masking up contaminants--------------------------

if 'yes' in input_file[10]:
        print '#-----------------------Masking up images-------------------------'
	#suff_mask=suff.split('.fits')[0]
	suff_mask=''
	mask_file=input_file[11]
	maps=ellipclass.pixel_maps(filters,mask_file,names,fit_box,box_opt,data_dir,suff_mask)
	maps.mask()

#------------------------------Computing ellipticity----------------------------------

if 'yes' in input_file[12]:
	print '#-----------------------Reading ellipticiy and P.A. from information file----------------'
	ellip_info=scipy.genfromtxt(input_file[0],comments='#',usecols=5,dtype='S')
	if ellip_info.size==1:
		ellip_info=np.array([ellip_info])
	ellip_dic={}

	for k in filters:
        	#ellip_dic[k]=scipy.genfromtxt(ellip_info[:,1][ellip_info[:,0]==k][0],comments='#',dtype='S')
		ellip_dic[k]=scipy.genfromtxt(ellip_info[filters==k][0],comments='#',dtype='S')

id_ellip=scipy.genfromtxt(input_file[3],comments='#',dtype='S',usecols=0)	

if id_ellip.size==1:
	id_ellip=np.array([id_ellip])

if 'yes' in input_file[13]:
        print '#----------------------Calculating ellipticity and P.A.--------------------'
	opt_ellip=1
	ellip_dic={}
	ellip_calc=ellipclass.ellipse(filters,opt_ellip,names,id_ellip,mag_zero,fit_box,ellip_dic,sky_dic,coord_file,suff,data_dir)
	ellip_dic=ellip_calc.task_ellip()
		
#--------------------------------Intensity profiles----------------------------------

if 'yes' in input_file[14]:

	print '#------------------------Computing intensity profiles-------------------'
	if 'yes' in input_file[13]:
		opt_ellip=3
	
	if 'no' in input_file[13]:
		opt_ellip=4

	if 'yes' in input_file[15]:
		el_thres=float(input_file[16])
		for k in filters:
			ellip_dic[k][:,1][ellip_dic[k][:,1].astype(float)>el_thres]=el_thres
			
	ellipse_file={}
	for k in filters:
		ellipse_file[k]=np.array([coord_file[:,0],coord_file[:,1],coord_file[:,2],ellip_dic[k][:,1],ellip_dic[k][:,2]]).T
	
	profile_calc=ellipclass.ellipse(filters,opt_ellip,names,id_ellip,mag_zero,fit_box,ellipse_file,sky_dic,ellipse_file,suff,data_dir)
	profile_calc.task_ellip()
	
#-----------------------------------PSF image-------------------------------------

if 'no' in input_file[17]:
	conv_opt='0'
        psf_name='none'
	psf_dic={}
	for k in filters:
		psf_dic[k]=psf_name

if 'yes' in input_file[17]:
	print '#----------------------------Setting PSF image-----------------------------'
	conv_opt='1'
        psf_opt=input_file[18]
	psf_class=ellipclass.psf_data(filters,input_file,data_dir)

	if '1' in psf_opt:
                psf_names={}
		psf_dic={}
		psf_dic=psf_class.psf_arrange()
                opt_ellip=5
                id_psf=np.array(['psf_1'])
                sky_empty={}
                ellipse_psf={}
                psf_suff='.fits'
                x0_psf=-99

		for k in psf_dic.keys():
			psf_names[k]=psf_dic[k].split('.fits')[0]
                        ellipse_psf[k]=0.05

                profile_psf=ellipclass.ellipse(filters,opt_ellip,psf_names,id_psf,mag_zero,x0_psf,ellipse_psf,sky_empty,x0_psf,psf_suff,data_dir)
		profile_psf.task_ellip()

##-------------------------------Dynamical models fitting-----------------------------------

if '1' in input_file[24]:
	if '1' in input_file[19]:
	
		opt_mof,opt_king,opt_wilson=0,0,0
	
		if 'yes' in input_file[20]:
			opt_mof=1
		if 'yes' in input_file[21]:
			opt_king=1
		if 'yes' in input_file[22]:
			opt_wilson=1
		
		lib_path=input_file[23]
	
		print '#-------------------------Setting up data for fitting--------------------------------'
		ellipclass.fitting_radius(names,data_dir,sky_dic,scale,opt_mof,opt_king,opt_wilson,lib_path)
		
	if ('1' in input_file[19]) or ('2' in input_file[19]):
	
		os.chdir(data_dir)
	
	        if 'yes' in input_file[20]:
			for k in filters:
				fact=1./(arc_pix[k])**2
				mof_file=open(os.path.join(cur_dir,'moffat_pars_'+k+'.dat'),'w')
				mof_file.write('#ID\t chi\t Npts\t rd\t rd_err_down\t rd_err_up\t gamma\t gamma_err_down\t gamma_err_up\t mu_0\n')
				list_moff=sorted(glob.glob('moffat_input_object_'+k+'*'))
				for item in list_moff:
					iter_str=item.split('moffat_input_object_'+k+'_')[-1].split('.dat')[0]
					print('#------------Fitting Moffat-EFF model to object '+iter_str+' in the '+k+'-band-----------------')
					#os.system(os.path.join(cur_dir,'./isothermal')+' < '+item)
					os.system(os.path.join(iso_path,'isothermal')+' < '+item)
					pars=np.genfromtxt('moffat_pars_object_'+k+'_'+iter_str+'.dat',comments='#')
					Io=np.genfromtxt('object_'+k+'_'+iter_str+'.dat',comments='#',usecols=1)
					mu_o=-2.5*np.log10(Io[0]*fact)+np.float(mag_zero[k])
					mof_file.write(iter_str+'\t'+str('%9.4f' %(pars[0]))+'\t'+str('%d' %(pars[1]/scale[k]))+'\t')
					mof_file.write(str('%9.4f' %(pars[2]))+'\t'+str('%9.4f' %(pars[2]-pars[4]))+'\t'+str('%9.4f' %(pars[5]-pars[2]))+'\t')
					mof_file.write(str('%9.4f' %(pars[3]))+'\t'+str('%9.4f' %(pars[3]-pars[6]))+'\t'+str('%9.4f' %(pars[7]-pars[3]))+'\t'+str('%9.4f' %(mu_o))+'\n')
				mof_file.close()
	
	        if 'yes' in input_file[21]:
			for k in filters:
				king_file=open(os.path.join(cur_dir,'king_pars_'+k+'.dat'),'w')
				king_file.write('#ID\t chi\t Npts\t ro\t ro_err_down\t ro_err_up\t W0\t W0_err_down\t W0_err_up\t mu_o\n')
				list_king=sorted(glob.glob('king_input_object_'+k+'*'))
				for item in list_king:
					iter_str=item.split('king_input_object_'+k+'_')[-1].split('.dat')[0]
					print('#------------Fitting King model to object '+iter_str+' in the '+k+'-band-----------------')
					#os.system('../isothermal < '+item)
					#os.system(os.path.join(cur_dir,'./isothermal')+' < '+item)
					os.system(os.path.join(iso_path,'isothermal')+' < '+item)
					pars=np.genfromtxt('king_pars_object_'+k+'_'+iter_str+'.dat',comments='#')
					Io=np.genfromtxt('object_'+k+'_'+iter_str+'.dat',comments='#',usecols=1)
					mu_o=-2.5*np.log10(Io[0]*fact)+np.float(mag_zero[k])
					king_file.write(iter_str+'\t'+str('%9.4f' %(pars[0]))+'\t'+str('%d' %(pars[1]/scale[k]))+'\t')
					king_file.write(str('%9.4f' %(pars[2]))+'\t'+str('%9.4f' %(pars[2]-pars[4]))+'\t'+str('%9.4f' %(pars[5]-pars[2]))+'\t')
					king_file.write(str('%9.4f' %(pars[3]))+'\t'+str('%9.4f' %(pars[3]-pars[6]))+'\t'+str('%9.4f' %(pars[7]-pars[3]))+'\t'+str('%9.4f' %(mu_o))+'\n')
				king_file.close()
	        
	        if 'yes' in input_file[22]:
			for k in filters:
				wilson_file=open(os.path.join(cur_dir,'wilson_pars_'+k+'.dat'),'w')
				wilson_file.write('#ID\t chi\t Npts\t ro\t ro_err_down\t ro_err_up\t W0\t W0_err_down\t W0_err_up\n')
				list_wilson=sorted(glob.glob('wilson_input_object_'+k+'*'))
				for item in list_wilson:
					iter_str=item.split('wilson_input_object_'+k+'_')[-1].split('.dat')[0]
					print('#------------Fitting Wilson model to object '+iter_str+' in the '+k+'-band-----------------')
					#os.system('../isothermal < '+item)
					os.system(os.path.join(iso_path,'isothermal')+' < '+item)
					pars=np.genfromtxt('wilson_pars_object_'+k+'_'+iter_str+'.dat',comments='#')
					Io=np.genfromtxt('object_'+k+'_'+iter_str+'.dat',comments='#',usecols=1)
					mu_o=-2.5*np.log10(Io[0]*fact)+np.float(mag_zero[k])
					wilson_file.write(iter_str+'\t'+str('%9.4f' %(pars[0]))+'\t'+str('%d' %(pars[1]/scale[k]))+'\t')
					wilson_file.write(str('%9.4f' %(pars[2]))+'\t'+str('%9.4f' %(pars[2]-pars[4]))+'\t'+str('%9.4f' %(pars[5]-pars[2]))+'\t')
					wilson_file.write(str('%9.4f' %(pars[3]))+'\t'+str('%9.4f' %(pars[3]-pars[6]))+'\t'+str('%9.4f' %(pars[7]-pars[3]))+'\t'+str('%9.4f' %(mu_o))+'\n')
				wilson_file.close()

	os.chdir(cur_dir)

#-------------------------------Plotting profiles--------------------------------

#if '3' in input_file[24]:
mod_fit=[]
if 'yes' in input_file[20]:
	mod_fit.append('moffat')
if 'yes' in input_file[21]:
	mod_fit.append('king')
if 'yes' in input_file[22]:
	mod_fit.append('wilson')

if 'yes' in input_file[25]:
	
	os.system('mkdir plots_profiles')
	plots_dir=os.path.join(cur_dir,'plots_profiles')
	ds9_path=input_file[26]
	resu.mosaics(filters,ds9_name,ds9_path,plots_dir,data_dir,names,scale)
	resu.plots(filters,plots_dir,mod_fit,data_dir,mag_zero,sky_dic,scale,arc_pix)

if ('yes' in input_file[27]) & (len(mod_fit)>=1):
	
	print('#---------------------Computing derived parameters------------------')
	der.derived_pars(filters,cur_dir,mod_fit,m_l,M_sun)	

