import numpy as np
import scipy
import pylab
from astropy.io import fits
import astropy.convolution 
import sys
import os
import scipy.interpolate
import traceback

try:
	nombreinput=sys.argv[1]
except IndexError:
        sys.exit('Syntax: %s inputfile see README for instructions' % (sys.argv[0]))

def interp_isothermal(par1,par2,mu,opt,pc_scale):
	W0=par1
        r0=par2*pc_scale
        W0_str=str('%4.1f' %(W0))
	if (opt==1):
                file_name='/Users/boliviacuevas/Documents/nprofit/nprofit_library/king_W0_'+W0_str+'.dat'
        if (opt==2):
                file_name='/Users/boliviacuevas/Documents/nprofit/nprofit_library/wilson_W0_'+W0_str+'.dat'
        file_name=file_name.replace(" ","")
        if os.path.exists(file_name)==True:
                iso_file=np.genfromtxt(file_name)
		inter=scipy.interpolate.interp1d(iso_file[:,0]*r0/pc_scale,iso_file[:,3])
		
	return inter,iso_file[-1,0]*r0

def models(x,par1,par2,mu,opt,inter):
	if opt==0:
		gamma=par1
		rd=par2
		return mu*(1+(x/rd)**2)**(-gamma/2.)
	if (opt==1) | (opt==2):
		return inter(x)/np.max(inter(x))*mu

ofile=open(sys.argv[1],'r')
pc_scale=np.float(ofile.readline().split('#')[0])
bg_val=np.float(ofile.readline().split('#')[0])
bg_rms=np.float(ofile.readline().split('#')[0])
im_size_x=np.int(ofile.readline().split('#')[0])
im_size_y=np.int(ofile.readline().split('#')[0])
psf_file=ofile.readline().split('#')[0].split('\t')[0]
hdu_psf=fits.open(psf_file)
psf=hdu_psf[0].data
lists=np.genfromtxt(sys.argv[1],comments='#',skip_header=6,dtype='S')

bg_arr=np.abs(scipy.random.normal(bg_val,bg_rms,im_size_x*im_size_y))
im_tmp=bg_arr.reshape(im_size_x,im_size_y)

x_vec=np.arange(0,im_size_x,1)
y_vec=np.arange(0,im_size_y,1)

reg_file=open('mock_sample.reg','w')
reg_file.write('# Region file format: DS9 version 4.0\n')
reg_file.write('#global color=green font="helvetica 10 normal"')
reg_file.write(' select=1 highlite=1 edit=1 move=1 delete=1 include=1 fixed=0 source\n')
reg_file.write('#fk5\n')

print('######Starting the mock sample simulation########')

for row in np.arange(0,len(lists[:,0]),1):
	ids=lists[row,0]
	x_c=np.int(lists[row,1])
	y_c=np.int(lists[row,2])
	reg_file.write('circle('+str('%3i' %(x_c))+','+str('%3i' %(y_c))+',10)')
	reg_file.write(' # color=green text={'+ids+'}\n')
	model=lists[row,3]
	mu=np.float(lists[row,4])
	par1=np.float(lists[row,5])
	par2=np.float(lists[row,6])
	ellip=np.float(lists[row,7])
        AR=np.sqrt(1-ellip**2)
	PA=np.float(lists[row,8])

	if 'moffat' in model:
		opt=0
		inter=0
	if 'king' in model:
		opt=1
		
	if 'wilson' in model:
		opt=2

	if (opt==1) | (opt==2):
		try:
			inter,rmax=interp_isothermal(par1,par2,mu,opt,pc_scale)
		except Exception,err:
			print('#####Error####')
			print('You have provided an incorrect path to the nprofit_library')
			print('Correct the path and try again')
			exit()

	x_mesh,y_mesh=np.meshgrid(y_vec+1-y_c,x_vec+1-x_c)
        x_pa=x_mesh*np.cos((90-PA)*np.pi/180)+y_mesh*np.sin((90-PA)*np.pi/180)
        y_pa=(-x_mesh*np.sin((90-PA)*np.pi/180)+y_mesh*np.cos((90-PA)*np.pi/180))/AR
	r=np.sqrt(x_pa**2+y_pa**2)

	if opt==0:
		mod_tmp=models(r,par1,par2,mu,opt,inter)
	else:
		r_mod=np.where(r<=rmax,r,rmax)
		mod_zero=np.zeros(im_size_x*im_size_y).reshape(im_size_x,im_size_y)
		mod_tmp=np.where(r_mod>rmax,mod_zero,models(r_mod,par1,par2,mu,opt,inter))

	im_tmp=im_tmp+mod_tmp

im_conv=astropy.convolution.convolve(im_tmp,psf,boundary='extend')

hdu=fits.PrimaryHDU(im_conv.T)

hdr = hdu.header
hdr['EXPTIME']=1.
#hdul=fits.HDUList([hdu])
hdul=fits.HDUList([hdu])
hdul.writeto('mock_sample.fits')
reg_file.close()
print('#######Run complete#########')
