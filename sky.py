import numpy as np
import scipy
import pylab as plt
import os
from pyraf import iraf

class coords_trans:

	def __init__(self,coords_name,filters,names):
	
		self.filters=filters
		self.coords_name=coords_name
		self.names=names

	def wcsxymatch(self):

		iraf.wcsxymatch.xcolumn=2
		iraf.wcsxymatch.ycolumn=3
		iraf.wcsxymatch.coords=self.coords_name
		iraf.verbose='no'
		input_wcs=self.names[self.filters[0]]
		ref_wcs=self.names[self.filters[0]]
		iraf.delete('coords_wcs.dat')
		output_wcs='coords_wcs.dat'
		iraf.wcsxymatch(input=input_wcs,referenc=ref_wcs,output=output_wcs)
		id_coords=scipy.genfromtxt(self.coords_name,usecols=0,dtype='S')
		x_wcs=np.round(scipy.genfromtxt(output_wcs,usecols=0).astype(float))
		y_wcs=np.round(scipy.genfromtxt(output_wcs,usecols=1).astype(float))
		iraf.delete('coords_wcs.dat')
		np.savetxt(output_wcs,np.array(zip(id_coords,x_wcs,y_wcs),dtype=[('id_coords','S16'),('x_wcs',int),('y_wcs',int)]),fmt="%10s %10d %10d")
		return output_wcs

class cut_image:

	def __init__(self,coords,filters,box_size,names,box_opt,cur_dir,data_dir):
	
		self.filters=filters
		self.id_coor=coords[:,0]
		self.x_coor=np.round(coords[:,1].astype(float)).astype(int)
		self.y_coor=np.round(coords[:,2].astype(float)).astype(int)
		self.box_opt=box_opt
		if self.box_opt==1:
			self.box_tmp=int(round(int(box_size)/2.))
		if self.box_opt==2:
			self.box_tmp_dic=box_size
                self.names=names
		self.suff='.fits'
                self.data_dir=data_dir
		self.cur_dir=cur_dir
		
	def boxes(self):

		for k in self.filters:

        		for i in np.arange(0,len(self.id_coor),1):
				if self.box_opt==1:
					self.box_size=self.box_tmp	
				if self.box_opt==2:
					self.box_size=int(round(int(self.box_tmp_dic[k][:,1][self.box_tmp_dic[k][:,0]==self.id_coor[i]][0])/2.))
				xmin=str(self.x_coor[i]-self.box_size)
               			ymin=str(self.y_coor[i]-self.box_size)
               		 	xmax=str(self.x_coor[i]+self.box_size)
               		 	ymax=str(self.y_coor[i]+self.box_size)
               		 	iraf.imcopy.verbose='no'
               		 	iraf.imcopy(self.names[k]+self.suff+'['+xmin+':'+xmax+','+ymin+':'+ymax+']',self.names[k]+'_'+self.id_coor[i]+self.suff)
                                os.system('mv '+self.names[k]+'_'+self.id_coor[i]+self.suff+' '+self.data_dir)

	def median_sky(self):
	
		sky_dic={}	
		for k in self.filters:

			sky_file=open('sky_'+k+'.dat','w')
			sky_file.write('#ID\t bg\t rms\n')
			id_sky=np.array([])
			sky_array=np.array([])
			sky_rms_array=np.array([])
			
        		for i in np.arange(0,len(self.id_coor),1):

				if self.box_opt==1:
					self.box_size=self.box_tmp	
				if self.box_opt==2:
					self.box_size=int(round(int(self.box_tmp_dic[k][:,1][self.box_tmp_dic[k][:,0]==self.id_coor[i]][0])/2.))
				iraf.imstat.fields='midpt'
                		sky=np.array([])
				sky_box=int(round(self.box_size*0.1))
				edge2=int(round(self.box_size+1-sky_box))
				edge3=int(self.box_size+1)
				box1='[1:'+str(sky_box)+',1:'+str(self.box_size*2)+']'	
				box2='[1:'+str(self.box_size*2)+',1:'+str(sky_box)+']'	
				box3='['+str(self.box_size*2-sky_box)+':'+str(self.box_size*2)+',1:'+str(self.box_size*2)+']'	
				box4='[1:'+str(self.box_size*2)+','+str(self.box_size*2-sky_box)+':'+str(self.box_size*2)+']'	
                		sky_1=iraf.imstat(self.data_dir+'/'+self.names[k]+'_'+self.id_coor[i]+self.suff+box1,Stdout=1)
                		sky_2=iraf.imstat(self.data_dir+'/'+self.names[k]+'_'+self.id_coor[i]+self.suff+box2,Stdout=1)
                		sky_3=iraf.imstat(self.data_dir+'/'+self.names[k]+'_'+self.id_coor[i]+self.suff+box3,Stdout=1)
                		sky_4=iraf.imstat(self.data_dir+'/'+self.names[k]+'_'+self.id_coor[i]+self.suff+box4,Stdout=1)
				iraf.imstat.fields='stddev'
                                sky_stddev=np.array([])
                		sky_1_stddev=iraf.imstat(self.data_dir+'/'+self.names[k]+'_'+self.id_coor[i]+self.suff+box1,Stdout=1)
                		sky_2_stddev=iraf.imstat(self.data_dir+'/'+self.names[k]+'_'+self.id_coor[i]+self.suff+box2,Stdout=1)
                		sky_3_stddev=iraf.imstat(self.data_dir+'/'+self.names[k]+'_'+self.id_coor[i]+self.suff+box3,Stdout=1)
                		sky_4_stddev=iraf.imstat(self.data_dir+'/'+self.names[k]+'_'+self.id_coor[i]+self.suff+box4,Stdout=1)

                		sky=np.append(sky,float(sky_1[1]))
                		sky=np.append(sky,float(sky_2[1]))
                		sky=np.append(sky,float(sky_3[1]))
                		sky=np.append(sky,float(sky_4[1]))
                		sky_stddev=np.append(sky_stddev,float(sky_1_stddev[1]))
                		sky_stddev=np.append(sky_stddev,float(sky_2_stddev[1]))
                		sky_stddev=np.append(sky_stddev,float(sky_3_stddev[1]))
                		sky_stddev=np.append(sky_stddev,float(sky_4_stddev[1]))
				sky_file.write(self.id_coor[i])
				sky_file.write('\t')
				sky_file.write(str(np.min(sky)))
				sky_file.write('\t')
				sky_file.write(str(np.min(sky_stddev)))
				sky_file.write('\n')
				id_sky=np.append(id_sky,self.id_coor[i])
				sky_array=np.append(sky_array,np.min(sky))
				sky_rms_array=np.append(sky_rms_array,np.min(sky_stddev))

			sky_dic[k]=np.array([id_sky,sky_array,sky_rms_array]).T
			if len(id_sky)==1:
				sky_dic[k]=np.array([id_sky,sky_array,sky_rms_array])
				
			sky_file.close()				
		return sky_dic

	def instat(self):
		iraf.task(instat="instat.cl")
		instat_dic={}
                for k in self.filters:

                        instat_file=open('sky_'+k+'.dat','w')
                        instat_file.write('#ID\t bg\t rms\n')
                        id_instat=np.array([])
                        bg_instat_array=np.array([])
                        rms_instat_array=np.array([])

                        for i in np.arange(0,len(self.id_coor),1):
				iraf.instat.verbose='no'
				
				iraf.instat(self.data_dir+'/'+self.names[k]+'_'+self.id_coor[i]+self.suff)
				bg=iraf.instat.bgvalue
				rms=iraf.instat.rmsvalue
				instat_file.write(self.id_coor[i])
				instat_file.write('\t')
				instat_file.write(str(bg))
				instat_file.write('\t')
				instat_file.write(str(rms))
				instat_file.write('\n')
                                id_instat=np.append(id_instat,self.id_coor[i])
                        	bg_instat_array=np.append(bg_instat_array,bg)
                        	rms_instat_array=np.append(rms_instat_array,rms)
			instat_dic[k]=np.array([id_instat,bg_instat_array,rms_instat_array]).T
                        instat_file.close()
		return instat_dic

class substract_sky:

	def __init__(self,sky_dic,box_size,names,box_opt,data_dir):
		self.sky_dic=sky_dic
		self.box_opt=box_opt
                if self.box_opt==1:
                        self.box_tmp=str(int(box_size)+1)
                if self.box_opt==2:
                        self.box_tmp_dic=box_size
                self.names=names
                self.suff='.fits'
                self.data_dir=data_dir

	def read_sky(self):
		for k in self.sky_dic.keys():		
			if self.box_opt==1:
				self.sky_box=self.box_tmp	
			if self.box_opt==2:
				self.sky_box=str(int(self.box_tmp_dic[k][:,1][self.box_tmp_dic[k][:,0]==i][0])+1)
			if np.size(self.sky_dic[k])<=3:
				i=self.sky_dic[k][0][0]
				sky_value=self.sky_dic[k][1][0]
				iraf.imexpr.dims=self.sky_box+','+self.sky_box
                        	iraf.imexpr.verbose='no'
                		iraf.imexpr('a-b',self.data_dir+'/'+self.names[k]+'_'+i+'_sky_sub'+self.suff,self.data_dir+'/'+self.names[k]+'_'+i+self.suff,sky_value)
			else:
				for i in self.sky_dic[k][:,0]: 
					sky_value=self.sky_dic[k][:,1][self.sky_dic[k][:,0]==i][0]
					iraf.imexpr.dims=self.sky_box+','+self.sky_box
                        		iraf.imexpr.verbose='no'
                			iraf.imexpr('a-b',self.data_dir+'/'+self.names[k]+'_'+i+'_sky_sub'+self.suff,self.data_dir+'/'+self.names[k]+'_'+i+self.suff,sky_value)
