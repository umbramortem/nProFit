import numpy as np
import scipy
import pylab as plt
import os
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.gridspec import GridSpec
from PIL import Image
import matplotlib.ticker as ticker
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import MultipleLocator

class mosaics:

	def __init__(self,filters,ds9,ds9_path,plots_dir,data_dir,names,scale):
	
		ds9_dic={'R':'-red','G':'-green','B':'-blue'}
		cmd_ds9=' -view colorbar no -scale log -scale mode zmax'
		files_path={}
		top_vec={}
		id_files={}
		for k in np.arange(0,len(filters),1):
			files_path[ds9[filters[k]]]=np.genfromtxt(os.path.join(data_dir,'list_intens_profiles_'+filters[k]+'.dat'),usecols=1,dtype='S')
			id_files[ds9[filters[k]]]=np.genfromtxt(os.path.join(data_dir,'list_intens_profiles_'+filters[k]+'.dat'),usecols=0,dtype='S')
	
			if files_path[ds9[filters[k]]].size==1:
				top_vec[ds9[filters[k]]]=np.array([files_path[ds9[filters[k]]]])
				id_files[ds9[filters[k]]]=np.array([id_files[ds9[filters[k]]]])
			else:	
				top_vec[ds9[filters[k]]]=files_path[ds9[filters[k]]]
	
		for i in np.arange(0,len(top_vec[ds9[filters[k]]]),1):
			cmd=''
			if filters.size==1:
				cmd_ds9_2=' '
				cmd=cmd+' '+os.path.join(data_dir,top_vec[ds9[filters[0]]][i].split('.dat')[0]+'.fits')+' '
			else:
				cmd_ds9_2=' -rgb '
				if filters.size==2:
					for k in ds9_dic.keys():
						if k not in ds9.keys():
							cmd=cmd+ds9_dic[k]
							cmd=cmd+' '+os.path.join(data_dir,top_vec[ds9[0]][i].split('.dat')[0]+'.fits')+' '
					for k in ds9.keys():
						cmd=cmd+ds9_dic[k]
						cmd=cmd+' '+os.path.join(data_dir,top_vec[ds9[k]][i].split('.dat')[0]+'.fits')+' '
			
				if filters.size==3:
					for k in ds9.keys():
						cmd=cmd+ds9_dic[ds9[k]]
						cmd=cmd+' '+os.path.join(data_dir,top_vec[ds9[k]][i].split('.dat')[0]+'.fits')+' '
			#top_vec[ds9[k]][i].split
			plots_path=os.path.join(plots_dir,'mosaic_'+id_files[ds9[filters[0]]][i]+'.png')
			
			cmd_end='-zoom to fit -height 450 -width 450 -saveimage png '+plots_path+' -exit'
                        if os.path.exists(plots_path)==False:
			        os.system(ds9_path+cmd_ds9+cmd_ds9_2+cmd+cmd_end)
                                img=Image.open(plots_path)
                                arr=np.array(img)
                                for j in np.arange(0,3,1):
                                        tmp=arr[:,:,j]
                                        arr_tmp=tmp[tmp!=255]
                                        n_size=int(np.sqrt(np.shape(arr_tmp)[0]))
                                        arr2=arr_tmp.reshape(n_size,n_size)
                                        if j==0:
                                                arr_3D=np.zeros((n_size,n_size,3))
                                        arr_3D[:,:,j]=arr2
                                new_im=Image.fromarray(arr_3D.astype(np.uint8))
                                new_im.save(plots_path)

class plots:
	
	def __init__(self,filters,plots_dir,models_list,data_dir,mag_zero,sky_dic,scale,arc_pix):
		
		pdf_pages=PdfPages(os.path.join(plots_dir,'profile_plots.pdf'))
		nums_mod=[0,4,6]
		ln_dic={}
		list_pars={}
		data_mod={}
		data_obs={}
		sky_rms={}
		rfit={}
		chi={}
		filter_colors=['blue','green','red']
		models_lab={'moffat':'Moffat-EFF','king':'King','wilson':'Wilson'}
		row_labels_dic={'king':[r'$r_0$',r'$W_0$'],'moffat':[r'$r_d$',r'$\gamma$'],'wilson':[r'$r_0$',r'$W_0$']}
		st=4
		filter_dic={}
		for i in np.arange(0,len(filters),1):
			filter_dic[filters[i]]=filter_colors[i]

		#print models_list
		for m in np.arange(0,len(models_list),1):
			ln_dic[models_list[m]]=nums_mod[m]
			list_pars[models_list[m]]={}
			rfit[models_list[m]]={}
			chi[models_list[m]]={}
			for k in filters:
				list_pars[models_list[m]][k]=np.genfromtxt(models_list[m]+'_pars_'+k+'.dat',comments='#',dtype='S')
				if np.size(list_pars[models_list[m]][k])<=10:
					rfit[models_list[m]][k]=np.float(list_pars[models_list[m]][k][2])
					chi[models_list[m]][k]=np.float(list_pars[models_list[m]][k][1])
				else:	
					rfit[models_list[m]][k]=list_pars[models_list[m]][k][:,2].astype(float)
					chi[models_list[m]][k]=list_pars[models_list[m]][k][:,1].astype(float)
		ln_dic['image']=2
		size_dic=len(ln_dic.keys())

		if np.size(list_pars[models_list[m]][k])<=10:
			ids=np.array([list_pars[models_list[0]][filters[0]][0]])
		else:
			ids=list_pars[models_list[0]][filters[0]][:,0]
                m=0

		#print list_pars	
		for i in np.arange(0,len(ids),1):
			min_int=1e6
			max_int=0.
			data_obs[ids[i]]={}
			data_mod[ids[i]]={}
			sky_rms[ids[i]]={}
			
			for k in filters:
				#print sky_dic[k]
				#print ids,'ids',len(ids)
				if len(ids)<=1:
					sky_rms[ids[i]][k]=np.float(sky_dic[k][2])
				else:
					sky_rms[ids[i]][k]=np.float(sky_dic[k][:,2][sky_dic[k][:,0]==ids[i]][0])
				data_obs[ids[i]][k]=np.genfromtxt(os.path.join(data_dir,'object_'+k+'_'+ids[i]+'.dat'),comments='#')
				rmax=data_obs[ids[i]][k][-1,0]/scale[k]
				data_mod[ids[i]][k]={}
				mask=(np.isnan(data_obs[ids[i]][k][:,1])==False) & (data_obs[ids[i]][k][:,1]>0)
				if np.min(data_obs[ids[i]][k][:,1])<min_int:
					min_int=np.min(data_obs[ids[i]][k][mask,1])
					min_filter=k
				if np.max(data_obs[ids[i]][k][:,1])>max_int:
					max_int=np.max(data_obs[ids[i]][k][mask,1])
					max_filter=k

                        fact_min=1./(arc_pix[min_filter][0]**2)
                        fact_max=1./(arc_pix[max_filter][0]**2)
			minmag=-2.5*np.log10(min_int*fact_min)+np.float(mag_zero[min_filter][0])
			maxmag=-2.5*np.log10(max_int*fact_max)+np.float(mag_zero[max_filter][0])
                        print minmag,maxmag
			
			if m>9.:
				m=0
			if m==0:	
				fig = plt.figure(figsize=(8.67,12.5), dpi=100)
                		fig.subplots_adjust(hspace=0., wspace=0.)

			for mod in ln_dic.keys(): 

				spec1=GridSpec(12,size_dic*2).new_subplotspec((m,ln_dic[mod]),colspan=2,rowspan=2)
                		spec2=GridSpec(12,size_dic*2).new_subplotspec((m+2,ln_dic[mod]),colspan=2)
				ax=fig.add_subplot(spec1)
                		bx=fig.add_subplot(spec2)
				bx.tick_params(which='both',direction='in',top='on',right='on',labelsize=10)
				ax.tick_params(which='both',direction='in',right='on',labelsize=10)
				bx.axhline(0.,ls='dashed',color='black',lw=1.5)

				list_table=[]
				if mod!='image':
					list_table_par1=[]
					list_table_par2=[]
					list_table_par3=[]
					col_labels=[]
					for k in filters:
						data_mod[ids[i]][k][mod]=np.genfromtxt(os.path.join(data_dir,mod+'_model_object_'+k+'_'+ids[i]+'.dat'),comments='#')
						sma_obs=data_obs[ids[i]][k][:,0]
						sma_pix=data_obs[ids[i]][k][:,0]/scale[k]
						if len(ids)==1:
							rfit_tmp=rfit[mod][k]
						else:
							rfit_tmp=rfit[mod][k][ids==ids[i]][0]
						mask_rfit=sma_pix<=rfit_tmp
						intens_obs=data_obs[ids[i]][k][:,1]
						int_err=data_obs[ids[i]][k][:,2]
						npix=data_obs[ids[i]][k][:,3]
						err=1.087*np.sqrt((sky_rms[ids[i]][k]/npix)**2+(int_err)**2)/intens_obs
						len_obs=len(sma_obs)
						fact=1./(arc_pix[k]**2)
						mag_obs=-2.5*np.log10(intens_obs*fact)+np.float(mag_zero[k][0])
						mag_mod=-2.5*np.log10(data_mod[ids[i]][k][mod][:len_obs]*fact)+np.float(mag_zero[k][0])
                                                res=mag_obs[:len(mag_mod)]-mag_mod
						ax.plot(sma_pix[~mask_rfit],mag_obs[~mask_rfit],color=filter_dic[k],ls='dashed')
						ax.errorbar(sma_pix[mask_rfit],mag_obs[mask_rfit],err[mask_rfit],marker='o',color=filter_dic[k],ls='none',mfc='white',ms=4,label=k)
						handles_leg, labels_leg = ax.get_legend_handles_labels()	
                                                ax.plot(sma_pix[:len(mag_mod)],mag_mod,color=filter_dic[k])
						ax.axvline(rfit_tmp,color=filter_dic[k],ls='dashed')
						ax.text(1.2,maxmag-0.2,models_lab[mod],fontsize=8)
						ax.yaxis.set_minor_locator(ticker.MultipleLocator(.5))
						ax.yaxis.set_major_locator(ticker.MultipleLocator(2))
						ax2=ax.twiny()
						ax2.tick_params(which='both',direction='in',labelsize=10)
						ax2.plot(sma_obs[:len(mag_mod)],mag_obs[:len(mag_mod)],color='none')
                                                #err_sub=err[:len(mag_mod)]
						#bx.errorbar(sma_pix[:len(mag_mod)][mask_rfit],res[mask_rfit]/4,err[:len(mag_mod)][mask_rfit],marker='o',color=filter_dic[k],ls='-',mfc='white',ms=4)
                                                mask_rfit2=sma_pix[:len(mag_mod)]<=rfit_tmp
                                                bx.plot(sma_pix[:len(mag_mod)],mag_obs[:len(mag_mod)],color='none')
                                                
						bx.errorbar(sma_pix[:len(mag_mod)][mask_rfit2],res[mask_rfit2]/4,err[:len(mag_mod)][mask_rfit2],marker='o',color=filter_dic[k],ls='-',mfc='white',ms=4)
						#bx.errorbar(sma_pix[:len(mag_mod)][~mask_rfit2],res[~mask_rfit2]/4,err[:len(mag_mod)][~mask_rfit2],color='none')
						bx.plot(sma_pix[:len(mag_mod)][~mask_rfit2],res[~mask_rfit2],marker='o',color='white',ms=1)
						bx.set_ylim(-0.25,0.25)
						bx.yaxis.set_minor_locator(ticker.MultipleLocator(.05))
						bx.yaxis.set_major_locator(ticker.MultipleLocator(.1))
						if ln_dic[mod]!=0:
							bx.set_yticks([])
                        				ax.set_yticks([])
							
						if ln_dic[mod]==0:
							ax.set_ylabel(r'$\mu$ (mag/arcsec$^2$)',size=12)
                                			bx.set_ylabel(r'$\Delta\mu$',size=12)

						if (m<9):
							bx.set_xticks([])
						ax2.set_xscale('log')
						if (m==0):
							ax2.xaxis.set_major_formatter(FormatStrFormatter('%d'))
							ax2.set_xlabel('SMA [pc]',size=10)
						if (m>0):
                                        		ax2.set_xticks([])
						if ((m==9) | (i==len(ids)-1)):
							bx.set_xlabel('SMA [pix]',size=10)
						if len(ids)<=1:
							list_table_par1.append(list_pars[mod][k][3][:st])
							list_table_par2.append(list_pars[mod][k][6][:st])
							list_table_par3.append(list_pars[mod][k][1][:st])
						else:
							list_table_par1.append(list_pars[mod][k][:,3][i][:st])
							list_table_par2.append(list_pars[mod][k][:,6][i][:st])
							list_table_par3.append(list_pars[mod][k][:,1][i][:st])
						col_labels.append(k)
					list_table.append(list_table_par1)
					list_table.append(list_table_par2)
					list_table.append(list_table_par3)
					row_labels=[row_labels_dic[mod][0],row_labels_dic[mod][1],r'$\chi^2$']
					the_table = ax.table(cellText=list_table,
                                                          colWidths = [0.2]*4,
                                                          rowLabels=row_labels,
                                                          colLabels=col_labels,loc='left',rowLoc='right',colLoc='right',
                                                          #bbox=Bbox([[1, 1], [3, 7]])
						          #bbox=[0.15,0.05, .2,.3])
						          #bbox=[.1,0.02, .5,.3])
						          bbox=[.1,0.02, .3,.3])
					the_table.auto_set_font_size(False)
                                        the_table.set_fontsize(4.5)
                                        #the_table.scale(1.4, 1.4)
					ax.set_xscale('log')
                                        ax.set_xticks([])
					#ax2.xaxis.set_major_formatter(FormatStrFormatter('%d'))
					bx.set_xscale('log')
                                        #ax.set_xticks([])
					bx.xaxis.set_major_formatter(FormatStrFormatter('%d'))
					ax.set_ylim(minmag+1,maxmag-1)
					#ax.invert_yaxis()		
					
					

				if mod=='image':
                        		img=Image.open(os.path.join(plots_dir,'mosaic_'+ids[i]+'.png'))
					size=int(img.size[0]/2.-1)
					circle1=plt.Circle((size,size),rfit_tmp/rmax*450.,color=filter_dic[k],fill=False,linewidth=2)
					ax.imshow(img,cmap='magma')
					ax.add_artist(circle1)
					bx.set_yticks([])
                                        ax.set_yticks([])
					bx.set_xticks([])
                                        ax.set_xticks([])
					bx.text(0.3, 0.85,'Object '+str(ids[i]),horizontalalignment='center',verticalalignment='center',transform = bx.transAxes,size=8)
					bx.legend(handles_leg,labels_leg,loc='lower left',handletextpad=0.0,fontsize=7)

				if ((m==9 and mod==ln_dic.keys()[-1]) or (i==len(ids)-1 and mod==ln_dic.keys()[-1])):
					pdf_pages.savefig(fig)
				
			m=m+3
		pdf_pages.close()
