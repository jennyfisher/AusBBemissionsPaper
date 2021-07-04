from netCDF4 import MFDataset
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap,cm
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as colors
from datetime import datetime,timedelta
import math
import numpy as np
from pyhdf.SD import SD, SDC

####################################
####################################
#######  FILL IN SETTINGS  #########
####################################
####################################

g = 9.8                            # m/simport matplotlib.colors as colors
NAv = 6.0221415e+23		   # Avogadro number
MWair = 28.9644                    # g/mole
psurf = 1013.25
H = (8.314*240)/(28.9644E-03*9.8)
pcol_const = (NAv* 10)/(MWair*g)     # scaling factor for turning vmr into pcol (using hPa levels)


Station = 'darwin'
Tracer = 'co'				
Inventory = 'GFED'
path = 'vaoab'
if Tracer == 'co':
MWtracer = 28.01
field = 'field34010'



#    		       Real	   	   ACCESS    
#                   lon      lat         J	I
#Wollongong      : 150.88   -34.41       44	176   
#Darwin          : 130.89   -12.43       62	166 
#Cape Grim 	 : 144.7    -40.7        40	173
#Cape Fergusson	 : 147.06   -19.28       57	174


for Year in ['2008,'2009','2010']:
	for Station in ['wollongong', 'Darwin', 'capegrim', 'capefergusson']:
		print('Processing '+Year+' '+Tracer+' '+Inventory+' at '+Station)

		####################################
		####################################
		###  SETTING UP OTHER VARIABLES  ###
		####################################
		####################################

		TCCON = False
		NDACC = False

		if Station == 'wollongong':
			NDACC = True
			I = 176
			J = 44
		if Station == 'darwin':
			TCCON = True
			I = 166
			J = 62
		if Station == 'capegrim':
			I = 173
			J = 40
		if Station == 'capefergusson':
			I = 174
			J = 57



		if NDACC:
			print("processing NDACC data")

			#Select hdf files
			dataset3 = SD('/home/574/mjd574/NDACC/'+Station+'/'+Tracer+'/fts_'+Tracer+'_'+Station[:4]+'_'+Year+'.hdf', SDC.READ)
			nd_times = dataset3.select('DATETIME')[:]
	

			#Convert nd time
			d0=datetime(2000,1,1,0)
			nd_datetime = []
			for time in nd_times:
				hrs=timedelta(hours=int(time*24))
				nd_datetime.append(d0+hrs)

			index = []
			index.append(np.where(np.array(nd_datetime) >= datetime(int(Year), 1, 1, 1))[0][0])
			index.append(np.where(np.array(nd_datetime) < datetime(int(Year)+1, 1, 1, 0))[-1][-1]+1)

			nd_datetimes = nd_datetime[index[0]:index[1]]


			nd_lvl = dataset3.select('ALTITUDE')[:]*1000
			nd_edges_temp = dataset3.select('ALTITUDE.BOUNDARIES')[:]*1000
			nd_pres= dataset3.select('PRESSURE_INDEPENDENT')[index[0]:index[1],:]
			nd_apr= dataset3.select(Tracer.upper()+'.COLUMN.PARTIAL_ABSORPTION.SOLAR_APRIORI')[index[0]:index[1],:]
			nd_avk= dataset3.select(Tracer.upper()+'.COLUMN_ABSORPTION.SOLAR_AVK')[index[0]:index[1],:]
			nd_pcol= dataset3.select(Tracer.upper()+'.COLUMN.PARTIAL_ABSORPTION.SOLAR')[index[0]:index[1],:]
			nd_xp = dataset3.select(Tracer.upper()+'.COLUMN_ABSORPTION.SOLAR')[index[0]:index[1]]

			#Create list of box's edges heights	
			nd_edges_temp = nd_edges_temp[0,:]
			nd_edges =[]
			for i  in np.arange(0,len(nd_lvl),1):
				nd_edges.append(nd_edges_temp[i])
			nd_edges.append(29.984505)

			#Create list of box's edge's pressures
			nd_edgePress =[]
			for i in np.arange(0,len(nd_lvl)+1,1):
				nd_edgePress.append(psurf*math.exp(-(nd_edges[i])/H))

			#Create list of box's center's pressures
			nd_edgePressdif=[]
			for i in np.arange(0,len(nd_lvl),1):
				nd_edgePressdif.append(nd_edgePress[i+1]-nd_edgePress[i])
	
			#Extract index for begining of each month + Create monthly averages of total columna, aprioris and avk
			index = []
			for n in np.arange(1,13,1):
				mytime = datetime(int(Year), n, 1, 0)
				index.append(np.where(np.array(nd_datetimes) > mytime)[0][0])
			index.append(len(nd_datetimes)-1)

			avg_nd_xp = []
			avg_apr = []
			avg_avk = []
			for n in np.arange(0,12,1):
				avg_apr.append(np.mean(nd_apr[index[n]:index[n+1],:], axis=(0)))
				avg_avk.append(np.mean(nd_avk[index[n]:index[n+1],:], axis=(0)))
				avg_nd_xp.append(np.mean(nd_xp[index[n]:index[n+1]], axis=(0)))
			avg_nd_xp = np.nan_to_num(np.array(avg_nd_xp))
			avg_apr = np.nan_to_num(np.array(avg_apr))
			avg_avk = np.nan_to_num(np.array(avg_avk))
		
		if TCCON:
			print("processing TCCON data")
			#select the file
			dataset3 = Dataset('Darwin_TCCON.nc')
			nd_times = dataset3.variables['time'][:]
			nd_aprtimes = dataset3.variables['prior_date'][:]

			#Convert nd time
			d0=datetime(1970,1,1,0)
			nd_datetimes = []
			nd_aprdatetimes = []
			for time in nd_times:
				hrs=timedelta(hours=int(time*24))
				nd_datetimes.append(d0+hrs)

			for time in nd_aprtimes:
				hrs=timedelta(hours=int(time*24))
				nd_aprdatetimes.append(d0+hrs)		

			#Find index limits of Year selected
			index = []
			index2 =[]
			startyear = datetime(int(Year),1,1,1)
			startyearb = datetime(int(Year),1,1,0)
			endyear = datetime(int(Year)+1,1,1,0)
			index.append(np.where(np.array(nd_datetimes) >= startyear)[0][0])	
			index.append(np.where(np.array(nd_datetimes) < endyear)[-1][-1] +1)	
			index2.append(np.where(np.array(nd_aprdatetimes) >= startyearb)[0][0])	
			index2.append(np.where(np.array(nd_aprdatetimes) < endyear)[-1][-1] +1)

			#Load variables for Year selected
			nd_datetimes = nd_datetimes[index[0]:index[1]]
			nd_aprdatetimes = nd_aprdatetimes[index2[0]:index2[1]]
			nd_vmr = dataset3.variables['xco_ppb'][index[0]:index[1]]
			nd_avk = dataset3.variables['ak_co'][::-1,:]
			nd_apr = dataset3.variables['prior_co'][index2[0]:index2[1],::-1]
			nd_lvl = dataset3.variables['prior_Height'][::-1]*1000+30

			nd_surf_p = dataset3.variables['pout_hPa'][index[0]:index[1]]
			nd_zen = dataset3.variables['asza_deg'][index[0]:index[1]]
			nd_zen_ref = dataset3.variables['ak_zenith'][:]
			nd_xp = nd_vmr * nd_surf_p * 2.12e13

			#Create list of box's edge's height
			nd_edges = []
			nd_edges.append([0])
			nd_edges.append(list(np.arange(500,71500,1000)))
			nd_edges = sum(nd_edges,[])[::-1]

			#Create list of box's edge's pressures
			nd_edgePress =[]
			for i in np.arange(0,len(nd_lvl)+1,1):
				nd_edgePress.append(psurf*math.exp(-(nd_edges[i])/H))

			#Create list of box's center's pressures
			nd_edgePressdif=[]
			for i in np.arange(0,len(nd_lvl),1):
				nd_edgePressdif.append(nd_edgePress[i+1]-nd_edgePress[i])

			#Create arrays of apriori and averaging kernel values for each time stamp
			nd_apr_temp=[]
			nd_avk_temp=[]
			for a in np.arange(0,len(nd_datetimes),1):
				b = (np.where(np.array(nd_aprdatetimes)<= nd_datetimes[a])[-1][-1])
				nd_apr_temp.append(nd_apr[b,:])
				c = np.where(nd_zen_ref == min(nd_zen_ref, key=lambda x:abs(x-nd_zen[a])))[0][0]
				nd_avk_temp.append(nd_avk[:,c])

			nd_apr = np.array(nd_apr_temp)
			nd_avk = np.array(nd_avk_temp)

			#Extract index for begining of each month + Create monthly averages of total columna, aprioris and avk
			index = []
			for n in np.arange(1,13,1):
				mytime = datetime(int(Year), n, 1, 0)
				index.append(np.where(np.array(nd_datetimes) > mytime)[0][0])

			index.append(len(nd_datetimes)-1)

			avg_nd_xp = []
			avg_apr = []
			avg_avk = []
			for n in np.arange(0,12,1):
				avg_apr.append(np.mean(nd_apr[index[n]:index[n+1],:], axis=(0)))
				avg_avk.append(np.mean(nd_avk[index[n]:index[n+1],:], axis=(0)))
				avg_nd_xp.append(np.mean(nd_xp[index[n]:index[n+1]], axis=(0)))

			avg_nd_xp =np.nan_to_num(np.array(avg_nd_xp))
			avg_apr =np.nan_to_num(np.array(avg_apr))
			avg_avk =np.nan_to_num(np.array(avg_avk))


		####################################
		####################################
		#######     ACCESS-UKCA     ########
		####################################
		####################################

		print("Processing ACCESS data")

		#Select ACCESS Files(s)
		dataset2 = MFDataset("/g/data1/p66/akl599/ACCESS/output/"+path+"/"+path+"a.pm-"+Year+"-*.nc")

		#extract ACCESS data
		ac_times=dataset2.variables['time'][:]
		ac_lvl = dataset2.variables['z1_hybrid_height'][:]
		ac_lon = np.arange(-180.9375,180.9375,1.875)
		ac_lat = np.arange(-90.625,90.625,1.25)

		#[time,lev,lat,lon_1]
		ac_vmr_temp=(dataset2.variables[field][:])* 1E9 * 28.9644 / MWtracer
		ac_vmr=np.concatenate((ac_vmr_temp[:,:,:,96:192],ac_vmr_temp[:,:,:,0:96]),axis=3)
		ac_vmr_surf=ac_vmr[:,0,:,:]  
		ac_vmr_surf_avg = np.mean(ac_vmr_surf, axis=(1,2))
		ac_vmr_surf_map = np.mean(ac_vmr_surf, axis=(0))
		ac_vmr_vert_avg = np.mean(ac_vmr, axis=(0,2,3))

		#convert ac time
		d0=datetime(1970,1,1,0)
		ac_datetimes = []
		for time in ac_times:
			hrs=timedelta(hours=int(time*24))
			ac_datetimes.append(d0+hrs)

		surf_stations=[]
		surf_stations.append(ac_datetimes)
		surf_stations.append(ac_vmr_surf[:,44,176])		# Extract surface vmr for wollongong
		surf_stations.append(ac_vmr_surf[:,62,166])		# Extract surface vmr for darwin
		surf_stations.append(ac_vmr_surf[:,39,173])		# Extract surface vmr for cape grim
		surf_stations.append(ac_vmr_surf[:,57,174])		# Extract surface vmr for capeferguson
		surf_stations=np.transpose(np.array(surf_stations))
		np.savetxt("Monthly_ACCESS_Stations_Surface_"+Tracer+"_"+Inventory+"_"+Year+"_"+Tracer.upper()+".csv",surf_stations, fmt="%s", delimiter=",", header="Month,ACCESS_"+Inventory+"_"+Tracer+"_Wollongong,ACCESS_"+Inventory+"_"+Tracer+"_Darwin,ACCESS_"+Inventory+"_"+Tracer+"_Cape_Grim,ACCESS_"+Inventory+"_"+Tracer+"_Cape_Ferguson")

		if NDACC or TCCON:
			station_ac_vmr = ac_vmr[:,:,J,I]
			rev_nd_lvl = nd_lvl[::-1]
			interp_station_ac_vmr = []

			mm = 12
	
			for i in np.arange(0,mm,1):
				interp_station_ac_vmr.append(np.interp(rev_nd_lvl[:],ac_lvl[:],station_ac_vmr[i,:]))
			interp_station_ac_vmr =np.array(interp_station_ac_vmr)[:,::-1]

			ac_pcol =[]
			for i in np.arange(0,mm,1):
				ac_pcol.append(2.12e13 * interp_station_ac_vmr[i,:] * nd_edgePressdif[:])
			ac_pcol =np.array(ac_pcol)

			ac_ap_pcol =[]
			for i in np.arange(0,mm,1):
				ac_ap_pcol.append(avg_apr[i,:] + avg_avk[i,:]*(ac_pcol[i,:]-avg_apr[i,:]))
			ac_ap_pcol =np.array(ac_ap_pcol)

			ac_tcol=np.sum(ac_ap_pcol, axis=(1))

			toprint = []
			toprint.append(ac_datetimes)
			toprint.append(ac_tcol)
			toprint.append(avg_nd_xp)
			toprint=np.array(toprint)	
			toprint=np.transpose(np.array(toprint))
			#np.savetxt("Monthly_"+Station+"_"+Year+"_"+Inventory+"_"+Tracer+"_ACCESS_vs_FTS.csv",toprint, delimiter=",", fmt="%s", header="Date_Time,"+Station+"_"+Inventory+"_"+Tracer+"_tcol,FTS_"+Station+"_"+Tracer+"_tcol")  

