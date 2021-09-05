from netCDF4 import MFDataset
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap,cm
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as colors
from datetime import datetime,timedelta
import math
import numpy as np
import sys
from pyhdf.SD import SD, SDC

####################################
####################################
#######  FILL IN SETTINGS  #########
####################################
####################################

NDAC_folder = '/NDACC/wollongong/CO/'
TCCON_file= '/TCCON/Darwin_TCCON.nc'

# Assuming GEOS_Chem files are in e.g. 'GC/rundirs/geos5_2x25_tropchem_GFED/2008/GFED_2008.nc' folders/files
# and station files are in e.g. 'GC/rundirs/geos5_2x25_tropchem_GFED/2008/Wollongong.nc' folder/files

GEOS_Chem_folder = '/GC/rundirs/'
			

g = 9.8                            
NAv = 6.0221415e+23		   # Avogadro number
MWair = 28.9644                    # g/mole
psurf = 1013.25
H = (8.314*240)/(28.9644E-03*9.8)
pcol_const = (NAv* 10)/(MWair*g)     # scaling factor for turning vmr into pcol (using hPa levels)

Tracer = 'co'
MWtracer = 28.01
gc_field = 'IJ_AVG_S__CO'
gc_field2 = 'BIOBSRCE__CO'

# 2x25		       Real	   	IDL/Phyton      
#                   lon      lat          I     J	
#Wollongong      : 150.88   -34.41        132   28    	  
#Darwin          : 130.89   -12.43        124   39    
#Cape Grim 	 : 144.7    -40.7         130   25      
#Cape Fergusson	 : 147.06   -19.28        131   35	




for Year in ['2008','2209','2010']:
	for Station in ['wollongong', 'Darwin', 'capegrim', 'capefergusson']:				
		for Inventory in ['FINN','GFED','QFED']:		   	

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
				I = 132
				J = 28
			if Station == 'darwin':
				TCCON = True
				I = 124
				J = 39
			if Station == 'capegrim':
				I = 130
				J = 25
				J = 35
			if Station == 'capefergusson':
				I = 131
				J = 35

	
			####################################
			####################################
			###########      FTS      ##########
			####################################
			####################################

			if NDACC:
				print("processing NDACC data")

				#Select hdf files
				dataset3 = SD(NDAC_folder+'fts_CO_woll_'+Year+'.hdf', SDC.READ)
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
				for i in np.arange(0,len(nd_lvl),1):
					nd_edges.append(nd_edges_temp[i])
				nd_edges.append(29.984505)

				nd_BXHGT=[]
				for i in np.arange(0,len(nd_lvl),1):
					nd_BXHGT.append(nd_edges[i]-nd_edges[i+1])

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
				dataset3 = Dataset(TCCON_file)
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
			#######      GEOS-Chem      ########
			####################################
			####################################

			#[time,lev,lat,lon]
			print("Processing GEOS-Chem data")
			#select GC File(s)
			dataset1 = MFDataset(GEOS_Chem_folder+'geos5_2x25_tropchem_'+Inventory+'/'+Year+'/'+Inventory+'_'+Year+'.nc')
			dataset4 = MFDataset(GEOS_Chem_folder+'geos5_2x25_tropchem_'+Inventory+'/'+Year+'/'+Station.capitalize()+'.nc')

			####################################
			#Monthly work
			####################################

			#Extract from Monthly files
			gc_times=dataset1.variables['time'][:]			# Extract GC time (TAU format - hours since 1985
			NAIR=dataset1.variables['BXHGHT_S__N_AIR_'][:] 		# Extract air density (molecules per m3)
			BXHT=dataset1.variables['BXHGHT_S__BXHEIGHT'][:] 	# Extract grid box height (metres)
			gc_area=dataset1.variables['DXYP__DXYP'][:]		# Extract grid column surface area
			gc_vmr=dataset1.variables[gc_field][:]			# Extract volume mixing ratios
			
			gc_bbemission=dataset1.variables[gc_field2][:]	# Extract fire emission (C atoms per m^2 per second)
			

			#Create long / lat grid
			gc_lon = np.arange(-181.25,181.25,2.5)			# Create grid edges for longitudes (res 2.5 deg)
			gc_lat = np.arange(-91,93,2)				# Create grid edges for latitude (res 2 deg)

			#convert gc time
			d0=datetime(1985,1,1,0)					# Tau 0 = 1/1/1985
			gc_datetimes = []					# Create list for new time array
			for time in gc_times:
				hrs=timedelta(hours=int(time))			# DateTime package extracts 'hours' since 1985 date format
				gc_datetimes.append(d0+hrs)			# Fill list with formated dates

			#Create box centres heights
			BXHT = np.mean(BXHT, axis=(0,2,3))			# Remove the time, long and lat dimensions
			gc_lvl = []						# Create list for level altitude
			gc_lvl.append((BXHT[0])/2)				# Fill first level
			for i in range(1,47):
				 gc_lvl.append((gc_lvl[i-1]+BXHT[i-1]/2)+(BXHT[i])/2)	# Fill in all other levels

			gc_vmr_surf = gc_vmr[:,0,:,:]				# Extract surface vmr
			gc_vmr_surf_avg = np.mean(gc_vmr_surf, axis=(1,2))	# Average to global surface vmr
			gc_vmr_surf_map = np.mean(gc_vmr_surf, axis=(0))	# Average to surface vmr / year
			gc_vmr_vert_avg = np.mean(gc_vmr, axis=(0,2,3))		# Average to global and yearly vertical profile

			NAIR = np.mean(NAIR, axis=(0,2,3))			# Remove the time, long and lat dimensions
			gc_molec_vert_avg = gc_vmr_vert_avg*1E-9*NAIR		# Convert to vmr to C/

			#Compare to ref model (here GFED4s)
			dataset_ref = MFDataset(GEOS_Chem_folder+'geos5_2x25_tropchem_GFED/'+Year+'/GFED_'+Year+'.nc')
			gc_vmr_ref=dataset_ref.variables[gc_field][:]
			gc_vmr_ref_surf = gc_vmr_ref[:,0,:,:]			
			gc_vmr_ref_surf_map = np.mean(gc_vmr_ref_surf, axis=(0))
			gc_vmr_ref_surf_station = gc_vmr_ref_surf[:,J,I]
			gc_vmr_surf_dif=(gc_vmr_surf_map[:,:]-gc_vmr_ref_surf_map[:,:])

			#Select stations from Monthly files

			surf_stations=[]
			surf_stations.append(gc_datetimes)
			surf_stations.append(gc_vmr_surf[:,28,132])		# Extract surface vmr for wollongong
			surf_stations.append(gc_vmr_surf[:,39,124])		# Extract surface vmr for darwin
			surf_stations.append(gc_vmr_surf[:,25,130])		# Extract surface vmr for cape grim
			surf_stations.append(gc_vmr_surf[:,35,131])		# Extract surface vmr for capeferguson
			surf_stations=np.transpose(np.array(surf_stations))
			np.savetxt("Monthly_GC_Stations_Surface_"+Tracer+"_"+Inventory+"_"+Year+"_"+Tracer.upper()+".csv",surf_stations, fmt="%s", delimiter=",", header="Month,Wollongong,Darwin,Cape_Grim,Cape_Ferguson")

			####################################
			#Station work
			####################################
			#extract GC station data into arrays
			station_vmr=dataset4.variables[Tracer.upper()][:]	# Extract station volume molar ratio (ppb Carbon atoms)

			station_times=dataset4.variables['Tau'][:]		# Extract statiom time (TAU format - hours since 1985)
			gc_airden=dataset4.variables['AIRDEN']


			#convert gc time = Dataset('2008.nc')
			d0=datetime(1985,1,1,0)					# Tau 0 = 1/1/1985
			station_datetimes = []					# Create list for new time array
			for time in station_times:
				hrs=timedelta(hours=int(time))			# DateTime package extracts 'hours' since 1985 date format
				station_datetimes.append(d0+hrs)		# Fill list with formated dates

			station_vmr_surf = station_vmr[:,0]
			surf_station=[]
			surf_station.append(station_datetimes)
			surf_station.append(station_vmr_surf)
			surf_station=np.transpose(np.array(surf_station))
			np.savetxt("GC_"+Station+"_Surface_"+Tracer+"_"+Inventory+"_"+Year+"_"+Tracer.upper()+".csv",surf_station, fmt="%s", delimiter=",", header="DateTime,"+Station+"_"+Tracer)

			if NDACC or TCCON:
				#Find matching FTS/GC time stamps and extrac GC data
				station_matched_vmr = []
				station_matched_airden = []
				for i in np.arange(0,len(nd_datetimes),1):
					for j in np.arange(0,len(station_datetimes),1):
						if nd_datetimes[i] == station_datetimes[j]:
							station_matched_vmr.append(station_vmr[j,:])
							station_matched_airden.append(gc_airden[j,:])
							break
				station_matched_vmr = np.array(station_matched_vmr)
				station_matched_airden = np.array(station_matched_airden)
				
				special_print=[]
				special_print.append(nd_datetimes[0])
				special_print.append(gc_lvl)
				special_print.append(station_matched_vmr[0,:])
				special_print.append(nd_lvl)
	
				#Interpolate GC data to FTS levels
				rev_nd_lvl = nd_lvl[::-1]
				interp_station_vmr = []
				for i in np.arange(0,len(nd_datetimes),1):
					interp_station_vmr.append(np.interp(rev_nd_lvl,gc_lvl,station_matched_vmr[i,:]))
				interp_station_vmr=np.array(interp_station_vmr)[:,::-1]

				special_print.append(interp_station_vmr[0,:])


				monthly_station_vmr=gc_vmr[:,:,J,I]
				interp_monthly_station_vmr=[]
				for i in np.arange(0,12,1):
					interp_monthly_station_vmr.append(np.interp(rev_nd_lvl,gc_lvl,monthly_station_vmr[i,:]))
				interp_monthly_station_vmr =np.array(interp_monthly_station_vmr)[:,::-1]

				#Calculate interpolated GC partial columns
				station_pcol =[]
				for i in np.arange(0,len(nd_datetimes),1):
					station_pcol.append(2.12e13 * interp_station_vmr[i,:] * nd_edgePressdif[:])
				station_pcol =np.array(station_pcol)						

				special_print.append(station_pcol[0,:])

				monthly_station_pcol =[]
				for i in np.arange(0,12,1):
					monthly_station_pcol.append(2.12e13 * interp_monthly_station_vmr[i,:] * nd_edgePressdif[:])
				monthly_station_pcol =np.array(monthly_station_pcol)					

				#Apply averaging kernels and aprioris on GC partial columns
				station_ap_pcol =[]
				for i in np.arange(0,len(nd_datetimes),1):
					station_ap_pcol.append(nd_apr[i,:] + nd_avk[i,:]*(station_pcol[i,:]-nd_apr[i,:]))
				station_ap_pcol =np.array(station_ap_pcol)					

				special_print.append(nd_apr[0,:])
				special_print.append(nd_avk[0,:])
				special_print.append(station_ap_pcol[0,:])
				special_print=np.transpose(special_print)
				np.savetxt("Special.csv",special_print, delimiter=",", fmt="%s", header="Date_Time,GC_levels,GC_vmr,ND_levels,GC_interp_vmr,GC_pcol,ND_apr,nd_avk,GC_smooth_pcol")

				monthly_station_ap_pcol =[]
				for i in np.arange(0,12,1):
					monthly_station_ap_pcol.append(avg_apr[i,:] + avg_avk[i,:]*(monthly_station_pcol[i,:]-avg_apr[i,:]))
				monthly_station_ap_pcol =np.array(monthly_station_ap_pcol)				

				#Sum partial columns into total columns
				station_tcol=np.sum(station_ap_pcol, axis=(1))
				monthly_station_tcol=np.sum(monthly_station_ap_pcol, axis=(1))	

				toprint = []
				toprint.append(nd_datetimes)
				toprint.append(station_tcol)
				toprint.append(nd_xp)
				toprint=np.transpose(np.array(toprint))
				np.savetxt(Station+"_"+Year+"_"+Inventory+"_"+Tracer+"_GC_vs_FTS.csv",toprint, delimiter=",", fmt="%s", header="Date_Time,"+Station+"_"+Inventory+"_GC_"+Tracer+"_tcol,FTS_"+Tracer+'_'+Station+"_tcol")	

				toprint2 = []
				toprint2.append(gc_datetimes)
				toprint2.append(monthly_station_tcol)
				toprint2.append(avg_nd_xp)
				toprint2=np.transpose(np.array(toprint2))
				np.savetxt("Monthly_"+Station+"_"+Year+"_"+Inventory+"_"+Tracer+"_GC_vs_FTS.csv",toprint2, delimiter=",", fmt="%s", header="Date_Time,"+Station+"_"+Inventory+"_GC_"+Tracer+"_tcol,FTS_"+Tracer+'_'+Station+"_tcol")	
	
