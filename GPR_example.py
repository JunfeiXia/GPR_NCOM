#import .mat file to format for GPy
#Junfei Xia edited based on John Lodise's code
#6-20-2010

import h5py 
import math 
import mat4py
from multiprocessing import Process, Queue
import multiprocessing as mp
from scipy import interpolate
import scipy
import io
from datetime import datetime,timedelta
# from netCDF4 import Dataset
import pyproj
import GPy

from GPy.kern import Kern
from GPy.core import Param
#import myKernel2D

import numpy as np
from matplotlib import pylab as pl
import matplotlib.pyplot as plt
from matplotlib import rc,rcParams
import matplotlib as mpl 

# lat0 = 36.
# lon0 = -4.
# NAD83=pyproj.Proj("+init=EPSG:2062",preserve_units = False) #Spain

################################Set up Coordinate system #################################
# lat0 = 28.2
# lon0 = -88.35
lat0 = 18.049999237060547
lon0 = -97.950012207031250
# NAD83=pyproj.Proj("EPSG:3453",preserve_units = False)#False = meters #Louisiana South (ftUS)
NAD83=pyproj.Proj('+proj=utm +zone=10 +ellps=WGS84',preserve_units = False)#False = meters #Louisiana South (ftUS)
x_ori,y_ori=NAD83(lon0,lat0) #define x,y origin using lat0,lon0
x_ori = x_ori
y_ori = y_ori
print(x_ori)

rearth = 6373.3 * 1000   # at lat0 = 29., in meters 


# x1,y1 = NAD83(-83,29)
# x2,y2 = NAD83(-83,28)
x1,y1= NAD83( -88.669711135515044, 28.454682423388576)
print(x1-x_ori,y1-x_ori)
# d = np.sqrt((x1-x2)**2 + (y1-y2)**2) # check to make sure we're using meters, should be roughly around 111,000m
# print(d)
###########################################################################

##########################################################################
#Import file with desired drifter velocity data

def getData(fname = 'Test_drifter_data.mat'):


    f = h5py.File(fname,'r')
    lon = f.get('DD_lon_select')
    lont = np.array(lon) 
	
    lat = f.get('DD_lat_select')
    latt = np.array(lat) 

    u_sel = f.get('u_select')
    uob = np.array(u_sel)

    v_sel = f.get('v_select')
    vob = np.array(v_sel)
	
    year = f.get('year')
    yr = np.array(year).astype(int) 
	
    mo = f.get('mo')
    mon = np.array(mo).astype(int)
	
    da = f.get('da')
    day= np.array(da).astype(int)
	
    hour = f.get('hr')
    hr= np.array(hour).astype(int)
	
    minute = f.get('mn')
    mn= np.array(minute).astype(int)
	
    time = f.get('py_time')
    tob = np.array(time)
	
    xob2_0,yob2_0=NAD83(lont,latt)
    
    xob2 = xob2_0
    yob2 = yob2_0
    
    
    yob2[np.where(np.isnan(vob))]=np.nan # avoid missing values 
    xob2[np.where(np.isnan(uob))]=np.nan
    yob2[np.where(np.isnan(latt))]=np.nan
    xob2[np.where(np.isnan(lont))]=np.nan
    xob = (xob2 - x_ori)/1000. # in km  # convert to km coordinate system with origin at lat0,lon0
    yob = (yob2 - y_ori)/1000. # in km
    
    #####define again for plotting purposes#######
    xob_plot = (xob2 - x_ori)/1000. # in km 
    yob_plot = (yob2 - y_ori)/1000. # in km
    tob_plot = tob[:,0]	
    ##############################

	
   # date_st = datetime(np.min(yr),np.min(mon),np.min(day),np.min(hr),np.min(mn))
    #date_end = datetime(np.max(yr),np.max(mon),np.max(day),np.max(hr),np.max(mn))

	
	
    return tob,yob,xob,latt,lont,vob,uob,yr,mon,day,hr,mn,xob_plot,yob_plot,tob_plot,fname
	

	
	##########################################################################################
def obsMatrix(tob,yob,xob,latt,lont,vob,uob):
# 	'''
#     Organize drifter data in vectors ready to use with GPy.
# 
#     Inputs:
#             tob: n by m matrix, with n position times for the m drifters
#             yob: n by m matrix with meridional coordinates of the drifters positions in km
#             xob: n by m matrix with zonal coordinates of the drifters positions in km
#             latt: n by m matrix with the latitude of the drifters positions
#             lont: n by m matrix with the longitude of the drifters positions
#             vob: n by m matrix with meridional velocity of the drifters in m/s
#             uob: n by m matrix with zonal velocity of the drifters in m/s
# 
#     Outputs:
#             X: (n x m) by 3 matrix, with time (hours), and with merdional and zonal positions (km)
# 	for each drifter position. 
#             LL: (n x m) by 3 matrix, with time (hours), and with latitude and longitude 
# 	for each drifter position. 
#             obs: (n x m) by 2 matrix, with merdional and zonal velocity components (m/s) 
#     '''
	tob = np.reshape(tob,[-1,1])
	yob = np.reshape(yob,[-1,1])
	xob = np.reshape(xob,[-1,1])
	latt = np.reshape(latt,[-1,1])
	lont = np.reshape(lont,[-1,1])
	uob = np.reshape(uob,[-1,1])
	vob = np.reshape(vob,[-1,1])
	validPoints = np.where( (~np.isnan(xob)) & (~np.isnan(yob)) & (~np.isnan(uob)) & (~np.isnan(vob)) )
	tob = tob[validPoints][:,None]
	xob = xob[validPoints][:,None]
	yob = yob[validPoints][:,None]
	lont = lont[validPoints][:,None]
	latt = latt[validPoints][:,None]
	uob = uob[validPoints][:,None]
	vob = vob[validPoints][:,None]

	# From here on, always use T,Y,X order
	X = np.concatenate([tob,yob,xob],axis=1)
	LL = np.concatenate([tob,latt,lont],axis=1)
	obs = np.concatenate([vob,uob],axis=1)
    
    

	return X, LL, obs
	#################################################################################################
	# Assuming that your X matrix 
	# X[:,0] has the observations times
	# X[:,1] has the meridional positions of the observations
	# X[:,2] has the zonal positions of the observations



def vel_field(X,LL,obs,tob):

	

	xmin = np.nanmin(X[:,2]) # lower limit in the zonal direction 
	xmax = np.nanmax(X[:,2])# lower limit in the zonal direction 

	ymin = np.nanmin(X[:,1]) # lower limit in the meridional direction 
	ymax = np.nanmax(X[:,1]) # lower limit in the meridional direction 

	tmin = np.nanmin(X[:,0]) # lower limit in the time direction 
	tmax = np.nanmax(X[:,0]) # lower limit in the timel direction
	
	print(xmin,xmax)
	print(ymin,ymax)
	print(tmin,tmax)


	# DEFINE COORDINATES of your desired grid 
	dx = 1 # space resolution of the grid. If positions in X are in km, 0.5 is 500 m resolution... used to be 0.05
	xg = np.arange(xmin,xmax+dx,dx)
	yg = np.arange(ymin,ymax+dx,dx) 
	print(x_ori)
	


	dt = .25 #  hours... for plotting, doesn't affect regression  
	tg = np.arange(tmin,tmax+dt,dt)

	# Create grid with size M = size(tg) * size(yg) * size(xg)

	YY,TT,XX = np.meshgrid(yg,tg,xg)
   	#  1st index vary with T
   	#  2nd index vary with Y
   	#  3rd index vary with X
	# reshape each dimension as a M by 1 vector
	TT = np.reshape(TT,[TT.size,1])
	YY = np.reshape(YY,[YY.size,1]) 
	XX = np.reshape(XX,[XX.size,1])

	# Merge dimensions, creating the (M by 3) Xg matrix to be used with GPy
	Xg =  np.concatenate([TT,YY,XX],axis=1)

	
	tc = tob[:,0]
	print(tc)
	
	
	 
	

	##############################################################
	# Now to use GPy, you only have to execute 3 lines for each velocity component
	# get covariance function for GPR
	# here I chose the squared exponential (RBF), but there are other options 
	k1 = GPy.kern.RBF(input_dim=3,ARD=True)     # covariance function
	# Other types of covariance function you could use:
	#k2 = GPy.kern.Matern32(input_dim=3,ARD=True)   
	#k3 = GPy.kern.Matern52(input_dim=3,ARD=True)   
	#k4 = GPy.kern.RatQuad(input_dim=3,ARD=True)   
	
	
	nres = 5 # restart optimization procedure 5 times
	
	k_2_Hyp = k1.copy() + k1.copy()  #creating sum of two squared exponentials as covariance function 

	model_u = GPy.models.GPRegression(X,obs[:,1][:,None],k_2_Hyp.copy())  
	model_u.optimize_restarts(messages=True,num_restarts=nres) 
	print(model_u[0],model_u[1],model_u[2],model_u[3],model_u[4],model_u[5],model_u[6],model_u[7],model_u[8]) #hyperparatmeters that make up model_u
	# order of hyperparameters: [variance 1 (m^2/s^2), time scale 1 (h), meridional(y) scale 1 (km), zonal(x) scale 1 (km), variance 2, time scale 2, meridional(y) scale 2, zonal(x) scale 2, gaussian noise(m^2/s^2)]
	U0,Ku0 = model_u.predict(Xg)
	#Ku is the diagonal of the posterior covariance. The square root of Ku are error estimates for U
	
	model_v = GPy.models.GPRegression(X,obs[:,0][:,None],k_2_Hyp.copy())
	model_v.optimize_restarts(messages=True,num_restarts=nres) 
	print(model_v[0],model_v[1],model_v[2],model_v[3],model_v[4],model_v[5],model_v[6],model_v[7],model_v[8])
	V0,Kv0 = model_v.predict(Xg) 
  		
#  	parallel coding to fill in grid faster/limited memory (uses 4 processors) 
#  	step = len(Xg[:,0])/4
#  	#print step
#  	Xg_1 = Xg[0:step,:]
#  	Xg_2 = Xg[step:2*step,:]
#  	Xg_3 = Xg[2*step:3*step,:]
#  	Xg_4 = Xg[3*step:,:]
#  	
#  	manager = mp.Manager()
#  	
#  	q_u_1 = manager.Queue()
#  	q_v_1 = manager.Queue()
#  	q_ku_1 = manager.Queue()
#  	q_kv_1 = manager.Queue()
#  	
#  	q_u_2 = manager.Queue()
#  	q_v_2 = manager.Queue()
#  	q_ku_2 = manager.Queue()
#  	q_kv_2 = manager.Queue()
#  	
#  	q_u_3 = manager.Queue()
#  	q_v_3 = manager.Queue()
#  	q_ku_3 = manager.Queue()
#  	q_kv_3 = manager.Queue()
#  	
#  	q_u_4 = manager.Queue()
#  	q_v_4 = manager.Queue()
#  	q_ku_4 = manager.Queue()
#  	q_kv_4 = manager.Queue()
#  	
#  	
# 
#  	p1 = Process(target=fill_vel_grid, args=(model_u,model_v,Xg_1,q_u_1,q_v_1,q_ku_1,q_kv_1))
# 	p1.start()
# 	p2 = Process(target=fill_vel_grid, args=(model_u,model_v,Xg_2,q_u_2,q_v_2,q_ku_2,q_kv_2))
# 	p2.start()
# 	p3 = Process(target=fill_vel_grid, args=(model_u,model_v,Xg_3,q_u_3,q_v_3,q_ku_3,q_kv_3))
# 	p3.start()
# 	p4 = Process(target=fill_vel_grid, args=(model_u,model_v,Xg_4,q_u_4,q_v_4,q_ku_4,q_kv_4))
# 	p4.start()
#  	
#  	#print q_u_1.get()[:] # q_ku_2.get().shape, q_ku_3.get().shape, q_u_4.get().shape
#  	
#  	U0 = np.concatenate((q_u_1.get(),q_u_2.get(),q_u_3.get(),q_u_4.get()),axis=None)
# 	V0 = np.concatenate((q_v_1.get(),q_v_2.get(),q_v_3.get(),q_v_4.get()),axis=None)
#  	Ku = np.concatenate((q_ku_1.get(),q_ku_2.get(),q_ku_3.get(),q_ku_4.get()),axis=None)
#  	Kv = np.concatenate((q_kv_1.get(),q_kv_2.get(),q_kv_3.get(),q_kv_4.get()),axis=None)
#  	 	
#  	print  	U0.shape, Ku.shape
#  	 	
# 	p1.join()
# 	p2.join()
# 	p3.join()
# 	p4.join() 
	
	    
	U = np.reshape(U0,[tg.size,yg.size,xg.size])  # reshape U into a size(tg),size(yg),size(xg) matrix
	Ku = np.reshape(Ku0,[tg.size,yg.size,xg.size])  # reshape posterior covariance into a size(tg),size(yg),size(xg) matrix
	#Ku is the diagonal of the posterior covariance. The square root of Ku are error estimates for U
	V = np.reshape(V0,[tg.size,yg.size,xg.size])  
	Kv = np.reshape(Kv0,[tg.size,yg.size,xg.size]) 
 	
 	
 	
	lat_degs_0 =( ((Xg[:,1])*1000.) + y_ori) # converting grid back to lat lon coords
	lon_degs_0 =( ((Xg[:,2])*1000.) + x_ori) # converting grid back to lat lon coords
	
	xg0 = np.reshape(Xg[:,2],[tg.size,yg.size,xg.size])  
	yg0 = np.reshape(Xg[:,1],[tg.size,yg.size,xg.size])
	
	lon_degs1, lat_degs1 = NAD83(lon_degs_0,lat_degs_0,inverse=True)
	
	print(lon_degs1)
	print(lat_degs1)
    
	lon_degs = np.reshape(lon_degs1,[tg.size,yg.size,xg.size])  
	lat_degs= np.reshape(lat_degs1,[tg.size,yg.size,xg.size]) 
 	
	return U,Ku,V,Kv,Xg,xg0,yg0,tg,model_u,model_v,lon_degs,lat_degs



	##################Parallell code#############
	
	
	
# def fill_vel_grid(model_u,model_v,Xg_0,q_u,q_v,q_ku,q_kv):
#   	
# 	Ku = np.zeros([len(Xg_0[:,0]),1])
# 	U = np.zeros([len(Xg_0[:,0]),1])
# 	Kv = np.zeros([len(Xg_0[:,0]),1])
# 	V = np.zeros([len(Xg_0[:,0]),1])
# 	
# 	ind = len(Xg_0[:,0])
# 	print ind
# 	
# 	intervals = ind/25000 + 1
# 	
# 	
#  	for j in range(intervals):
#  		if ind > 25000:
#  			
#  			U_0,Ku_0 = model_u.predict(Xg_0[ind-25005:ind,:])   
#  			V_0,Kv_0 = model_v.predict(Xg_0[ind-25005:ind,:])
#  			U[ind-25005:ind,:] = U_0
#  			Ku[ind-25005:ind,:] = Ku_0
#  			V[ind-25005:ind,:] = V_0
#  			Kv[ind-25005:ind,:] = Kv_0
#  			
#  			
# 
#  			print ind, len(Xg_0[ind-25000:ind])
#  			ind = ind - 25000
#  
#   		if ind < 25000 and ind > 0:
#   			U_0,Ku_0 = model_u.predict(Xg_0[0:ind+1,:])
#   			V_0,Kv_0 = model_v.predict(Xg_0[0:ind+1,:])
#   			U[0:ind+1,:] = U_0
#   			Ku[0:ind+1,:] = Ku_0
#   			V[0:ind+1,:] = V_0
#   			Kv[0:ind+1,:] = Kv_0
#  		
#  			print ind, len(Xg_0[0:ind])
#  			ind = ind - 25000
# 
# 
# 	
# 	q_u.put(U)
# 	q_v.put(V)
# 	q_ku.put(Ku)
# 	q_kv.put(Kv)
	

 	
 	
def SaveMat(fname,xg0,yg0,U,V,X,Ku,Kv,tg,obs,LL,yob,xob,tob,yr,mon,day,hr,mn,uob,vob,latt,lont,xob_plot,yob_plot,tob_plot,model_u,model_v,lon_degs,lat_degs):

	data = {'xg0':xg0,'yg0':yg0,'U':U,'V':V,'X':X,'Ku':Ku,'Kv':Kv,'tg':tg,'obs':obs,'LL':LL,'yob':yob,'xob':xob,'tob':tob,'yr':yr,'mon':mon,'day':day,'hr':hr,'mn':mn,'uob':uob,'vob':vob,'latt':latt,'lont':lont,'xob_plot':xob_plot,'yob_plot':yob_plot,'tob_plot':tob_plot,'model_u':model_u.param_array,'model_v':model_v.param_array,'lon_degs':lon_degs,'lat_degs':lat_degs}

	scipy.io.savemat(fname[0:-4] + '_DD_OI.mat', data)

# # run the code
# tob,yob,xob,latt,lont,vob,uob,yr,mon,day,hr,mn,xob_plot,yob_plot,tob_plot,fname = getData()
#
# X,LL,obs = obsMatrix(tob,yob,xob,latt,lont,vob,uob)
# U,Ku,V,Kv,Xg,xg0,yg0,tg,model_u,model_v,lon_degs,lat_degs = vel_field(X,LL,obs,tob)
#
# SaveMat(fname,xg0,yg0,U,V,X,Ku,Kv,tg,obs,LL,yob,xob,tob,yr,mon,day,hr,mn,uob,vob,latt,lont,xob_plot,yob_plot,tob_plot,model_u,model_v,lon_degs,lat_degs)

