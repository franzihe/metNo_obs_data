# %%
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import math
import functions as fct
from matplotlib.lines import Line2D
from datetime import datetime, timedelta


# %%
ds = xr.open_dataset('https://thredds.met.no/thredds/dodsC/remotesensingradiosonde/2023/02/andoya_20230131-20230228.nc', decode_times=True)

# %%
now = datetime.utcnow()
print('Now:', now)

# %%
time_id = '{}-{}-{}T{}'.format(now.year, now.month, now.day,now.hour)
# time_id = '2023-2-17T08'


# %%
sonde_time = ds.time.sel(time=time_id, method='nearest')#['time.year'].values
print('latest Radiosonde launch @: {}-{}-{}T{}-{} UTC'.format(sonde_time['time.year'].values, 
                                                          sonde_time['time.month'].values,
                                                          sonde_time['time.day'].values,
                                                          sonde_time['time.hour'].values,
                                                          sonde_time['time.minute'].values,),)

# %%
# 	Atmospheric Pressure
pres = ds.air_pressure.sel(time=time_id, method='nearest')             #'units': 'hPa'
pres= pres*100
pres = pres.assign_attrs(units='Pa', valid_max=1100.0*100)
pres = pres.drop('time')


# %%
ds = (ds.sel(time=time_id, method='nearest').drop('time')).assign_coords(time_from_launch=pres.values)
ds = ds.rename({'time_from_launch': 'pressure'})
ds = ds.dropna('pressure')

# %%
# 	Atmospheric Pressure
pres = ds.air_pressure              #'units': 'hPa'
pres= pres*100
pres = pres.assign_attrs(units='Pa', valid_max=1100.0*100)

# Altitude ('Derived from geopotential_height')
z = ds.altitude                  #'units': 'm'

# Temperature
temp = ds.air_temperature           #'units': 'K'
temp = temp-273.15
temp = temp.assign_attrs(units='Celsius', valid_max=400-273.15)
# 	Dewpoint Temperature
Td = ds.dew_point_temperature    #'units': 'K'
Td = Td-273.15
Td = Td.assign_attrs(units='Celsius', valid_max=400-273.15)
# 	Wind Direction
wd = ds.wind_from_direction     #'units': 'degree'
wd = np.deg2rad(wd)
wd = wd.assign_attrs(units='rad', valid_max='2pi')
# 	Wind Speed
ws = ds.wind_speed               #'units': 'm/s'

# calculate wind components
u, v = fct.get_wind_components(ws, wd)





# %%
fig, axsm = plt.subplots(1,2, figsize=(9,10),gridspec_kw={'width_ratios': [12, 1]})
axs = axsm.flatten()

# sounding
fct.plot_sounding_axes(axs[0])

# plot Temperature
linecolor_T = fct.linecolor_T
linewidth_T = fct.linewidth_T
axs[0].semilogy(temp + fct.skew(pres),pres, base=math.e, color =linecolor_T, linewidth = (linewidth_T+1.5))


# plot dewpoint
linecolor_Td = fct.linecolor_Td
linewidth_Td = fct.linewidth_Td
axs[0].semilogy(Td + fct.skew(pres), pres, base=math.e, color=linecolor_Td, linewidth = (linewidth_Td+1.5))

# legend
tT = r'Temperature'
lT = Line2D(range(10), range(10), linestyle='-', marker='', linewidth=(linewidth_T+1.5), color=linecolor_T)

tTd = r'Dew-point Temperature'
lTd = Line2D(range(10), range(10), linestyle='-', marker='', linewidth=(linewidth_Td+1.5), color=linecolor_Td)

axs[0].legend((lT, lTd,),(tT, tTd, ), loc='upper right',
           fontsize=12, handlelength=5)

# wind barbs
fct.plot_wind_axes(axs[1])
# # fct.plot_wind_barbs(axs[1],z[::15],pres[::15],u[::15],v[::15])
fct.plot_wind_barbs(z,pres,u,v)


# # plot labels for std heights
# for plvl in fct.plevs_std:
#     zlvl = fct.interp_height(z/1000.,pres,plvl)
#     fct.label_m(fct.Tmin+2.55,plvl, str(int(zlvl)), axs[0])
    
    

# # Adjust plot margins.
plt.tight_layout()
plt.subplots_adjust(left=0.03, bottom=0.03, right=0.97, top=0.97, wspace=0.02, hspace=0.12)

axs[0].set_title('{}-{}-{}T{}'.format(sonde_time['time.year'].values,sonde_time['time.month'].values, sonde_time['time.day'].values, sonde_time['time.hour'].values), 
                 fontsize = 18, color=fct.blue);


# %%



