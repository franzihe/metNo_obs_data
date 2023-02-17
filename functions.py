import numpy as np
import math
import matplotlib.pyplot as plt

### Define colorbar colors
champ = 255.
blue = np.array([1,74,159])/champ           # for the date

# adjusted functions from https://github.com/cwebster2/pyMeteo

## constants
# This defines the skew-angle of the T axis
skew_angle = 37.5 
"""This defines the skewness of the T axis"""

# These define the domain plotted.  T values are at @1000 mb (e.g. unskewed)
Tmin = -40.0
"""Sets the left boundary of the plot. Temperature at 1000 mb (C)"""
Tmax = 40.0
"""Sets the right boundary of the plot. Temperature at 1000 mb (C)"""
pbot = 100000.0
# pbot = 95000.0
"""Sets the bottom boundary of the plot. Pressure (Pa)"""
#ptop = 10000.0
ptop = 20000.0
"""Sets the top boundary of the plot. Pressure (Pa)"""

## Values below used for plotting 
dp = 5000.0
"""The delta pressure used in calculating adiabats"""
ptickbot = 100000.0
"""Lowest elevated pressure level to be labelled in the plot (Pa)"""
pticktop = 10000.0
"""Highest elevated pressure level to be labelled in the plot (Pa)"""
# tickdp = 10**4
tickdp = 20000
"""The spacing between pressure labels in the plot (Pa)"""
plevs = np.arange(pbot,ptop-1,-dp)
"""List of pressure levels to do calculations on"""
fontscalefactor = 5

## Lines to plot
isotherms    = np.arange(-150,51,10)              # in degrees C
"""List of isotherms to plot on the skew-t.  In degrees C"""
isobars      = np.arange(ptickbot,ptop-1,-10000)   # in Pa
"""List of isobars to plot on the skew-t.  In Pascals"""
dry_adiabats = np.arange(-40,210,10)              # in degrees C
"""List of dry adiabats to plot on the skew-t.  In degrees C"""
#moist_adiabats = np.concatenate((np.arange(-15.,10.1,5.),np.arange(12.,45.1,2.)))
moist_adiabats = np.arange(8,32,4)
"""List of moist adiabats to plot on the skew-t.  In degrees C"""
#mixing_ratios = [0.2,0.4,0.8,1,2,3,4,6,8,10,14,18,24,32,40]
mixing_ratios = [0.4,1,2,3,5,8,12,20]
"""List of water vapor mixing ratio lines to plot. In g/kg"""

## Linewidths
#lw_major = 0.6
lw_major = 2.3
"""Line width of 'major' lines.  E.g. Lines plotted at 10 C intervals or 50 mb intervals"""
#lw_minor = 0.42
lw_minor = lw_major
"""Line width of 'minor' lines.  E.g. Lines not plotted at the major intervals"""

## Linecolors
lc_major = 'grey'
"""Line color of 'major' lines.  E.g. Lines plotted at 10 C intervals or 50 mb intervals"""
lc_minor = 'lightgrey'
"""Line color of 'minor' lines.  E.g. Lines not plotted at the major intervals"""

## Skew-T line parameters
linecolor_T = 'black'
"""Line color of environmental temperature profile"""
linewidth_T = 2.5
"""Line width of environmental temperature profile"""

linecolor_Td = 'green'
"""Line color of environmental dew-point temperature profile"""
linewidth_Td = 2.5
"""Line width of environmental dew-point temperature profile"""

linecolor_Parcel_T = 'red'
"""Line color of lifted surface parcel temperature profile"""
linewidth_Parcel_T = 1.0
"""Line width of lifted surface parcel temperature profile"""

linecolor_Twb = 'blue'
"""Line color of environmental wet-bulb temperature profile"""
linewidth_Twb = 0.5
"""Line width of environmental wet-bulb temperature profile"""

linecolor_Tve = 'black'
"""Line color of environmental virtual temperature profile"""
linewidth_Tve = 0.7
"""Line width of environmental virtual temperature profile"""
linestyle_Tve = '--'
"""Line style of environmental virtual temperature profile"""

linecolor_Tvp = 'red'
"""Line color of lifted surface parcel virtual temperature profile"""
linewidth_Tvp = 0.7
"""Line width of lifted surface parcel virtual temperature profile"""
linestyle_Tvp = '--'
"""Line style of lifted surface parcel virtual temperature profile"""

#for plotted lines
pb_plot=105000
pt_plot=10000
pt_plot2=20000
dp_plot=1000
plevs_plot = np.arange(pb_plot,pt_plot-1,-dp_plot)
plevs_plot2 = np.arange(pb_plot,pt_plot2-1,-dp_plot)
plevs_std = [100000,85000,70000,50000,40000,30000,25000,20000,]#15000]

#TODO: enforce square w/ aspect ratio of plot
#Domain of the hodograph
umin = -22.5
umax = 27.5
vmin = -12.5
vmax = 27.5

# constants
p00 = 100000.   # reference pressure
T00 = 273.15

L = 2.501e6    # latent heat of vaporization
Rd = 287.04         # gas constant dry air
Rv = 461.5              # gas constant water vapor
epsilon = Rd/Rv
cpd = 1005.7             # what about cpd vs cpv
cp = 1005.7              # what about cpd vs cpv
cv = 718.

kappa = (cp-cv)/cp
kappa_d = Rd/cp
g = 9.81


def get_wind_components(wind_speed, wind_direction):
    r"""Calculate the U, V wind vector components from the speed and direction.
    Parameters
    ----------
    speed : `pint.Quantity`
        Wind speed (magnitude)
    wind_direction : `pint.Quantity`
        Wind direction, specified as the direction from which the wind is
        blowing (0-2 pi radians or 0-360 degrees), with 360 degrees being North.
    Returns
    -------
    u, v : tuple of `pint.Quantity`
        The wind components in the X (East-West) and Y (North-South)
        directions, respectively.
        
    -------
    Taken from metpy.calc https://github.com/Unidata/MetPy/blob/ba3f58a377001831283c7a78ea7813245e8c3abe/src/metpy/calc/basic.py
    """
    u = -wind_speed * np.sin(wind_direction)
    v = -wind_speed * np.cos(wind_direction)
    return u, v

def skew(p):
    """Puts the skew in skew-T
    :parameter p: pressure level of the point.
    This calculates the skew of the T axis for a point to plot.  
    This assumes a logarithmic y axes and uses the variable
    :py:data:skew_angle to determine the skew.  This is the 
    magic that turns a cartesian plot into a Skew-T/Log-p plot.
    """
    return skew_angle * np.log(p00/p)

def label(x, y, s, c, r, axes):
	axes.text(x,y*100,s, verticalalignment='center', horizontalalignment='center', weight='bold', fontsize=2*fontscalefactor, color=c, rotation=r)


def met_T(theta,p):
    """Convert Potential Temperature :math:`\\theta` to Temperature
    :parameter theta: Potential temperature (K)
    :parameter p: Pressure (Pa)
    :returns: Temperature (K)
    """
    return theta * (p00/p)**-kappa_d

def draw_isotherms(axes):
    """Plot isotherms on axes
    :parameter axes: The axes to draw on
    :type axes: :py:class:`matplotlib.axes`
    This function draws isotherms every 10 C.
    """
    col_isotherm = np.array([255,26,26])/255.
    for T in isotherms:
        if (T % 10 == 0):
#           axes.semilogy(T + skew(plevs_plot), plevs_plot, basey=math.e, color = lc_major, linewidth= lw_major)
           axes.semilogy(T + skew(plevs_plot), plevs_plot, base=math.e, color = col_isotherm, linewidth= lw_major)
        else:
#           axes.semilogy(T + skew(plevs_plot), plevs_plot, basey=math.e, color = lc_minor, linewidth= lw_minor)
           axes.semilogy(T + skew(plevs_plot), plevs_plot, base=math.e, color = col_isotherm, linewidth= lw_minor)
    for T in np.arange(-40, 40, 10):
#        label(T+skew(87500),875, str(T), 'red', 90.-skew_angle, axes)
        label(T+skew(85000),875, str(T),  col_isotherm, 90.-skew_angle, axes)
    
#    for T in np.arange(-80, -40, 10):
 #       label(T+skew(17500),175, str(T), 'red', 90.-skew_angle, axes)
    p = [325,425,550,700]
    T = np.arange(-80, -40, 10)
    for i in range(0,4):
        label(-40+skew(87500),p[i], str(T[i]),  col_isotherm, 90.-skew_angle, axes)
        
        
def draw_isobars(axes):
    """Plot isobars on axes
    :parameter axes: The axes to draw on
    :type axes: :py:class:`matplotlib.axes`
    This function draws isobars at intervals of 2*dp.
    """
    col_isobar = np.array([189,190,193])/255.
    for i in isobars:
        if (i % 5000 == 0):
#            axes.plot([Tmin, Tmax], [i,i], color = lc_major, linewidth = lw_major)
            axes.plot([Tmin, Tmax], [i,i], color = col_isobar, linewidth = lw_major)
        else:
#            axes.plot([Tmin, Tmax], [i,i], color = lc_minor, linewidth = lw_minor)
            axes.plot([Tmin, Tmax], [i,i], color = col_isobar, linewidth = lw_minor)
    for i in np.arange(1000,150,-50):
#        label(-10-((1000-i)*.025),i,str(i),'black',0, axes)
        label(-10-((1000-i)*.025),i,str(i),col_isobar,0, axes)
        
def draw_dry_adiabat(axes):
    """Plot dry adiabats on axes
    :parameter axes: The axes to draw on
    :type axes: :py:class:`matplotlib.axes`
    This function calculates dry adiabats
    and plots these lines.  Adiabats are calculated 
    every 10 K
    """
    col_dry_adiabat = np.array([61,171,226])/255.
    for T in dry_adiabats:
        dry_adiabat = met_T(T+T00,plevs_plot) - T00 + skew(plevs_plot)
        if (T % 10 == 0):
#            axes.semilogy(dry_adiabat, plevs_plot, basey=math.e, color = lc_major, linewidth = lw_major)
            axes.semilogy(dry_adiabat, plevs_plot, base=math.e, color = col_dry_adiabat, linewidth = lw_major)
        else:
#            axes.semilogy(dry_adiabat, plevs_plot, basey=math.e, color = lc_minor, linewidth = lw_minor)
            axes.semilogy(dry_adiabat, plevs_plot, base=math.e, color = col_dry_adiabat, linewidth = lw_minor)
            
    for T in np.arange(-20, 120, 10):
        p = (600. - 3.5*T)*100.
        x = met_T(T+T00,p) -T00 + skew(p)
        x1 = met_T(T+T00,p+.5*dp_plot) -T00 + skew(p+.5*dp_plot)
        x2 = met_T(T+T00,p-.5*dp_plot) -T00 + skew(p-.5*dp_plot)
        dx = x2-x1
        theta = math.atan2(-dp_plot,-dx) * 180/math.pi +37
#        label(x,p/100,str(T),'black',theta, axes)
        label(x,p/100,str(T),col_dry_adiabat,theta, axes)

def es(T):
	#return 611.2 * np.exp((Lv(T)/Rv)*((1./T00)-(1./T)))
        return 611.2 * np.exp(17.67*(T-T00)/(T-29.65))


def Lv(T):
#TODO: Temp dependance
	return L

def w_vs(T,pd):
	return epsilon * (es(T)/pd)

def dTdz_moist(T,p):
  pd = p - es(T)
  num = 1. + ((Lv(T) * w_vs(T,pd))/(Rd*T))
  den = 1. + ((Lv(T)**2 * w_vs(T,pd))/(cpd * Rv * T**2))
  return (-g/cpd)*(num/den)

def dTdp_moist(T,p):
	return dTdz_moist(T,p) * -((Rd*T)/(p*g))
     
def draw_moist_adiabat(axes):
    """Plot moist adiabats on axes
    :parameter axes: The axes to draw on
    :type axes: :py:class:`matplotlib.axes`
    This function calculates moist adiabats
    and plots these lines.  Adiabats are calculated for
    values of T at 1000mb from -15 to 45 C every 5 C between
    -15 and 10 C and every 2.5 C between 12.5 and 45 C.
    """
    col_moist_adiabat = np.array([202,137,194])/255.
    ps_blo = [p for p in plevs_plot if p > 100000]
    ps_blo.reverse()
    ps = [p for p in plevs_plot2 if p < 100000]
    for T in moist_adiabats:
        T_1000 = T = T + T00
        moist_adiabat = []
        # work backwards from 1000mb
        for p in ps_blo:
            T += dTdp_moist(T,p) * dp_plot
            moist_adiabat.append(T - T00 + skew(p))
        #reverse list order
        moist_adiabat.reverse()
        # insert 1000mb point
        T = T_1000
        moist_adiabat.append(T - T00)
        # work forwards from 1000mb
        for p in ps:
            T -= dTdp_moist(T,p) * dp_plot
            moist_adiabat.append(T - T00 + skew(p))
            # draw labels
#            if (p == 22000):
            if (p == 41000):
                if (T_1000 >= T00 and T_1000 <= 30+T00):
                    label(T-T00+skew(p),p/100,str(int(T_1000-T00)),col_moist_adiabat, 0, axes)
        if (int(T_1000 - T00) % 5 == 0):            
#            axes.semilogy(moist_adiabat, plevs_plot2, basey=math.e, color = lc_major, linewidth = lw_major)
            axes.semilogy(moist_adiabat, plevs_plot2, base=math.e, color = col_moist_adiabat, linewidth = lw_major, linestyle= '--')
        else:
 #           axes.semilogy(moist_adiabat, plevs_plot2, basey=math.e, color = lc_minor, linewidth = lw_minor)
            axes.semilogy(moist_adiabat, plevs_plot2, base=math.e, color = col_moist_adiabat, linewidth = lw_major, linestyle= '--')
#        label(T+skew(45000),410,str(T),col_moist_adiabat, 0, axes)

def TMR(W, p):
  # Computes temperature on mixing ratio w at pressure p.
  # TMR in C, w in g/kg dry air, p in millibars.
  # TODO: change this to something else?
  x = np.log10(W * p / (622 + W))
  TMR = 10 ** (0.0498646455 * x + 2.4082965) - 280.23475 + 38.9114 * ((10 ** (0.0915 * x) - 1.2035) ** 2)
  return TMR


def draw_water_mix_ratio(axes):
    """Plot lines of constant water vapor mixing ratio on axes
    :parameter axes: The axes to draw on
    :type axes: :py:class:`matplotlib.axes`
    This function calculates isolines of constant water vapor
    mixing ratio and plots these lines.  Values of w calculated
    are given by the list variable w.
    """
    col_water_mix_ratio = np.array([104,208,121])/255.
    #TODO: put w and the top plevel for plotting somewhere configurable
    ps = [p for p in plevs if p>=20000 and p<=105000]
    for W in mixing_ratios:
        water_mix = []
        for p in ps:
            T = TMR(W,p/100.) 
            water_mix.append(T + skew(p))
 #       axes.semilogy(water_mix, ps, basey=math.e, color = 'grey', linestyle = '--', linewidth = .5)
        axes.semilogy(water_mix, ps, base=math.e, color = col_water_mix_ratio, linestyle = '--', linewidth = lw_major)

        # Label the isoline
#        T = TMR(W,1075.)
        T = TMR(W,625.)
#        label(T+skew(107500.), 1075, str(W), 'black', -15, axes)
        label(T+skew(62500.), 625, str(W), col_water_mix_ratio, -15, axes)
        
def remove_tick_labels(axes):
  axes.tick_params(axis='x', top='off', bottom='off', which='both')#, labelsize=0)
  axes.tick_params(axis='y', left='off', right='off', which='both')#, labelsize=0)
  for xlabel_i in axes.get_xticklabels():
    xlabel_i.set_visible(False)
    xlabel_i.set_fontsize(0.0)
  for xlabel_i in axes.get_yticklabels():
    xlabel_i.set_fontsize(0.0)
    xlabel_i.set_visible(False)
    
# Puts the skew in skew-T
def skew(p):
    """Puts the skew in skew-T
    :parameter p: pressure level of the point.
    This calculates the skew of the T axis for a point to plot.  
    This assumes a logarithmic y axes and uses the variable
    :py:data:skew_angle to determine the skew.  This is the 
    magic that turns a cartesian plot into a Skew-T/Log-p plot.
    """
    return skew_angle * np.log(p00/p)

def plot_sounding_axes(axes):
  """Plots Skew-T/Log-P axes
  This will plot isotherms, isobars, dry and moist adiabats, 
  lines of constant water vapor mixing ratio, labels and 
  setup the y axes to be reversed.
  :paramter axes: The axes to draw on
  """
  draw_isotherms(axes)
  draw_isobars(axes)
  draw_dry_adiabat(axes)
  draw_moist_adiabat(axes)
  draw_water_mix_ratio(axes)
  remove_tick_labels(axes)
  axes.axis([Tmin, Tmax, pbot, ptop])
  axes.set_ylim(axes.get_ylim()[::1])
  

def draw_wind_line(axes):
	wind_line = []
	for p in plevs_plot:
		wind_line.append(0)
	axes.semilogy(wind_line, plevs_plot, color='black', linewidth=.5)
	# plot circles at certain levels?
	for p in plevs_std:
		axes.semilogy(0,p, color='black', markersize=3, marker='.')
  
def plot_wind_axes(axes):
  # plot wind barbs
#   draw_wind_line(axes)
  axes.set_axis_off()
  axes.axis([-0.5,0.5,pbot,ptop])
  
def plot_wind_barbs(z, p, u, v):
   
   pp = p.assign_coords(pressure=z.values).rename({'pressure':'height'})
   u2 = u.assign_coords(pressure=z.values).rename({'pressure':'height'})
   v2 = v.assign_coords(pressure=z.values).rename({'pressure':'height'})
   
   for i in range(0,12*1000,200):
    # print(i)
    if (pp.sel(height=i,method='nearest') > pt_plot):
      #   print(i, pp.sel(height=i,method='nearest').values)
        plt.barbs(0,pp.sel(height=i,method='nearest'),
                  u2.sel(height=i,method='nearest'),
                  v2.sel(height=i,method='nearest'), 
                  length=6, linewidth=.5, pivot='middle')
   
   # for i in np.arange(100000,15000,-100):
   #  # print(i)
   #  if (p.sel(pressure=i,method='nearest') > pt_plot):
   #    #   print(i, p.sel(pressure=i,method='nearest').values)
   #      plt.barbs(0,
   #                p.sel(pressure=i,method='nearest'),
   #                u.sel(pressure=i,method='nearest'),
   #                v.sel(pressure=i,method='nearest'), 
   #                length=6, linewidth=.5, pivot='middle')
   #  for i in np.arange(0,len(z)):
   #      if (p[i] > pt_plot):
   #          plt.barbs(0,p[i],u[i],v[i], length=6, linewidth=.5, pivot='middle')
            
   #          # plt.barbs(Tmin+4,p[i],u[i],v[i], length=6, linewidth=.5)
            
   # # if (u is not None and v is not None):
   # #    #draw_wind_line(axes)
   # #    for i in np.arange(0,len(z),2):
   # #      if (p[i] > pt_plot):
   # #          plt.barbs(Tmin+4,pres[i],u[i],v[i], length=8, linewidth=2.)
            
def interp_height(z, p, plvl):
   #interpolates height to a pressure level

   nlevs = len(p)

   # check bounds
   if plvl > p[0]:
      return 0
   if plvl < p[nlevs-1]:
      return -1

   z0 = 0
   while p[z0] > plvl:
      z0 = z0 + 1 

   if p[z0] == plvl:
      return z[z0]

   z1 = nlevs-1
   while p[z1] < plvl:
      z1 = z1 - 1

   #interpolate to height.  
   # Code adapted from NSHARP 95 John Hart NSSFC KCMO. 

   zdiff = z[z1] - z[z0]
   pdiff = np.log( p[z0] / p[z1])
   pdist = np.log( p[z0] / plvl)
   height = z[z0] + (( pdist / pdiff) * zdiff)

   return height


def label_m(x, y, s, axes):
	axes.text(x,y,s, verticalalignment='center', horizontalalignment='right', fontsize=3*fontscalefactor)