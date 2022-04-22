########### updated 22/04/2022 ###########

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import savefig
import matplotlib.gridspec as gridspec
from matplotlib.patches import Polygon
import astropy.units as u
from astropy.time import Time
#from datetime import datetime
from astropy.wcs import WCS
import astroplan as ap
import astroplan.plots as app
from astroplan.plots import plot_sky
import matplotlib.dates as mdates
from regions import RectangleSkyRegion, RectanglePixelRegion
from astroplan import FixedTarget
from astropy.coordinates import SkyCoord
import os


#########################################################################################################
######################################## Planning observations ##########################################
#########################################################################################################

#set plots parameters
plotpar = {'axes.labelsize': 25,
           'xtick.labelsize': 20,
           'ytick.labelsize': 20,
           'text.usetex': True}
plt.rcParams.update(plotpar)

###################################
#planning observations with astropy
###################################
class cha(object):
    def __init__(self):
        self.starid   = None
        self.site     = None
        self.time     = None
        self.radius   = None
        self.savefig  = None
        self.by       = None

    @classmethod
    #getting observing sites
    def get_sites(self):
        # make sure stuff is up to date
        #download_IERS_A()
        from astropy.coordinates import EarthLocation
        lista = EarthLocation.get_site_names()
        return lista
    
    
    @classmethod
    #finding chart
    def fch_obs(self, starid, radius, by=None, ra=None, dec=None, savefig=None):
        
        if not savefig:
            savefig = 'no'
        else:
            savefig = savefig
        if not by:
            by = 'id'
        else:
            by = 'coord'      
        if not ra:
            ra = []
        else:
            ra = ra
        if not dec:
            dec = []
        else:
            dec = dec
        
        if by == 'id':
            target = ap.FixedTarget.from_name(starid)
        if by == 'coord':
            target_coord = SkyCoord(ra=ra*u.deg, dec=dec*u.deg)
            target = FixedTarget(coord=target_coord, name=starid)
        
        #Plotting
        fig = plt.figure(figsize=(7,7))
        axn, hdu = app.plot_finder_image(target, fov_radius=5*u.arcmin)
        
        wcs = WCS(hdu.header)
        wcs_im = wcs.celestial
        #transform = axn.get_transform('icrs')
        #making the rectangle
        rectangle_sky   = RectangleSkyRegion(center=target.coord, width=0.3 * u.arcmin, 
                                           height=0.3 * u.arcmin, angle=0 * u.deg)
        rectangle_topix = rectangle_sky.to_pixel(wcs_im)
        rectangle_pix   = RectanglePixelRegion(rectangle_topix.center, 
                                               rectangle_topix.width, rectangle_topix.height, 
                                               rectangle_topix.angle)
        rectangle = Polygon(rectangle_topix.corners, facecolor='None', edgecolor='r', ls='-')
        axn.add_artist(rectangle)
        if savefig == 'yes':
            plt.savefig('FC_'+starid+'.pdf', bbox_inches='tight')
    
    
    @classmethod
    #sky position
    def sky_position(self, starid, site, time, by=None, ra=None, dec=None, savefig=None):
        
        if not savefig:
            savefig = 'no'
        else:
            savefig = savefig
        if not by:
            by = 'id'
        else:
            by = 'coord'      
        if not ra:
            ra = []
        else:
            ra = ra
        if not dec:
            dec = []
        else:
            dec = dec
        
        if by == 'id':
            target = ap.FixedTarget.from_name(starid)
        if by == 'coord':
            target_coord = SkyCoord(ra=ra*u.deg, dec=dec*u.deg)
            target = FixedTarget(coord=target_coord, name=starid)
            
        site      = ap.Observer.at_site(site)
        
        fig             = plt.figure(figsize=(7,7))
        observe_time    = Time(time)
        #observe_time    = Time(datetime.utcnow())
        sunset_tonight  = site.sun_set_time(observe_time, which= 'previous')#, horizon=-18*u.deg)
        sunrise_tonight = site.sun_rise_time(sunset_tonight, which='next')#, horizon=-18*u.deg)
        observe_times   =  sunset_tonight + (sunrise_tonight - sunset_tonight)*np.linspace(0, 1, 20)
        plot_sky(target, site, observe_times)
        #plt.title('%s'%(star))
        #plt.legend(loc='best')
        plt.tight_layout()
        if savefig == 'yes':
            plt.savefig('SKP_'+starid+'.pdf', bbox_inches='tight')
    
    
    @classmethod
    #moon separation
    def moon_sep(self, starid, site, time, by=None, ra=None, dec=None, savefig=None):
        
        if not savefig:
            savefig = 'no'
        else:
            savefig = savefig
        if not by:
            by = 'id'
        else:
            by = 'coord'      
        if not ra:
            ra = []
        else:
            ra = ra
        if not dec:
            dec = []
        else:
            dec = dec
        
        if by == 'id':
            target = ap.FixedTarget.from_name(starid)
        if by == 'coord':
            target_coord = SkyCoord(ra=ra*u.deg, dec=dec*u.deg)
            target = FixedTarget(coord=target_coord, name=starid)
            
        site      = ap.Observer.at_site(site)
        
        observe_time    = Time(time)
        sunset_tonight  = site.sun_set_time(observe_time, which= 'nearest')#, horizon=-18*u.deg)
        sunrise_tonight = site.sun_rise_time(sunset_tonight, which='nearest')#, horizon=-18*u.deg)
        
        if sunset_tonight.ymdhms.hour <= 23.0 and sunset_tonight.ymdhms.hour >= 15.0:
            time_grid = Time(time) + np.arange(sunset_tonight.ymdhms.hour-24, 
                                           sunrise_tonight.ymdhms.hour+1, 0.1)*u.hour
            
        if sunset_tonight.ymdhms.hour >= 0.0 and sunset_tonight.ymdhms.hour <= 13.0:
            time_grid = Time(time) + np.arange(sunset_tonight.ymdhms.hour, 
                                           sunrise_tonight.ymdhms.hour+1, 0.1)*u.hour
        
        moon      = site.moon_altaz(time_grid) # the moon :)
        
        fig             = plt.figure(figsize=(8,5))
        gs  = gridspec.GridSpec(1, 1)
        ax  = plt.subplot(gs[0])
        hours = mdates.HourLocator()
        hoursFmt = mdates.DateFormatter('%H')
        ax.xaxis.set_major_locator(hours)
        ax.xaxis.set_major_formatter(hoursFmt)
        ax.plot_date(time_grid.plot_date, moon.separation(target.coord).deg, fmt='-',
                  label=target.name)
        ax.axhline(27, color='r', linestyle='dashed', label='Minimum Lunar Separation')
        ax.legend(loc='best')
        ax.set_ylabel('Lunar distance [deg]')
        ax.set_xlabel('Time [UTC]')
        plt.tight_layout()
        if savefig == 'yes':
            plt.savefig('Moon_'+starid+'.pdf', bbox_inches='tight')
    
    
    @classmethod
    #object visibility
    def plan_one_obs(self, starid, site, time, by=None, ra=None, dec=None, savefig=None):
        
        if not savefig:
            savefig = 'no'
        else:
            savefig = savefig
        if not by:
            by = 'id'
        else:
            by = 'coord'      
        if not ra:
            ra = []
        else:
            ra = ra
        if not dec:
            dec = []
        else:
            dec = dec
        
        if by == 'id':
            target = ap.FixedTarget.from_name(starid)
        if by == 'coord':
            target_coord = SkyCoord(ra=ra*u.deg, dec=dec*u.deg)
            target = FixedTarget(coord=target_coord, name=starid)
            
        site      = ap.Observer.at_site(site)
        
        observe_time    = Time(time)
        sunset_tonight  = site.sun_set_time(observe_time, which= 'nearest')#, horizon=-18*u.deg)
        sunrise_tonight = site.sun_rise_time(sunset_tonight, which='nearest')#, horizon=-18*u.deg)
        
        if sunset_tonight.ymdhms.hour <= 23.0 and sunset_tonight.ymdhms.hour >= 15.0:
            time_grid = Time(time) + np.arange(sunset_tonight.ymdhms.hour-24, 
                                           sunrise_tonight.ymdhms.hour+1, 0.1)*u.hour
            
        if sunset_tonight.ymdhms.hour >= 0.0 and sunset_tonight.ymdhms.hour <= 13.0:
            time_grid = Time(time) + np.arange(sunset_tonight.ymdhms.hour, 
                                           sunrise_tonight.ymdhms.hour+1, 0.1)*u.hour
        
        moon      = site.moon_altaz(time_grid) # the moon :)
        
        ####plotting####
        #visibility
        
        fig = plt.figure(figsize=(10,7))
        gs  = gridspec.GridSpec(1, 1)
        #gs.update(left=0.13, right=0.97, wspace=0.2, hspace=0.25, top=0.96, bottom=0.1)
        
        ax0  = plt.subplot(gs[0,:])
        app.plot_altitude(target, site, time_grid, brightness_shading=True, 
                          airmass_yaxis=True, ax=ax0)
        ax0.plot(time_grid.plot_date, moon.alt, color = 'gray', ls = ':', label = 'Moon')
        ax0.legend(loc='best')
        
        plt.tight_layout()
        if savefig == 'yes':
            plt.savefig('Visibility_'+starid+'.pdf', bbox_inches='tight')
        
    @classmethod
    #all in one
    def plan_all_obs(self, starid, site, time, radius, by=None, ra=None, dec=None):
        
        if not by:
            by = 'id'
        else:
            by = 'coord'      
        if not ra:
            ra = []
        else:
            ra = ra
        if not dec:
            dec = []
        else:
            dec = dec
        
        if by == 'id':
            #visibility
            self.plan_one_obs(starid, site, time, savefig='yes')
            #finding chart
            self.fch_obs(starid, radius, savefig='yes')
            #sky position
            self.sky_position(starid, site, time, savefig='yes')
            #lunar separation
            self.moon_sep(starid, site, time, savefig='yes')
        
        if by == 'coord':
            #visibility
            self.plan_one_obs(starid, site, time, by=by, ra=ra, dec=dec, savefig='yes')
            #finding chart
            self.fch_obs(starid, radius, by=by, ra=ra, dec=dec, savefig='yes')
            #sky position
            self.sky_position(starid, site, time, by=by, ra=ra, dec=dec, savefig='yes')
            #lunar separation
            self.moon_sep(starid, site, time, by=by, ra=ra, dec=dec, savefig='yes')
        
        #merging pdfs
        from PyPDF2 import PdfFileReader, PdfFileWriter, PdfFileMerger
        
        figs = ['Visibility_'+starid+'.pdf', 'FC_'+starid+'.pdf', 'SKP_'+starid+'.pdf', 'Moon_'+starid+'.pdf']
        mergedObject = PdfFileMerger(strict=False)
        for fileNumber in figs:
            mergedObject.append(PdfFileReader(fileNumber, 'rb'), import_bookmarks=False)
        mergedObject.write('plan_obs_'+starid+'.pdf')
        
        #removing files
        for i in figs:
            if os.path.isfile(i):
                os.remove(i)
    
    @classmethod
    #multiple observations
    def plan_mult_obs(self, starids, site, time, by=None, ra=None, dec=None):
        
        if not by:
            by = 'id'
        else:
            by = 'coord'      
        if not ra:
            ra = []
        else:
            ra = ra
        if not dec:
            dec = []
        else:
            dec = dec
        
        if by == 'id':
            targets   = [ap.FixedTarget.from_name(target) for target in starid]
        if by == 'coord':
            target_coord = SkyCoord(ra=ra*u.deg, dec=dec*u.deg)
            targets      = [FixedTarget(coord=target_coord[i], name=starid[i]) 
                            for i in range(len(target_coord))]
        
        site      = ap.Observer.at_site(site)
        
        #visibility
        observe_time    = Time(time)
        sunset_tonight  = site.sun_set_time(observe_time, which= 'nearest')#, horizon=-18*u.deg)
        sunrise_tonight = site.sun_rise_time(sunset_tonight, which='nearest')#, horizon=-18*u.deg)
        
        if sunset_tonight.ymdhms.hour <= 23.0 and sunset_tonight.ymdhms.hour >= 15.0:
            time_grid = Time(time) + np.arange(sunset_tonight.ymdhms.hour-24, 
                                           sunrise_tonight.ymdhms.hour+1, 0.1)*u.hour
            
        if sunset_tonight.ymdhms.hour >= 0.0 and sunset_tonight.ymdhms.hour <= 13.0:
            time_grid = Time(time) + np.arange(sunset_tonight.ymdhms.hour, 
                                           sunrise_tonight.ymdhms.hour+1, 0.1)*u.hour
        
        moon      = site.moon_altaz(time_grid) # the moon :)
        
        ####plotting####
        #visibility
        
        fig = plt.figure(figsize=(13, 17))
        gs  = gridspec.GridSpec(2, 1)
        #gs.update(left=0.13, right=0.97, wspace=0.2, hspace=0.25, top=0.96, bottom=0.1)
        
        ax0  = plt.subplot(gs[0])
        app.plot_altitude(targets, site, time_grid, brightness_shading=True, 
                          airmass_yaxis=True, ax=ax0)
        ax0.plot(time_grid.plot_date, moon.alt, color = 'gray', ls = ':', label = 'Moon')
        ax0.legend(bbox_to_anchor=(1.08,1), loc="upper left", fontsize=12)
        
        
        #sky position
        ax1  = plt.subplot(gs[1])
        hours = mdates.HourLocator()
        hoursFmt = mdates.DateFormatter('%H')
        ax1.xaxis.set_major_locator(hours)
        ax1.xaxis.set_major_formatter(hoursFmt)
        for target in targets:
            ax1.plot_date(time_grid.plot_date, moon.separation(target.coord).deg, fmt='-', 
                          label=target.name)
        ax1.axhline(27, color='r', linestyle='dashed', label='Minimum Lunar Separation')
        ax1.legend(loc='best')
        ax1.set_ylabel('Lunar distance [deg]')
        ax1.set_xlabel('Time [UTC]')
        ax1.legend(bbox_to_anchor=(1.08,1), loc="upper left", fontsize=12)
        
        plt.savefig('plan_mult_obs.pdf', bbox_inches='tight')


#############
#useful tools
#############
class chatool(object):
    def __init__(self):
        self.v          = None
        self.snr        = None
        self.grating    = None
        self.slit       = None
        self.seeing     = None
        self.bin        = None

    @classmethod
    #Exposure time calculator for the Goodman spectrograph of the SOAR telescope. 
    #Based on the calculation done by David Sanmartin (http://www.lna.br/soar/NSO/ETC/index.html)
    #Information about grating can be found at 
    #https://noirlab.edu/science/programs/ctio/instruments/goodman-high-throughput-spectrograph/overview
    def texp_goodman(self, v, snr, seeing, grating, slit):
        from scipy import special
        
        #grating and dispersion
        if grating == 400:
            dispersion = 1
        if grating == 600:
            dispersion = 0.65
        if grating == 930:
            dispersion = 0.42
        if grating == 2400:
            dispersion = 0.12
        if grating == 2100:
            dispersion = 0.15
        if grating == 1800:
            dispersion = 0.19
        
        #calculatin Texp
        flux   = 1.0                                                  #by default in units e-/seg/A
        Tslit  = special.erf((slit/2)/(np.sqrt(2)*seeing/2.35482))
        counts = flux*dispersion*Tslit*10**(0.4*(18.8 - v))           #(e-/pixel/sec)
        Texp   = snr**2/counts
        
        return Texp
    
    @classmethod
    #based on MIKE manual: http://www.ucolick.org/~rab/MIKE/usersguide.html
    #see user's manual: http://www.lco.cl/?epkb_post_type_1=the-mike-magellan-inamori-kyocera-echelle-users-guide
    def texp_mike(self, v, snr, seeing, slit, bin=None):
        import scipy.stats
        
        if not bin:
            bin = 'no'
        else:
            bin = 'yes'
            
        if bin == 'no':
            #no spectral binning, only use with <=0.5"
            light_in_slit = scipy.stats.norm(0, seeing).cdf(slit) - scipy.stats.norm(0, seeing).cdf(-slit)
            Texp          = (((snr/50.0)**2)*2500.0)/(1*0.05*10.0**(0.4*(18.4-v))*light_in_slit)
        
        if bin == 'yes':
            #two-pixel spectral binning, use with >0.5" 
            light_in_slit = scipy.stats.norm(0, seeing).cdf(slit) - scipy.stats.norm(0, seeing).cdf(-slit)
            Texp          = (((snr/20.0)**2)*400.0)/(1*0.1*10.0**(0.4*(18.4-v))*light_in_slit)
        
        return Texp
    
    @classmethod
    #from deg to hmsdms
    def de_to_se(self, ra, dec):
        import pandas as pd
        
        c = SkyCoord(ra=ra*u.degree, dec=dec*u.degree)
        d = c.to_string('hmsdms', sep=':')#, format='unicode')
        
        raf = []
        decf = []
        for i in d:
            coor = i.split(sep=' ')
            raf.append(coor[0])
            decf.append(coor[1])
            #print (coor[0])
        result = pd.DataFrame({'ra':raf, 'dec':decf})
        result.to_csv('coord_hmsdms.csv', index=False)
        
    @classmethod
    #scaling Texp
    def esc_texp(self, v1, texp1, v2):
        
        deltV = v1-v2 
        Texpf = texp1/pow(2.5, deltV)
        
        return Texpf


