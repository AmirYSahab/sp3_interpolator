'''
Created on Sep 24, 2020

@author: amir
'''
from datetime import datetime
import scipy.constants

class settings():

    def __init__(self, dataEditing_mode=1):
        self.filter = {
            'ws'    : 'ws/champ_2003_297_307/EKF/',# folder of observations/precise orbir/ephemeris/etc.
            'begin' : self.dt([2003,10,29,0,0,0]),# beginning time of OD (UTC)
            'end'   : self.dt([2003,10,29,1,30,0]),# end time of OD (UTC)
            'type'  : 'kalman',# Type of the filtering algorithm
            'obsType': 'Graphic',# Observation Type
            'Dynmod' : 0,# Set State Vector paremeters that will be estimatedcknowledge receipt of your application for the position of Doctoral candidate
            'measurementNoise' : 0.1,
            'obsSISRE'  : 1.5,
            'stepSize' : 30, # Filter output interval in second
            'timeUpdateThreshold' : 1,
            'enableSaveFilterOutput' :1, # Update threshold for the
            'ephemeris' : 'precise', #'brodcast'
            'EphemerisInterval' : 15*60, # eohemeris modelin interval
            'sigmascale' : 2
            }
        # self.filter.obsType = 'Graphic';                            # Observation Type
        # self.filter.obsType = 'Code';                          # Observation Type
        # self.filter.obsType = 'navsol';                        # Observation Type
        
        # To estimate atmospheric drag, solar radiation presure, empirical
        # accelertaions and markov process corelletion time
        # Set enableDynModParam to "0" or "1" or "2"
        # 0 : do not estimate
        # 1 : estimate atmospheric drag coefficient, solar radiation presure
        #     coefficient and empirical accelertaions
        # 2 : estimate atmospheric drag coefficient, solar radiation presure
        #     coefficient, empirical accelertaions and markov process corelletion
        #     time
        self.init = {
            'CD'              : 2,# initial value for drag acceleration coeff
            'CR'              : 1.5,# initial value for solar radiation acceleration coeff
            'Psrad'           : 4.56e-6,# Solar radiation pressure at 1 AU, [N/m^2], IERS 96
            'empacc'          : [1e-6,1e-6,1e-6],# in RTN
            'corelTime'       : 600,# initial value for correlation time of marcov process
            'posXYZ'          : [0,0,0],
            'velXYZ'          : [0,0,0],#[0.05,0.05,0.05],
            'atmDragCoeff'    : 0.001,
            'solarRadCoeff'   : 0.001,
            'empAccellRTN'    : [1e-9,1e-9,1e-9 ],#[1000e-9 1000e-9 1000e-9]; # m/s2
            'corelTime'       : 1e-2,
            'recClcBias'      : 1,
            'ambiguityBias'   : 0.01,
            }
        
        self.sat = {
            'a2m'   : (1.22/522),# cross sectıonal area to mass (m2/kg)
            'mass'  : 522,# satellite mass (kg)
            'incId' : 4,# satellite orbital parameter representing inclination
            'nmax'  : 70,# maximum degree  of Stoke's coefficients
            }
        # Set Time System Parameters
        self.TimeRefSystem = {
            'UT1_UTC' : -0.3652860,# in second
            'TAI_UTC' : 32,# in second
            'PolarMotion_xp' : 0.220270,# in second
            'PolarMotion_yp' : 0.242220,# in second
            'MJD2000' : 2451545.0-2400000.5,  # Modified Julian date of J2000
            }
        # Set filter time propagation parameters
        # propagation step size in second
        # Set statistical parameters of measurements and auxiliary parameters
        # Set Initial statistical parameters

        # Set parameters required for data editing
        # Set the data editing mode.The value of 1 for filtSet.dataEditing.mode
        # indicates recursive outlier detection. The value of 2 is for robust filtering
        
        # if recursive outlier detection mode is selected consider the following settings
        self.dataEditing = {}
        if dataEditing_mode == 1:
            self.dataEditing = {'outlierFactor' : 5000,# Outlier Factor
                           'elevationThreshold' : 10, # elevation threshold in degree
                           }
            # if robust mode is selected consider the following settings
        elif dataEditing_mode == 2:
            # Set the level of significance (los) value for the chi-square
            # distribution used to detect faulty measurements. The los can be
            # set to one of following values ;
            # los=> (.995 .990 .975 .950 .900 .100 .050 .025 .010 .005 )
            self.dataEditig = {'chiSquare_loss' : 0.95}
        
        # Minumum number of observation
        self.dataEditing['minNumObs'] = 4
        # GPS recever Antenna Offset from the base location in radial, alongTrack
        # crossTrack directions
        self.dataEditing['AntOffFromSatBaseLoc'] = [-0.4306 , -1.488, 0]
        
        # For Unscented Kalman Filter; sigma-vector parameters
        self.ukfParams = {
        'kappa' : 0,
        'alfa' : 0.1,
        'beta' : 2,
        }
        
        self.constants = {'speed_of_light' : scipy.constants.physical_constants['speed of light in vacuum'][0],
                            'lambda_L1'  : (scipy.constants.physical_constants['speed of light in vacuum'][0]/1575.42e6), # wavelength of L1 carrier
                            'f1' : 1575.42e6, #Hz;
                            'f2' : 1227.6e6, #Hz
                            'earth_Radius' : 6378136.46,# Earth's mean radius(m)
                            'earth_GM' : 3986004.415e8,# Earth's gravity constant, [m^3/s^2],EGM2008
                            'Sun_GM' : 1.32712438e20,# Sun's gravity constant [m^3/s^2]; IAU 1976
                            'Moon_GM' : 4902799059741.11,# Moon's garvity constant
                            'earth_w' : 7.2921158553e-5,# Earth angular velocity, NIMA 1997
                            'AU' : 149597870000.0,# Astronomical unit [m], IAU 1976

                          }

    
    def dt(self, t):
        dstring = '{0:04}-{1:02}-{2:02} {3:02}:{4:02}:0{5}'\
                        .format(t[0],t[1],t[2],t[3],t[4],t[5])
        return datetime.strptime(dstring,"%Y-%m-%d %H:%M:%S")

if __name__ == '__main__':
    settings = settings(dataEditing_mode=1)
    print('done')