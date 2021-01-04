'''
Created on Dec 4, 2020

@author: amir
'''
from astropy import coordinates as coord
from astropy import units as u
from astropy.time import Time

class transform():
    def __init__(self):
        pass
        # position of satellite in GCRS or J20000 ECI:
        
        
    def ECI2ECEF(self,epoch:str, coord_:tuple):
        t = Time(epoch)
        cartrep = coord.CartesianRepresentation(x=coord_[0], 
                                                y=coord_[1],
                                                z=coord_[2], unit=u.meter)
        gcrs = coord.GCRS(cartrep, obstime=t)
        itrs = gcrs.transform_to(coord.ITRS(obstime=t))
        loc = coord.EarthLocation(*itrs.cartesian.xyz )
        print(loc.lat, loc.lon, loc.height)
        
        
if __name__ == '__main__':
    tr = transform(epoch= '2017-09-27 12:22:00', coord_=(5713846.540659178, 3298890.8383577876, 0))
    