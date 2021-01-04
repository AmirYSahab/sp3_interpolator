'''
Created on Aug 18, 2020

@author: amir
'''
import datetime, calendar
import os
from astropy.time import Time


class time():
    pass

class gnss_time:
    def __init__(self,input):
        self.time = time()

        self.input = input
        self.script_path = os.path.dirname(os.path.realpath(__file__))
        self.get_astro_time()

    def get_astro_time(self):
        if len(self.input) == 3:
            dstring = '{:04}-{:02}-{:02} 00:00:00'.format(self.input[0],self.input[1],self.input[2])
        else:
            dstring = '{:4}-{:02}-{:02} {:02}:{:02}:{:02}'.format(self.input[0],self.input[1],self.input[2],self.input[3],self.input[4],self.input[5])
        datetimeformat = "%Y-%m-%d %H:%M:%S"
        epoch = datetime.datetime.strptime(dstring,datetimeformat)
        self.time.astroTime = Time(epoch)
        self.time.gpsWeek = self.gps_week(self.time.astroTime.gps)
        self.time.weekDay = self.week_day()
    def gps_week(self, time_in_secs):
        return int(time_in_secs/(3600*24*7))
    
    def week_day(self):
        a =  self.time.astroTime.gps/(3600*24*7) - int(self.time.astroTime.gps/(3600*24*7))
        return int(a*7)
        
    def doy(self):
        return int(self.time.astroTime.yday.split(':')[1])
    
    def jdy(self):
        return self.time.astroTime.jd
    
    def mjd(self):
        return self.time.astroTime.mjd
    
    def iso(self):
        return self.time.astroTime.iso
    
    def isot(self):
        return self.time.astroTime.isot
    
    def gps_secods(self):
        return self.time.astroTime.gps
    
    def time_string(self):
        return str(self.time.astroTime.datetime64)
    
if __name__ == '__main__':
    gnss_time((2017,6,16))
    
    