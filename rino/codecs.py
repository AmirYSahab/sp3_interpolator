'''
Created on Aug 21, 2020

@author: amir
'''
import datetime, time
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import os
try:
    import astropy
    from astropy.time import Time
except:
    command = 'python3 -m pip install astropy'
    os.system(command)
    import astropy
    from astropy.time import Time

class read_rinex:
    def __init__(self,rinex, return_type = 'tensor'):
        '''
        return_type = 'raw' or 'xDataFrame' or 'tensor'
        '''
        self.rinex = rinex 
        l = rinex[0].split()
        self.version = float(l[0])
        self.decode(return_type = return_type)
        
    def decode(self, return_type):
        if self.version == 2.0:
            self.__20__(self.rinex, return_type = return_type)
        
    def __20__(self, rinex,  return_type):
        for i,r in enumerate(rinex):
            if r.split()[0] == 'END': 
                header_end_idx = i
                break
        
        self.header = self.header_v20_(rinex[0:i])
        self.observations = self.body_v20_(rinex[i+1:], return_type = return_type)
        
    def body_v20_(self, rinex, return_type):
        n = self.header['n_obs']
        obs_flags = self.header['obs_list']# columns
        n_lines = len(rinex)
        self.body = {}
        self.b = []
        i = 0
        SATS = [None]*100
        D = []

        i = 0
        while True:
            if i>=len(rinex): break
            tline = rinex[i].split()
            obs_time = self.time_line(tline)
            self.body[obs_time] = {}
            n_obs = int(tline[7])
            sats = [int(sat) for sat in tline[8:]]
            for sv in sats:
                SATS[sv-1] = sv
                
            l = rinex[i+1]
            l = [float(x) for x in l.split()]
            if len(l) == len(obs_flags)+2:
                    obs_flags.insert(1, 'sig_{}'.format(obs_flags[0]))
                    obs_flags.insert(3, 'sig_{}'.format(obs_flags[2]))
            
            data = {}
                
            for idx, line in enumerate(rinex[i+1:i+1+n_obs]):
                #data[sats[idx]] = {flag:float(x) for flag,x in zip(obs_flags,line.split())}
                data[sats[idx]] = [float(x) for flag,x in zip(obs_flags,line.split())]
            self.body[obs_time] = data
            
            i = i + n_obs+1
        #self.body2 = pd.DataFrame(self.body)#, columns = obs_flags)
        '''
        for t in self.body2:
            print(self.body2[t][1]['L1'])
        '''
        allSats = [sv for sv in SATS if sv != None]
        sv = np.arange(max(allSats)+1)#lines
        #_ = [None]*len(obs_flags)
        #_ = [None]*(max(allSats)+1)
        obs_times = list(self.body.keys()) # depth
        
        if return_type == 'raw':
            observations = [None]*(len(obs_times))

            for idx,(epoch,body) in enumerate(self.body.items()):
                thislayer=[[None]*len(obs_flags)]*len(sv)
                for sat,item in body.items():
                    thislayer[sat-1] = item            
                observations[idx] = thislayer
            return {
                'observations' : observations,
                'epochs' : obs_times,
                'flags' : obs_flags,
                'sv' : sv
                }, self.header
                
        elif return_type == 'xDataFrame':
            
            observations = [None]*(len(obs_times))

            for idx,(epoch,body) in enumerate(self.body.items()):
                thislayer=[[None]*len(obs_flags)]*len(sv)
                for sat,item in body.items():
                    thislayer[sat-1] = item            
                observations[idx] = thislayer
            return xr.Dataset({'data':(['epoch','sv', 'flag'],observations)},
                                       coords = {'epoch' : obs_times,
                                                 'flag' : obs_flags,
                                                 'sv' : sv},
                                        attrs = {'epoch':'date-time','sv':sv, 'flag': obs_flags})
        elif return_type == 'tensor':
            tensor = np.zeros((len(sv),len(obs_flags),len(obs_times)))
            for idx,(epoch,body) in enumerate(self.body.items()):
                for svn, obs in body.items():
                    tensor[svn,:,idx]=obs
#             for i,sat in enumerate(allSats):
#                 fig, axs = plt.subplots(7, 1)
#                 for j,f in enumerate(obs_flags):
#                     d = tensor[i,j,:]
#                     axs[j].plot(d, np.arange(len(d)),'.')
#                 plt.show()

            return {'tensor':tensor, 
                    'rows':sv,
                     'columns':obs_flags,
                     'depth': obs_times,
                     'help':' rows are the svn, columns, the observation type and depth are the epoch'}
            
            
            
            #a = list(observations['data'].loc[observations.epoch[0]])
    def time_line(self, l):
        ll = [int(x) for x in l[:5]]
        ll.append(float(l[5]))
        ll[5] = float(l[5])
        
        if int(l[0]) > 80:
            dstring = '19{0:02}-{1:02}-{2:02} {3:02}:{4:02}:0{5}'\
                        .format(int(ll[0]),int(ll[1]),int(ll[2]),int(ll[3]),int(ll[4]),ll[5])
        else:
            dstring = '20{0:02}-{1:02}-{2:02} {3:02}:{4:02}:{5}'\
                        .format(int(ll[0]),int(ll[1]),int(ll[2]),int(ll[3]),int(ll[4]),ll[5])
                
        return Time(dstring)#,"%Y-%m-%d %H:%M:%S.%f")
    
    def obs_lines(self, lines):
        pass
    
    def header_v20_(self, rinex):
        return {
            'date' :                self.date(rinex[1]),#l2
            'system':               rinex[2].split()[0],#l3
            'approx_pos' :          [float(x) for x in rinex[6].split()[:3]],#l6
            'anyenna_Delat_H-E-N' : [float(x) for x in rinex[7].split()[:3]],#l7
            'Wavelength' :          rinex[8].split()[0],#l8_0
            'wlFactor' :            rinex[8].split()[1],#l8_1
            'L1/2' :                rinex[8].split()[2],#l8_2
            'n_obs' :               int(rinex[9].split()[0]),#l9_0    
            'obs_list' :            rinex[9].split()[1 : 1 + int(rinex[9].split()[0])],#l9
            'interval' :            float(rinex[10].split()[0]),#l10
            'firstObsTime' :        self.firstObsTime(rinex[11]),#l11
            'Voltage SNR is mapped to signal strength':rinex[12].split()[-2],#l12
            'SNRV' :                rinex[13].split()[1:-1],#l13
            'sig' :                 rinex[14].split()[2:-1],#l14
            }
        
    def firstObsTime(self, line):
        l = line.split('/')[0].split()
        dstring = '{}-{}-{} {}:{}:{}'.format(l[0],l[1],l[2],l[3],l[4],l[5])
        #return datetime.datetime.strptime(dstring,"%Y-%m-%d %H:%M:%S.%f")
        return Time(dstring) 
    
    def date(self,line): 
        dline = line.split('/')[0].split() 
        dstring = '{} {}'.format(dline[-3],dline[-2])
        #return datetime.datetime.strptime(dstring,"%Y-%m-%d %H:%M:%S")
        return Time(dstring)
    
    def system(self, line):
        l = line.split()

    def __30__(self, rinex):
        pass
    
    def __301__(self, rinex):
        pass
    
    def __304__(self, rinex):
        pass
    
class read_sp3:
    def __init__(self,sp3, return_type = 'tensor'):
        self.sp3 = sp3 
        self.data = self.decode(return_type = return_type)

    def decode(self, return_type):
        
        for first_line_Number,l in enumerate(self.sp3):
            if l[0]=='*': 
                l = self.sp3[first_line_Number+1]
                break
        obs_flags = ['x','y','z','dt']
        #times = [datetime.datetime.strptime('{}-{}-{} {}:{}:{}'.format(int(l[3:7]),int(l[8:10]),int(l[11:13]),int(l[14:16]),int(l[17:19]),float(l[20:31])),"%Y-%m-%d %H:%M:%S.%f") for l in self.sp3 if l[0]=='*']
        times = [Time('{}-{}-{} {}:{}:{}'.format(int(l[3:7]),int(l[8:10]),int(l[11:13]),int(l[14:16]),int(l[17:19]),float(l[20:31]))) for l in self.sp3 if l[0]=='*']
        svns = np.unique([int(i[2:4]) for i in self.sp3 if i[0]=='P' or i[0]=='p'])
        allSats = np.arange(max(svns)+1)
        n_depth = len(times)
        n_width = len(l.split()[2:])
        n_height = len(allSats)
        prns = []
        tensor = np.zeros((n_height,n_width,n_depth))
        epoch = -1
        for l in self.sp3:
            if l[0]=='P':
                sat = int(l[2:4])
                #print(sat,epoch)
                prns.append(sat)
                d = l[5:].split()
                tensor[sat,:,epoch] = [float(d[0]),float(d[1]),float(d[2]),float(d[3])]
            elif l[0]=='*':
                epoch += 1

        if return_type == 'tensor' or return_type == 'raw':
            return {'tensor':tensor, 
                    'rows':allSats,
                     'columns':obs_flags,
                     'depth': times,
                     'help':' rows are the svn, columns, the observation type and depth are the epoch'}
            
        elif return_type == 'xDataFrame':
            return xr.Dataset({'data':(['epoch','sv', 'flag'],tensor)},
                                       coords = {'epoch' : times,
                                                 'flag' : obs_flags,
                                                 'sv' : svns},
                                        attrs = {'epoch':'date-time','sv':svns, 'flag': obs_flags})
                
            