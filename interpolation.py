'''
Created on Oct 29, 2020

@author: amir
'''
import os,sys
try:
    command = '"/usr/bin/python3.8" "/home/amir/.p2/pool/plugins/org.python.pydev.core_7.7.0.202008021154/pysrc/setup_cython.py" build_ext --inplace'
    os.system(command)
except:
    pass

import rino.rinexReader as rnxrdr
from rino.rinexReader import __rinex__,__sp3__
import datetime
import numpy as np
from astropy.time import TimeDelta
from astropy.time import Time
from _dbus_bindings import Dictionary
import matplotlib.pyplot as plt
import multiprocessing
import warnings
from functools import partial
import multiprocessing as mp
import time

warnings.simplefilter('ignore', np.RankWarning) # Ignore polyfit warnings


global val_epochs, val_epochs_secs, output, window, fit_epochs_secs, rinex_epochs, epochs

class __sp3Interp__():
    def __init__(self):
        self.output = {'rows':[],'columns':[],'depth':[],'tensor':[]}
        self.general_output = {'rows':[],'columns':[],'depth':[],'tensor':[]}
        self.parallel_output = {}
    
    def merge_days(self,day1,day2):
        '''
        Merges tensors and epochs(depth) of day1 and day2 for a single system dataset 
        returns the corresponding dictionary of the merged data
        inputs:
            day1 : the tensor dict of first day
            day2 : the tensor dict of second day (to be piled on the other day)
        output:
            tensor dict of the piled data
        '''
        r,w,c = day1['tensor'].shape 
        c1 = len(day1['depth'])
        c2 = len(day2['depth'])
        c = c1+c2
        
        epochs = day1['depth']
        
        epochs.extend(day2['depth'])
        
        output = {'tensor':np.zeros((r,w,c)),
                  'rows':day1['rows'],
                  'columns':day1['columns'],
                  'depth': epochs,
                  'help': day1['help']}
        
        output['tensor'][:,:,:c1] = day1['tensor']
        output['tensor'][:,:,c1:] = day2['tensor'] 
         
        return output
    
    def read_days(self, year, month, day, n_days = 1, email = 'e189904@metu.edu.tr'):
        '''
        downloads and reads the rinex and sp3 files fo the days of interest and merges them
        returns corresponding dictionary of the sp3 and rinex
        The dictionaries iclude keys as below:
            1- 'tensor': 3d tensor of the data with 
            {1st dimention as rows/flags, 
            2nd dimention as columns/satellite prns 
            and 3rd dimention as depth/epochs}
            2- 'rows': flags for 
        '''
        
        # read rinex
        dtime = datetime.datetime.strptime('{}-{}-{}'.format(year,month,day),"%Y-%m-%d")
        
        rinex = __rinex__()
        gps_sp3 = __sp3__()
        
        rnxData = [None]*n_days
        gps_sp3Data = [None]*n_days
        glonass_sp3Data = [None]*n_days
        
        for day_ in range(n_days):
            d = dtime + datetime.timedelta(days = day_)
            
            # read rinex
            url = rinex.generate_url(leo = 'champ', date = (d.year,d.month,d.day))
            rinex.download(url=url, email = email)
            
            file_path = '{}/{}'.format(rinex.download_path,url.split(sep='/')[-1])
            rinex.unzip(file_path)
            
            rinex.hatanaka2rinex()
            
            rinex.read(rinex.path)
            
            rnxData[day_] = rinex.data.observations
            
            rinex_interval = rinex.data.header['interval']
            
            # read sp3 for gps
            gpsurl = gps_sp3.generate_url('gps', date = (d.year,d.month,d.day))
            
            gps_sp3.download(gpsurl, email=email)
            
            file_path = '{}/{}'.format(gps_sp3.download_path,gpsurl.split(sep='/')[-1])
            
            gps_sp3.unzip(file_path)
            
            gps_sp3Data[day_] = gps_sp3.read(gps_sp3.path, return_type='tensor')
            
            # read sp3 for glonass
            glonass_sp3 = __sp3__()
            
            glonassurl = glonass_sp3.generate_url('glonass', date = (d.year,d.month,d.day))
            
            glonass_sp3.download(glonassurl, email=email)
            
            file_path = '{}/{}'.format(glonass_sp3.download_path,glonassurl.split(sep='/')[-1])
            
            glonass_sp3.unzip(file_path)
            
            glonass_sp3Data[day_] = glonass_sp3.read(glonass_sp3.path, return_type='tensor')
            
        for i,d_rnx in enumerate(rnxData):
            
            d_gps_sp3 = gps_sp3Data[i]
            
            d_glonass_sp3 = glonass_sp3Data[i]
            
            if i == 0:
                rinex_ = d_rnx
                gps_sp3_ = d_gps_sp3
                glonass_sp3_ = d_glonass_sp3
            else:
                rinex_ = self.merge_days(rinex_, d_rnx)
                gps_sp3_ = self.merge_days(gps_sp3_, d_gps_sp3)
                glonass_sp3_ = self.merge_days(glonass_sp3_, d_glonass_sp3)
    
        return rinex_,gps_sp3_, glonass_sp3_, rinex_interval
    
    def get_epochs(self, yesterday,today,tomorrow, init_arg, end_arg):
        '''
        inputs: 
            the data Dictionaries of the current date/dates 
            (if the date of interest is more than one day the today dictonary needs to include all days of interest)
            dates may be constructed using function pile dates
            - yesterday: is the tensor dict of the day prior to days/day of interest
            - today: is the tensor dict of the day/days of interest
                     today mey include the tenosr dicts of multiple days
            -tomrrow: is the tensor dict of a day after the day of interest or in case of multiple days of interest\
                      it includes the tensor dict of the day after the last day of interest   
        output (list):
            - the piled epochs of the date/s of interest including 3 hours befor and after the edges of the epochs
            - the number of added epochs to the beginnings (d_y) and end (d_t) of the array
        '''
        
        interpolation_epochs = yesterday['depth'][init_arg:]
        
        t = tomorrow['depth'][:end_arg]
        
        interpolation_epochs.extend(today['depth'])
        
        interpolation_epochs.extend(t)
        
        s_y = len(yesterday['depth'])
        
        s_t = len(t)
        
        return interpolation_epochs, s_y, s_t
    
    def get_previous_and_next_day(self, sp3, system):
        '''
        gets the sp3 and rinex of nex and previous dates
        inpute:
            sp3 tensor dictionary
            system: gps or glonass
        '''
        # get previous day (one day before first date)
        sp3_ = __sp3__()
        
        e = gps_sp3['depth'][0].iso.split(' ')[0].split('-')
        
        year = int(e[0])
        month = int(e[1])
        day = int(e[2])
        
        y = (datetime.date(year,month,day) - datetime.timedelta(1))
        
        gpsurl = sp3_.generate_url(system, date = (y.year,y.month,y.day))
        
        sp3_.download(gpsurl, email=email)
        
        file_path = '{}/{}'.format(sp3_.download_path,gpsurl.split(sep='/')[-1])
        
        sp3_.unzip(file_path)
        
        yesterday = sp3_.read(sp3_.path, return_type='tensor')
        
        # get one day later (one day after last date)
        e = gps_sp3['depth'][-1].iso.split(' ')[0].split('-')
        
        year = int(e[0])
        month = int(e[1])
        day = int(e[2])
        
        t = (datetime.date(year,month,day) + datetime.timedelta(1))
        
        gpsurl = sp3_.generate_url(system, date = (t.year,t.month,t.day))
        
        sp3_.download(gpsurl, email=email)
        
        file_path = '{}/{}'.format(sp3_.download_path,gpsurl.split(sep='/')[-1])
        
        sp3_.unzip(file_path)
        
        tomorrow = sp3_.read(sp3_.path, return_type='tensor')
        
        yesterday_epochs = [t.mjd for t in yesterday['depth']]
        tomorrow_epochs = [t.mjd for t in tomorrow['depth']]
        
        dt = TimeDelta(3*3600,format='sec')
        init_epoch =(sp3['depth'][0] -dt).mjd 
        end_epoch = (sp3['depth'][-1] + dt).mjd
        
        # get arguments of 3 hurs prior of initial epoch of current day
        init_arg = np.where(yesterday_epochs == init_epoch)[0][0]
        
        # get arguments of 3 hurs after of last epoch of current day
        end_arg = np.where(tomorrow_epochs == end_epoch)[0][0]
        
        # interpolatiopn epochs
        interpolation_epochs,s_y,s_t = self.get_epochs(yesterday, sp3, tomorrow, init_arg, end_arg)
        
        # interpolation tensor (including 3 hurs befor and after
        r,w,c = sp3['tensor'].shape
        c = len(interpolation_epochs)
        interpolation_tensor = np.zeros((r,w,c))
        
        interpolation_tensor[:,:,:s_y-init_arg] = yesterday['tensor'][:,:,init_arg:]
        interpolation_tensor[:,:,s_y-init_arg:-s_t] = sp3['tensor'][:,:,:]
        interpolation_tensor[:,:,-s_t:] = tomorrow['tensor'][:,:,:end_arg]
        interpolation_dict = {'tensor': interpolation_tensor, 
                             'depth': interpolation_epochs,
                             'rows': sp3['rows'],
                             'columns': sp3['columns'],
                             'help': sp3['help']}
        
        return interpolation_dict
    
    def prepare_sp3_for_interpolation(self, gps_sp3, glonass_sp3):
        start_epoch = rinex['depth'][0]
        
        if not start_epoch in gps_sp3['depth']:
            mjds = [t.mjd for t in gps_sp3['depth']]
            a = np.ones_like(mjds)*start_epoch.mjd
            a -= mjds
            a = [abs(i) for i in a]
            arg = np.argmin(a)
            if arg != 0:
                arg-=1
        else:
            arg = np.where(gps_sp3['depth'] == start_epoch)[0][0]
            
        gps_sp3['tensor'] = gps_sp3['tensor'][:,:,arg:]
        gps_sp3['depth'] = gps_sp3['depth'][arg:]
        
        if not start_epoch in glonass_sp3['depth']:
            mjds = [t.mjd for t in glonass_sp3['depth']]
            a = np.ones_like(mjds)*start_epoch.mjd
            a -= mjds
            a = [abs(i) for i in a]
            arg = np.argmin(a)
            if arg != 0:
                arg-=1
        else:
            arg = np.where(glonass_sp3['depth'] == start_epoch)[0][0]
            
        glonass_sp3['tensor'] = glonass_sp3['tensor'][:,:,arg:]
        glonass_sp3['depth'] = glonass_sp3['depth'][arg:]
    
        # get one day ago and next day **gps**
        gps_interpolation_dict = self.get_previous_and_next_day(gps_sp3, system='gps')
        glonass_interpolation_dict = self.get_previous_and_next_day(glonass_sp3, system='glonass')
        
        return gps_interpolation_dict, glonass_interpolation_dict
    
    def prepare_parallel_loop(self, sp3):
        n_cpu = multiprocessing.cpu_count()
        n_p1 = len(sp3['columns']) # number of columns in sp3 file (x,y,z,dt)
        if n_p1 >= n_cpu:
            n_p1 = n_cpu
            n_p2 = 1
        else:
            n_p2 = n_cpu-n_p1
        return n_p1,n_p2
        
    def get_bad_epochs(self, T):
        bad_epoch_ind = np.array(np.where(T >= 999999.99)[0])
        #print('bad epochs: ', bad_epoch_ind)
        if bad_epoch_ind.size == 0: return []
        return bad_epoch_ind
    
    def remove_bad_epochs(self, T,bad_epochs):
        return np.delete(T,bad_epochs)
    
    def find_nearest(self, array, value):
        array = np.asarray(array)
        idx = (np.abs(array - value)).argmin()
        return array[idx],idx

    def __epoch_interpolate__(self, epoch):
        '''
        interpolates sp3 form start epoch to end epoch
        suitable to use in offline processes.
        if you need to apply a real time algorithm __interpolate_rt__ function is suggested
        
        inputs: 
            sp3_: sp3 tensor dictionary
            rinex: rinex tensor dictionary
            rinex_interval: the proper epoch interval (generaly should be selected as the interval of the rinex observations)
            
        outputs:
            interpolated sp3 tensor dictionary
            
        Note: 
            Tensor dictionary include 4 keys as below:
                'tensor': include a 3d matrix by axises indicating 
                    (for sp3: satellite prns as first dimention, coordinate flags [x,y,z or dt] as second dimention, and epochs as third dim)
                    (for rinex: satellite prns as first dimention, observation flags [L1,L2,C1,C2 .etc.] as second dimention, and epochs as third dim)
                'rows': prns
                'columns': flags (sp3: [x,y,z or dt], rinex: [L1,L2,C1,C2 .etc.]
                'depth': epochs
        '''
        print(epoch)
        #self.parallel_output[epoch.value] = {}
        
        if epoch in self.output['depth']: return
        
        dt_indice = np.where(np.array(self.sp3['columns'])=='dt')[0][0]# indice of dt data in sp3
        
        #repochs = Time(Time(rinex['depth']), format = 'mjd')
        epochs = Time(Time(self.sp3['depth']), format = 'mjd')
        
        val,arg = self.find_nearest(epochs.value, value = epoch.mjd)
        if epoch.mjd-val<0:# epoch is before the earest epoch in sp3
            a,b = 9,8#arguments for epoch range selection
            c,d = arg-1,arg#validity epochs for interpolation
        else:
            # epoch is after the earest epoch in sp3
            a,b = 8,9#arguments for epoch selection
            c,d = arg,arg+1#validity epochs for interpolation
    
        val_epoch0,val_epochn = self.sp3['depth'][c],self.sp3['depth'][d]
        # mjds in seconds (interpolant time axis)
        interp_epochs = np.multiply(epochs[arg-a:arg+b].value,3600*24)
        # mjds in seconds (interpolation time axis)
        
        repochs = self.rinex['depth']
        intersection = np.intersect1d(np.where(self.rinex['depth']<=val_epochn), np.where(self.rinex['depth']>=val_epoch0))
        output_epochs = Time(rinex['depth'])[intersection]
        #output_epochs = np.arange(val_epoch0.mjd*3600*24,val_epochn.mjd*3600*24+rinex_interval,rinex_interval)
        
        # reinitialize epochs to prevent ill conditioning
        t_ = np.subtract(interp_epochs,interp_epochs[0])# reinitialize to prevent ill conditioning
        T_ = Time(Time(output_epochs),format = 'mjd')
        T = np.multiply(T_.value,3600*24)
        T = np.subtract(T,interp_epochs[0])
        
        output_tensor = np.zeros((len(self.sp3['rows']),len(self.sp3['columns']),len(output_epochs)))
        self.output = {'tensor':output_tensor, 'columns' : self.sp3['columns'], 'rows':self.sp3['rows'], 'depth':output_epochs}
        
        if not epoch.mjd*3600*24 in T_.value*3600*24:
            raise Warning('correct the output time axis based on rinex time interval')
        
        self.pins = {'sp3':self.sp3,'t':t_,'ind':dt_indice,'arg':arg,'a':a,'b':b,'T':T, 'output':self.output}
        #     list_ = get_interpolation(sp3_['rows'][4], pins)
        #self.get_interpolation(sp3_['rows'][1])
        #get_interpolation_ = partial(self.get_interpolation,pins=self.pins)
        #results = pool.map(howmany_within_range_rowonly, [row for row in data])
        
        
#         pool = mp.Pool(mp.cpu_count())
#         start = time.time()
#         pool.map(self.get_interpolation,[sat for sat,_ in enumerate(self.sp3['rows'])])
#         print(time.time()-start)
#         start = time.time()
        [self.get_interpolation(sat) for sat in self.sp3['rows']]
#         print(time.time()-start)
        
        self.parallel_output[epoch.value]= self.output
        
        #self.output = {'rows':[],'columns':[],'depth':[],'tensor':[]}
        #self.parallel_output[epoch.value]['tensor'] = self.output['tensor']
        
#         self.general_output['depth'].extend(list(self.output['depth'].value))
#         try:
#             #np.concatenate(self.general_output['tensor'],self.output['tensor'], axis = None)
#             self.general_output['tensor'] = np.dstack((self.general_output['tensor'],self.output['tensor']))
#         except:
#             self.general_output['tensor'] = self.output['tensor']
         
    def get_interpolation(self, sat):
        if not all(self.pins['sp3']['tensor'][sat,0,:]==0): 
            bad_indices = self.get_bad_epochs(self.pins['sp3']['tensor'][sat,self.pins['ind'],self.pins['arg']-self.pins['a']:self.pins['arg']+self.pins['b']])
            t = self.remove_bad_epochs(self.pins['t'], bad_indices)
            #self.output = np.zeros((len(self.pins['sp3']['columns']),len(t)))
            
            if not len(t)<16:
                #print(bad_indices)
                for i,xyz in enumerate(self.pins['sp3']['columns']):
                    #print(xyz, len(t), sat)
                    X = self.remove_bad_epochs(self.pins['sp3']['tensor'][sat,i,self.pins['arg']-self.pins['a']:self.pins['arg']+self.pins['b']], bad_indices)
                    # get interpolarion coefficients
                    poly = np.polyfit(t, X, 16 )
                    # interpolated time series
                    self.output['tensor'][sat,i,:] = np.polyval(poly, self.pins['T'])[:]
          
    def __interpolate__(self, sp3_,rinex, epoch, check_muerical_stability = True, save_figures = False, return_interp_rmse = False):
        '''
        interpolates sp3 form start epoch to end epoch
        suitable to use in offline processes.
        if you need to apply a real time algorithm __interpolate_rt__ function is suggested
        
        inputs: 
            sp3_: sp3 tensor dictionary
            rinex: rinex tensor dictionary
            rinex_interval: the proper epoch interval (generaly should be selected as the interval of the rinex observations)
            
        outputs:
            interpolated sp3 tensor dictionary
            
        Note: 
            Tensor dictionary include 4 keys as below:
                'tensor': include a 3d matrix by axises indicating 
                    (for sp3: satellite prns as first dimention, coordinate flags [x,y,z or dt] as second dimention, and epochs as third dim)
                    (for rinex: satellite prns as first dimention, observation flags [L1,L2,C1,C2 .etc.] as second dimention, and epochs as third dim)
                'rows': prns
                'columns': flags (sp3: [x,y,z or dt], rinex: [L1,L2,C1,C2 .etc.]
                'depth': epochs
        '''
        
        dt_indice = np.where(np.array(sp3_['columns'])=='dt')[0][0]# indice of dt data in sp3
        repochs = Time(Time(rinex['depth']), format = 'mjd')
        epochs = Time(Time(sp3_['depth']), format = 'mjd')
        val,arg = self.find_nearest(epochs.value, value = epoch.mjd)
        
        if epoch.mjd-val<0:# epoch is before the earest epoch in sp3
            a,b = 9,8#arguments for epoch range selection
            c,d = arg-1,arg#validity epochs for interpolation
        else:# epoch is after the earest epoch in sp3
            a,b = 8,9#arguments for epoch selection
            c,d = arg,arg+1#validity epochs for interpolation
    
        val_epoch0,val_epochn = sp3_['depth'][c],sp3_['depth'][d]
        # mjds in seconds (interpolant time axis)
        interp_epochs = np.multiply(epochs[arg-a:arg+b].value,3600*24)
        # mjds in seconds (interpolation time axis)
        
        repochs = rinex['depth']
        intersection = np.intersect1d(np.where(rinex['depth']<=val_epochn), np.where(rinex['depth']>=val_epoch0))
        output_epochs = Time(rinex['depth'])[intersection]
        
        #output_epochs = np.arange(val_epoch0.mjd*3600*24,val_epochn.mjd*3600*24+rinex_interval,rinex_interval)
        # reinitialize epochs to prevent ill conditioning
        t_ = np.subtract(interp_epochs,interp_epochs[0])# reinitialize to prevent ill conditioning
        T_ = Time(Time(output_epochs),format = 'mjd')
        T = np.multiply(T_.value,3600*24)
        T = np.subtract(T,interp_epochs[0])
        
        output_tensor = np.zeros((len(sp3_['rows']),len(sp3_['columns']),len(output_epochs)))
        output = {'tensor':output_tensor, 'columns' : sp3_['columns'], 'rows':sp3_['rows'], 'depth':output_epochs}
        
        if not epoch.mjd*3600*24 in T_.value*3600*24: raise Warning('correct the output time axis based on rinex time interval')
        
        RMSE = {}
        for sat in sp3_['rows']:
            if all(sp3_['tensor'][sat,0,:])==0: continue
            RMSE[sat] = {}
            bad_indices = self.get_bad_epochs(sp3_['tensor'][sat,dt_indice,arg-a:arg+b])
            
            t = self.remove_bad_epochs(t_, bad_indices)
            
            if len(t)<16:continue
            
            print(bad_indices)
            for i,xyz in enumerate(sp3_['columns']):
                print(xyz, len(t), sat)
                X = self.remove_bad_epochs(sp3_['tensor'][sat,i,arg-a:arg+b], bad_indices)
                # get interpolarion coefficients
                poly = np.polyfit(t, X, 16 )
                
                # test coefficients
                if return_interp_rmse:RMSE[sat][xyz] = self.__rmse__(np.polyval(poly, t),X)
                    
                # check the residuals based on machine epsilon (optional)
                if check_muerical_stability:
                    diff = np.subtract(np.polyval(poly, t),X)
                    if max(diff)> np.finfo(np.float).eps*1e11:
                        raise Warning('poly fit is not done properly, \n max(residuals) = {} \n machine epsilon = {}'.format(max(diff),np.finfo(np.float).eps))
                
                # interpolated time series
                output['tensor'][sat,i,:] = np.polyval(poly, T)
        
        if save_figures:self.__save_figures__(output)
        
        if return_interp_rmse: 
            return output, RMSE
        else:
            return output
    
    def __save_figures__(self, input_tensor_dict):
        cwd = os.path.dirname(os.path.realpath(__file__))
        t_ax = input_tensor_dict['depth']     
           
        for sat in input_tensor_dict['rows']:
            if all(input_tensor_dict['tensor'][sat,0,:]==0):continue
            for j,xyz in enumerate(input_tensor_dict['columns']):
                plt.figure()
                plt.xlabel('epochs')
                plt.ylabel(xyz)
                plt.plot(Time(input_tensor_dict['depth'],format = 'mjd').value,input_tensor_dict['tensor'][sat,j,:])
                plt.savefig('{}/output_figures/{}_{}.png'.format(cwd,sat,xyz), bbox_inches = 'tight',dpi=600)
        
    def __rmse__(self, predictions, targets):
        return np.sqrt(((predictions - targets) ** 2).mean())
    
    def paral(self,epoch):
        if not epoch in self.output['depth']:
            self.__epoch_interpolate__(epoch)
    
    def __interpolate_offline__(self, sp3,rinex):
        '''
        interpolates sp3 form start epoch to end epoch
        suitable to use in offline processes.
        if you need to apply a real time algorithm __interpolate_rt__ function is suggested
        
        inputs: 
            sp3_: sp3 tensor dictionary
            rinex: rinex tensor dictionary
            rinex_interval: the proper epoch interval (generaly should be selected as the interval of the rinex observations)
            
        outputs:
            interpolated sp3 tensor dictionary
            
        Note: 
            Tensor dictionary include 4 keys as below:
                'tensor': include a 3d matrix by axises indicating 
                    (for sp3: satellite prns as first dimention, coordinate flags [x,y,z or dt] as second dimention, and epochs as third dim)
                    (for rinex: satellite prns as first dimention, observation flags [L1,L2,C1,C2 .etc.] as second dimention, and epochs as third dim)
                'rows': prns
                'columns': flags (sp3: [x,y,z or dt], rinex: [L1,L2,C1,C2 .etc.]
                'depth': epochs
        '''
        self.rinex = rinex
        self.sp3 = sp3
#         import time
#         start = time.time()
#         [self.paral(epoch) for epoch in rinex['depth']]
#         print(time.time()-start)
        pool = mp.Pool(mp.cpu_count())
        start = time.time()
        pool.map(self.paral,rinex['depth'])
        print(time.time()-start)
        
if __name__ == '__main__':
    sp3 = __sp3Interp__()
    # day of process (if multiple days it is the first day and n_days gives the number of days) 
    year = 2002
    month = 2
    day = 1
    n_days = 1
    
    email = 'e189904@metu.edu.tr'
    rinex,gps_sp3, glonass_sp3, rinex_interval = sp3.read_days(year, month, day, n_days, email)
    gps_interpolation_dict, glonass_interpolation_dict = sp3.prepare_sp3_for_interpolation(gps_sp3, glonass_sp3 )
    sp3.rinex = rinex
    sp3.sp3 = gps_interpolation_dict
    #sp3.__epoch_interpolate__(rinex['depth'][100])
    #sp3.sp3 = glonass_interpolation_dict
    #sp3.__epoch_interpolate__(rinex['depth'][100])
    sp3.__interpolate_offline__(gps_interpolation_dict, rinex)
    sp3.__interpolate_offline__(glonass_interpolation_dict, rinex)