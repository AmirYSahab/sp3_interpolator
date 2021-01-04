import sched, time
from rino import rinexReader
s = sched.scheduler(time.time, time.sleep)

# def do_something(sc): 
#     print("Doing stuff...")
#     # do your stuff
#     s.enter(2, 1, do_something, (sc,))
#     print(time.time())
# 
# s.enter(5, 1, do_something, (s,))
# s.run()


class cast(object):

    def __init__(self, interval, data):
        '''
        Constructor
        '''
        self.data = data['observations']
        self.epochs = data['epochs']
        self.flags = data['flags']
        self.svn = data['sv']
        self.interval = interval
        
        self.epoch_no = 0
        #self.do_something(0)
    
    def do_something(self, sc): 
        print("Doing stuff...")
    # do your stuff
        s.enter(1, 1, self.do_something, (sc,))
        print(time.time())
        
    def stream(self, sc): 
        epoch = self.epochs[self.epoch_no]
        print(epoch)
        print(self.svn)
        print(self.flags)
        print(self.data[self.epoch_no])
        self.epoch_no+=1

        
        # do your stuff
        s.enter(1, 1, self.stream, (sc,))
        
    
    
if __name__ == '__main__':
    rinex = rinexReader.__rinex__()
    url = rinex.generate_url(leo = 'champ', date = (2002,2,1))
    rinex.download(url)
    file_path = '{}/{}'.format(rinex.download_path,url.split(sep='/')[-1])
    rinex.unzip(file_path)
    rinex.hatanaka2rinex(rinex.hpath)
    rinex.read(rinex.path, return_raw = True)
    cast = cast(rinex.header['interval'],rinex.data)
#     
    s.enter(0, 1, cast.stream, (s,))
    data = s.run()
    print(data)

    