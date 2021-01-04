'''
Created on Aug 17, 2020

@author: amir
'''
import os
import time
# gunzip
print('installing gunzip')
cmd0 = 'sudo apt update'
os.system(cmd0)
cmd1 = 'sudo apt install gzip'
os.system(cmd1)
print('changing the owner of prerequisits')
script_path = os.path.dirname(os.path.realpath(__file__))
CRNX_path = '{}/RNXCMP/bin/CRX2RNX'.format(script_path)
cmd = 'chmod +x CRX_path'
print('changing the owner of prerequisits')
time.sleep(10)
#chmod +x RNXCMP/bin/CRX2RNX