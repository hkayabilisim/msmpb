from msmpb_util import *
import pandas as pd
import json
import matplotlib.pyplot as plt
import urllib.request
import tarfile
import hashlib

def md5(file):
    md5_object = hashlib.md5()
    block_size = 128 * md5_object.block_size
    a_file = open(file, 'rb')
    chunk = a_file.read(block_size)
    while chunk:
        md5_object.update(chunk)
        chunk = a_file.read(block_size)
    md5_hash = md5_object.hexdigest()

    return md5_hash

url = 'http://timewarping.org/waterbox.tgz'
downloaded_filename = 'waterbox.tgz'

if md5(downloaded_filename)=="e3778af8fc53587660e97487449357f4":
    print("%s already downloaded"%downloaded_filename)
else:
    print("Downloading and extracting %s" % url)
    urllib.request.urlretrieve(url, downloaded_filename)
    tar = tarfile.open(downloaded_filename, "r:gz")
    tar.extractall()
    tar.close()

experiment = {'data': {'name': 'Box44', 'N': 65850, 'A': 88},
       'trial': [ 
            {'nu':4, 'M': 90, 'a0': 9, 'msmout': [], 'pmeout': []}, 
            {'nu':4, 'M': 84, 'a0': 9, 'msmout': [], 'pmeout': []},
            {'nu':4, 'M': 80, 'a0': 9, 'msmout': [], 'pmeout': []},
            {'nu':4, 'M': 72, 'a0': 9, 'msmout': [], 'pmeout': []},
            {'nu':4, 'M': 70, 'a0': 9, 'msmout': [], 'pmeout': []},  
            {'nu':6, 'M': 60, 'a0': 9, 'msmout': [], 'pmeout': []}, 
            {'nu':6, 'M': 56, 'a0': 9, 'msmout': [], 'pmeout': []},
            {'nu':6, 'M': 54, 'a0': 9, 'msmout': [], 'pmeout': []},   
            {'nu':6, 'M': 50, 'a0': 9, 'msmout': [], 'pmeout': []}, 
            {'nu':6, 'M': 48, 'a0': 9, 'msmout': [], 'pmeout': []},                        
            {'nu':4, 'M': 54, 'a0':12, 'msmout': [], 'pmeout': []}, 
            {'nu':4, 'M': 50, 'a0':12, 'msmout': [], 'pmeout': []},   
            {'nu':4, 'M': 48, 'a0':12, 'msmout': [], 'pmeout': []}, 
            {'nu':4, 'M': 42, 'a0':12, 'msmout': [], 'pmeout': []}, 
            {'nu':4, 'M': 40, 'a0':12, 'msmout': [], 'pmeout': []},   
            {'nu':6, 'M': 42, 'a0':12, 'msmout': [], 'pmeout': []},  
            {'nu':6, 'M': 40, 'a0':12, 'msmout': [], 'pmeout': []}, 
            {'nu':6, 'M': 36, 'a0':12, 'msmout': [], 'pmeout': []},   
            {'nu':6, 'M': 32, 'a0':12, 'msmout': [], 'pmeout': []}, 
            {'nu':6, 'M': 30, 'a0':12, 'msmout': [], 'pmeout': []}]}

for t in experiment['trial']:
    for snapshot in np.arange(1000,101000,1000):
        bincoor='waterbox/%s/waterbox.%s.coor' % (snapshot,snapshot)
        potfile='waterbox/%s/pme.pot' % (snapshot)
        #bincoor='/Users/hkaya/Downloads/waterbox/output/waterbox.%s.coor' % (snapshot)
        #potfile='/Users/hkaya/Google Drive/msm_periodic_paper/outputs/%s/pme.pot' % (snapshot)
        out = runMSM(experiment['data'],t['M'],t['a0'],t['nu'],bincoor=bincoor,potfile=potfile)
        t['msmout'].append(out)
        out = runPME(experiment['data'],t['M'],t['a0'],t['nu'],bincoor=bincoor,potfile=potfile)
        t['pmeout'].append(out)
    
with open('../../results/Figure5.json', 'w') as outfile:
    json.dump(experiment, outfile)
