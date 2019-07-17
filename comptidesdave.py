#!/usr/bin/env python
import this
from obspy.core import read, UTCDateTime, Stream
from obspy.clients.fdsn.client import Client
from obspy.signal.polarization import particle_motion_odr
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import hilbert
import matplotlib as mpl
#Set font parameters using matplotlib
mpl.rc('font',family='serif')
mpl.rc('font',serif='Times') 
mpl.rc('text', usetex=True)
mpl.rc('font',size=18)


debug = True
sta = 'COR'
chans = '00'
presloc='30'
net = "IU"
stime= UTCDateTime('2019-121T00:00:00.0')
etime= stime + 3.*24.*60.*60.

client=Client()
inv = client.get_stations(network="IU", station=sta, starttime=stime, endtime=etime, channel="LH*", level='response')

if debug:
    print(inv)

ctime = stime
st = Stream()

#while ctime <= etime:
#    string = '/tr1/telemetry_days/IU_' + sta + '/' + str(ctime.year) + '/' + \
#                str(ctime.year) + '_' + str(ctime.julday).zfill(3) + '/*'
#    st += read(string + chans + '_LH*')
#    st += read(string + presloc + '_LDO*')
#    ctime += 24.*60.*60.
st = read('cordata')
#st += client.get_waveforms(net, sta, "00","LH*", stime, etime)
#st += client.get_waveforms(net, sta, "30","LDO", stime, etime)
st.detrend('constant')
st.merge(fill_value=0)
st.decimate(5)
st.decimate(2)
st.decimate(6)
st.decimate(6)
# Convert to velocity
st.attach_response(inv)

# We now have the data and the metadata so we should rotate it and do particle motion

st.rotate(method="->ZNE", inventory=inv)
st.select(channel="LH*").remove_response(output='ACC')
st.filter('bandpass',freqmin=0.03/1000., freqmax=0.04/1000.)
st.taper(0.05)
if debug:
    print(st)
    


fig = plt.figure(1,figsize=(12,12))
plt.subplots_adjust(hspace=0.001)
for idx, chan in enumerate(['LHZ', 'LHN','LHE']):
    stT = st.select(channel= chan)
    t = np.arange(len(stT[0].data))*360./(24.*60.*60.)
    ax1 = plt.subplot(3, 1, idx+1)
    for tr in stT:
        ax1.plot(t, tr.data*10**9, label=tr.stats.location)
    ax1.text(1., .9*max(tr.data*10**9), chan)
    ax2 = ax1.twinx()
    ax2.plot(t, st.select(channel="LDO")[0].data,label='LDO',color='g', alpha=.5)
    ax2.plot(t, np.imag(hilbert(st.select(channel="LDO")[0].data)),label='Hilbert LDO',color='b', alpha=.5)
    ax2.set_ylabel('Pressure')

    if idx == 0:
        plt.title(sta + ' ' + str(stime.year) + ' ' + str(stime.julday).zfill(3))
    #plt.legend()
    if idx+ 1 < 3:
        plt.subplots_adjust(hspace=0.001)

    plt.xlim((min(t),max(t)))
handles, labels = ax1.get_legend_handles_labels()
handles2, labels2 = ax2.get_legend_handles_labels()
handles += handles2
labels += labels2
ax = plt.gca()
fig.legend(handles, labels, loc = 'lower center', ncol = 5, fontsize = 17)
plt.xlabel('Time (days)')  
#plt.show()
plt.savefig(sta + '_' + str(stime.year) + '_' + str(stime.julday).zfill(3) + '.pdf', format='PDF',dpi=400)

            
