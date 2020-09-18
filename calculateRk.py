import numpy as np
import sys
import fct3Dsift as FS
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
import matplotlib as mpl
import time
import os


def GetDistanceFromBorder(brainMask,k):
    """
    k:single keypoint
    """
    xyz=tuple(np.int32(k[FS.kI.XYZ]))
    s=k[FS.kI.scale]
    Idraw=np.zeros(brainMask.shape)

    r=1
    FS.drawSphere(Idraw,xyz[0],xyz[1],xyz[2],r)
    if brainMask[xyz]: #inside the brain
        

        while np.all(np.logical_and(Idraw,brainMask)==Idraw) and r < 3*s and r<20:
            r+=1
            FS.drawSphere(Idraw,xyz[0],xyz[1],xyz[2],r)
        
        if r<20 and r<3*s:
            dst=r-1
        else:
            dst=20

    else:
         while not(np.any(np.logical_and(Idraw,brainMask))):
            r+=1
            FS.drawSphere(Idraw,xyz[0],xyz[1],xyz[2],r)
         dst=-(r-1)
    return dst

def GetRbStats(mKey,brainMask):
    brainSize=brainMask.shape
    mK=mKey
    stats=np.zeros((mK.shape[0],4))
    for i in range(mK.shape[0]):
        k=mK[i,:]
        sigma=k[FS.kI.scale]
        s=int((k[FS.kI.scale]))
        j=np.log2(s/1.6)*3
        sIPlus1=int((1.6*np.power(2,(j+1)/3)))
        xyz=tuple(np.int32(k[FS.kI.XYZ]))

        I0=np.zeros(brainSize)
        I1=np.zeros(brainSize)

        FS.drawSphere(I0,xyz[0],xyz[1],xyz[2],2*s)
        FS.drawSphere(I1,xyz[0],xyz[1],xyz[2],2*sIPlus1)
        sum0t=np.sum(I0)
        sum1t=np.sum(I1)

        I0g=gaussian_filter(I0,s)
        I1g=gaussian_filter(I1,sIPlus1)

        sum0m=np.sum(I0g[brainMask])
        sum1m=np.sum(I1g[brainMask])
        rk=((sum0m/sum0t)+(sum1m/sum1t))/2
        distFromBorder=GetDistanceFromBorder(brainMask,k)
        stats[i,0]=rk
        stats[i,1]=distFromBorder/sigma
        stats[i,2]=distFromBorder
        stats[i,3]=sigma
    return stats

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]


dImg=os.path.join(os.getcwd(),'imageExample') #skull-stripped image
pVol=os.path.join(dImg,'OAS1_0001.hdr')
pKey=os.path.join(dImg,'OAS1_0001.key')
[v,h]=FS.ReadImage(pVol)
mK=FS.ReadKeypoints(pKey)
mK=FS.FilterKeyByRotation(mK) #keypoints with the same rotation have the same Rk
brainSize=v.shape

brainMask=v>0

NumberOfKeypointsToTest=20

stats=GetRbStats(mK[0:NumberOfKeypointsToTest,:],brainMask)
a1=stats[stats[:,2]!=20,:] #distance from border of 20 or more were not calculated to fasten the process
plt.scatter(a1[:,1],a1[:,0],c=a1[:,3],s=(a1[:,3]**1.7)*1,alpha=0.3,cmap='tab20')
cbar=plt.colorbar()
cbar.set_label(r'$\sigma$', rotation=0)
#cbar.ax.get_yaxis().labelpad = 15
plt.title('$R_b$ vs. distance factor for skull-stripped brain keypoints')
plt.xlabel(r'distance factor (distance to border/$\sigma$)')
plt.ylabel('$R_b$')
plt.show()