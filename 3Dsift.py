# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 13:42:16 2020

@author: Etenne Pepyn
"""


# -*- coding: utf-8 -*-
"""
Created on Tue Oct 22 16:16:02 2019
import sys
sys.argv=['',r"S:\HCP_NoSkullStrip_T1w\Keypoints_VoxelCoordinates\100307_T1w_NoSkullStrip.key"]
execfile(r'S:\keySkullStripping\Python\visualizeFeatures.py')
@author: Etenne Pepyn
"""
import os
from pathlib import Path
import numpy as np 
import nibabel as nib
import re
from scipy.ndimage import gaussian_filter
import pandas
import subprocess
import scipy.ndimage as ndi
from pyflann import *

class kI:
    """
    column indexes for data in 3D sift keypoint files
    """
    scale=3
    XYZ=slice(0,3,1)
    descriptor=slice(17,81,1)
    flag=16
    noFlag=[x for x in range(81) if x != 16] 
    
#********** Read write utilities ***********************************************
def ReadImage(pFile):
    """
    Read volume image

    Parameters
    ----------
    pFile : path of the volume file

    Returns
    -------
    intensity array and header

    """
    pFile=str(Path(pFile))
    img=nib.load(pFile)
    arr=np.squeeze(img.get_fdata())
    h=img.header
    return[arr,h]    
    
def SaveImage(pFile,arr,header):
    """
    Saves a volume image

    Parameters
    ----------
    pFile : path of the volume image
    arr : intensity array
    header 

    Returns
    -------
    None.

    """
    im2=nib.Nifti1Image(arr,affine=None,header=header)
    nib.save(im2,pFile)
    
def ReadKeypoints(pFile):
    """
    Reads a keypoint file

    Parameters
    ----------
    pFile : path of the keypoint file

    Returns
    -------
    mKey : matrix of keypoint

    """
    pFile=str(Path(pFile))
    file= open (pFile,'r')
    i=0
    end=0
    while end==0:
        r=file.readline()
        if 'Scale-space' in r:
            end=1
        i=i+1
    file.close()
    a=np.arange(0,i)
    listC=list(a)
    df=pandas.read_csv(pFile,sep='\t',header=None,index_col=False,skiprows=listC)
    if len(df.columns)==82: #happens when there are tab at the end of a line (shouldnt be there)
        df=df.drop(81,axis=1)

    mKey=df.to_numpy(dtype=np.float64)

    return mKey

def GetResolutionHeaderFromKeyFile(pFile):
    """
    Get resolution and header information from keypoint file

    Parameters
    ----------
    pFile : keypoint file path

    Returns
    -------
    resolution array of shape (3,) and header as a string

    """
    rReso=re.compile('(\d+ \d+ \d+)')
    
    file= open (pFile,'r')
    #skI
    header=''
    end=0
    while end==0:
        r=file.readline()
        header=header+r
        if 'resolution' in r or 'Resolution' in r:
            resoString=rReso.findall(r)
            if resoString ==[]:
                print('resolution not found in'+pFile)
            else:
                resoString=resoString[0]
            resolution=[int(i) for i in resoString.split()]
        if 'Scale-space' in r:
            end=1
    if r[:5]!='Scale':
        print('ERROR: keypoint format not supported in'+pFile)
    file.close()
    if not('resolution' in locals()):
        resolution=np.array([0,0,0])
    return [resolution,header]

def WriteKeyFile(path,mKey,header='default'):
    """
    Write a key file.

    Parameters
    ----------
    path : path of the key file
    mKey : matrix of keypoint to write
    header : string, when provided no header it writes a token header to keep
    format, avoid this situation.

    Returns
    -------
    None.

    """
    if header=='default':
        header="# featExtract 1.1\n# Extraction Voxel Resolution (ijk) : 176 208 176\
        \n#Extraction Voxel Size (mm)  (ijk) : 1.000000 1.000000 1.000000\
        \n#Feature Coordinate Space: millimeters (gto_xyz)\nFeatures: " +str(mKey.shape[0])+"\
        \nScale-space location[x y z scale] orientation[o11 o12 o13 o21 o22 o23 o31 o32 o32] 2nd moment eigenvalues[e1 e2 e3] info flag[i1] descriptor[d1 .. d64]\n"
    else:
        p=re.compile('Features: \d+')
        header=p.sub('Features: '+str(mKey.shape[0])+' ',header)
    fW=open(path,'w',newline='\n')
    fW.write(header)
    for i in range(mKey.shape[0]):
        for j in range(mKey.shape[1]):
            if j>=16:
                n=str(int(mKey[i,j]))
                fW.write(n+'\t')
            else:
                n=mKey[i,j]
                fW.write(format(n,'.6f')+'\t')
            
        fW.write('\n')
    fW.close()


def ExtractAllKeys(dVolume,dDst,pExecutable=r"S:\75mmHCP\featExtract.exe",volumeID='(\d{6})[_.]'):
    """
    Extract keypoints from all image volumes in dVolume.

    Parameters
    ----------
    dVolume : volume directory
    dDst : destination directory
    pExecutable : path of the featExtract executable
    volumeID : regex volume ID, used to shorten names

    Returns
    -------
    Keypoint files located in dDst

    """
    rVolumeID=re.compile(volumeID)
    allF=os.listdir(dVolume)
    print('extracting keypoints')
    for f in allF:
        num=rVolumeID.findall(f)[0]
        print(num)
        list_files = subprocess.run([pExecutable,os.path.join(dVolume,f),os.path.join(dDst,num+'.key')])
        print("The exit code was: %d" % list_files.returncode)
        


#******************* Filtering keypoint images **********************************
def FilterKeyByScale(mKey,minScale,maxScale):
    """
    Returns a subset of keypoint with scale between min and maxScale from the matrix of keypoint mKey. 

    Parameters
    ----------
    mKey : n keypoints in a matrix of size (n,81)
    minScale : float
    maxScale : float

    Returns
    -------
    k1 : subset of keypoints

    """
    mKey1=mKey[mKey[:,kI.scale]>minScale,:]
    mKey1=mKey1[mKey1[:,kI.scale]<maxScale,:]
    return mKey1
    
def FilterKeyByClass(mKey,classValue):
    """
    Returns keys with flags equal to classValue

    """
    return mKey[np.int32(mKey[:,kI.flag])&15==classValue,:]

def FilterKeyByRotation(mKey,rotation=False):
    """
    Returns keys with the specified rotation
    """
    if rotation==False:
        return mKey[~np.int32(mKey[:,kI.flag])&32==32,:]
    else:
        return mKey[np.int32(mKey[:,kI.flag])&32==32,:]

def FilterKeysWithSingleMask(mKey,mask,distanceRatio=0):
    """
    Filter a set of key to conserve only those positionned inside mask.
    If the distanceRatio parameter is used, 
    keypoints less than distanceRatio*scale from the mask's edge will be 
    removed from the resulting set. 
    
    In this case, the distance is an approximation. For each keypoint, 
    we verify that a cube of side 2*distanceRatio*scale center on keypoint
    is totally inside the mask. If it is not, the keypoint is rejected. 
    

    Parameters
    ----------
    mKey : n keypoints in a matrix of size (n,81)
    mask : binary matrix of the size of the image used to extract keypoints in mKey.
        Contains a mask of the brain.
    distanceRatio : should vary between 0 and 2 most of the time.

    Returns
    -------
    mKey2 : keypoints present in mask

    """
    mKey2=np.zeros(mKey.shape)
    prevXYZ=np.array([1,2,3])
    transfered=0
    for i in range(mKey.shape[0]):
        XYZ=np.int32(mKey[i,kI.XYZ])
        if ~np.allclose(XYZ,prevXYZ): #used to skip keypoints with different rotation
            if distanceRatio==0:
                if mask[tuple(XYZ)]==True:
                    mKey2[i,:]=mKey[i,:]
                    transfered=1
                else:
                    transfered=0
            else:
                c=distanceRatio*int(mKey[i,kI.scale])
                if np.all(mask[XYZ[0]-c:XYZ[0]+c,XYZ[1]-c:XYZ[1]+c,XYZ[2]-c:XYZ[2]+c]):
                    mKey2[i,:]=mKey[i,:]         
                    transfered=1
                else:
                    transfered=0
        elif transfered==1:
            mKey2[i,:]=mKey[i,:]
        prevXYZ=XYZ
    mKey2=mKey2[~np.all(mKey2==0,axis=1)]
    return mKey2
        
def FilterKeyWithShiftedMasks(mKey,mask,distanceRatio=1):   
    """
    Filter a set of key to conserve only those positionned inside mask, or 
    at a specified distance from the mask's border. The distance is
    distanceRatio*scale.
    
    Very similar to FilterKeysWithSingleMask, but always apply a distance ratio.
    
    The distance ratio is better applied in this function, using many 
    erroded masks, one for each scale. By erroding masks with a sphere element
    of radius a, we generate a new mask where every point will be at least a
    away from the original mask's border. 

    Parameters
    ----------
    mKey : n keypoints in a matrix of size (n,81)
    mask : binary matrix of the size of the image used to extract keypoints in mKey.
        Contains a mask of the brain.
    distanceRatio : should vary between 0 and 2 most of the time.

    Returns
    -------
    mKey2: keypoints present in shifted masks

    """
    maxScale=9
    scales=np.zeros(maxScale)
    for i in range(maxScale):
        scales[i]=1.6*np.power(2,(i/3))
    
    #GenerateMasks
    lMask=[]
    for i in  range(maxScale):         
        r=int(np.round(distanceRatio*scales[i]))
        d=1+2*r
        sphereElement=np.zeros((d,d,d))
        drawSphere(sphereElement,r,r,r,r)
        lMask.append(ndi.binary_erosion(mask,structure=sphereElement))
    mKey2=_FilterKeysWithScaleMasks(mKey,lMask,scales)
    return mKey2
                
def _FilterKeysWithScaleMasks(mKey,scaleMasks,scales):
    #used by FilterKeyFileWithShiftedMasks
    mKey2=np.zeros(mKey.shape)
    prevXYZ=np.array([1,2,3])
    transfered=0
    for i in range(mKey.shape[0]):
        XYZ=np.int32(mKey[i,kI.XYZ])
        if ~np.allclose(XYZ,prevXYZ): #used to skip keypoints with different rotation
            scale=mKey[i,kI.scale]
            j=np.argmin(np.abs(scales-scale)) #closest scale mask
            sMask=scaleMasks[j]
            if sMask[tuple(XYZ)]==True:
                mKey2[i,:]=mKey[i,:]
                transfered=1
            else:
                transfered=0
        elif transfered==1:
            mKey2[i,:]=mKey[i,:]
        prevXYZ=XYZ
    mKey2=mKey2[~np.all(mKey2==0,axis=1)]
    return mKey2

#******************* analysing and modifying keypoints
    
def CompareKeyImages(mKey1,mKey2):
    """
    Returns the number of exact keypoints match between 2 set of keypoints

    Parameters
    ----------
    mKey1 : n keypoints in a matrix of size (n,81)
    mKey2 : n keypoints in a matrix of size (n,81)

    Returns
    -------
    s : number of exact keypoints

    """
    s=0
    if mKey1.shape[0]>mKey2.shape[0]:
        for i in range(mKey2.shape[0]):
            if np.sum(np.all(mKey1[:,kI.noFlag]-mKey2[i,kI.noFlag]<1e-005,axis=1))==1:
                s+=1
    else:
        for i in range(mKey1.shape[0]):
            if np.sum(np.all(mKey2[:,kI.noFlag]-mKey1[i,kI.noFlag]<1e-005,axis=1))==1:
                s+=1
    return s     


def SubstractKeyImages(mKeyPositive,mKeyNegative):
    """
    Set difference mKeyPositive-mKeyNegative

    Parameters
    ----------
    mKeyPositive : n keypoints in a matrix of size (n,81)
    mKeyNegative : n keypoints in a matrix of size (n,81)

    Returns
    -------
    out : resulting set

    """
    rest=np.copy(mKeyPositive)
    for i in range(rest.shape[0]):
        if np.sum(np.all(mKeyNegative[:,kI.noFlag]-rest[i,kI.noFlag]<1e-005,axis=1))==1:
            rest[i,:]=0
    out=rest[~np.all(rest==0,axis=1)]
    return out

def FindMatchBetweenBrainKeypoints(mKeySkullStrippedKeypoint,mKeyOriginalKeypoint,maxDistance=5000):
    """
    Produce a list of matching keypoint between 2 keypoints list. This is meant
    as an example of how to use the FLANN library. Those match were used to
    visualize which keypoints were equivalent in an original image and a 
    skull-stripped one. 
    Files contain only non-rotated keypoints. At maximum 5000 distance (of descriptors),
    the maximum difference of position is 4 pixel (measured manuallyon patient number 100206 from HCP)

    Parameters
    ----------
    pSkullStrippedKeypoint : path of skull-stripped keypoints
    pOriginalKeypoint : path of original keypoints
    maxDistance : maximum euclidean distance between descriptors to be considered
    a match.

    Returns
    -------
    matching keypoints`

    """

    kSS=mKeySkullStrippedKeypoint
    kOri=mKeyOriginalKeypoint
    flannOri = FLANN()
    paramsOri = flannOri.build_index(kOri[:,kI.descriptor], algorithm="kdtree",trees=8);
    outputMat=np.zeros((0,2))
    for i in range(kSS.shape[0]):
        idx1, dist1 = flannOri.nn_index(kSS[i,kI.descriptor],1, checks=paramsOri["checks"])
        idx=int(idx1[0])
        dist=dist1[0]
        if dist<maxDistance:
            a=np.array([[i,idx]])
            outputMat=np.concatenate((outputMat,a),axis=0)
        # outputMat[i,0]=dist
        # outputMat[i,1]=np.sqrt(np.sum((kSS[i,kI.XYZ]-kOri[idx,kI.XYZ])**2))
        
    # pandas.DataFrame(outputMat).to_csv(pOut, header=None, index=None)
    return np.int64(outputMat)


#*********** Function to use other funciton more effectively with files *********
    
def ApplyFuncToKeyFile(fct,args):
    """
    Apply a function to a single keypoint file

    example: ApplyFuncToKeyFile(FilterKeyByRotation,[r"S:\funWithKey\100307ori.key",True])
    Parameters
    ----------
    fct : function name
    args : arguments of the function, in the same order, with the keypoint path
    instead of the keypoint matrix

    Returns
    -------
    Modifies the function

    """
    k=ReadKeypoints(args[0])
    r,h=GetResolutionHeaderFromKeyFile(args[0])
    k1=fct(k,*args[1:])
    WriteKeyFile(args[0],k1,h)
    
def ApplyFuncToKeyInDirectory(fct,args,dDst=None):
    """
    Apply function to all keypoint file in a directory. Either modifies the original
    or create a modified version in dDst
    
    example: 
    ApplyFuncToKeyInDirectory(FilterKeyByScale,[r'S:\funWithKey\ori',3,6],r'S:\funWithKey\destination')

    Parameters
    ----------
    fct : function name
    args : arguments of the function, in the same order, with the keypoint directory path
    instead of the keypoint matrix
    dDst : destination directory

    Returns
    -------
    None.

    """
    srcFolder=args[0]
    if dDst==None:
        dDst=srcFolder
    
    allF=os.listdir(srcFolder)
    
    for f in allF:
        if os.path.splitext(f)[1]=='.key': 
            k=ReadKeypoints(os.path.join(srcFolder,f))
            r,h=GetResolutionHeaderFromKeyFile(os.path.join(srcFolder,f))
            k1=fct(k,*args[1:])
            WriteKeyFile(os.path.join(dDst,f),k1,h)
            
def ApplyFuncToDirectoryPair(fct,args,dDst,nameFormat='(.+)\.'):
    """
    Execute functions needing 2 file input. First directory need to contain .key  files.
    Second directory can contain either .key files or masks as medical skull-stripped images.
    
    The matching files from the second directory are found using nameFormat as a
    Regular expression for finding a captured group containing the shared name
    between the 2 file types. This allows for slight difference between the 2 name
    formats. The current Regex will return the name of the file before the '.'.
    
    
    
    example:
    ApplyFuncToDirectoryPair(FilterKeyWithShiftedMasks,[pOriginal,pMask,1.5],dst)

    Parameters
    ----------
    fct : function name
    args : arguments of the function, in the same order, with the keypoint directory path
    instead of the keypoint matrix and a directory path instead of a keypoint matrix 
    or a mask
    dDst: destination directory,
    nameFormat : regular expression for finding the shared name of the 2 file
    types. The default is '(.+)\.'.

    Returns
    -------
    A keypoint file in dDst for each keypoint file in the first directory

    """
    dMaster=args[0]
    dSlave=args[1]
    rName=re.compile(nameFormat)
    allMasterFile=os.listdir(dMaster)
    allSlaveFile=os.listdir(dSlave)
    
    if not os.path.exists(dDst):
        os.makedirs(dDst)
        
    for nameMaster in allMasterFile:
        if os.path.splitext(nameMaster)[1]=='.key': 
            name=rName.findall(nameMaster)[0]
            slaveName=[x for x in allSlaveFile if name in x][0]
            pSlave=os.path.join(dSlave,slaveName)
            arg1=ReadKeypoints(os.path.join(dMaster,nameMaster))
            r,h=GetResolutionHeaderFromKeyFile(os.path.join(dMaster,nameMaster))
            if Path(pSlave).suffix=='.key':
                arg2=ReadKeypoints(pSlave)
            else: #we assume we want a mask
                arg2=ReadImage(pSlave)[0]>0
                    
            k1=fct(arg1,arg2,*args[2:])
            WriteKeyFile(os.path.join(dDst,name+'.key'),k1,h)
        
        
#************ DRAW SPHERE ************************
def fill(array2D,y0):
    for x in range(array2D.shape[0]):
           ind=np.argwhere(array2D[x,:])
           if ind.any():
               d=int(y0-ind[0])
               array2D[x,y0-d:y0+d+1]=True
    
def drawCircle(array, x0, y0, radius):
    #mid-point circle drawing algorithm
    f = 1 - radius
    ddf_x = 1
    ddf_y = -2 * radius
    x = 0
    y = radius
    array[x0, y0 + radius]=True
    array[x0, y0 - radius]=True
    array[x0 + radius, y0]=True
    array[x0 - radius, y0]=True

    while x < y:
        if f >= 0: 
            y -= 1
            ddf_y += 2
            f += ddf_y
        x += 1
        ddf_x += 2
        f += ddf_x    
        array[x0 + x, y0 + y ]=True
        array[x0 - x, y0 + y ]=True
        array[x0 + x, y0 - y ]=True
        array[x0 - x, y0 - y ]=True
        array[x0 + y, y0 + x ]=True
        array[x0 - y, y0 + x ]=True
        array[x0 + y, y0 - x ]=True
        array[x0 - y, y0 - x ]=True
def drawSphere(array,x0,y0,z0,radius):
    r=0
    drawCircle(array[:,:,z0],x0,y0,radius)
    fill(array[:,:,z0],y0)
    for i in range(radius,0,-1):  
        
        while radius+0.5>np.sqrt(r**2 +i**2):
            r+=1
        r-=1
        drawCircle(array[:,:,z0+i],x0,y0,r)
        drawCircle(array[:,:,z0-i],x0,y0,r)
        fill(array[:,:,z0+i],y0)
        fill(array[:,:,z0-i],y0)