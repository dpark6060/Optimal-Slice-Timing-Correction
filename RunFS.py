#!/usr/bin/env python
import glob
import os
import sys, traceback
import scipy.signal as sig
import numpy as np
import nibabel as nb
from multiprocessing import Pool
import subprocess as sp
import re

class WorkingFile:
    
    def __init__(self):
        self.FileBase=''
        self.FileName=''
        self.RootDir=''
        self.OutputFile=''
        self.TxtOutput=''
        self.JobDir=''
        self.ErrorOutput=''
        self.StdOutput=''
        self.Fs=0.0
        self.Int=0
        self.Fnew=0.0
        self.Cores=0
        self.mem=2.0
        self.SubmitFiles=[]
        self.InputNii=''
        self.force=False
        self.startSlice=0
        self.Direction=1
        self.TR=0
        self.UseOrder=False
        self.OrderFile=''
        self.Ref=0
        self.UserOutput=''
        



def usage():
    print'USAGE: ./RunFS.py -in <InputFile> -tr <TR> -i <Int> [-out <OutputFile> ] [-s <StartSlice>] [-d <Direction>] [-o <SliceORderFile>] [-r <RefSlice>] [-c <nCores>] [-m <Memory>] [-Force]'
    print '<>:\t User Defined Input'
    print '[]:\t Optional Input\n'
    
    print '-in:\t The input image to perform STC on\n'
    print '-tr:\t Set the TR of the original fMRI data in seconds\n'
    
    print '-i:\t Set the sampling interleave acquisition parameter. This number refers to how far away the next sampled slice is'
    print('\t0: No STC, only run filtering')
    print('\t1: Sequential acquisition - The next slice to be acquired from any slice s is slice s+1')
    print('\t2: Even-Odd interleave (Every Other) - The next slice to be acquired from any slice s is slice s+2')
    print('\tn: Skip n slices interleave - The next slice to be acquired from any slice s is slice s+n\n')

    print '-out:\t Set the Output file.  If you just specify a file name, it will be saved in the same directory as the original Nii file.\n'

    print '-s:\t Set the starting slice - The slice that was acquired first in the sequence.  Default is slice 1, the bottom most slice.'
    print('\tThis starts the interleave from that slice.  If your interleave parameter is "1"')
    print('\t and your starting slice is "3", your slice acquisition sequence will be modled as:')
    print('\t3\n\t5\n\t7\n\t9...\n')
    
    print('-d:\t value 1 or -1.  Set the direction of slice acquisition.')
    print '\t1: implies ascending slice acquisition from starting slice: (3,5,7,9...)'
    print '\t-1: implies descending slice acquisition from starting slice: (7,5,3,...)\n'
    
    print('-o:\t Slice Order File.  If present, all interelave parameters are ignored, and slices are shifted using the slice order file')
    print('\t we refer to the bottom slice in the image as slice 1, not slice 0')
    
    print '-r:\t Set the Reference slice'
    print('\tThis is the slice the data is aligned to.  Default is the first slice\n')
    print('-c:\t the number of cores free on your computer for parallel processing.  The default is one less than the number of cores your computer has.')
    print('\tIt is highly recommended that you use parallel processing to speed up this operation\n')
    print('-m:\t amount of free memory you have, in GB, so the system knows how much load it can give each core. The default is 4\n')
    print('-Force:\t Forces the program to proceed even if it estimates a memory error, or if it detects that STC has already been run.  Use -Force if you need to re-run with new parameters')
    print('\tWARNING: THIS CAN POTENTIALLY CRASH YOUR COMPUTER\n\n\n')

    return()


def MakeT(Int,Ref,NumSlices,Tr,StartSlice,Direction):
# Creates timing sequence for the slices in a given image.
# Int - the interpolation value (1 - sequential, 2 - Odd/Even, 3 - 1,4,7,... 2,5,8,... 3,6,9..., etc)
# Ref - the slice number to align to. 0 indicates the first slice.
# NumSlices - the total number of slices
# Tr - the TR of the data, or the sampling rate.

# Calculate dt: the time interval between any two slices acquired sequentially
    dt=float(Tr)/float(NumSlices)
    IntSeq=[]
 
# Create slice acquisition order sequence
    for i in range(0,Direction*int(Int),Direction*1):
        IntSeq.extend(list(np.mod(range(i+StartSlice,int(Direction*NumSlices+StartSlice),int(Direction*Int)),NumSlices)))    
    IntSeq=np.array(IntSeq)

# Initialize slice timing array
    TimeList=range(len(IntSeq))
    trash,TimeList=zip(*sorted(zip(IntSeq,TimeList)))
    
    TimeList=np.array(TimeList)
    TimeList=(TimeList)*dt
    TimeList=TimeList-TimeList[Ref]
# Zero the acquisition time around the slice of interest    
    return(IntSeq,TimeList)

def FiltShift(Inputs):
    
    Zimg=Inputs[0]
    tShift=float(Inputs[1])    
    Fs=float(Inputs[2])
    Fnew=float(Inputs[3])
    if len(Inputs)>=5:
        stdout=Inputs[4]
    else:
        stdout=''
    if len(Inputs)>=6:
        stderror=Inputs[5]
    else:
        stderror=''
# The main filtershift function, takes a single vector ZTshift
# ZTshift: [Slice Number, Time Shift]
# where Slice Number is the index of the slice to be shifted, and Time Shift is the amount to shift.  This can be positive or negative.

    if not stdout=='':
        sys.stdout=open(stdout,'w')
    if not stderror=='':
        sys.stderror=open(stderror,'w')
        
# Set sampling information  
    Tr= 1./Fs
    Foriginal = Fs # Hz
    #Fnew = np.ceil( #Hz
    StopGain = 60 # -1*dB
    TranWidth = 0.08 # Hz
    print 'TR:{}'.format(Tr)
    print 'Foriginal:{}'.format(Fs)
    print 'Fnew:{}'.format(Fnew)
    print 'tShift:{}'.format(tShift)

# The slice is extracted from the 4D global image and reshaped into a 3D volume, representing the 2D slice and its time series.
    #dbg.write(time.ctime()+': Extracting slice {}\n'.format(Z))
    #dbg.flush()
    base,fname=os.path.split(Zimg)
    
    ext='.exe'
    fnameStrip=fname
    while ext:
        fnameStrip,ext=os.path.splitext(fnameStrip)
        
# Create lowpass Filter
    CutOffFreq=0.2
    SamplingRate=Fnew
    Nyq = SamplingRate/2.0   
    N, beta = sig.kaiserord(StopGain,TranWidth/Nyq)        
    BPF = sig.firwin(N, CutOffFreq, window=('kaiser', beta), pass_zero=True, scale=True, nyq=Nyq)
    
    Sig=np.loadtxt(Zimg)
# The time axis must be the last dimension.
    tdim=len(list(np.shape(Sig)))-1    
    
# Length Time Sig (LTS) is the length of the time series associated with the slice signal.
    LTS=Sig.shape[-1]

# Padding is added to the front end end of the slice, where half the signal is mirrored on each end
# FR - Front Range: the range of indicies to be padded on the front (beginning) of the signal
# FR starts at 1 and ends at LTS/2. (0 is the first index)
# BR - Back Range: the range of indicies to be padded on the back (end) of the signal
# BR starts at LTS-1 and goes to LST/2
    FR=np.array(range(int(round(LTS/2.)),0,-1))
    BR=np.array(range(LTS-1,int(FR.shape[-1])-1,-1))

# Pad the signal, with conditions for each dimension so that all data shapes can be accomodated. 
    if tdim==0:
    # One dimensional case (Single Vector)
    
# The signal to be padded on the front is the values of Sig at indecies FR
# The signal to be padded to the     
        FrontPad=(Sig[FR])  
        BackPad=Sig[BR]

# The length of the padding is stored as LFP
        LFP=FrontPad.shape[-1]
    

    elif tdim==1:
    
        FrontPad=Sig[:,FR]
        BackPad=Sig[:,BR]
    

        LFP=FrontPad.shape[-1]

    elif tdim==2:
    
        FrontPad=Sig[:,:,FR]
        BackPad=Sig[:,:,BR]

        LFP=FrontPad.shape[-1]

    elif tdim==3:
    
        FrontPad=Sig[:,:,:,FR]
        BackPad=Sig[:,:,:,BR]

        LFP=FrontPad.shape[-1]

    else:
        print('Bad Array Dimensions for Padding')

# The padding is added to the signal along the time axis
    Sig=np.concatenate((FrontPad,Sig,BackPad),-1)
    #dbg.write('{}: Slice {}, Done\n'.format(time.ctime(),Z))        
    #dbg.flush()
# Upsampling/interpolation paramaters are calculated
# S: Upsampling Factor, or the number of samples to be added between existing samples
# Dm: the dimensions that the upsampled signal will be
# SS: the dimensions of the padded signal
    S=int(round(Fnew/Foriginal))
    Dm=list(np.shape(Sig))
    SS=Dm
    tdim=len(Dm)-1
    Dm[tdim]=Dm[tdim]*S

# Initialize upsampled signal as zeros
    ZeroSig=np.zeros(Dm)

    k=0

# Create Zero Padded Signal    
    if tdim==0:

    
# Assign every S samples in ZeroSig to values in Sig    
        ZeroSig[::S]=Sig
    elif tdim==1:       
        ZeroSig[:,::S]=Sig
    elif tdim==2:        
        ZeroSig[:,:,::S]=Sig                 
    elif tdim==3:       
        ZeroSig[:,:,:,::S]=Sig
    else:
        print("Bad Array Dimensions")  

# Cleanup Sig as it's no longer needed   
    del Sig    

# Filter the Zero padded signal with the designed filter

    
    ZeroSig = sig.filtfilt(BPF, [1], ZeroSig*float(S), axis=-1,padtype='even', padlen=0)
    

# Calculate new frequency and time parameters for the upsampled signal
    tdim=len(SS)-1
    Fs=Fnew
    Ts=1/Fs

# Initialize a variable the length of the padded signal at the original frequency
    Sig=np.zeros(SS)
    
# Calculate the number of indicies to shift when resampling    
    shift=round(tShift*Fs)
    print 'shift:{}'.format(shift)
    print 'tdim:{}'.format(tdim)
# Shift the Signal
    #dbg.write('{}: Slice {}, Shifting...\n'.format(time.ctime(),Z))        
    #dbg.flush()
    if tdim==0:
    # One Dimensional Case (1D vector)
    
        if shift>0:
    # If the shift is larger than zero
    # Extend the Upsampled signal by repeating the values at the beginning of the signal by the shift amount
    # then resample the signal, starting at index 0, every S indecies, to one S of the end
            Rep=np.tile(ZeroSig[0],shift)
            ZeroSig=np.append(Rep,ZeroSig,-1)
            Sig=ZeroSig[range(0,ZeroSig.shape[-1]-S,S)]       
    
        else:
    # If the Shift is less than zero
    # Extend the Upsampled signal by repeating the values at the end of the signal by the shift amount
    # Then resample the signal, starting at index shift, every s indicies, to the end
            Rep=np.tile(ZeroSig[-1],abs(shift))
            ZeroSig=np.append(ZeroSig,Rep,-1)
            Sig=ZeroSig[range(int(abs(shift)),ZeroSig.shape[-1]-1,S)]
    
    # Crop the signal to remove the padding preformed earlier    
        Sig=Sig[LFP:LFP+LTS]
    
    elif tdim==1:   
        if shift>0:
            Rep=np.tile(np.expand_dims(ZeroSig[:,0],axis=-1),[1,shift])
            #Rep=np.expand_dims(Rep,axis=-1)

            ZeroSig=np.append(Rep,ZeroSig,-1)
            Sig=ZeroSig[:,range(0,ZeroSig.shape[-1]-S,S)]       
        
        else:
            Rep=np.tile(np.expand_dims(ZeroSig[:,-1],axis=-1),[1,abs(shift)])
            #Rep=np.expand_dims(Rep,axis=-1)

            ZeroSig=np.append(ZeroSig,Rep,-1)
            Sig=ZeroSig[:,range(int(abs(shift)),ZeroSig.shape[-1]-1,S)]
            
        Sig=Sig[:,LFP:LFP+LTS]

    elif tdim==2: 

        if shift>0:
        
            Rep=np.tile(np.expand_dims(ZeroSig[:,:,0],axis=-1),[1,1,shift])
            #Rep=np.expand_dims(Rep,axis=-1)
            ZeroSig=np.append(Rep,ZeroSig,-1)
            Sig=ZeroSig[:,:,range(0,ZeroSig.shape[-1]-S,S)]       
        
        else:
            Rep=np.tile(np.expand_dims(ZeroSig[:,:,-1],axis=-1),[1,1,abs(shift)])
            #Rep=np.expand_dims(Rep,axis=-1)

            ZeroSig=np.append(ZeroSig,Rep,-1)
            Sig=ZeroSig[:,:,range(int(abs(shift)),ZeroSig.shape[-1]-1,S)]
        
        Sig=Sig[:,:,LFP:LFP+LTS]

       
    elif tdim==3:   
        if shift>0:
        
            Rep=np.tile(np.expand_dims(ZeroSig[:,:,:,0],axis=-1),[1,1,1,shift])
            #Rep=np.expand_dims(Rep,axis=-1)
            ZeroSig=np.append(Rep,ZeroSig,-1)
            Sig=ZeroSig[:,:,:,range(0,ZeroSig.shape[-1]-S,S)]       
        
        else:
            Rep=np.tile(np.expand_dims(ZeroSig[:,:,:,-1],axis=-1),[1,1,1,shift])
    
            #Rep=np.expand_dims(Rep,axis=-1)
            ZeroSig=np.append(ZeroSig,Rep,-1)
            Sig=ZeroSig[:,:,:,range(int(abs(shift)),ZeroSig.shape[-1]-1,S)]
        
        Sig=Sig[:,:,:,LFP:LFP+LTS]

    else:
    
        print("Bad Array Dimensions")  

    savefile=os.path.join(base,'FS_'+fname)
    np.savetxt(savefile,Sig)
    os.system('chmod 777 {}'.format(savefile))  

    return()

def CreateJobsFromFile(f):

    SubmitFiles=[]
    
    Tr=1.0/f.Fs
    RefSlice=0
    
    Nii=nb.load(f.InputNii)
    NiiData=Nii.get_data()
    per=np.percentile(NiiData,3)
    print per
    lx,ly,lz,lt=NiiData.shape
    
    if f.UseOrder:
        try:
            SliceOrder=np.loadtxt(f.OrderFile)
            SliceOrder=SliceOrder-1
        except:
            print('Error Loading Slice Order Flie {}'.format(f.OrderFile))
        
        if not len(SliceOrder) ==lz:
            print('Slice Order File does not match number of slices in {}'.format(f.InputNii))
            quit()
            
        shift=(TR-TR*1.0/lz)/lz
        SliceShift=SliceOrder*shift
        SliceShift=SliceShift-SliceShift[Start]
    
    else:
        
        SliceOrder, SliceShift = MakeT(f.Int, f.Ref, lz, Tr,f.startSlice, f.Direction)
    
    print('Using Slice Order:\n{}'.format(SliceOrder+1))
    
    GBperCore=f.mem/f.Cores*0.9
    BytePerCore=GBperCore*np.power(1024,3)    
    Nelements=lt*2*(f.Fnew/f.Fs)
    MemTotal=Nelements*31.48

    
    
    
    while MemTotal>BytePerCore:
        f.Cores=f.Cores-1
        if f.Cores<1:
            print ('Need More memory, there WILL (probably) be a memory error.  use -Force to force')
            if f.force==False:
                quit()
            
        GBperCore=f.mem/f.Cores*.9
        BytePerElement=np.dtype(float).itemsize
        BytePerCore=GBperCore*np.power(1024,3)
        
    
    txtLen=int(np.floor(BytePerCore)/(MemTotal))
    ParameterFile=os.path.join(f.TxtOutput,'Parameters.txt')
    Rerun=True
    
    if os.path.exists(ParameterFile):
        parameters=open(ParameterFile,'r')
        line=parameters.read()
        groups=re.search('txtLen:(.*)\n',line)
        OldtxtLen=int(groups.groups()[0])
        groups=re.search('Nelements:(.*)\n',line)
        OldNelements=int(np.round(float(groups.groups()[0])))
        if (OldtxtLen==txtLen and OldNelements==Nelements) and f.force==False:
            Rerun=False
 
    SliceLen=glob.glob(os.path.join(f.TxtOutput,'SliceLen_*.txt'))
    
    
    if (not SliceLen) or Rerun:
        
        
        if Rerun:
            cmd='rm {}/*.txt'.format(f.TxtOutput)
            pr=sp.Popen(cmd,shell=True)
            pr.wait()
            
        for z in range(lz):
            Extract=np.squeeze(NiiData[:,:,z,:])
            Extract=np.reshape(Extract,(-1,lt))
            LexOrig=len(Extract)
            ind=np.where(Extract.min(1)>per)
            Extract=Extract[ind]
            lentxt=os.path.join(f.TxtOutput,'SliceLen_{}.txt'.format(z))
            indtxt=os.path.join(f.TxtOutput,'SliceInd_{}.txt'.format(z))
            
            np.savetxt(lentxt,[LexOrig])
            np.savetxt(indtxt,ind[0])
              
            for ii,inc in enumerate(range(0,len(Extract),txtLen)):
                
                if inc+txtLen>=len(Extract):
                    txtdat=Extract[inc:]
                else:
                    txtdat=Extract[inc:inc+txtLen]
                
                TxtName=os.path.join(f.TxtOutput,'Image_{}_{}.txt'.format(z,ii))
                
                np.savetxt(TxtName,txtdat)
            
                SubmitFiles.append([TxtName,'{}'.format(SliceShift[z]),'{}'.format(f.Fs),'{}'.format(f.Fnew),os.path.join(f.StdOutput,'Image_{}_{}.txt'.format(z,ii)),os.path.join(f.ErrorOutput,'Image_{}_{}.txt'.format(z,ii))])
        
        ParamOut=open(ParameterFile,'w')
        ParamOut.write('txtLen:{}\n'.format(txtLen))
        ParamOut.write('Nelements:{}\n'.format(Nelements))
        ParamOut.close()
           
    else:
        print('{} Text Files Exist, Recreating Jobs'.format(f.FileBase))
        
        
        for z in range(lz):            
            txtfiles=glob.glob(os.path.join(f.TxtOutput,'{}_{}_*.txt'.format('Image',z)))
            for ii,ft in enumerate(txtfiles):
                
                
                TxtName=os.path.join(f.TxtOutput,'{}_{}_{}.txt'.format('Image',z,ii))
                shPath=os.path.join(f.JobDir,'{}_slice{}_{}.sh'.format(f.FileBase,z,ii))              
                if os.path.isfile(shPath):
                    os.remove(shPath)
                
                Output=os.path.join(f.TxtOutput,'FS_{}_{}_{}.txt'.format('Image',z,ii))
                if not os.path.exists(Output):
                
                    SubmitFiles.append([TxtName,'{}'.format(SliceShift[z]),'{}'.format(f.Fs),'{}'.format(f.Fnew),os.path.join(f.StdOutput,'Image_{}_{}.txt'.format(z,ii)),os.path.join(f.ErrorOutput,'Image_{}_{}.txt'.format(z,ii))])

    np.savetxt(os.path.join(f.TxtOutput,'SliceOrder.txt'),SliceOrder)
    np.savetxt(os.path.join(f.TxtOutput,'SliceShift.txt'), SliceShift)
    return(SubmitFiles)

def SetUpWorkingPaths(f):
    InputNii=f.InputNii
    
    if f.UserOutput=='':    
        RootDir,FileName=os.path.split(InputNii)
        
    else:
        RootDir,FileName=os.path.split(f.UserOutput)
        if RootDir=='':
            RootDir,trash=os.path.split(InputNii)


    trash='.exe'
    FileBase=FileName
    while trash:
        FileBase,trash=os.path.splitext(FileBase)
    
    
    
    
    
    OutputFile=os.path.join(RootDir,'{}_FS.nii.gz'.format(FileBase))
    
    
    
    n=1
    
    while os.path.exists(OutputFile):
        OutputFile=os.path.join(RootDir,'{}_FS{}.nii.gz'.format(FileBase,n))
        n=n+1
    
    WrkDir=os.path.join(RootDir,'{}_Work'.format(FileBase))
        
    TxtOutput=os.path.join(WrkDir,'WorkFiles')
    JobDir=os.path.join(WrkDir,'ShFiles')
    ErrorOutput=os.path.join(WrkDir,'StdError')
    StdOutput=os.path.join(WrkDir,'StdOutput')
    
    
    if not os.path.exists(WrkDir):
        os.makedirs(WrkDir)
    
    if not os.path.exists(TxtOutput):
        os.makedirs(TxtOutput)
    
    if not os.path.exists(JobDir):
        os.makedirs(JobDir)
        
    if not os.path.exists(ErrorOutput):
        os.makedirs(ErrorOutput)
    
    if not os.path.exists(StdOutput):
        os.makedirs(StdOutput)
        
    
    f.FileBase=FileBase
    f.FileName=FileName
    f.OutputFile=OutputFile
    f.TxtOutput=TxtOutput
    f.JobDir=JobDir
    f.ErrorOutput=ErrorOutput
    f.StdOutput=StdOutput
    f.WrkDir=WrkDir
    f.InputNii=InputNii
    f.RootDir=RootDir
    return(f)

def CompileNii(f):
    OriginalFile=f.InputNii
    
    UnCorNii=os.path.join(OriginalFile)
    UnCorNii=nb.load(UnCorNii)
    lx,ly,lz,lt=UnCorNii.get_data().shape
    
    BaseDir=f.RootDir
    fname=f.FileBase
    subID=fname

    SplitBase=subID
    
    StackList=np.zeros((lx,ly,1,lt))
    
    txtDir=f.TxtOutput
    
    for z in range(lz):
        zstack=np.zeros((lx*ly,lt))
        ztxts=os.path.join(txtDir,'FS_{}_{}_*.txt'.format('Image',z))
        nzchunks=len(glob.glob(ztxts))
        ind=np.loadtxt(os.path.join(txtDir,'SliceInd_{}.txt'.format(z)))
        ind=np.array([int(n) for n in ind])
        
        zchunks=[]
        
        
        for ii in range(nzchunks):
            #print(os.path.join(txtDir,'FS_{}_{}_{}.txt'.format('Image',z,ii)))
            newdata=np.loadtxt(os.path.join(txtDir,'FS_{}_{}_{}.txt'.format('Image',z,ii)))
            
            if newdata.ndim==1:
                newdata=np.expand_dims(newdata,0)
                
            zchunks.extend(newdata)
            
        if nzchunks>0:
            zstack[ind,:]=zchunks
            
        zstack=np.reshape(zstack,(lx,ly,1,lt))
        StackList=np.concatenate((StackList,zstack),2)
        
    StackList=StackList[:,:,1:,:]
    
    img=nb.Nifti1Image(StackList,UnCorNii.get_affine(),UnCorNii.get_header())
    savefile=os.path.join(BaseDir,'{}_FS'.format(subID))
    nb.loadsave.save(img,savefile)
    return()

def main(WorkingFile):
    
    if WorkingFile.Cores==0:
        WorkingFile.Cores=multiprocessing.cpu_count()
        if WorkingFile.Cores>2:
            WorkingFile.Cores=WorkingFile.Cores-1

    WorkingFile = SetUpWorkingPaths(WorkingFile)
    
    x,y,z,t=nb.load(WorkingFile.InputNii).get_shape()
    Fnew=np.ceil(z*WorkingFile.Fs)    
    WorkingFile.Fnew=Fnew
    
    if os.path.exists(WorkingFile.OutputFile):
        print('FS arleady run on this file:\nOutput Located at:\n{}'.format(WorkingFile.OutputFile))
        quit()
    
    WorkingFile.SubmitFiles=CreateJobsFromFile(WorkingFile)

    NtxtFiles=len(glob.glob(os.path.join(WorkingFile.TxtOutput,'Image*.txt')))
    ProcessedFiles=len(glob.glob(os.path.join(WorkingFile.TxtOutput,'FS_Image*.txt')))
    if not ProcessedFiles==NtxtFiles:
        p=Pool(WorkingFile.Cores)
        p.map(FiltShift,WorkingFile.SubmitFiles)
        p.close()
        p.join()
    else:
        print 'Completed for {}'.format(WorkingFile.FileBase)
    
    ProcessedFiles=len(glob.glob(os.path.join(WorkingFile.TxtOutput,'FS_{}*.txt'.format('Image'))))
    
    while ProcessedFiles<NtxtFiles:
        print('Waiting for Files to Finish')
        ProcessedFiles=len(glob.glob(os.path.join(WorkingFile.TxtOutput,'FS_{}*.txt'.format('Image'))))
        time.sleep(60)
        
    CompileNii(WorkingFile)
    
    cmd='gzip -f {}/{}_FS.nii'.format(WorkingFile.RootDir,WorkingFile.FileBase)
    pr=sp.Popen(cmd,shell=True)
    pr.wait()

if __name__ == '__main__':
    
    #-in <InputFile> -tr <TR> -i <Int> [-s <StartSlice>] [-d <Direction>] [-o <SliceORderFile>] [-r <RefSlice>] [-c <nCores>] [-m <Memory>] [-Force]'
    if len(sys.argv)>=3:

        try:
            f=WorkingFile()
            f.Fs=0
            f.Int=0
            f.Ref=0
            f.Cores=0
            f.mem=2
            f.force=False
            f.startSlice=0
            f.Direction=1
            f.OrderFile=''
            f.UseOrder=False
            f.TR=0
                        
            arg=[]
            arg=sys.argv[1:]
        
            while arg:
                i=str(arg.pop(0))
                i=i.rstrip()
                
                if i=='-tr':
                    f.TR=float(arg.pop(0))
                    f.Fs=1./f.TR
                    
                elif i=='-in':
                    f.InputNii=arg.pop(0)
            
                    if not os.path.exists(f.InputNii):
                        print ('Invalid File: {}'.format(f.InputNii))
                        quit()
                
                elif i=='-s':
                    f.startSlice=int(arg.pop(0))-1
                    
                elif i=='-d':
                    f.Direction=int(arg.pop(0))
                    if not (f.Direction==1 or f.Direction ==-1):
                        print ('Invalid value for Direction (-d).  Must be 1 or -1')
                        quit()
                
                elif i=='-o':
                    f.UseOrder=True
                    f.OrderFile=arg.pop(0)

                elif i=='-i':
                    f.Int=int(arg.pop(0))

                elif i=='-r':
                    f.Ref=int(arg.pop(0))-1

                elif i=='-c':
                    f.Cores=int(arg.pop(0))

                elif i=='-m':
                    f.mem=float(arg.pop(0))
                    
                elif i=='-Force':
                    f.force=True
                    
                elif i=='-out':
                    f.UserOutput=arg.pop(0)
                    
                else:
                    print('Unrecognized Argument "{}"'.format(i))

            main(f)

        except:
            
            print traceback.print_exc(file=sys.stdout)            
          
    else:
        usage()
                    
