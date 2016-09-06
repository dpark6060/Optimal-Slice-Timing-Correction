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
        self.UseTiming=False
        self.TimingFile=''
        



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
    
    print('-o:\t Slice Order File.  This file is the order in which each slice was acquired.  If present, all interelave parameters are ignored, and slices are shifted using the slice order file')
    print('\t we refer to the bottom slice in the image as slice 1, not slice 0')

    print('-t:\t Slice Timing File.  If present, all interelave parameters are ignored, and slices are shifted using the slice order file')
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
    
    # Sets up jobs for the image based on slices and memory size
    
    SubmitFiles=[]
    
    # Initialize the TR and Reference slice    
    Tr=1.0/f.Fs
    RefSlice=0
    
    # Load Nifti Data
    Nii=nb.load(f.InputNii)
    NiiData=Nii.get_data()
    
    # Calculate the third percentile of data - anything below this will be condisered noise/ outside the brain
    # and will not be processed to save time.  In the final image, any voxels that fall into this category will
    # be set to zero.
    
    per=np.percentile(NiiData,3)
    print per
    lx,ly,lz,lt=NiiData.shape
    
    # If we are using a slice order file
    if f.UseOrder:
        try:
            SliceOrder=np.loadtxt(f.OrderFile)
            SliceOrder=SliceOrder-1
        except:
            print('Error Loading Slice Order Flie {}'.format(f.OrderFile))
        
        if not len(SliceOrder) ==lz:
            print('Slice Order File does not match number of slices in {}'.format(f.InputNii))
            quit()
        
        # Load the file and calculate the exact time of acquision based of the TR and the number of slices
        shift=(TR-TR*1.0/lz)/lz
        SliceShift=SliceOrder*shift
        SliceShift=SliceShift-SliceShift[Start]
    
    
    #TODO Add Slice Timing File Here
    
    # Otherwise just create a slice timing file using the start slice and intereleave parameter
    else:
        
        SliceOrder, SliceShift = MakeT(f.Int, f.Ref, lz, Tr,f.startSlice, f.Direction)
    
    print('Using Slice Order:\n{}'.format(SliceOrder+1))
    
    
    # Calculate the size that each job can be given the availible memory and number of cores.
    GBperCore=f.mem/f.Cores*0.9 # Availible memory (GB) divided by the number of cores - reduced by 0.9 - a conservitive measure on how much memory each job can have
    BytePerCore=GBperCore*np.power(1024,3)  # Calculate the number of bytes this represents per job
    Nelements=lt*2*(f.Fnew/f.Fs)    # Calculate the number of elements there will be in one voxel's time series after upsampling and padding
    MemTotal=Nelements*31.48        # 31.48 is a calculated value for how many bytes a single array element takes up in numpy - an estimate of how much memory a single voxel's timeseries will take up

    
    
    # If the amount of memory required for ONE voxel is more than we can spare given the number of cores we're using, reduce the number of cores by one.
    while MemTotal>BytePerCore: 
        f.Cores=f.Cores-1
        if f.Cores<1:
            print ('Need More memory, there WILL (probably) be a memory error.  use -Force to force')
            print ('THIS IS NOT RECOMMENDED, THIS WILL CRASH YOUR COMPUTER')
            
            # If for some reason you know better than my little program, force the computer to try to process it anyways.
            if f.force==False:
                quit()
         
        # Recalculate the bytes per core with the new core count   
        GBperCore=f.mem/f.Cores*.9
        BytePerElement=np.dtype(float).itemsize
        BytePerCore=GBperCore*np.power(1024,3)
        
    
    # Now determine just HOW many voxels we can send to each core at a time without going over our alloted memory
    # The bytes per core divided by the total memory required for one voxel, rounded down
    txtLen=int(np.floor(BytePerCore)/(MemTotal))
    
    # Saved in a parameters file that is referenced later if you try to rerun the analysis - if nothing in the parameters file has changed, the analysis does not re-run.
    ParameterFile=os.path.join(f.TxtOutput,'Parameters.txt')
    Rerun=True
    
    
    # If the parameters file exists already, load it and read some values
    if os.path.exists(ParameterFile):
        parameters=open(ParameterFile,'r')
        line=parameters.read()
        groups=re.search('txtLen:(.*)\n',line)
        OldtxtLen=int(groups.groups()[0])
        groups=re.search('Nelements:(.*)\n',line)
        OldNelements=int(np.round(float(groups.groups()[0])))
        
        # If we're still using the same number of voxels and elemets pre voxel, nothing has changed...
        if (OldtxtLen==txtLen and OldNelements==Nelements) and f.force==False:
            #...So do not rerurn.
            Rerun=False
 
    SliceLen=glob.glob(os.path.join(f.TxtOutput,'SliceLen_*.txt'))
    
    # If there are no previous sliceLen files (glob returns empty), or we're forcing a rerun
    if (not SliceLen) or Rerun:
        
        # Remove all previous logs if we're rerunning
        if Rerun:
            cmd='rm {}/*.txt'.format(f.TxtOutput)
            pr=sp.Popen(cmd,shell=True)
            pr.wait()
        
        #Extract NII data slice by slice and reshape into a 2-D array
        for z in range(lz):
            Extract=np.squeeze(NiiData[:,:,z,:])
            Extract=np.reshape(Extract,(-1,lt))
            LexOrig=len(Extract)
            
            # Any time series where the min along the time axis is still above the bottom 3% intensity value of the data, consider this signal
            ind=np.where(Extract.min(1)>per)
            # Extract those signals
            Extract=Extract[ind]
            
            # Create two text files, one that stores the length of the slice (number of voxels) before removing noise voxels, and another that stores the xyz coordinates of these voxels.
            lentxt=os.path.join(f.TxtOutput,'SliceLen_{}.txt'.format(z))
            indtxt=os.path.join(f.TxtOutput,'SliceInd_{}.txt'.format(z))
            
            # Save the values
            np.savetxt(lentxt,[LexOrig])
            np.savetxt(indtxt,ind[0])
            
            # Extract voxel time series' and save as text files for parallelization - remember, can only store txtLen voxels' worth of time series.
            for ii,inc in enumerate(range(0,len(Extract),txtLen)):
                
                # the current index plus the increment is larger than the file, just take to the end.
                if inc+txtLen>=len(Extract):
                    txtdat=Extract[inc:]
                else:
                    txtdat=Extract[inc:inc+txtLen]
                
                # Save the voxel time series as a text file indicating the slice number and the portion of that slice
                TxtName=os.path.join(f.TxtOutput,'Image_{}_{}.txt'.format(z,ii))                
                np.savetxt(TxtName,txtdat)
                
                # Create a job file containing the timing shift of the slice that this job is working on, the image file it's working on, and the output directory.
                SubmitFiles.append([TxtName,'{}'.format(SliceShift[z]),'{}'.format(f.Fs),'{}'.format(f.Fnew),os.path.join(f.StdOutput,'Image_{}_{}.txt'.format(z,ii)),os.path.join(f.ErrorOutput,'Image_{}_{}.txt'.format(z,ii))])
        
        # Write more values to the parameter file
        ParamOut=open(ParameterFile,'w')
        ParamOut.write('txtLen:{}\n'.format(txtLen))
        ParamOut.write('Nelements:{}\n'.format(Nelements))
        ParamOut.close()
    
    # But if the text files exist already, simply recreate the jobs.      
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
    
    # Save the slice order and slice shift files for reference.
    np.savetxt(os.path.join(f.TxtOutput,'SliceOrder.txt'),SliceOrder)
    np.savetxt(os.path.join(f.TxtOutput,'SliceShift.txt'), SliceShift)
    
    # Return a list of jobs to submit.
    return(SubmitFiles)

def SetUpWorkingPaths(f):
    # This creates a working directory structured for the code to work on.  
    InputNii=f.InputNii
    
    # If the user didn't specify an output directory, simply work in the path of the input NII file
    if f.UserOutput=='':    
        RootDir,FileName=os.path.split(InputNii)
    
    # Otherwise create and set up a working directory where they user specified
    else:
        RootDir,FileName=os.path.split(f.UserOutput)
        if RootDir=='':
            RootDir,trash=os.path.split(InputNii)

    # Strip the file name of its extension
    trash='.exe'
    FileBase=FileName
    while trash:
        FileBase,trash=os.path.splitext(FileBase)
    
    # Save the final output in the specified output directory with the suffix "_FS"
    OutputFile=os.path.join(RootDir,'{}_FS.nii.gz'.format(FileBase))
    
    # If that output exists, I guess you're rerunning, so let's just add a number to the end so it doesn't save over
    n=1
    while os.path.exists(OutputFile):
        OutputFile=os.path.join(RootDir,'{}_FS{}.nii.gz'.format(FileBase,n))
        n=n+1
    
    # Create the work directory and sub-directory structure
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
        
    # Set internal variables
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
    # Reconstruct the processed Time series text files into a NIFTI file
    
    # First load the original file for dimension information
    ## TODO This could be stored in the f structure maybe
    OriginalFile=f.InputNii    
    UnCorNii=os.path.join(OriginalFile)
    UnCorNii=nb.load(UnCorNii)
    lx,ly,lz,lt=UnCorNii.get_data().shape
    
    # Load directory information for output file structure
    BaseDir=f.RootDir
    fname=f.FileBase
    subID=fname    
    SplitBase=subID
    txtDir=f.TxtOutput
    
    # Reconstruction will occur by slice, each slice being this size and shape.
    StackList=np.zeros((lx,ly,1,lt))    
    
    # Loop through and reconstruct each slice
    for z in range(lz):
        
        #The Original time series matrix was lx*ly by lt in size, so create a variable to store that
        zstack=np.zeros((lx*ly,lt))
        
        # Find all the text files that belong to this slice
        ztxts=os.path.join(txtDir,'FS_{}_{}_*.txt'.format('Image',z))
        nzchunks=len(glob.glob(ztxts))
        
        # Load the SliceInd files, which tell us which index each row in the FS text file corresponds to in the Original time series file
        ind=np.loadtxt(os.path.join(txtDir,'SliceInd_{}.txt'.format(z)))
        ind=np.array([int(n) for n in ind])
        
        # The voxels we filtered aren't every voxel in the slice, so we'll store them in a separate variable that will be less than or equal to lx*ly by lt in size
        zchunks=[]
        
        # Now loop through the text files associated with the slice
        for ii in range(nzchunks):
            #print(os.path.join(txtDir,'FS_{}_{}_{}.txt'.format('Image',z,ii)))
            newdata=np.loadtxt(os.path.join(txtDir,'FS_{}_{}_{}.txt'.format('Image',z,ii)))
            
            # Resize the loaded data so it can be appended to a 2d matrix
            if newdata.ndim==1:
                newdata=np.expand_dims(newdata,0)
            
            # Add the new data do the zchunks 
            zchunks.extend(newdata)
        
        # Assuming that there were ANY timeseries text files for this slice, just assign the "zchunks" rows to the "zstack" rows as indicated by "ind"
        if nzchunks>0:
            zstack[ind,:]=zchunks
        
        # Reshape the processed time series back to the proper slice dimensions
        zstack=np.reshape(zstack,(lx,ly,1,lt))
        # Add this slice to the stack of slices
        StackList=np.concatenate((StackList,zstack),2)
    
    # Apparently my code attaches an extra slice of zeros, so just crop it off.
    StackList=StackList[:,:,1:,:]
    
    # Save the image using the affine and header from the original Nifti file
    img=nb.Nifti1Image(StackList,UnCorNii.get_affine(),UnCorNii.get_header())
    savefile=os.path.join(BaseDir,'{}_FS'.format(subID))
    nb.loadsave.save(img,savefile)
    return()

def main(WorkingFile):
    
    # If the user didn't set the cores (not recommended) we have to guess - take the core count and as long as you have more than two, take one fewer than the max number of cores.
    if WorkingFile.Cores==0:
        WorkingFile.Cores=multiprocessing.cpu_count()
        if WorkingFile.Cores>2:
            WorkingFile.Cores=WorkingFile.Cores-1
    
    # Set up a working directory based on user specifications (or lack there of)
    WorkingFile = SetUpWorkingPaths(WorkingFile)
    
    # Get the dimensions of the Nifti Image and Calculate the new frequency needed to sample the image
    x,y,z,t=nb.load(WorkingFile.InputNii).get_shape()
    Fnew=np.ceil(z*WorkingFile.Fs)    
    WorkingFile.Fnew=Fnew
    
    # If there's already an output file, we're quitting the program.  It takes too long to run.
    if os.path.exists(WorkingFile.OutputFile):
        print('FS arleady run on this file:\nOutput Located at:\n{}'.format(WorkingFile.OutputFile))
        quit()
    
    # Create jobs to submit to parallel processing
    WorkingFile.SubmitFiles=CreateJobsFromFile(WorkingFile)
    
    # This is so that if you error out mid-image, you won't need to rerun the ENTIRE thing, it's smart enough to see that there are slices that are already processed
    # Get a count of the number of text files that need to be be processed
    NtxtFiles=len(glob.glob(os.path.join(WorkingFile.TxtOutput,'Image*.txt')))
    # And a count of the number of files that ARE processed
    ProcessedFiles=len(glob.glob(os.path.join(WorkingFile.TxtOutput,'FS_Image*.txt')))
    # If they're equal, we're done...if not, create a parallel pool of the jobs created
    if not ProcessedFiles==NtxtFiles:
        p=Pool(WorkingFile.Cores)
        p.map(FiltShift,WorkingFile.SubmitFiles)
        p.close()
        p.join()
    else:
        print 'Completed for {}'.format(WorkingFile.FileBase)
    
    # Check for Completion condition - while the processed files are fewer than the original files, keep waiting for the jobs to run
    ProcessedFiles=len(glob.glob(os.path.join(WorkingFile.TxtOutput,'FS_{}*.txt'.format('Image'))))
    while ProcessedFiles<NtxtFiles:
        print('Waiting for Files to Finish')
        ProcessedFiles=len(glob.glob(os.path.join(WorkingFile.TxtOutput,'FS_{}*.txt'.format('Image'))))
        time.sleep(60)
    
    # Once the jobs are completed, Reconstruct the text files into a Nifti Image
    CompileNii(WorkingFile)
    
    # Gzip that if it isn't already
    cmd='gzip -f {}/{}_FS.nii'.format(WorkingFile.RootDir,WorkingFile.FileBase)
    pr=sp.Popen(cmd,shell=True)
    pr.wait()
    
    # That's it, we're done.  Go Home.

if __name__ == '__main__':
    
    # I'm not commenting this, it's just a crude parser.
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
                    
