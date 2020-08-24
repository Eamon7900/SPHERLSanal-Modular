import csv
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import os
import sys
import math
import scipy.interpolate as interpolate
from scipy.signal import argrelextrema
from scipy import signal
from scipy import stats
from scipy import interpolate

#A class to store and process the data from one ModelDump
class ModelDump:
    fileName = "" #The fileName of this model dump 
    time = 0 #Time for this ModelDump
    phase = 0 #Phase of this ModelDump
   

    def __init__(self, fileName, time, grid):
        self.fileName = fileName
        self.time = time
        self.grid = grid
        
    def calcPhase(self, period):
        
        self.phase = (self.time % period)/period

#A class to store and process a set of ModelDumps
class DataSet:
    #Define instance variables:
    fileRanges = [] #Array of Dictionaries containing start, end and baseName for each filerange
    data = [] #Array of ModelDumps containing their respective data
    times = [] #Array of the times of all model dumps this DataSet contains 
    dataDir = "data/"
    profileExt = "_pro.txt"
    period = 0

    #Disects a filerange into dictionary containining a start, and end, and a basefilename

    @staticmethod
    def disectFileRange(fileRange):
        #get base file name
        parts0=fileRange.split('[',1)
        if len(parts0)<2:# maybe a single file?
            parts3=fileRange.rsplit('_t',1)
            if parts3[1].isdigit():
                start=int(parts3[1])
                end=start+1
                baseFileName=parts3[0]+'_t'
            else:
                print(disectFileRange.__name__+": unrecognized file type "+fileRange\
                    +" expecting something with an _tXXXXXXXX suffix, where X's denote digits.")
                return {} 

        else:# a range of files
            baseFileName=parts0[0]
            
            #get interation range
            parts1=parts0[1].split(']')
            parts2=parts1[0].split('-')
            start=int(parts2[0])
            if parts2[1]=="*":
                end=99999999
            else:
                end=int(parts2[1])        
        return {
            "start" : start,
            "end" : end,
            "baseName" : baseFileName
        } 

    def __init__(self, dataPath, fileRanges):
        if not dataPath == "":
            self.dataDir = dataPath + "/"
        print(self.dataDir) #TEST
        for fileRange in fileRanges:
            self.fileRanges.append(DataSet.disectFileRange(fileRange))
        for fileRange in self.fileRanges:
            for i in range(fileRange["start"], fileRange["end"]):
                #Fill zeroes for each file in range
                fileName = self.dataDir + fileRange["baseName"] + str(i).rjust(8,'0') + self.profileExt
                if(not os.path.isfile(fileName)): continue
                with open(fileName) as dataFile:
                    # Generate the time array:
                    time= dataFile.readline()
                    time = time.strip().split()
                    time = np.float(time[1])
                    self.times.append(time)

                    #Read in the data from the grid of each dump file
                    data = pd.read_csv(dataFile, delimiter='\s+').to_dict(orient='list')
                    #curGrid = [] #Each grid is an array of dictionaries (one for each row)
                    #for row in data:
                    #    curGrid.append(row)
                    self.data.append(ModelDump(fileName, time, data))

        #print("Read in time series of: " + str(self.times))
        self.times = np.array(self.times)
        
        # calculate the period using autocorrelation:
        seriesTime = 20 * 86400 # always calculate over last 20 days of the time series
        periodResolution = 0.0001 * 24 * 3600

        nsteps = np.int(seriesTime / periodResolution)#
        nsteps = np.min((10**5, nsteps))

        maxTime = max(self.times)
        minTime = maxTime - seriesTime

        #use the last U0 in each model dump with a stable time
        StableTime = []
        Ur0_short = []
        for modelDump in self.data:
            if modelDump.time > minTime:
                Ur0_short.append(modelDump.grid["U0"][-1]) 
                StableTime.append(modelDump.time)

        print("Stable times for period autocorrelation: " + str(StableTime))
        print("Stable Ur0 for period autocorrelation: " + str(Ur0_short))

        print("\n Beginning autocorrelation")
        minTime = np.max((StableTime[0],minTime))
        print("Minimum Time:" + str(minTime))

        InterpedTime = np.linspace(minTime, maxTime, nsteps)
        print("Interpolated time series: " + str(InterpedTime))

        Interped = interpolate.interp1d(StableTime, Ur0_short)
        InterpedUr0 = Interped(InterpedTime[1:])
        print("Interpolated Ur0: " + str(InterpedUr0))


        Correlation = np.correlate(InterpedUr0,InterpedUr0,mode='full')
        print("Correlation: " + str(Correlation))

        timestep = (np.max(StableTime) - np.min(StableTime)) / nsteps
        peaks = signal.argrelextrema(Correlation[nsteps:nsteps+np.int(nsteps/4)],np.greater)
        print("Peaks: " + str(peaks))

        period  = np.mean(np.diff(peaks))*timestep/86400
        PeriodError = np.std(np.diff(peaks))*timestep/86400

        print("calculated period of " + str(period) + " from data set") #DEBUG 

        self.period = period
        self.pererr = PeriodError

        
        