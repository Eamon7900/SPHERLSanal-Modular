import csv
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import os
import sys

#A class to store and process the data from one ModelDump
class ModelDump:
    fileName = "" #The fileName of this model dump 
    time = 0 #Time for this ModelDump
    grid = [] #Array of dictionaries containing data for each zone

    def __init__(self, fileName, time, grid):
        self.fileName = fileName
        self.time = time
        self.grid = grid

#A class to store and process a set of ModelDumps
class DataSet:
    #Define instance variables:
    fileRanges = [] #Array of Dictionaries containing start, end and baseName for each filerange
    data = [] #Array of ModelDumps containing their respective data
    times = [] #Array of the times of all model dumps this DataSet contains 
    dataDir = "data/"
    profileExt = "_pro.txt"
    
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
                print(disectFileName.__name__+": unrecognized file type "+fileName\
                    +" expecting something with an _tXXXXXXXX suffix, where X's denote digits.")
                return {} 

        else:# a range of files
            baseFileName=parts0[0]
            
            #get interation range
            parts1=parts0[1].split(']')
            parts2=parts1[0].split('-')
            start=int(parts2[0])
            if parts2[1]=="*":
                end=sys.maxsize
            else:
                end=int(parts2[1])        
        return {
            "start" : start,
            "end" : end,
            "baseName" : baseFileName
        } 

    def __init__(self,fileRanges):
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
                    print("TimeArray: " + str(self.times)) #TEST

                    #Read in the data from the grid of each dump file
                    data = pd.read_csv(dataFile, delimiter='\s+').to_dict(orient='records')
                    curGrid = [] #Each grid is an array of dictionaries (one for each row)
                    for row in data:
                        curGrid.append(row)
                    self.data.append(ModelDump(fileName, time, curGrid))

        self.times = np.array(self.times)
        # calculate the period using autocorrelation: