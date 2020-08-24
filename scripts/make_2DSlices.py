#!/usr/bin/env python
#!/usr/bin/python2

# usage: plot_2DSlices.py [-h] [-m] [-b] [-r] [-c] fileName <eosFileName>

# Plots 2D slices from xml configuration file specified by "fileName". A
# description of the xml configuration file can be found in the reference file
# "plot_2DSlices_reference.xml" found with this script


#import getopt
import argparse
import sys
import os
from math import *
import glob
import numpy as np
import time
import disect_filename
import paths
import xml.etree.ElementTree as xml
import re

nNumCoords=7
colors=['r','g','b','c','m','y']
defaultXML="make_2DSlices_reference.xml"
class File2DSlice:
  def load(self,fileName):
    '''sets:
    fileName, file name of the 2D slice
    planeType, type of the 2D slice ("rt","rp", "tp")
    eosFile, file name of the equition of state file, if using a gamma-law gas it is None
    gamma, value of gamma for a gamma-law gass, if using an equation of state table it is None
    coordinateNames, Names of the coordinates
    coordinates, values of the coordinates
    dataNames, names of the data columns
    data, the data columns
    '''
    
    self.fileName=fileName
    if not os.access(fileName,os.F_OK|os.R_OK):
      print("error opening ",fileName," for reading")
      return
    f=open(fileName,'r')
    
    #set time step index
    temp=fileName[fileName.find("_t")+2:]
    self.index=temp[:temp.find("_2D")]
    
    #set type of plane
    self.planeType=f.readline().strip()
    
    #read time
    line=f.readline()
    words=line.split()
    self.time=float(words[1])
    
    #read gamma/eos
    if words[2]!="0":
      self.eosFile=words[3]
      self.gamma=None
    else:
      self.eosFile=None
      self.gamma=float(words[3])
    
    #read coordinates
    self.coordinateNames=[]
    self.coordinates=[]
    for i in range(nNumCoords):
      line=f.readline()
      words=line.split()
      self.coordinates.append([])
      self.coordinateNames.append(words[0])
      for j in range(1,len(words)):
        self.coordinates[i].append(float(words[j]))
    
    #read in data names
    line=f.readline()
    words=line.split()
    self.dataNames=words
    
    #read in data
    self.data=[]
    line=f.readline()
    nCount=0#indicates the variable
    for word in line.split():
      self.data.append([])
      if word!="-":
        self.data[nCount].append(float(word))
      nCount+=1
    for line in f:
      nCount=0
      for word in line.split():
        if word!="-":
          self.data[nCount].append(float(word))
        nCount+=1
def isFloat(testFloat):
  '''try to convert to a float, if an exception is thrown it isn't a float otherwise it is'''
  try:
    float(testFloat)
  except:
    return False
  else:
    return True
def isInt(testInt):
  '''try to convert to a integer, if an exception is thrown it isn't an integer otherwise it is'''
  try:
    int(testInt)
  except:
    return False
  else:
    return True
def main():
  #set parser options
  parser=argparse.ArgumentParser(description="Plots 2D slices from xml configuration file "\
    +"specified by \"fileName\". A description of the xml configuration file can be found in the"\
    +" reference file \"plot_2DSlices_reference.xml\" found with this script")
  parser.add_argument('fileName',action="store",type=str,help="Name of the xml configuration file")
  parser.add_argument('-p',dest='plot',action="store_true",default=False, help="Only plot existing slices, do not remake")
  parser.add_argument('-r',action="store_true",default=False,help="Remove distributed binary files")
  parser.add_argument('-c',action="store_true",default=False,help="Show codes as they will be "\
    +"executed, mostly a debugging option")
  parser.add_argument('eosFile', action="store", nargs='?', const=1, default="", type=str, help="Name of eosFile if not specified in binary file")
  
  #parse arguments
  parsed=parser.parse_args()
  XML=parsed.fileName; 
  #get xml settings
  if XML[0:1] == "." or XML[0:1] == "/" :  #If an absolute path to the XML is given, use it.
    print("Plotting 2D slices with XML: " + XML)
    settings=parseXMLFile(XML);
  else: #if an absolute path is not given, assume the XML file is in the scripts directory
    print("Plotting 2D Slices with XML: " + paths.scriptPath + XML)
    settings=parseXMLFile(paths.scriptPath + XML)
  #create slices
  createSlices(settings,parsed)
def parseXMLFile(fileName):
  #get root element
  tree=xml.parse(fileName)
  root=tree.getroot()
  settings={}
  
  if root==None:
    print("Require a \"figure\" root node")
  
  #get output file and extension
  if root.get("outputFile")!=None and root.get("outputFile")!="":
    outputFile=root.get("outputFile")
    [path,ext]=os.path.splitext(outputFile)
    supportedFileTypes=["png", "pdf", "ps", "eps", "svg"]
    if ext[1:] not in supportedFileTypes:
      print("File type \""+ext[1:]+"\" not suported. Supported types are ",supportedFileTypes,\
        " please choose one of those")
      quit()
    settings['outputFile']=path
    settings['fileFormat']=ext[1:]
  else:#ues defaults
    settings['outputFile']="2Dslice"
    settings['fileFormat']="png"
  
  #get figure width
  if root.get("width")!=None and root.get("width")!="":
    if isFloat(root.get("width")):
      settings['figWidth']=float(root.get("width"))
    else:
      print("Expecting a float for \"width\" attribute of node \"figure\"")
      quit()
  else:
    settings['figWidth']=20.0
    
  #get figure height
  if root.get("height")!=None and root.get("height")!="":
    if isFloat(root.get("height")):
      settings['figHeight']=float(root.get("height"))
    else:
      print("Expecting a float for \"height\" attribute of node \"figure\"")
      quit()
  else:
    settings['figHeight']=10.0
  
  #set location of left of the plot area
  settings['figLeft']=0.05
  if root.get("figLeft")!=None and root.get("figLeft")!="":
    settings['figLeft']=float(root.get("figLeft"))
  
  #set location of top of the plot area
  settings['figTop']=0.95
  if root.get("figTop")!=None and root.get("figTop")!="":
    settings['figTop']=float(root.get("figTop"))
  
  #set location of bottom of the plot area
  settings['figBottom']=0.05
  if root.get("figBottom")!=None and root.get("figBottom")!="":
    settings['figBottom']=float(root.get("figBottom"))
  
  #set font size
  settings['fontSize']=12.0
  if root.get("fontSize")!=None and root.get("fontSize")!="":
    settings['fontSize']=float(root.get("fontSize"))
  
  #set location of right side of the figure
  settings['figRight']=0.95
  if root.get("figRight")!=None and root.get("figRight")!="":
    settings['figRight']=float(root.get("figRight"))
  
  #figure horizontal spacing
  settings['figSpace']=0.0
  if root.get("figSpace")!=None and root.get("figSpace")!="":
    settings['figSpace']=float(root.get("figSpace"))
  
  #get figure dpi
  if root.get("dpi")!=None and root.get("dpi")!="":
    if isInt(root.get("dpi")):
      settings['figDpi']=int(root.get("dpi"))
    else:
      print("Expecting an integer for \"dpi\" attribute of node \"figure\"")
      quit()
  else:
    settings['figDpi']=90
  
  #get starting index for output file names
  if root.get("startIndex")!=None and root.get("startIndex")!="":
    if isInt(root.get("startIndex")):
      settings['startIndex']=int(root.get("startIndex"))
    else:
      print("Expecting an integer for \"startIndex\" attribute of node \"figure\"")
      quit()
  else:
    settings['startIndex']=0
  
  #get figure title
  if root.get("title")!=None:
    settings['title']=root.get("title")
  else:
    settings['title']=""
  
  #get baseFileName
  inputFileNameElement=root.findall("inputFileName")[0]
  if root.findtext("inputFileName")==None or root.findtext("inputFileName")=="":
    print("Requires an \"inputFileName\" node")
    quit()
  settings['inputFileName']=root.findtext("inputFileName")
  
  #get file frequency
  settings["fileFrequency"]=1
  if inputFileNameElement.get("frequency")!=None:
    settings["fileFrequency"]=int(inputFileNameElement.get("frequency"))
    
  
  # get plane elements
  planeElements=root.findall("plane")
  if planeElements==None:
    print("Requires at least one plane node")
    quit()
  planes=[]
  for planeElement in planeElements:
    
    plane={}
    #get plane type
    allowedPlaneTypes=["rt","rp","tp"]
    if planeElement.get("type") not in allowedPlaneTypes:
      print("Requires a plane \"type\" attribute from one of the following",allowedPlaneTypes)
      quit()
    plane['planeType']=planeElement.get("type")
    
    #get plane index
    if not isInt(planeElement.get("index")):
      print("Requires a plane \"index\" attribute that indicates the location of the plane in the "\
        +"model, should be an integer indicating the zone in the direction perpendicular to the"\
        +" plane.")
      quit()
    plane['planeIndex']=int(planeElement.get("index"))
    
    #get grid
    plane['grid']=None#use default
    if planeElement.get("grid")!=None:
      if planeElement.get("grid").lower() in ["both","major"]:
        plane['grid']=planeElement.get("grid").lower()
    
    #get x-axis
    xaxisElement=planeElement.find("xaxis")
    if xaxisElement==None:
      print("Requires an \"xaxis\" node")
      quit()
   
   #get xmin
    if xaxisElement.get("min")!=None:
      if not isFloat(xaxisElement.get("min")):
        print("Expecting a float for xaxis \"min\" attribute")
        quit()
      plane['xMin']=float(xaxisElement.get("min"))
    else:#use default
      plane['xMin']=None
    
    #get xmax
    if xaxisElement.get("max")!=None:
      if not isFloat(xaxisElement.get("max")):
        print("Expecting a float for xaxis \"max\" attribute")
        quit()
      plane['xMax']=float(xaxisElement.get("max"))
    else:#use default
      plane['xMax']=None
      
    #get xaxis formula
    if xaxisElement.text==None or xaxisElement.text=="":
      print("Must specify an xaxis varible")
      quit()
    plane['xFormula']=xaxisElement.text
    
    #get xlabel
    if xaxisElement.get("label")!=None:
      plane['xLabel']=xaxisElement.get("label")
    else:#use default
      plane['xLabel']="x axis label"
    
    #get x minorTics
    plane['xminortics']=False#use default
    if xaxisElement.get("minortics")!=None:
      if xaxisElement.get("minortics").lower() in ["true","1","t","yes","y"]:
        plane['xminortics']=True
    
    #get y-axis
    yaxisElement=planeElement.find("yaxis")
    if yaxisElement==None:
      print("Requires an \"yaxis\" node")
      quit()
   
   #get ymin
    if yaxisElement.get("min")!=None:
      if not isFloat(yaxisElement.get("min")):
        print("Expecting a float for yaxis \"min\" attribute")
        quit()
      plane['yMin']=float(yaxisElement.get("min"))
    else:#use default
      plane['yMin']=None
    
    #get xmax
    if yaxisElement.get("max")!=None:
      if not isFloat(yaxisElement.get("max")):
        print("Expecting a float for yaxis \"max\" attribute")
        quit()
      plane['yMax']=float(yaxisElement.get("max"))
    else:#use default
      plane['yMax']=None
      
    #get xaxis formula
    if yaxisElement.text==None or yaxisElement.text=="":
      print("Must specify an yaxis varible")
      quit()
    plane['yFormula']=yaxisElement.text
    
    #get ylabel
    if yaxisElement.get("label")!=None:
      plane['yLabel']=yaxisElement.get("label")
    else:#use default
      plane['yLabel']="y axis label"
    
    #get y minorTics
    plane['yminortics']=False#use default
    if yaxisElement.get("minortics")!=None:
      if yaxisElement.get("minortics").lower() in ["true","1","t","yes","y"]:
        plane['yminortics']=True
    
    #get scalor
    scalorElement=planeElement.find("scalor")
    if scalorElement==None:
      plane['scalorFormula']=None
      plane['scalorMin']=sys.float_info.min
      plane['scalorMax']=sys.float_info.max
      plane['scalorPallet']="jet"
      plane['scalorLabel']=None
    else:
      
      #get scalorMin
      if scalorElement.get("min")!=None:
        if not isFloat(scalorElement.get("min")):
          print("Expecting a float for scalor \"min\" attribute")
          quit()
        plane['scalorMin']=float(scalorElement.get("min"))
      else:#use default
        plane['scalorMin']=None
      
      #get scalorMax
      if scalorElement.get("max")!=None:
        if not isFloat(scalorElement.get("max")):
          print("Expecting a float for scalor \"max\" attribute")
          quit()
        plane['scalorMax']=float(scalorElement.get("max"))
      else:#use default
        plane['scalorMax']=None
      
      #get scalor formula
      if scalorElement.text==None or scalorElement.text=="":
        print("Must specify an scalor varible")
        quit()
      plane['scalorFormula']=scalorElement.text
      
      #get scalor label
      if scalorElement.get("label")!=None:
        plane['scalorLabel']=scalorElement.get("label")
      else:#use default
        plane['scalorLabel']="scalor label"
      
      #get scalor pallet
      if scalorElement.get("pallet")!=None:
        plane['scalorPallet']=scalorElement.get("pallet")
        if plane['scalorPallet']=="stellar":
          if scalorElement.get("palletFocus")!=None:
            if not isFloat(scalorElement.get("palletFocus")):
              print("Expecting a float for scalor \"palletFocus\" attribute")
              quit()
            plane['palletFocus']=float(scalorElement.get("palletFocus"))
          else:#use default
            print("Scalor pallet \"stellar\" also needs the attribute \"palletFocus\" to be set."\
              " This indicates which value of the scale the pallet should set to white.")
      else:#use default
        plane['scalorPallet']="jet"
    
    #read in vectors
    vectorElements=planeElement.findall("vector")
    vectors=[]
    for vectorElement in vectorElements:
      
      vector={}
      
      #get xfrequency
      vector['xfrequency']=1
      if isInt(vectorElement.get("xfrequency")):
        vector['xfrequency']=int(vectorElement.get("xfrequency"))
      
      #get yfrequency
      vector['yfrequency']=1
      if isInt(vectorElement.get("yfrequency")):
        vector['yfrequency']=int(vectorElement.get("yfrequency"))
      
      #get scale
      vector['scale']=1.0
      if isFloat(vectorElement.get("scale")):
        vector['scale']=float(vectorElement.get("scale"))
        
      #get vector thickness
      vector['thickness']=1.0
      if isFloat(vectorElement.get("thickness")):
        vector['thickness']=float(vectorElement.get("thickness"))
      
      #get color
      vector['color']='k'
      if vectorElement.get("color")!=None and vectorElement.get("color")!="":
        vector['color']=vectorElement.get("color")
      
      #get label
      labelElement=vectorElement.find("label")
      if labelElement!=None:
        
        #get label xposition
        if labelElement.get("xpos")!=None:
          if isFloat(labelElement.get("xpos")):
            vector['labelXPos']=labelElement.get("xpos")
          else:
            print("Expecting a float for \"xPos\" in vector label of vector",nCount)
        else:#default
          vector['labelXPos']=0.1
          
        #get label yposition
        if labelElement.get("ypos")!=None:
          if isFloat(labelElement.get("ypos")):
            vector['labelYPos']=labelElement.get("ypos")
          else:
            print("Expecting a float for \"ypos\" in vector label of vector",nCount)
        else:#default
          vector['labelYPos']=0.92
      
        #get text
        if labelElement.text!=None:
          vector['label']=labelElement.text
        else:#use default
          vector['label']=""
      else:
        vector['label']=""
        vector['labelYPos']=0.0
        vector['labelXPos']=0.0
      
      #get xposition formula
      if vectorElement.findtext("xposition")==None or vectorElement.findtext("xposition")=="":
        print("Must specify an xposition varible for vector")
        quit()
      vector['xposition']=vectorElement.findtext("xposition")
        
      #get yposition formula
      if vectorElement.findtext("yposition")==None or vectorElement.findtext("yposition")=="":
        print("Must specify an yposition varible for vector")
        quit()
      vector['yposition']=vectorElement.findtext("yposition")
      
      #get xcomponent formula
      if vectorElement.findtext("xcomponent")==None or vectorElement.findtext("xcomponent")=="":
        print("Must specify an xcomponent varible for vector")
        quit()
      vector['xcomponent']=vectorElement.findtext("xcomponent")
      
      #get ycomponent formula
      if vectorElement.findtext("ycomponent")==None or vectorElement.findtext("ycomponent")=="":
        print("Must specify an ycomponent varible for vector")
        quit()
      vector['ycomponent']=vectorElement.findtext("ycomponent")
      
      #add vector settings to list of vectors
      vectors.append(vector)
    
    #add list of vectors to plane
    plane['vectors']=vectors
    
    #add plane to list of planes
    planes.append(plane)
  
  #add list of planes to settings
  settings['planes']=planes
  return settings
def createSlices(settings,parsed):
  #get base file name
  [start,end,baseFileName]=disect_filename.disectFileName(settings['inputFileName'])
  
  nCount=0
  for plane in settings['planes']:
    
    #make sure that all combined binary files have 2D slices made
    if plane['planeType']=="rt":
      nPlaneID=0
      planeID="k"
    if plane['planeType']=="tp":
      nPlaneID=1
      planeID="i"
    if plane['planeType']=="rp":
      nPlaneID=2
      planeID="j"

    if(parsed.eosFile!=""):
      cmd = 'mk2DSlice' + ' ' + "\"" + settings['inputFileName'] + "\"" + ' ' + parsed.eosFile + ' ' + str(nPlaneID) + ' ' + str(plane['planeIndex']);
    else:
      cmd = 'mk2DSlice' + ' ' + "\"" + settings['inputFileName'] + "\"" + str(nPlaneID) + ' ' + str(plane['planeIndex']);
    print(cmd)
    os.system(cmd);

    #get and sort files
if __name__ == "__main__":
  main()
  