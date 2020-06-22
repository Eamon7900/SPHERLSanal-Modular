#ifndef MAIN_H
#define MAIN_H

/** @file
  
  Header file for \ref main.cpp.
  
*/

//#include "fftw++.h"
#include "../../config.h"
#ifdef FFTW_ENABLE
  #include <fftw3.h>
#endif
#ifdef HDF_ENABLE
  #include "mfhdf.h"
#endif
#include "exception2.h"

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <exception>
#include <sys/stat.h>
#include <cmath>
#include <csignal>
#include <fenv.h>
#include <limits>
#include <vector>

//functions
void convertDistBinToAscii(std::string sFileNameBase);
void combineBinFiles(std::string sFileNameBase);
void convertCollBinToAscii(std::string sFileName);
void convertCollAsciiToBin(std::string sFileName);
void makeRadialProFromColBin(std::string sFileName);
void printHelp();
bool bFileExists(std::string strFilename);
void fpSignalHandler(int nSig);
void make2DSlice(std::string sFileName,int nPlane,int nPlaneIndex);
void convertBinToLNA(std::string sFileName);
void computeFourierTrans(std::string sInFileName,std::string sOutFileName);
#endif
struct watchzone{
  std::vector<double> vecdT;//2
  std::vector<double> vecdU_ip1half;//3
  std::vector<double> vecdU_im1half;//4
  std::vector<double> vecdU0_ip1half;//5
  std::vector<double> vecdU0_im1half;//6
  std::vector<double> vecdQ0;//7
  std::vector<double> vecdV_jp1half;//8
  std::vector<double> vecdV_jm1half;//9
  std::vector<double> vecdQ1;//10
  std::vector<double> vecdW_kp1half;//11
  std::vector<double> vecdW_km1half;//12
  std::vector<double> vecdQ2;//13
  std::vector<double> vecdR_ip1half;//14
  std::vector<double> vecdR_im1half;//15
  std::vector<double> vecdDensity;//16
  std::vector<double> vecdDensityAve;//17
  std::vector<double> vecdE;//18
  std::vector<double> vecdP;//19
  std::vector<double> vecdTemp;//20
  std::vector<double> vecdDelM_r_t0;//21
  std::vector<double> vecdDelM_r;//22
  std::vector<double> vecdErrorDelM_r;//23
};

#ifdef HDF_ENABLE
void convertBinToHDF4(std::string sFileName);/**<
  converts a collected binary file to an hdf file
*/
#endif
