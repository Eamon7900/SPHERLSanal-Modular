#include <math.h>
#include <iomanip>
#include <unistd.h>
#include <string.h>
#include <vector>
#include <glob.h>

#include "binfile.h"
#include "datafile.h"

DataFile::DataFile(string fileName) 
: BinaryFile(fileName) {}

double DataFile::dCalRhoAve3D(double ****dGrid,int nI,int nStartY,int nEndY,int nStartZ,int nEndZ){
  int j;
  int k;
  double dSum=0.0;
  double dVolume=0.0;
  double dRFactor=0.33333333333333333*(pow(dGrid[this->nR][nI][0][0],3.0)
    -pow(dGrid[this->nR][nI-1][0][0],3.0));
  double dVolumeTemp;
  double dDeltaCosThetaIJK;
  double dDeltaPhi;
  for(j=nStartY;j<nEndY;j++){
    dDeltaCosThetaIJK=cos(dGrid[this->nTheta][0][j-1][0])-cos(dGrid[this->nTheta][0][j][0]);
    for(k=nStartZ;k<nEndZ;k++){
      dDeltaPhi=dGrid[this->nPhi][0][0][k]-dGrid[this->nPhi][0][0][k-1];
      dVolumeTemp=dRFactor*dDeltaCosThetaIJK*dDeltaPhi;
      dSum+=dVolumeTemp*dGrid[this->nD][nI][j][k];
      dVolume+=dVolumeTemp;
    }
  }
  return dSum/dVolume;
}

double DataFile::dCalRhoAve2D(double ****dGrid,int nI,int nStartY,int nEndY,int nStartZ,int nEndZ){
  int j;
  int k;
  double dSum=0.0;
  double dVolume=0.0;
  double dRFactor=0.33333333333333333*(pow(dGrid[this->nR][nI][0][0],3.0)
    -pow(dGrid[this->nR][nI-1][0][0],3.0));
  double dDeltaCosThetaIJK;
  double dVolumeTemp;
  for(j=nStartY;j<nEndY;j++){
    dDeltaCosThetaIJK=cos(dGrid[this->nTheta][0][j-1][0])-cos(dGrid[this->nTheta][0][j][0]);
    for(k=nStartZ;k<nEndZ;k++){
      dVolumeTemp=dRFactor*dDeltaCosThetaIJK;
      dSum+=dVolumeTemp*dGrid[this->nD][nI][j][k];
      dVolume+=dVolumeTemp;
    }
  }
  return dSum/dVolume;
}



