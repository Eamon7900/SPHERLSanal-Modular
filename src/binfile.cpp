#include <math.h>
#include <iomanip>
#include <unistd.h>
#include <string.h>

#include "binfile.h"

BinaryFile::BinaryFile(string fileName) {
    this->fileName = fileName;
}

double BinaryFile::dCalRhoAve3D(double ****dGrid,int nI,int nStartY,int nEndY,int nStartZ,int nEndZ){
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

double BinaryFile::dCalRhoAve2D(double ****dGrid,int nI,int nStartY,int nEndY,int nStartZ,int nEndZ){
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

void BinaryFile::setExeDir(){
  char buff[1024];
  ssize_t len = readlink("/proc/self/exe", buff, sizeof(buff)-1);
  if (len != -1) {
    buff[len] = '\0';
    this->sExeDir=std::string(buff);
    //find the first "/" from the end
    unsigned pos=this->sExeDir.find_last_of("/");
    
    //keep from the begging to the location of the last "/" to remove the name
    //of the executable
    this->sExeDir=this->sExeDir.substr(0,pos);
    
    //check to see if the last directory is "bin" if so remove that also
    //as installed versions put the exe's into the bin directory and sExeDir
    //should point the top level directory.
    pos=this->sExeDir.find_last_of("/");
    std::string sBin=this->sExeDir.substr(pos+1,3);
    
    //if installed remove bin directory
    if(sBin.compare("bin")==0){
      this->sExeDir=this->sExeDir.substr(0,pos);
    }
    
  } else {
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": error determining executable path"<<std::endl;
    throw exception2(ssTemp.str(),OUTPUT);
  }
}


/*static*/ bool BinaryFile::bFileExists(std::string strFilename){
  std::ifstream ifTest;
  ifTest.open(strFilename.c_str(),std::ios::in);
  if(!ifTest){
    return false;//doesn't exsist
  }
  else{
    ifTest.close();
    return true;//does exsist
  }
}