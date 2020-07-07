#include <math.h>
#include <iomanip>
#include <unistd.h>
#include <glob.h>
#include <string.h>
#include <regex>
#include <vector>

#include "eos.h"
#include "exception2.h"
#include "binfile.h"
#include "paths.h"
#include "datafile.h"
#include "argparser.h"

void makeRadialProFromColBin(DataFile& bin);

//Creates a radial profile from a ColBin file.
int main(int argc, char* argv[]){
  string eosFile="";
  bool extraInfo=false;
  /*Parse arguments using switch since there are only 4 cases:
  mkRadPro -e <fileRange> <eosFile>: argc = 4
  mkRadPro -e <fileRange> : argc = 3
  mkRadPro <fileRange> <eosFile>: argc = 3
  mkRadPro <fileRange> : argc = 2;
  */

  string fileRange;
  switch(argc){
    case 4:
      fileRange=argv[2];
      eosFile=argv[3];
      break;
    case 3:
      if(strcmp(argv[1],"-e") == 0){ //no eosFile
        extraInfo=true;
        fileRange=argv[2];
      } else {
        fileRange=argv[1]; 
        eosFile=argv[2];
      }
      break;
   case 2:
     fileRange = argv[1];
     break;
  }

  ArgParser argParser(fileRange);
  std::vector<std::string> filesInRange = argParser.getFilesInRange(); 

  if(eosFile!="")
    std::cout << "creating radial profiles using eos file: " << eosFile << std::endl;
  //extract the numerical range of binary files from fileRange
  std::cout << filesInRange.empty() << std::endl;
  //Run through every file in range to check if it exists
  while(!filesInRange.empty()){
      DataFile curBin(filesInRange.back());
      filesInRange.pop_back();

      curBin.sEOSFile = eosFile;
      curBin.bExtraInfoInProfile = extraInfo;
      std::cout << "Creating rdial profile for: " << curBin.sFileName << std::endl;
      makeRadialProFromColBin(curBin);
  }

  return 0;
}

void makeRadialProFromColBin(DataFile& bin){//updated
  //open input file
  std::string sExtension=bin.sFileName.substr(bin.sFileName.size()-4,1);
  if(sExtension.compare(".")==0){//if there is an extension remove it
    bin.sFileName=bin.sFileName.substr(0,bin.sFileName.size()-4);
  }
  if(bin.sFileName.size()==0){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": no input file specified\n";
    throw exception2(ssTemp.str(),INPUT);
  }

  std::ifstream ifFile; 
  ifFile.open(bin.sFileName.c_str(),std::ios::binary);
  if(!ifFile.good()){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": input file \""
      <<bin.sFileName<<"\" didn't open properly\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //check that it is a binary file
  char cTemp;
  ifFile.read((char*)(&cTemp),sizeof(char));
  if(cTemp!='b'){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": input file \""
      <<bin.sFileName<<"\" isn't a binary file.\n";
    throw exception2(ssTemp.str(),INPUT); 
  }
  
  //check that it is the correct version
  int nTemp;
  ifFile.read((char*)(&nTemp),sizeof(int));
  if(nTemp!=bin.nDumpFileVersion){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": inpput file \""
      <<bin.sFileName<<"\" version \""<<nTemp
      <<"\" isn't the supported version \""<<bin.nDumpFileVersion<<"\".\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //read in time
  double dTime;;
  ifFile.read((char*)(&dTime),sizeof(double));
  
  //read in the time step index
  int nTimeStepIndex;
  ifFile.read((char*)(&nTimeStepIndex),sizeof(int));
  
  //read in the time step
  double dTimeStep1;
  ifFile.read((char*)(&dTimeStep1),sizeof(double));
  
  //read in the time step
  double dTimeStep2;
  ifFile.read((char*)(&dTimeStep2),sizeof(double));
  
  //read in alpha
  double dAlpha;
  ifFile.read((char*)(&dAlpha),sizeof(double));
  
  //read in gamma law
  int nGammaLaw;
  ifFile.read((char*)(&nGammaLaw),sizeof(int));
  
  double dGamma;
  std::string sEOSTable;
  eos* eosTable;
  if(nGammaLaw==0){
    ifFile.read((char*)(&dGamma),sizeof(double));
    eosTable = new eos();
  }
  else{
    char *cBuffer=new char[nGammaLaw+1];
    ifFile.read(cBuffer,nGammaLaw*sizeof(char));
    cBuffer[nGammaLaw]='\0';
    sEOSTable=cBuffer;
    delete [] cBuffer;
    
    if(bin.sEOSFile!=""){//overwrite sEOSTable if sEOSFile is set
      sEOSTable=bin.sEOSFile;
    }
    
    //set the exe directory
    //test to see if sEOSTable is relative to the execuatable directory
    std::string eosFileName;
    if (sEOSTable.substr(0,1)!="/" && sEOSTable.substr(0,2)!="./" && sEOSTable.substr(0,1)!="~"){
      //if absolute path not specified, assume EOS file is in eos folder. 
      eosFileName=bin.sRootDir+PATHS::EOS+sEOSTable;
    }
    else{
      eosFileName=sEOSTable;
    }

    eosTable = new eos(eosFileName); 
    eosTable->readBin(); 
  }
  
  //read in artificial viscosity
  double dA;
  ifFile.read((char*)(&dA),sizeof(double));
  
  //read in artificial viscosity threshold
  double dAVThreshold;
  ifFile.read((char*)(&dAVThreshold),sizeof(double));
  
  //read in global grid size
  int nSizeGlobe[3];
  ifFile.read((char*)(nSizeGlobe),3*sizeof(int));
  
  //read in periodicity
  int nPeriodic[3];
  ifFile.read((char*)(nPeriodic),3*sizeof(int));
  
  //read in number of 1D Zones
  int nNum1DZones;
  ifFile.read((char*)(&nNum1DZones),sizeof(int));
  
  //read in number of ghostcells
  int nNumGhostCells;
  ifFile.read((char*)(&nNumGhostCells),sizeof(int));
  
  //read in number of grid variables
  int nNumVars;
  ifFile.read((char*)(&nNumVars),sizeof(int));
  
  //get variable info, and set grid sizes
  int **nSize=new int*[nNumVars];
  int **nVarInfo=new int*[nNumVars];
  int l;
  for(int n=0;n<nNumVars;n++){
    nSize[n]=new int[3];
    nVarInfo[n]=new int[3];
    ifFile.read((char*)(nVarInfo[n]),(4)*sizeof(int));
    for(l=0;l<3;l++){
      if(nSizeGlobe[l]==1){
        nVarInfo[n][l]=-1;
      }
      if(nVarInfo[n][l]==-1){//variable not defined in direction l
        nSize[n][l]=1;
      }
      else if(nVarInfo[n][l]==1&&l==0){//interface variable
        nSize[n][l]=nSizeGlobe[l]+1;
      }
      else{
        nSize[n][l]=nSizeGlobe[l];
      }
    }
  }
  
  //figure out number of dimensions
  int nNumDims=0;
  if(nSizeGlobe[0]>1){
    nNumDims++;
  }
  if(nSizeGlobe[1]>1){
    nNumDims++;
  }
  if(nSizeGlobe[2]>1){
    nNumDims++;
  }
  
  //set variable indices
  int nNumIntVars=0;
  if(nGammaLaw==0){//using gamma law gas
    if(nNumDims==1){
      nNumIntVars=5;
      bin.nM=0;
      bin.nDM=1;
      bin.nR=2;
      bin.nD=3;
      bin.nU=4;
      bin.nU0=5;
      bin.nE=6;
      bin.nP=nNumVars+0;
      bin.nQ=nNumVars+1;
      bin.nKEP=nNumVars+2;
      bin.nC=nNumVars+3;
      bin.nKETot=nNumVars+4;
      bin.nV=-1;
      bin.nW=-1;
      bin.nT=-1;
      bin.nTheta=-1;
      bin.nPhi=-1;
      bin.nKappa=-1;
      bin.nGamma=-1;
    }
    else if(nNumDims==2){
      nNumIntVars=5;
      bin.nM=0;
      bin.nTheta=1;
      bin.nDM=2;
      bin.nR=3;
      bin.nD=4;
      bin.nU=5;
      bin.nU0=6;
      bin.nV=7;
      bin.nE=8;
      bin.nP=nNumVars+0;
      bin.nQ=nNumVars+1;
      bin.nKEP=nNumVars+2;
      bin.nC=nNumVars+3;
      bin.nKETot=nNumVars+4;
      bin.nW=-1;
      bin.nT=-1;
      bin.nPhi=-1;
      bin.nKappa=-1;
      bin.nGamma=-1;
    }
    else if(nNumDims==3){
      nNumIntVars=5;
      bin.nM=0;
      bin.nTheta=1;
      bin.nPhi=2;
      bin.nDM=3;
      bin.nR=4;
      bin.nD=5;
      bin.nU=6;
      bin.nU0=7;
      bin.nV=8;
      bin.nW=9;
      bin.nE=10;
      bin.nP=nNumVars+0;
      bin.nQ=nNumVars+1;
      bin.nKEP=nNumVars+2;
      bin.nC=nNumVars+3;
      bin.nKETot=nNumVars+4;
      bin.nT=-1;
      bin.nKappa=-1;
      bin.nGamma=-1;
    }
  }
  else{//using a tabulated equation of state
    if(nNumDims==1){
      nNumIntVars=11;
      bin.nM=0;
      bin.nDM=1;
      bin.nR=2;
      bin.nD=3;
      bin.nU=4;
      bin.nU0=5;
      bin.nT=6;
      bin.nE=nNumVars+0;
      bin.nQ=nNumVars+1;
      bin.nP=nNumVars+2;
      bin.nKappa=nNumVars+3;
      bin.nGamma=nNumVars+4;
      bin.nL_rad=nNumVars+5;
      bin.nL_con=nNumVars+6;
      bin.nF_con=nNumVars+7;
      bin.nKEP=nNumVars+8;
      bin.nC=nNumVars+9;
      bin.nKETot=nNumVars+10;
      bin.nV=-1;
      bin.nW=-1;
      bin.nTheta=-1;
      bin.nPhi=-1;
    }
    else if(nNumDims==2){
      nNumIntVars=11;
      bin.nM=0;
      bin.nTheta=1;
      bin.nDM=2;
      bin.nR=3;
      bin.nD=4;
      bin.nU=5;
      bin.nU0=6;
      bin.nV=7;
      bin.nT=8;
      bin.nE=nNumVars+0;
      bin.nQ=nNumVars+1;
      bin.nP=nNumVars+2;
      bin.nKappa=nNumVars+3;
      bin.nGamma=nNumVars+4;
      bin.nL_rad=nNumVars+5;
      bin.nL_con=nNumVars+6;
      bin.nF_con=nNumVars+7;
      bin.nKEP=nNumVars+8;
      bin.nC=nNumVars+9;
      bin.nKETot=nNumVars+10;
      bin.nPhi=-1;
      bin.nW=-1;
    }
    else if(nNumDims==3){
      nNumIntVars=11;
      bin.nM=0;
      bin.nTheta=1;
      bin.nPhi=2;
      bin.nDM=3;
      bin.nR=4;
      bin.nD=5;
      bin.nU=6;
      bin.nU0=7;
      bin.nV=8;
      bin.nW=9;
      bin.nT=10;
      bin.nE=nNumVars+0;
      bin.nQ=nNumVars+1;
      bin.nP=nNumVars+2;
      bin.nKappa=nNumVars+3;
      bin.nGamma=nNumVars+4;
      bin.nL_rad=nNumVars+5;
      bin.nL_con=nNumVars+6;
      bin.nF_con=nNumVars+7;
      bin.nKEP=nNumVars+8;
      bin.nC=nNumVars+9;
      bin.nKETot=nNumVars+10;
    }
  }
  
  //read in grid
  double ****dGrid=new double***[nNumVars];
  int nGhostCellsX;
  int nGhostCellsY;
  int nGhostCellsZ;
  int nSizeX1;
  int nSizeX2;
  int nStartY;
  int nEndY;
  int nSizeY;
  int nStartZ;
  int nEndZ;
  int nSizeZ;
  int i;
  int j;
  int k;
  double *dTempArray;
  double dSum;
  int nCount;
  for(int n=0;n<nNumVars;n++){
    
    nGhostCellsX=1;
    if(nVarInfo[n][0]==-1){
      nGhostCellsX=0;
    }
    nGhostCellsY=1;
    if(nVarInfo[n][1]==-1){
      nGhostCellsY=0;
    }
    nGhostCellsZ=1;
    if(nVarInfo[n][2]==-1){
      nGhostCellsZ=0;
    }
    
    //make some space to hold the variables
    dGrid[n]=new double**[nSize[n][0]+nGhostCellsX*2*nNumGhostCells];
    
    //read in 1D part of the grid
    nSizeX1=nGhostCellsX*(nNum1DZones+nNumGhostCells);//may be need to +1 if only one proc and variable in interface centered
    if (nVarInfo[n][0]==1&&nPeriodic[0]==0){
      nSizeX1=nGhostCellsX*(nNum1DZones+1+nNumGhostCells);
    }
    nSizeY=1;
    nSizeZ=1;
    for(i=0;i<nSizeX1;i++){
      dGrid[n][i]=new double*[nSizeY];
      for(j=0;j<nSizeY;j++){
        dGrid[n][i][j]=new double[nSizeZ];
        ifFile.read((char*)(dGrid[n][i][j]),nSizeZ*sizeof(double));
      }
    }
    
    //read in the rest of the grid
    nSizeX2=nSize[n][0]+nGhostCellsX*2*nNumGhostCells;
    nSizeY=nSize[n][1]+nGhostCellsY*2*nNumGhostCells;//assume y and z are always periodic
    nSizeZ=nSize[n][2]+nGhostCellsZ*2*nNumGhostCells;
    for(i=nSizeX1;i<nSizeX2;i++){
      dGrid[n][i]=new double*[nSizeY];
      for(j=0;j<nSizeY;j++){
        dGrid[n][i][j]=new double[nSizeZ];
        ifFile.read((char*)(dGrid[n][i][j]),nSizeZ*sizeof(double));
      }
    }
  }
  ifFile.close();
  
  //radialize the grid
  double **dMax=new double*[nNumVars+nNumIntVars];
  double **dMin=new double*[nNumVars+nNumIntVars];
  double **dAve=new double*[nNumVars+nNumIntVars];
  int **nMaxJIndex=new int*[nNumVars+nNumIntVars];
  int **nMaxKIndex=new int*[nNumVars+nNumIntVars];
  int **nMinJIndex=new int*[nNumVars+nNumIntVars];
  int **nMinKIndex=new int*[nNumVars+nNumIntVars];
  double dMaxTemp;
  double dMinTemp;
  for(int n=0;n<nNumVars;n++){
    
    nGhostCellsX=1;
    if(nVarInfo[n][0]==-1){
      nGhostCellsX=0;
    }
    nGhostCellsY=1;
    if(nVarInfo[n][1]==-1){
      nGhostCellsY=0;
    }
    nGhostCellsZ=1;
    if(nVarInfo[n][2]==-1){
      nGhostCellsZ=0;
    }
    
    //make some space to hold max,min and average
    dMax[n]=new double[nSize[n][0]+nGhostCellsX*2*nNumGhostCells];
    dMin[n]=new double[nSize[n][0]+nGhostCellsX*2*nNumGhostCells];
    dAve[n]=new double[nSize[n][0]+nGhostCellsX*2*nNumGhostCells];
    
    nMaxJIndex[n]=new int[nSize[n][0]+nGhostCellsX*2*nNumGhostCells];
    nMaxKIndex[n]=new int[nSize[n][0]+nGhostCellsX*2*nNumGhostCells];
    nMinJIndex[n]=new int[nSize[n][0]+nGhostCellsX*2*nNumGhostCells];
    nMinKIndex[n]=new int[nSize[n][0]+nGhostCellsX*2*nNumGhostCells];
    
    //read in 1D part of the grid
    nSizeX1=nGhostCellsX*(nNum1DZones+nNumGhostCells);//may be need to +1 if only one proc and variable in interface centred
    if (nVarInfo[n][0]==1&&nPeriodic[0]==0){
      nSizeX1=nGhostCellsX*(nNum1DZones+1+nNumGhostCells);
    }
    nSizeY=1;
    nSizeZ=1;
    for(i=0;i<nSizeX1;i++){
      //find average max, and min
      dMax[n][i]=dGrid[n][i][0][0];
      dMin[n][i]=dGrid[n][i][0][0];
      dAve[n][i]=dGrid[n][i][0][0];
      nMaxJIndex[n][i]=0;
      nMaxKIndex[n][i]=0;
      nMinJIndex[n][i]=0;
      nMinKIndex[n][i]=0;
    }
    
    //read in the rest of the grid
    nSizeX2=nSize[n][0]+nGhostCellsX*2*nNumGhostCells;
    nSizeY=nSize[n][1]+nGhostCellsY*2*nNumGhostCells;
    nSizeZ=nSize[n][2]+nGhostCellsZ*2*nNumGhostCells;
    nStartY=nGhostCellsY*nNumGhostCells;
    nEndY=nSize[n][1]+nStartY;
    nStartZ=nGhostCellsZ*nNumGhostCells;
    nEndZ=nSize[n][2]+nStartZ;
    for(i=nSizeX1;i<nSizeX2;i++){
      dMaxTemp=-1.0*std::numeric_limits<double>::max();
      dMinTemp=std::numeric_limits<double>::max();
      dSum=0.0;
      nCount=0;
      for(j=nStartY;j<nEndY;j++){
        for(k=nStartZ;k<nEndZ;k++){
          //find average max, and min
          if(dGrid[n][i][j][k]>dMaxTemp){
            dMaxTemp=dGrid[n][i][j][k];
            nMaxJIndex[n][i]=j;
            nMaxKIndex[n][i]=k;
          }
          if(dGrid[n][i][j][k]<dMinTemp){
            dMinTemp=dGrid[n][i][j][k];
            nMinJIndex[n][i]=j;
            nMinKIndex[n][i]=k;
          }
          dSum+=dGrid[n][i][j][k];
          nCount++;
        }
      }
      dMax[n][i]=dMaxTemp;
      dMin[n][i]=dMinTemp;
      dAve[n][i]=dSum/double(nCount);
    }
  }
  
  //allocate space for internal variables
  for(int n=nNumVars;n<nNumVars+nNumIntVars;n++){
    dMax[n]=new double[nSize[bin.nD][0]+nGhostCellsX*2*nNumGhostCells];
    dMin[n]=new double[nSize[bin.nD][0]+nGhostCellsX*2*nNumGhostCells];
    dAve[n]=new double[nSize[bin.nD][0]+nGhostCellsX*2*nNumGhostCells];
    nMaxJIndex[n]=new int[nSize[bin.nD][0]+nGhostCellsX*2*nNumGhostCells];
    nMaxKIndex[n]=new int[nSize[bin.nD][0]+nGhostCellsX*2*nNumGhostCells];
    nMinJIndex[n]=new int[nSize[bin.nD][0]+nGhostCellsX*2*nNumGhostCells];
    nMinKIndex[n]=new int[nSize[bin.nD][0]+nGhostCellsX*2*nNumGhostCells];
  }
  
  //calculate filling factor for upflow
  nGhostCellsX=1;
  if(nVarInfo[bin.nU][0]==-1){
    nGhostCellsX=0;
  }
  nGhostCellsY=1;
  if(nVarInfo[bin.nU][1]==-1){
    nGhostCellsY=0;
  }
  nGhostCellsZ=1;
  if(nVarInfo[bin.nU][2]==-1){
    nGhostCellsZ=0;
  }
  nSizeX1=nGhostCellsX*(nNum1DZones+nNumGhostCells);
  if (nVarInfo[bin.nU][0]==1&&nPeriodic[0]==0){
    nSizeX1=nGhostCellsX*(nNum1DZones+1+nNumGhostCells);
  }
  nSizeX2=nSize[bin.nU][0]+nGhostCellsX*2*nNumGhostCells;
  nSizeY=nSize[bin.nU][1]+nGhostCellsY*2*nNumGhostCells;
  nSizeZ=nSize[bin.nU][2]+nGhostCellsZ*2*nNumGhostCells;
  nStartY=nGhostCellsY*nNumGhostCells;
  nEndY=nSize[bin.nU][1]+nStartY;
  nStartZ=nGhostCellsZ*nNumGhostCells;
  nEndZ=nSize[bin.nU][2]+nStartZ;
  double *dUpFlowFillingFactor=new double[nSize[bin.nU][0]+nGhostCellsX*2
    *nNumGhostCells];
  for(int i=0;i<nSizeX1;i++){
    dUpFlowFillingFactor[i]=0.0;//no-upflow in 1D part
  }
  int nCountUp=0;
  int nCountTotal=0;
  for(int i=nSizeX1;i<nSizeX2;i++){
    for(int j=nStartY;j<nEndY;j++){
      for(int k=nStartZ;k<nEndZ;k++){
        
        //if there is an up-flow (removing flow due to pulsation)
        if(dGrid[bin.nU][i][j][k]-dGrid[bin.nU0][i][0][0]>0.0){
          nCountUp++;
        }
        nCountTotal++;
      }
    }
    
    //record filling factor, and reset counts
    dUpFlowFillingFactor[i]=double(nCountUp)/double(nCountTotal);
    nCountUp=0;
    nCountTotal=0;
  }
  
  if(nGammaLaw!=0){//set P,E,kappa,gamma, Q, L_rad and L_con, KE, C, <rho>
    
    //allocate space
    nGhostCellsX=1;
    if(nVarInfo[bin.nD][0]==-1){/*all internal variables are centred quantities, will be the same as 
      the density*/
      nGhostCellsX=0;
    }
    nGhostCellsY=1;
    if(nVarInfo[bin.nD][1]==-1){
      nGhostCellsY=0;
    }
    nGhostCellsZ=1;
    if(nVarInfo[bin.nD][2]==-1){
      nGhostCellsZ=0;
    }
    
    double dF_con;
    double dP_i;
    double dE_i;
    double dKappa_i;
    double dGamma_i;
    double dCp_i;
    double dCp_ip1;
    double dCp_ip1half;
    double dRho_ip1half;
    double dTAve_i;
    double dTAve_ip1;
    double dTAve_ip1half;
    double dT_ip1halfjk;
    double dP_ip1;
    double dE_ip1;
    double dKappa_ip1;
    double dKappa_ip1half;
    double dGamma_ip1;
    double dQ;
    double dQ0=0.0;
    double dQ1=0.0;
    double dQ2=0.0;
    double dMaxE;
    double dMinE;
    double dSumE;
    double dMaxP;
    double dMinP;
    double dSumP;
    double dMaxF_con;
    double dMinF_con;
    double dSumF_con;
    double dMaxKappa;
    double dMinKappa;
    double dSumKappa;
    double dMaxGamma;
    double dMinGamma;
    double dSumGamma;
    double dMaxQ;
    double dMinQ;
    double dSumQ;
    double dC;
    double dMaxC;
    double dMinC;
    double dSumC;
    double dRSq_i;
    double dA_ip1half;
    double dA_im1half;
    double dTheta_j;
    double dTheta_jp1half;
    double dTheta_jm1half;
    double dA_jp1half;
    double dA_jm1half;
    double dA_j;
    double dASq=dA*dA;
    double dDVDtThreshold;
    double dDVDt_mthreshold;
    double dDVDt;
    double dLSum_rad=0.0;
    double dLSum_con=0.0;
    double dAreaSum=0.0;
    double dR_i;
    double dRSq_ip1half;
    double dT4_ip1;
    double dT4_i;
    double dArea;
    double dArea1;
    double dArea2;
    double dPhi_kp1half;
    double dPhi_km1half;
    double dU_i;
    double dU_ijk;
    double dV_ijk;
    double dW_ijk;
    double dVSq_ijk;
    double dKETotSum;
    double dDMSum;
    double dDMTemp;
    
    //set 1D part of the grid
    nSizeX1=nGhostCellsX*(nNum1DZones+nNumGhostCells);/*maybe need to +1 if only one proc and 
      variable is interface centered*/
    if (nVarInfo[bin.nD][0]==1&&nPeriodic[0]==0){
      nSizeX1=nGhostCellsX*(nNum1DZones+1+nNumGhostCells);
    }
    nSizeY=1;
    nSizeZ=1;
    for(i=0;i<nSizeX1;i++){//find average max, and min in 1D region
      
      //get P,E,Kappa,Gamma
      eosTable->getPEKappaGamma(dGrid[bin.nT][i][0][0],dGrid[bin.nD][i][0][0],dP_i,dE_i,dKappa_i,dGamma_i);
      eosTable->getPEKappaGamma(dGrid[bin.nT][i+1][0][0],dGrid[bin.nD][i+1][0][0],dP_ip1,dE_ip1,dKappa_ip1
        ,dGamma_ip1);
      
      //calculate Q
      dR_i=(dGrid[bin.nR][i+1][0][0]+dGrid[bin.nR][i][0][0])*0.5;
      dRSq_i=dR_i*dR_i;
      dA_ip1half=dGrid[bin.nR][i+1][0][0]*dGrid[bin.nR][i+1][0][0];
      dA_im1half=dGrid[bin.nR][i][0][0]*dGrid[bin.nR][i][0][0];
      dC=sqrt(dGamma_i*dP_i/dGrid[bin.nD][i][0][0]);
      dDVDtThreshold=dAVThreshold*dC;
      dDVDt=(dA_ip1half*dGrid[bin.nU][i+1][0][0]
        -dA_im1half*dGrid[bin.nU][i][0][0])/dRSq_i;
      if(dDVDt<-1.0*dDVDtThreshold){//being compressed
        dDVDt_mthreshold=dDVDt+dDVDtThreshold;
        dQ=dASq*dGrid[bin.nD][i][0][0]*dDVDt_mthreshold*dDVDt_mthreshold;
      }
      else{
        dQ=0.0;
      }
      
      //calculate luminosity from cell and add to sum
      dLSum_rad=0.0;
      dArea=dA_ip1half*4.0*bin.dPi;
      dT4_i=pow(dGrid[bin.nT][i][0][0],4);
      dT4_ip1=pow(dGrid[bin.nT][i+1][0][0],4);
      dKappa_ip1half=(dT4_i/dKappa_i+dT4_ip1/dKappa_ip1)
        /(dT4_i+dT4_ip1);
      dAve[bin.nL_rad][i]=-16.0*bin.dPi*bin.dSigma*dA_ip1half*dKappa_ip1half
        /3.0*(dT4_ip1-dT4_i)/(dGrid[bin.nDM][i][0][0]+dGrid[bin.nDM][i+1][0][0])
        *2.0*dArea;
      dAve[bin.nL_con][i]=0.0;
      dU_i=(dGrid[bin.nU0][i+1][0][0]+dGrid[bin.nU0][i][0][0])*0.5;
      dAve[bin.nKEP][i]=0.5*dGrid[bin.nDM][i][0][0]*dU_i*dU_i;
      dAve[bin.nKETot][i]=dAve[bin.nKEP][i];//nothing but pulsation energy in 1D
      dMax[bin.nE][i]=dE_i;
      dMin[bin.nE][i]=dE_i;
      dAve[bin.nE][i]=dE_i;
      dMax[bin.nP][i]=dP_i;
      dMin[bin.nP][i]=dP_i;
      dAve[bin.nP][i]=dP_i;
      dMax[bin.nF_con][i]=0.0;
      dMin[bin.nF_con][i]=0.0;
      dAve[bin.nF_con][i]=0.0;
      dMax[bin.nKappa][i]=dKappa_i;
      dMin[bin.nKappa][i]=dKappa_i;
      dAve[bin.nKappa][i]=dKappa_i;
      dMax[bin.nGamma][i]=dGamma_i;
      dMin[bin.nGamma][i]=dGamma_i;
      dAve[bin.nGamma][i]=dGamma_i;
      dMax[bin.nQ][i]=dQ;
      dMin[bin.nQ][i]=dQ;
      dAve[bin.nQ][i]=dQ;
      dMax[bin.nC][i]=dC;
      dMin[bin.nC][i]=dC;
      dAve[bin.nC][i]=dC;
      nMaxJIndex[bin.nP][i]=0;
      nMaxKIndex[bin.nP][i]=0;
      nMinJIndex[bin.nP][i]=0;
      nMinKIndex[bin.nP][i]=0;
      nMaxJIndex[bin.nE][i]=0;
      nMaxKIndex[bin.nE][i]=0;
      nMinJIndex[bin.nE][i]=0;
      nMinKIndex[bin.nE][i]=0;
      nMaxJIndex[bin.nKappa][i]=0;
      nMaxKIndex[bin.nKappa][i]=0;
      nMinJIndex[bin.nKappa][i]=0;
      nMinKIndex[bin.nKappa][i]=0;
      nMaxJIndex[bin.nGamma][i]=0;
      nMaxKIndex[bin.nGamma][i]=0;
      nMinJIndex[bin.nGamma][i]=0;
      nMinKIndex[bin.nGamma][i]=0;
      nMaxJIndex[bin.nQ][i]=0;
      nMaxKIndex[bin.nQ][i]=0;
      nMinJIndex[bin.nQ][i]=0;
      nMinKIndex[bin.nQ][i]=0;
      nMaxJIndex[bin.nC][i]=0;
      nMaxKIndex[bin.nC][i]=0;
      nMinJIndex[bin.nC][i]=0;
      nMinKIndex[bin.nC][i]=0;
      nMaxJIndex[bin.nL_rad][i]=0;
      nMaxKIndex[bin.nL_rad][i]=0;
      nMinJIndex[bin.nL_rad][i]=0;
      nMinKIndex[bin.nL_rad][i]=0;
      nMaxJIndex[bin.nKEP][i]=0;
      nMaxKIndex[bin.nKEP][i]=0;
      nMinJIndex[bin.nKEP][i]=0;
      nMinKIndex[bin.nKEP][i]=0;
      nMaxJIndex[bin.nKETot][i]=0;
      nMaxKIndex[bin.nKETot][i]=0;
      nMinJIndex[bin.nKETot][i]=0;
      nMinKIndex[bin.nKETot][i]=0;
    }
    
    //set the rest of the grid
    nSizeX2=nSize[bin.nD][0]+nGhostCellsX*2*nNumGhostCells;
    nSizeY=nSize[bin.nD][1]+nGhostCellsY*2*nNumGhostCells;
    nSizeZ=nSize[bin.nD][2]+nGhostCellsZ*2*nNumGhostCells;
    nStartY=nGhostCellsY*nNumGhostCells;
    nEndY=nSize[bin.nD][1]+nStartY;
    nStartZ=nGhostCellsZ*nNumGhostCells;
    nEndZ=nSize[bin.nD][2]+nStartZ;
    for(i=nSizeX1;i<nSizeX2;i++){
      dMaxE=-1.0*std::numeric_limits<double>::max();
      dMinE=std::numeric_limits<double>::max();
      dSumE=0.0;
      dMaxP=-1.0*std::numeric_limits<double>::max();
      dMinP=std::numeric_limits<double>::max();
      dSumP=0.0;
      dMaxF_con=-1.0*std::numeric_limits<double>::max();
      dMinF_con=std::numeric_limits<double>::max();
      dSumF_con=0.0;
      dMaxKappa=-1.0*std::numeric_limits<double>::max();
      dMinKappa=std::numeric_limits<double>::max();
      dSumKappa=0.0;
      dMaxGamma=-1.0*std::numeric_limits<double>::max();
      dMinGamma=std::numeric_limits<double>::max();
      dSumGamma=0.0;
      dMaxQ=-1.0*std::numeric_limits<double>::max();
      dMinQ=std::numeric_limits<double>::max();
      dSumQ=0.0;
      dMaxC=-1.0*std::numeric_limits<double>::max();
      dMinC=std::numeric_limits<double>::max();
      dSumC=0.0;
      nCount=0;
      dR_i=(dGrid[bin.nR][i+1][0][0]+dGrid[bin.nR][i][0][0])*0.5;
      dRSq_ip1half=dGrid[bin.nR][i+1][0][0]*dGrid[bin.nR][i+1][0][0];
      dRSq_i=dR_i*dR_i;
      dA_ip1half=dGrid[bin.nR][i+1][0][0]*dGrid[bin.nR][i+1][0][0];
      dA_im1half=dGrid[bin.nR][i][0][0]*dGrid[bin.nR][i][0][0];
      dLSum_rad=0.0;
      dLSum_con=0.0;
      dAreaSum=0.0;
      dArea1=dRSq_ip1half*4.0*bin.dPi;
      dKETotSum=0.0;
      dDMSum=0.0;
      dDMTemp=dGrid[bin.nDM][i][0][0];//in 1D these are the same
      
      //start debug LCon
      /*
      double bDebugLCon=false;
      std::ofstream ofDebugLCon;
      int nWidthOutputField=25;
      int nWidthIntOutputField=12;
      if(i==125 || i==128){
        bDebugLCon=true;
        std::stringstream ssDebugFile;
        ssDebugFile<<"debugLCon_i"<<i<<".txt";
        std::cout<<ssDebugFile.str().c_str()<<std::endl;
        ofDebugLCon.open(ssDebugFile.str().c_str());
        
        //set double output precision
        ofDebugLCon.precision(16);
        ofDebugLCon.unsetf(std::ios::fixed);
        ofDebugLCon.setf(std::ios::scientific);
        
        ofDebugLCon
          <<std::setw(nWidthIntOutputField)<<"j(0)"
          <<std::setw(nWidthIntOutputField)<<"k(1)"
          <<std::setw(nWidthOutputField)<<"dCp_ip1half(2)"
          <<std::setw(nWidthOutputField)<<"dRho_ip1half(3)"
          <<std::setw(nWidthOutputField)<<"dT_ip1halfjk(4)"
          <<std::setw(nWidthOutputField)<<"dTAve_ip1half(5)"
          <<std::setw(nWidthOutputField)<<"dU_ip1halfjk(6)"
          <<std::setw(nWidthOutputField)<<"dU0_ip1half(7)"
          <<std::setw(nWidthOutputField)<<"dA_ip1halfjk(8)"
          <<std::setw(nWidthOutputField)<<"dLCon_ip1halfjk(9)"<<std::endl;
      }*/
      //end debug LCon
      
      for(j=nStartY;j<nEndY;j++){
        
        if(nNumDims>1){
          dTheta_jp1half=dGrid[bin.nTheta][0][j][0];
          if(j==0){
            dTheta_jm1half=dGrid[bin.nTheta][0][j][0]-(dGrid[bin.nTheta][0][j+1][0]-dGrid[bin.nTheta][0][j][0]);
          }
          else{
            dTheta_jm1half=dGrid[bin.nTheta][0][j-1][0];
          }
          dTheta_j=(dTheta_jp1half+dTheta_jm1half)*0.5;
          dA_jp1half=sin(dTheta_jp1half);
          dA_jm1half=sin(dTheta_jm1half);
          dA_j=sin(dTheta_j);
          dArea2=dArea1/2.0*(cos(dTheta_jm1half)-cos(dTheta_jp1half));
          dArea=dArea2;
          dDMTemp=(cos(dTheta_jm1half)-cos(dTheta_jp1half))*dGrid[bin.nD][i][j][0];
        }
        
        for(k=nStartZ;k<nEndZ;k++){
          if(nNumDims>2){
            dPhi_kp1half=dGrid[bin.nPhi][0][0][k];
            if(k==0){
              dPhi_km1half=dGrid[bin.nPhi][0][0][k]-(dGrid[bin.nPhi][0][0][k+1]-dGrid[bin.nPhi][0][0][k]);
            }
            else{
              dPhi_km1half=dGrid[bin.nPhi][0][0][k-1];
            }
            dArea=dArea2/(2.0*bin.dPi)*(dPhi_kp1half-dPhi_km1half);
            dDMTemp=(cos(dTheta_jm1half)-cos(dTheta_jp1half))
              *(dPhi_kp1half-dPhi_km1half)*dGrid[bin.nD][i][j][k];
          }
          
          //get P,E,Kappa,Gamma, calculate luminosity from cell and add to sum
          eosTable->getPEKappaGammaCp(dGrid[bin.nT][i][j][k],dGrid[bin.nD][i][j][k],dP_i,dE_i,dKappa_i
            ,dGamma_i,dCp_i);
          if(i<nSizeX2-1){
            eosTable->getPEKappaGammaCp(dGrid[bin.nT][i+1][j][k],dGrid[bin.nD][i+1][j][k],dP_ip1,dE_ip1
              ,dKappa_ip1,dGamma_ip1,dCp_ip1);
            dT4_i=pow(dGrid[bin.nT][i][j][k],4);
            dT4_ip1=pow(dGrid[bin.nT][i+1][j][k],4);
            dKappa_ip1half=(dT4_i/dKappa_i+dT4_ip1/dKappa_ip1)
              /(dT4_i+dT4_ip1);
            dLSum_rad=dLSum_rad-16.0*bin.dPi*bin.dSigma*dRSq_ip1half*dKappa_ip1half/3.0
              *(dT4_ip1-dT4_i)/(dGrid[bin.nDM][i][0][0]+dGrid[bin.nDM][i+1][0][0])*2.0*dArea;
            dRho_ip1half=(dGrid[bin.nD][i][j][k]+dGrid[bin.nD][i+1][j][k])*0.5;
            dTAve_ip1half=(dAve[bin.nT][i]+dAve[bin.nT][i+1])*0.5;
            dT_ip1halfjk=(dGrid[bin.nT][i][j][k]+dGrid[bin.nT][i+1][j][k])*0.5;
            dCp_ip1half=(dCp_i+dCp_ip1)*0.5;
            
            //start debug LCon
            /*
            if (bDebugLCon){
              double dLCon=dCp_ip1half*dRho_ip1half*(dT_ip1halfjk-dTAve_ip1half)
                *(dGrid[bin.nU][i+1][j][k]-dGrid[bin.nU0][i+1][0][0])*dArea;
              ofDebugLCon
                <<std::setw(nWidthIntOutputField)<<j
                <<std::setw(nWidthIntOutputField)<<k
                <<std::setw(nWidthOutputField)<<dCp_ip1half
                <<std::setw(nWidthOutputField)<<dRho_ip1half
                <<std::setw(nWidthOutputField)<<dT_ip1halfjk
                <<std::setw(nWidthOutputField)<<dTAve_ip1half
                <<std::setw(nWidthOutputField)<<dGrid[bin.nU][i+1][j][k]
                <<std::setw(nWidthOutputField)<<dGrid[bin.nU0][i+1][0][0]
                <<std::setw(nWidthOutputField)<<dArea
                <<std::setw(nWidthOutputField)<<dLCon<<std::endl;
            }*/
            //end debug LCon
            
            dLSum_con=dLSum_con+dCp_ip1half*dRho_ip1half*(dT_ip1halfjk-dTAve_ip1half)
              *(dGrid[bin.nU][i+1][j][k]-dGrid[bin.nU0][i+1][0][0])*dArea;
            dF_con=dCp_ip1half*dRho_ip1half*(dT_ip1halfjk-dTAve_ip1half)
              *(dGrid[bin.nU][i+1][j][k]-dGrid[bin.nU0][i+1][0][0]);
          }
          else{
            //use surface boundary condition
            dLSum_rad=dLSum_rad+bin.dSigma*dArea*pow(pow(2.0,0.25)*dGrid[bin.nT][i][j][k],4);
            dLSum_con=0.0;
          }
          dAreaSum+=dArea;
          
          //calculate total kinetic energy
          if(i==nSizeX1){
            dU_ijk=(dGrid[bin.nU][i+1][j][k]+dGrid[bin.nU][i][0][0])*0.5;
          }
          else{
            dU_ijk=(dGrid[bin.nU][i+1][j][k]+dGrid[bin.nU][i][j][k])*0.5;
          }
          dVSq_ijk=dU_ijk*dU_ijk;
          if(nNumDims>1){
            dV_ijk=(dGrid[bin.nV][i][j+1][k]+dGrid[bin.nV][i][j][k])*0.5;
            dVSq_ijk+=dV_ijk*dV_ijk;
          }
          if(nNumDims>2){
            dW_ijk=(dGrid[bin.nW][i][j][k+1]+dGrid[bin.nW][i][j][k])*0.5;
            dVSq_ijk+=dW_ijk*dW_ijk;
          }
          dKETotSum+=dDMTemp*dVSq_ijk;
          dDMSum+=dDMTemp;
          
          if(dF_con>dMaxF_con){
            dMaxF_con=dF_con;
            nMaxJIndex[bin.nF_con][i]=j;
            nMaxKIndex[bin.nF_con][i]=k;
          }
          if(dF_con<dMinF_con){
            dMinF_con=dF_con;
            nMinJIndex[bin.nF_con][i]=j;
            nMinKIndex[bin.nF_con][i]=k;
          }
          dSumF_con+=dF_con;
          
          //calculate Q
          dC=sqrt(dGamma_i*dP_i/dGrid[bin.nD][i][j][k]);
          dDVDtThreshold=dAVThreshold*dC;
          if(nNumDims>=1){
            if(i==nSizeX1){
              dDVDt=(dA_ip1half*dGrid[bin.nU][i+1][j][k]
                -dA_im1half*dGrid[bin.nU][i][0][0])/dRSq_i;
            }
            else{
              dDVDt=(dA_ip1half*dGrid[bin.nU][i+1][j][k]
                -dA_im1half*dGrid[bin.nU][i][j][k])/dRSq_i;
            }
            if(dDVDt<-1.0*dDVDtThreshold){//being compressed
              dDVDt_mthreshold=dDVDt+dDVDtThreshold;
              dQ0=dASq*dGrid[bin.nD][i][j][k]*dDVDt_mthreshold*dDVDt_mthreshold;
            }
            else{
              dQ0=0.0;
            }
          }
          else if(nNumDims>=2){
            if(j==0){
              dDVDt=(dA_jp1half*dGrid[bin.nV][i][j][k]
                -dA_jm1half*dGrid[bin.nV][i][nSizeY-1][k])/dA_j;
            }
            else{
              dDVDt=(dA_jp1half*dGrid[bin.nV][i][j][k]
                -dA_jm1half*dGrid[bin.nV][i][j-1][k])/dA_j;
            }
            if(dDVDt<-1.0*dDVDtThreshold){//being compressed
              dDVDt_mthreshold=dDVDt+dDVDtThreshold;
              dQ1=dASq*dGrid[bin.nD][i][j][k]
                *dDVDt_mthreshold*dDVDt_mthreshold;
            }
            else{
              dQ1=0.0;
            }
          }
          else if(nNumDims==3){
            if(k==0){
              dDVDt=(dGrid[bin.nW][i][j][k]-dGrid[bin.nW][i][j][nSizeZ-1]);
            }
            else{
              dDVDt=(dGrid[bin.nW][i][j][k]-dGrid[bin.nW][i][j][k-1]);
            }
            if(dDVDt<-1.0*dDVDtThreshold){//being compressed
              dDVDt_mthreshold=dDVDt+dDVDtThreshold;
              dQ2=dASq*dGrid[bin.nD][i][j][k]*dDVDt_mthreshold*dDVDt_mthreshold;
            }
            else{
              dQ2=0.0;
            }
          }
          dQ=dQ0+dQ1+dQ2;
          
          if(dP_i>dMaxP){
            dMaxP=dP_i;
            nMaxJIndex[bin.nP][i]=j;
            nMaxKIndex[bin.nP][i]=k;
          }
          if(dP_i<dMinP){
            dMinP=dP_i;
            nMinJIndex[bin.nP][i]=j;
            nMinKIndex[bin.nP][i]=k;
          }
          dSumP+=dP_i;
          if(dE_i>dMaxE){
            dMaxE=dE_i;
            nMaxJIndex[bin.nE][i]=j;
            nMaxKIndex[bin.nE][i]=k;
          }
          if(dE_i<dMinE){
            dMinE=dE_i;
            nMinJIndex[bin.nE][i]=j;
            nMinKIndex[bin.nE][i]=k;
          }
          dSumE+=dE_i;
          if(dKappa_i>dMaxKappa){
            dMaxKappa=dKappa_i;
            nMaxJIndex[bin.nKappa][i]=j;
            nMaxKIndex[bin.nKappa][i]=k;
          }
          if(dKappa_i<dMinKappa){
            dMinKappa=dKappa_i;
            nMinJIndex[bin.nKappa][i]=j;
            nMinKIndex[bin.nKappa][i]=k;
          }
          dSumKappa+=dKappa_i;
          if(dGamma_i>dMaxGamma){
            dMaxGamma=dGamma_i;
            nMaxJIndex[bin.nGamma][i]=j;
            nMaxKIndex[bin.nGamma][i]=k;
          }
          if(dGamma_i<dMinGamma){
            dMinGamma=dGamma_i;
            nMinJIndex[bin.nGamma][i]=j;
            nMinKIndex[bin.nGamma][i]=k;
          }
          dSumGamma+=dGamma_i;
          if(dQ>dMaxQ){
            dMaxQ=dQ;
            nMaxJIndex[bin.nQ][i]=j;
            nMaxKIndex[bin.nQ][i]=k;
          }
          if(dQ<dMinQ){
            dMinQ=dQ;
            nMinJIndex[bin.nQ][i]=j;
            nMinKIndex[bin.nQ][i]=k;
          }
          dSumQ+=dQ;
          if(dC>dMaxC){
            dMaxC=dC;
            nMaxJIndex[bin.nC][i]=j;
            nMaxKIndex[bin.nC][i]=k;
          }
          if(dC<dMinC){
            dMinC=dC;
            nMinJIndex[bin.nC][i]=j;
            nMinKIndex[bin.nC][i]=k;
          }
          dSumC+=dC;
          nMaxJIndex[bin.nL_rad][i]=0;
          nMaxKIndex[bin.nL_rad][i]=0;
          nMinJIndex[bin.nL_rad][i]=0;
          nMinKIndex[bin.nL_rad][i]=0;
          nMaxJIndex[bin.nKEP][i]=0;
          nMaxKIndex[bin.nKEP][i]=0;
          nMinJIndex[bin.nKEP][i]=0;
          nMinKIndex[bin.nKETot][i]=0;
          nMaxJIndex[bin.nKETot][i]=0;
          nMaxKIndex[bin.nKETot][i]=0;
          nMinJIndex[bin.nKETot][i]=0;
          nMinKIndex[bin.nKETot][i]=0;
          nCount++;
        }
      }
      if(nNumDims==3){
        dAve[bin.nD][i]=bin.dCalRhoAve3D(dGrid,i,nStartY,nEndY,nStartZ,nEndZ);
      }
      else if(nNumDims==2){
        dAve[bin.nD][i]=bin.dCalRhoAve2D(dGrid,i,nStartY,nEndY,nStartZ,nEndZ);
      }
      else{
        dAve[bin.nD][i]=dGrid[bin.nD][i][0][0];
      }
      
      //start debug LCon
      /*
      if (bDebugLCon){
        ofDebugLCon
          <<std::setw(nWidthOutputField)<<dLSum_con
          <<std::setw(nWidthOutputField)<<dAreaSum
          <<std::setw(nWidthOutputField)<<dPi
          <<std::setw(nWidthOutputField)<<dRSq_ip1half
          <<std::endl;
        ofDebugLCon.close();
        bDebugLCon=false;
      }*/
      //end debug LCon
      
      dAve[bin.nL_rad][i]=dLSum_rad/dAreaSum*4.0*bin.dPi*dRSq_ip1half;
      dAve[bin.nL_con][i]=dLSum_con/dAreaSum*4.0*bin.dPi*dRSq_ip1half;
      dU_i=(dGrid[bin.nU0][i+1][0][0]+dGrid[bin.nU0][i][0][0])*0.5;
      dAve[bin.nKEP][i]=0.5*dGrid[bin.nDM][i][0][0]*dU_i*dU_i;
      dAve[bin.nKETot][i]=0.5*dKETotSum/dDMSum*dGrid[bin.nDM][i][0][0];
      dMax[bin.nP][i]=dMaxP;
      dMin[bin.nP][i]=dMinP;
      dAve[bin.nP][i]=dSumP/double(nCount);
      dMax[bin.nF_con][i]=dMaxF_con;
      dMin[bin.nF_con][i]=dMinF_con;
      dAve[bin.nF_con][i]=dSumF_con/double(nCount);
      dMax[bin.nE][i]=dMaxE;
      dMin[bin.nE][i]=dMinE;
      dAve[bin.nE][i]=dSumE/double(nCount);
      dMax[bin.nKappa][i]=dMaxKappa;
      dMin[bin.nKappa][i]=dMinKappa;
      dAve[bin.nKappa][i]=dSumKappa/double(nCount);
      dMax[bin.nGamma][i]=dMaxGamma;
      dMin[bin.nGamma][i]=dMinGamma;
      dAve[bin.nGamma][i]=dSumGamma/double(nCount);
      dMax[bin.nQ][i]=dMaxQ;
      dMin[bin.nQ][i]=dMinQ;
      dAve[bin.nQ][i]=dSumQ/double(nCount);
      dMax[bin.nC][i]=dMaxC;
      dMin[bin.nC][i]=dMinC;
      dAve[bin.nC][i]=dSumC/double(nCount);
    }
  }
  else{//set P, Q, KE, and C
    
    //allocate space
    nGhostCellsX=1;
    if(nVarInfo[bin.nD][0]==-1){/*all internal variables are cetnered quantities, will be the same as 
      the density*/
      nGhostCellsX=0;
    }
    nGhostCellsY=1;
    if(nVarInfo[bin.nD][1]==-1){
      nGhostCellsY=0;
    }
    nGhostCellsZ=1;
    if(nVarInfo[bin.nD][2]==-1){
      nGhostCellsZ=0;
    }
    double dP;
    double dQ;
    double dMaxP;
    double dMinP;
    double dSumP;
    double dMaxQ;
    double dMinQ;
    double dSumQ;
    double dMaxC;
    double dMinC;
    double dSumC;
    double dQ0=0.0;
    double dQ1=0.0;
    double dQ2=0.0;
    double dRSq_i;
    double dA_ip1half;
    double dA_im1half;
    double dTheta_j;
    double dTheta_jp1half;
    double dTheta_jm1half;
    double dA_jp1half;
    double dA_jm1half;
    double dA_j;
    double dASq=dA*dA;
    double dC;
    double dDVDtThreshold;
    double dDVDt_mthreshold;
    double dDVDt;
    double dU_i;
    
    //set 1D part of the grid
    nSizeX1=nGhostCellsX*(nNum1DZones+nNumGhostCells);//may be need to +1 if only one proc and variable in interface centered
    if (nVarInfo[bin.nD][0]==1&&nPeriodic[0]==0){
      nSizeX1=nGhostCellsX*(nNum1DZones+1+nNumGhostCells);
    }
    nSizeY=1;
    nSizeZ=1;
    for(i=0;i<nSizeX1;i++){//find average max, and min in 1D region
      
      //calculate P
      dP=dGrid[bin.nD][i][0][0]*(dGamma-1.0)*dGrid[bin.nE][i][0][0];
      
      //calculate Q
      dRSq_i=(dGrid[bin.nR][i+1][0][0]+dGrid[bin.nR][i][0][0])*0.5;
      dA_ip1half=dGrid[bin.nR][i+1][0][0]*dGrid[bin.nR][i+1][0][0];
      dA_im1half=dGrid[bin.nR][i][0][0]*dGrid[bin.nR][i][0][0];
      dC=sqrt(dGamma*dP/dGrid[bin.nD][i][j][k]);
      dDVDtThreshold=dAVThreshold*dC;
      dDVDt=(dA_ip1half*dGrid[bin.nU][i+1][j][k]
        -dA_im1half*dGrid[bin.nU][i][j][k])/dRSq_i;
      if(dDVDt<-1.0*dDVDtThreshold){//being compressed
        dDVDt_mthreshold=dDVDt+dDVDtThreshold;
        dQ=dASq*dGrid[bin.nD][i][j][k]*dDVDt_mthreshold*dDVDt_mthreshold;
      }
      else{
        dQ=0.0;
      }
      dU_i=(dGrid[bin.nU0][i+1][0][0]+dGrid[bin.nU0][i][0][0])*0.5;
      dAve[bin.nKEP][i]=0.5*dGrid[bin.nDM][i][0][0]*dU_i*dU_i;
      dMax[bin.nP][i]=dP;
      dMin[bin.nP][i]=dP;
      dAve[bin.nP][i]=dP;
      dMax[bin.nQ][i]=dQ;
      dMin[bin.nQ][i]=dQ;
      dAve[bin.nQ][i]=dQ;
      dMax[bin.nC][i]=dC;
      dMin[bin.nC][i]=dC;
      dAve[bin.nC][i]=dC;
      nMaxJIndex[bin.nP][i]=0;
      nMaxKIndex[bin.nP][i]=0;
      nMinJIndex[bin.nP][i]=0;
      nMinKIndex[bin.nP][i]=0;
      nMaxJIndex[bin.nQ][i]=0;
      nMaxKIndex[bin.nQ][i]=0;
      nMinJIndex[bin.nQ][i]=0;
      nMinKIndex[bin.nQ][i]=0;
      nMaxJIndex[bin.nC][i]=0;
      nMaxKIndex[bin.nC][i]=0;
      nMinJIndex[bin.nC][i]=0;
      nMinKIndex[bin.nC][i]=0;
      nMaxJIndex[bin.nKEP][i]=0;
      nMaxKIndex[bin.nKEP][i]=0;
      nMinJIndex[bin.nKEP][i]=0;
      nMinKIndex[bin.nKEP][i]=0;
      
    }
    
    //set the rest of the grid
    nSizeX2=nSize[bin.nD][0]+nGhostCellsX*2*nNumGhostCells;
    nSizeY=nSize[bin.nD][1]+nGhostCellsY*2*nNumGhostCells;
    nSizeZ=nSize[bin.nD][2]+nGhostCellsZ*2*nNumGhostCells;
    nStartY=nGhostCellsY*nNumGhostCells;
    nEndY=nSize[bin.nD][1]+nStartY;
    nStartZ=nGhostCellsZ*nNumGhostCells;
    nEndZ=nSize[bin.nD][2]+nStartZ;
    for(i=nSizeX1;i<nSizeX2;i++){
      dMaxP=-1.0*std::numeric_limits<double>::max();
      dMinP=std::numeric_limits<double>::max();
      dSumP=0.0;
      dMaxQ=-1.0*std::numeric_limits<double>::max();
      dMinQ=std::numeric_limits<double>::max();
      dSumQ=0.0;
      dMaxC=-1.0*std::numeric_limits<double>::max();
      dMinC=std::numeric_limits<double>::max();
      dSumC=0.0;
      nCount=0;
      for(j=nStartY;j<nEndY;j++){
        dTheta_jp1half=dGrid[bin.nTheta][0][j][0];
        if(j==0){
          dTheta_jm1half=dGrid[bin.nTheta][0][j][0]-(dGrid[bin.nTheta][0][j+1][0]-dGrid[bin.nTheta][0][j][0]);
        }
        dTheta_j=(dTheta_jp1half+dTheta_jm1half)*0.5;
        dA_jp1half=sin(dTheta_jp1half);
        dA_jm1half=sin(dTheta_jm1half);;
        dA_j=sin(dTheta_j);
        for(k=nStartZ;k<nEndZ;k++){
          
          //calculate P
          dP=dGrid[bin.nD][i][j][k]*(dGamma-1.0)*dGrid[bin.nE][i][j][k];
        
          //calculate Q
          dC=sqrt(dGamma*dP/dGrid[bin.nD][i][j][k]);
          dDVDtThreshold=dAVThreshold*dC;
          if(nNumDims>=1){
            dDVDt=(dA_ip1half*dGrid[bin.nU][i+1][j][k]
              -dA_im1half*dGrid[bin.nU][i][j][k])/dRSq_i;
            if(dDVDt<-1.0*dDVDtThreshold){//being compressed
              dDVDt_mthreshold=dDVDt+dDVDtThreshold;
              dQ0=dASq*dGrid[bin.nD][i][j][k]*dDVDt_mthreshold*dDVDt_mthreshold;
            }
            else{
              dQ0=0.0;
            }
          }
          else if(nNumDims>=2){
            if(j==0){
              dDVDt=(dA_jp1half*dGrid[bin.nV][i][j][k]
                -dA_jm1half*dGrid[bin.nV][i][nSizeY-1][k])/dA_j;
            }
            else{
              dDVDt=(dA_jp1half*dGrid[bin.nV][i][j][k]
                -dA_jm1half*dGrid[bin.nV][i][j-1][k])/dA_j;
            }
            if(dDVDt<-1.0*dDVDtThreshold){//being compressed
              dDVDt_mthreshold=dDVDt+dDVDtThreshold;
              dQ1=dASq*dGrid[bin.nD][i][j][k]
                *dDVDt_mthreshold*dDVDt_mthreshold;
            }
            else{
              dQ1=0.0;
            }
          }
          else if(nNumDims==3){
            if(k==0){
              dDVDt=(dGrid[bin.nW][i][j][k]-dGrid[bin.nW][i][j][nSizeZ-1]);
            }
            else{
              dDVDt=(dGrid[bin.nW][i][j][k]-dGrid[bin.nW][i][j][k-1]);
            }
            if(dDVDt<-1.0*dDVDtThreshold){//being compressed
              dDVDt_mthreshold=dDVDt+dDVDtThreshold;
              dQ2=dASq*dGrid[bin.nD][i][j][k]*dDVDt_mthreshold*dDVDt_mthreshold;
            }
            else{
              dQ2=0.0;
            }
          }
          dQ=dQ0+dQ1+dQ2;
          
          if(dP>dMaxP){
            dMaxP=dP;
            nMaxJIndex[bin.nP][i]=j;
            nMaxKIndex[bin.nP][i]=k;
          }
          if(dP<dMinP){
            dMinP=dP;
            nMinJIndex[bin.nP][i]=j;
            nMinKIndex[bin.nP][i]=k;
          }
          dSumP+=dP;
          if(dQ>dMaxQ){
            dMaxQ=dQ;
            nMaxJIndex[bin.nQ][i]=j;
            nMaxKIndex[bin.nQ][i]=k;
          }
          if(dQ<dMinQ){
            dMinQ=dQ;
            nMinJIndex[bin.nQ][i]=j;
            nMinKIndex[bin.nQ][i]=k;
          }
          dSumQ+=dQ;
          if(dC>dMaxC){
            dMaxC=dC;
            nMaxJIndex[bin.nC][i]=j;
            nMaxKIndex[bin.nC][i]=k;
          }
          if(dC<dMinC){
            dMinC=dC;
            nMinJIndex[bin.nC][i]=j;
            nMinKIndex[bin.nC][i]=k;
          }
          dSumC+=dC;
          
          nCount++;
        }
      }
      dU_i=(dGrid[bin.nU0][i+1][0][0]+dGrid[bin.nU0][i][0][0])*0.5;
      if(nNumDims==3){
        dAve[bin.nD][i]=bin.dCalRhoAve3D(dGrid,i,nStartY,nEndY,nStartZ,nEndZ);
      }
      else if(nNumDims==2){
        dAve[bin.nD][i]=bin.dCalRhoAve2D(dGrid,i,nStartY,nEndY,nStartZ,nEndZ);
      }
      else{
        dAve[bin.nD][i]=dGrid[bin.nD][i][0][0];
      }
      dAve[bin.nKEP][i]=0.5*dGrid[bin.nDM][i][0][0]*dU_i*dU_i;
      dMax[bin.nP][i]=dMaxP;
      dMin[bin.nP][i]=dMinP;
      dAve[bin.nP][i]=dSumP/double(nCount);
      dMax[bin.nQ][i]=dMaxQ;
      dMin[bin.nQ][i]=dMinQ;
      dAve[bin.nQ][i]=dSumQ/double(nCount);
      dMax[bin.nC][i]=dMaxC;
      dMin[bin.nC][i]=dMinC;
      dAve[bin.nC][i]=dSumC/double(nCount);
      nMaxJIndex[bin.nKEP][i]=0;
      nMaxKIndex[bin.nKEP][i]=0;
      nMinJIndex[bin.nKEP][i]=0;
      nMinKIndex[bin.nKEP][i]=0;
    }
  }
  
  //open output file
  std::string sFileNameOut=bin.sFileName+"_pro.txt";
  std::ofstream ofFile;
  ofFile.open(sFileNameOut.c_str());
  if(!ofFile.good()){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": output file \""
      <<sFileNameOut<<" didn't open properly\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //set double output precision
  ofFile.precision(bin.nPrecisionAscii);
  if(bin.bScientific){
    ofFile.unsetf(std::ios::fixed);
    ofFile.setf(std::ios::scientific);
  }
  else{
    ofFile.unsetf(std::ios::scientific);
    ofFile.setf(std::ios::fixed);
  }
  
  //write out header to profile file
  int nWidthOutputField=25;
  int nWidthIntOutputField=12;
  ofFile<<"time= "<<dTime<<" [s]"<<std::endl;
  ofFile<<std::setw(nWidthIntOutputField)<<"Zone_#(1)"
    <<std::setw(nWidthOutputField)<<"M_r_im1half[g](2)"
    <<std::setw(nWidthOutputField)<<"DM_r[g](3)"
    <<std::setw(nWidthOutputField)<<"ErrorDM_r(4)"
    <<std::setw(nWidthOutputField)<<"R_im1half[cm](5)"
    <<std::setw(nWidthOutputField)<<"R_ip1half[cm](6)"
    <<std::setw(nWidthOutputField)<<"<D>[g/cm^3](7)"
    <<std::setw(nWidthOutputField)<<"D_max[g/cm^3](8)"
    <<std::setw(nWidthIntOutputField)<<"D_max_j(9)"
    <<std::setw(nWidthIntOutputField)<<"D_max_k(10)"
    <<std::setw(nWidthOutputField)<<"D_min[g/cm^3](11)"
    <<std::setw(nWidthIntOutputField)<<"D_min_j(12)"
    <<std::setw(nWidthIntOutputField)<<"D_min_k(13)"
    <<std::setw(nWidthOutputField)<<"U_ave_im1half[cm/s](14)"
    <<std::setw(nWidthOutputField)<<"U_max_im1half[cm/s](15)"
    <<std::setw(nWidthIntOutputField)<<"U_max_j(16)"
    <<std::setw(nWidthIntOutputField)<<"U_max_k(17)"
    <<std::setw(nWidthOutputField)<<"U_min_im1half[cm/s](18)"
    <<std::setw(nWidthIntOutputField)<<"U_min_j(19)"
    <<std::setw(nWidthIntOutputField)<<"U_min_k(20)"
    <<std::setw(nWidthOutputField)<<"U0[cm/s](21)"
    <<std::setw(nWidthOutputField)<<"V_ave[cm/s](22)"
    <<std::setw(nWidthOutputField)<<"V_max[cm/s](23)"
    <<std::setw(nWidthIntOutputField)<<"V_max_j(24)"
    <<std::setw(nWidthIntOutputField)<<"V_max_k(25)"
    <<std::setw(nWidthOutputField)<<"V_min[cm/s](26)"
    <<std::setw(nWidthIntOutputField)<<"V_min_j(27)"
    <<std::setw(nWidthIntOutputField)<<"V_min_k(28)"
    <<std::setw(nWidthOutputField)<<"W_ave[cm/s](29)"
    <<std::setw(nWidthOutputField)<<"W_max[cm/s](30)"
    <<std::setw(nWidthIntOutputField)<<"W_max_j(31)"
    <<std::setw(nWidthIntOutputField)<<"W_max_k(32)"
    <<std::setw(nWidthOutputField)<<"W_min[cm/s](33)"
    <<std::setw(nWidthIntOutputField)<<"W_min_j(34)"
    <<std::setw(nWidthIntOutputField)<<"W_min_k(35)"
    <<std::setw(nWidthOutputField)<<"Q[dynes/cm^2](36)"
    <<std::setw(nWidthOutputField)<<"Q_max[dynes/cm^2](37)"
    <<std::setw(nWidthIntOutputField)<<"Q_max_j(38)"
    <<std::setw(nWidthIntOutputField)<<"Q_max_k(39)"
    <<std::setw(nWidthOutputField)<<"Q_min[dynes/cm^2](40)"
    <<std::setw(nWidthIntOutputField)<<"Q_min_j(41)"
    <<std::setw(nWidthIntOutputField)<<"Q_min_k(42)"
    <<std::setw(nWidthOutputField)<<"E_ave[erg/g](43)"
    <<std::setw(nWidthOutputField)<<"E_max[erg/g](44)"
    <<std::setw(nWidthIntOutputField)<<"E_max_j(45)"
    <<std::setw(nWidthIntOutputField)<<"E_max_k(46)"
    <<std::setw(nWidthOutputField)<<"E_min[erg/g](47)"
    <<std::setw(nWidthIntOutputField)<<"E_min_j(48)"
    <<std::setw(nWidthIntOutputField)<<"E_min_k(49)"
    <<std::setw(nWidthOutputField)<<"T_ave[K](50)"
    <<std::setw(nWidthOutputField)<<"T_max[K](51)"
    <<std::setw(nWidthIntOutputField)<<"T_max_j(52)"
    <<std::setw(nWidthIntOutputField)<<"T_max_k(53)"
    <<std::setw(nWidthOutputField)<<"T_min[K](54)"
    <<std::setw(nWidthIntOutputField)<<"T_min_j(55)"
    <<std::setw(nWidthIntOutputField)<<"T_min_k(56)"
    <<std::setw(nWidthOutputField)<<"Kap_ave[cm^2/g](57)"
    <<std::setw(nWidthOutputField)<<"Kap_min[cm^2/g](58)"
    <<std::setw(nWidthOutputField)<<"Kap_max[cm^2/g](59)"
    <<std::setw(nWidthOutputField)<<"L_rd_im1half[L_sun](60)"
    <<std::setw(nWidthOutputField)<<"L_cv_im1half[L_sun](61)"
    <<std::setw(nWidthOutputField)<<"F_cv_ave[L_sun/cm^2](62)"
    <<std::setw(nWidthOutputField)<<"F_cv_max[L_sun/cm^2](63)"
    <<std::setw(nWidthOutputField)<<"F_cv_min[L_sun/cm^2](64)"
    <<std::setw(nWidthOutputField)<<"KEP[ergs](65)"
    <<std::setw(nWidthOutputField)<<"KETot[ergs](66)"
    <<std::setw(nWidthOutputField)<<"P_ave[dynes/cm^2](67)"
    <<std::setw(nWidthOutputField)<<"P_min[dynes/cm^2](68)"
    <<std::setw(nWidthOutputField)<<"P_max[dynes/cm^2](69)"
    <<std::setw(nWidthOutputField)<<"Gam_ave[na](70)"
    <<std::setw(nWidthOutputField)<<"Gam_min[na](71)"
    <<std::setw(nWidthOutputField)<<"Gam_max[na](72)"
    <<std::setw(nWidthOutputField)<<"C_ave[cm/s](73)"
    <<std::setw(nWidthOutputField)<<"C_max[cm/s](74)"
    <<std::setw(nWidthIntOutputField)<<"C_max_j(75)"
    <<std::setw(nWidthIntOutputField)<<"C_max_k(76)"
    <<std::setw(nWidthOutputField)<<"C_min[cm/s](77)"
    <<std::setw(nWidthIntOutputField)<<"C_min_j(78)"
    <<std::setw(nWidthIntOutputField)<<"C_min_k(79)"
    <<std::setw(nWidthOutputField)<<"UpFillFac(80)";
  if(bin.bExtraInfoInProfile){
    ofFile<<std::setw(nWidthOutputField)<<"DlnPDlnT(81)"
      <<std::setw(nWidthOutputField)<<"DlnPDlnRho(82)"
      <<std::setw(nWidthOutputField)<<"DEDT(83)";
  }
  ofFile<<std::endl;
  
  //write out profile
  double dErrorDM_r;
  for(int i=0;i<nSizeGlobe[0]+2*nNumGhostCells;i++){
    
    //calculate mass error
    dErrorDM_r=(4.0/3.0*bin.dPi*dAve[bin.nD][i]*(pow(dGrid[bin.nR][i+1][0][0],3.0)
      -pow(dGrid[bin.nR][i][0][0],3.0))-dGrid[bin.nDM][i][0][0])/dGrid[bin.nDM][i][0][0];
    double dDlnPDlnT;
    double dDlnPDlnRho;
    double dDEDT;
    eosTable->getDlnPDlnTDlnPDlnPDEDT(dAve[bin.nT][i],dAve[bin.nD][i],dDlnPDlnT,dDlnPDlnRho,dDEDT);
    
    nMaxJIndex[bin.nP][i]=0;
    nMaxKIndex[bin.nP][i]=0;
    nMinJIndex[bin.nP][i]=0;
    nMinKIndex[bin.nP][i]=0;
    
    ofFile<<std::setw(nWidthIntOutputField)<<i//1
      <<std::setw(nWidthOutputField)<<dGrid[bin.nM][i][0][0]//2
      <<std::setw(nWidthOutputField)<<dGrid[bin.nDM][i][0][0]//3
      <<std::setw(nWidthOutputField)<<dErrorDM_r//4
      <<std::setw(nWidthOutputField)<<dGrid[bin.nR][i][0][0]//5
      <<std::setw(nWidthOutputField)<<dGrid[bin.nR][i+1][0][0]//6
      <<std::setw(nWidthOutputField)<<dAve[bin.nD][i]//7
      <<std::setw(nWidthOutputField)<<dMax[bin.nD][i]//8
      <<std::setw(nWidthIntOutputField)<<nMaxJIndex[bin.nD][i]//9
      <<std::setw(nWidthIntOutputField)<<nMaxKIndex[bin.nD][i]//10
      <<std::setw(nWidthOutputField)<<dMin[bin.nD][i]//11
      <<std::setw(nWidthIntOutputField)<<nMinJIndex[bin.nD][i]//12
      <<std::setw(nWidthIntOutputField)<<nMinKIndex[bin.nD][i]//13
      <<std::setw(nWidthOutputField)<<dAve[bin.nU][i]//14
      <<std::setw(nWidthOutputField)<<dMax[bin.nU][i]//15
      <<std::setw(nWidthIntOutputField)<<nMaxJIndex[bin.nU][i]//16
      <<std::setw(nWidthIntOutputField)<<nMaxKIndex[bin.nU][i]//17
      <<std::setw(nWidthOutputField)<<dMin[bin.nU][i]//18
      <<std::setw(nWidthIntOutputField)<<nMinJIndex[bin.nU][i]//19
      <<std::setw(nWidthIntOutputField)<<nMinKIndex[bin.nU][i]//20
      <<std::setw(nWidthOutputField)<<dGrid[bin.nU0][i][0][0];//21
    if(nNumDims>1){
      ofFile
        <<std::setw(nWidthOutputField)<<dAve[bin.nV][i]//22
        <<std::setw(nWidthOutputField)<<dMax[bin.nV][i]//23
        <<std::setw(nWidthIntOutputField)<<nMaxJIndex[bin.nV][i]//24
        <<std::setw(nWidthIntOutputField)<<nMaxKIndex[bin.nV][i]//25
        <<std::setw(nWidthOutputField)<<dMin[bin.nV][i]//26
        <<std::setw(nWidthIntOutputField)<<nMinJIndex[bin.nV][i]//27
        <<std::setw(nWidthIntOutputField)<<nMinKIndex[bin.nV][i];//28
    }
    else{
      ofFile
        <<std::setw(nWidthOutputField)<<"-"//22
        <<std::setw(nWidthOutputField)<<"-"//23
        <<std::setw(nWidthIntOutputField)<<"-"//24
        <<std::setw(nWidthIntOutputField)<<"-"//25
        <<std::setw(nWidthOutputField)<<"-"//26
        <<std::setw(nWidthIntOutputField)<<"-"//27
        <<std::setw(nWidthIntOutputField)<<"-";//28
    }
    if(nNumDims>2){
      ofFile
        <<std::setw(nWidthOutputField)<<dAve[bin.nW][i]//29
        <<std::setw(nWidthOutputField)<<dMax[bin.nW][i]//30
        <<std::setw(nWidthIntOutputField)<<nMaxJIndex[bin.nW][i]//31
        <<std::setw(nWidthIntOutputField)<<nMaxKIndex[bin.nW][i]//32
        <<std::setw(nWidthOutputField)<<dMin[bin.nW][i]//33
        <<std::setw(nWidthIntOutputField)<<nMinJIndex[bin.nW][i]//34
        <<std::setw(nWidthIntOutputField)<<nMinKIndex[bin.nW][i];//35
    }
    else{
      ofFile
        <<std::setw(nWidthOutputField)<<"-"//29
        <<std::setw(nWidthOutputField)<<"-"//30
        <<std::setw(nWidthIntOutputField)<<"-"//31
        <<std::setw(nWidthIntOutputField)<<"-"//32
        <<std::setw(nWidthOutputField)<<"-"//33
        <<std::setw(nWidthIntOutputField)<<"-"//34
        <<std::setw(nWidthIntOutputField)<<"-";//35
    }
    ofFile
      <<std::setw(nWidthOutputField)<<dAve[bin.nQ][i]//36
      <<std::setw(nWidthOutputField)<<dMax[bin.nQ][i]//37
      <<std::setw(nWidthIntOutputField)<<nMaxJIndex[bin.nQ][i]//38
      <<std::setw(nWidthIntOutputField)<<nMaxKIndex[bin.nQ][i]//39
      <<std::setw(nWidthOutputField)<<dMin[bin.nQ][i]//40
      <<std::setw(nWidthIntOutputField)<<nMinJIndex[bin.nQ][i]//41
      <<std::setw(nWidthIntOutputField)<<nMinKIndex[bin.nQ][i];//42
    ofFile
      <<std::setw(nWidthOutputField)<<dAve[bin.nE][i]//43
      <<std::setw(nWidthOutputField)<<dMax[bin.nE][i]//44
      <<std::setw(nWidthIntOutputField)<<nMaxJIndex[bin.nE][i]//45
      <<std::setw(nWidthIntOutputField)<<nMaxKIndex[bin.nE][i]//46
      <<std::setw(nWidthOutputField)<<dMin[bin.nE][i]//47
      <<std::setw(nWidthIntOutputField)<<nMinJIndex[bin.nE][i]//48
      <<std::setw(nWidthIntOutputField)<<nMinKIndex[bin.nE][i];//49
    if(nGammaLaw==0){
      ofFile
        <<std::setw(nWidthOutputField)<<"-"//50
        <<std::setw(nWidthOutputField)<<"-"//51
        <<std::setw(nWidthIntOutputField)<<"-"//52
        <<std::setw(nWidthIntOutputField)<<"-"//53
        <<std::setw(nWidthOutputField)<<"-"//54
        <<std::setw(nWidthIntOutputField)<<"-"//55
        <<std::setw(nWidthIntOutputField)<<"-"//56
        <<std::setw(nWidthOutputField)<<"-"//57
        <<std::setw(nWidthOutputField)<<"-"//58
        <<std::setw(nWidthOutputField)<<"-"//59
        <<std::setw(nWidthOutputField)<<"-"//60
        <<std::setw(nWidthOutputField)<<"-";//61
    }
    else{
      ofFile
        <<std::setw(nWidthOutputField)<<dAve[bin.nT][i]//50
        <<std::setw(nWidthOutputField)<<dMax[bin.nT][i]//51
        <<std::setw(nWidthIntOutputField)<<nMaxJIndex[bin.nT][i]//52
        <<std::setw(nWidthIntOutputField)<<nMaxKIndex[bin.nT][i]//53
        <<std::setw(nWidthOutputField)<<dMin[bin.nT][i]//54
        <<std::setw(nWidthIntOutputField)<<nMinJIndex[bin.nT][i]//55
        <<std::setw(nWidthIntOutputField)<<nMinKIndex[bin.nT][i]//56
        <<std::setw(nWidthOutputField)<<dAve[bin.nKappa][i]//57
        <<std::setw(nWidthOutputField)<<dMax[bin.nKappa][i]//58
        <<std::setw(nWidthOutputField)<<dMin[bin.nKappa][i];//59
      if(i!=0){
        ofFile
          <<std::setw(nWidthOutputField)<<dAve[bin.nL_rad][i-1]/bin.dLSun//60
          <<std::setw(nWidthOutputField)<<dAve[bin.nL_con][i-1]/bin.dLSun//61
          <<std::setw(nWidthOutputField)<<dAve[bin.nF_con][i-1]/bin.dLSun//62
          <<std::setw(nWidthOutputField)<<dMax[bin.nF_con][i-1]/bin.dLSun//63
          <<std::setw(nWidthOutputField)<<dMin[bin.nF_con][i-1]/bin.dLSun;//64
      }
      else{//luminosity not defined at inner interface
        ofFile
          <<std::setw(nWidthOutputField)<<"-"//60
          <<std::setw(nWidthOutputField)<<"-"//61
          <<std::setw(nWidthOutputField)<<"-"//62
          <<std::setw(nWidthOutputField)<<"-"//63
          <<std::setw(nWidthOutputField)<<"-";//64
      }
    }
    ofFile
      <<std::setw(nWidthOutputField)<<dAve[bin.nKEP][i]//65
      <<std::setw(nWidthOutputField)<<dAve[bin.nKETot][i]//66
      <<std::setw(nWidthOutputField)<<dAve[bin.nP][i]//67
      <<std::setw(nWidthOutputField)<<dMax[bin.nP][i]//68
      <<std::setw(nWidthOutputField)<<dMin[bin.nP][i]//69
      <<std::setw(nWidthOutputField)<<dAve[bin.nGamma][i]//70
      <<std::setw(nWidthOutputField)<<dMax[bin.nGamma][i]//71
      <<std::setw(nWidthOutputField)<<dMin[bin.nGamma][i]//72
      <<std::setw(nWidthOutputField)<<dAve[bin.nC][i]//73
      <<std::setw(nWidthOutputField)<<dMax[bin.nC][i]//74
      <<std::setw(nWidthIntOutputField)<<nMaxJIndex[bin.nC][i]//75
      <<std::setw(nWidthIntOutputField)<<nMaxKIndex[bin.nC][i]//76
      <<std::setw(nWidthOutputField)<<dMin[bin.nC][i]//77
      <<std::setw(nWidthIntOutputField)<<nMinJIndex[bin.nC][i]//78
      <<std::setw(nWidthIntOutputField)<<nMinKIndex[bin.nC][i]//79
      <<std::setw(nWidthOutputField)<<dUpFlowFillingFactor[i];//80
    if(bin.bExtraInfoInProfile){
      ofFile<<std::setw(nWidthOutputField)<<dDlnPDlnT//81
        <<std::setw(nWidthOutputField)<<dDlnPDlnRho//82
        <<std::setw(nWidthOutputField)<<dDEDT;//83
    }
    ofFile<<std::endl;
  }
  ofFile
    <<std::setw(nWidthIntOutputField)<<(nSizeGlobe[0]+2*nNumGhostCells)//1
    <<std::setw(nWidthOutputField)<<dGrid[bin.nM][nSizeGlobe[0]+2*nNumGhostCells][0][0]//2
    <<std::setw(nWidthOutputField)<<"-"//3
    <<std::setw(nWidthOutputField)<<"-"//4
    <<std::setw(nWidthOutputField)<<dGrid[bin.nR][nSizeGlobe[0]+2*nNumGhostCells][0][0]//5
    <<std::setw(nWidthOutputField)<<"-"//6
    <<std::setw(nWidthOutputField)<<"-"//7
    <<std::setw(nWidthOutputField)<<"-"//8
    <<std::setw(nWidthIntOutputField)<<"-"//9
    <<std::setw(nWidthIntOutputField)<<"-"//10
    <<std::setw(nWidthOutputField)<<"-"//11
    <<std::setw(nWidthIntOutputField)<<"-"//12
    <<std::setw(nWidthIntOutputField)<<"-"//13
    <<std::setw(nWidthOutputField)<<dAve[bin.nU][nSizeGlobe[0]+2*nNumGhostCells]//14
    <<std::setw(nWidthOutputField)<<dMax[bin.nU][nSizeGlobe[0]+2*nNumGhostCells]//15
    <<std::setw(nWidthIntOutputField)<<nMaxJIndex[bin.nU][i]//16
    <<std::setw(nWidthIntOutputField)<<nMaxKIndex[bin.nU][i]//17
    <<std::setw(nWidthOutputField)<<dMin[bin.nU][nSizeGlobe[0]+2*nNumGhostCells]//18
    <<std::setw(nWidthIntOutputField)<<nMinJIndex[bin.nU][i]//19
    <<std::setw(nWidthIntOutputField)<<nMinKIndex[bin.nU][i]//20
    <<std::setw(nWidthOutputField)<<dAve[bin.nU0][nSizeGlobe[0]+2*nNumGhostCells]//21
    <<std::setw(nWidthOutputField)<<"-"//22
    <<std::setw(nWidthOutputField)<<"-"//23
    <<std::setw(nWidthIntOutputField)<<"-"//24
    <<std::setw(nWidthIntOutputField)<<"-"//25
    <<std::setw(nWidthOutputField)<<"-"//26
    <<std::setw(nWidthIntOutputField)<<"-"//27
    <<std::setw(nWidthIntOutputField)<<"-"//28
    <<std::setw(nWidthOutputField)<<"-"//29
    <<std::setw(nWidthOutputField)<<"-"//30
    <<std::setw(nWidthIntOutputField)<<"-"//31
    <<std::setw(nWidthIntOutputField)<<"-"//32
    <<std::setw(nWidthOutputField)<<"-"//33
    <<std::setw(nWidthIntOutputField)<<"-"//34
    <<std::setw(nWidthIntOutputField)<<"-"//35
    <<std::setw(nWidthOutputField)<<"-"//36
    <<std::setw(nWidthOutputField)<<"-"//37
    <<std::setw(nWidthIntOutputField)<<"-"//38
    <<std::setw(nWidthIntOutputField)<<"-"//39
    <<std::setw(nWidthOutputField)<<"-"//40
    <<std::setw(nWidthIntOutputField)<<"-"//41
    <<std::setw(nWidthIntOutputField)<<"-"//42
    <<std::setw(nWidthOutputField)<<"-"//43
    <<std::setw(nWidthOutputField)<<"-"//44
    <<std::setw(nWidthIntOutputField)<<"-"//45
    <<std::setw(nWidthIntOutputField)<<"-"//46
    <<std::setw(nWidthOutputField)<<"-"//47
    <<std::setw(nWidthIntOutputField)<<"-"//48
    <<std::setw(nWidthIntOutputField)<<"-"//49
    <<std::setw(nWidthOutputField)<<"-"//50
    <<std::setw(nWidthOutputField)<<"-"//51
    <<std::setw(nWidthIntOutputField)<<"-"//52
    <<std::setw(nWidthIntOutputField)<<"-"//53
    <<std::setw(nWidthOutputField)<<"-"//54
    <<std::setw(nWidthIntOutputField)<<"-"//55
    <<std::setw(nWidthIntOutputField)<<"-"//56
    <<std::setw(nWidthOutputField)<<"-"//57
    <<std::setw(nWidthOutputField)<<"-"//58
    <<std::setw(nWidthOutputField)<<"-"//59
    <<std::setw(nWidthOutputField)<<dAve[bin.nL_rad][nSizeGlobe[0]+2*nNumGhostCells-1]/bin.dLSun//60
    <<std::setw(nWidthOutputField)<<dAve[bin.nL_con][nSizeGlobe[0]+2*nNumGhostCells-1]/bin.dLSun//61
    <<std::setw(nWidthOutputField)<<dAve[bin.nF_con][nSizeGlobe[0]+2*nNumGhostCells-1]/bin.dLSun//62
    <<std::setw(nWidthOutputField)<<dMax[bin.nF_con][nSizeGlobe[0]+2*nNumGhostCells-1]/bin.dLSun//63
    <<std::setw(nWidthOutputField)<<dMin[bin.nF_con][nSizeGlobe[0]+2*nNumGhostCells-1]/bin.dLSun//64
    <<std::setw(nWidthOutputField)<<"-"//65
    <<std::setw(nWidthOutputField)<<"-"//66
    <<std::setw(nWidthOutputField)<<"-"//67
    <<std::setw(nWidthOutputField)<<"-"//68
    <<std::setw(nWidthOutputField)<<"-"//69
    <<std::setw(nWidthOutputField)<<"-"//70
    <<std::setw(nWidthOutputField)<<"-"//71
    <<std::setw(nWidthOutputField)<<"-"//72
    <<std::setw(nWidthOutputField)<<"-"//73
    <<std::setw(nWidthOutputField)<<"-"//74
    <<std::setw(nWidthIntOutputField)<<"-"//75
    <<std::setw(nWidthIntOutputField)<<"-"//76
    <<std::setw(nWidthOutputField)<<"-"//77
    <<std::setw(nWidthIntOutputField)<<"-"//78
    <<std::setw(nWidthIntOutputField)<<"-"//79
    <<std::setw(nWidthOutputField)<<dUpFlowFillingFactor[nSizeGlobe[0]+2*nNumGhostCells];//80
  if(bin.bExtraInfoInProfile){
    ofFile<<std::setw(nWidthOutputField)<<"-"//81
      <<std::setw(nWidthOutputField)<<"-"//82
      <<std::setw(nWidthOutputField)<<"-";//83
  }
  ofFile<<std::endl;
  
  ofFile.close();
  
  //delete grid
  for(int n=0;n<nNumVars;n++){
    nGhostCellsX=1;
    if(nVarInfo[n][0]==-1){
      nGhostCellsX=0;
    }
    nGhostCellsY=1;
    if(nVarInfo[n][1]==-1){
      nGhostCellsY=0;
    }
    nGhostCellsZ=1;
    if(nVarInfo[n][2]==-1){
      nGhostCellsZ=0;
    }
      
    //delete 1D part of the grid
    nSizeX1=nGhostCellsX*(nNum1DZones+nNumGhostCells);
    if (nVarInfo[n][0]==1&&nPeriodic[0]==0){
      nSizeX1=nGhostCellsX*(nNum1DZones+1+nNumGhostCells);
    }
    nSizeY=1;
    nSizeZ=1;
    for(i=0;i<nSizeX1;i++){
      for(j=0;j<nSizeY;j++){
        delete [] dGrid[n][i][j];
      }
      delete [] dGrid[n][i];
    }
    
    //read in the rest of the grid
    nSizeX2=nSize[n][0]+nGhostCellsX*2*nNumGhostCells;
    nSizeY=nSize[n][1]+nGhostCellsY*2*nNumGhostCells;//assume y and z are always periodic
    nSizeZ=nSize[n][2]+nGhostCellsZ*2*nNumGhostCells;
    for(i=nSizeX1;i<nSizeX2;i++){
      for(j=0;j<nSizeY;j++){
        delete [] dGrid[n][i][j];
      }
      delete [] dGrid[n][i];
    }
    delete [] dGrid[n];
    
    delete [] dMax[n];
    delete [] dMin[n];
    delete [] dAve[n];
    delete [] nMaxJIndex[n];
    delete [] nMaxKIndex[n];
    delete [] nMinJIndex[n];
    delete [] nMinKIndex[n];
  }
  delete [] dMax;
  delete [] dMin;
  delete [] dAve;
  delete [] nMaxJIndex;
  delete [] nMaxKIndex;
  delete [] nMinJIndex;
  delete [] nMinKIndex;
  delete [] dGrid;
  
  delete [] dUpFlowFillingFactor;
  for(int n=0;n<nNumVars;n++){
    delete [] nSize[n];
    delete [] nVarInfo[n];
  }
  delete [] nSize;
  delete [] nVarInfo;
}
