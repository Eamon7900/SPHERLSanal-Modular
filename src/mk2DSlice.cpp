#include <math.h>
#include <iomanip>
#include <unistd.h>
#include <string.h>

#include "binfile.h"
#include "datafile.h"
#include "eos.h"
int main(int argc, char** argv){
}
void make2DSlice(DataFile *bin, int nPlane,int nPlaneIndex){//updated
  
  int nWidthOutputField=25;
  
  //open input file
  if(bin->sFileName.size()==0){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": no input file specified\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  std::ifstream ifFile;
  ifFile.open(bin->sFileName.c_str(),std::ios::binary);
  if(!ifFile.good()){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": input file \""
      <<bin->sFileName<<"\" didn't open properly\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //check that it is a binary file
  char cTemp;
  ifFile.read((char*)(&cTemp),sizeof(char));
  if(cTemp!='b'){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": input file \""
      <<bin->sFileName<<"\" isn't a binary file.\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //check that it is the correct version
  int nTemp;
  ifFile.read((char*)(&nTemp),sizeof(int));
  if(nTemp!=bin->nDumpFileVersion){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": inpput file \""
      <<bin->sFileName<<"\" version \""<<nTemp
      <<"\" isn't the supported version \""<<bin->nDumpFileVersion<<"\".\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //read in binary input file and write out in averaged profiles
  
  //read in time
  double dTime;
  ifFile.read((char*)(&dTime),sizeof(double));
  
  //read in the time step index
  int nTimeStepIndex;
  ifFile.read((char*)(&nTimeStepIndex),sizeof(int));
  
  //read in timestep
  double dTimeStep1;
  ifFile.read((char*)(&dTimeStep1),sizeof(double));
  
  //read in timestep
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
  eos *eosTable;
  if(nGammaLaw==0){
    ifFile.read((char*)(&dGamma),sizeof(double));
  }
  else{
    char *cBuffer=new char[nGammaLaw+1];
    ifFile.read(cBuffer,nGammaLaw*sizeof(char));
    cBuffer[nGammaLaw]='\0';
    sEOSTable=cBuffer;
    delete [] cBuffer;
    if(bin->sEOSFile!=""){//overwrite sEOSTable if sEOSFile is set
      sEOSTable=bin->sEOSFile;
    }
    eosTable = new eos(sEOSTable);
    eosTable->readBin();
  }
  
  //read in artificial viscosity
  double dA;
  ifFile.read((char*)(&dA),sizeof(double));
  double dASq=dA*dA;
  
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
  for(int n=0;n<nNumVars;n++){
    nSize[n]=new int[3];
    nVarInfo[n]=new int[3];
    ifFile.read((char*)(nVarInfo[n]),(4)*sizeof(int));
    for(int l=0;l<3;l++){
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
  
  //set global grid sizes
  int nSize0=nSizeGlobe[0]+1+2*nNumGhostCells;
  if(nSizeGlobe[0]==1){
    nSize0=nSizeGlobe[0]+1;//don't need ghost cells if grid not defined in direction l
  }
  int nSize1=nSizeGlobe[1]+1+2*nNumGhostCells;
  if(nSizeGlobe[1]==1){
    nSize1=nSizeGlobe[1];//don't need ghost cells if grid not defined in direction l
  }
  int nSize2=nSizeGlobe[2]+1+2*nNumGhostCells;
  if(nSizeGlobe[2]==1){
    nSize2=nSizeGlobe[2]+1;//don't need ghost cells if grid not defined in direction l
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
  
  //read in grid
  double ****dGrid=new double***[nNumVars];
  for(int n=0;n<nNumVars;n++){
    
    int nGhostCellsX=1;
    if(nVarInfo[n][0]==-1){
      nGhostCellsX=0;
    }
    int nGhostCellsY=1;
    if(nVarInfo[n][1]==-1){
      nGhostCellsY=0;
    }
    int nGhostCellsZ=1;
    if(nVarInfo[n][2]==-1){
      nGhostCellsZ=0;
    }
    
    dGrid[n]=new double**[nSize[n][0]+nGhostCellsX*2*nNumGhostCells];
    
    //read in 1D part of the grid
    int nSizeX1=nGhostCellsX*(nNum1DZones+nNumGhostCells);//may be need to +1 if only one proc and variable in interface centered
    if (nVarInfo[n][0]==1&&nPeriodic[0]==0){
      nSizeX1=nGhostCellsX*(nNum1DZones+1+nNumGhostCells);
    }
    int nSizeY=1;
    int nSizeZ=1;
    for(int i=0;i<nSizeX1;i++){
      dGrid[n][i]=new double*[nSizeY];
      for(int j=0;j<nSizeY;j++){
        dGrid[n][i][j]=new double[nSizeZ];
        ifFile.read((char*)(dGrid[n][i][j]),nSizeZ*sizeof(double));
      }
    }
    
    //read in the rest of the grid
    int nSizeX2=nSize[n][0]+nGhostCellsX*2*nNumGhostCells;
    nSizeY=nSize[n][1]+nGhostCellsY*2*nNumGhostCells;
    nSizeZ=nSize[n][2]+nGhostCellsZ*2*nNumGhostCells;
    for(int i=nSizeX1;i<nSizeX2;i++){
      dGrid[n][i]=new double*[nSizeY];
      for(int j=0;j<nSizeY;j++){
        if(ifFile.eof()){
          std::stringstream ssTemp;
          ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
            <<": reached end of file sooner than expected\n";
          throw exception2(ssTemp.str(),INPUT);
        }
        dGrid[n][i][j]=new double[nSizeZ];
        ifFile.read((char*)(dGrid[n][i][j]),(nSizeZ)*sizeof(double));
      }
    }
  }
  ifFile.close();
  
  //set variable indices
  if(nGammaLaw==0){//using gamma law gas
    if(nNumDims==1){
      bin->nM=0;
      bin->nTheta=-1;
      bin->nPhi=-1;
      bin->nDM=1;
      bin->nR=2;
      bin->nD=3;
      bin->nU=4;
      bin->nU0=5;
      bin->nV=-1;
      bin->nW=-1;
      bin->nE=6;
      bin->nT=-1;
    }
    else if(nNumDims==2){
      bin->nM=0;
      bin->nTheta=1;
      bin->nPhi=-1;
      bin->nDM=2;
      bin->nR=3;
      bin->nD=4;
      bin->nU=5;
      bin->nU0=6;
      bin->nV=7;
      bin->nW=-1;
      bin->nE=8;
      bin->nT=-1;
    }
    else if(nNumDims==3){
      bin->nM=0;
      bin->nTheta=1;
      bin->nPhi=2;
      bin->nDM=3;
      bin->nR=4;
      bin->nD=5;
      bin->nU=6;
      bin->nU0=7;
      bin->nV=8;
      bin->nW=9;
      bin->nE=10;
      bin->nT=-1;
    }
  }
  else{//using a tabulated equaiton of state
    if(nNumDims==1){
      bin->nM=0;
      bin->nTheta=-1;
      bin->nPhi=-1;
      bin->nDM=1;
      bin->nR=2;
      bin->nD=3;
      bin->nU=4;
      bin->nU0=5;
      bin->nV=-1;
      bin->nW=-1;
      bin->nE=-1;
      bin->nT=6;
    }
    else if(nNumDims==2){
      bin->nM=0;
      bin->nTheta=1;
      bin->nPhi=-1;
      bin->nDM=2;
      bin->nR=3;
      bin->nD=4;
      bin->nU=5;
      bin->nU0=6;
      bin->nV=7;
      bin->nW=-1;
      bin->nE=-1;
      bin->nT=8;
    }
    else if(nNumDims==3){
      bin->nM=0;
      bin->nTheta=1;
      bin->nPhi=2;
      bin->nDM=3;
      bin->nR=4;
      bin->nD=5;
      bin->nU=6;
      bin->nU0=7;
      bin->nV=8;
      bin->nW=9;
      bin->nE=-1;
      bin->nT=10;
    }
  }
  
  //check that slice index is within the model
  if (nPlane==0){//r-theta
    if(nPlaneIndex>nSize[bin->nD][2]+2*nNumGhostCells-1){
      std::stringstream ssTemp;
      ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": r-theta plane index "<<nPlaneIndex
        <<" is larger than size of input model in phi-direction of "
        <<(nSize[bin->nD][2]+2*nNumGhostCells-1)<<std::endl;
      throw exception2(ssTemp.str(),INPUT);
    }
    else if(nPlaneIndex>0 && nNumDims==2){
      std::stringstream ssTemp;
      ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": r-theta plane index "<<nPlaneIndex
        <<" is larger than size of input model in phi-direction of "
        <<(nSize[bin->nD][2]-1)<<std::endl;
      throw exception2(ssTemp.str(),INPUT);
    }
  }
  else if(nPlane==1){//theta-phi
    if(nPlaneIndex>nSize[bin->nD][0]+2*nNumGhostCells-1){
      std::stringstream ssTemp;
      ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": theta-phi plane index "<<nPlaneIndex
        <<" is larger than size of input model in r-direction of "
        <<(nSize[bin->nD][0]+2*nNumGhostCells-1)<<std::endl;
      throw exception2(ssTemp.str(),INPUT);
    }
    if(nPlaneIndex<nNum1DZones-1){
      std::stringstream ssTemp;
      ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": theta-phi plane index "<<nPlaneIndex
        <<" is inside the 1D region of input model, 3D region begins at radial zone "<<nNum1DZones
        <<std::endl;
      throw exception2(ssTemp.str(),INPUT);
    }
  }
  else if(nPlane==2){//r-phi
    if(nPlaneIndex>nSize[bin->nD][1]+2*nNumGhostCells-1){
      std::stringstream ssTemp;
      ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": r-phi plane index "<<nPlaneIndex
        <<" is larger than size of input model in theta-direction of "
        <<(nSize[bin->nD][1]+2*nNumGhostCells-1)<<std::endl;
      throw exception2(ssTemp.str(),INPUT);
    }
  }
  
  //open output file
  std::stringstream fileNameOut;
  if(nPlane==0){//r-theta plane
    if(nPlaneIndex>=0&&nPlaneIndex<nSize1){
      fileNameOut<<bin->sFileName<<"_2Dk="<<nPlaneIndex<<".txt";
    }
    else{
      std::stringstream ssTemp;
      ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": nPlaneIndex="<<nPlaneIndex
        <<" is outside theta zoning\n";
      throw exception2(ssTemp.str(),SYNTAX);
    }
  }
  else if(nPlane==1){//theta-phi plane
    if(nPlaneIndex>nNum1DZones+nNumGhostCells&&nPlaneIndex<nSize0){
      fileNameOut<<bin->sFileName<<"_2Di="<<nPlaneIndex<<".txt";
    }
    else{
      std::stringstream ssTemp;
      ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": nPlaneIndex="<<nPlaneIndex
        <<" is inside 1D radial only region or outside of radial zoning\n"
        <<" can not make theta-phi plane slice\n";
      throw exception2(ssTemp.str(),SYNTAX);
    }
  }
  else if(nPlane==2){//phi-r plane
    if(nPlaneIndex>=0&&nPlaneIndex<nSize2){
      fileNameOut<<bin->sFileName<<"_2Dj="<<nPlaneIndex<<".txt";
    }
    else{
      std::stringstream ssTemp;
      ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": nPlaneIndex="<<nPlaneIndex
        <<" is outside phi zoning\n";
      throw exception2(ssTemp.str(),SYNTAX);
    }
  }
  else{
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": unknown plane type\n";
    throw exception2(ssTemp.str(),SYNTAX);
  }
  std::ofstream ofFile;
  ofFile.open(fileNameOut.str().c_str());
  if(!ofFile.good()){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": output file \""
      <<fileNameOut.str()<<" didn't open properly\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //set double output precision
  ofFile.precision(bin->nPrecisionAscii);
  if(bin->bScientific){
    ofFile.unsetf(std::ios::fixed);
    ofFile.setf(std::ios::scientific);
  }
  else{
    ofFile.unsetf(std::ios::scientific);
    ofFile.setf(std::ios::fixed);
  }
  
  //write header to file
  if(nPlane==0){//r-theta plane
    ofFile<<"rt"<<std::endl;
  }
  else if(nPlane==1){//theta-phi plane
    ofFile<<"tp"<<std::endl;
  }
  else if(nPlane==2){//phi-r plane
    ofFile<<"rp"<<std::endl;
  }
  ofFile<<"t[s] "<<dTime<<" ";
  if(nGammaLaw==0){
    ofFile<<"0 "<<dGamma;
  }
  else{
    ofFile<<sEOSTable.size()<<" "<<sEOSTable;
  }
  ofFile<<std::endl;
  
  //write out slice
  double dP;
  double dE;
  double dKappa;
  double dQ;
  double dQ0;
  double dQ1;
  double dQ2;
  double dR_i;
  double dRSq_i;
  double dA_ip1half;
  double dA_im1half;
  double dRSq_ip1half;
  double dC;
  double dDVDtThreshold;
  double dDVDt;
  double dDVDt_mthreshold;
  double dTheta_jp1half;
  double dTheta_jm1half;
  double dTheta_j;
  double dA_jp1half;
  double dA_jm1half;
  double dA_j;
  double dCp;
  if(nPlane==0){//r-theta
    
    //write out cooridnate variables
    //i
    ofFile<<"im1half(0) ";
    for( int i=0;i<nSize[bin->nR][0]+2*nNumGhostCells;i++){
      ofFile<<i<<" ";
    }
    ofFile<<std::endl;
    
    //M_r
    ofFile<<"M_r_im1half[g](1) ";
    for( int i=0;i<nSize[bin->nR][0]+2*nNumGhostCells;i++){
      ofFile<<dGrid[bin->nM][i][0][0]<<" ";
    }
    ofFile<<std::endl;
    
    //R
    ofFile<<"R_im1half[cm](2) ";
    for( int i=0;i<nSize[bin->nR][0]+2*nNumGhostCells;i++){
      ofFile<<dGrid[bin->nR][i][0][0]<<" ";
    }
    ofFile<<std::endl;
    
    //j
    ofFile<<"jm1half(3) ";
    for( int j=0;j<nSize[bin->nD][1]+1+2*nNumGhostCells;j++){//add an extra since it is an interface
      ofFile<<j<<" ";
    }
    ofFile<<std::endl;
    
    //theta
    int nInterFaceY=0;
    if(nPeriodic[1]==0){
      nInterFaceY=1;
    }
    ofFile<<"theta_jm1half[rad](4) ";
    if(nPeriodic[1]==1){/*write out inner interface if periodic, if not periodic, it is already 
      included*/
      double dInnerTheta=dGrid[bin->nTheta][0][0][0]-(dGrid[bin->nTheta][0][nSize[bin->nD][1]+nNumGhostCells-2][0]
        -dGrid[bin->nTheta][0][nSize[bin->nD][1]+nNumGhostCells-3][0]);
      ofFile<<dInnerTheta<<" ";
     
    }
    for( int j=0;j<nSize[bin->nD][1]+nInterFaceY+2*nNumGhostCells;j++){
      ofFile<<dGrid[bin->nTheta][0][j][0]<<" ";
    }
    ofFile<<std::endl;
    
    //k
    if(nNumDims>2){
      ofFile<<"km1half(5) "<<nPlaneIndex<<" "<<nPlaneIndex+1<<std::endl;
    }
    else{
      ofFile<<"km1half(5) "<<0<<" "<<1<<std::endl;
    }
    
    //phi
    if(nNumDims>2){
      ofFile<<"phi_km1half[rad](6) "<<dGrid[bin->nPhi][0][0][nPlaneIndex]<<" ";
      if(nPlaneIndex==nSize[bin->nPhi][2]+2*nNumGhostCells-1){//if in last zone, need to do something to get outter phi
        ofFile<<dGrid[bin->nPhi][0][0][nPlaneIndex]+(dGrid[bin->nPhi][0][0][nPlaneIndex]
          -dGrid[bin->nPhi][0][0][nPlaneIndex-1])<<std::endl;
      }
      else{
        ofFile<<dGrid[bin->nPhi][0][0][nPlaneIndex+1]<<std::endl;
      }
    }
    else{
      ofFile<<"phi_km1half[rad](6) "<<0.0<<" "<<0.0<<std::endl;
    }
    
    //write out header for 2D data
    ofFile
      <<std::setw(nWidthOutputField)<<"U_im1halfjk[cm/s](7)"
      <<std::setw(nWidthOutputField)<<"U0_im1half[cm/s](8)"
      <<std::setw(nWidthOutputField)<<"V_ijm1halfk[cm/s](9)"
      <<std::setw(nWidthOutputField)<<"W_ijkm1half[cm/s](10)"
      <<std::setw(nWidthOutputField)<<"D_ijk[g/cm^3](11)"<<std::setw(nWidthOutputField)
      <<"D_rel_dif_hor_ave(12)"
      <<std::setw(nWidthOutputField)<<"E_ijk[erg/g](13)"<<std::setw(nWidthOutputField)
      <<"E_rel_dif_hor_ave(14)"
      <<std::setw(nWidthOutputField)<<"T_ijk[K](15)"<<std::setw(nWidthOutputField)
      <<"T_rel_dif_hor_ave(16)"
      <<std::setw(nWidthOutputField)<<"P_ijk[dynes/cm^2](17)"<<std::setw(nWidthOutputField)
      <<"P_rel_dif_hor_ave(18)"
      <<std::setw(nWidthOutputField)<<"Q_ijk[dynes/cm^2](19)"<<std::setw(nWidthOutputField)
      <<"Q_rel_dif_hor_ave(20)"
      <<std::setw(nWidthOutputField)<<"Kap_ij[cm^2/g](21)"<<std::setw(nWidthOutputField)
      <<"Kap_rel_dif_hor_ave(22)"
      <<std::setw(nWidthOutputField)<<"Gam_ijk[na](23)"<<std::setw(nWidthOutputField)
      <<"Gam_rel_dif_hor_ave(24)"<<std::setw(nWidthOutputField)<<"C_P[na](25)"
      <<std::endl;
    
    //copy 1D region to 2D grid
    int nSizeX1=nNum1DZones+nNumGhostCells;
    int nSizeY2=nSize[bin->nD][1]+2*nNumGhostCells;
    int nSizeZ=nSize[bin->nD][2]+2*nNumGhostCells;
    for(int i=0;i<nSizeX1;i++){
      dR_i=(dGrid[bin->nR][i+1][0][0]+dGrid[bin->nR][i][0][0])*0.5;
      dRSq_ip1half=dGrid[bin->nR][i+1][0][0]*dGrid[bin->nR][i+1][0][0];
      dRSq_i=dR_i*dR_i;
      dA_ip1half=dGrid[bin->nR][i+1][0][0]*dGrid[bin->nR][i+1][0][0];
      dA_im1half=dGrid[bin->nR][i][0][0]*dGrid[bin->nR][i][0][0];
      for(int j=0;j<nSizeY2;j++){
        
        ofFile
          <<std::setw(nWidthOutputField)<<dGrid[bin->nU][i][0][0]//0
          <<std::setw(nWidthOutputField)<<dGrid[bin->nU0][i][0][0]//1
          <<std::setw(nWidthOutputField)<<0.0;//2
        if(bin->nW==-1){
          ofFile<<std::setw(nWidthOutputField)<<"-";//3
        }
        else{
          ofFile<<std::setw(nWidthOutputField)<<0.0;//3
        }
        ofFile<<std::setw(nWidthOutputField)<<dGrid[bin->nD][i][0][0]
          <<std::setw(nWidthOutputField)<<0.0;//4,5
        
        if(nGammaLaw!=0){//set P,E,kappa,gamma, Q, and L
          
          //get P, E, kappa, and gamma
         eosTable->getPEKappaGammaCp(dGrid[bin->nT][i][0][0],dGrid[bin->nD][i][0][0],dP,dE,dKappa,dGamma,dCp);
          
          //calculate Q
          dC=sqrt(dGamma*dP/dGrid[bin->nD][i][0][0]);
          dDVDtThreshold=dAVThreshold*dC;
          dDVDt=(dA_ip1half*dGrid[bin->nU][i+1][0][0]
            -dA_im1half*dGrid[bin->nU][i][0][0])/dRSq_i;
          if(dDVDt<-1.0*dDVDtThreshold){//being compressed
            dDVDt_mthreshold=dDVDt+dDVDtThreshold;
            dQ=dASq*dGrid[bin->nD][i][0][0]*dDVDt_mthreshold*dDVDt_mthreshold;
          }
          else{
            dQ=0.0;
          }
          
          //print them out
          ofFile
            <<std::setw(nWidthOutputField)<<dE<<std::setw(nWidthOutputField)<<0.0//6,7
            <<std::setw(nWidthOutputField)<<dGrid[bin->nT][i][0][0]<<std::setw(nWidthOutputField)<<0.0//8,9
            <<std::setw(nWidthOutputField)<<dP<<std::setw(nWidthOutputField)<<0.0//10,11
            <<std::setw(nWidthOutputField)<<dQ<<std::setw(nWidthOutputField)<<0.0//12,13
            <<std::setw(nWidthOutputField)<<dKappa<<std::setw(nWidthOutputField)<<0.0//14,15
            <<std::setw(nWidthOutputField)<<dGamma<<std::setw(nWidthOutputField)<<0.0//16,17
            <<std::setw(nWidthOutputField)<<dCp<<std::setw(nWidthOutputField);//18
          
        }
        else{//set P and Q
          
          dP=dGrid[bin->nD][i][0][0]*(dGamma-1.0)*dGrid[bin->nE][i][0][0];
          
          
          //calculate Q
          dR_i=(dGrid[bin->nR][i+1][0][0]+dGrid[bin->nR][i][0][0])*0.5;
          dRSq_i=dR_i*dR_i;
          dA_ip1half=dGrid[bin->nR][i+1][0][0]*dGrid[bin->nR][i+1][0][0];
          dA_im1half=dGrid[bin->nR][i][0][0]*dGrid[bin->nR][i][0][0];
          dC=sqrt(dGamma*dP/dGrid[bin->nD][i][0][0]);
          dDVDtThreshold=dAVThreshold*dC;
          dDVDt=(dA_ip1half*dGrid[bin->nU][i+1][0][0]
            -dA_im1half*dGrid[bin->nU][i][0][0])/dRSq_i;
          if(dDVDt<-1.0*dDVDtThreshold){//being compressed
            dDVDt_mthreshold=dDVDt+dDVDtThreshold;
            dQ=dASq*dGrid[bin->nD][i][0][0]*dDVDt_mthreshold*dDVDt_mthreshold;
          }
          else{
            dQ=0.0;
          }
          
          //print them out
          ofFile
            <<std::setw(nWidthOutputField)<<dGrid[bin->nE][i][0][0]<<std::setw(nWidthOutputField)<<0.0//6,7
            <<std::setw(nWidthOutputField)<<"-"<<std::setw(nWidthOutputField)<<"-"//8,9
            <<std::setw(nWidthOutputField)<<dP<<std::setw(nWidthOutputField)<<0.0//10,11
            <<std::setw(nWidthOutputField)<<dQ<<std::setw(nWidthOutputField)<<0.0//12,13
            <<std::setw(nWidthOutputField)<<"-"<<std::setw(nWidthOutputField)<<"-"//14,15
            <<std::setw(nWidthOutputField)<<dGamma<<std::setw(nWidthOutputField)<<0.0//16,17
            <<std::setw(nWidthOutputField)<<"-";//18
        }
        ofFile<<std::endl;
      }
    ofFile
      <<std::setw(nWidthOutputField)<<"-"//0
      <<std::setw(nWidthOutputField)<<"-"//1
      <<std::setw(nWidthOutputField)<<dGrid[bin->nV][i][0][0]//2
      <<std::setw(nWidthOutputField)<<"-"//3
      <<std::setw(nWidthOutputField)<<"-"<<std::setw(nWidthOutputField)//4
      <<"-"//5
      <<std::setw(nWidthOutputField)<<"-"<<std::setw(nWidthOutputField)//6
      <<"-"//7
      <<std::setw(nWidthOutputField)<<"-"<<std::setw(nWidthOutputField)//8
      <<"-"//9
      <<std::setw(nWidthOutputField)<<"-"<<std::setw(nWidthOutputField)//10
      <<"-"//11
      <<std::setw(nWidthOutputField)<<"-"<<std::setw(nWidthOutputField)//12
      <<"-"//13
      <<std::setw(nWidthOutputField)<<"-"<<std::setw(nWidthOutputField)//14
      <<"-"//15
      <<std::setw(nWidthOutputField)<<"-"<<std::setw(nWidthOutputField)//16
      <<"-"//17
      <<std::setw(nWidthOutputField)<<"-"//18
      <<std::endl;
    }
    
    //3D region
    int nSizeX2=nSize[bin->nD][0]+2*nNumGhostCells;
    for(int i=nSizeX1;i<nSizeX2;i++){
      
      dR_i=(dGrid[bin->nR][i+1][0][0]+dGrid[bin->nR][i][0][0])*0.5;
      dRSq_ip1half=dGrid[bin->nR][i+1][0][0]*dGrid[bin->nR][i+1][0][0];
      dRSq_i=dR_i*dR_i;
      dA_ip1half=dGrid[bin->nR][i+1][0][0]*dGrid[bin->nR][i+1][0][0];
      dA_im1half=dGrid[bin->nR][i][0][0]*dGrid[bin->nR][i][0][0];
      
      //calculate horizontal average of various quantities
      double dDAve=0.0;
      double dEAve=0.0;
      double dTAve=0.0;
      double dPAve=0.0;
      double dQAve=0.0;
      double dKappaAve=0.0;
      double dGammaAve=0.0;
      int nCount=0;
      if(nGammaLaw!=0){//set P,E,kappa,gamma, Q, and L
        
        for(int j=0;j<nSizeY2;j++){
          
          dTheta_jp1half=dGrid[bin->nTheta][0][j][0];
          if(j==0){
            dTheta_jm1half=dGrid[bin->nTheta][0][j][0]-(dGrid[bin->nTheta][0][j+1][0]-dGrid[bin->nTheta][0][j][0]);
          }
          else{
            dTheta_jm1half=dGrid[bin->nTheta][0][j-1][0];
          }
          dTheta_j=(dTheta_jp1half+dTheta_jm1half)*0.5;
          dA_jp1half=sin(dTheta_jp1half);
          dA_jm1half=sin(dTheta_jm1half);
          dA_j=sin(dTheta_j);
          
          //get P, E, kappa, and gamma
         eosTable->getPEKappaGammaCp(dGrid[bin->nT][i][j][nPlaneIndex],dGrid[bin->nD][i][j][nPlaneIndex],dP
            ,dE,dKappa,dGamma,dCp);
          
          //calculate Q
          dC=sqrt(dGamma*dP/dGrid[bin->nD][i][j][nPlaneIndex]);
          dDVDtThreshold=dAVThreshold*dC;
          
          //Q0
          if(nNumDims>=1){
            if(i==nSizeX1){
              dDVDt=(dA_ip1half*dGrid[bin->nU][i+1][j][nPlaneIndex]
                -dA_im1half*dGrid[bin->nU][i][0][0])/dRSq_i;
            }
            else{
              dDVDt=(dA_ip1half*dGrid[bin->nU][i+1][j][nPlaneIndex]
                -dA_im1half*dGrid[bin->nU][i][j][nPlaneIndex])/dRSq_i;
            }
            if(dDVDt<-1.0*dDVDtThreshold){//being compressed
              dDVDt_mthreshold=dDVDt+dDVDtThreshold;
              dQ0=dASq*dGrid[bin->nD][i][j][nPlaneIndex]*dDVDt_mthreshold*dDVDt_mthreshold;
            }
            else{
              dQ0=0.0;
            }
          }
          
          //Q1
          if(nNumDims>=2){
            if(j==0){
              dDVDt=(dA_jp1half*dGrid[bin->nV][i][j][nPlaneIndex]
                -dA_jm1half*dGrid[bin->nV][i][nSizeY2-1][nPlaneIndex])/dA_j;
            }
            else{
              dDVDt=(dA_jp1half*dGrid[bin->nV][i][j][nPlaneIndex]
                -dA_jm1half*dGrid[bin->nV][i][j-1][nPlaneIndex])/dA_j;
            }
            if(dDVDt<-1.0*dDVDtThreshold){//being compressed
              dDVDt_mthreshold=dDVDt+dDVDtThreshold;
              dQ1=dASq*dGrid[bin->nD][i][j][nPlaneIndex]
                *dDVDt_mthreshold*dDVDt_mthreshold;
            }
            else{
              dQ1=0.0;
            }
          }
          
          //Q2
          if(nNumDims==3){
            if(nPlaneIndex==0){
              dDVDt=(dGrid[bin->nW][i][j][nPlaneIndex]-dGrid[bin->nW][i][j][nSizeZ-1]);
            }
            else{
              dDVDt=(dGrid[bin->nW][i][j][nPlaneIndex]-dGrid[bin->nW][i][j][nPlaneIndex-1]);
            }
            if(dDVDt<-1.0*dDVDtThreshold){//being compressed
              dDVDt_mthreshold=dDVDt+dDVDtThreshold;
              dQ2=dASq*dGrid[bin->nD][i][j][nPlaneIndex]*dDVDt_mthreshold*dDVDt_mthreshold;
            }
            else{
              dQ2=0.0;
            }
          }
          dQ=dQ0+dQ1+dQ2;
          
          dDAve+=dGrid[bin->nD][i][j][nPlaneIndex];
          dEAve+=dE;
          dTAve+=dGrid[bin->nT][i][j][nPlaneIndex];
          dPAve+=dP;
          dQAve+=dQ;
          dKappaAve+=dKappa;
          dGammaAve+=dGamma;
          nCount++;
        }
        dDAve=dDAve/double(nCount);
        dEAve=dEAve/double(nCount);
        dTAve=dTAve/double(nCount);
        dPAve=dPAve/double(nCount);
        dQAve=dQAve/double(nCount);
        dKappaAve=dKappaAve/double(nCount);
        dGammaAve=dGammaAve/double(nCount);
      }
      else{//set P and Q
        
        for(int j=0;j<nSizeY2;j++){
          
          if(nNumDims>1){
            dTheta_jp1half=dGrid[bin->nTheta][0][j][0];
            if(j==0){
              dTheta_jm1half=dGrid[bin->nTheta][0][j][0]-(dGrid[bin->nTheta][0][j+1][0]-dGrid[bin->nTheta][0][j][0]);
            }
            else{
              dTheta_jm1half=dGrid[bin->nTheta][0][j-1][0];
            }
            dTheta_j=(dTheta_jp1half+dTheta_jm1half)*0.5;
            dA_jp1half=sin(dTheta_jp1half);
            dA_jm1half=sin(dTheta_jm1half);
            dA_j=sin(dTheta_j);
          }
          
          //get P
          dP=dGrid[bin->nD][i][j][nPlaneIndex]*(dGamma-1.0)*dGrid[bin->nE][i][j][nPlaneIndex];
          
          //calculate Q
          dC=sqrt(dGamma*dP/dGrid[bin->nD][i][j][nPlaneIndex]);
          dDVDtThreshold=dAVThreshold*dC;
          
          //Q0
          if(nNumDims>=1){
            if(i==nSizeX1){
              dDVDt=(dA_ip1half*dGrid[bin->nU][i+1][j][nPlaneIndex]
                -dA_im1half*dGrid[bin->nU][i][0][0])/dRSq_i;
            }
            else{
              dDVDt=(dA_ip1half*dGrid[bin->nU][i+1][j][nPlaneIndex]
                -dA_im1half*dGrid[bin->nU][i][j][nPlaneIndex])/dRSq_i;
            }
            if(dDVDt<-1.0*dDVDtThreshold){//being compressed
              dDVDt_mthreshold=dDVDt+dDVDtThreshold;
              dQ0=dASq*dGrid[bin->nD][i][j][nPlaneIndex]*dDVDt_mthreshold*dDVDt_mthreshold;
            }
            else{
              dQ0=0.0;
            }
          }
          
          //Q1
          if(nNumDims>=2){
            if(j==0){
              dDVDt=(dA_jp1half*dGrid[bin->nV][i][j][nPlaneIndex]
                -dA_jm1half*dGrid[bin->nV][i][nSizeY2-1][nPlaneIndex])/dA_j;
            }
            else{
              dDVDt=(dA_jp1half*dGrid[bin->nV][i][j][nPlaneIndex]
                -dA_jm1half*dGrid[bin->nV][i][j-1][nPlaneIndex])/dA_j;
            }
            if(dDVDt<-1.0*dDVDtThreshold){//being compressed
              dDVDt_mthreshold=dDVDt+dDVDtThreshold;
              dQ1=dASq*dGrid[bin->nD][i][j][nPlaneIndex]
                *dDVDt_mthreshold*dDVDt_mthreshold;
            }
            else{
              dQ1=0.0;
            }
          }
          
          //Q2
          if(nNumDims==3){
            if(nPlaneIndex==0){
              dDVDt=(dGrid[bin->nW][i][j][nPlaneIndex]-dGrid[bin->nW][i][j][nSizeZ-1]);
            }
            else{
              dDVDt=(dGrid[bin->nW][i][j][nPlaneIndex]-dGrid[bin->nW][i][j][nPlaneIndex-1]);
            }
            if(dDVDt<-1.0*dDVDtThreshold){//being compressed
              dDVDt_mthreshold=dDVDt+dDVDtThreshold;
              dQ2=dASq*dGrid[bin->nD][i][j][nPlaneIndex]*dDVDt_mthreshold*dDVDt_mthreshold;
            }
            else{
              dQ2=0.0;
            }
          }
          dQ=dQ0+dQ1+dQ2;
          
          dDAve+=dGrid[bin->nD][i][j][nPlaneIndex];
          dEAve+=dGrid[bin->nE][i][j][nPlaneIndex];
          dPAve+=dP;
          dQAve+=dQ;
          nCount++;
        }
        dDAve=dDAve/double(nCount);
        dEAve=dEAve/double(nCount);
        dPAve=dPAve/double(nCount);
        dQAve=dQAve/double(nCount);
      }
      
      for(int j=0;j<nSizeY2;j++){
        dTheta_jp1half=dGrid[bin->nTheta][0][j][0];
        if(j==0){
          dTheta_jm1half=dGrid[bin->nTheta][0][j][0]-(dGrid[bin->nTheta][0][j+1][0]-dGrid[bin->nTheta][0][j][0]);
        }
        else{
          dTheta_jm1half=dGrid[bin->nTheta][0][j-1][0];
        }
        dTheta_j=(dTheta_jp1half+dTheta_jm1half)*0.5;
        dA_jp1half=sin(dTheta_jp1half);
        dA_jm1half=sin(dTheta_jm1half);
        dA_j=sin(dTheta_j);
        
        if(i==nSizeX1){//first zone inside 3D region, has 1D u_im1half
          ofFile<<std::setw(nWidthOutputField)<<dGrid[bin->nU][i][0][0];//0
        }
        else{
          ofFile<<std::setw(nWidthOutputField)<<dGrid[bin->nU][i][j][nPlaneIndex];//0
        }
        ofFile
          <<std::setw(nWidthOutputField)<<dGrid[bin->nU0][i][0][0]//1
          <<std::setw(nWidthOutputField)<<dGrid[bin->nV][i][j][nPlaneIndex];//2
        if (bin->nW!=-1){
          ofFile
            <<std::setw(nWidthOutputField)<<dGrid[bin->nW][i][j][nPlaneIndex];//3
        }
        else{
          ofFile
            <<std::setw(nWidthOutputField)<<"-";//3
        }
        ofFile
          <<std::setw(nWidthOutputField)<<dGrid[bin->nD][i][j][nPlaneIndex]//4
          <<std::setw(nWidthOutputField)<<(dGrid[bin->nD][i][j][nPlaneIndex]-dDAve)/dDAve;//5
        if(nGammaLaw!=0){
          //get P, E, kappa, and gamma
         eosTable->getPEKappaGammaCp(dGrid[bin->nT][i][j][nPlaneIndex],dGrid[bin->nD][i][j][nPlaneIndex],dP
            ,dE,dKappa,dGamma,dCp);
          
          //calculate Q, uses dGamma
          dC=sqrt(dGamma*dP/dGrid[bin->nD][i][j][nPlaneIndex]);
          dDVDtThreshold=dAVThreshold*dC;
          
          //Q0
          if(nNumDims>=1){
            if(i==nSizeX1){
              dDVDt=(dA_ip1half*dGrid[bin->nU][i+1][j][nPlaneIndex]
                -dA_im1half*dGrid[bin->nU][i][0][0])/dRSq_i;
            }
            else{
              dDVDt=(dA_ip1half*dGrid[bin->nU][i+1][j][nPlaneIndex]
                -dA_im1half*dGrid[bin->nU][i][j][nPlaneIndex])/dRSq_i;
            }
            if(dDVDt<-1.0*dDVDtThreshold){//being compressed
              dDVDt_mthreshold=dDVDt+dDVDtThreshold;
              dQ0=dASq*dGrid[bin->nD][i][j][nPlaneIndex]*dDVDt_mthreshold*dDVDt_mthreshold;
            }
            else{
              dQ0=0.0;
            }
          }
          
          //Q1
          if(nNumDims>=2){
            if(j==0){
              dDVDt=(dA_jp1half*dGrid[bin->nV][i][j][nPlaneIndex]
                -dA_jm1half*dGrid[bin->nV][i][nSizeY2-1][nPlaneIndex])/dA_j;
            }
            else{
              dDVDt=(dA_jp1half*dGrid[bin->nV][i][j][nPlaneIndex]
                -dA_jm1half*dGrid[bin->nV][i][j-1][nPlaneIndex])/dA_j;
            }
            if(dDVDt<-1.0*dDVDtThreshold){//being compressed
              dDVDt_mthreshold=dDVDt+dDVDtThreshold;
              dQ1=dASq*dGrid[bin->nD][i][j][nPlaneIndex]
                *dDVDt_mthreshold*dDVDt_mthreshold;
            }
            else{
              dQ1=0.0;
            }
          }
          
          //Q1
          if(nNumDims==3){
            if(nPlaneIndex==0){
              dDVDt=(dGrid[bin->nW][i][j][nPlaneIndex]-dGrid[bin->nW][i][j][nSizeZ-1]);
            }
            else{
              dDVDt=(dGrid[bin->nW][i][j][nPlaneIndex]-dGrid[bin->nW][i][j][nPlaneIndex-1]);
            }
            if(dDVDt<-1.0*dDVDtThreshold){//being compressed
              dDVDt_mthreshold=dDVDt+dDVDtThreshold;
              dQ2=dASq*dGrid[bin->nD][i][j][nPlaneIndex]*dDVDt_mthreshold*dDVDt_mthreshold;
            }
            else{
              dQ2=0.0;
            }
          }
          dQ=dQ0+dQ1+dQ2;
          
          ofFile
            <<std::setw(nWidthOutputField)<<dE<<std::setw(nWidthOutputField)<<(dE-dEAve)/dEAve//6,7
            <<std::setw(nWidthOutputField)<<dGrid[bin->nT][i][j][nPlaneIndex]//8
            <<std::setw(nWidthOutputField)<<(dGrid[bin->nT][i][j][nPlaneIndex]-dTAve)/dTAve//9
            <<std::setw(nWidthOutputField)<<dP<<std::setw(nWidthOutputField)<<(dP-dPAve)/dPAve;//10,11
          if(dQAve==0.0){/*if QAve is zero, this can only be the case if all Q's are zero in the 
            horizontal zone also*/
            ofFile
              <<std::setw(nWidthOutputField)<<dQ<<std::setw(nWidthOutputField)<<0.0;//12,13
          }
          else{
            ofFile
              <<std::setw(nWidthOutputField)<<dQ<<std::setw(nWidthOutputField)<<(dQ-dQAve)/dQAve;//12,13
          }
          ofFile
            <<std::setw(nWidthOutputField)<<dKappa<<std::setw(nWidthOutputField)//14
            <<(dKappa-dKappaAve)/dKappaAve//15
            <<std::setw(nWidthOutputField)<<dGamma<<std::setw(nWidthOutputField)//16
            <<(dGamma-dGammaAve)/dGammaAve<<std::setw(nWidthOutputField)<<dCp<<std::endl;//17,18
        }
        else{
          
          //get P
          dP=dGrid[bin->nD][i][j][nPlaneIndex]*(dGamma-1.0)*dGrid[bin->nE][i][j][nPlaneIndex];
          
          //calculate Q
          dC=sqrt(dGamma*dP/dGrid[bin->nD][i][j][nPlaneIndex]);
          dDVDtThreshold=dAVThreshold*dC;
          
          //Q0
          if(nNumDims>=1){
            if(i==nSizeX1){
              dDVDt=(dA_ip1half*dGrid[bin->nU][i+1][j][nPlaneIndex]
                -dA_im1half*dGrid[bin->nU][i][0][0])/dRSq_i;
            }
            else{
              dDVDt=(dA_ip1half*dGrid[bin->nU][i+1][j][nPlaneIndex]
                -dA_im1half*dGrid[bin->nU][i][j][nPlaneIndex])/dRSq_i;
            }
            if(dDVDt<-1.0*dDVDtThreshold){//being compressed
              dDVDt_mthreshold=dDVDt+dDVDtThreshold;
              dQ0=dASq*dGrid[bin->nD][i][j][nPlaneIndex]*dDVDt_mthreshold*dDVDt_mthreshold;
            }
            else{
              dQ0=0.0;
            }
          }
          
          //Q1
          if(nNumDims>=2){
            if(j==0){
              dDVDt=(dA_jp1half*dGrid[bin->nV][i][j][nPlaneIndex]
                -dA_jm1half*dGrid[bin->nV][i][nSizeY2-1][nPlaneIndex])/dA_j;
            }
            else{
              dDVDt=(dA_jp1half*dGrid[bin->nV][i][j][nPlaneIndex]
                -dA_jm1half*dGrid[bin->nV][i][j-1][nPlaneIndex])/dA_j;
            }
            if(dDVDt<-1.0*dDVDtThreshold){//being compressed
              dDVDt_mthreshold=dDVDt+dDVDtThreshold;
              dQ1=dASq*dGrid[bin->nD][i][j][nPlaneIndex]
                *dDVDt_mthreshold*dDVDt_mthreshold;
            }
            else{
              dQ1=0.0;
            }
          }
          
          //Q2
          if(nNumDims==3){
            if(nPlaneIndex==0){
              dDVDt=(dGrid[bin->nW][i][j][nPlaneIndex]-dGrid[bin->nW][i][j][nSizeZ-1]);
            }
            else{
              dDVDt=(dGrid[bin->nW][i][j][nPlaneIndex]-dGrid[bin->nW][i][j][nPlaneIndex-1]);
            }
            if(dDVDt<-1.0*dDVDtThreshold){//being compressed
              dDVDt_mthreshold=dDVDt+dDVDtThreshold;
              dQ2=dASq*dGrid[bin->nD][i][j][nPlaneIndex]*dDVDt_mthreshold*dDVDt_mthreshold;
            }
            else{
              dQ2=0.0;
            }
          }
          dQ=dQ0+dQ1+dQ2;
          
          ofFile
            <<std::setw(nWidthOutputField)<<dE<<std::setw(nWidthOutputField)<<(dE-dEAve)/dEAve//6,7
            <<std::setw(nWidthOutputField)<<"-"<<std::setw(nWidthOutputField)<<"-"//8,9
            <<std::setw(nWidthOutputField)<<dP<<std::setw(nWidthOutputField)<<(dP-dPAve)/dPAve;//10,11
          if(dQAve==0.0){/*if QAve is zero, this can only be the case if all Q's are zero in the 
            horizontal zone also*/
            ofFile
              <<std::setw(nWidthOutputField)<<dQ<<std::setw(nWidthOutputField)<<0.0;//12,13
          }
          else{
            ofFile
              <<std::setw(nWidthOutputField)<<dQ<<std::setw(nWidthOutputField)<<(dQ-dQAve)/dQAve;//12,13
          }
          ofFile
            <<std::setw(nWidthOutputField)<<"-"<<std::setw(nWidthOutputField)<<"-"//14,15
            <<std::setw(nWidthOutputField)<<dGamma<<std::setw(nWidthOutputField)<<0.0//16,17
            <<std::setw(nWidthOutputField)<<dCp//18
            <<std::endl;
        }
      }
      ofFile
        <<std::setw(nWidthOutputField)<<"-"//0
        <<std::setw(nWidthOutputField)<<"-"//1
        <<std::setw(nWidthOutputField)<<dGrid[bin->nV][i][0][nPlaneIndex]//2
        <<std::setw(nWidthOutputField)<<"-"//3
        <<std::setw(nWidthOutputField)<<"-"<<std::setw(nWidthOutputField)//4
        <<"-"//5
        <<std::setw(nWidthOutputField)<<"-"<<std::setw(nWidthOutputField)//6
        <<"-"//7
        <<std::setw(nWidthOutputField)<<"-"<<std::setw(nWidthOutputField)//8
        <<"-"//9
        <<std::setw(nWidthOutputField)<<"-"<<std::setw(nWidthOutputField)//10
        <<"-"//11
        <<std::setw(nWidthOutputField)<<"-"<<std::setw(nWidthOutputField)//12
        <<"-"//13
        <<std::setw(nWidthOutputField)<<"-"<<std::setw(nWidthOutputField)//14
        <<"-"//15
        <<std::setw(nWidthOutputField)<<"-"<<std::setw(nWidthOutputField)//16
        <<"-"<<std::setw(nWidthOutputField)<<"-"//17,18
        <<std::endl;
    }
    for(int j=0;j<nSizeY2;j++){
      ofFile
        <<std::setw(nWidthOutputField)<<dGrid[bin->nU][nSizeX2][j][nPlaneIndex]//0
        <<std::setw(nWidthOutputField)<<dGrid[bin->nU0][nSizeX2][0][0]//1
        <<std::setw(nWidthOutputField)<<"-"//2
        <<std::setw(nWidthOutputField)<<"-"//3
        <<std::setw(nWidthOutputField)<<"-"<<std::setw(nWidthOutputField)//4
        <<"-"//5
        <<std::setw(nWidthOutputField)<<"-"<<std::setw(nWidthOutputField)//6
        <<"-"//7
        <<std::setw(nWidthOutputField)<<"-"<<std::setw(nWidthOutputField)//8
        <<"-"//9
        <<std::setw(nWidthOutputField)<<"-"<<std::setw(nWidthOutputField)//10
        <<"-"//11
        <<std::setw(nWidthOutputField)<<"-"<<std::setw(nWidthOutputField)//12
        <<"-"//13
        <<std::setw(nWidthOutputField)<<"-"<<std::setw(nWidthOutputField)//14
        <<"-"//15
        <<std::setw(nWidthOutputField)<<"-"<<std::setw(nWidthOutputField)//16
        <<"-"<<std::setw(nWidthOutputField)<<"-"<<std::endl;//17,18
    }
  }
  if(nPlane==1){//theta-phi
    
    //i
    ofFile<<"im1half[g](0) "<<nPlaneIndex<<" "
      <<nPlaneIndex+1<<std::endl;
    
    //M_r
    ofFile<<"M_r_im1half[g](1) "<<dGrid[bin->nM][nPlaneIndex][0][0]<<" "
      <<dGrid[bin->nM][nPlaneIndex+1][0][0]<<std::endl;
    
    //R
    ofFile<<"R_im1half[cm](2) "<<dGrid[bin->nR][nPlaneIndex][0][0]<<" "
      <<dGrid[bin->nR][nPlaneIndex+1][0][0]<<std::endl;
    
    //j
    ofFile<<"jm1half(3) ";
    for( int j=0;j<nSize[bin->nD][1]+1+2*nNumGhostCells;j++){//add an extra since it is an interface
      ofFile<<j<<" ";
    }
    ofFile<<std::endl;
    
    //write out theta
    int nInterFaceY=0;
    if(nPeriodic[1]==0){
      nInterFaceY=1;
    }
    ofFile<<"theta_jm1half[rad](4) ";
    if(nPeriodic[1]==1){/*write out inner interface if periodic, if not periodic, it is already 
      included*/
      double dInnerTheta=dGrid[bin->nTheta][0][0][0]-(dGrid[bin->nTheta][0][nSize[bin->nD][1]+nNumGhostCells-2][0]
        -dGrid[bin->nTheta][0][nSize[bin->nD][1]+nNumGhostCells-3][0]);
      ofFile<<dInnerTheta<<" ";
     
    }
    for( int j=0;j<nSize[bin->nD][1]+nInterFaceY+2*nNumGhostCells;j++){
      ofFile<<dGrid[bin->nTheta][0][j][0]<<" ";
    }
    ofFile<<std::endl;
    
    //k
    ofFile<<"km1half(5) ";
    for( int k=0;k<nSize[bin->nD][2]+1+2*nNumGhostCells;k++){//add an extra since it is an interface
      ofFile<<k<<" ";
    }
    ofFile<<std::endl;
    
    //write out phi
    int nInterFaceZ=0;
    if(nPeriodic[2]==0){
      nInterFaceZ=1;
    }
    ofFile<<"phi_km1half[rad](6) ";
    if(nPeriodic[2]==1){/*write out inner interface if periodic, if not periodic, it is already
      included*/
      double dInnerTheta=dGrid[bin->nPhi][0][0][0]-(dGrid[bin->nPhi][0][0][nSize[bin->nD][2]+nNumGhostCells-2]
        -dGrid[bin->nPhi][0][0][nSize[bin->nD][2]+nNumGhostCells-3]);
      ofFile<<dInnerTheta<<" ";
     
    }
    for( int k=0;k<nSize[bin->nD][2]+nInterFaceZ+2*nNumGhostCells;k++){
      ofFile<<dGrid[bin->nPhi][0][0][k]<<" ";
    }
    ofFile<<std::endl;
    
    ofFile
      <<std::setw(nWidthOutputField)<<"U_im1halfjk[cm/s](7)"
      <<std::setw(nWidthOutputField)<<"U0_im1half[cm/s](8)"
      <<std::setw(nWidthOutputField)<<"V_ijm1halfk[cm/s](9)"
      <<std::setw(nWidthOutputField)<<"W_ijkm1half[cm/s](10)"
      <<std::setw(nWidthOutputField)<<"D_ijk[g/cm^3](11)"<<std::setw(nWidthOutputField)
      <<"D_rel_dif_hor_ave(12)"
      <<std::setw(nWidthOutputField)<<"E_ijk[erg/g](13)"<<std::setw(nWidthOutputField)
      <<"E_rel_dif_hor_ave(14)"
      <<std::setw(nWidthOutputField)<<"T_ijk[K](15)"<<std::setw(nWidthOutputField)
      <<"T_rel_dif_hor_ave(16)"
      <<std::setw(nWidthOutputField)<<"P_ijk[dynes/cm^2](17)"<<std::setw(nWidthOutputField)
      <<"P_rel_dif_hor_ave(18)"
      <<std::setw(nWidthOutputField)<<"Q_ijk[dynes/cm^2](19)"<<std::setw(nWidthOutputField)
      <<"Q_rel_dif_hor_ave(20)"
      <<std::setw(nWidthOutputField)<<"Kap_ij[cm^2/g](21)"<<std::setw(nWidthOutputField)
      <<"Kap_rel_dif_hor_ave(22)"
      <<std::setw(nWidthOutputField)<<"Gam_ijk[na](23)"<<std::setw(nWidthOutputField)
      <<"Gam_rel_dif_hor_ave(24)"
      <<std::setw(nWidthOutputField)<<"C_P[na](24)"
      <<std::endl;
    
    //3D region
    int nSizeX1=nNum1DZones+nNumGhostCells;
    int nSizeY=nSize[bin->nD][1]+2*nNumGhostCells;
    int nSizeZ=nSize[bin->nD][2]+2*nNumGhostCells;
    int i=nPlaneIndex;
    double dDAve=0.0;
    double dEAve=0.0;
    double dTAve=0.0;
    double dPAve=0.0;
    double dQAve=0.0;
    double dKappaAve=0.0;
    double dGammaAve=0.0;
    int nCount=0;
    double dCp;
    dR_i=(dGrid[bin->nR][i+1][0][0]+dGrid[bin->nR][i][0][0])*0.5;
    dRSq_ip1half=dGrid[bin->nR][i+1][0][0]*dGrid[bin->nR][i+1][0][0];
    dRSq_i=dR_i*dR_i;
    dA_ip1half=dGrid[bin->nR][i+1][0][0]*dGrid[bin->nR][i+1][0][0];
    dA_im1half=dGrid[bin->nR][i][0][0]*dGrid[bin->nR][i][0][0];
    
    //calculate horizontal average of various quantities
    if(nGammaLaw!=0){//set P,E,kappa,gamma, Q, and L
      
      for(int j=0;j<nSizeY;j++){
      
        //calculate some theta areas
        dTheta_jp1half=dGrid[bin->nTheta][0][j][0];
        if(j==0){
          dTheta_jm1half=dGrid[bin->nTheta][0][j][0]-(dGrid[bin->nTheta][0][j+1][0]-dGrid[bin->nTheta][0][j][0]);
        }
        else{
          dTheta_jm1half=dGrid[bin->nTheta][0][j-1][0];
        }
        dTheta_j=(dTheta_jp1half+dTheta_jm1half)*0.5;
        dA_jp1half=sin(dTheta_jp1half);
        dA_jm1half=sin(dTheta_jm1half);
        dA_j=sin(dTheta_j);
        
        for(int k=0;k<nSizeZ;k++){
          
          //get P, E, kappa, and gamma
         eosTable->getPEKappaGammaCp(dGrid[bin->nT][i][j][k],dGrid[bin->nD][i][j][k],dP,dE,dKappa,dGamma,dCp);
          
          //calculate Q
          dC=sqrt(dGamma*dP/dGrid[bin->nD][i][j][k]);
          dDVDtThreshold=dAVThreshold*dC;
          
          //Q0
          if(nNumDims>=1){
            if(i==nSizeX1){
              dDVDt=(dA_ip1half*dGrid[bin->nU][i+1][j][k]
                -dA_im1half*dGrid[bin->nU][i][0][0])/dRSq_i;
            }
            else{
              dDVDt=(dA_ip1half*dGrid[bin->nU][i+1][j][k]
                -dA_im1half*dGrid[bin->nU][i][j][k])/dRSq_i;
            }
            if(dDVDt<-1.0*dDVDtThreshold){//being compressed
              dDVDt_mthreshold=dDVDt+dDVDtThreshold;
              dQ0=dASq*dGrid[bin->nD][i][j][k]*dDVDt_mthreshold*dDVDt_mthreshold;
            }
            else{
              dQ0=0.0;
            }
          }
          
          //Q1
          if(nNumDims>=2){
            if(j==0){
              dDVDt=(dA_jp1half*dGrid[bin->nV][i][j][k]
                -dA_jm1half*dGrid[bin->nV][i][nSizeY-1][k])/dA_j;
            }
            else{
              dDVDt=(dA_jp1half*dGrid[bin->nV][i][j][k]
                -dA_jm1half*dGrid[bin->nV][i][j-1][k])/dA_j;
            }
            if(dDVDt<-1.0*dDVDtThreshold){//being compressed
              dDVDt_mthreshold=dDVDt+dDVDtThreshold;
              dQ1=dASq*dGrid[bin->nD][i][j][k]
                *dDVDt_mthreshold*dDVDt_mthreshold;
            }
            else{
              dQ1=0.0;
            }
          }
          
          //Q2
          if(nNumDims==3){
            if(k==0){
              dDVDt=(dGrid[bin->nW][i][j][k]-dGrid[bin->nW][i][j][nSizeZ-1]);
            }
            else{
              dDVDt=(dGrid[bin->nW][i][j][k]-dGrid[bin->nW][i][j][k-1]);
            }
            if(dDVDt<-1.0*dDVDtThreshold){//being compressed
              dDVDt_mthreshold=dDVDt+dDVDtThreshold;
              dQ2=dASq*dGrid[bin->nD][i][j][k]*dDVDt_mthreshold*dDVDt_mthreshold;
            }
            else{
              dQ2=0.0;
            }
          }
          dQ=dQ0+dQ1+dQ2;
          
          dDAve+=dGrid[bin->nD][i][j][k];
          dEAve+=dE;
          dTAve+=dGrid[bin->nT][i][j][k];
          dPAve+=dP;
          dQAve+=dQ;
          dKappaAve+=dKappa;
          dGammaAve+=dGamma;
          nCount++;
        }
      }
      dDAve=dDAve/double(nCount);
      dEAve=dEAve/double(nCount);
      dTAve=dTAve/double(nCount);
      dPAve=dPAve/double(nCount);
      dQAve=dQAve/double(nCount);
      dKappaAve=dKappaAve/double(nCount);
      dGammaAve=dGammaAve/double(nCount);
    }
    else{//set P and Q
      
      for(int j=0;j<nSizeY;j++){
      
        //calculate some theta areas
        dTheta_jp1half=dGrid[bin->nTheta][0][j][0];
        if(j==0){
          dTheta_jm1half=dGrid[bin->nTheta][0][j][0]-(dGrid[bin->nTheta][0][j+1][0]-dGrid[bin->nTheta][0][j][0]);
        }
        else{
          dTheta_jm1half=dGrid[bin->nTheta][0][j-1][0];
        }
        dTheta_j=(dTheta_jp1half+dTheta_jm1half)*0.5;
        dA_jp1half=sin(dTheta_jp1half);
        dA_jm1half=sin(dTheta_jm1half);
        dA_j=sin(dTheta_j);
        
        for(int k=0;k<nSizeZ;k++){
          
          //get P
          dP=dGrid[bin->nD][i][j][k]*(dGamma-1.0)*dGrid[bin->nE][i][j][k];
          
          //calculate Q
          dC=sqrt(dGamma*dP/dGrid[bin->nD][i][j][k]);
          dDVDtThreshold=dAVThreshold*dC;
          
          
          //Q0
          if(nNumDims>=1){
            if(i==nSizeX1){
              dDVDt=(dA_ip1half*dGrid[bin->nU][i+1][j][k]
                -dA_im1half*dGrid[bin->nU][i][0][0])/dRSq_i;
            }
            else{
              dDVDt=(dA_ip1half*dGrid[bin->nU][i+1][j][k]
                -dA_im1half*dGrid[bin->nU][i][j][k])/dRSq_i;
            }
            if(dDVDt<-1.0*dDVDtThreshold){//being compressed
              dDVDt_mthreshold=dDVDt+dDVDtThreshold;
              dQ0=dASq*dGrid[bin->nD][i][j][k]*dDVDt_mthreshold*dDVDt_mthreshold;
            }
            else{
              dQ0=0.0;
            }
          }
          
          //Q1
          if(nNumDims>=2){
            if(j==0){
              dDVDt=(dA_jp1half*dGrid[bin->nV][i][j][k]
                -dA_jm1half*dGrid[bin->nV][i][nSizeY-1][k])/dA_j;
            }
            else{
              dDVDt=(dA_jp1half*dGrid[bin->nV][i][j][k]
                -dA_jm1half*dGrid[bin->nV][i][j-1][k])/dA_j;
            }
            if(dDVDt<-1.0*dDVDtThreshold){//being compressed
              dDVDt_mthreshold=dDVDt+dDVDtThreshold;
              dQ1=dASq*dGrid[bin->nD][i][j][k]
                *dDVDt_mthreshold*dDVDt_mthreshold;
            }
            else{
              dQ1=0.0;
            }
          }
          
          //Q2
          if(nNumDims==3){
            if(k==0){
              dDVDt=(dGrid[bin->nW][i][j][k]-dGrid[bin->nW][i][j][nSizeZ-1]);
            }
            else{
              dDVDt=(dGrid[bin->nW][i][j][k]-dGrid[bin->nW][i][j][k-1]);
            }
            if(dDVDt<-1.0*dDVDtThreshold){//being compressed
              dDVDt_mthreshold=dDVDt+dDVDtThreshold;
              dQ2=dASq*dGrid[bin->nD][i][j][k]*dDVDt_mthreshold*dDVDt_mthreshold;
            }
            else{
              dQ2=0.0;
            }
          }
          dQ=dQ0+dQ1+dQ2;
          
          dDAve+=dGrid[bin->nD][i][j][k];
          dEAve+=dGrid[bin->nE][i][j][k];
          dPAve+=dP;
          dQAve+=dQ;
          nCount++;
        }
      }
      dDAve=dDAve/double(nCount);
      dEAve=dEAve/double(nCount);
      dPAve=dPAve/double(nCount);
      dQAve=dQAve/double(nCount);
    }
    
    for(int j=0;j<nSizeY;j++){
      
      //calculate some theta areas
      dTheta_jp1half=dGrid[bin->nTheta][0][j][0];
      if(j==0){
        dTheta_jm1half=dGrid[bin->nTheta][0][j][0]-(dGrid[bin->nTheta][0][j+1][0]-dGrid[bin->nTheta][0][j][0]);
      }
      else{
        dTheta_jm1half=dGrid[bin->nTheta][0][j-1][0];
      }
      dTheta_j=(dTheta_jp1half+dTheta_jm1half)*0.5;
      dA_jp1half=sin(dTheta_jp1half);
      dA_jm1half=sin(dTheta_jm1half);
      dA_j=sin(dTheta_j);
      
      for(int k=0;k<nSizeZ;k++){
        
        ofFile
          <<std::setw(nWidthOutputField)<<dGrid[bin->nU][i][j][k]//0
          <<std::setw(nWidthOutputField)<<dGrid[bin->nU0][i][0][0]//1
          <<std::setw(nWidthOutputField)<<dGrid[bin->nV][i][j][k];//2
        if(bin->nW!=-1){
          ofFile
            <<std::setw(nWidthOutputField)<<dGrid[bin->nW][i][j][k];//3
        }
        else{
          ofFile
            <<std::setw(nWidthOutputField)<<"-";//3
        }
        ofFile
          <<std::setw(nWidthOutputField)<<dGrid[bin->nD][i][j][k]//4
          <<std::setw(nWidthOutputField)<<(dGrid[bin->nD][i][j][k]-dDAve)/dDAve;//5
        if(nGammaLaw!=0){
          
          //get P,E,Kappa,Gamma, calculate luminosity from cell and add to sum
         eosTable->getPEKappaGammaCp(dGrid[bin->nT][i][j][k],dGrid[bin->nD][i][j][k],dP,dE,dKappa
            ,dGamma,dCp);
          
          //calculate Q
          dC=sqrt(dGamma*dP/dGrid[bin->nD][i][j][k]);
          dDVDtThreshold=dAVThreshold*dC;
          if(nNumDims>=1){
            if(i==nSizeX1){
              dDVDt=(dA_ip1half*dGrid[bin->nU][i+1][j][k]
                -dA_im1half*dGrid[bin->nU][i][0][0])/dRSq_i;
            }
            else{
              dDVDt=(dA_ip1half*dGrid[bin->nU][i+1][j][k]
                -dA_im1half*dGrid[bin->nU][i][j][k])/dRSq_i;
            }
            if(dDVDt<-1.0*dDVDtThreshold){//being compressed
              dDVDt_mthreshold=dDVDt+dDVDtThreshold;
              dQ0=dASq*dGrid[bin->nD][i][j][k]*dDVDt_mthreshold*dDVDt_mthreshold;
            }
            else{
              dQ0=0.0;
            }
          }
          else if(nNumDims>=2){
            if(j==0){
              dDVDt=(dA_jp1half*dGrid[bin->nV][i][j][k]
                -dA_jm1half*dGrid[bin->nV][i][nSizeY-1][k])/dA_j;
            }
            else{
              dDVDt=(dA_jp1half*dGrid[bin->nV][i][j][k]
                -dA_jm1half*dGrid[bin->nV][i][j-1][k])/dA_j;
            }
            if(dDVDt<-1.0*dDVDtThreshold){//being compressed
              dDVDt_mthreshold=dDVDt+dDVDtThreshold;
              dQ1=dASq*dGrid[bin->nD][i][j][k]
                *dDVDt_mthreshold*dDVDt_mthreshold;
            }
            else{
              dQ1=0.0;
            }
          }
          else if(nNumDims==3){
            if(k==0){
              dDVDt=(dGrid[bin->nW][i][j][k]-dGrid[bin->nW][i][j][nSizeZ-1]);
            }
            else{
              dDVDt=(dGrid[bin->nW][i][j][k]-dGrid[bin->nW][i][j][k-1]);
            }
            if(dDVDt<-1.0*dDVDtThreshold){//being compressed
              dDVDt_mthreshold=dDVDt+dDVDtThreshold;
              dQ2=dASq*dGrid[bin->nD][i][j][k]*dDVDt_mthreshold*dDVDt_mthreshold;
            }
            else{
              dQ2=0.0;
            }
          }
          dQ=dQ0+dQ1+dQ2;
          
          ofFile
            <<std::setw(nWidthOutputField)<<dE<<std::setw(nWidthOutputField)<<(dE-dEAve)/dEAve//6,7
            <<std::setw(nWidthOutputField)<<dGrid[bin->nT][i][j][k]//8
            <<std::setw(nWidthOutputField)<<(dGrid[bin->nT][i][j][k]-dTAve)/dTAve//9
            <<std::setw(nWidthOutputField)<<dP<<std::setw(nWidthOutputField)<<(dP-dPAve)/dPAve;//10,11
          if(dQAve==0.0){/*if QAve is zero, this can only be the case if all Q's are zero in the 
            horizontal zone also*/
            ofFile
              <<std::setw(nWidthOutputField)<<dQ<<std::setw(nWidthOutputField)<<0.0;//12,13
          }
          else{
            ofFile
              <<std::setw(nWidthOutputField)<<dQ<<std::setw(nWidthOutputField)<<(dQ-dQAve)/dQAve;//12,13
          }
          ofFile
            <<std::setw(nWidthOutputField)<<dKappa<<std::setw(nWidthOutputField)//14
            <<(dKappa-dKappaAve)/dKappaAve//15
            <<std::setw(nWidthOutputField)<<dGamma<<std::setw(nWidthOutputField)//16
            <<(dGamma-dGammaAve)/dGammaAve<<std::setw(nWidthOutputField)//17
            <<dCp//18
            <<std::endl;
        }
        else{
          
          //get P
          dP=dGrid[bin->nD][i][j][k]*(dGamma-1.0)*dGrid[bin->nE][i][j][k];
          
          //calculate Q
          dC=sqrt(dGamma*dP/dGrid[bin->nD][i][j][k]);
          dDVDtThreshold=dAVThreshold*dC;
          
          //Q0
          if(nNumDims>=1){
            if(i==nSizeX1){
              dDVDt=(dA_ip1half*dGrid[bin->nU][i+1][j][k]
                -dA_im1half*dGrid[bin->nU][i][0][0])/dRSq_i;
            }
            else{
              dDVDt=(dA_ip1half*dGrid[bin->nU][i+1][j][k]
                -dA_im1half*dGrid[bin->nU][i][j][k])/dRSq_i;
            }
            if(dDVDt<-1.0*dDVDtThreshold){//being compressed
              dDVDt_mthreshold=dDVDt+dDVDtThreshold;
              dQ0=dASq*dGrid[bin->nD][i][j][k]*dDVDt_mthreshold*dDVDt_mthreshold;
            }
            else{
              dQ0=0.0;
            }
          }
          
          //Q1
          if(nNumDims>=2){
            if(j==0){
              dDVDt=(dA_jp1half*dGrid[bin->nV][i][j][k]
                -dA_jm1half*dGrid[bin->nV][i][nSizeY-1][k])/dA_j;
            }
            else{
              dDVDt=(dA_jp1half*dGrid[bin->nV][i][j][k]
                -dA_jm1half*dGrid[bin->nV][i][j-1][k])/dA_j;
            }
            if(dDVDt<-1.0*dDVDtThreshold){//being compressed
              dDVDt_mthreshold=dDVDt+dDVDtThreshold;
              dQ1=dASq*dGrid[bin->nD][i][j][k]
                *dDVDt_mthreshold*dDVDt_mthreshold;
            }
            else{
              dQ1=0.0;
            }
          }
          
          //Q2
          if(nNumDims==3){
            if(k==0){
              dDVDt=(dGrid[bin->nW][i][j][k]-dGrid[bin->nW][i][j][nSizeZ-1]);
            }
            else{
              dDVDt=(dGrid[bin->nW][i][j][k]-dGrid[bin->nW][i][j][k-1]);
            }
            if(dDVDt<-1.0*dDVDtThreshold){//being compressed
              dDVDt_mthreshold=dDVDt+dDVDtThreshold;
              dQ2=dASq*dGrid[bin->nD][i][j][k]*dDVDt_mthreshold*dDVDt_mthreshold;
            }
            else{
              dQ2=0.0;
            }
          }
          dQ=dQ0+dQ1+dQ2;
          
          ofFile
            <<std::setw(nWidthOutputField)<<dE<<std::setw(nWidthOutputField)<<(dE-dEAve)/dEAve//6,7
            <<std::setw(nWidthOutputField)<<"-"<<std::setw(nWidthOutputField)<<"-"//8,9
            <<std::setw(nWidthOutputField)<<dP<<std::setw(nWidthOutputField)<<(dP-dPAve)/dPAve;//10,11
          if(dQAve==0.0){/*if QAve is zero, this can only be the case if all Q's are zero in the 
            horizontal zone also*/
            ofFile
              <<std::setw(nWidthOutputField)<<dQ<<std::setw(nWidthOutputField)<<0.0;//12,13
          }
          else{
            ofFile
              <<std::setw(nWidthOutputField)<<dQ<<std::setw(nWidthOutputField)<<(dQ-dQAve)/dQAve;//12,13
          }
          ofFile
            <<std::setw(nWidthOutputField)<<"-"<<std::setw(nWidthOutputField)<<"-"//14,15
            <<std::setw(nWidthOutputField)<<dGamma<<std::setw(nWidthOutputField)<<0.0//16,17
            <<std::setw(nWidthOutputField)<<"-"//18
            <<std::endl;
        }
      }
      ofFile
        <<std::setw(nWidthOutputField)<<"-"//0
        <<std::setw(nWidthOutputField)<<"-"//1
        <<std::setw(nWidthOutputField)<<"-";//2
      if(bin->nW!=-1){
        ofFile
          <<std::setw(nWidthOutputField)<<dGrid[bin->nW][i][j][0];//3
      }
      else{
        ofFile
          <<std::setw(nWidthOutputField)<<"-";//3
      }
      ofFile
        <<std::setw(nWidthOutputField)<<"-"<<std::setw(nWidthOutputField)//4
        <<"-"//5
        <<std::setw(nWidthOutputField)<<"-"<<std::setw(nWidthOutputField)//6
        <<"-"//7
        <<std::setw(nWidthOutputField)<<"-"<<std::setw(nWidthOutputField)//8
        <<"-"//9
        <<std::setw(nWidthOutputField)<<"-"<<std::setw(nWidthOutputField)//10
        <<"-"//11
        <<std::setw(nWidthOutputField)<<"-"<<std::setw(nWidthOutputField)//12
        <<"-"//13
        <<std::setw(nWidthOutputField)<<"-"<<std::setw(nWidthOutputField)//14
        <<"-"//15
        <<std::setw(nWidthOutputField)<<"-"<<std::setw(nWidthOutputField)//16
        <<"-"<<std::setw(nWidthOutputField)<<"-"//17,18
        <<std::endl;
    }
    for(int k=0;k<nSizeZ;k++){
      ofFile
        <<std::setw(nWidthOutputField)<<"-"//0
        <<std::setw(nWidthOutputField)<<"-"//1
        <<std::setw(nWidthOutputField)<<dGrid[bin->nV][i][0][k]//2
        <<std::setw(nWidthOutputField)<<"-"//3
        <<std::setw(nWidthOutputField)<<"-"<<std::setw(nWidthOutputField)//4
        <<"-"//5
        <<std::setw(nWidthOutputField)<<"-"<<std::setw(nWidthOutputField)//6
        <<"-"//7
        <<std::setw(nWidthOutputField)<<"-"<<std::setw(nWidthOutputField)//8
        <<"-"//9
        <<std::setw(nWidthOutputField)<<"-"<<std::setw(nWidthOutputField)//10
        <<"-"//11
        <<std::setw(nWidthOutputField)<<"-"<<std::setw(nWidthOutputField)//12
        <<"-"//13
        <<std::setw(nWidthOutputField)<<"-"<<std::setw(nWidthOutputField)//14
        <<"-"//15
        <<std::setw(nWidthOutputField)<<"-"<<std::setw(nWidthOutputField)//16
        <<"-"<<std::setw(nWidthOutputField)<<"-"//17,18
        <<std::endl;
    }
  }
  if(nPlane==2){//r-phi
     
    //i
    ofFile<<"im1half(0) ";
    for( int i=0;i<nSize[bin->nR][0]+2*nNumGhostCells;i++){
      ofFile<<i<<" ";
    }
    ofFile<<std::endl;
    
    //M_r
    ofFile<<"M_r_im1half[g](1) ";
    for( int i=0;i<nSize[bin->nR][0]+2*nNumGhostCells;i++){
      ofFile<<dGrid[bin->nM][i][0][0]<<" ";
    }
    ofFile<<std::endl;
    
    //R
    ofFile<<"R_im1half[cm](2) ";
    for( int i=0;i<nSize[bin->nR][0]+2*nNumGhostCells;i++){
      ofFile<<dGrid[bin->nR][i][0][0]<<" ";
    }
    ofFile<<std::endl;
    
    //j
    if(nNumDims>2){
      ofFile<<"jm1half(3) "<<nPlaneIndex<<" "<<nPlaneIndex+1<<std::endl;
    }
    else{
      ofFile<<"jm1half(3) "<<0<<" "<<1<<std::endl;
    }
    
    //theta
    ofFile<<"theta_jm1half[rad](4) "<<dGrid[bin->nTheta][0][nPlaneIndex][0]<<" ";
    if(nPlaneIndex==nSize[bin->nTheta][1]+2*nNumGhostCells-1){//if in last zone, need to do something to get outter phi
      ofFile<<dGrid[bin->nTheta][0][nPlaneIndex][0]+(dGrid[bin->nTheta][0][nPlaneIndex][0]
        -dGrid[bin->nTheta][0][nPlaneIndex-1][0])<<std::endl;
    }
    else{
      ofFile<<dGrid[bin->nTheta][0][nPlaneIndex+1][0]<<std::endl;
    }
    
    //k
    ofFile<<"km1half(5) ";
    for( int k=0;k<nSize[bin->nD][2]+1+2*nNumGhostCells;k++){//add an extra since it is an interface
      ofFile<<k<<" ";
    }
    ofFile<<std::endl;
    
    //phi
    int nInterFaceZ=0;
    if(nPeriodic[2]==0){
      nInterFaceZ=1;
    }
    ofFile<<"phi_km1half[rad](6) ";
    if(nPeriodic[2]==1){/*write out inner interface if periodic, if not periodic, it is already
      included*/
      double dInnerTheta=dGrid[bin->nPhi][0][0][0]-(dGrid[bin->nPhi][0][0][nSize[bin->nD][2]+nNumGhostCells-2]
        -dGrid[bin->nPhi][0][0][nSize[bin->nD][2]+nNumGhostCells-3]);
      ofFile<<dInnerTheta<<" ";
     
    }
    for( int k=0;k<nSize[bin->nD][2]+nInterFaceZ+2*nNumGhostCells;k++){
      ofFile<<dGrid[bin->nPhi][0][0][k]<<" ";
    }
    ofFile<<std::endl;
    
    ofFile
      <<std::setw(nWidthOutputField)<<"U_im1halfjk[cm/s](7)"
      <<std::setw(nWidthOutputField)<<"U0_im1half[cm/s](8)"
      <<std::setw(nWidthOutputField)<<"V_ijm1halfk[cm/s](9)"
      <<std::setw(nWidthOutputField)<<"W_ijkm1half[cm/s](10)"
      <<std::setw(nWidthOutputField)<<"D_ijk[g/cm^3](11)"<<std::setw(nWidthOutputField)
      <<"D_rel_dif_hor_ave(12)"
      <<std::setw(nWidthOutputField)<<"E_ijk[erg/g](13)"<<std::setw(nWidthOutputField)
      <<"E_rel_dif_hor_ave(14)"
      <<std::setw(nWidthOutputField)<<"T_ijk[K](15)"<<std::setw(nWidthOutputField)
      <<"T_rel_dif_hor_ave(16)"
      <<std::setw(nWidthOutputField)<<"P_ijk[dynes/cm^2](17)"<<std::setw(nWidthOutputField)
      <<"P_rel_dif_hor_ave(18)"
      <<std::setw(nWidthOutputField)<<"Q_ijk[dynes/cm^2](19)"<<std::setw(nWidthOutputField)
      <<"Q_rel_dif_hor_ave(20)"
      <<std::setw(nWidthOutputField)<<"Kap_ij[cm^2/g](21)"<<std::setw(nWidthOutputField)
      <<"Kap_rel_dif_hor_ave(22)"
      <<std::setw(nWidthOutputField)<<"Gam_ijk[na](23)"<<std::setw(nWidthOutputField)
      <<"Gam_rel_dif_hor_ave(24)"<<std::setw(nWidthOutputField)<<"C_P[na](24)"
      <<std::endl;
    
    //1D region
    int nSizeX1=nNum1DZones+nNumGhostCells;
    int nSizeY2=nSize[bin->nD][1]+2*nNumGhostCells;
    int nSizeZ=nSize[bin->nD][2]+2*nNumGhostCells;
    for(int i=0;i<nSizeX1;i++){
      dR_i=(dGrid[bin->nR][i+1][0][0]+dGrid[bin->nR][i][0][0])*0.5;
      dRSq_ip1half=dGrid[bin->nR][i+1][0][0]*dGrid[bin->nR][i+1][0][0];
      dRSq_i=dR_i*dR_i;
      dA_ip1half=dGrid[bin->nR][i+1][0][0]*dGrid[bin->nR][i+1][0][0];
      dA_im1half=dGrid[bin->nR][i][0][0]*dGrid[bin->nR][i][0][0];
      for(int k=0;k<nSizeZ;k++){
        
        ofFile
          <<std::setw(nWidthOutputField)<<dGrid[bin->nU][i][0][0]//0
          <<std::setw(nWidthOutputField)<<dGrid[bin->nU0][i][0][0]//1
          <<std::setw(nWidthOutputField)<<0.0//2
          <<std::setw(nWidthOutputField)<<0.0//3
          <<std::setw(nWidthOutputField)<<dGrid[bin->nD][i][0][0]<<std::setw(nWidthOutputField)<<0.0;//4,5
        
        if(nGammaLaw!=0){//set P,E,kappa,gamma, Q, and L,Cp
          
          //get P, E, kappa, and gamma
         eosTable->getPEKappaGammaCp(dGrid[bin->nT][i][0][0],dGrid[bin->nD][i][0][0],dP
            ,dE,dKappa,dGamma,dCp);
          
          //calculate Q
          dC=sqrt(dGamma*dP/dGrid[bin->nD][i][0][0]);
          dDVDtThreshold=dAVThreshold*dC;
          dDVDt=(dA_ip1half*dGrid[bin->nU][i+1][0][0]
            -dA_im1half*dGrid[bin->nU][i][0][0])/dRSq_i;
          if(dDVDt<-1.0*dDVDtThreshold){//being compressed
            dDVDt_mthreshold=dDVDt+dDVDtThreshold;
            dQ=dASq*dGrid[bin->nD][i][0][0]*dDVDt_mthreshold*dDVDt_mthreshold;
          }
          else{
            dQ=0.0;
          }
          
          //print them out
          ofFile
            <<std::setw(nWidthOutputField)<<dE<<std::setw(nWidthOutputField)<<0.0//6,7
            <<std::setw(nWidthOutputField)<<dGrid[bin->nT][i][0][0]<<std::setw(nWidthOutputField)<<0.0//8,9
            <<std::setw(nWidthOutputField)<<dP<<std::setw(nWidthOutputField)<<0.0//10,11
            <<std::setw(nWidthOutputField)<<dQ<<std::setw(nWidthOutputField)<<0.0//12,13
            <<std::setw(nWidthOutputField)<<dKappa<<std::setw(nWidthOutputField)<<0.0//14,15
            <<std::setw(nWidthOutputField)<<dGamma<<std::setw(nWidthOutputField)<<0.0//16,17
            <<std::setw(nWidthOutputField)<<dCp;//18
          
        }
        else{//set P and Q
          
          dP=dGrid[bin->nD][i][0][0]*(dGamma-1.0)*dGrid[bin->nE][i][0][0];
          
          
          //calculate Q
          dR_i=(dGrid[bin->nR][i+1][0][0]+dGrid[bin->nR][i][0][0])*0.5;
          dRSq_i=dR_i*dR_i;
          dA_ip1half=dGrid[bin->nR][i+1][0][0]*dGrid[bin->nR][i+1][0][0];
          dA_im1half=dGrid[bin->nR][i][0][0]*dGrid[bin->nR][i][0][0];
          dC=sqrt(dGamma*dP/dGrid[bin->nD][i][0][0]);
          dDVDtThreshold=dAVThreshold*dC;
          dDVDt=(dA_ip1half*dGrid[bin->nU][i+1][0][0]
            -dA_im1half*dGrid[bin->nU][i][0][0])/dRSq_i;
          if(dDVDt<-1.0*dDVDtThreshold){//being compressed
            dDVDt_mthreshold=dDVDt+dDVDtThreshold;
            dQ=dASq*dGrid[bin->nD][i][0][0]*dDVDt_mthreshold*dDVDt_mthreshold;
          }
          else{
            dQ=0.0;
          }
          
          //print them out
          ofFile
            <<std::setw(nWidthOutputField)<<dGrid[bin->nE][i][0][0]<<std::setw(nWidthOutputField)<<0.0//6,7
            <<std::setw(nWidthOutputField)<<"-"<<std::setw(nWidthOutputField)<<"-"//8,9
            <<std::setw(nWidthOutputField)<<dP<<std::setw(nWidthOutputField)<<0.0//10,11
            <<std::setw(nWidthOutputField)<<dQ<<std::setw(nWidthOutputField)<<0.0//12,13
            <<std::setw(nWidthOutputField)<<"-"<<std::setw(nWidthOutputField)<<"-"//14,15
            <<std::setw(nWidthOutputField)<<dGamma<<std::setw(nWidthOutputField)<<0.0//16,17
            <<std::setw(nWidthOutputField)<<"-";//18
        }
        ofFile<<std::endl;
      }
      ofFile
        <<std::setw(nWidthOutputField)<<"-"//0
        <<std::setw(nWidthOutputField)<<"-"//1
        <<std::setw(nWidthOutputField)<<"-";//2
      if(bin->nW!=-1){
        ofFile
          <<std::setw(nWidthOutputField)<<0.0;//3
      }
      else{
        ofFile
          <<std::setw(nWidthOutputField)<<"-";//3
      }
      ofFile
        <<std::setw(nWidthOutputField)<<"-"<<std::setw(nWidthOutputField)//4
        <<"-"//5
        <<std::setw(nWidthOutputField)<<"-"<<std::setw(nWidthOutputField)//6
        <<"-"//7
        <<std::setw(nWidthOutputField)<<"-"<<std::setw(nWidthOutputField)//8
        <<"-"//9
        <<std::setw(nWidthOutputField)<<"-"<<std::setw(nWidthOutputField)//10
        <<"-"//11
        <<std::setw(nWidthOutputField)<<"-"<<std::setw(nWidthOutputField)//12
        <<"-"//13
        <<std::setw(nWidthOutputField)<<"-"<<std::setw(nWidthOutputField)//14
        <<"-"//15
        <<std::setw(nWidthOutputField)<<"-"<<std::setw(nWidthOutputField)//16
        <<"-"<<std::setw(nWidthOutputField)<<"-"//17,18
        <<std::endl;
    }
    
    //3D region
    int nSizeX2=nSize[bin->nD][0]+2*nNumGhostCells;
    for(int i=nSizeX1;i<nSizeX2;i++){
      
      dR_i=(dGrid[bin->nR][i+1][0][0]+dGrid[bin->nR][i][0][0])*0.5;
      dRSq_ip1half=dGrid[bin->nR][i+1][0][0]*dGrid[bin->nR][i+1][0][0];
      dRSq_i=dR_i*dR_i;
      dA_ip1half=dGrid[bin->nR][i+1][0][0]*dGrid[bin->nR][i+1][0][0];
      dA_im1half=dGrid[bin->nR][i][0][0]*dGrid[bin->nR][i][0][0];
      
      //calculate horizontal average of various quantities
      double dDAve=0.0;
      double dEAve=0.0;
      double dTAve=0.0;
      double dPAve=0.0;
      double dQAve=0.0;
      double dKappaAve=0.0;
      double dGammaAve=0.0;
      int nCount=0;
      
      //calculate some theta areas
      if(nNumDims>1){
        dTheta_jp1half=dGrid[bin->nTheta][0][nPlaneIndex][0];
        if(nPlaneIndex==0){
          dTheta_jm1half=dGrid[bin->nTheta][0][nPlaneIndex][0]-(dGrid[bin->nTheta][0][nPlaneIndex+1][0]
            -dGrid[bin->nTheta][0][nPlaneIndex][0]);
        }
        else{
          dTheta_jm1half=dGrid[bin->nTheta][0][nPlaneIndex-1][0];
        }
        dTheta_j=(dTheta_jp1half+dTheta_jm1half)*0.5;
        dA_jp1half=sin(dTheta_jp1half);
        dA_jm1half=sin(dTheta_jm1half);
        dA_j=sin(dTheta_j);
      }
      
      if(nGammaLaw!=0){//set P,E,kappa,gamma, Q, and L
        
        for(int k=0;k<nSizeZ;k++){
          
          //get P, E, kappa, and gamma
         eosTable->getPEKappaGammaCp(dGrid[bin->nT][i][nPlaneIndex][k],dGrid[bin->nD][i][nPlaneIndex][k],dP
            ,dE,dKappa,dGamma,dCp);
          
          //calculate Q
          dC=sqrt(dGamma*dP/dGrid[bin->nD][i][nPlaneIndex][k]);
          dDVDtThreshold=dAVThreshold*dC;
          
          //Q0
          if(nNumDims>=1){
            if(i==nSizeX1){
              dDVDt=(dA_ip1half*dGrid[bin->nU][i+1][nPlaneIndex][k]
                -dA_im1half*dGrid[bin->nU][i][0][0])/dRSq_i;
            }
            else{
              dDVDt=(dA_ip1half*dGrid[bin->nU][i+1][nPlaneIndex][k]
                -dA_im1half*dGrid[bin->nU][i][nPlaneIndex][k])/dRSq_i;
            }
            if(dDVDt<-1.0*dDVDtThreshold){//being compressed
              dDVDt_mthreshold=dDVDt+dDVDtThreshold;
              dQ0=dASq*dGrid[bin->nD][i][nPlaneIndex][k]*dDVDt_mthreshold*dDVDt_mthreshold;
            }
            else{
              dQ0=0.0;
            }
          }
          
          //Q1
          if(nNumDims>=2){
            if(nPlaneIndex==0){
              dDVDt=(dA_jp1half*dGrid[bin->nV][i][nPlaneIndex][k]
                -dA_jm1half*dGrid[bin->nV][i][nSizeY2-1][k])/dA_j;
            }
            else{
              dDVDt=(dA_jp1half*dGrid[bin->nV][i][nPlaneIndex][k]
                -dA_jm1half*dGrid[bin->nV][i][nPlaneIndex-1][k])/dA_j;
            }
            if(dDVDt<-1.0*dDVDtThreshold){//being compressed
              dDVDt_mthreshold=dDVDt+dDVDtThreshold;
              dQ1=dASq*dGrid[bin->nD][i][nPlaneIndex][k]
                *dDVDt_mthreshold*dDVDt_mthreshold;
            }
            else{
              dQ1=0.0;
            }
          }
          
          //Q2
          if(nNumDims==3){
            if(k==0){
              dDVDt=(dGrid[bin->nW][i][nPlaneIndex][k]-dGrid[bin->nW][i][nPlaneIndex][nSizeZ-1]);
            }
            else{
              dDVDt=(dGrid[bin->nW][i][nPlaneIndex][k]-dGrid[bin->nW][i][nPlaneIndex][k-1]);
            }
            if(dDVDt<-1.0*dDVDtThreshold){//being compressed
              dDVDt_mthreshold=dDVDt+dDVDtThreshold;
              dQ2=dASq*dGrid[bin->nD][i][nPlaneIndex][k]*dDVDt_mthreshold*dDVDt_mthreshold;
            }
            else{
              dQ2=0.0;
            }
          }
          dQ=dQ0+dQ1+dQ2;
          
          dDAve+=dGrid[bin->nD][i][nPlaneIndex][k];
          dEAve+=dE;
          dTAve+=dGrid[bin->nT][i][nPlaneIndex][k];
          dPAve+=dP;
          dQAve+=dQ;
          dKappaAve+=dKappa;
          dGammaAve+=dGamma;
          nCount++;
        }
        dDAve=dDAve/double(nCount);
        dEAve=dEAve/double(nCount);
        dTAve=dTAve/double(nCount);
        dPAve=dPAve/double(nCount);
        dQAve=dQAve/double(nCount);
        dKappaAve=dKappaAve/double(nCount);
        dGammaAve=dGammaAve/double(nCount);
      }
      else{//set P and Q
        
        for(int k=0;k<nSizeZ;k++){
          
          if(nNumDims>1){
            dTheta_jp1half=dGrid[bin->nTheta][0][nPlaneIndex][0];
            if(nPlaneIndex==0){
              dTheta_jm1half=dGrid[bin->nTheta][0][nPlaneIndex][0]-(dGrid[bin->nTheta][0][nPlaneIndex+1][0]
                -dGrid[bin->nTheta][0][nPlaneIndex][0]);
            }
            else{
              dTheta_jm1half=dGrid[bin->nTheta][0][nPlaneIndex-1][0];
            }
            dTheta_j=(dTheta_jp1half+dTheta_jm1half)*0.5;
            dA_jp1half=sin(dTheta_jp1half);
            dA_jm1half=sin(dTheta_jm1half);
            dA_j=sin(dTheta_j);
          }
          
          //get P
          dP=dGrid[bin->nD][i][nPlaneIndex][k]*(dGamma-1.0)*dGrid[bin->nE][i][nPlaneIndex][k];
          
          //calculate Q
          dC=sqrt(dGamma*dP/dGrid[bin->nD][i][nPlaneIndex][k]);
          dDVDtThreshold=dAVThreshold*dC;
          
          //Q0
          if(nNumDims>=1){
            if(i==nSizeX1){
              dDVDt=(dA_ip1half*dGrid[bin->nU][i+1][nPlaneIndex][k]
                -dA_im1half*dGrid[bin->nU][i][0][0])/dRSq_i;
            }
            else{
              dDVDt=(dA_ip1half*dGrid[bin->nU][i+1][nPlaneIndex][k]
                -dA_im1half*dGrid[bin->nU][i][nPlaneIndex][k])/dRSq_i;
            }
            if(dDVDt<-1.0*dDVDtThreshold){//being compressed
              dDVDt_mthreshold=dDVDt+dDVDtThreshold;
              dQ0=dASq*dGrid[bin->nD][i][nPlaneIndex][k]*dDVDt_mthreshold*dDVDt_mthreshold;
            }
            else{
              dQ0=0.0;
            }
          }
          
          //Q1
          if(nNumDims>=2){
            if(nPlaneIndex==0){
              dDVDt=(dA_jp1half*dGrid[bin->nV][i][nPlaneIndex][k]
                -dA_jm1half*dGrid[bin->nV][i][nSizeY2-1][k])/dA_j;
            }
            else{
              dDVDt=(dA_jp1half*dGrid[bin->nV][i][nPlaneIndex][k]
                -dA_jm1half*dGrid[bin->nV][i][nPlaneIndex-1][k])/dA_j;
            }
            if(dDVDt<-1.0*dDVDtThreshold){//being compressed
              dDVDt_mthreshold=dDVDt+dDVDtThreshold;
              dQ1=dASq*dGrid[bin->nD][i][nPlaneIndex][k]
                *dDVDt_mthreshold*dDVDt_mthreshold;
            }
            else{
              dQ1=0.0;
            }
          }
          
          //Q2
          if(nNumDims==3){
            if(k==0){
              dDVDt=(dGrid[bin->nW][i][nPlaneIndex][k]-dGrid[bin->nW][i][nPlaneIndex][nSizeZ-1]);
            }
            else{
              dDVDt=(dGrid[bin->nW][i][nPlaneIndex][k]-dGrid[bin->nW][i][nPlaneIndex][k-1]);
            }
            if(dDVDt<-1.0*dDVDtThreshold){//being compressed
              dDVDt_mthreshold=dDVDt+dDVDtThreshold;
              dQ2=dASq*dGrid[bin->nD][i][nPlaneIndex][k]*dDVDt_mthreshold*dDVDt_mthreshold;
            }
            else{
              dQ2=0.0;
            }
          }
          dQ=dQ0+dQ1+dQ2;
          
          dDAve+=dGrid[bin->nD][i][nPlaneIndex][k];
          dEAve+=dGrid[bin->nE][i][nPlaneIndex][k];
          dPAve+=dP;
          dQAve+=dQ;
          nCount++;
        }
        dDAve=dDAve/double(nCount);
        dEAve=dEAve/double(nCount);
        dPAve=dPAve/double(nCount);
        dQAve=dQAve/double(nCount);
      }
      
      for(int k=0;k<nSizeZ;k++){

        if(nNumDims>1){
          dTheta_jp1half=dGrid[bin->nTheta][0][nPlaneIndex][0];
          if(nPlaneIndex==0){
            dTheta_jm1half=dGrid[bin->nTheta][0][nPlaneIndex][0]-(dGrid[bin->nTheta][0][nPlaneIndex+1][0]
              -dGrid[bin->nTheta][0][nPlaneIndex][0]);
          }
          else{
            dTheta_jm1half=dGrid[bin->nTheta][0][nPlaneIndex-1][0];
          }
          dTheta_j=(dTheta_jp1half+dTheta_jm1half)*0.5;
          dA_jp1half=sin(dTheta_jp1half);
          dA_jm1half=sin(dTheta_jm1half);
          dA_j=sin(dTheta_j);
        }
        
        if(i==nSizeX1){//first zone inside 3D region, has 1D u_im1half
          ofFile<<std::setw(nWidthOutputField)<<dGrid[bin->nU][i][0][0];//0
        }
        else{
          ofFile<<std::setw(nWidthOutputField)<<dGrid[bin->nU][i][nPlaneIndex][k];//0
        }
        ofFile
          <<std::setw(nWidthOutputField)<<dGrid[bin->nU0][i][0][0]//1
          <<std::setw(nWidthOutputField)<<dGrid[bin->nV][i][nPlaneIndex][k];//2
        if(bin->nW!=-1){
          ofFile
            <<std::setw(nWidthOutputField)<<dGrid[bin->nW][i][nPlaneIndex][k];//3
        }
        else{
          ofFile
            <<std::setw(nWidthOutputField)<<"-";//3
        }
        ofFile
          <<std::setw(nWidthOutputField)<<dGrid[bin->nD][i][nPlaneIndex][k]//4
          <<std::setw(nWidthOutputField)<<(dGrid[bin->nD][i][nPlaneIndex][k]-dDAve)/dDAve;//5
        if(nGammaLaw!=0){
          //get P, E, kappa, and gamma,Cp
         eosTable->getPEKappaGammaCp(dGrid[bin->nT][i][nPlaneIndex][k],dGrid[bin->nD][i][nPlaneIndex][k],dP
            ,dE,dKappa,dGamma,dCp);
          
          //calculate Q
          dC=sqrt(dGamma*dP/dGrid[bin->nD][i][nPlaneIndex][k]);
          dDVDtThreshold=dAVThreshold*dC;
          
          //Q0
          if(nNumDims>=1){
            if(i==nSizeX1){
              dDVDt=(dA_ip1half*dGrid[bin->nU][i+1][nPlaneIndex][k]
                -dA_im1half*dGrid[bin->nU][i][0][0])/dRSq_i;
            }
            else{
              dDVDt=(dA_ip1half*dGrid[bin->nU][i+1][nPlaneIndex][k]
                -dA_im1half*dGrid[bin->nU][i][nPlaneIndex][k])/dRSq_i;
            }
            if(dDVDt<-1.0*dDVDtThreshold){//being compressed
              dDVDt_mthreshold=dDVDt+dDVDtThreshold;
              dQ0=dASq*dGrid[bin->nD][i][nPlaneIndex][k]*dDVDt_mthreshold*dDVDt_mthreshold;
            }
            else{
              dQ0=0.0;
            }
          }
          
          //Q1
          if(nNumDims>=2){
            if(nPlaneIndex==0){
              dDVDt=(dA_jp1half*dGrid[bin->nV][i][nPlaneIndex][k]
                -dA_jm1half*dGrid[bin->nV][i][nSizeY2-1][k])/dA_j;
            }
            else{
              dDVDt=(dA_jp1half*dGrid[bin->nV][i][nPlaneIndex][k]
                -dA_jm1half*dGrid[bin->nV][i][nPlaneIndex-1][k])/dA_j;
            }
            if(dDVDt<-1.0*dDVDtThreshold){//being compressed
              dDVDt_mthreshold=dDVDt+dDVDtThreshold;
              dQ1=dASq*dGrid[bin->nD][i][nPlaneIndex][k]
                *dDVDt_mthreshold*dDVDt_mthreshold;
            }
            else{
              dQ1=0.0;
            }
          }
          
          //Q2
          if(nNumDims==3){
            if(k==0){
              dDVDt=(dGrid[bin->nW][i][nPlaneIndex][k]-dGrid[bin->nW][i][nPlaneIndex][nSizeZ-1]);
            }
            else{
              dDVDt=(dGrid[bin->nW][i][nPlaneIndex][k]-dGrid[bin->nW][i][nPlaneIndex][k-1]);
            }
            if(dDVDt<-1.0*dDVDtThreshold){//being compressed
              dDVDt_mthreshold=dDVDt+dDVDtThreshold;
              dQ2=dASq*dGrid[bin->nD][i][nPlaneIndex][k]*dDVDt_mthreshold*dDVDt_mthreshold;
            }
            else{
              dQ2=0.0;
            }
          }
          dQ=dQ0+dQ1+dQ2;
          
          ofFile
            <<std::setw(nWidthOutputField)<<dE<<std::setw(nWidthOutputField)<<(dE-dEAve)/dEAve//6,7
            <<std::setw(nWidthOutputField)<<dGrid[bin->nT][i][nPlaneIndex][k]//8
            <<std::setw(nWidthOutputField)<<(dGrid[bin->nT][i][nPlaneIndex][k]-dTAve)/dTAve//9
            <<std::setw(nWidthOutputField)<<dP<<std::setw(nWidthOutputField)<<(dP-dPAve)/dPAve;//10,11
          if(dQAve==0.0){/*if QAve is zero, this can only be the case if all Q's are zero in the 
            horizontal zone also*/
            ofFile
              <<std::setw(nWidthOutputField)<<dQ<<std::setw(nWidthOutputField)<<0.0;//12,13
          }
          else{
            ofFile
              <<std::setw(nWidthOutputField)<<dQ<<std::setw(nWidthOutputField)<<(dQ-dQAve)/dQAve;//12,13
          }
          ofFile
            <<std::setw(nWidthOutputField)<<dKappa<<std::setw(nWidthOutputField)//14
            <<(dKappa-dKappaAve)/dKappaAve//15
            <<std::setw(nWidthOutputField)<<dGamma<<std::setw(nWidthOutputField)//16
            <<(dGamma-dGammaAve)/dGammaAve<<std::setw(nWidthOutputField)<<dCp<<std::endl;//17,18
        }
        else{
          
          //get P
          dP=dGrid[bin->nD][i][nPlaneIndex][k]*(dGamma-1.0)*dGrid[bin->nE][i][nPlaneIndex][k];
          
          //calculate Q
          dC=sqrt(dGamma*dP/dGrid[bin->nD][i][nPlaneIndex][k]);
          dDVDtThreshold=dAVThreshold*dC;
          
          //Q0
          if(nNumDims>=1){
            if(i==nSizeX1){
              dDVDt=(dA_ip1half*dGrid[bin->nU][i+1][nPlaneIndex][k]
                -dA_im1half*dGrid[bin->nU][i][0][0])/dRSq_i;
            }
            else{
              dDVDt=(dA_ip1half*dGrid[bin->nU][i+1][nPlaneIndex][k]
                -dA_im1half*dGrid[bin->nU][i][nPlaneIndex][k])/dRSq_i;
            }
            if(dDVDt<-1.0*dDVDtThreshold){//being compressed
              dDVDt_mthreshold=dDVDt+dDVDtThreshold;
              dQ0=dASq*dGrid[bin->nD][i][nPlaneIndex][k]*dDVDt_mthreshold*dDVDt_mthreshold;
            }
            else{
              dQ0=0.0;
            }
          }
          
          //Q1
          if(nNumDims>=2){
            if(nPlaneIndex==0){
              dDVDt=(dA_jp1half*dGrid[bin->nV][i][nPlaneIndex][k]
                -dA_jm1half*dGrid[bin->nV][i][nSizeY2-1][nPlaneIndex])/dA_j;
            }
            else{
              dDVDt=(dA_jp1half*dGrid[bin->nV][i][nPlaneIndex][k]
                -dA_jm1half*dGrid[bin->nV][i][nPlaneIndex-1][k])/dA_j;
            }
            if(dDVDt<-1.0*dDVDtThreshold){//being compressed
              dDVDt_mthreshold=dDVDt+dDVDtThreshold;
              dQ1=dASq*dGrid[bin->nD][i][nPlaneIndex][k]
                *dDVDt_mthreshold*dDVDt_mthreshold;
            }
            else{
              dQ1=0.0;
            }
          }
          
          //Q2
          if(nNumDims==3){
            if(k==0){
              dDVDt=(dGrid[bin->nW][i][nPlaneIndex][k]-dGrid[bin->nW][i][nPlaneIndex][nSizeZ-1]);
            }
            else{
              dDVDt=(dGrid[bin->nW][i][nPlaneIndex][k]-dGrid[bin->nW][i][nPlaneIndex][k-1]);
            }
            if(dDVDt<-1.0*dDVDtThreshold){//being compressed
              dDVDt_mthreshold=dDVDt+dDVDtThreshold;
              dQ2=dASq*dGrid[bin->nD][i][nPlaneIndex][k]*dDVDt_mthreshold*dDVDt_mthreshold;
            }
            else{
              dQ2=0.0;
            }
          }
          dQ=dQ0+dQ1+dQ2;
          
          ofFile
            <<std::setw(nWidthOutputField)<<dE<<std::setw(nWidthOutputField)<<(dE-dEAve)/dEAve//6,7
            <<std::setw(nWidthOutputField)<<"-"<<std::setw(nWidthOutputField)<<"-"//8,9
            <<std::setw(nWidthOutputField)<<dP<<std::setw(nWidthOutputField)<<(dP-dPAve)/dPAve;//10,11
          if(dQAve==0.0){/*if QAve is zero, this can only be the case if all Q's are zero in the 
            horizontal zone also*/
            ofFile
              <<std::setw(nWidthOutputField)<<dQ<<std::setw(nWidthOutputField)<<0.0;//12,13
          }
          else{
            ofFile
              <<std::setw(nWidthOutputField)<<dQ<<std::setw(nWidthOutputField)<<(dQ-dQAve)/dQAve;//12,13
          }
          ofFile
            <<std::setw(nWidthOutputField)<<"-"<<std::setw(nWidthOutputField)<<"-"//14,15
            <<std::setw(nWidthOutputField)<<"-"<<std::setw(nWidthOutputField)<<"-"//16,17
            <<std::setw(nWidthOutputField)<<"-"<<std::endl;//18
        }
      }
      ofFile
        <<std::setw(nWidthOutputField)<<"-"//0
        <<std::setw(nWidthOutputField)<<"-"//1
        <<std::setw(nWidthOutputField)<<"-";//2
      if(bin->nW!=-1){
        ofFile
          <<std::setw(nWidthOutputField)<<dGrid[bin->nW][i][nPlaneIndex][0];//3
      }
      else{
        ofFile
          <<std::setw(nWidthOutputField)<<"-";//3
      }
      ofFile
        <<std::setw(nWidthOutputField)<<"-"<<std::setw(nWidthOutputField)//4
        <<"-"//5
        <<std::setw(nWidthOutputField)<<"-"<<std::setw(nWidthOutputField)//6
        <<"-"//7
        <<std::setw(nWidthOutputField)<<"-"<<std::setw(nWidthOutputField)//8
        <<"-"//9
        <<std::setw(nWidthOutputField)<<"-"<<std::setw(nWidthOutputField)//10
        <<"-"//11
        <<std::setw(nWidthOutputField)<<"-"<<std::setw(nWidthOutputField)//12
        <<"-"//13
        <<std::setw(nWidthOutputField)<<"-"<<std::setw(nWidthOutputField)//14
        <<"-"//15
        <<std::setw(nWidthOutputField)<<"-"<<std::setw(nWidthOutputField)//16
        <<"-"<<std::setw(nWidthOutputField)<<"-"//17,18
        <<std::endl;
    }
    for(int k=0;k<nSizeZ;k++){
      ofFile
        <<std::setw(nWidthOutputField)<<dGrid[bin->nU][nSizeX2][nPlaneIndex][k]//0
        <<std::setw(nWidthOutputField)<<dGrid[bin->nU0][nSizeX2][0][0]//1
        <<std::setw(nWidthOutputField)<<"-"//2
        <<std::setw(nWidthOutputField)<<"-"//3
        <<std::setw(nWidthOutputField)<<"-"<<std::setw(nWidthOutputField)//4
        <<"-"//5
        <<std::setw(nWidthOutputField)<<"-"<<std::setw(nWidthOutputField)//6
        <<"-"//7
        <<std::setw(nWidthOutputField)<<"-"<<std::setw(nWidthOutputField)//8
        <<"-"//9
        <<std::setw(nWidthOutputField)<<"-"<<std::setw(nWidthOutputField)//10
        <<"-"//11
        <<std::setw(nWidthOutputField)<<"-"<<std::setw(nWidthOutputField)//12
        <<"-"//13
        <<std::setw(nWidthOutputField)<<"-"<<std::setw(nWidthOutputField)//14
        <<"-"//15
        <<std::setw(nWidthOutputField)<<"-"<<std::setw(nWidthOutputField)//16
        <<"-"<<std::setw(nWidthOutputField)<<"-"//17,18
        <<std::endl;
    }
  }
  
  ofFile.close();
  for(int n=0;n<nNumVars;n++){
    delete [] nSize[n];
    delete [] nVarInfo[n];
  }
  delete [] nSize;
  delete [] nVarInfo;
}