
#include <string.h>
#include <vector>
#include <fstream>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <glob.h>

#include "binfile.h"
#include "exception2.h"

#ifdef _WIN32
    #include <Windows.h>
#elif defined __linux
    #include <sstream>
    #include <unistd.h>
#elif defined __APPLE__
    #include <mach-o/dyld.h>
#endif

#define MAX_PATH_LEN 2048

static std::string getExeDir_WIN32();
static std::string getExeDir_linux();
static std::string getExeDir_APPLE();

BinaryFile::BinaryFile(){
    setRootDir();
}

BinaryFile::BinaryFile(std::string filename){
    sFileName = filename;
    setRootDir();
}

bool BinaryFile::bFileExists(std::string strFilename){
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

void BinaryFile::setRootDir(){
    std::string exePath;
    #ifdef _WIN32
        exePath = getExeDir_WIN32();
    #elif defined __linux__
        exePath = getExeDir_linux();
    #elif defined __APPLE__
        exePath = getExeDir_APPLE();
    #else
        std::cout << "Unable to detect executable path on this platform. Put any nessecary files in working directory to continue." << std::endl;
        sRootDir=".";
        return;
    #endif
    //find the first "/" from the end
    unsigned pos=exePath.find_last_of("/");
    
    //keep from the begging to the location of the last "/" to remove the name
    //of the executable
    exePath=exePath.substr(0,pos);
    
    //check to see if the last directory is "bin" if so remove that also
    //as installed versions put the exe's into the bin directory and exePath
    //should point the top level directory.
    pos=exePath.find_last_of("/");
    std::string sBin=exePath.substr(pos+1,3);
    
    //if installed remove bin directory
    if(sBin.compare("bin")==0){
        exePath=exePath.substr(0,pos);
    }
    sRootDir = exePath;
}

static std::string getExeDir_WIN32(){
    #ifdef _WIN32
        char* exePath[2048];
        GetModuleFileName(NULL, exePath, 2048);
        std::string sExePath(exePath);
        return sExePath;
    #endif
    return NULL;
}

static std::string getExeDir_linux(){
    #ifdef __linux__
        char buff[MAX_PATH_LEN];
        std::string exePath;
        ssize_t len = readlink("/proc/self/exe", buff, sizeof(buff)-1);
        if (len != -1) {
            buff[len] = '\0';
            exePath=std::string(buff);
        } else {
            std::stringstream ssTemp;
            ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
            <<": error determining executable path"<<std::endl;
            throw exception2(ssTemp.str(),OUTPUT);
        }
        return exePath;
    #endif
    return NULL;
}

static std::string getExeDir_APPLE(){
    #ifdef __APPLE__
        char path[MAX_PATH_LEN];
        uint32_t pathLen = MAX_PATH_LEN;
        _NSGetExecutablePath(path, &pathLen);
        return std::string(path);
    #endif
    return NULL;
}

