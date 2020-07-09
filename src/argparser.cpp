#include <string>
#include <vector>
#include <iostream>
#include <glob.h>
#include <limits.h>

#include "argparser.h"

ArgParser::ArgParser(std::string fileRange){
  std::string range;   
  for(int i = 0; i < fileRange.length(); i++){
    if(fileRange[i] == '['){
      fileRange = fileRange.substr(0, i); //Leave just the base filename
      ++i;
      while(fileRange[i] != ']'){
        range += fileRange[i++];  
      }
    }
    
  }

  if(range.length() == 0){ //If there isn't a '[' in the filerange, assume only one file specifed.
    baseFileName = fileRange;
    singleFile = true;
    return;
  } 

  std::string sLower;
  sLower.reserve(range.length() / 2);
  std::string sUpper;
  sUpper.reserve(range.length() / 2);

  //Split the range into two substrings:
  int i = 0;
  while(range[i] != '-'){
    sLower+=range[i];
    i++;
  }
  i++;
  for(;i < range.length(); i++){
    sUpper+=range[i];
  }
  
  //convert lower and upper bounds of range into integers
  range_l = stoi(sLower);
  if(sUpper == "*"){
    range_u = INT_MAX;
  }else {
    range_u= stoi(sUpper);
  }
  baseFileName = fileRange;
}

std::string ArgParser::getBaseFileName(){
  return baseFileName;
}

std::vector<std::string> ArgParser::getFilesInRange(){
  if(singleFile){
    std::vector<std::string> file;
    file.push_back(baseFileName);
    return file;
  }

  std::string rangeGlob = "[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]";
  std::string filesExistGlob = baseFileName.c_str() + rangeGlob;


  glob_t allFiles;
  glob(filesExistGlob.c_str(), NULL, NULL, &allFiles);


  std::vector<std::string> matches;
  for(int i = 0; i < allFiles.gl_pathc; i++){
    std::string fname = allFiles.gl_pathv[i];
    int fnum = atoi(fname.substr(baseFileName.length(), fname.length()).c_str());  
    if(range_l <= fnum && fnum <= range_u){
      matches.push_back(fname);
    }
  }
  return matches;
}