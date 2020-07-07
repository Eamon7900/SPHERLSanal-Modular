#ifndef BINFILE_H
#define BINFILE_H

#include <string>
#include <vector>
#include <iostream>
#ifdef _WIN32
    #include <Windows.h>
#elif defined __linux
    #include <sstream>
    #include <unistd.h>
#elif defined __APPLE__
    #include <mach-o/dyld.h>
#endif
class BinaryFile {
    private:
        void setRootDir();
    protected:
        BinaryFile(); //Default constructor for eos that is part of binary DataFile (has no file name)
        BinaryFile(std::string fileName);
    public:
        std::string sFileName; //The name of the binary file
        std::string sRootDir; //The root installation directory of SPHERLSanal-Modular
        bool bFileExists(std::string strFilename);
};

#endif