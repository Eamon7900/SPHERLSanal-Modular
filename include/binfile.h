#ifndef BINFILE_H
#define BINFILE_H

#include <string>

class BinaryFile {
    protected:
        BinaryFile(); //Default constructor for eos that is part of binary DataFile (has no file name)
        BinaryFile(std::string fileName);
    public:
        std::string sFileName; //The name of the binary file
        std::string sExeName; //The path of the executable which is operating on the binary file

};

#endif