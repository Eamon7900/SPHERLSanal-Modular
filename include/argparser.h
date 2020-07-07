#include <string>

class ArgParser {
    private:
        int range_l;
        int range_u;
        std::string baseFileName;
    public:
        ArgParser(std::string fileRange);
        std::string getBaseFileName();
        std::vector<std::string> getFilesInRange();
};