#include "fileExists.h"

//' Check if a file exists, if it does, do nothing, otherwise throw an exception
//'
//' @param name Filename / file path
//' @export
void fileExists (const std::string& name) {
    if ( access( name.c_str(), F_OK ) == -1 ) {
        std::ostringstream stringStream;
        stringStream << "File does not exist: " << name;
        throw std::string(stringStream.str());
    }
}

bool fileExistsBool(const std::string& name) {
    if (FILE *file = fopen(name.c_str(), "r")) {
        fclose(file);
        return true;
    } else {
        return false;
    }
}


void fileExists(const std::string& name, bool exitIfNotExist, bool exitIfExist) {
    bool existStatus = fileExistsBool(name);
    if(not existStatus) {
        if(exitIfNotExist) {
            std::ostringstream stringStream;
            stringStream << "File does not exist: " << name;
            throw std::string(stringStream.str());
        }
    }
    else {
        if(exitIfExist) {
            std::ostringstream stringStream;
            stringStream << "File already exists: " << name;
            throw std::string(stringStream.str());
        }
    }
}

bool fileExists(const std::string& name, bool rmIfExist) {
    bool existStatus = fileExistsBool(name);
    if(existStatus) {
        if(rmIfExist) {
            if( remove(name.c_str()) != 0 ) {
                std::ostringstream stringStream;
                stringStream << "Failed to remove file:  " << name;
                throw std::string(stringStream.str());
            }
        }
    }
    return existStatus;
}


