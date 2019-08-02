//
//  Etc.hpp
//  BNC
//
//  Created by Chung-Kyun HAN on 1/8/19.
//  Copyright © 2019 Chung-Kyun HAN. All rights reserved.
//

#ifndef Etc_hpp
#define Etc_hpp

#include <stdio.h>
#include <cstdlib>
#include <fstream>


double get_random_number();

void createCSV(std::string fpath, char *header);
void appendRow(std::string fpath, char *row);
void appendRows(std::string fpath, char **rows, int numRows);

class FilePathOrganizer {
public:
    std::string logPath, solPathCSV, solPathTXT, tempPath;
    //
    FilePathOrganizer(const char * argv[]);
    FilePathOrganizer(std::string, std::string, std::string);
    ~FilePathOrganizer();
    //
    void appendLog(char *);
private:
    void initLogFile();
};


#endif /* Etc_hpp */
