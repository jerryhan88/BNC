//
//  Etc.hpp
//  BNC
//
//  Created by Chung-Kyun HAN on 1/8/19.
//  Copyright © 2019 Chung-Kyun HAN. All rights reserved.
//

#ifndef Etc_hpp
#define Etc_hpp

#include <iostream>
#include <stdio.h>
#include <cstdlib>
#include <fstream>
#include <ctime>
#include <chrono>

#include <boost/filesystem.hpp>

namespace fs = boost::filesystem;

double get_random_number();

void createCSV(std::string fpath, char *header);
void appendRow(std::string fpath, char *row);
void appendRows(std::string fpath, char **rows, int numRows);

class FilePathOrganizer {
public:
    std::string logPath, solPathCSV, solPathTXT,
                lpPath, ilpPath;
    //
    FilePathOrganizer(const char * argv[]);
    FilePathOrganizer(std::string, std::string, std::string);
    FilePathOrganizer(const fs::path &appr_dpath, const std::string &postfix);
    ~FilePathOrganizer() {}
    //
    void appendLog(char *);
private:
    void initLogFile();
};

class TimeTracker {
public:
    std::clock_t c_start;
    std::chrono::high_resolution_clock::time_point w_start;
    //
    TimeTracker() {
        c_start = std::clock();
        w_start = std::chrono::high_resolution_clock::now();
    }
    ~TimeTracker() {}
    //
    double get_elipsedTimeCPU();
    double get_elipsedTimeWall();
};

#endif /* Etc_hpp */
