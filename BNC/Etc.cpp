//
//  Etc.cpp
//  BNC
//
//  Created by Chung-Kyun HAN on 1/8/19.
//  Copyright © 2019 Chung-Kyun HAN. All rights reserved.
//

#include "Etc.hpp"

#define BUF_SIZE 1024


double get_random_number() {
    return ((double) rand()) / RAND_MAX;
}

void createCSV(std::string fpath, char *header) {
    std::fstream fout;
    fout.open(fpath, std::ios::out);
    fout << header << "\n";
    fout.close();
}

void appendRow(std::string fpath, char *row) {
    std::fstream fout;
    fout.open(fpath, std::ios::out | std::ios::app);
    fout << row << "\n";
    fout.close();
}

void appendRows(std::string fpath, char **rows, int numRows) {
    std::fstream fout;
    fout.open(fpath, std::ios::out | std::ios::app);
    for (int i = 0; i < numRows; i++)
        fout << rows[i] << "\n";
    fout.close();
}

FilePathOrganizer::FilePathOrganizer(const char * argv[]) {
    this->logPath = argv[1];
    this->solPathCSV = argv[2];
    this->solPathTXT = argv[3];
    initLogFile();
}

FilePathOrganizer::FilePathOrganizer(std::string logPath, std::string solPathCSV, std::string solPathTXT) {
    this->logPath = logPath;
    this->solPathCSV = solPathCSV;
    this->solPathTXT = solPathTXT;
}

FilePathOrganizer::FilePathOrganizer(const fs::path &appr_dpath, const std::string &postfix) {
    std::string appr_name = appr_dpath.leaf().string();
    //
    fs::path logPath(appr_dpath);
    logPath.append(appr_name + "_" + postfix + ".log");
    this->logPath = logPath.string();
    //
    fs::path solPathCSV(appr_dpath);
    solPathCSV.append(appr_name + "_" + postfix + ".csv");
    this->solPathCSV = solPathCSV.string();
    //
    fs::path solPathTXT(appr_dpath);
    solPathTXT.append(appr_name + "_" + postfix + ".txt");
    this->solPathTXT = solPathTXT.string();
    //
    fs::path lpPath(appr_dpath);
    lpPath.append(appr_name + "_" + postfix + ".lp");
    this->lpPath = lpPath.string();
    //
    fs::path ilpPath(appr_dpath);
    ilpPath.append(appr_name + "_" + postfix + ".ilp");
    this->ilpPath = ilpPath.string();
}

void FilePathOrganizer::initLogFile() {
    char header[BUF_SIZE];
    sprintf(header, "nid,numIter,upperBound,SEC0,SEC1,SEC2,note");
    createCSV(logPath, header);
}

void FilePathOrganizer::appendLog(char *row) {
    appendRow(logPath, row);
}


double TimeTracker::get_elipsedTimeCPU() {
    std::clock_t c_end = std::clock();
    return (c_end-c_start) / (double) CLOCKS_PER_SEC;
}

double TimeTracker::get_elipsedTimeWall() {
    std::chrono::high_resolution_clock::time_point w_end = std::chrono::high_resolution_clock::now();
    return std::chrono::duration<double, std::milli>(w_end-w_start).count() / 1000.0;
}
