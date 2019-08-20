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
}

void appendRow(std::string fpath, char *row) {
    std::fstream fout;
    fout.open(fpath, std::ios::out | std::ios::app);
    fout << row << "\n";
}

void appendRows(std::string fpath, char **rows, int numRows) {
    std::fstream fout;
    fout.open(fpath, std::ios::out | std::ios::app);
    for (int i = 0; i < numRows; i++)
        fout << rows[i] << "\n";
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

void FilePathOrganizer::initLogFile() {
    char header[BUF_SIZE];
    sprintf(header, "nid,numIter,upperBound,SEC1,SEC2,note");
    createCSV(logPath, header);
}

void FilePathOrganizer::appendLog(char *row) {
    appendRow(logPath, row);
}

//void csvExample(const PathOrganizer &po) {
//    char header[BUF_SIZE];
//    sprintf(header, "c%d,c%d,c%d,c%d", 1, 2, 3, 4);
//
//    createCSV(po.solPathCSV, header);
//    char row[BUF_SIZE];
//    sprintf(row, "1, 2, 3, 4");
//    appendRow(po.solPathCSV, row);
//    char *rows[5];
//    for (int i = 0; i < 5; i++) {
//        rows[i] = row;
//    }
//    appendRows(po.solPathCSV, rows, 5);
//}


FilePathOrganizer::~FilePathOrganizer() {
    
}
