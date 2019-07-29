//
//  main.cpp
//  BNC
//
//  Created by Chung-Kyun HAN on 20/7/19.
//  Copyright © 2019 Chung-Kyun HAN. All rights reserved.
//

//#include <cassert>
//#include <cstdlib>
//#include <cmath>
//#include <sstream>

//
#include <fstream>
#include "PDPTW_ILP.hpp"

#define BUF_SIZE 1024

void createCSV(string fpath, char *header) {
    fstream fout;
    fout.open(fpath, ios::out);
    fout << header << "\n";
}

void appendRow(string fpath, char *row) {
    fstream fout;
    fout.open(fpath, ios::out | ios::app);
    fout << row << "\n";
}


void appendRows(string fpath, char **rows, int numRows) {
    fstream fout;
    fout.open(fpath, ios::out | ios::app);
    for (int i = 0; i < numRows; i++)
        fout << rows[i] << "\n";
}

void csvExample() {
    char header[BUF_SIZE];
    sprintf(header, "c%d,c%d,c%d,c%d", 1, 2, 3, 4);
    
    string fpath = "/Users/ckhan/workspace/BNC/BNC/test.csv";
    createCSV(fpath, header);
    char row[BUF_SIZE];
    sprintf(row, "1, 2, 3, 4");
    appendRow(fpath, row);
    char *rows[5];
    for (int i = 0; i < 5; i++) {
        rows[i] = row;
    }
    appendRows(fpath, rows, 5);
}

int main(int argc, const char * argv[]) {
//    run_example();
    run_example_PDPTW();
    
//    csvExample();
    
    return 0;
}
