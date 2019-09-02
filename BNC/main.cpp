//
//  main.cpp
//  BNC
//
//  Created by Chung-Kyun HAN on 20/7/19.
//  Copyright © 2019 Chung-Kyun HAN. All rights reserved.
//

//
#include <iostream>

#include <nlohmann/json.hpp>
#include <boost/filesystem.hpp>
#include <map>
#include <vector>

#include "Problem.hpp"
#include "Etc.hpp"
#include "BnB/Tree.hpp"
#include "MathematicalModel.hpp"


#define BUF_SIZE 2048

namespace fs = boost::filesystem;

struct path_leaf_string
{
    std::string operator()(const fs::directory_entry& entry) const
    {
        return entry.path().leaf().string();
    }
};

std::vector<fs::path> get_probFiles(fs::path prob_dpath) {
    std::vector<fs::path> probFiles;
    std::vector<std::string> prob_fileNames;
    fs::directory_iterator start(prob_dpath);
    fs::directory_iterator end;
    std::transform(start, end, std::back_inserter(prob_fileNames), path_leaf_string());
    std::sort(prob_fileNames.begin(), prob_fileNames.end());
    for (std::string fn: prob_fileNames) {
        fs::path prob_fpath(prob_dpath);
        probFiles.push_back(prob_fpath.append(fn).string());
    }
    return probFiles;
}

void write_solution(Problem *prob, FilePathOrganizer &fpo, TimeTracker &tt,
                    MathematicalModel *mm,
                    double gap) {
    char row[BUF_SIZE];
    std::fstream fout_csv;
    fout_csv.open(fpo.solPathCSV, std::ios::out);
    fout_csv << "objV,Gap,eliCpuTime,eliWallTime" << "\n";
    sprintf(row, "%f,%f,%f,%f",
                mm->grbModel->get(GRB_DoubleAttr_ObjVal),
                gap,
                tt.get_elipsedTimeCPU(),
                tt.get_elipsedTimeWall());
    fout_csv << row << "\n";
    fout_csv.close();
    //
    double **x_ij = new double *[(*prob).N.size()];
    double *u_i = new double[(*prob).N.size()];
    for (int i: (*prob).N) {
        x_ij[i] = new double[(*prob).N.size()];
    }
    mm->get_x_ij(x_ij);
    mm->get_u_i(u_i);
    std::map<int, int> _route;
    for (int i: (*prob).N) {
        for (int j: (*prob).N) {
            if (x_ij[i][j] > 0.5) {
                _route[i] = j;
            }
        }
    }
    std::fstream fout_txt;
    fout_txt.open(fpo.solPathTXT, std::ios::out);
    std::vector<int> seq;
    int i = (*prob).o;
    fout_txt << "Visiting sequence" << "\n";
    
    while (true) {
        fout_txt << i << "-";
        seq.push_back(i);
        i = _route[i];
        if (i == (*prob).d) {
            fout_txt << i << "\n\n";
            break;
        }
    }
    seq.push_back(i);
    
    
    fout_txt << "Arrival time" << "\n";
    for (int i: seq) {
        fout_txt << i << ": " << u_i[i] << "\n";
    }
    fout_txt.close();
    //
    
    for (int i: (*prob).N) {
        delete [] x_ij[i];
    }
    delete [] x_ij;
    delete [] u_i;
}


int main(int argc, const char * argv[]) {
    
    
    fs::path prob_dpath(argv[1]);
    std::string appr_name(argv[2]);
    fs::path appr_dpath(prob_dpath);
    appr_dpath.remove_leaf();
    appr_dpath.append(appr_name);
    if (fs::status(appr_dpath).type() == fs::status_unknown ||
            fs::status(appr_dpath).type() == fs::file_not_found) {
        fs::create_directory(appr_dpath);
    }
    //
    std::vector<fs::path> probFiles = get_probFiles(prob_dpath);
    for (fs::path prob_fpath: probFiles) {
        Problem *prob = Problem::read_json(prob_fpath.string());
        std::string postfix = prob->problemName;
        FilePathOrganizer fpo(appr_dpath, postfix);
        TimeTracker tt;
        double gap;
        if (appr_name == "ILP") {
            MathematicalModel mm('I', prob);
            (*mm.grbModel).set(GRB_StringParam_LogFile, fpo.logPath);
            (*mm.grbModel).optimize();
            if ((*mm.grbModel).get(GRB_IntAttr_Status) == INFEASIBLE_MODEL) {
                (*mm.grbModel).computeIIS();
                (*mm.grbModel).write(fpo.ilpPath);
            }
             gap = (*mm.grbModel).get(GRB_DoubleAttr_MIPGap);
            write_solution(prob, fpo, tt, &mm, gap);
        } else {
            assert(appr_name == "BNC");
            
            
            MathematicalModel mm('I', prob);
            (*mm.grbModel).set(GRB_StringParam_LogFile, fpo.logPath);
            
            (*mm.grbModel).set(GRB_IntParam_PreCrush, 1);
            
            
            subtourelim cb = subtourelim(mm.x_ij, 12);
            (*mm.grbModel).setCallback(&cb);
            
            (*mm.grbModel).optimize();
            
            
            
//            Tree tree(prob, &fpo);
//            tree.run_BnB();
//            MathematicalModel *mm = tree.get_incumbentMM();
//            gap = 0.0;
//            write_solution(prob, fpo, tt, mm, gap);
        }
    }
    return 0;
}
