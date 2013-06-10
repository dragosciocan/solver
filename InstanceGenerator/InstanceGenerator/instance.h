//
//  instance.h
//  InstanceGenerator
//
//  Created by Dragos Ciocan on 2/5/13.
//  Copyright (c) 2013 Dragos Ciocan. All rights reserved.
//

#ifndef InstanceGenerator_instance_h
#define InstanceGenerator_instance_h

#include <ext/hash_map>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include "global_problem.h"

#include <ilcplex/ilocplex.h>
ILOSTLBEGIN

using namespace std;

namespace distributed_solver {
class Instance {
    int num_advertisers_;
    int num_impressions_;
    int num_slots_; // Number of ad slots per impression.
    long double bid_sparsity_; // Percentage of impressions advertiser bids for.
    int num_shards_;
    long double width_;
    long double epsilon_;
    vector<__gnu_cxx::hash_map<int, long double> > bids_matrix_; // Matrix of bids of advertisers for impressions.
    vector<__gnu_cxx::hash_map<int, long double> > transpose_bids_matrix_;
    
    vector<long double> budgets_;
    vector<long double> slacks_;
    vector<long double> avg_slacks_;
    vector<__gnu_cxx::hash_map<int, long double> >* primal_sol_;
    vector<__gnu_cxx::hash_map<int, long double> > avg_primal_sol_;
    
    // Multiplicative weights related vars.
    int iteration_count_;
    vector<long double> weights_;
    GlobalProblem global_problem_;
    
public:
    Instance(int num_advertisers, int num_impressions, int num_slots, long double bid_sparsity);
    long double max_bid_;
    
    // Generation and output functions.
    void GenerateInstance();
    void WriteInstanceToCSV(std::string file_name_handle);
    void GenerateAndWriteInstance(std::string file_name_handle);
    void SetBudgets(long double factor);
    
    // Creates current global problem.
    void RunMultiplicativeWeights(long double num_iterations, long double numerical_accuracy_tolerance);
    void CreateGlobalProblem();
    void ComputeWeightedBudget();
    void UpdateGlobalProblem();
    void UpdateWeights();
    void CalculateInstanceWidth();
    void CalculateSlacks();
    void UpdateAvgPrimal(int t);
    void UpdateAvgSlacks(int t);
    void ReportWorstInfeasibility(int t);
    void BuildPrimals();
    void ReportWeightStats();
    void ComputeCPLEXRevenue();

 private:
    void VerifySolution();
    long double CalculateGlobalMWProblemOpt();
};
}

#endif
