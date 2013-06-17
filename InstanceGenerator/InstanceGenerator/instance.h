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
#include "allocation_mw.h"

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
    long double scaling_factor_;
    vector<__gnu_cxx::hash_map<int, long double> > bids_matrix_; // Matrix of bids of advertisers for impressions.
    vector<__gnu_cxx::hash_map<int, long double> > transpose_bids_matrix_;
    vector<long double> budgets_;
    vector<__gnu_cxx::hash_map<int, long double> >* primal_sol_;
    vector<__gnu_cxx::hash_map<int, long double> >* avg_primal_sol_;
    
    vector<__gnu_cxx::hash_map<int, pair<long double, long double> > >* solution_;
    
    // Multiplicative weights related vars.
    int iteration_count_;
    long double numerical_accuracy_tolerance_;
    
public:
    Instance(int num_advertisers, int num_impressions, int num_slots, long double bid_sparsity, long double epsilon,
             long double scaling_factor, long double numerical_accuracy_toleranc);
    long double max_bid_;
    
    // Generation and output functions.
    void GenerateInstance();
    void WriteInstanceToCSV(std::string file_name_handle);
    void GenerateAndWriteInstance(std::string file_name_handle);
    void SetBudgets();
    
    // Creates current global problem.
    void RunMultiplicativeWeights(long double num_iterations, long double numerical_accuracy_tolerance);
    static void UpdateAvgPrimal(int t,
                                vector<__gnu_cxx::hash_map<int, long double> >* sol,
                                vector<__gnu_cxx::hash_map<int, long double> >* avg_sol);
    void BuildPrimals();
    static void ResetPrimal(vector<__gnu_cxx::hash_map<int, long double> >* sol);
    static void ResetCurrentPrimal(vector<__gnu_cxx::hash_map<int, pair<long double, long double> > >* sol);
};
}

#endif
