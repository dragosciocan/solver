 //
//  instance.cpp
//  InstanceGenerator
//
//  Created by Dragos Ciocan on 2/5/13.
//  Copyright (c) 2013 Dragos Ciocan. All rights reserved.
//

#include <cmath>
#include "convex_hull.h"
#include "instance.h"
#include "naive_mw.h"

namespace distributed_solver {
    typedef IloArray<IloNumArray>    NumMatrix;
    typedef IloArray<IloNumVarArray> NumVarMatrix;
    
    Instance::Instance(int num_advertisers, int num_impressions, int num_slots, long double bid_sparsity,
                       long double epsilon, long double scaling_factor, long double numerical_accuracy_tolerance) {
        epsilon_ = epsilon;
        scaling_factor_ = scaling_factor;
        iteration_count_ = 0;
        numerical_accuracy_tolerance_ = numerical_accuracy_tolerance;
        
        num_advertisers_ = num_advertisers;
        num_impressions_ = num_impressions;
        num_slots_ = num_slots;
        bid_sparsity_ = bid_sparsity;
        num_shards_ = 10;
        
        budgets_ = vector<long double>(num_advertisers_);
        SetBudgets();
    }
    
    void Instance::GenerateInstance() {        
        transpose_bids_matrix_.reserve(num_impressions_);
        for (int i = 0; i < num_impressions_; ++i) {
            transpose_bids_matrix_.push_back(*new __gnu_cxx::hash_map<int, long double>());
        }
        
        srand(1);
        max_bid_ = 0;
        for (int j = 0; j < num_advertisers_; ++j) {
            // Generate bids for advertiser j.
            __gnu_cxx::hash_map<int, long double> bid_row;
            for (int i = 0; i < (bid_sparsity_ * num_impressions_); ++i) {
                int index = rand() % num_impressions_;
                long double bid = (long double) (rand() + 1) / ((long double) RAND_MAX);
                if (max_bid_ < bid) {
                    max_bid_ = bid;
                }
                bid_row[index] = bid;
                transpose_bids_matrix_[index][j] = bid;
            }
            bids_matrix_.push_back(bid_row);
        }
        cout << "Generated instance \n";
    }
    
    void Instance::WriteInstanceToCSV(string file_name_handle) {
        // Open file.
        string file_name;
        // cout << file_name_handle + "\n";
        // file_name = file_name_handle + to_string(num_advertisers_) + "x" + to_string(num_impressions_) + "x" + to_string(num_slots_) + "x" + to_string(bid_sparsity_) + ".csv";
        // __gnu_cxx::ofstream file;
        // cout << "Writing instance " + file_name + "\n";
        // file.open(file_name);
        // Go through bids and write them to file shards.
        string bid_row;
        
        for (int j = 0; j < num_advertisers_; ++j) {
            // Get bid vector for each advertiser, and write it to a row of the csv.
            /*
            for (__gnu_cxx::hash_map<int, long double>::iterator iter = bids_matrix_[j].begin();
                 iter != bids_matrix_[j].end();
                 ++iter) {
                bid_row = bid_row + to_string(iter->first) + "," + to_string(iter->second) + ",";
            }
             */
            bid_row = bid_row + "\n";
            // file << bid_row;
            bid_row = "";
        }
        
        // Close file
        // file.close();
    }
    
    void Instance::GenerateAndWriteInstance(string file_name_handle) {
        srand(1);
        string file_name;
        int num_advertisers_per_shard = num_advertisers_ / num_shards_;
        string bid_row;
        for (int k = 0; k < num_shards_; ++k) {
            // file_name = file_name_handle + to_string(num_advertisers_) + "x" + to_string(num_impressions_) + "x" + to_string(num_slots_) + "x" + to_string(bid_sparsity_)  + ".csv@" + to_string(k);
            // __gnu_cxx::ofstream file;
            cout << "Writing instance " + file_name + "\n";
            // file.open(file_name);
            for (int j = 0; j < num_advertisers_per_shard; ++j) {
                // Generate bids for advertiser j.
                __gnu_cxx::hash_map<int, long double> bid_map;
                for (int i = 0; i < (bid_sparsity_ * num_impressions_); ++i) {
                    int index = rand() % num_impressions_;
                    long double bid = (long double) (rand() + 1) / ((long double) RAND_MAX);
                    bid_map[index] = bid;
                }
                for (__gnu_cxx::hash_map<int, long double>::iterator iter = bid_map.begin(); iter != bid_map.end(); ++iter) {
                    // bid_row = bid_row + to_string(iter->first) + "," + to_string(iter->second) + ",";
                }
                bid_row = bid_row + "\n";
                // file << bid_row;
                bid_row = "";
            }
            // file.close();
        }
    }
    
    void Instance::RunMultiplicativeWeights(long double num_iterations, long double numerical_accuracy_tolerance) {
        BuildPrimals();
        
        bool mw_algorithm = true;
        
        if (mw_algorithm) {
            AllocationMW alloc_mw = AllocationMW(num_advertisers_, num_impressions_, num_slots_,
                                                 bid_sparsity_, max_bid_, epsilon_, numerical_accuracy_tolerance_,
                                                 primal_sol_, avg_primal_sol_,
                                                 &bids_matrix_, &transpose_bids_matrix_,
                                                 &budgets_, solution_);
            alloc_mw.RunAllocationMW(num_iterations);
        } else {
            NaiveMW naive_mw = NaiveMW(num_impressions_, num_advertisers_,
                                       max_bid_, epsilon_, width_, bid_sparsity_,
                                       primal_sol_, avg_primal_sol_,
                                       &bids_matrix_, &transpose_bids_matrix_,
                                       &budgets_);
            naive_mw.RunNaiveMW(num_iterations);
        }
    }
    
    void Instance::SetBudgets() {
        long double average_bid = 0.5; // Need to change this manually depending on how bids are drawn.
        for (int j = 0; j < num_advertisers_; ++j) {
            budgets_[j] = average_bid * (num_impressions_ / num_advertisers_) * scaling_factor_;
        }
    }
    
    void Instance::UpdateAvgPrimal(int t,
                                   vector<__gnu_cxx::hash_map<int, long double> >* sol,
                                   vector<__gnu_cxx::hash_map<int, long double> >* avg_sol) {
        long double new_value = 0;
        for (int i = 0; i < (*avg_sol).size(); ++i) {
            for (__gnu_cxx::hash_map<int, long double>::iterator iter = (*avg_sol)[i].begin();
                 iter != (*avg_sol)[i].end();
                 ++iter) {
                if ((*sol)[i].find(iter->first) != (*sol)[i].end()) {
                    new_value = (*sol)[i].find(iter->first)->second;
                } else {
                    new_value = 0;
                }
                iter->second = (long double)(t - 1) / t * iter->second + (long double)1 / t * new_value;
            }
        }
    }
    
    void Instance::BuildPrimals() {
        primal_sol_ = new vector<__gnu_cxx::hash_map<int, long double> >();
        avg_primal_sol_ = new vector<__gnu_cxx::hash_map<int, long double> >();
        for (int j = 0; j < num_advertisers_; ++j) {
            __gnu_cxx::hash_map<int, long double> row;
            
            __gnu_cxx::hash_map<int, pair<long double, long double> > row_2;
            
            for (__gnu_cxx::hash_map<int, long double>::iterator iter = bids_matrix_[j].begin();
                 iter != bids_matrix_[j].end();
                 ++iter) {
                row[iter->first] = 0;
                row_2[iter->first] = make_pair(0.0, 0.0);
            }
            primal_sol_->push_back(row);
            avg_primal_sol_->push_back(row);
            solution_->push_back(row_2);
        }
    }
    
    void Instance::ResetPrimal(vector<__gnu_cxx::hash_map<int, long double> >* sol) {
        for (int i = 0; i < (*sol).size(); ++i) {
            for (__gnu_cxx::hash_map<int, long double>::iterator iter = (*sol)[i].begin(); iter != (*sol)[i].end(); ++iter) {
                iter->second = 0.0;
            }
        }
    }
    
    void Instance::ResetCurrentPrimal(vector<__gnu_cxx::hash_map<int, pair<long double, long double> > >* sol) {
        for (int i = 0; i < (*sol).size(); ++i) {
            for (__gnu_cxx::hash_map<int, pair<long double, long double> >::iterator iter =
                 (*sol)[i].begin();
                 iter != (*sol)[i].end();
                 ++iter) {
                iter->second.first = 0.0;
            }
        }

    }
}
