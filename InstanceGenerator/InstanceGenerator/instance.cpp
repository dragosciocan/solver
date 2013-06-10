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
    
    Instance::Instance(int num_advertisers, int num_impressions, int num_slots, long double bid_sparsity) :
    global_problem_(num_impressions, max_bid_, bid_sparsity_ * num_impressions_) {
        epsilon_ = 0.5;
        iteration_count_ = 0;
        
        num_advertisers_ = num_advertisers;
        num_impressions_ = num_impressions;
        num_slots_ = num_slots;
        bid_sparsity_ = bid_sparsity;
        num_shards_ = 10;
        
        budgets_ = vector<long double>(num_advertisers_);
        weights_ = vector<long double>(num_advertisers_);
        slacks_ = vector<long double>(num_advertisers_);
        avg_slacks_ = vector<long double>(num_advertisers_);
        for (int i = 0; i < num_advertisers_; ++i) {
            weights_[i] = 1;
            // Perturbation against degeneracy
            // weights_[i] = 1 + ((long double) (rand() + 1) / ((long double) RAND_MAX)) / 10000;
            // cout << "Setting weight = ";
            // cout << weights_[i];
            // cout << "\n";
            slacks_[i] = 0;
        }
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
        SetBudgets(0.25);
        CalculateInstanceWidth();
        BuildPrimals();
        
        //ComputeCPLEXRevenue();
        
        bool mw_algorithm = false;
        
        if (mw_algorithm) {
        CreateGlobalProblem();
        global_problem_.InitializeBudgetAllocation();
            
        for (int t = 1; t <= num_iterations; ++t) {
            cout << "Entering iteration ";
            cout << t;
            cout << "\n";
                        
            // Get primal solution.
            global_problem_.ConstructPrimal(primal_sol_, t);
                
            // VerifySolution();
            
            // Calculate slacks, update averages and recalculate weights.
            cout << "Running MW update \n";
            CalculateSlacks();
            UpdateAvgPrimal(t);
            UpdateAvgSlacks(t);
            ReportWorstInfeasibility(t);
            UpdateWeights();
            UpdateGlobalProblem();
            ReportWeightStats();
        }
        } else {
            NaiveMW naive_mw = NaiveMW(num_impressions_, num_advertisers_, max_bid_, epsilon_, width_, bid_sparsity_,
                                       budgets_);
            naive_mw.RunNaiveMW(&bids_matrix_, primal_sol_, &avg_primal_sol_);
        }
    }
    
    void Instance::ComputeWeightedBudget() {
        // Compute budget.
        global_problem_.budget_ = 0;
        for (int i = 0; i < num_advertisers_; ++i) {
            global_problem_.budget_ = global_problem_.budget_ + weights_[i] * budgets_[i];
        }
    }
    
    void Instance::CreateGlobalProblem() {
        ComputeWeightedBudget();
        for (int i = 0; i < global_problem_.num_partitions_; ++i) {
            // Construct subproblems.
            vector<pair<long double, long double> >* coefficients;
            coefficients = new vector<pair<long double, long double> >();
            vector<int>* advertiser_index = new vector<int>();
            
            // Go through bids matrix and identify all bids for impression i.
            int subproblem_size = 0;
            /*for (int j = 0; j < num_advertisers_; ++j) {
                    if (bids_matrix_[j].find(i) != bids_matrix_[j].end()) {
                        coefficients->push_back(make_pair(bids_matrix_[j].find(i)->second,
                                                               bids_matrix_[j].find(i)->second * weights_[j]));
                    advertiser_index->push_back(j);
                    ++subproblem_size;
                }
            }*/
            
            for (__gnu_cxx::hash_map<int, long double>::iterator iter = transpose_bids_matrix_[i].begin();
                 iter != transpose_bids_matrix_[i].end();
                 ++iter) {
                coefficients->push_back(make_pair(iter->second,
                                                       weights_[iter->first]));
                advertiser_index->push_back(iter->first);
                ++subproblem_size;
            }
            
            global_problem_.subproblems_.push_back(Subproblem(subproblem_size, coefficients, advertiser_index));
        }
    }
    
    void Instance::UpdateGlobalProblem() {
        // Update weights.
        for (int i = 0; i < global_problem_.num_partitions_; ++i) {
            for (int j = 0; j < global_problem_.subproblems_[i].num_vars_; ++j) {
                global_problem_.subproblems_[i].constraints_[j].coefficient_ = global_problem_.subproblems_[i].constraints_[j].price_ * weights_[global_problem_.subproblems_[i].advertiser_index_->at(j)];
                global_problem_.subproblems_[i].constraints_[j].weight_ = weights_[global_problem_.subproblems_[i].advertiser_index_->at(j)];
                global_problem_.subproblems_[i].constraints_[j].is_active_ = true;
            }
        }
        ComputeWeightedBudget();
    }
    
    void Instance::CalculateInstanceWidth() {
        width_ = max_bid_ * (num_impressions_ * bid_sparsity_);
        for (int i = 0; i < num_advertisers_; ++i) {
            if (budgets_[i] > width_) {
                width_ = budgets_[i];
            }
        }
    }
    
    void Instance::CalculateSlacks() {
        for (int j = 0; j < num_advertisers_; ++j) {
            slacks_[j] = (-1) * budgets_[j];
            for (__gnu_cxx::hash_map<int, long double>::iterator iter = (*primal_sol_)[j].begin(); iter != (*primal_sol_)[j].end(); ++iter) {
                slacks_[j] = slacks_[j] + iter->second * bids_matrix_[j].find(iter->first)->second;
            }
        }
    }
    
    void Instance::UpdateWeights() {
        for (int j = 0; j < num_advertisers_; ++j) {
            long double tmp = slacks_[j]/width_;
            if (abs(tmp) > 1) {
                cout << "Slack normalization error\n";
            }
            if (tmp >= 0) {
                weights_[j] = weights_[j] * pow((1 + epsilon_), tmp);
            }
            else {
                weights_[j] = weights_[j] * pow((1 - epsilon_), -tmp);
            }
            // cout << "weight " << j << " set to " << weights_[j] << "\n";
        }
    }
    
    void Instance::SetBudgets(long double scaling_factor) {
        long double average_bid = 0.5; // Need to change this manually depending on how bids are drawn.
        for (int j = 0; j < num_advertisers_; ++j) {
            budgets_[j] = average_bid * (num_impressions_ / num_advertisers_) * scaling_factor;
        }
    }
    
    void Instance::UpdateAvgPrimal(int t) {
        long double new_value = 0;
        for (int i = 0; i < num_advertisers_; ++i) {
            for (__gnu_cxx::hash_map<int, long double>::iterator iter = avg_primal_sol_[i].begin(); iter != avg_primal_sol_[i].end(); ++iter) {
                if ((*primal_sol_)[i].find(iter->first) != (*primal_sol_)[i].end()) {
                    new_value = (*primal_sol_)[i].find(iter->first)->second;
                } else {
                    new_value = 0;
                }
                iter->second = (long double)(t - 1) / t * iter->second + (long double)1 / t * new_value;
            }
        }
    }
    
    void Instance::UpdateAvgSlacks(int t) {
        for (int i = 0; i < num_advertisers_; ++i) {
            avg_slacks_[i] = (long double) (t - 1) / t * avg_slacks_[i] + (long double) 1 / t * slacks_[i];
        }
    }
    
    void Instance::ReportWorstInfeasibility(int t) {
        long double max_infeasibility = 0.0;
        int max_infeasibility_index = -1;
        for (int i = 0; i < num_advertisers_; ++i) {
            if ((avg_slacks_[i] > 0) and ((avg_slacks_[i] / budgets_[i]) > max_infeasibility)){
                max_infeasibility = avg_slacks_[i] / budgets_[i];
                max_infeasibility_index = i;
            }
        }
        cout << "At iteration ";
        cout << t;
        cout << ", max infeasiblity was ";
        cout << max_infeasibility;
        cout << " on constraint ";
        cout << max_infeasibility_index;
        cout << "\n";
        
        ofstream report_file;
        report_file.open("/Users/ciocan/Documents/Google/data/infeasibility_report.csv", ios::app);
        report_file << t;
        report_file << ", ";
        report_file << max_infeasibility;
        report_file << ", ";
        report_file << max_infeasibility_index;
        report_file << "\n";
        report_file.close();
    }
    
    void Instance::BuildPrimals() {
        primal_sol_ = new vector<__gnu_cxx::hash_map<int, long double> >();
        for (int j = 0; j < num_advertisers_; ++j) {
            __gnu_cxx::hash_map<int, long double> row;
            for (__gnu_cxx::hash_map<int, long double>::iterator iter = bids_matrix_[j].begin(); iter != bids_matrix_[j].end(); ++iter) {
                row[iter->first] = 0;
            }
            // cout << j;
            primal_sol_->push_back(row);
            avg_primal_sol_.push_back(row);
        }
    
    }
    
    void Instance::VerifySolution() {
        IloEnv env;
        IloInt i, j;
        IloModel model(env);
        IloInt nbImpressions = num_impressions_;
        IloInt nbAdvertisers = num_advertisers_;
        
        NumVarMatrix x(env, nbImpressions);
        // NumVarMatrix y(env, nbImpressions);
        
        for(i = 0; i < nbImpressions; i++) {
            x[i] = IloNumVarArray(env, nbAdvertisers, 0.0, 1.0, ILOFLOAT);
            // y[i] = IloNumVarArray(env, nbAdvertisers, 0.0, 1.0, ILOFLOAT);
        }
        
        // Add impression constraints.
        for(i = 0; i < nbImpressions; i++) {
            model.add(IloSum(x[i]) <= 1.0);
        }
        
        // Add weighted contraint.
        IloExpr weighted_constraint(env);
        for(j = 0; j < nbImpressions; j++) {      // demand must meet supply
            for(i = 0; i < nbAdvertisers; i++) {
                weighted_constraint += ((double) transpose_bids_matrix_[j][i]) * ((double) weights_[i]) * x[j][i];
            }
        }
        model.add(weighted_constraint <= ((double) global_problem_.budget_));
        weighted_constraint.end();
        
        IloExpr obj_exp(env);
        for(i = 0; i < nbImpressions; i++) {
            for(j = 0; j < nbAdvertisers; j++) {
                obj_exp += ((double) transpose_bids_matrix_[i][j]) * x[i][j];
            }
        }
        
        model.add(IloMaximize(env, obj_exp));
        obj_exp.end();
    
        IloCplex cplex(env);
        cplex.setOut(env.getNullStream());
        cplex.extract(model);
        
        // Optimize the problem and obtain solution.
        if ( !cplex.solve() ) {
            env.error() << "Failed to optimize LP" << endl;
            throw(-1);
        }
        
        IloNumArray vals(env);
        
        long double sum_b = 0;
        for (int a = 0; a < num_advertisers_; ++a) {
            //slacks_[a] = 0;
            for (i = 0; i < nbImpressions; i++) {
                if (cplex.getValue(x[i][a]) > 0) {
                    //slacks_[a] = slacks_[a] + cplex.getValue(x[i][a]) * bids_matrix_[a].find(i)->second;
                    sum_b += cplex.getValue(x[i][a]) * bids_matrix_[a].find(i)->second * weights_[a];
                }
            }
        }
        
        cout << "Cplex buget allocation is = ";
        cout << sum_b;
        cout << "\n";
        
        if (cplex.getObjValue() == CalculateGlobalMWProblemOpt()) {
            cout << "Solution checks \n";
        } else {
            cout << "Solution does not check, Cplex opt is ";
            cout << cplex.getObjValue();
            cout << "\n";
        }
        env.end();
    }
    
    long double Instance::CalculateGlobalMWProblemOpt() {
        long double MW_revenue = 0;
        for (int a = 0; a < (*primal_sol_).size(); ++a) {
            for (__gnu_cxx::hash_map<int, long double>::iterator iter = (*primal_sol_)[a].begin();
                 iter != (*primal_sol_)[a].end();
                 ++iter) {
                MW_revenue += iter->second * bids_matrix_[a][iter->first];
            }
        }
        cout << "MW rev ";
        cout << MW_revenue;
        cout << "\n";
        return MW_revenue;
    }
    
    void Instance::ReportWeightStats() {
        long double min_weight = 100000;
        long double max_weight = 0;
        for (int i = 0; i < num_advertisers_; ++i) {
            if (weights_[i] < min_weight) {
                min_weight = weights_[i];
            }
            if (weights_[i] > max_weight) {
                max_weight = weights_[i];
            }
        }
        cout << "min weight = ";
        cout << min_weight;
        cout << ", max weight = ";
        cout << max_weight;
        cout << "\n";
    }
    
    void Instance::ComputeCPLEXRevenue() {
        IloEnv env;
        IloInt i, j;
        IloModel model(env);
        IloInt nbImpressions = num_impressions_;
        IloInt nbAdvertisers = num_advertisers_;
        
        NumVarMatrix x(env, nbImpressions);
        
        for(i = 0; i < nbImpressions; i++) {
            x[i] = IloNumVarArray(env, nbAdvertisers, 0.0, 1.0, ILOFLOAT);
        }
        
        // Add impression constraints.
        for(i = 0; i < nbImpressions; i++) {
            model.add(IloSum(x[i]) <= 1.0);
        }
        
        // Add weighted contraint.
        for(j = 0; j < nbAdvertisers; j++) {
            IloExpr curr_adv_constraint(env);
            for (__gnu_cxx::hash_map<int, long double>::iterator iter = bids_matrix_[j].begin();
                 iter != bids_matrix_[j].end();
                 ++iter) {
                curr_adv_constraint += ((double) iter->second) * ((double) weights_[j]) * x[iter->first][j];
            }
            model.add(curr_adv_constraint <= ((double) budgets_[j]));
        }
        
        IloExpr obj_exp(env);
        for(i = 0; i < nbImpressions; i++) {
            for(j = 0; j < nbAdvertisers; j++) {
                obj_exp += ((double) transpose_bids_matrix_[i][j]) * x[i][j];
            }
        }
        
        model.add(IloMaximize(env, obj_exp));
        obj_exp.end();
        
        IloCplex cplex(env);
        cplex.setOut(env.getNullStream());
        cplex.extract(model);
        
        // Optimize the problem and obtain solution.
        if ( !cplex.solve() ) {
            env.error() << "Failed to optimize LP" << endl;
            throw(-1);
        }
        
        cplex.exportModel("/Users/ciocan/Documents/Google/data/example.lp");
        
        cout << "CPLEX opt is " << cplex.getObjValue() << "\n";
        
        env.end();    
    }
}
