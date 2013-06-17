//
//  naive_mw.cpp
//  InstanceGenerator
//
//  Created by Dragos Ciocan on 5/26/13.
//  Copyright (c) 2013 Dragos Ciocan. All rights reserved.
//

#include <algorithm>
#include <cmath>

#include "naive_mw.h"

namespace distributed_solver {
    NaiveMW::NaiveMW(int num_impressions, int num_advertisers, long double max_bid,
                     long double epsilon, long double width, long double bid_sparsity,
                    vector<__gnu_cxx::hash_map<int, long double> >* primal_sol,
                     vector<__gnu_cxx::hash_map<int, long double> >* avg_primal_sol,
                     vector<__gnu_cxx::hash_map<int, long double> >* bids_matrix,
                     vector<__gnu_cxx::hash_map<int, long double> >* transpose_bids_matrix,
                     vector<long double>* budgets) {
        num_impressions_ = num_impressions;
        num_advertisers_ = num_advertisers;
        max_bid_ = max_bid;
        epsilon_ = epsilon;
        bid_sparsity_ = bid_sparsity;
        width_ = fmax(width, bid_sparsity * num_advertisers_);
        
        for (int i = 0; i < num_impressions_; ++i) {
            impression_weights_.push_back(1.0);
            impression_slacks_.push_back(0.0);
            impression_avg_slacks_.push_back(0.0);
        }
        budgets_ = budgets;
        primal_sol_ = primal_sol;
        avg_primal_sol_ = avg_primal_sol;
        bids_matrix_ = bids_matrix;
        transpose_bids_matrix_ = transpose_bids_matrix;
        
        for (int i = 0; i < num_advertisers_; ++i) {
            advertiser_weights_.push_back(1.0);
            advertiser_slacks_.push_back(0.0);
            advertiser_avg_slacks_.push_back(0.0);
        }
    }
     
    void NaiveMW::CreateRatios() {
        ratios_.clear();
        int index = 0;
        for (int a = 0; a < (*bids_matrix_).size(); ++a) {
            for (__gnu_cxx::hash_map<int, long double>::iterator iter = (*bids_matrix_)[a].begin();
                 iter != (*bids_matrix_)[a].end(); ++iter) {
                ratios_.push_back(KnapsackRatio());
                ratios_[index].ratio = iter->second /
                                       (iter->second * advertiser_weights_[a] + impression_weights_[iter->first]);
                ratios_[index].coefficient = iter->second * advertiser_weights_[a] + impression_weights_[iter->first];
                ratios_[index].advertiser_index = a;
                ratios_[index].impression_index = iter->first;
                index++;
            }
        }
    }
    
    void NaiveMW::NaiveMWPrimal() {
        int adv_index = -1;
        int imp_index = -1;
        
        // Update ratios;
        for (int j = 0; j < ratios_.size(); ++j) {
            adv_index = ratios_[j].advertiser_index;
            imp_index = ratios_[j].impression_index;
            ratios_[j].ratio = (*bids_matrix_)[adv_index][imp_index] /
                               (advertiser_weights_[adv_index] * (*bids_matrix_)[adv_index][imp_index] +
                                impression_weights_[imp_index]);
            ratios_[j].coefficient = advertiser_weights_[adv_index] * (*bids_matrix_)[adv_index][imp_index] +
                                     impression_weights_[imp_index];
        }
        
        // Sort all knapsack ratios (price/coefficient).
        sort(ratios_.begin(), ratios_.end());
        
        // Allocate primal.
        long double remaining_budget = weighted_budget_;
        int current_index = 0;
        while ((remaining_budget > 0) && (current_index < ratios_.size())) {
            adv_index = ratios_[current_index].advertiser_index;
            imp_index = ratios_[current_index].impression_index;
            
            (*primal_sol_)[adv_index][imp_index] = fmin(1,
                                                       remaining_budget / ratios_[current_index].coefficient);
            //cout << "setting " << imp_index << ", " << adv_index << " to " << (*primal_sol)[adv_index][imp_index] << "\n";
            remaining_budget -= fmin(ratios_[current_index].coefficient,
                                     remaining_budget);
            current_index++;
        }
    }
    
    void NaiveMW::UpdateWeightedBudget() {
        weighted_budget_ = 0;
        for (int a = 0; a < num_advertisers_; ++a) {
            weighted_budget_ += advertiser_weights_[a] * (*budgets_)[a];
        }
        for (int i = 0; i < num_impressions_; ++i) {
            weighted_budget_ += impression_weights_[i];
        }
        cout << "set weighted budget to " << weighted_budget_ << "\n";
    }
    
    void NaiveMW::UpdateAvgSlacks(int t) {
        for (int i = 0; i < num_advertisers_; ++i) {
            advertiser_avg_slacks_[i] = (long double) (t - 1) / t * advertiser_avg_slacks_[i] +
                                        (long double) 1 / t * advertiser_slacks_[i];
        }
        for (int i = 0; i < num_impressions_; ++i) {
            impression_avg_slacks_[i] = (long double) (t - 1) / t * impression_avg_slacks_[i] +
                                        (long double) 1 / t * impression_slacks_[i];
        }
    }
    
    void NaiveMW::ReportWorstInfeasibility(int t) {
        long double max_infeasibility = 0.0;
        int max_infeasibility_index = -1;
        int max_infeasibility_type = -1;
        for (int i = 0; i < num_advertisers_; ++i) {
            if ((advertiser_avg_slacks_[i] > 0) and
                ((advertiser_avg_slacks_[i] / (*budgets_)[i]) > max_infeasibility)){
                max_infeasibility = advertiser_avg_slacks_[i] / (*budgets_)[i];
                max_infeasibility_index = i;
                max_infeasibility_type = 0;
            }
        }
        for (int i = 0; i < num_impressions_; ++i) {
            if ((impression_avg_slacks_[i] > 0) and
                (impression_avg_slacks_[i] > max_infeasibility)){
                max_infeasibility = impression_avg_slacks_[i];
                max_infeasibility_index = i;
                max_infeasibility_type = 1;
            }
        }
        cout << "At iteration ";
        cout << t;
        cout << ", max infeasiblity was ";
        cout << max_infeasibility;
        cout << " on constraint ";
        cout << max_infeasibility_index;
        cout << " of type ";
        if (max_infeasibility_type == 0) {
            cout << "advertiser";
        }
        if (max_infeasibility_type == 1) {
            cout << "impression";
        }
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
    
    void NaiveMW::UpdateAdvertiserWeights() {
        for (int j = 0; j < num_advertisers_; ++j) {
            long double tmp = advertiser_slacks_[j] / width_;
            if (abs(tmp) > 1) {
                cout << "Slack normalization error on advertiser weight\n";
            }
            if (tmp >= 0) {
                advertiser_weights_[j] = advertiser_weights_[j] * pow((1 + epsilon_), tmp);
            }
            else {
                advertiser_weights_[j] = advertiser_weights_[j] * pow((1 - epsilon_), -tmp);
            }
        }
    }
    
    void NaiveMW::UpdateImpressionWeights() {
        for (int j = 0; j < num_impressions_; ++j) {
            long double tmp = impression_slacks_[j] / width_;
            if (abs(tmp) > 1) {
                cout << "Slack normalization error on impression weight\n";
            }
            if (tmp >= 0) {
                impression_weights_[j] = impression_weights_[j] * pow((1 + epsilon_), tmp);
            }
            else {
                impression_weights_[j] = impression_weights_[j] * pow((1 - epsilon_), -tmp);
            }
        }
    }
    
    void NaiveMW::UpdateAdvertiserSlacks() {
        for (int j = 0; j < num_advertisers_; ++j) {
            advertiser_slacks_[j] = (-1) * (*budgets_)[j];
            for (__gnu_cxx::hash_map<int, long double>::iterator iter = (*primal_sol_)[j].begin();
                 iter != (*primal_sol_)[j].end();
                 ++iter) {
                advertiser_slacks_[j] += iter->second * (*bids_matrix_)[j].find(iter->first)->second;
            }
        }
    }
    
    void NaiveMW::UpdateImpressionSlacks() {
        for (int i = 0; i < num_impressions_; ++i) {
            impression_slacks_[i] = -1.0;
        }
        
        for (int j = 0; j < num_advertisers_; ++j) {
            for (__gnu_cxx::hash_map<int, long double>::iterator iter = (*primal_sol_)[j].begin();
                 iter != (*primal_sol_)[j].end();
                 ++iter) {
                impression_slacks_[iter->first] += iter->second;
            }
        }
    }

    
    void NaiveMW::RunNaiveMW(int num_iterations) {
        CreateRatios();
        for (int t = 1; t < num_iterations; ++t) {
            if (t == 10000) {
                epsilon_ = epsilon_ / 5;
            }
            
            if (t == 20000) {
                epsilon_ = epsilon_ / 2;
            }
            /*
            if (t == 10000) {
                epsilon_ = epsilon_ / 2;
            }
            */
            Instance::ResetPrimal(primal_sol_);
            UpdateWeightedBudget();
            NaiveMWPrimal();
            UpdateAdvertiserSlacks();
            UpdateImpressionSlacks();
            UpdateAvgSlacks(t);
            Instance::UpdateAvgPrimal(t, primal_sol_, avg_primal_sol_);
            UpdateAdvertiserWeights();
            UpdateImpressionWeights();
            ReportWorstInfeasibility(t);
            CalculateRevenue();
        }
    }
    
    void NaiveMW::CalculateRevenue() {
        long double MW_revenue = 0;
        for (int a = 0; a < (*primal_sol_).size(); ++a) {
            for (__gnu_cxx::hash_map<int, long double>::iterator iter = (*primal_sol_)[a].begin();
                 iter != (*primal_sol_)[a].end();
                 ++iter) {
                MW_revenue += iter->second * (*bids_matrix_)[a][iter->first];
            }
        }
        cout << "MW rev ";
        cout << MW_revenue;
        cout << "\n";
    }
}