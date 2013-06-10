//
//  global_problem.cpp
//  InstanceGenerator
//
//  Created by Dragos Ciocan on 2/14/13.
//  Copyright (c) 2013 Dragos Ciocan. All rights reserved.
//

#include <algorithm>
#include <cmath>
#include "global_problem.h"

namespace distributed_solver {
    GlobalProblem::GlobalProblem(int num_partitions, long double max_bid, long double advertiser_indegree) {
        num_partitions_ = num_partitions;
        budget_ = 0;
    }
    
    void GlobalProblem::FindOptimalBudgetAllocation() {
        long double remaining_budget = budget_;
        
        // Create list of slopes.
        vector<Slope> slopes;
        for (int i = 0; i < num_partitions_; ++i) {
            // Reset budget allocation.
            budget_allocation_[i].second = 0;
            
            if (subproblems_[i].envelope_points_.size() > 1) {
                // -1 because at last envelope return on budget is 0.
                for (int j = 0; j < subproblems_[i].envelope_points_.size() - 1; ++j) {
                    slopes.push_back(Slope(subproblems_[i].envelope_points_[j].first, i, j));
                }
            }
        }
        // Sort list of slopes.
        sort(slopes.begin(), slopes.end(), compare_Slope());
        
        // Allocate budgets in decreasing order of slopes.
        int sp_index;
        for (int k = 0; k < slopes.size(); ++k) {
            sp_index = slopes[k].subproblem_index_;
            long double allocation_increase = min(remaining_budget,
                                                       subproblems_[sp_index].budget_cutoffs_[slopes[k].region_index_ + 1] -
                                                       subproblems_[sp_index].budget_cutoffs_[slopes[k].region_index_]);
            if ((allocation_increase > 0) && (remaining_budget > 0)) {
                budget_allocation_[sp_index] = make_pair(slopes[k].region_index_, budget_allocation_[sp_index].second + allocation_increase);
            }
            remaining_budget = remaining_budget - allocation_increase;
        }
    }
    
    void GlobalProblem::ConstructPrimal(vector<__gnu_cxx::hash_map<int, long double> >* primal_sol, int iteration) {
        // Reset primal solution.
        GlobalProblem::ResetPrimal(primal_sol);
        
        // Find optimal budget allocations to problems.
        cout << "Solving subproblems \n";
        for (int i = 0; i < num_partitions_; ++i) {
            //subproblems_[i].SolveSubproblem(iteration, i);
            subproblems_[i].SolveSubproblemConvexHull(iteration, i);
        }
        cout << "Find optimal budget allocation \n";
        FindOptimalBudgetAllocation();
        
        // Calculate primal solution for each subproblem.
        cout << "Constructing primal \n";
        long double dual_val = 0;
        primal_assignment_test_ = 0;
        for (int i = 0; i < num_partitions_; ++i) {
            if (budget_allocation_[i].second > 0) {
                long double u = subproblems_[i].envelope_points_[budget_allocation_[i].first].first;
                long double v = subproblems_[i].envelope_points_[budget_allocation_[i].first].second;
                dual_val += u * budget_allocation_[i].second + v;
                ConstructSubproblemPrimal(primal_sol, i, budget_allocation_[i].second, budget_allocation_[i].first);
            }
        }
        cout << "Dual Value = ";
        cout << dual_val;
        cout << "\n";
    }
    
    void GlobalProblem::ConstructSubproblemPrimal(vector<__gnu_cxx::hash_map<int, long double> >* primal_sol,
                                                  int subproblem_index, long double budget_allocation, int opt_region) {
        // Figure out opt u, v.
        long double u = subproblems_[subproblem_index].envelope_points_[opt_region].first;
        long double v = subproblems_[subproblem_index].envelope_points_[opt_region].second;
        
        long double allocation_value = 0;
         
        // If optimum is u = 0, optimal allocation is greedy wrt price.
        if (u == 0) {
            int max_price_index;
            long double max_price = 0;
            for (int i = 0; i < subproblems_[subproblem_index].num_vars_; ++i) {
                if ((max_price < subproblems_[subproblem_index].constraints_[i].price_) and (subproblems_[subproblem_index].constraints_[i].is_active_)) {
                    max_price_index = i;
                    max_price = subproblems_[subproblem_index].constraints_[i].price_;
                }
            }
            
            (*primal_sol)[subproblems_[subproblem_index].advertiser_index_->at(max_price_index)][subproblem_index] = 1;
            
            allocation_value = (*primal_sol)[subproblems_[subproblem_index].advertiser_index_->at(max_price_index)][subproblem_index] * subproblems_[subproblem_index].constraints_[max_price_index].price_;
            
            primal_assignment_test_ += allocation_value;
        }
        
        // If optimum is v = 0, optimal allocation is greedy wrt price/weight ratio.
        if (v == 0) {
            int max_ratio_index;
            long double max_ratio = 0;
            for (int i = 0; i < subproblems_[subproblem_index].num_vars_; ++i) {
                if (max_ratio < (subproblems_[subproblem_index].constraints_[i].price_ /
                                 subproblems_[subproblem_index].constraints_[i].coefficient_) and
                    (subproblems_[subproblem_index].constraints_[i].is_active_)) {
                    max_ratio_index = i;
                    max_ratio = (subproblems_[subproblem_index].constraints_[i].price_ /
                                 subproblems_[subproblem_index].constraints_[i].coefficient_);
                }
            }
             
            (*primal_sol)[subproblems_[subproblem_index].advertiser_index_->at(max_ratio_index)][subproblem_index] =budget_allocation / subproblems_[subproblem_index].constraints_[max_ratio_index].coefficient_;
           
            allocation_value = (*primal_sol)[subproblems_[subproblem_index].advertiser_index_->at(max_ratio_index)][subproblem_index] * subproblems_[subproblem_index].constraints_[max_ratio_index].price_;
            primal_assignment_test_ += allocation_value;
        }
        
        // If opt is u, v > 0, the optimum whp only has two positive allocations, which are the
        // solutions to a system of 2 equations.
        if ((u > 0) and (v > 0)) {
            vector<int> tight_constraint_indices;
            for (int i = 0; i < subproblems_[subproblem_index].num_vars_; ++i) {
                if (subproblems_[subproblem_index].constraints_[i].is_active_) {
                    long double slack = subproblems_[subproblem_index].constraints_[i].price_ -
                    (u * subproblems_[subproblem_index].constraints_[i].coefficient_ + v);
                    if (slack < 0) { slack = (-1) * slack;}
                    if (slack < 0.000000000000000001) {
                        tight_constraint_indices.push_back(i);
                    }
                }
            }
            if (tight_constraint_indices.size() > 2) {
                cout << "ERROR, PERTURB PRICES \n";
            }
            if (tight_constraint_indices.size() == 1) {
                (*primal_sol)[subproblems_[subproblem_index].advertiser_index_->
                              at(tight_constraint_indices[0])]
                            [subproblem_index] = fmin(budget_allocation / subproblems_[subproblem_index].constraints_[tight_constraint_indices[0]].coefficient_, 1);
                allocation_value = (*primal_sol)[subproblems_[subproblem_index].advertiser_index_-> at(tight_constraint_indices[0])][subproblem_index] * subproblems_[subproblem_index].constraints_[tight_constraint_indices[0]].coefficient_;
                primal_assignment_test_ += allocation_value;
            }
            if (tight_constraint_indices.size() == 2) {
                int first_index = tight_constraint_indices[0];
                int second_index = tight_constraint_indices[1];
                
                long double x_1 = (budget_allocation -
                              subproblems_[subproblem_index].constraints_[second_index].coefficient_) / (subproblems_[subproblem_index].constraints_[first_index].coefficient_ - subproblems_[subproblem_index].constraints_[second_index].coefficient_);
                (*primal_sol)[subproblems_[subproblem_index].advertiser_index_->at(first_index)][subproblem_index] = x_1;
                (*primal_sol)[subproblems_[subproblem_index].advertiser_index_->at(second_index)][subproblem_index] = 1 - x_1;
                
                allocation_value = x_1 * subproblems_[subproblem_index].constraints_[first_index].price_ + (1-x_1)*subproblems_[subproblem_index].constraints_[second_index].price_;
                primal_assignment_test_ += allocation_value;
                
            }
        }
        if (((allocation_value - (u * budget_allocation + v)) > 0.0000001) ||
            ((allocation_value - (u * budget_allocation + v)) < -0.0000001)) {
            cout << "**************Error at problem************ " << subproblem_index << " equal to " <<
                         allocation_value << " - " << (u * budget_allocation + v) << "\n";
        }
    }
    
    void GlobalProblem::ResetPrimal(vector<__gnu_cxx::hash_map<int, long double> >* primal_sol) {
        for (int i = 0; i < (*primal_sol).size(); ++i) {
            for (__gnu_cxx::hash_map<int, long double>::iterator iter = (*primal_sol)[i].begin(); iter != (*primal_sol)[i].end(); ++iter) {
                iter->second = 0;
            }
        }
    }
    
    void GlobalProblem::InitializeBudgetAllocation() {
        for (int i = 0; i < num_partitions_; ++i) {
            budget_allocation_.push_back(make_pair(0, 0.0));
        }
    }
    
    void GlobalProblem::SetTolerance(long double numerical_accuracy_tolerance) {
        numerical_accuracy_tolerance_ = numerical_accuracy_tolerance;
    }
    
    Slope::Slope(long double slope, int subproblem_index, int region_index) {
        slope_ = slope;
        subproblem_index_ = subproblem_index;
        region_index_ = region_index;
    }
}