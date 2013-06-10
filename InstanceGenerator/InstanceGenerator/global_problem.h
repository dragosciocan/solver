
//
//  global_problem.h
//  InstanceGenerator
//
//  Created by Dragos Ciocan on 2/14/13.
//  Copyright (c) 2013 Dragos Ciocan. All rights reserved.
//

#ifndef InstanceGenerator_global_problem_h
#define InstanceGenerator_global_problem_h

#include <cmath>
#include <iostream>
#include <vector>

#include "subproblem.h"

using namespace std;

namespace distributed_solver {
    class GlobalProblem {
    public:
        int num_partitions_;
        long double budget_;
        vector<pair<int, long double> > budget_allocation_;
        vector<Subproblem> subproblems_;
        vector<long double> slacks_;
        
        GlobalProblem(int num_partitions, long double max_bid, long double advertiser_indegree);
        void InitializeInstance();
        void InitializeBudgetAllocation();
        void ConstructPrimal(vector<__gnu_cxx::hash_map<int, long double> >* primal_sol, int iteration);
        
    private:
        void FindOptimalBudgetAllocation();
        void ConstructSubproblemPrimal(vector<__gnu_cxx::hash_map<int, long double> >* primal_sol,
                                       int subproblem_index, long double budget_allocation, int opt_region);
        void ResetPrimal(vector<__gnu_cxx::hash_map<int, long double> >* primal_sol);
        
        long double numerical_accuracy_tolerance_;
        void SetTolerance(long double numerical_accuracy_tolerance);
        
        long double primal_assignment_test_;
    };
    
    class Slope {
    public:
        long double slope_;
        int subproblem_index_;
        int region_index_;
        Slope(long double slope, int subproblem_index, int region_index);
    };
    
    struct compare_Slope
    {
        bool operator() (const Slope & lhs, const Slope & rhs) {
            return lhs.slope_ > rhs.slope_;
        }
    };

}

#endif
