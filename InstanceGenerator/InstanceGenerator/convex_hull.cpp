//
//  convex_hull.cpp
//  InstanceGenerator
//
//  Created by Dragos Ciocan on 5/21/13.
//  Copyright (c) 2013 Dragos Ciocan. All rights reserved.
//

#include "convex_hull.h"


// 2D cross product of OA and OB vectors, i.e. z-component of their 3D cross product.
// Returns a positive value, if OAB makes a counter-clockwise turn,
// negative for clockwise turn, and zero if the points are collinear.
coord2_t cross(const Point &O, const Point &A, const Point &B)
{
    return (A.x - O.x) * (B.y - O.y) - (A.y - O.y) * (B.x - O.x);
}

// Returns a list of points on the convex hull in counter-clockwise order.
// Note: the last point in the returned list is the same as the first one.
vector<Point> convex_hull(vector<Point> P)
{
    int n = P.size(), k = 0;
    vector<Point> H(2*n);
    
    // Sort points lexicographically
    sort(P.begin(), P.end());
    
    //cout << "printing points within hull function \n";
    //for (int p = 0; p < P.size(); ++p) {
    //    cout << P[p].x << ", " << P[p].y << "\n";
    //}
    
    
    // Build lower hull
    for (int i = 0; i < n; i++) {
        //cout << "k = " << k << "\n";
        //cout << "product is " << cross(H[k-2], H[k-1], P[i]) << "\n";
        while (k >= 2 && cross(H[k-2], H[k-1], P[i]) <= 0.000000000001) k--;
        // H[k++] = P[i];
        H[k] = P[i];
        k++;
        //cout << "set at index " << k - 1 << " equal to " << P[i].x << ", " << P[i].y << "\n";
    }
    /*
    int poop = k;
    cout << "printing lower hull \n";
    for (int p = 0; p <= poop; ++p) {
        cout << H[p].x << ", " << H[p].y << "\n";
    }
    */
    // Build upper hull
    /*for (int i = n-2, t = k+1; i >= 0; i--) {
        while (k >= t && cross(H[k-2], H[k-1], P[i]) <= 0) k--;
        H[k++] = P[i];
    }
    */
    
    H.resize(k);
    // original line H.resize(k-1);
    /*
    cout << "printing hull \n";
    for (int p = 0; p < H.size(); ++p) {
        cout << H[p].x << ", " << H[p].y << "\n";
    }
    */ 
    /*
    cout << "printing higher hull \n";
    for (int p = poop + 1; p < H.size(); ++p) {
        cout << H[p].x << ", " << H[p].y << "\n";
    }
    */
    return H;
}