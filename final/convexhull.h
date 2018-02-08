#ifndef CONVEXHULL_H
#define CONVEXHULL_H

#include <bits/stdc++.h>
#include "fundamentals.h"

using namespace std;

vector<point> graham_scan(const vector<point>&);
vector<point> jarvis_march(const vector<point>&);
pair<point, point> upper_bridge(const vector<point>&, long double);
vector<point> upper_hull(const vector<point>&, point, point);
vector<point> kirkpatrick_seidel(const vector<point>&);

#include "convexhull.cpp"

#endif
