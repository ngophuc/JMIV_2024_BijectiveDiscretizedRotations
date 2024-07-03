#ifndef UtilityFunction_H
#define UtilityFunction_H

#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <assert.h>

#include "DGtal/helpers/StdDefs.h"
#include "DGtal/base/Common.h"

#include "Tree.h"
//#include "CircularDoublyLinkedList.h"
//#include "DoublyLinkedList.h"

using namespace DGtal;
using namespace Z2i;
using namespace std;

//Compute the norme of a 2D vector
int normVector(const vector<int>& v);

//Generate Pythagore bijective angles in Z2
vector<vector<int> > generateTripletPythagoreAngles(int r=20);

//Exact comparison between two triplets encoded hinges angles
int compareCosinusHingeAngles(const vector<int>& a1, const vector<int>& a2);
int compareHingeAngles(const vector<int>& a1, const vector<int>& a2);

//Exact comparison of a pythagore and an hinge angles
int compareCosinusHingePythagoreAngles(const vector<int>& triplet, const vector<int>& pythagore);
int compareTripletPythagoreAngles(const vector<int>& triplet, const vector<int>& pythagore);

//Other functions for the tests
double computeAngle (vector<int> triplet);
double calculCosinus(long int p, long int q, long int k);
double angle_from_sin_cos(double sinx, double cosx);
#endif // UtilityFunction_H
