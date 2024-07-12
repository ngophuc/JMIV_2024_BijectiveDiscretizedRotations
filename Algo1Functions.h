#ifndef Algo1Function_H
#define Algo1Function_H

#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <assert.h>

#include "DGtal/helpers/StdDefs.h"
#include "DGtal/base/Common.h"

#include "Tree.h"
#include "CircularDoublyLinkedList.h"
#include "UtilityFunctions.h"

using namespace DGtal;
using namespace Z2i;
using namespace std;

//Step 1: generate euclidean ball of radius rho (Equation 12)
vector<Point> genBall(int r=10);
//Step 2: generate triplet (Equation 14)
vector<vector<int> > generateTripletHingeAngleFull(int r=10);
vector<vector<int> > generateTripletHingeAngle(int r=10);
//Step 3.1 : sort triplets (Equation 20)
vector<vector<int> > sortTriplet(const vector<vector<int> >& triplets);
//Step 3.2: remove double triplets (Property 2)
vector<vector<int> > removeDoubleTriplet(const vector<vector<int> >& triplets);
//Step 4: build the circle graph from sorted hinge angles 
node* createCircleGraph(const vector<vector<int> >& pqk_trie);
//Step 5+6: build delta and sort hinge angles in function of delta
vector<pair<int, node*> > sortGraphNode(node* circleGraph);

//Watershed tree construction (Algorithm 1: lines 1 to 4)
vector<nodeTree*> createFirstTree(node* circleGraph, int mu);

//Check belonging Pythagore angle for intervals of the tree
nodeTree* findNodeOfPythagoreAngle(const vector<int>& angle, const vector<nodeTree*>& tree, int mu);
#endif // Algo1Function_H
