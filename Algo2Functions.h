#ifndef Algo2Function_H
#define Algo2Function_H

#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <assert.h>

#include "DGtal/helpers/StdDefs.h"
#include "DGtal/base/Common.h"

#include "Tree.h"
#include "DoublyLinkedList.h"
#include "UtilityFunctions.h"

using namespace DGtal;
using namespace Z2i;
using namespace std;

vector<Point> genBall(int r=10);
vector<vector<Point> > generatePoints(int mu);
vector<vector<vector<int> > > sortTripletNew(int r=10);
vector<vector<int> > sortTriplet(int r=10);

vector<vector<int> > generateTripletHingeAngleFull(int mu=10);
vector<vector<int> > generateTripletHingeAngle(int mu=10);

vector<nodeTree*> findPythagoreAngleForSegments(vector<nodeTree*>& tree, int mu);
vector<vector<int> > findPythagoreAngleForSegments(const vector<int>& tLeft, const vector<int>& tRight, const vector<vector<int> >& bAngles); 

bool bijectiveBound(nodeTree** newNode, const vector<vector<int> >& pAngles);
int insertNewTripletsToTreeBijective(vector<nodeTree*>& tree, node** head, int idNewNode, const vector<vector<int> >& triplets);

bool injectiveVerif(const vector<pair<vector<int>, bool > >& vecInj, const vector<int>& triplet);
vector<pair<vector<int>, bool > > genVecInjectiveAngle(const vector<int>& triplet);
#endif // Algo2Function_H
