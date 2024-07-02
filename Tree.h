#ifndef Tree_h
#define Tree_h

#include <vector>
#include <cassert>

using namespace std;

typedef struct nodeTree{
  int idNode;
  int idParent;
  bool isBijective=false;
  nodeTree* p_parent;
  vector<int> tripletLeft;
  vector<int> tripletRight;
  int alpha=-1;
  int omega=-1;
  vector<vector<int> > pAngles;//List of bijective angles belonging to the segment
} nodeTree;

nodeTree* createNodeTree (int id, vector<int> tf, vector<int> tr);

nodeTree* addParent(nodeTree* node, int idF, nodeTree* pF);

nodeTree* updateOmega(nodeTree* node, int o);
nodeTree* updateAlpha(nodeTree* node, int a);
nodeTree* updateBijectivity(nodeTree* node, bool b);

#endif //Tree_h
