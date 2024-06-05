#ifndef Tree_h
#define Tree_h

#include <vector>
#include <cassert>

using namespace std;

typedef struct nodeTree{
  int idNode;
  int idParent;
  bool isBijective;
  int level;
  nodeTree* p_parent;
  vector<int> tripletLeft;
  vector<int> tripletRight;
  int alpha;
  int omega;
} nodeTree;

nodeTree* createNodeTree (int id, vector<int> tf, vector<int> tr);

nodeTree* addParent(nodeTree* node, int idF, nodeTree* pF);

nodeTree* updateOmega(nodeTree* node, int o);
nodeTree* updateAlpha(nodeTree* node, int a);
nodeTree* updateBijectivity(nodeTree* node, bool b);

void updateLevelBottomUp(vector<nodeTree*> tree);

#endif //Tree_h
