#include "Tree.h"

nodeTree* createNodeTree (int id, vector<int> tf, vector<int> tr) {
  nodeTree* aNode = new nodeTree();
  // assign data to the node
  aNode->idNode = id;
  aNode->tripletLeft = tf;
  aNode->tripletRight = tr;
  aNode->isBijective = false;
  return aNode;
}

nodeTree* addParent(nodeTree* node, int idF, nodeTree* pF) {
  node->idParent = idF;
  node->p_parent = pF;
  return node;
}

nodeTree* updateOmega(nodeTree* node, int o){
  node->omega = o;
  return node;
}

nodeTree* updateAlpha(nodeTree* node, int a){
  node->alpha = a;
  return node;
}

nodeTree* updateBijectivity(nodeTree* node, bool b){
  node->isBijective = b;
  return node;
}

void updateLevelBottomUp(vector<nodeTree*> tree){
  nodeTree* root = tree.back();
  root->level = 0;
  assert(root->idParent==root->idNode); //&& root->alpha=0
  for(int it=tree.size()-2; it>=0; it--) { //update level from root to leaves (back to front)
    nodeTree* aNode = tree.at(it);
    int idParent = aNode->idParent;
    aNode->level = tree.at(idParent)->level + 1;//increase in the level from its parent
  }
}