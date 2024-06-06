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
