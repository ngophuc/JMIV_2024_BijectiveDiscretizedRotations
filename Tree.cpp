#include "Tree.h"

nodeTree* createNodeTree (int id, vector<int> tf, vector<int> tr) {
  nodeTree* aNode = new nodeTree();
  // assign data to the node
  aNode->idNode = id;
  aNode->tripletLeft = tf;
  aNode->tripletRight = tr;
  aNode->isBijective = false;
  aNode->alpha = 0;
  aNode->omega = -1;
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

//By Eq. 48
nodeTree* updateAlpha(nodeTree* node){
  vector<int> tL = node->tripletLeft;
  vector<int> tR = node->tripletRight;
  
  int n1 = tL.at(0)*tL.at(0)+tL.at(1)*tL.at(1);
  int n2 = tR.at(0)*tR.at(0)+tR.at(1)*tR.at(1);
  node->alpha = n1>n2 ? n1 : n2;
  return node;
}

nodeTree* updateBijectivity(nodeTree* node, bool b){
  node->isBijective = b;
  return node;
}

nodeTree* updateInjectivity(nodeTree* node, bool b){
  node->isInjective = b;
  return node;
}
