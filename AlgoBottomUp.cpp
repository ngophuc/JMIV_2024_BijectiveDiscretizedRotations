#include "UtilityFunctions.h"
#include "Algo1Functions.h"

using namespace DGtal;
using namespace Z2i;

string displayVector2(std::vector<int> v) {
  return "(" + std::to_string(v[0]) + ", "+ std::to_string(v[1]) + ", "+ std::to_string(v[2]) + ")" ;
}

int main(int , char**) {
  int r=5;//square of radius
  //Section 9: Discrete rotation tree construction: Bottom-up
  
  //Section 9.1: Initialization
  //Step 1+2: generate euclidean ball of radius rho (Eq 12) and generate triplet (Eq 14)
  vector<vector<int> > triplets = generateTripletHingeAngle(r);
  //Add virtual nodes of 0 in triplets
  triplets.insert(triplets.begin(), {0,0,1});
  //Step 3: sort triplets (Eq 16) and discard non prime triplets (Property 2)
  vector<vector<int> > sortedTriplets = sortTriplet(triplets);
  vector<vector<int> > vecTriplet = removeDoubleTriplet(sortedTriplets);
  //Step 4: hinge angles is encoded by (sorted) triplets already
  //Step 5: build the circle graph from sorted hinge angles
  node* circleGraph = createCircleGraph(vecTriplet);
  //printf("Created circular doubly linked list : \n");
  //displayList(circleGraph);
  //Step 5+6: build delta and sort hinge angles in function of delta
  vector<pair<int, node*> > sortedNormNode = sortGraphNode(circleGraph);
  /*
   //test remove nodes in the graph
   for(auto p : sortedNormNode) {
   printf("After delete %d : \n", p.second->data.idNode);
   deleteNode(&circleGraph, p.second);
   displayList(circleGraph);
   }
   */
  
  //section 9.2 Watershed tree construction (Algorithm 1)
  //lines 1 to 4
  vector<nodeTree*> tree = createFirstTree(circleGraph, r);
  //cout<<"Display tree from circleGraph"<<endl;
  //displayTree(circleGraph);
  int idNode;
  //lines 5 to 19
  //Each iteration, a hinge angle is removed from the circular graph
  for(auto aNode : sortedNormNode) {
    idNode = aNode.second->data.idNode;
    if(aNode!=sortedNormNode.back()) {//Not remove the last element (keep one node in the tree)
      //line 6: processing nodes by the sorted norm in sortedNormNode
      node* aNodeCircleGraph = aNode.second; //node h to be removed from the circular graph
      int normDeletedNode = aNode.first; //norm of node to be removed
      //line 7: set of segments being adjacent to h
      node* aNodePrev = aNodeCircleGraph->prev;
      nodeTree* aNodeTree = aNodeCircleGraph->data.p_segment;//node of tree
      nodeTree* aNodeTreePrev = aNodePrev->data.p_segment;//prev node of node in the tree
      //line 10-18: there exist two distinct segments
      //line 11: create a new segment by merging two segments being adjacent to h
      nodeTree* aNodeParent = createNodeTree (tree.size(), aNodeTreePrev->tripletLeft,aNodeTree->tripletRight);//new node
      //line 12: update omega of the new node in the lower level
      updateOmega(aNodeParent, normDeletedNode-1);
      //line 13: add new node to the tree
      tree.push_back(aNodeParent);
      //line 14: update tree by adding the new segment
      updatePSegment(aNodePrev, aNodeParent);
      //lines 15 to 18: connect the two existing segments to the parent
      //line 18: make connections aNodeTree, aNodeTreePrev -> aNodeParent
      addParent(aNodeTree, aNodeParent->idNode, aNodeParent);
      addParent(aNodeTreePrev, aNodeParent->idNode, aNodeParent);
      //line 16: delete aNodeCircleGraph in the circular graph
      deleteNode(&circleGraph, aNodeCircleGraph);
      //line 17: update alpha of two segments after the fusion
      updateAlpha(aNodeTree, normDeletedNode);
      updateAlpha(aNodeTreePrev, normDeletedNode);
      /*
       //display circular graph
       printf("----- Display circular graph, after delete %d : \n", idNode);
       displayList(circleGraph);
       //display tree
       cout<<"Display tree from vector"<<endl;
       for(auto e : tree)
       cout<<"id="<<e->idNode<<": "<<displayVector2(e->tripletLeft)<<" --> "<<displayVector2(e->tripletRight)<<" >>> idParent ="<<e->idParent<<" avec (alpha="<<e->alpha<<", omega="<<e->omega<<")"<<endl;
       */
    }
    //lines 8-9: one segment
    //update alpha for last node (rood of the tree)
    updateAlpha(tree.back(), 0);
    tree.back()->idParent=tree.back()->idNode;//root has itself as parent
    /* CHECK
     //display last level
     //display circular graph
     printf("----- Display circular graph(last) \n");
     displayList(circleGraph);
     //display tree
     cout<<"Display tree from vector (last)"<<endl;
     for(auto e : tree)
     cout<<"id="<<e->idNode<<": "<<displayVector2(e->tripletLeft)<<" --> "<<displayVector2(e->tripletRight)<<" >>> idParent="<<e->idParent<<" avec (alpha="<<e->alpha<<", omega="<<e->omega<<")"<<endl;
     */
  }
  
  //lines 20 to 26: unbinary the tree
  //Line 20: go through the tree from the root to the leaves
  //Att: the root (alpha=0 and omega=0: it=tree.size()-1) is a particular note that we keep anyway
  for(int it=tree.size()-2; it>=0; it--) {
    //Line 21: for each node
    nodeTree* aNode = tree.at(it);
    nodeTree* aNodeParent = aNode->p_parent;//Except the Root has no parent, all other node has a parent
    //Line 21: check wherether the node's parent is removed or kept
    if(aNodeParent->alpha > aNodeParent->omega) {//the parent node needs to be removed
      //Lines 25 to 27: the node takes the grandparent as its parent
      nodeTree* aNodeGrandParent = aNodeParent->p_parent;
      addParent(aNode, aNodeGrandParent->idNode, aNodeGrandParent);
    }
  }
  
  /* CHECK */
  cout<<"------- Display full tree from vector"<<endl;
  for(auto e : tree){
    if(e->alpha > e->omega) {
      cout<<"e->id ="<<e->idNode<<" ---> remove..."<<endl;
    }
   cout<<"id="<<e->idNode<<": "<<displayVector2(e->tripletLeft)<<" --> "<<displayVector2(e->tripletRight)<<" >>> idParent="<<e->idParent<<"with (alpha="<<e->alpha<<", omega="<<e->omega<<")"<<endl;
  }
   
  //Clean tree: retrive only node having alpha < omega
  vector<nodeTree*> clean_tree;
  for(auto e : tree) {
    if(e->alpha <= e->omega) {//Keep only nodes having alpha <= omega
      clean_tree.push_back(e);
    }
  }
  
  cout<<"------- Display clean tree from vector"<<endl;
  for(auto e : clean_tree){
    cout<<"id="<<e->idNode<<": "<<displayVector2(e->tripletLeft)<<" --> "<<displayVector2(e->tripletRight)<<" >>> idParent="<<e->idParent<<" with (alpha="<<e->alpha<<", omega="<<e->omega<<")"<<endl;
  }
  
  //Section 10.5.2: bijective discrete rotations
  vector<vector<int> > pAngles = generateTripletPythagoreAngles(r);
  for(auto angle : pAngles) {
    nodeTree* node = findNodeOfPythagoreAngle(angle, clean_tree, r);
    if(node!=NULL) {
      updateBijectivity(node, true);
      nodeTree* nodeParent = node->p_parent;
      while (nodeParent->alpha!=0 && nodeParent->isBijective!=true) {
        updateBijectivity(nodeParent, true);
        nodeParent = nodeParent->p_parent;
      }
      if(nodeParent->alpha==0)//bijectivity of the root
        updateBijectivity(nodeParent, true);
    }
  }
  
  cout<<"------- Display clear tree with bijectivity from vector"<<endl;
  for(auto e : clean_tree){
    cout<<"id="<<e->idNode<<": "<<displayVector2(e->tripletLeft)<<" --> "<<displayVector2(e->tripletRight)<<" >>> idParent="<<e->idParent<<" with (alpha="<<e->alpha<<", omega="<<e->omega<<") and bijectiviy="<<e->isBijective<<endl;
  }
  
  return 0;
}
