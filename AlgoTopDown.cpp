#include "UtilityFunctions.h"
#include "Algo2Functions.h"

using namespace DGtal;
using namespace Z2i;

string displayVector2(std::vector<int> v) {
  return "(" + std::to_string(v[0]) + ", "+ std::to_string(v[1]) + ", "+ std::to_string(v[2]) + ")" ;
}

int main(int , char**) {
  //int r=5;//square of radius
  int mu=5;
  //Section 9: Discrete rotation tree construction: Top-down
  
  //Generate sorted triplet B_mu sort without removing double angles
  //vector<vector<vector<int> > > vecTriplet = sortTripletNew(mu);

  //Generate sorted triplet B_mu sort with removing double angles
  vector<vector<int> > vecSortedTriplet = sortTriplet(mu);
  //Get triplets of C_rho from B_mu (Fonction Build: Line 1)
  vector<vector<vector<int> > > vecTriplet (mu+1, vector<vector<int> >{});
  for(size_t it=0; it<vecSortedTriplet.size(); it++) {
    int p = vecSortedTriplet.at(it).at(0);
    int q = vecSortedTriplet.at(it).at(1);
    int k = vecSortedTriplet.at(it).at(2);
    int rho = p*p+q*q;
    vecTriplet.at(rho).push_back(vecSortedTriplet.at(it));
  }

  //A iteration S_{rho-1} -> S_{rho}
  vector<nodeTree*> tree;
  //Level 0: 0 -> 2*pi
  nodeTree* root = createNodeTree(0,vector<int> {0,0,0}, vector<int> {0,0,0});
  updateOmega(root, 0);//Omega of the root is 0
  updateBijectivity(root, true);
  updateInjectivity(root, true);//Bijective -> Injective
  tree.push_back(root);
  //Level 1: rho = 1
  vector<vector<int> > tripletL1 = vecTriplet.at(1);//Get triplet of level 1    
  //Add virtual node of 0 in tripletL1
  tripletL1.insert(tripletL1.begin(), {0,0,1});
  //Left node [last triplet -> first triplet]
  int lastIdNode = 1;
  //Middle nodes
  //First node that in loop
  nodeTree* treeNode = createNodeTree(lastIdNode,tripletL1.back(), tripletL1.front());
  updateAlpha(treeNode);
  addParent(treeNode, root->idNode, root);
  updateBijectivity(treeNode, true);
  tree.push_back(treeNode);
  lastIdNode += 1;
  for(size_t it=0; it<tripletL1.size()-1; it++) {
    nodeTree* treeNode = createNodeTree(lastIdNode+it,tripletL1.at(it), tripletL1.at(it+1));
    updateAlpha(treeNode);
    addParent(treeNode, root->idNode, root);
    updateBijectivity(treeNode, true);
    updateInjectivity(treeNode, true);//Bijective -> Injective
    tree.push_back(treeNode);
  }
  lastIdNode = tree.back()->idNode+1;
  
  //Associate Bij angles to each segment of the fist level
  findPythagoreAngleForSegments(tree, mu);
  
  //Create a graph/list of level 1
  node* currentLevel = NULL;
  for(size_t it=1; it<tree.size(); it++) {
    insertEnd(&currentLevel, tree.at(it)->idNode, tree.at(it));
  }

  //Loop for level 2 to mu (Algo 2: Lines 3-4)
  for(int rho=2; rho<=mu; rho++) {
    //Get triplets at level rho
    vector<vector<int> > triplet_rho = vecTriplet.at(rho);
    cout<<"------ Triplets at L="<<rho<<": "<<endl;
    for(size_t it=0; it<triplet_rho.size(); it++) {
      printf("[%d, %d, %d] (%3.2f) \n", triplet_rho.at(it).at(0), triplet_rho.at(it).at(1), triplet_rho.at(it).at(2), computeAngle(triplet_rho.at(it))*180/M_PI);
    }
    cout<<"------"<< endl;
    
    //Version with two lists in //
    if(triplet_rho.size()!=0) {
      //Insert the triplets to the tree
      lastIdNode = insertNewTripletsToTreeBijective(tree, &currentLevel, lastIdNode, triplet_rho);
      //Clean tree: remove intermediate nodes having alpha > omega
      for(size_t it=tree.size()-1; it>0; it--) {//First node is the root
        nodeTree* aNode = tree.at(it);
        assert(aNode->idNode != aNode->idParent);
        if(aNode->omega>0 && aNode->alpha > aNode->omega)//intermediate node
          tree.erase (tree.begin()+it);
      }
      
      //Display tree
      cout<<"Display clean tree at level "<<rho<<endl;
      for(auto e : tree)
        cout<<"id="<<e->idNode<<": "<<displayVector2(e->tripletLeft)<<" --> "<<displayVector2(e->tripletRight)<<" >>> idParent="<<e->idParent<<" with (alpha="<<e->alpha<<", omega="<<e->omega<<")"<<endl;
    }
    
    // Check for injectivity
    for(size_t it=0; it<triplet_rho.size(); it++) {
      vector<pair<vector<int>, bool > > vecInj = genVecInjectiveAngle(triplet_rho.at(it));
      
      //Map to currentLevel
      //Fist node (loop node) is bijective, therefore injective
      updateInjectivity(currentLevel->data.p_segment, true);
      //Handling the next nodes
      node* head = currentLevel->next;
      while (head != NULL) {
        nodeTree* current_seg = head->data.p_segment;
        nodeTree* parent_seg = current_seg->p_parent;
        vector<int> left_triplet = current_seg->tripletLeft;
        vector<int> right_triplet = current_seg->tripletRight;
        bool isInj = false;
        if(current_seg->isBijective == true)
          isInj = true; //bijective -> injective
        else {
          if(parent_seg->isInjective == false){//parent non inj -> children non inj
            isInj = false;
          }
          else { //only if the parent is injective then verify the injectivity of the children
            isInj = injectiveVerif(vecInj, right_triplet);
          }
        }
        if(current_seg->alpha == rho) //update if it is first appear in the tree
          updateInjectivity(current_seg, isInj);
        //Move to the next segment
        head = head->next;
      }
    }
  }
  
  //Update omega (Algo 2: Lines 5-6)
  for(size_t it=0; it<tree.size(); it++) {
    nodeTree* aNode = tree.at(it);
    if(aNode->omega==-1 /*&& aNode->alpha == mu*/)
      updateOmega(aNode, mu);
  }

  cout<<"------- Display clear tree with bijectivity from vector"<<endl;
  for(auto e : tree){
    //cout<<"id="<<e->idNode<<": "<<displayVector2(e->tripletLeft)<<" --> "<<displayVector2(e->tripletRight)<<" >>> idParent="<<e->idParent<<" with (alpha="<<e->alpha<<", omega="<<e->omega<<") and bijectiviy="<<e->isBijective<<endl;
    if(e->omega!=0 && e->alpha!=0) {
      cout<<displayVector2(e->tripletLeft)<<" --> "<<displayVector2(e->tripletRight)<<" >>> (alpha="<<e->alpha<<", omega="<<e->omega<<") "<<computeAngle(e->tripletLeft)*180/M_PI<<"; "<<computeAngle(e->tripletRight)*180/M_PI<<" et bij="<<e->isBijective<<", inj="<<e->isInjective<<endl;
    }
  }

  return 0;
}
