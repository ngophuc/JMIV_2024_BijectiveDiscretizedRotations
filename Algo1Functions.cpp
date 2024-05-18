#include "Algo1Functions.h"

//Step 1: generate euclidean ball of radius rho (Equation 12)
vector<Point> genBall(int r) { //m=mu**2, r=rho**2
  vector<Point> res;
  int x=0, y=0;
  while(y*y<=r) {
    if(y!=0) {
      res.push_back(Point(0,y));
      res.push_back(Point(0,-y));
    }
    else {
      res.push_back(Point(0,0));
    }
    x=1;
    while(x*x+y*y<=r) {
      if(y!=0) {
        res.push_back(Point(x,y));
        res.push_back(Point(-x,y));
      }
      res.push_back(Point(x,-y));
      res.push_back(Point(-x,-y));
      x++;
    }
    y++;
  }
  return res;
}

//Step 2: generate triplet (Equation 14)
vector<vector<int> > generateTriplet(int r) {
  vector<vector<int> > res;
  vector<Point> pts = genBall(r);
  for(size_t it=0; it<pts.size(); it++) {
    int p = pts.at(it)[0];//px
    int q = pts.at(it)[1];//py
    int pq_2 = p*p+q*q;
    int k=0;
    while (4*k*k+4*k+1 <= pq_2) {
      vector<int> v = {p, q, k};
      res.push_back(v);
      k++;
    }
    k=-1;
    while (4*k*k+4*k+1 <= pq_2) {
      vector<int> v = {p, q, k};
      res.push_back(v);
      k--;
    }
  }
  return res;
}

//Step 3 : sort triplets
vector<vector<int> > sortTriplet(const vector<vector<int> >& triplets) {
  //vector<vector<int> > pqk_trie = generateTriplet(r);
  vector<vector<int> > pqk_trie;
  for(size_t it=0; it<triplets.size(); it++)
    pqk_trie.push_back(triplets.at(it));
  
  std::sort(begin(pqk_trie), end(pqk_trie), [](vector<int> t1, vector<int> t2) {
    return compareHingeAngles(t1, t2)<0;
  });
  /* CHECK
  cout<<"After sorting"<<endl;
  for(auto p : pqk_trie) {
    cout<<"("<<p[0]<<","<<p[1]<<","<<p[2]<<") ==> "<<computeAngle(p)*180/M_PI<<endl;
  }
  */
  return pqk_trie;
}
                       
//Step 4: remove double triplets
vector<vector<int> > removeDoubleTriplet(const vector<vector<int> >& pqk_trie) {
  size_t it=0, it_next=1, it_min;
  double norme1, norme2, angle1, angle2;
  vector<vector<int> > pqk_non_double;
  while(it<pqk_trie.size() && it_next<pqk_trie.size()-1) {
    it_next = it+1;
    vector<int> t1 = pqk_trie.at(it);
    vector<int> t2 = pqk_trie.at(it_next);
    angle1 = computeAngle(t1);
    angle2 = computeAngle(t2);
    it_min = it;
    while(fabs(angle1-angle2)<1e-6 && it_next<pqk_trie.size()-1) {
      norme1 = t1[0]*t1[0]+t1[1]*t1[1];
      norme2 = t2[0]*t2[0]+t2[1]*t2[1];
      if(norme2<norme1)
        it_min = it_next;
      it_next++;
      t1 = pqk_trie.at(it_min);
      t2 = pqk_trie.at(it_next);
      angle1 = computeAngle(t1);
      angle2 = computeAngle(t2);
    }
    pqk_non_double.push_back(pqk_trie.at(it_min));
    it=it_next;
  }
  /* CHECK
  cout<<"Après supp"<<endl;
  for(auto p : pqk_non_double) {
    cout<<"("<<p[0]<<","<<p[1]<<","<<p[2]<<") ==> "<<computeAngle(p)*180/M_PI<<endl;
  }
  */
  return pqk_non_double;
}

//Step 5+6: build delta and sort hinge angles in function of delta
vector<pair<int, node*> > sortGraphNode(node* circleGraph) {
  vector<pair<int, node*> > sortedNormNode;
  //Step 5: build delta function
  node* tmpNode = circleGraph;
  node* lastNode = circleGraph->prev;
  while (tmpNode != lastNode) {
    vector<int> triplet = tmpNode->data.triplet;
    int n = normVector(triplet);
    sortedNormNode.push_back(make_pair(n, tmpNode));
    tmpNode = tmpNode->next;
  }
  //Add last elt
  vector<int> triplet = lastNode->data.triplet;
  int n = normVector(triplet);
  sortedNormNode.push_back(make_pair(n, lastNode));
  /* CHECK
   cout<<"Affiche vector pair"<<endl;
   for (auto e : sortedNormNode) {
   cout<<"norme="<<e.first<<"--> node:"<<displayNodeData(e.second)<<endl;
   }
  */
  
  //Step 6: sort hinge angles in function of delta (in decreasing order)
  std::sort(sortedNormNode.begin(), sortedNormNode.end(), [](const pair<int, node*>& t1, const pair<int, node*>& t2) {
    return t1.first>t2.first;
  });
  /* CHECK
   cout<<"Après le trie"<<endl;
   for (auto e : sortedNormNode) {
   cout<<"norme="<<e.first<<"--> node:"<<displayNodeData(e.second)<<endl;
   }
  */
  return sortedNormNode;
}

node* createCircleGraph(const vector<vector<int> >& pqk_trie) {
  node* circleGraph = NULL;
  for(size_t it=0; it<pqk_trie.size(); it++) {
    insertEnd(&circleGraph, it, pqk_trie.at(it));
  }
  return circleGraph;
}

//Algorithm 1: Watershed tree construction (lines 1 to 4)
vector<nodeTree*> createFirstTree(node* circleGraph, int mu) {
   vector<nodeTree*> tree;
   node* tmpNode = circleGraph;
   //all nodes of the tree execept the last node
   while (circleGraph!=NULL && tmpNode->next != circleGraph) {
     nodeTree* node = createNodeTree(tmpNode->data.idNode, tmpNode->data.triplet, tmpNode->next->data.triplet);
     updatePSegment(tmpNode,node); //update pointeur segment of the node
     //line 4: update omega at m (first level)
     updateOmega(node, mu);
     tree.push_back(node);
     tmpNode = tmpNode->next;
   }
   //last node (circular graph)
   if(tmpNode!=NULL) {
     nodeTree* node = createNodeTree(tmpNode->data.idNode, tmpNode->data.triplet, tmpNode->next->data.triplet);
     updatePSegment(tmpNode,node);
     //line 4: update omega at m (first level)
     updateOmega(node, mu);
     tree.push_back(node);
   }
   return tree;
}

/*___________________________________________________________*/
/* Check belonging Pythagore angle for intervals of the tree */
/*___________________________________________________________*/
 nodeTree* findNodeOfPythagoreAngle(const vector<int>& angle, const vector<nodeTree*>& tree) {
   int idLeft = 0;
   int idRight = tree.size();
   int idMid = 0;
   while (idLeft <= idRight) {
     idMid = (idLeft + idRight)/ 2;
     nodeTree* node = tree.at(idMid);
     vector<int> tL = node->tripletLeft;
     vector<int> tR = node->tripletRight;
     if(node->omega==5) {
       int cL = compareTripletPythagoreAngles(tL, angle);
       int cR = compareTripletPythagoreAngles(tR, angle);
       if(cL==-1 && cR==1) //tL < angle < tR
         return node;
       if(cL==1) //angle < tL => ignore right half
         idRight = idMid - 1;
       if(cR==-1) //angle > tR => ignore left half
         idLeft = idMid + 1;
     }
     else
       return NULL;
   }
   return NULL;
 }
