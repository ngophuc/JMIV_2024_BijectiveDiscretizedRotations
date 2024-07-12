#include "Algo2Functions.h"

//Generate euclidean ball of radius rho (Equation 12)
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

//Generate points in B_mu, then k for each (p,q) and then sort them
vector<vector<Point> > generatePoints(int mu) {
  vector<Point> pts = genBall(mu);
  vector<vector<Point> > vecPoint (mu+1, vector<Point>{});
  for(size_t it=0; it<pts.size(); it++) {
    int p = pts.at(it)[0];
    int q = pts.at(it)[1];
    int rho = p*p+q*q;
    vecPoint.at(rho).push_back(pts.at(it));
  }
  return vecPoint;
}

//Sort the triplets
vector<vector<vector<int> > > sortTripletNew(int r) {
  //Gen points
  vector<vector<Point> > vecPoint = generatePoints(r);
  for(size_t it=0; it<vecPoint.size(); it++) {
    for(auto p : vecPoint.at(it))
      cout<<p<<", ";
    cout<<endl;
  }
  //Gen (p,q,k) for each point
  vector<vector<vector<int> > > vecPqk (r+1, vector<vector<int> >{});
  for(size_t it=0; it<vecPoint.size(); it++) {
    if(vecPoint.at(it).size()!=0) {
      int p = vecPoint.at(it).front()[0];
      int q = vecPoint.at(it).front()[1];
      int pq_2 = p*p+q*q;
      vector<vector<int> > res;
      for(size_t it_bis=0; it_bis<vecPoint.at(it).size(); it_bis++) {
        p = vecPoint.at(it).at(it_bis)[0];
        q = vecPoint.at(it).at(it_bis)[1];
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
      //sort res then push to the final vector
      std::sort(begin(res), end(res), [](vector<int> t1, vector<int> t2) {
        return compareHingeAngles(t1, t2)<0; });
      
      //Display sorted vector
      for(size_t it=0; it<res.size(); it++)
        printf("[%d, %d, %d] (%3.2f) ", res.at(it).at(0), res.at(it).at(1),res.at(it).at(2), computeAngle(res.at(it))*180/M_PI);
      cout<<endl;
      
      for(size_t it=0; it<res.size(); it++)
        vecPqk.at(pq_2).push_back(res.at(it));
    }
  }
  return vecPqk;
}

//For angles in [0, 2PI)
vector<vector<int> > generateTripletHingeAngleFull(int mu) {
  vector<vector<int> > res;
  //res.push_back({0, 0, 0})
  vector<Point> pts = genBall(mu);//generatePointDisk(r);
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

//For angles in (0,PI/4)
vector<vector<int> > generateTripletHingeAngle(int r) {
  vector<vector<int> > res;
  //res.push_back({0, 0, 0})
  vector<Point> pts = genBall(r);//generatePointDisk(r);
  for(size_t it=0; it<pts.size(); it++) {
    int p = pts.at(it)[0];//px
    int q = pts.at(it)[1];//py
    if(p>=0 && q>=0 && p<=q) { //0 <= p <= q
      int pq_2 = p*p+q*q;
      int k=0;
      while (4*k*k+4*k+1 <= pq_2) {
        vector<int> v = {p, q, k};
        if(computeAngle(v)<=M_PI_4)
          res.push_back(v);
        k++;
      }
      k=-1;
      while (4*k*k+4*k+1 <= pq_2) {
        vector<int> v = {p, q, k};
        if(computeAngle(v)<=M_PI_4)
          res.push_back(v);
        k--;
      }
    }
  }
  return res;
}

vector<vector<int> > sortTriplet(int r) {
  vector<vector<int> > pqk_trie = generateTripletHingeAngleFull(r);
  //trie
  std::sort(begin(pqk_trie), end(pqk_trie), [](vector<int> t1, vector<int> t2) {
    return compareHingeAngles(t1, t2)<0;
  });

  size_t it=0, it_next=1, it_min;
  double norme1, norme2, angle1, angle2;
  vector<vector<int> > pqk_non_double;
  while(it<pqk_trie.size() && it_next<pqk_trie.size()-1) {
    it_next = it+1;
    //cout<<"it="<<it<<" vs it_next="<<it_next<<endl;
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
  //Push the last angle
  pqk_non_double.push_back(pqk_trie.back());

  return pqk_non_double;
}

vector<nodeTree*> findPythagoreAngleForSegments(vector<nodeTree*>& tree, int mu) {
  vector<vector<int> > bAngles = generateTripletPythagoreAnglesFull(mu);
  std::sort(begin(bAngles), end(bAngles), [](vector<int> t1, vector<int> t2) {
    double a1 = angle_from_sin_cos(double(t1[1])/double(t1[2]), double(t1[0])/double(t1[2]));
    double a2 = angle_from_sin_cos(double(t2[1])/double(t2[2]), double(t2[0])/double(t2[2]));
    return a1<a2;
  });
  
  //Ignore the first segment corresponding to the root, start it_tree=1
  size_t it_tree=1, it_pythagore=0;
  //For fist angles for the loop segment
  vector<int> tL = tree.at(it_tree)->tripletLeft;
  vector<int> tR = tree.at(it_tree)->tripletRight;
  vector<int> pAngle = bAngles.at(it_pythagore);
  while (compareTripletPythagoreAngles(tR, pAngle)==1) { //tR>pAngle
    //cout<<" -> angle pythagore="<<angle_from_sin_cos( double(pAngle[1])/pAngle[2], double(pAngle[0])/pAngle[2] )*180/M_PI<<" -> hA left="<<computeAngle(tL)*180/M_PI<<" --> hA right="<<computeAngle(tR)*180/M_PI<<endl;
    tree.at(it_tree)->pAngles.push_back(pAngle);//add the pA to the node
    it_pythagore++;//Move to the next angle
    pAngle = bAngles.at(it_pythagore);
  }
  it_tree++;//Move to the next segment
  
  //Other segments
  while (it_tree<tree.size() && it_pythagore<bAngles.size()) {
    tL = tree.at(it_tree)->tripletLeft;
    tR = tree.at(it_tree)->tripletRight;
    pAngle = bAngles.at(it_pythagore);
    int cL = compareTripletPythagoreAngles(tL, pAngle);
    int cR = compareTripletPythagoreAngles(tR, pAngle);
    if(cL==1) {//tL > pA
      it_pythagore++;//Move to the next angle
    }
    else {//tL <= pA
      if(cR==1) { //tL < pA < tR
        tree.at(it_tree)->pAngles.push_back(pAngle);//add the pA to the node
        it_pythagore++;//Move to the next angle
      }
      else {// pA >= tR
        it_tree++;//Move to the next segment
      }
    }
  }
  
  //Last angles for the loop segment !
  while (it_pythagore<bAngles.size()) {
    tL = tree.at(1)->tripletLeft;
    tR = tree.at(1)->tripletRight;
    pAngle = bAngles.at(it_pythagore);
    tree.at(1)->pAngles.push_back(pAngle);//add the pA to the node
    it_pythagore++;//Move to the next angle
  }
  
  return tree;
}

vector<vector<int> > findPythagoreAngleForSegments(const vector<int>& tLeft, const vector<int>& tRight, const vector<vector<int> >& bAngles) {
  vector<vector<int> > res;
  if(bAngles.size()==0)
    return res;//Empty vector: non bij seg
  
  if(compareHingeAngles(tLeft, tRight)==1) { //tL>tR: case of loop node
    size_t it_pythagore = 0;
    vector<int> pAngle = bAngles.at(it_pythagore);
    while (it_pythagore<bAngles.size() && compareTripletPythagoreAngles(tRight, pAngle)==1) { //tR>pAngle
      res.push_back(pAngle);//add the pA to the node
      it_pythagore++;//Move to the next angle
      if(it_pythagore<bAngles.size())
         pAngle = bAngles.at(it_pythagore);
    }
    while (it_pythagore<bAngles.size() && compareTripletPythagoreAngles(tLeft, pAngle)==1) { //tL>pAngle
      it_pythagore++;//Move to the next angle
      if(it_pythagore<bAngles.size())
         pAngle = bAngles.at(it_pythagore);
    }
    while (it_pythagore<bAngles.size() && compareTripletPythagoreAngles(tLeft, pAngle)==-1) { //tL<pAngle
      res.push_back(pAngle);//add the pA to the node
      it_pythagore++;//Move to the next angle
      if(it_pythagore<bAngles.size())
        pAngle = bAngles.at(it_pythagore);
    }
  }
  else {//Other cases: tL<tR
    size_t it_pythagore = 0;
    vector<int> pAngle = bAngles.at(it_pythagore);
    while (it_pythagore<bAngles.size() && compareTripletPythagoreAngles(tLeft, pAngle)==1) { //tL>pAngle
      it_pythagore++;//Move to the next angle
      if(it_pythagore<bAngles.size())
        pAngle = bAngles.at(it_pythagore);
    }
    
    while (it_pythagore<bAngles.size() && compareTripletPythagoreAngles(tRight, pAngle)==1) { //tL<pAngle<tR
      res.push_back(pAngle);//add the pA to the node
      it_pythagore++;//Move to the next angle
      if(it_pythagore<bAngles.size())
        pAngle = bAngles.at(it_pythagore);
    }
  }
  return res;
}

bool bijectiveBound(nodeTree** newNode, const vector<vector<int> >& pAngles) {
  vector<int> tL = (*newNode)->tripletLeft;
  vector<int> tR = (*newNode)->tripletRight;
  vector<vector<int> > res = findPythagoreAngleForSegments(tL, tR, pAngles);

  //Add bij angles to the new node
  for(size_t it=0; it<res.size(); it++) {
    (*newNode)->pAngles.push_back(res.at(it));
    vector<int> pAngle = res.at(it);
  }

  return res.size()!=0;
}

int insertNewTripletsToTreeBijective(vector<nodeTree*>& tree, node** head, int idNewNode, const vector<vector<int> >& triplets) {
  bool isBij = false;
  int id_triplet = 0;
  node* temp = *head;
  nodeTree* current_seg=temp->data.p_segment; //=> loop node
  //Processing of the first segment
  vector<int> triplet = triplets.at(id_triplet);
  
  while (compareHingeAngles(triplet, current_seg->tripletRight)==-1) {//triplet < firstNode->tripletRight
    int normeTriplet = triplet[0]*triplet[0]+triplet[1]*triplet[1];
    //Insert triplet here
    vector<int> tripletLeft = current_seg->tripletLeft;
    vector<int> tripletRight = current_seg->tripletRight;
    vector<vector<int> > bAngles = current_seg->pAngles;
    //Insert node to the tree
    nodeTree* aNewNodeLeft = createNodeTree(idNewNode, tripletLeft, triplet);
    updateAlpha(aNewNodeLeft, normeTriplet);//Update alpha left node
    addParent(aNewNodeLeft, current_seg->idNode, current_seg);
    idNewNode++;
    nodeTree* aNewNodeRight = createNodeTree(idNewNode, triplet, tripletRight);
    updateAlpha(aNewNodeRight, normeTriplet);//Update alpha of right node
    addParent(aNewNodeRight, current_seg->idNode, current_seg);
    idNewNode++;
    //Check bijective from pAngles of the parent node
    isBij = bijectiveBound(&aNewNodeLeft, bAngles);
    updateBijectivity(aNewNodeLeft, isBij);
    isBij = bijectiveBound(&aNewNodeRight, bAngles);
    updateBijectivity(aNewNodeRight, isBij);
    tree.push_back(aNewNodeLeft);
    tree.push_back(aNewNodeRight);
    updateOmega(current_seg, normeTriplet-1); //Update omega of the parent
    //Test Omega
    if(current_seg->omega>0 && current_seg->omega<current_seg->alpha) {
      cout<<"----------> Node to be removed id="<<current_seg->idNode<<endl;
      addParent(aNewNodeLeft, current_seg->idParent, current_seg->p_parent);
      addParent(aNewNodeRight, current_seg->idParent, current_seg->p_parent);
    }
    //Insert node to the current level
    node* aNewNodeL = insertAfter(temp, aNewNodeLeft->idNode, aNewNodeLeft);
    node* aNewNodeR = insertAfter(aNewNodeL, aNewNodeRight->idNode, aNewNodeRight);
    deleteNode(head, temp);
    //Move to the next triplet if exists
    id_triplet++;
    if(id_triplet<triplets.size()) {
      triplet = triplets.at(id_triplet);
      //Update current seg with the new right segment
      temp = aNewNodeR;
      current_seg=temp->data.p_segment;
    }
  }

  //Processing of normal cases + right segment
  temp = temp->next;
  current_seg=temp->data.p_segment;
  
  while (id_triplet<triplets.size()) {
    vector<int> triplet = triplets.at(id_triplet);
    int normeTriplet = triplet[0]*triplet[0]+triplet[1]*triplet[1];
    if(compareHingeAngles(triplet, current_seg->tripletRight)==-1) {//triplet < current_seg->tripletRight
      //Insert new node in the segment

      //Insert triplet here
      vector<int> tripletLeft = current_seg->tripletLeft;
      vector<int> tripletRight = current_seg->tripletRight;
      vector<vector<int> > bAngles = current_seg->pAngles;
      for(size_t it=0; it<bAngles.size(); it++) {
        vector<int> pAngle = bAngles.at(it);
      }

      assert(compareHingeAngles(tripletLeft, triplet)==-1 && compareHingeAngles(triplet,tripletRight)==-1); //tripletLeft<triplet<tripletRight
      //Insert node to the tree
      nodeTree* aNewNodeLeft = createNodeTree(idNewNode, tripletLeft, triplet);
      updateAlpha(aNewNodeLeft, normeTriplet);//Update alpha left node
      addParent(aNewNodeLeft, current_seg->idNode, current_seg);
      idNewNode++;
      nodeTree* aNewNodeRight = createNodeTree(idNewNode, triplet, tripletRight);
      updateAlpha(aNewNodeRight, normeTriplet);//Update alpha of right node
      addParent(aNewNodeRight, current_seg->idNode, current_seg);
      idNewNode++;
      //Check bijective from pAngles of the parent node !
      isBij = bijectiveBound(&aNewNodeLeft, bAngles);
      updateBijectivity(aNewNodeLeft, isBij);
      isBij = bijectiveBound(&aNewNodeRight, bAngles);
      updateBijectivity(aNewNodeRight, isBij);
      //Push new nodes to the tree
      tree.push_back(aNewNodeLeft);
      tree.push_back(aNewNodeRight);
      updateOmega(current_seg, normeTriplet-1); //Update omega of the parent
      //Test Omega
      if(current_seg->omega>0 && current_seg->omega<current_seg->alpha) {
        cout<<"----------> Node to be removed id="<<current_seg->idNode<<endl;
        addParent(aNewNodeLeft, current_seg->idParent, current_seg->p_parent);
        addParent(aNewNodeRight, current_seg->idParent, current_seg->p_parent);
      }
      //Insert node to the current level
      node* aNewNodeL = insertAfter(temp, aNewNodeLeft->idNode, aNewNodeLeft);
      node* aNewNodeR = insertAfter(aNewNodeL, aNewNodeRight->idNode, aNewNodeRight);
      deleteNode(head, temp);

      //Move to the next triplet if exists
      id_triplet++;
      if(id_triplet<triplets.size()) {
        triplet = triplets.at(id_triplet);
        //Update current seg with the new right segment
        temp = aNewNodeR;
        current_seg=temp->data.p_segment;
      }
    }
    else {
      if(compareHingeAngles(triplet, current_seg->tripletRight)==0) {//triplet == current_seg->tripletRight
        id_triplet++;
        triplet = triplets.at(id_triplet);
        temp = temp->next;
        current_seg=temp->data.p_segment;
      }
      else{
        if(temp->next != NULL) {//triplet > current_seg->tripletRight and not the last segment
          temp = temp->next;
          current_seg=temp->data.p_segment;
          triplet = triplets.at(id_triplet);
        }
        else {//After the last segment
          assert(compareHingeAngles(triplet, current_seg->tripletRight)==1);//triplet > current_seg->tripletRight
          //Insert new node on the right
          //Insert triplet here
          current_seg=(*head)->data.p_segment; //=> loop node
          vector<int> tripletLeft = current_seg->tripletLeft;
          vector<int> tripletRight = current_seg->tripletRight;
          vector<vector<int> > bAngles = current_seg->pAngles;
          
          for(size_t it=0; it<bAngles.size(); it++) {
            vector<int> pAngle = bAngles.at(it);
          }
          
          //Insert node to the tree
          nodeTree* aNewNodeLeft = createNodeTree(idNewNode, tripletLeft, triplet);
          updateAlpha(aNewNodeLeft, normeTriplet);//Update alpha left node
          addParent(aNewNodeLeft, current_seg->idNode, current_seg);
          idNewNode++;
          nodeTree* aNewNodeRight = createNodeTree(idNewNode, triplet, tripletRight);
          updateAlpha(aNewNodeRight, normeTriplet);//Update alpha of right node
          addParent(aNewNodeRight, current_seg->idNode, current_seg);
          idNewNode++;
          //Check bijective from pAngles of the parent node !
          isBij = bijectiveBound(&aNewNodeLeft, bAngles);
          updateBijectivity(aNewNodeLeft, isBij);
          isBij = bijectiveBound(&aNewNodeRight, bAngles);
          updateBijectivity(aNewNodeRight, isBij);
          
          tree.push_back(aNewNodeLeft);
          tree.push_back(aNewNodeRight);
          updateOmega(current_seg, normeTriplet-1); //Update omega of the parent
          //Test Omega
          if(current_seg->omega>0 && current_seg->omega<current_seg->alpha) {
            cout<<"----------> Node to be removed id="<<current_seg->idNode<<endl;
            addParent(aNewNodeLeft, current_seg->idParent, current_seg->p_parent);
            addParent(aNewNodeRight, current_seg->idParent, current_seg->p_parent);
          }
          //Insert node to the current level
          node* aNewNodeL = insertAfter(temp, aNewNodeLeft->idNode, aNewNodeLeft);
          node* aNewNodeR = insertAfter(*head, aNewNodeRight->idNode, aNewNodeRight);
          deleteNode(head, *head);
          //Move to the next triplet if exists
          id_triplet++;
          if(id_triplet<triplets.size()) {
            triplet = triplets.at(id_triplet);
            //Update current seg with the new right segment
            temp = aNewNodeR;
            current_seg=temp->data.p_segment;
          }
        }
      }
    }
  }
  return idNewNode;
}
