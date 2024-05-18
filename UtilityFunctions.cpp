#include "UtilityFunctions.h"

//Compute the norme of a 2D vector
int normVector(const vector<int>& v) {
  return v[0]*v[0]+v[1]*v[1];
}

/*______________________________________________________________*/
/*  Exact comparison between two triplets encoded hinges angles */
/*______________________________________________________________*/
long long int _4lambda2(int p, int q, int k) {
  return 4*(p*p+q*q)-(4*k*k+4*k+1);
}

int sign(int p, int q, int k) {
  long long int left = 2*q;
  long long int right = -1*p*(2*k+1);
  long long int left2 = q*q*_4lambda2(p,q,k);
  long long int right2 = right*right;
  
  //different signs
  if(left>=0 && right<0)
    return 1;
  if(left<=0 && right>0)
    return -1;
  if(left<0 && right>=0)
    return -1;
  if(left>0 && right<=0)
    return 1;
  
  // both null
  if(left==0 && right==0)
    return 0;
  
  //same sign
  if(left2==right2)
    return 0;
  
  if(left2>right2) {
    if(left>0) //both are positive
      return 1;
    else //both are negative
      return -1;
  }
  if(left2<right2) {
    if(left>0) //both are positive
      return -1;
    else //both are negative
      return 1;
  }
}

int signCosinus(int p, int q, int k) {
  return sign(p, q, k);
}

int signSinus(int p, int q, int k) {
  return sign(-1*q, p, k);
}

int is_Zero(int p, int q, int k) {
  return ((signSinus(p, q, k) == 0 && signCosinus(p, q, k) == 1));
}

int is_Pi(int p, int q, int k) {
  return ((signSinus(p, q, k) == 0 && signCosinus(p, q, k) == -1));
}

int is_Between_Zero_And_Pi(int p, int q, int k) {
  return (signSinus(p, q, k) > 0);
}

int is_Between_Pi_And_2Pi(int p, int q, int k) {
  return (signSinus(p, q, k) < 0);
}

int sign_RHS(int p, int q, int k, int r, int s, int l) {
  // returns
  // -1 if RHS negative
  // 0 if RHS null
  // 1 if RHS positive
  
  // one of the terms is null
  if(s==0 || q==0){
    if (s==0 && q==0)
      return 0;
    if (s > 0 || q < 0)
      return 1;
    return -1;
  }
  
  // no term is null
  if(s>0 && q<0)
    return 1;
  if(s<0 && q>0)
    return -1;
  
  // both terms have the same sign
  long long int RHS_pq_2=(s*(p*p + q*q))*(s*(p*p + q*q))*_4lambda2(r,s,l);
  long long int RHS_rs_2=(q*(r*r + s*s))*(q*(r*r + s*s))*_4lambda2(p,q,k);
  
  int RHS_pq_2_gt_RHS_rs_2 = 0;
  if(RHS_pq_2 > RHS_rs_2)
    RHS_pq_2_gt_RHS_rs_2 = 1;
  if(RHS_pq_2 < RHS_rs_2)
    RHS_pq_2_gt_RHS_rs_2 = -1;
  if(s>0)
    return RHS_pq_2_gt_RHS_rs_2;
  return -1*RHS_pq_2_gt_RHS_rs_2;
}

int compareCosinusHingeAngles(const vector<int>& a1, const vector<int>& a2) {
  // returns
  // -1 if cos a1 < cos a2
  // 0 if cos a1 == cos a2
  // 1 if cos a1 > cos a2
  
  //vector<long long int> alpha1; //(p,q,k)
  int p = a1[0];
  int q = a1[1];
  int k = a1[2];
  
  //vector<long long int> alpha2; //(r,s,l)
  int r = a2[0];
  int s = a2[1];
  int l = a2[2];
  
  int p2 = p*p;
  int q2 = q*q;
  long long int r2 = r*r;
  long long int s2 = s*s;
  long long int A = (r2 + s2)*p*(2*k + 1);
  long long int B = (p2 + q2)*r*(2*l + 1);
  long long int C = (r2 + s2)*q;
  long long int D = (p2 + q2)*s;
  
  long long int LHS_3_9 = A-B;
  int sign_RHS_3_9 = sign_RHS(p,q,k,r,s,l);
  int sign_LHS_3_9 = 0;
  if(LHS_3_9>0)
    sign_LHS_3_9=1;
  if(LHS_3_9<0)
    sign_LHS_3_9=-1;
  
  // LHS and RHS both null
  if(sign_LHS_3_9==0 && sign_RHS_3_9==0)
    return 0;
  
  // LHS and RHS of different signs, non both null
  if(sign_LHS_3_9<=0 && sign_RHS_3_9>0)
    return -1;
  if(sign_LHS_3_9>0 && sign_RHS_3_9<=0)
    return 1;
  if(sign_LHS_3_9<0 && sign_RHS_3_9>=0)
    return -1;
  if(sign_LHS_3_9>=0 && sign_RHS_3_9<0)
    return 1;
  
  // LHS and RHS of same sign, non null
  long long int LHS_3_11 = (A-B)*(A-B) - (C*C*_4lambda2(p,q,k) + D*D*_4lambda2(r,s,l));
  long long int RHS_3_11_partial = -1*C*D;
  int sign_3_9 = sign_LHS_3_9;
  
  // LHS_3_11 and RHS_3_11 both null
  if(LHS_3_11==0 && RHS_3_11_partial==0)
    return 0;
  
  // LHS_3_11 and RHS_3_11 of different signs, non both null
  if(LHS_3_11>0 && RHS_3_11_partial<=0)
    return sign_3_9;
  if(LHS_3_11<=0 && RHS_3_11_partial>0)
    return -1*sign_3_9;
  if(LHS_3_11<0 && RHS_3_11_partial>=0)
    return -1*sign_3_9;
  if(LHS_3_11>=0 && RHS_3_11_partial<0)
    return sign_3_9;
  
  // LHS_3_11 and RHS_3_11 of same sign, non null
  long long int LHS_3_12 = LHS_3_11*LHS_3_11;
  long long int RHS_3_12 = 4*_4lambda2(p,q,k)*_4lambda2(r,s,l)*RHS_3_11_partial*RHS_3_11_partial;
  
  int sign_3_11 = 1;
  if (LHS_3_11<0)
    sign_3_11 = -1;
  
  if (LHS_3_12 == RHS_3_12)
    return 0;
  if (LHS_3_12 > RHS_3_12)
    return sign_3_9*sign_3_11;
  if (LHS_3_12 < RHS_3_12)
    return -1*sign_3_9*sign_3_11;
}

// Compare triplets (Theorem 3.9 PHD Yohan, page 32)
int compareHingeAngles(const vector<int>& a1, const vector<int>& a2) {
  // returns
  // -1 if a1 < a2
  // 0 if a1 == a2
  // 1 if a1 > a2
  
  //vector<int> alpha1; //(p,q,k)
  int p = a1[0];
  int q = a1[1];
  int k = a1[2];
  
  //vector<int> alpha2; //(r,s,l)
  int r = a2[0];
  int s = a2[1];
  int l = a2[2];
  
  int a1_is_Zero = is_Zero(p,q,k);
  int a2_is_Zero = is_Zero(r,s,l);
  int a1_is_Pi = is_Pi(p,q,k);
  int a2_is_Pi = is_Pi(r,s,l);
  int a1_is_Between_Zero_And_Pi = is_Between_Zero_And_Pi(p,q,k);
  int a2_is_Between_Zero_And_Pi = is_Between_Zero_And_Pi(r,s,l);
  int a1_is_Between_Pi_And_2Pi = is_Between_Pi_And_2Pi(p,q,k);
  int a2_is_Between_Pi_And_2Pi = is_Between_Pi_And_2Pi(r,s,l);
  
  // one of the angles is 0
  if (a1_is_Zero || a2_is_Zero)
  {
    if (a1_is_Zero && a2_is_Zero)
      return 0;
    if (a1_is_Zero)
      return -1;
    return 1;
  }
  // no angle is 0 ; one of the angles is pi
  if (a1_is_Pi || a2_is_Pi)
  {
    if (a1_is_Pi && a2_is_Pi)
      return 0;
    if (a1_is_Between_Zero_And_Pi || a2_is_Between_Pi_And_2Pi)
      return -1;
    return 1;
  }
  // no angle is 0 ; no angle is pi ; one angle in each half space
  if ((a1_is_Between_Zero_And_Pi && a2_is_Between_Pi_And_2Pi) || (a1_is_Between_Pi_And_2Pi && a2_is_Between_Zero_And_Pi)){
    if (a1_is_Between_Pi_And_2Pi && a2_is_Between_Zero_And_Pi)
      return 1;
    return -1;
  }
  
  // no angle is 0 ; no angle is pi ; both angles in the same half space
  int compare_cos_a1_and_cos_2 = compareCosinusHingeAngles(a1,a2);
  if (a1_is_Between_Zero_And_Pi)
    return -1*compare_cos_a1_and_cos_2;
  return compare_cos_a1_and_cos_2;
  
}

/*______________________________________________________*/
/*  Exact comparison of a pythagore and an hinge angles */
/*______________________________________________________*/
int compareCosinusHingePythagoreAngles(const vector<int>& triplet, const vector<int>& pythagore) {
  // returns
  // -1 if cos a1 < cos a2
  // 0 if cos a1 == cos a2
  // 1 if cos a1 > cos a2
  
  int p=triplet[0];//alpha or a1
  int q=triplet[1];
  int k=triplet[2];
  int a=pythagore[0];//theta or a2
  int b=pythagore[1];
  int c=pythagore[2];
  
  long long int _4lambda2 = 4*(p*p + q*q) - (4*k*k+4*k+1);
  long long int LHS_3_6 = c*p*(2*k + 1) - 2*a*(p*p + q*q);
  long long int RHS_3_6_partiel = -2*c*q;
  if(LHS_3_6>0 && RHS_3_6_partiel<0) {
    return 1; //cos alpha > cos theta
  }
  else {
    if(LHS_3_6<0 && RHS_3_6_partiel>0) {
      return -1; //cos alpha < cos theta
    }
    else {
      //LHS_3_6_2 and RHS_3_6_2 of same sign, non null (3.6 -> 3.7)
      long long int LHS_3_6_2 = LHS_3_6*LHS_3_6;
      long long int RHS_3_6_2 = _4lambda2*c*c*q*q;
      if(LHS_3_6_2>RHS_3_6_2) {
        if(LHS_3_6>0) //both are positive in 3.6
          return 1;
        else //both are negative in 3.6
          return -1;
      }
      else {//if(LHS_3_6_2>RHS_3_6_2) {
        if(LHS_3_6>0) //both are positive in 3.6
          return -1;
        else //both are negative in 3.6
          return 1;
        //}
      }
    }
  }
}

int compareTripletPythagoreAngles(const vector<int>& triplet, const vector<int>& pythagore) {
  // triplet (alpha) vs pythagore (theta)
  // returns
  // -1 if alpha < theta
  // 0 if alpha == theta
  // 1 if alpha > theta
  
  int p=triplet[0];//alpha or a1
  int q=triplet[1];
  int k=triplet[2];
  int a=pythagore[0];//theta or a2
  int b=pythagore[1];
  int c=pythagore[2];
  
  int a1_is_Zero = is_Zero(p,q,k);
  int a2_is_Zero = a==1 ? 1 : 0;
  int a1_is_Pi = is_Pi(p,q,k);
  int a2_is_Pi = a==-1 ? 1 : 0;
  int a1_is_Between_Zero_And_Pi = is_Between_Zero_And_Pi(p,q,k);
  int a2_is_Between_Zero_And_Pi = b>0 ? 1 : 0;
  int a1_is_Between_Pi_And_2Pi = is_Between_Pi_And_2Pi(p,q,k);
  int a2_is_Between_Pi_And_2Pi = b<0 ? 1 : 0;
  
  // one of the angles is 0
  if (a1_is_Zero || a2_is_Zero) {
    if (a1_is_Zero && a2_is_Zero)
      return 0;
    if (a1_is_Zero)
      return -1;
    return 1;
  }
  // no angle is 0 ; one of the angles is pi
  if (a1_is_Pi || a2_is_Pi) {
    if (a1_is_Pi && a2_is_Pi)
      return 0;
    if (a1_is_Between_Zero_And_Pi || a2_is_Between_Pi_And_2Pi)
      return -1;
    return 1;
  }
  
  // no angle is 0 ; no angle is pi ; one angle in each half space
  if ((a1_is_Between_Zero_And_Pi && a2_is_Between_Pi_And_2Pi) || (a1_is_Between_Pi_And_2Pi && a2_is_Between_Zero_And_Pi)) {
    if (a1_is_Between_Pi_And_2Pi && a2_is_Between_Zero_And_Pi)
      return 1;
    return -1;
  }
  
  // no angle is 0 ; no angle is pi ; both angles in the same half space
  // compare by cos (alpha) - cos (theta) : Eq 3.6
  int compare_cos_a1_and_cos_2 = compareCosinusHingePythagoreAngles(triplet, pythagore);
  if (a1_is_Between_Zero_And_Pi)
    return -1*compare_cos_a1_and_cos_2;
  return compare_cos_a1_and_cos_2;
}

/*__________________________________________*/
/*  Other functions on angles for the tests */
/*__________________________________________*/
double computeAngle (vector<int> triplet) {
  int p = triplet[0];
  int q = triplet[1];
  int k = triplet[2];
  int p2 = p*p;
  int q2 = q*q;
  int k2 = k*k;
  double lambda = sqrt(p2 + q2 - (k+0.5)*(k+0.5));
  double cos_alpha = (q*lambda+p*(k+0.5))/(p2 + q2);//(p*lambda + q*(k+0.5)) / (p*p + q*q);
  double sin_alpha = (p*lambda-q*(k+0.5))/(p2 + q2);//(p*(k+0.5) - q*lambda) / (p*p + q*q);
  double alpha = angle_from_sin_cos(sin_alpha, cos_alpha);//acos(cos_alpha);
  return alpha;
}

double calculCosinus(long int p, long int q, long int k) {
  double lambda = sqrt(p*p+q*q - (k+0.5)*(k+0.5));
  return (q*lambda+p*(k+0.5))/(p*p+q*q);
}

double angle_from_sin_cos(double sinx, double cosx) //result in -pi to +pi range
{
  double ang_from_cos = acos(cosx);
  double ang_from_sin = asin(sinx);
  double sin2 = sinx*sinx;
  if(sinx<0)
  {
    ang_from_cos = -ang_from_cos;
    if(cosx<0) //both negative
      ang_from_sin = -M_PI -ang_from_sin;
  }
  else if(cosx<0)
    ang_from_sin = M_PI - ang_from_sin;
  //now favor the computation coming from the
  //smaller of sinx and cosx, as the smaller
  //the input value, the smaller the error
  double angle = (1.0-sin2)*ang_from_sin + sin2*ang_from_cos;
  //report angle into 0 to 2PI
  if(angle<0)
    angle = angle + 2*M_PI;
  
  return angle;
}

/*___________________________________________*/
/* Generate Pythagore bijective angles in Z2 */
/*___________________________________________*/
vector<vector<int> > generateTripletPythagoreAngles(int r) {
  vector<vector<int> > pythagoreAngle;
  //p=0
  vector<int> v = {0, 1, 1}; //PI/2
  pythagoreAngle.push_back(v);
  v = {1, 0, 1}; //0
  pythagoreAngle.push_back(v);
  v = {-1, 0, 1}; //PI
  pythagoreAngle.push_back(v);
  v = {0, -1, 1}; //3PI/2
  pythagoreAngle.push_back(v);
  //p!=0
  int p=1, a, b, c;
  while(p*p<4*r) {
    a = 2*p*(p+1);//cos = a/c
    b = 2*p+1;//sin = a/c
    c = 2*p*p+2*p+1;
    v = {a, b, c};//a,b,c
    pythagoreAngle.push_back(v);
    v = {b, a, c};//b,a, c
    pythagoreAngle.push_back(v);
    
    v = {-a, b, c};
    pythagoreAngle.push_back(v);
    v = {b, -a, c};
    pythagoreAngle.push_back(v);
    
    v = {a, -b, c};
    pythagoreAngle.push_back(v);
    v = {-b, a, c};
    pythagoreAngle.push_back(v);
    
    v = {-a, -b, c};
    pythagoreAngle.push_back(v);
    v = {-b, -a, c};
    pythagoreAngle.push_back(v);
    p++;
  }
  /*
   vector<double> pAngle;
   for(size_t it=0; it<pythagoreAngle.size(); it++) {
   int a = pythagoreAngle.at(it)[0];
   int b = pythagoreAngle.at(it)[1];
   int c = pythagoreAngle.at(it)[2];
   
   double angle = angle_from_sin_cos( double(b)/c, double(a)/c );
   pAngle.push_back(angle);
   }
   sort(pAngle.begin(), pAngle.end());
   cout << "Sorted \n";
   for (auto x : pAngle)
   cout << x << " " <<endl;
   */
  return pythagoreAngle;
}

