#ifndef AROYA_NMF_H
#define AROYA_NMF_H

#include<vector>
#include<cmath>
using namespace std;

#include"eigen3/Eigen/Dense"
using namespace Eigen;

#include"KMeans.h"

class NMF{
    public:
      NMF();
      ~NMF();
      void setAdjData(const vector<vector<double>> &);
      void setK(const int &);
      void run();
      vector<int> getLabel();
      AroyaKMeans kmeans;
	  void setLearnRate(const double&);

    private:
      int rows, cols, K;
      MatrixXd A, U, V;
      int *label;
	  double learnRate;
};

#endif
#include"NMF.cpp"