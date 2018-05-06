#ifndef AROYA_LOUVAIN_H
#define AROYA_LOUVAIN_H

#include<vector>
using namespace std;

class Louvain{
    public:
      Louvain();
      ~Louvain();
      void setAdjData(const vector<vector<double>> &);
      void run();
      vector<int> getLabel();

    private:
      double **W, *Degrees, *WeightInnerCommunity, *WeightCommunity;
      int rows, cols, m;
      int *label;
      vector<vector<int>> communityNodes;
      void combineCommunity(const int &CommunityA, const int &CommunityB);
      double modularityGain(const int&newNode,const int&Community);
};

#endif
#include"Louvain.cpp"