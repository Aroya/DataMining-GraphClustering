#ifndef AROYA_NMF_C
#define AROYA_NMF_C
#include"NMF.h"


NMF::NMF(){
	learnRate = 0.1;
    K = 1;
    label = nullptr;
}
NMF::~NMF(){
    if(label!=nullptr)
        free(label);
}

void NMF::setK(const int &t) { K = t; }

void NMF::setAdjData(const vector<vector<double>> &newData)
{
    rows = newData.size();
    assert(rows > 0);
    cols = newData[0].size();
    A.resize(rows, cols);
    int i, j;
    for (i = 0; i < rows;i++)
    {
        for (j = 0; j < cols;j++)
        {
            A(i, j) = newData[i][j];
        }
    }
    U.resize(rows, K);
    V.resize(rows, K);
    U.setRandom();
    V.setRandom();
    // if(label!=nullptr)
    //     free(label);
    // label = (int *)malloc(sizeof(int) * rows);
}
void NMF::run(){
    MatrixXd U1, U2, V1, V2;
    int i, j;
    double Loss = 0;
	double Loss_before = -1;
    while (Loss > Loss_before)
    {
        U1 = A * V;
        U2 = U * V.transpose() * V;
        V1 = A.transpose() * U;
        V2 = V * U.transpose() * U;
        for (i = 0; i < rows;i++)
        {
            for (j = 0; j < K;j++)
            {
                U(i, j) += learnRate*U(i, j) * U1(i, j) / U2(i, j);
                V(i, j) += learnRate*V(i, j) * V1(i, j) / V2(i, j);
            }
        }
        Loss_before = Loss;
        Loss = (A - U * V.transpose()).norm();
		printf("Now Loss:%f\n", Loss);
    }
    vector<vector<double>> transfer;
    for (i = 0; i < rows;i++)
    {
        transfer.push_back(vector<double>());
        for (j = 0; j < K;j++)
            transfer[i].push_back(U(i, j));
    }
    kmeans.setClusters(K);
    kmeans.setData(transfer);
    kmeans.setBord(0.01);
    kmeans.run();
}
vector<int> NMF::getLabel() { return kmeans.getFlag(); }
void NMF::setLearnRate(const double&t) { learnRate = t; }
#endif