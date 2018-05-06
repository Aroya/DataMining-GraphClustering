#ifndef AROYA_LOUVAIN_C
#define AROYA_LOUVAIN_C

#include"Louvain.h"

Louvain::Louvain(){
    W = nullptr;
    rows = 0;
    cols = 0;
    label = nullptr;
    Degrees = nullptr;
    WeightCommunity = nullptr;
    WeightInnerCommunity = nullptr;
}
Louvain::~Louvain(){
    if(W!=nullptr){
        for (int i = 0; i < rows;i++)
            free(W[i]);
        free(W);
    }
    //label
    if(label!=nullptr)
        free(label);
    //Degree Array
    if(Degrees!=nullptr)
        free(Degrees);
    if(WeightInnerCommunity!=nullptr)
        free(WeightInnerCommunity);
    if(WeightCommunity!=nullptr)
        free(WeightCommunity);
}
void Louvain::setAdjData(const vector<vector<double>>&newData){
    int i, j;
    double sum;
    //space process
    if(W!=nullptr){
        for (i = 0; i < rows;i++)
            free(W[i]);
        free(W);
    }
    rows = newData.size();
    cols = newData[0].size();
    W = (double **)malloc(sizeof(double *) * rows);
    //label
    if(label!=nullptr)
        free(label);
    //Degree Array
    if(Degrees!=nullptr)
        free(Degrees);
    if(!communityNodes.empty()){
        communityNodes = vector<vector<int>>();
    }
    if(WeightInnerCommunity!=nullptr)
        free(WeightInnerCommunity);
    if(WeightCommunity!=nullptr)
        free(WeightCommunity);
    label = (int *)malloc(sizeof(int) * rows);
    Degrees = (double *)malloc(sizeof(double) * rows);
    for (i = 0; i < rows;i++)
        communityNodes.push_back(vector<int>());
    WeightInnerCommunity = (double *)malloc(sizeof(double) * rows);
    WeightCommunity = (double *)malloc(sizeof(double) * rows);
    m = 0;
    for (i = 0; i < rows;i++){
        
        W[i] = (double *)malloc(sizeof(double) * cols);
        //save data into the space
        sum = 0;
        for (j = 0; j < cols;j++)
        {
            sum += (W[i][j] = newData[i][j]);
        }
        m += sum;
        Degrees[i] = sum;

        //Community Informations
        label[i] = i;
        communityNodes[i].push_back(i);
        WeightInnerCommunity[i] = 0;
        WeightCommunity[i] = sum;
    }
}

double Louvain::modularityGain(const int&newNode,const int&Community)
{
    double gain, tempDBL;
    int index, i, j, length;
    //Sum(IN)=Sum(Node to other nodes in community)+Sum(IN before)
    tempDBL = 0;
    length = communityNodes[Community].size();
    for (i = 0; i < length;i++){
        tempDBL += W[newNode][communityNodes[Community][i]];
    }
    double Sum_IN_ = tempDBL + WeightInnerCommunity[Community];

    //Sum(TOT)=Degree(NewNode)+Degree(Community before)
    double Sum_TOT_ = Degrees[newNode] + WeightCommunity[Community];
    
    //k(newNode)=Degree(newNode)
    double k_I = Degrees[newNode];

    //k(newNode;IN)=Sum(Node to other nodes in Community)
    double k_I_IN = tempDBL;

    //m has been already storaged

    //caculate gain
    gain = (Sum_IN_ + k_I_IN) / m - pow((Sum_TOT_ + k_I) / m, 2) / 4
    - Sum_IN_ / 2 / m - pow(Sum_TOT_ / m, 2) / 4 - pow(k_I / m, 2) / 4;
    return gain;
}
void Louvain::run(){
    //when combine 2 community
    //gain = Sum(i in Community1;i to Community2)
    int i, j, k, length, thisLabel, targetCommunity, maxINT;
    double tempDBL, maxDBL;
    for (i = 0; i < rows;i++){
        thisLabel = label[i];
        //judge near points
        tempDBL = 0;
        maxINT = -1;
		maxDBL = 0;
        for (j = 0; j < rows;j++){
            //find different community points
            if(W[i][j]>0&&thisLabel!=label[j])
            {
                targetCommunity = label[j];
                //caculate the gain if combines two community
                //Node → Community
                length = communityNodes[thisLabel].size();
                tempDBL = 0;
                for (k = 0; k < length;k++)
                {
                    tempDBL += modularityGain(communityNodes[thisLabel][k], targetCommunity);
                }
                if(tempDBL>maxDBL)
                {
                    maxDBL = tempDBL;
                    maxINT = j;
                }
            }
        }
        //maximum gain is positive, than combine 2 community
        if(maxDBL>0){
            //combineCommunity
            combineCommunity(label[i], label[maxINT]);
            i = 0;
        }
    }

    //normalize label
    maxINT = -1;
    bool *alive = (bool *)malloc(sizeof(bool) * rows);
    for (i = 0; i < rows;i++)
    {
        alive[i] = false;
    }
    for (i = 0; i < rows;i++)
    {
        alive[label[i]] = true;
    }
    int *newMap = (int *)malloc(sizeof(int) * rows);
    j = 0;
    for (i = 0; i < rows;i++){
        if(alive[i])
            newMap[i] = j++;
    }
    for (i = 0; i < rows;i++)
    {
        label[i] = newMap[label[i]];
    }
    free(alive);
    free(newMap);
}

void Louvain::combineCommunity(const int &CommunityA, const int &CommunityB){
    //combine B to A
    double tempDBL;
    int lengthA = communityNodes[CommunityA].size();
    int lengthB = communityNodes[CommunityB].size();
    int i, j;
    tempDBL = 0;
    for (i = 0; i < lengthB;i++){
        for (j = 0; j < lengthA;j++){
            tempDBL += W[communityNodes[CommunityA][j]][communityNodes[CommunityB][i]];
        }
        //change the community sign
		communityNodes[CommunityA].push_back(communityNodes[CommunityB][i]);
        label[communityNodes[CommunityB][i]] = CommunityA;
    }
    WeightInnerCommunity[CommunityA] += tempDBL + WeightInnerCommunity[CommunityB];
    WeightCommunity[CommunityA] += WeightCommunity[CommunityB];
}

vector<int> Louvain::getLabel(){
    vector<int>temp;
    for (int i = 0; i < rows;i++)temp.push_back(label[i]);
    return temp;
}

#endif