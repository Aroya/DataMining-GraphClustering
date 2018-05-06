#ifndef AROYA_KMEANS_H
#define AROYA_KMEANS_H

#include<iostream>
#include<vector>
#include<cmath>
#include<fstream>
using namespace std;

//#define AROYA_DEBUG

class AroyaKMeans {
public:
	AroyaKMeans();
	void setClusters(const int&clusters);					//���þ������
	void setData(const vector<vector<double>>&yourData);	//��reader������Ϣ
	void run();												//���о���
	void setBord(const double&newBord);						//������������
	void writeFile(const char*fileName,const bool&withData=false);	//�������д���ļ�
	vector<int>getFlag();
private:
	vector<vector<double>>data;		//data
	int*cluster;					//ָ��data���ڵ�centre
	vector<vector<double>>centre;	//����
	double **distance;				//data�����ĵľ���
	int rows, columns, clusters;
	double bord;					//Ĭ��5%���µı仯ʱ������
};


#endif
#include"KMeans.cpp"