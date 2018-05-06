/*
this code need Intel MKL to run arma
*/
#include<iostream>
using namespace std;
#include"lib/Reader.h"
#include"lib/ReaderHelper.h"
#include"lib/Spectral.h"
#include"lib/Louvain.h"
#include"lib/NMF.h"
char *dataset = "../GraphClustering/dataset/wisconsin/wisconsin_adj.txt";
//char *label = "../GraphClustering/dataset/texas/texas_label.txt";
char *output= "../GraphClustering/output.txt";
int main() {
	AroyaReader reader;
	reader.setSplit('\t');
	reader.read(dataset);
	AroyaReaderHelper readerHelper;
	readerHelper.insertAll(reader, false);
	
	//AroyaReader labelReader;
	//labelReader.read(label);
	//AroyaReaderHelper labelHelper;
	//labelHelper.insertAll(labelReader);
	//vector<vector<double>> labelTemp = labelHelper.getData();
	//vector<int> trueLabel;
	//for (int i = 0; i < labelTemp.size();i++)
	//	trueLabel.push_back(int(labelTemp[i][0]));

	Spectral spectral;
	spectral.setAdjData(readerHelper.getData());
	spectral.setSpectralType(NormalizedCut);
	spectral.setKEigenVectors(2);
	spectral.myans.setBord(0.01);
	vector<int> thisAns;
	thisAns = spectral.getFlag();
	spectral.myans.setClusters(6);
	spectral.run();
	writeFile(spectral.getFlag(), output);

	//Louvain louvain;
	//louvain.setAdjData(readerHelper.getData());
	//louvain.run();
	//vector<int> louvainAns = louvain.getLabel();
	//writeFile(louvainAns,output);

	//NMF nmf;
	//nmf.setK(5);
	//nmf.setAdjData(readerHelper.getData());
	//nmf.setLearnRate(0.01);
	//nmf.run();
	//writeFile(nmf.getLabel(),output);


	//system("pause");
	return 0;
}