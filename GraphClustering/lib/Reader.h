#ifndef	AROYA_READER_H
#define AROYA_READER_H
#include<string>
#include<vector>
#include<fstream>
#include<iostream>
#include<sstream>

using namespace std;
void test();

class AroyaReader {
private:
	vector<vector<string>>data;		//ȫ����string�ݴ�
	stringstream internalSst;		//ת��ʹ��
	int rows, columns;
	char split;						//default:csv:,
public:
	AroyaReader();
	void read(const char*fileName);
	string getStringData(const int&rows, const int&columns);
	double getDoubleData(const int&rows, const int&columns);
	int findTable(const char*tableName);
	int getRows();
	int getColumns();
	void setSplit(const char&);
};

#endif
#include"Reader.cpp"