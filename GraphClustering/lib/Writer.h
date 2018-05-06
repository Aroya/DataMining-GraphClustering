#ifndef AROYA_WRITER_H
#define AROYA_WRITER_H
#include<vector>
using namespace std;
void writeFile(const vector<vector<double>>&data,
	const vector<int>&flag,
	const char*dfilename);
void writeFile(const vector<int>&flag,
	const char*dfilename);
#endif // AROYA_WRITER
#include"Writer.cpp"