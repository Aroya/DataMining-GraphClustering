#ifndef AROYA_READER_HELPER_C
#define AROYA_READER_HELPER_C
#include"ReaderHelper.h"

AroyaReaderHelper::AroyaReaderHelper() {
	//nothing to do yet
}

//insert table
void AroyaReaderHelper::insert(AroyaReader&reader, const char*tableName, const char*myTableName) {
	if (myTableName == nullptr)myTableName = tableName;
	int t = reader.findTable(tableName);	//find table index
	int length = reader.getRows();			//get data size
	int myPosition = findTable(myTableName);//find myTableName index

	if (myPosition < 0) {					//not found create new table
		myPosition = buffer.size();			//point to end of buffer
		buffer.push_back(empty);			//new buffer
		table.push_back(myTableName);		//new table
	}

	for (int i = 1; i < length; i++) {		//push to buffer
		buffer[myPosition].push_back(reader.getDoubleData(i, t));
	}
}

void AroyaReaderHelper::insert(AroyaReader&reader, const int&col,const bool&exist_table) {
	//table����
	if(exist_table) table.push_back(reader.getStringData(0, col));
	else table.push_back("N/A");
	int i, l = reader.getRows(), j = buffer.size();
	//���д���
	buffer.push_back(empty);
	for (i = 1 - int(!exist_table); i < l; i++) {
		buffer[j].push_back(reader.getDoubleData(i, col));
	}
}

void AroyaReaderHelper::insertAll(AroyaReader&reader, const bool&exist_table) {
	int cols = reader.getColumns();
	for (int i = 0; i < cols; i++) {
		insert(reader, i, exist_table);
	}
}
int AroyaReaderHelper::findTable(const char*tableName) {
	int i, j;
	i = table.size();
	for (j = 0; j < i; j++) {
		if (table[j] == tableName)return j;
	}
	return -1;
}

void AroyaReaderHelper::writeFile(const char*fileName) {
	vector<vector<double>>ans = getData();
	ofstream fout;
	fout.open(fileName);
	if (fout.is_open()) {
		int i, j, rows, cols;
		bool notFirst;
		cols = table.size();
		rows = ans.size();
		notFirst = false;
		for (i = 0; i < cols; i++) {
			if (notFirst) {
				fout << ',' << table[i];
			}
			else {
				fout << table[i];
				notFirst = true;
			}
		}
		fout << endl;
		for (i = 0; i < rows; i++) {
			notFirst = false;
			for (j = 0; j < cols; j++) {
				if (notFirst) {
					fout << ',' << ans[i][j];
				}
				else {
					fout << ans[i][j];
					notFirst = true;
				}
			}
			fout << endl;
		}
	}
	else {
		printf("AroyaReaderHelper::writeFile()cannot write %s\n",fileName);
	}
}

int AroyaReaderHelper::getRows() {
	return buffer[0].size();
}
int AroyaReaderHelper::getColumns() {
	return table.size();
}

vector<vector<double>> AroyaReaderHelper::getData() {
	vector<vector<double>>temp;
	vector<double>item;
	int i, j, rows, columns;
	rows = buffer[0].size();
	columns = buffer.size();
	for (i = 0; i < rows; i++) {
		item.clear();
		for (j = 0; j < columns; j++) {
			item.push_back(buffer[j][i]);
		}
		temp.push_back(item);
	}
	return temp;
}

double** AroyaReaderHelper::getDataPointer() {
	double**temp;
	int i, j, rows, columns;
	rows = buffer[0].size();
	columns = buffer.size();
	temp = new double*[rows];
	for (i = 0; i < rows; i++) {
		temp[i] = new double[columns];
		for (j = 0; j < columns; j++) {
			temp[i][j] = buffer[j][i];
		}
	}
	return temp;
}

void AroyaReaderHelper::saveTable(const char*fileName) {
	ofstream fout;
	fout.open(fileName);
	if (fout.is_open()) {
		int i, j = getColumns();
		fout << table[0];
		for (i = 1; i < j; i++) {
			fout << ',' << table[i];
		}
		fout << endl;
		fout.close();
	}
	else {
		printf("AroyaReaderHelper::saveTable() cannot open %s\n", fileName);
	}
}
//����Ϊ��¼��table˳��
void AroyaReaderHelper::transformTable(const char*fileName) {
	ifstream fin;
	fin.open(fileName);
	if (fin.is_open()) {
		AroyaReader reader;
		reader.read(fileName);

		int i, j, k, l, p;
		j = reader.getColumns();
		vector<string>table_;
		for (i = 0; i < j; i++)table_.push_back(reader.getStringData(0, i));
		
		l = getRows();
		vector<vector<double>>buffer_;
		for (i = 0; i < j; i++) {
			buffer_.push_back(empty);
			for (k = 0; k < l; k++) {
				buffer_[i].push_back(0);
			}
		}
		for (i = 0; i < j; i++) {
			p = findTable(table_[i].c_str());
			if (p > -1) {
				for (k = 0; k < l; k++) {
					buffer_[i][k] = buffer[p][k];
				}
			}
		}
		buffer = buffer_;
		table = table_;
	}
	else {
		printf("AroyaReaderHelper::transformTable() cannot open %s\n", fileName);
	}
}
#endif