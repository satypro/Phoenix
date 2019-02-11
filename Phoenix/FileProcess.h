#ifndef FILE_READ_H
#define FILE_READ_H

#include <stdio.h>
#include <string.h>
#include <vector>
#include <iostream>
#include "data_type.h"
#include <pcl/io/pcd_io.h>


using namespace std;

/*****************************************************************************************

Class Name		:	FileRead
Purpose			:	read the input file and store data(x,y,z coordinates and intensity values) in the vectors data structure
Created by		:	Beena
Created date		:
Modification		:
Modified Function	:

*****************************************************************************************/
class FileProcess
{
public:
	typedef pcl::PointCloud<pcl::PointXYZ>::Ptr CloudType;

	FileProcess() {}
	~FileProcess() {}


	string getFileType(std::string sFilePath);
	bool read(string sFilePath, CloudType &cloudData, vector <float>&intensity);
	bool readLabels(string filename, vector <int>&labels);
	void displayLasInfo();

private:
	bool lasFileRead(const string sFilePath, CloudType &cloudData, vector <float>&intensity);
	bool plyFileRead(const string sFilePath, CloudType &cloudData, vector <float>&intensity);
	bool readOriginalPtCloud(string filename, CloudType &cloudData, vector <float>&intensity);
	bool readLasFileHeaderInfo(FILE *fp);
	bool offFileRead(string filename, CloudType &cloudData, vector <float>&intensity);

	bool xyzFileRead(string filename, CloudType &cloudData, vector <float>&intensity);
	SLasHeaderBlock _stPubBlk;

};


#endif
