#pragma once
#include <pcl/io/pcd_io.h>
#include <string>
#include <stdio.h>
#include <vector>
#include "SearchNeighbour.h"
#include "SearchNeighbourFactory.h"
#include "SearchStructure.h"
#include "ProcessRequest.h"

class Processing
{
public:
	typedef pcl::PointCloud<pcl::PointXYZ> cl;
	Processing();
	~Processing();
	void buildStructure(ProcessRequest request);
	void processCloud(ProcessRequest request);
	void processPoint(ProcessRequest request, int index);
	bool fileRead(std::string cloudFile, std::string labelFile, std::string optimalFile);
	pcl::PointCloud<pcl::PointXYZ>::Ptr getCloud();

private:
	pcl::PointCloud<pcl::PointXYZ>::Ptr _inCloud;
	void scale();
	std::string _filePath, _FileExtension;
	std::vector <float> _intensity;
	std::vector <int> optimalScale;
	std::vector <int> labels;
	SearchNeighbour* search;
};