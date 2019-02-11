#pragma once
#include <stdio.h>
#include <vector>
#include "OctreeSerial.h"
#include "data_type.h"
#include "SearchNeighbour.h"
#include "SearchStructure.h"
#include "SearchNeighbourFactory.h"
#include "ProcessRequest.h"

using namespace std;
using namespace Eigen;

class Classifier
{
public:
	Classifier();
	~Classifier();
	typedef Eigen::Matrix<double, 3, 1> ColumnVector;
	typedef Eigen::Matrix<double, 1, 3> RowVector;
	typedef pcl::PointXYZ PointT;
	typedef pcl::octree::OctreeContainerPointIndices LeafContainerT;
	typedef pcl::octree::OctreeContainerEmpty BranchContainerT;

	virtual void setInputCloud(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud) = 0;
	virtual void setOctree(OctreeXYZ<PointT, LeafContainerT, BranchContainerT > * octree) = 0;
	virtual void setParams(float rmin, float rmax, float rmaxpt, float epsi, float scale) = 0;
	virtual void setSearch(SearchNeighbour* search)
	{
		this->search = search;
	}

	virtual void process(ProcessRequest processRequest, vector <probfeatnode> *probval, vector <tensorType> *aggregated_tensor)
	{
		this->request = request;
	}

	virtual void processPoint(ProcessRequest request, probfeatnode *_probval, tensorType *_accum, int index)
	{
		this->request = request;
	}

	virtual void prob_measure(vector <probfeatnode> *probval) 
	{
	}

	virtual void prob_measure(vector <probfeatnode> *probval, vector <tensorType> *aggregated_tensor) 
	{
	}

	virtual void prob_measure(vector <probfeatnode> *probval, vector <tensorType> *aggregated_tensor, vector<int> optimalScale) 
	{
	}

	virtual void calculatePartialDerivative(float radius) = 0;
	virtual void calculateSecondDerivative(float radius) = 0;
	virtual void calculateThirdDerivative(float radius) = 0;

	SearchNeighbour* search;
	ProcessRequest request;
};