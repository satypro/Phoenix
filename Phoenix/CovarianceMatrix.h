#ifndef COVARIANCEMATRIX_H
#define COVARIANCEMATRIX_H

#include <stdio.h>
#include <vector>
#include "data_type.h"
#include <pcl/point_types.h>
#include <pcl/point_cloud.h>
#include <pcl/io/pcd_io.h>
#include "OctreeSerial.h"
#include <pcl/octree/octree_impl.h>
#include <pcl/octree/octree_iterator.h>
#include <pcl/octree/octree_pointcloud_pointvector.h>
#include <pcl/octree/octree_container.h>
#include <iostream>
#include <vector>
#include <ctime>
#include "SearchNeighbour.h"
#include "SearchStructure.h"
#include "SearchNeighbourFactory.h"
#include "Classifier.h"

using namespace std;
using namespace Eigen;

class CovarianceMatrix : public Classifier
{
public:
	typedef Eigen::Matrix<double, 3, 1> ColumnVector;
	typedef Eigen::Matrix<double, 1, 3> RowVector;
	typedef pcl::PointXYZ PointT;
	typedef pcl::octree::OctreeContainerPointIndices LeafContainerT;
	typedef pcl::octree::OctreeContainerEmpty BranchContainerT;

	CovarianceMatrix() 
	{
		search = SearchNeighbourFactory::GetNeighbourSearchDataStructure(SearchStructure::OctTree);
	}

	void setInputCloud(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud);

	void setOctree(OctreeXYZ<PointT, LeafContainerT, BranchContainerT > * octree);

	void setParams(float rmin, float rmax, float rmaxpt, float epsi, float scale);

	void prob_measure(vector <probfeatnode> *probval, vector <tensorType> *aggregated_tensor);
	//New Prob measure
	void prob_measure(probfeatnode *probval, tensorType *tensor, int index);

	void calculatePartialDerivative(float radius) {
		//Method not Implemented
	}

	void calculateSecondDerivative(float radius) {
		//Method not Implemented
	}

	void calculateThirdDerivative(float radius) {
		//Method not Implemented
	}

private:

	bool setAccumzero();

	bool voteAccum(float radius, vector<probfeatnode>&probfeatVal);

	bool probMeasure(vector<probfeatnode>&probfeatVal);

	// New Method
	bool voteAccum(float radius, probfeatnode probVal, int index);
	// New Method
	bool probMeasure(probfeatnode &probfeatVal, int index);

	void writeGlyphVars(float radius);

	void computeSaliencyVals(glyphVars &glyphtemp, int idx);

	void glyphAnalysis(glyphVars &glyphtemp);

	bool eigendecomposition(glyphVars &glyphtemp, size_t idx);

	float _rmin, _rmax, _rmaxpt, _epsi, _scale;

	float _sigma;

	pcl::PointCloud<pcl::PointXYZ>::Ptr _inCloud;

	OctreeXYZ<PointT, LeafContainerT, BranchContainerT > * _octree;

	vector < tensorType	> _accum;

	vector<glyphVars> _glyphData;

	vector < tensorType	> averaged_tensors;

	vector < tensorType	> averaged_accum;

	float max_cl, min_cl;

	SearchNeighbour* search;

	tensorType _n_accum;
	glyphVars _n_glyphData;
	tensorType _n_avg_acum;
	tensorType _n_avg_tensor;
};
#endif