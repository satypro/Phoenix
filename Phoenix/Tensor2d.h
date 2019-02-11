#ifndef TENSOR2D_H
#define TENSOR2D_H

#include <stdio.h>
#include <vector>
#include "OctreeSerial.h"
#include "data_type.h"
#include "SearchNeighbour.h"
#include "SearchStructure.h"
#include "SearchNeighbourFactory.h"
#include "Classifier.h"

using namespace std;
using namespace Eigen;


class Tensor2d : public Classifier
{
public:
	Tensor2d() {}
	typedef Eigen::Matrix<double, 3, 1> ColumnVector;
	typedef Eigen::Matrix<double, 1, 3> RowVector;
	typedef pcl::PointXYZ PointT;
	typedef pcl::octree::OctreeContainerPointIndices LeafContainerT;
	typedef pcl::octree::OctreeContainerEmpty BranchContainerT;

	void setInputCloud(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud);

	void setOctree(OctreeXYZ<PointT, LeafContainerT, BranchContainerT > * octree);

	void setParams(float rmin, float rmax, float rmaxpt, float epsi, float scale);

	void prob_measure(vector <probfeatnode> *probval);

	void calculatePartialDerivative(float radius);

	void calculateSecondDerivative(float radius) {
		//Method not Implemented
	}

	void calculateThirdDerivative(float radius) {
		//Method not Implemented
	}

private:
	void getnormalizeTensor(size_t idx, float weight);

	bool setAccumzero();

	bool structureTensor2D(float radius);

	bool probMeasure(vector<probfeatnode>&probfeatVal, float radius);

	void writeGlyphVars(float radius);

	void computeSaliencyVals(glyphVars &glyphtemp, int i, float radius);

	void glyphAnalysis(glyphVars &glyphtemp);

	bool eigendecomposition(glyphVars &glyphtemp, size_t idx);

	void average_eigenvalues();

	float _rmin, _rmax, _rmaxpt, _epsi, _scale;

	float _sigma;

	pcl::PointCloud<pcl::PointXYZ>::Ptr _inCloud;

	OctreeXYZ<PointT, LeafContainerT, BranchContainerT > * _octree;

	vector < tensorType	> _accum;

	vector<glyphVars> _glyphData;
	vector < tensorType	> averaged_tensors;

	vector <float> Ix;
	vector <float> Iy;

	double average_mu1;
	double average_mu2;
	double max_mu1, max_mu2, min_mu1, min_mu2;
};

#endif // TENSOR2D_H
