#include "DiffusedNormalVoting.h"
#include <algorithm>
#include <vector>
#include <pcl/common/common.h>
#include <pcl/common/eigen.h>
#include <pcl/common/centroid.h>
#include "eig3.h"

float current_radius = 0.0;
static int fNdx = 0;

void DiffusedNormalVoting::setInputCloud(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud)
{
	_inCloud = cloud;
}

void DiffusedNormalVoting::setOctree(OctreeXYZ<PointT, LeafContainerT, BranchContainerT > * octree)
{
	_octree = octree;
}

void DiffusedNormalVoting::setParams(float rmin, float rmax, float rmaxpt, float epsi, float scale)
{
	_rmin = rmin;
	_rmax = rmax;
	_rmaxpt = rmaxpt;
	_epsi = epsi;
	_scale = scale;

	_sigma = 1.0;
	_delta = 0.16;
}

bool DiffusedNormalVoting::setAccumzero()
{
	if (_inCloud->points.size() == 0)
		return false;

	_accum.resize(_inCloud->points.size());
	averaged_tensors.resize(_inCloud->points.size());
	averaged_accum.resize(_inCloud->points.size());

	for (size_t i = 0; i < _inCloud->points.size(); i++)
	{
		for (int j = 0; j < 3; j++)
		{
			_accum[i].evec0[j] = 0.0;
			_accum[i].evec1[j] = 0.0;
			_accum[i].evec2[j] = 0.0;

			averaged_tensors[i].evec0[j] = 0.0;
			averaged_tensors[i].evec1[j] = 0.0;
			averaged_tensors[i].evec2[j] = 0.0;

			averaged_accum[i].evec0[j] = 0.0;
			averaged_accum[i].evec1[j] = 0.0;
			averaged_accum[i].evec2[j] = 0.0;
		}
	}

	return true;
}

tensorType DiffusedNormalVoting::compute3DBallVote(ColumnVector V, float *weight)
{
	double norm, coeff, s, t;
	Matrix3d vv(3, 3), voteTemp(3, 3);
	tensorType ballTensorVote;
	Matrix3d I(3, 3);

	voteTemp.setZero();
	vv.setZero();
	I.setIdentity();

	s = V.norm(); // if s is zero

	t = (s*s) / (_sigma * _sigma);

	coeff = exp(-1.0 * t);

	*weight += coeff;

	if (V.norm() != 0.0)
		V = V.normalized();

	norm = V.transpose() * V; //if norm is zero
	norm = abs(norm);

	if (norm != 0.0)
	{
		vv = V * V.transpose();
		vv = vv / norm;  // if norm is zero
	}


	voteTemp = coeff * (I - vv);

	ballTensorVote.evec0[0] = voteTemp(0, 0);
	ballTensorVote.evec0[1] = voteTemp(0, 1);
	ballTensorVote.evec0[2] = voteTemp(0, 2);

	ballTensorVote.evec1[0] = voteTemp(1, 0);
	ballTensorVote.evec1[1] = voteTemp(1, 1);
	ballTensorVote.evec1[2] = voteTemp(1, 2);

	ballTensorVote.evec2[0] = voteTemp(2, 0);
	ballTensorVote.evec2[1] = voteTemp(2, 1);
	ballTensorVote.evec2[2] = voteTemp(2, 2);

	return ballTensorVote;
}

bool DiffusedNormalVoting::addBallVote(size_t idx, tensorType ballTensorVote, float lambdaN)
{

	for (int j = 0; j < 3; j++)
	{
		_accum[idx].evec0[j] = _accum[idx].evec0[j] + lambdaN * ballTensorVote.evec0[j];
		_accum[idx].evec1[j] = _accum[idx].evec1[j] + lambdaN * ballTensorVote.evec1[j];
		_accum[idx].evec2[j] = _accum[idx].evec2[j] + lambdaN * ballTensorVote.evec2[j];
	}

	return true;

}

bool DiffusedNormalVoting::makeVector(pcl::PointXYZ p_i, pcl::PointXYZ p_j, ColumnVector *V)
{
	ColumnVector temp;

	temp(0, 0) = double((p_j.x - p_i.x));
	temp(1, 0) = double((p_j.y - p_i.y));
	temp(2, 0) = double((p_j.z - p_i.z));

	double len = sqrt(temp(0, 0) * temp(0, 0) + temp(1, 0) * temp(1, 0) + temp(2, 0) * temp(2, 0));

	if (len == 0.0)
		return false;

	(*V)(0, 0) = temp(0, 0);
	(*V)(1, 0) = temp(1, 0);
	(*V)(2, 0) = temp(2, 0);

	return true;
}

void DiffusedNormalVoting::getnormalizeTensor(size_t idx, float weight)
{
	if (weight == 0.0)
		return;

	for (int j = 0; j < 3; j++)
	{
		_accum[idx].evec0[j] = _accum[idx].evec0[j] / weight;
		_accum[idx].evec1[j] = _accum[idx].evec1[j] / weight;
		_accum[idx].evec2[j] = _accum[idx].evec2[j] / weight;
	}

	for (int j = 0; j < 3; j++)
	{
		_n_accum.evec0[j] = _n_accum.evec0[j] / weight;
		_n_accum.evec1[j] = _n_accum.evec1[j] / weight;
		_n_accum.evec2[j] = _n_accum.evec2[j] / weight;
	}
}

bool DiffusedNormalVoting::voteAccum(float radius, vector<probfeatnode>&probfeatVal)
{
	std::vector<int> pointIdxRadiusSearch;
	std::vector<float> pointRadiusSquaredDistance;
	pcl::PointXYZ searchPoint, tempPoint;
	tensorType ballTensorVote;
	ColumnVector clVec;
	float weight;

	int min_neighbour = 1000;
	int max_neighbour = 0;
	int sum_neighbour = 0;

	for (size_t i = 0; i < _inCloud->points.size(); i++)
	{
		searchPoint = _inCloud->points[i];

		//        probfeatVal[i].neighborhoodSize = 0;

		if (search->radiusSearch(searchPoint, radius, pointIdxRadiusSearch, pointRadiusSquaredDistance) > 0)
		{
			weight = 0.0;

			if (pointIdxRadiusSearch.size() < min_neighbour)
				min_neighbour = pointIdxRadiusSearch.size();
			if (pointIdxRadiusSearch.size() > max_neighbour)
				max_neighbour = pointIdxRadiusSearch.size();
			sum_neighbour += pointIdxRadiusSearch.size();
			//            probfeatVal[i].neighborhoodSize = pointIdxRadiusSearch.size ();

			for (size_t j = 0; j < pointIdxRadiusSearch.size(); j++)
			{
				tempPoint = _inCloud->points[pointIdxRadiusSearch[j]];

				if (makeVector(searchPoint, tempPoint, &clVec))
				{
					ballTensorVote = compute3DBallVote(clVec, &weight);
					addBallVote(i, ballTensorVote, 1.0);   // initially consider eigenvalue = 1.0

				}
			}

			getnormalizeTensor(i, weight);

			pointIdxRadiusSearch.clear();
			pointRadiusSquaredDistance.clear();
		}
		else
			min_neighbour = 0;
	}

	std::cout << "Min neighbour = " << min_neighbour << " Max neighbour = " << max_neighbour << " average = " << sum_neighbour / _inCloud->points.size() << std::endl;

	return true;
}

bool DiffusedNormalVoting::voteAccum(float radius, probfeatnode probVal, int index)
{
	std::vector<int> pointIdxRadiusSearch;
	std::vector<float> pointRadiusSquaredDistance;
	pcl::PointXYZ searchPoint, tempPoint;
	tensorType ballTensorVote;
	ColumnVector clVec;
	float weight;

	int min_neighbour = 1000;
	int max_neighbour = 0;
	int sum_neighbour = 0;

	searchPoint = _inCloud->points[index];

	if (search->radiusSearch(searchPoint, radius, pointIdxRadiusSearch, pointRadiusSquaredDistance) > 0)
	{
		weight = 0.0;

		if (pointIdxRadiusSearch.size() < min_neighbour)
			min_neighbour = pointIdxRadiusSearch.size();
		if (pointIdxRadiusSearch.size() > max_neighbour)
			max_neighbour = pointIdxRadiusSearch.size();
		
		sum_neighbour += pointIdxRadiusSearch.size();

		for (size_t j = 0; j < pointIdxRadiusSearch.size(); j++)
		{
			tempPoint = _inCloud->points[pointIdxRadiusSearch[j]];

			if (makeVector(searchPoint, tempPoint, &clVec))
			{
					ballTensorVote = compute3DBallVote(clVec, &weight);
					addBallVote(index, ballTensorVote, 1.0);
			}
		}

		getnormalizeTensor(index, weight);

		pointIdxRadiusSearch.clear();
		pointRadiusSquaredDistance.clear();
	}
	else
	{
		min_neighbour = 0;
	}

	std::cout << " Min neighbour = "  << min_neighbour 
			  << " Max neighbour = "  << max_neighbour 
			  << " average = " << sum_neighbour / _inCloud->points.size() 
			  << std::endl;

	return true;
}

void DiffusedNormalVoting::processPoint(ProcessRequest request, probfeatnode *_probval, tensorType *_accum,  int index)
{
	_rmin = request.rmin;
	_rmax = request.rmax;
	_rmaxpt = request.rmaxpt;
	_epsi = request.epsi;
	_scale = request.scale;

	_sigma = 1.0;
	_delta = 0.16;

	prob_measure(_probval, _accum, index);
}

void DiffusedNormalVoting::process(ProcessRequest request, 
	vector <probfeatnode> *_probval, 
	vector <tensorType> *_aggregated_tensor)
{
	//vector <probfeatnode> _probval;
	//vector <tensorType> _accum; // aggregated tensor

	_rmin = request.rmin;
	_rmax = request.rmax;
	_rmaxpt = request.rmaxpt;
	_epsi = request.epsi;
	_scale = request.scale;

	_sigma = 1.0;
	_delta = 0.16;

	prob_measure(_probval, _aggregated_tensor);
}

void DiffusedNormalVoting::prob_measure(vector <probfeatnode> *probval, vector<tensorType> *aggregated_tensor)
{
	double radius, dDeltaRadius;

	if (_scale == 0.0 || _rmin == 0.0 || _rmax == 0.0 || _rmin >= _rmax || _inCloud->points.size() == 0)
	{
		cout << "invalid configuration parameters for classification module" << endl;
		return;
	}

	dDeltaRadius = (_rmax - _rmin) / (_scale - 1.0);
	radius = _rmin;

	vector<probfeatnode> probfeatVal;

	probfeatVal.resize(_inCloud->points.size());

	for (size_t i = 0; i < _inCloud->points.size(); i++)
	{
		probfeatVal[i].prob[0] = 0;
		probfeatVal[i].prob[1] = 0;
		probfeatVal[i].prob[2] = 0;

		probfeatVal[i].featStrength[0] = 0;
		probfeatVal[i].featStrength[1] = 0;
		probfeatVal[i].featStrength[2] = 0;

		probfeatVal[i].csclcp[0] = 0;
		probfeatVal[i].csclcp[1] = 0;
		probfeatVal[i].csclcp[2] = 0;
	}

	while (radius <= _rmax)
	{
		cout << "radius " << radius << " start (scale = " << _scale << endl;
		current_radius = radius;
		_sigma = radius;
		setAccumzero();
		voteAccum(radius, probfeatVal);

		for (size_t idx = 0; idx < _inCloud->points.size(); idx++)
		{
			averaged_accum[idx].evec0[0] += _accum[idx].evec0[0];
			averaged_accum[idx].evec0[1] += _accum[idx].evec0[1];
			averaged_accum[idx].evec0[2] += _accum[idx].evec0[2];

			averaged_accum[idx].evec1[0] += _accum[idx].evec1[0];
			averaged_accum[idx].evec1[1] += _accum[idx].evec1[1];
			averaged_accum[idx].evec1[2] += _accum[idx].evec1[2];

			averaged_accum[idx].evec2[0] += _accum[idx].evec2[0];
			averaged_accum[idx].evec2[1] += _accum[idx].evec2[1];
			averaged_accum[idx].evec2[2] += _accum[idx].evec2[2];
		}

		//probMeasure(probfeatVal,radius);
		radius += dDeltaRadius;
		cout << "radius " << radius - dDeltaRadius << " endl " << endl;
	}
	_scale = 1.0; //SATENDRA COMMENT
	// Single tensor Implementation
	for (size_t i = 0; i < _inCloud->points.size(); i++)
	{
		probfeatVal[i].csclcp[0] = 0;
		probfeatVal[i].csclcp[1] = 0;
		probfeatVal[i].csclcp[2] = 0;
	}

	_accum.swap(averaged_accum);
	probMeasure(probfeatVal, _rmax); // SATENDRA COMMENT

//   ground_density(probfeatVal);
	//writeGlyphVars(radius - dDeltaRadius);  SATENDRA COMMENT

	for (size_t i = 0; i < _inCloud->points.size(); i++)
	{
		probfeatVal[i].prob[0] = probfeatVal[i].prob[0] / _scale;
		probfeatVal[i].prob[1] = probfeatVal[i].prob[1] / _scale;
		probfeatVal[i].prob[2] = probfeatVal[i].prob[2] / _scale;

		probfeatVal[i].featStrength[0] = probfeatVal[i].featStrength[0] / _scale;
		probfeatVal[i].featStrength[1] = probfeatVal[i].featStrength[1] / _scale;
		probfeatVal[i].featStrength[2] = probfeatVal[i].featStrength[2] / _scale;

		probfeatVal[i].csclcp[0] = probfeatVal[i].csclcp[0] / _scale;
		probfeatVal[i].csclcp[1] = probfeatVal[i].csclcp[1] / _scale;
		probfeatVal[i].csclcp[2] = probfeatVal[i].csclcp[2] / _scale;

		//        probfeatVal[i].sum_eigen = probfeatVal[i].sum_eigen/_scale;
		probfeatVal[i].planarity = probfeatVal[i].planarity / _scale;
		probfeatVal[i].anisotropy = probfeatVal[i].anisotropy / _scale;
		probfeatVal[i].sphericity = probfeatVal[i].sphericity / _scale;
		probfeatVal[i].linearity = probfeatVal[i].linearity / _scale;
		probfeatVal[i].omnivariance = probfeatVal[i].omnivariance / _scale;
		probfeatVal[i].eigenentropy = probfeatVal[i].eigenentropy / _scale;
	}

	probval->resize(_inCloud->points.size());
	for (size_t i = 0; i < _inCloud->points.size(); i++)
		(*probval)[i] = probfeatVal[i];

	aggregated_tensor->resize(_inCloud->points.size());
	for (size_t i = 0; i < _inCloud->points.size(); i++)
		(*aggregated_tensor)[i] = _accum[i];

	probfeatVal.clear();
	return;
}

void DiffusedNormalVoting::prob_measure(probfeatnode *probval, tensorType *tensor, int index)
{
	double radius, dDeltaRadius;

	if (_scale == 0.0 || _rmin == 0.0 || _rmax == 0.0 || _rmin >= _rmax || _inCloud->points.size() == 0)
	{
		cout << "invalid configuration parameters for classification module" << endl;
		return;
	}

	probfeatnode probfeatVal;

	dDeltaRadius = (_rmax - _rmin) / (_scale - 1.0);
	radius = _rmin;

	probfeatVal.prob[0] = 0;
	probfeatVal.prob[1] = 0;
	probfeatVal.prob[2] = 0;
	probfeatVal.featStrength[0] = 0;
	probfeatVal.featStrength[1] = 0;
	probfeatVal.featStrength[2] = 0;
	probfeatVal.csclcp[0] = 0;
	probfeatVal.csclcp[1] = 0;
	probfeatVal.csclcp[2] = 0;

	while (radius <= _rmax)
	{
		cout << "radius " << radius << " start (scale = " << _scale << endl;
		current_radius = radius;
		_sigma = radius;
		setAccumzero();
		voteAccum(radius, probfeatVal, index);
		
		_n_avg_acum.evec0[0] += _n_accum.evec0[0];
		_n_avg_acum.evec0[1] += _n_accum.evec0[0];
		_n_avg_acum.evec0[2] += _n_accum.evec0[0];
		
		_n_avg_acum.evec1[0] += _n_accum.evec1[0];
		_n_avg_acum.evec1[1] += _n_accum.evec1[1];
		_n_avg_acum.evec1[2] += _n_accum.evec1[2];

		_n_avg_acum.evec2[0] += _n_accum.evec2[0];
		_n_avg_acum.evec2[1] += _n_accum.evec2[1];
		_n_avg_acum.evec2[2] += _n_accum.evec2[2];

		probMeasure(probfeatVal, radius, index);
		radius += dDeltaRadius;
		cout << "radius " << radius - dDeltaRadius << " endl " << endl;
	}

	//	SATENDRA COMMENT
	_scale = 1.0; 
	//	Single tensor Implementation
	probfeatVal.csclcp[0] = 0;
	probfeatVal.csclcp[1] = 0;
	probfeatVal.csclcp[2] = 0;

	_accum.swap(averaged_accum);
	
	probMeasure(probfeatVal, radius, index); // SATENDRA COMMENT

	// ground_density(probfeatVal);
	// writeGlyphVars(radius - dDeltaRadius);  SATENDRA COMMENT

	probfeatVal.prob[0] = probfeatVal.prob[0] / _scale;
	probfeatVal.prob[1] = probfeatVal.prob[1] / _scale;
	probfeatVal.prob[2] = probfeatVal.prob[2] / _scale;

	probfeatVal.featStrength[0] = probfeatVal.featStrength[0] / _scale;
	probfeatVal.featStrength[1] = probfeatVal.featStrength[1] / _scale;
	probfeatVal.featStrength[2] = probfeatVal.featStrength[2] / _scale;

	probfeatVal.csclcp[0] = probfeatVal.csclcp[0] / _scale;
	probfeatVal.csclcp[1] = probfeatVal.csclcp[1] / _scale;
	probfeatVal.csclcp[2] = probfeatVal.csclcp[2] / _scale;

	probfeatVal.planarity = probfeatVal.planarity / _scale;
	probfeatVal.anisotropy = probfeatVal.anisotropy / _scale;
	probfeatVal.sphericity = probfeatVal.sphericity / _scale;
	probfeatVal.linearity = probfeatVal.linearity / _scale;
	probfeatVal.omnivariance = probfeatVal.omnivariance / _scale;
	probfeatVal.eigenentropy = probfeatVal.eigenentropy / _scale;

	(*probval) = probfeatVal;
	(*tensor) = _n_accum;
	
	return;
}

bool DiffusedNormalVoting::probMeasure(probfeatnode &probfeatVal, float radius, int index)
{
		eigendecomposition(_n_glyphData, index);
		computeSaliencyVals(_n_glyphData);
		glyphAnalysis(_n_glyphData);

		glyphVars glyphtemp = _n_glyphData;

		if (glyphtemp.evals[2] == 0.0 && glyphtemp.evals[1] == 0.0 && glyphtemp.evals[0] == 0.0)
		{
			probfeatVal.prob[0] = probfeatVal.prob[0] + 1;
		}
		else
		{
			if (glyphtemp.evals[2] >= _epsi * glyphtemp.evals[0]) //ev0>ev1>ev2
				probfeatVal.prob[0] = probfeatVal.prob[0] + 1;

			if (glyphtemp.evals[1] < _epsi * glyphtemp.evals[0]) //ev0>ev1>ev2
				probfeatVal.prob[1] = probfeatVal.prob[1] + 1;

			if (glyphtemp.evals[2] < _epsi * glyphtemp.evals[0])  //ev0>ev1>ev2
				probfeatVal.prob[2] = probfeatVal.prob[2] + 1;


			probfeatVal.featStrength[0] += ((glyphtemp.evals[2] * glyphtemp.evals[1]) / (glyphtemp.evals[0] * glyphtemp.evals[0]));
			probfeatVal.featStrength[1] += ((glyphtemp.evals[2] * (glyphtemp.evals[0] - glyphtemp.evals[1])) / (glyphtemp.evals[0] * glyphtemp.evals[1]));
			probfeatVal.featStrength[2] += glyphtemp.evals[2] / (glyphtemp.evals[0] + glyphtemp.evals[1] + glyphtemp.evals[2]);

		}

		probfeatVal.csclcp[0] += _n_glyphData.csclcp[0];
		probfeatVal.csclcp[1] += _n_glyphData.csclcp[1];
		probfeatVal.csclcp[2] += _n_glyphData.csclcp[2];

		probfeatVal.sum_eigen = _n_glyphData.evals[0] + _n_glyphData.evals[1] + _n_glyphData.evals[2];
		probfeatVal.evecs[0] = _n_glyphData.evecs[0];
		probfeatVal.evecs[1] = _n_glyphData.evecs[1];
		probfeatVal.evecs[2] = _n_glyphData.evecs[2];
		probfeatVal.evecs[3] = _n_glyphData.evecs[3];
		probfeatVal.evecs[4] = _n_glyphData.evecs[4];
		probfeatVal.evecs[5] = _n_glyphData.evecs[5];
		probfeatVal.evecs[6] = _n_glyphData.evecs[6];
		probfeatVal.evecs[7] = _n_glyphData.evecs[7];
		probfeatVal.evecs[8] = _n_glyphData.evecs[8];

		if (_n_glyphData.evals[0] != 0)
		{
			float len = glyphtemp.evals[2] + glyphtemp.evals[1] + glyphtemp.evals[0];
			float lamda0 = glyphtemp.evals[0] / len;
			float lamda1 = glyphtemp.evals[1] / len;
			float lamda2 = glyphtemp.evals[2] / len;

			probfeatVal.planarity += (_n_glyphData.evals[0] - _n_glyphData.evals[1]) / _n_glyphData.evals[0];
			probfeatVal.anisotropy += (_n_glyphData.evals[1] - _n_glyphData.evals[2]) / _n_glyphData.evals[0];
			probfeatVal.sphericity += (_n_glyphData.evals[0] - _n_glyphData.evals[2]) / _n_glyphData.evals[0];

			probfeatVal.linearity += (_n_glyphData.evals[0] - _n_glyphData.evals[1]) / _n_glyphData.evals[0];
			probfeatVal.omnivariance += cbrt(_n_glyphData.evals[1] * _n_glyphData.evals[2] * _n_glyphData.evals[0]);
			probfeatVal.eigenentropy += -1 * (_n_glyphData.evals[0] * log(_n_glyphData.evals[0]) + _n_glyphData.evals[1] * log(_n_glyphData.evals[1]) + _n_glyphData.evals[2] * log(_n_glyphData.evals[2]));

			probfeatVal.eigenvalues[0] = lamda0;
			probfeatVal.eigenvalues[1] = lamda1;
			probfeatVal.eigenvalues[2] = lamda2;
		}
		else
		{
			probfeatVal.planarity = 0;
			probfeatVal.anisotropy = 0;
			probfeatVal.sphericity = 0;
			probfeatVal.linearity = 0;
			probfeatVal.omnivariance = 0;
			probfeatVal.eigenentropy = 0;
		}

	cout << "in probMeasure() scale = " << _scale << endl;

	return true;
}

bool DiffusedNormalVoting::probMeasure(vector<probfeatnode>&probfeatVal, float radius)
{
	if (probfeatVal.size() != _inCloud->points.size())
		return false;

	_glyphData.clear();
	_glyphData.resize(_inCloud->points.size());

	for (size_t i = 0; i < _inCloud->points.size(); i++)
	{
		eigendecomposition(_glyphData[i], i);
		computeSaliencyVals(_glyphData[i]);
		glyphAnalysis(_glyphData[i]);

		glyphVars glyphtemp = _glyphData[i];

		if (glyphtemp.evals[2] == 0.0 && glyphtemp.evals[1] == 0.0 && glyphtemp.evals[0] == 0.0)
		{
			probfeatVal[i].prob[0] = probfeatVal[i].prob[0] + 1;
		}
		else
		{
			if (glyphtemp.evals[2] >= _epsi * glyphtemp.evals[0]) //ev0>ev1>ev2
				probfeatVal[i].prob[0] = probfeatVal[i].prob[0] + 1;

			if (glyphtemp.evals[1] < _epsi * glyphtemp.evals[0]) //ev0>ev1>ev2
				probfeatVal[i].prob[1] = probfeatVal[i].prob[1] + 1;

			if (glyphtemp.evals[2] < _epsi * glyphtemp.evals[0])  //ev0>ev1>ev2
				probfeatVal[i].prob[2] = probfeatVal[i].prob[2] + 1;


			probfeatVal[i].featStrength[0] += ((glyphtemp.evals[2] * glyphtemp.evals[1]) / (glyphtemp.evals[0] * glyphtemp.evals[0]));
			probfeatVal[i].featStrength[1] += ((glyphtemp.evals[2] * (glyphtemp.evals[0] - glyphtemp.evals[1])) / (glyphtemp.evals[0] * glyphtemp.evals[1]));
			probfeatVal[i].featStrength[2] += glyphtemp.evals[2] / (glyphtemp.evals[0] + glyphtemp.evals[1] + glyphtemp.evals[2]);

		}

		probfeatVal[i].csclcp[0] += _glyphData[i].csclcp[0];
		probfeatVal[i].csclcp[1] += _glyphData[i].csclcp[1];
		probfeatVal[i].csclcp[2] += _glyphData[i].csclcp[2];

		probfeatVal[i].sum_eigen = _glyphData[i].evals[0] + _glyphData[i].evals[1] + _glyphData[i].evals[2];
		probfeatVal[i].evecs[0] = _glyphData[i].evecs[0];
		probfeatVal[i].evecs[1] = _glyphData[i].evecs[1];
		probfeatVal[i].evecs[2] = _glyphData[i].evecs[2];
		probfeatVal[i].evecs[3] = _glyphData[i].evecs[3];
		probfeatVal[i].evecs[4] = _glyphData[i].evecs[4];
		probfeatVal[i].evecs[5] = _glyphData[i].evecs[5];
		probfeatVal[i].evecs[6] = _glyphData[i].evecs[6];
		probfeatVal[i].evecs[7] = _glyphData[i].evecs[7];
		probfeatVal[i].evecs[8] = _glyphData[i].evecs[8];

		//        if(radius == _rmin)
		//        {
		//            probfeatVal[i].evecs[3] = _glyphData[i].evecs[6];
		//            probfeatVal[i].evecs[4] = _glyphData[i].evecs[7];
		//            probfeatVal[i].evecs[5] = _glyphData[i].evecs[8];
		//        } else if(radius == _rmax) {
		//            probfeatVal[i].evecs[6] = _glyphData[i].evecs[6];
		//            probfeatVal[i].evecs[7] = _glyphData[i].evecs[7];
		//            probfeatVal[i].evecs[8] = _glyphData[i].evecs[8];
		//        }

		if (_glyphData[i].evals[0] != 0)
		{
			float len = glyphtemp.evals[2] + glyphtemp.evals[1] + glyphtemp.evals[0];
			float lamda0 = glyphtemp.evals[0] / len;
			float lamda1 = glyphtemp.evals[1] / len;
			float lamda2 = glyphtemp.evals[2] / len;

			probfeatVal[i].planarity += (_glyphData[i].evals[0] - _glyphData[i].evals[1]) / _glyphData[i].evals[0];
			probfeatVal[i].anisotropy += (_glyphData[i].evals[1] - _glyphData[i].evals[2]) / _glyphData[i].evals[0];
			probfeatVal[i].sphericity += (_glyphData[i].evals[0] - _glyphData[i].evals[2]) / _glyphData[i].evals[0];

			probfeatVal[i].linearity += (_glyphData[i].evals[0] - _glyphData[i].evals[1]) / _glyphData[i].evals[0];
			probfeatVal[i].omnivariance += cbrt(_glyphData[i].evals[1] * _glyphData[i].evals[2] * _glyphData[i].evals[0]);
			probfeatVal[i].eigenentropy += -1 * (_glyphData[i].evals[0] * log(_glyphData[i].evals[0]) + _glyphData[i].evals[1] * log(_glyphData[i].evals[1]) + _glyphData[i].evals[2] * log(_glyphData[i].evals[2]));


			probfeatVal[i].eigenvalues[0] = lamda0;
			probfeatVal[i].eigenvalues[1] = lamda1;
			probfeatVal[i].eigenvalues[2] = lamda2;
		}
		else
		{
			probfeatVal[i].planarity = 0;
			probfeatVal[i].anisotropy = 0;
			probfeatVal[i].sphericity = 0;
			probfeatVal[i].linearity = 0;
			probfeatVal[i].omnivariance = 0;
			probfeatVal[i].eigenentropy = 0;
		}
	}

	cout << "in probMeasure() scale = " << _scale << endl;

	// Uncomment to write eigenvalues per scale to file

//    char* filenames[6] = {"1.txt","2.txt","3.txt","4.txt","5.txt", "6.txt"};
//    ofstream fout(filenames[fNdx++]);
//    fout<< _glyphData.size()<<endl;
//    for(int i=0;i<probfeatVal.size(); i++) {

//        fout<< i<< " " << _glyphData[i].evals[0] << " "<< _glyphData[i].evals[1] << " "<< _glyphData[i].evals[2] << " " <<
//               _glyphData[i].evecs[0] << " " << _glyphData[i].evecs[1] << " " << _glyphData[i].evecs[2] << " " <<
//                _glyphData[i].evecs[3] << " " << _glyphData[i].evecs[4] << " " << _glyphData[i].evecs[5] << " " <<
//                _glyphData[i].evecs[6] << " " << _glyphData[i].evecs[7] << " " << _glyphData[i].evecs[8] <<endl;
//    }
//    fout.close();

	return true;

}

void DiffusedNormalVoting::writeGlyphVars(float radius)
{
	if (_glyphData.size() == 0)
		return;

	ofstream fout("glyphvars.txt");
	fout << _inCloud->points.size() << "\n";

	fout << radius << "\n";

	for (size_t i = 0; i < _glyphData.size(); i++)
	{
		fout << _inCloud->points[i].x << " " << _inCloud->points[i].y << " " << _inCloud->points[i].z << " ";
		fout << _glyphData[i].evals[0] << " " << _glyphData[i].evals[1] << " " << _glyphData[i].evals[2] << " ";
		fout << _glyphData[i].evecs[0] << " " << _glyphData[i].evecs[1] << " " << _glyphData[i].evecs[2] << " ";
		fout << _glyphData[i].evecs[3] << " " << _glyphData[i].evecs[4] << " " << _glyphData[i].evecs[5] << " ";
		fout << _glyphData[i].evecs[6] << " " << _glyphData[i].evecs[7] << " " << _glyphData[i].evecs[8] << " ";
		fout << _glyphData[i].uv[0] << " " << _glyphData[i].uv[1] << " ";
		fout << _glyphData[i].abc[0] << " " << _glyphData[i].abc[1] << " " << _glyphData[i].abc[2] << " ";
		fout << _glyphData[i].csclcp[0] << " " << _glyphData[i].csclcp[1] << " " << _glyphData[i].csclcp[2] << "\n";
	}

	fout.close();
	_glyphData.clear();

}

void DiffusedNormalVoting::computeSaliencyVals(glyphVars &glyphtemp)
{
	float len = glyphtemp.evals[2] + glyphtemp.evals[1] + glyphtemp.evals[0];

	float cl = 0.0, cp = 0.0, cs = 0.0;

	if (len != 0.0)
	{
		cl = (glyphtemp.evals[0] - glyphtemp.evals[1]) / len; //ev0>ev1>ev2
		cp = (2 * (glyphtemp.evals[1] - glyphtemp.evals[2])) / len;//(2.0 * (eigen_values[1] - eigen_values[0]));
		cs = 1 - (cl + cp); //1.0 - cl - cp;
	}

	glyphtemp.csclcp[0] = cs;
	glyphtemp.csclcp[1] = cl;
	glyphtemp.csclcp[2] = cp;
}

void DiffusedNormalVoting::glyphAnalysis(glyphVars &glyphtemp)
{

	double eps = 1e-4;
	double evals[3], uv[2], abc[3];

	double norm = sqrt(glyphtemp.evals[2] * glyphtemp.evals[2] + glyphtemp.evals[1] * glyphtemp.evals[1] + glyphtemp.evals[0] * glyphtemp.evals[0]);

	//    if(norm != 0)
	//    {
	//        glyphtemp.evals[2] = glyphtemp.evals[2]/norm;  //normalized the eigenvalues for superquadric glyph such that sqrt(lamda^2 + lambad^1 +lambda^0) = 1
	//        glyphtemp.evals[1] = glyphtemp.evals[1]/norm;
	//        glyphtemp.evals[0] = glyphtemp.evals[0]/norm;

	//    }

	evals[0] = glyphtemp.evals[2];   //ev0>ev1>ev2    //evals0>evals>1>evals2
	evals[1] = glyphtemp.evals[1];
	evals[2] = glyphtemp.evals[0];

	//tenGlyphBqdUvEval(uv, evals);
	//tenGlyphBqdAbcUv(abc, uv, 0.0);  // 3.0 for superquadric glyph, 0.0 for ellpsoid

	//norm=ELL_3V_LEN(evals);

	if (norm < eps)
	{
		double weight = norm / eps;
		abc[0] = weight * abc[0] + (1 - weight);
		abc[1] = weight * abc[1] + (1 - weight);
		abc[2] = weight * abc[2] + (1 - weight);
	}

	glyphtemp.uv[0] = uv[0];
	glyphtemp.uv[1] = uv[1];

	glyphtemp.abc[0] = abc[0];
	glyphtemp.abc[1] = abc[1];
	glyphtemp.abc[2] = abc[2];
}

void DiffusedNormalVoting::getdiffusionvelocity(Eigen::Vector3f evals, metaVelData *diffVel)
{
	evals[0] = -evals[0] / _delta;   //ev2>ev1>ev0
	evals[1] = -evals[1] / _delta;
	evals[2] = -evals[2] / _delta;

	evals[0] = exp(evals[0]);
	evals[1] = exp(evals[1]);
	evals[2] = exp(evals[2]);


	(*diffVel).vel[0] = evals[2];    //ev0>ev1>ev2 ->vel2>vel1>vel0
	(*diffVel).vel[1] = evals[1];
	(*diffVel).vel[2] = evals[0];

	(*diffVel).index[0] = 2;
	(*diffVel).index[1] = 1;
	(*diffVel).index[2] = 0;

	return;
}

bool DiffusedNormalVoting::eigendecomposition(glyphVars &glyphtemp, size_t idx)
{

	float A[3][3], V[3][3], d[3];
	Eigen::Vector3f eigen_values;
	metaVelData diffVel;

	for (int i = 0; i < 3; i++)
	{
		A[i][0] = _accum[idx].evec0[i];
		A[i][1] = _accum[idx].evec1[i];
		A[i][2] = _accum[idx].evec2[i];
	}

	eigen_decomposition(A, V, d); //d[2] > d[1] > d[0]

	eigen_values[0] = d[0];
	eigen_values[1] = d[1];
	eigen_values[2] = d[2];

	float len = d[0] + d[1] + d[2];
	// lamda0>lambda1>lambda2
	float lamda0 = d[2] / len;
	float lamda1 = d[1] / len;
	float lamda2 = d[0] / len;
	Eigen::MatrixXf e0(3, 1);
	e0 << V[0][2], V[1][2], V[2][2];
	e0.normalize();
	Eigen::MatrixXf e1(3, 1);
	e1 << V[0][1], V[1][1], V[2][1];
	e1.normalize();
	Eigen::MatrixXf e2(3, 1);
	e2 << V[0][0], V[1][0], V[2][0];
	e2.normalize();
	Eigen::Matrix3f T;
	T << 0, 0, 0, 0, 0, 0, 0, 0, 0;

	T += lamda0 * e0 * e0.transpose();
	T += lamda1 * e1 * e1.transpose();
	T += lamda2 * e2 * e2.transpose();

	averaged_tensors[idx].evec0[0] += T(0, 0);
	averaged_tensors[idx].evec0[1] += T(0, 1);
	averaged_tensors[idx].evec0[2] += T(0, 2);

	averaged_tensors[idx].evec1[0] += T(1, 0);
	averaged_tensors[idx].evec1[1] += T(1, 1);
	averaged_tensors[idx].evec1[2] += T(1, 2);

	averaged_tensors[idx].evec2[0] += T(2, 0);
	averaged_tensors[idx].evec2[1] += T(2, 1);
	averaged_tensors[idx].evec2[2] += T(2, 2);

	getdiffusionvelocity(eigen_values, &diffVel);

	glyphtemp.evals[2] = diffVel.vel[0];       //evals0>evals>1>evals2  //vel2>vel1>vel0
	glyphtemp.evals[1] = diffVel.vel[1];
	glyphtemp.evals[0] = diffVel.vel[2];

	glyphtemp.evecs[0] = V[0][diffVel.index[2]];
	glyphtemp.evecs[1] = V[1][diffVel.index[2]];
	glyphtemp.evecs[2] = V[2][diffVel.index[2]];

	glyphtemp.evecs[3] = V[0][diffVel.index[1]];
	glyphtemp.evecs[4] = V[1][diffVel.index[1]];
	glyphtemp.evecs[5] = V[2][diffVel.index[1]];

	glyphtemp.evecs[6] = V[0][diffVel.index[0]];
	glyphtemp.evecs[7] = V[1][diffVel.index[0]];
	glyphtemp.evecs[8] = V[2][diffVel.index[0]];

	return true;

}

//for each point, if it's local neighborhood has higher number of point-type features than
//line-type and surface-type, set all points to point-type features
void DiffusedNormalVoting::filterPointTypes(float radius, std::vector<probfeatnode>&probfeatVal, vector <tensorType> *aggregated_tensor)
{
	vector <tensorType> aggregated_tensor_tmp;
	aggregated_tensor_tmp.resize(_inCloud->points.size());
	for (size_t i = 0; i < _inCloud->points.size(); i++)
	{
		aggregated_tensor_tmp[i] = (*aggregated_tensor)[i];
	}

	std::vector<int> pointIdxRadiusSearch;
	std::vector<float> pointRadiusSquaredDistance;
	pcl::PointXYZ searchPoint;

	for (size_t i = 0; i < _inCloud->points.size(); i++)
	{
		searchPoint = _inCloud->points[i];

		if (search->radiusSearch(searchPoint, radius, pointIdxRadiusSearch, pointRadiusSquaredDistance) > 0)
		{
			int countPointType = 0;
			float sum_cl = 0.0;
			float sum_cp = 0.0;
			float sum_cs = 0.0;

			for (size_t j = 0; j < pointIdxRadiusSearch.size(); j++)
			{
				if (probfeatVal[pointIdxRadiusSearch[j]].csclcp[0] > probfeatVal[pointIdxRadiusSearch[j]].csclcp[2] && _glyphData[pointIdxRadiusSearch[j]].csclcp[0] > _glyphData[pointIdxRadiusSearch[j]].csclcp[1])
				{
					countPointType++;
					sum_cl += probfeatVal[pointIdxRadiusSearch[j]].csclcp[1];
					sum_cp += probfeatVal[pointIdxRadiusSearch[j]].csclcp[2];
					sum_cs += probfeatVal[pointIdxRadiusSearch[j]].csclcp[0];
				}
			}

			if (pointIdxRadiusSearch.size() - countPointType <= countPointType)
			{
				tensorType averageTensor;
				averageTensor.evec0[0] = 0;
				averageTensor.evec0[1] = 0;
				averageTensor.evec0[2] = 0;
				averageTensor.evec1[0] = 0;
				averageTensor.evec1[1] = 0;
				averageTensor.evec1[2] = 0;
				averageTensor.evec2[0] = 0;
				averageTensor.evec2[1] = 0;
				averageTensor.evec2[2] = 0;
				for (size_t j = 0; j < pointIdxRadiusSearch.size(); j++)
				{
					averageTensor.evec0[0] += aggregated_tensor_tmp[pointIdxRadiusSearch[j]].evec0[0];
					averageTensor.evec0[1] += aggregated_tensor_tmp[pointIdxRadiusSearch[j]].evec0[1];
					averageTensor.evec0[2] += aggregated_tensor_tmp[pointIdxRadiusSearch[j]].evec0[2];

					averageTensor.evec1[0] += aggregated_tensor_tmp[pointIdxRadiusSearch[j]].evec1[0];
					averageTensor.evec1[1] += aggregated_tensor_tmp[pointIdxRadiusSearch[j]].evec1[1];
					averageTensor.evec1[2] += aggregated_tensor_tmp[pointIdxRadiusSearch[j]].evec1[2];

					averageTensor.evec2[0] += aggregated_tensor_tmp[pointIdxRadiusSearch[j]].evec2[0];
					averageTensor.evec2[1] += aggregated_tensor_tmp[pointIdxRadiusSearch[j]].evec2[1];
					averageTensor.evec2[2] += aggregated_tensor_tmp[pointIdxRadiusSearch[j]].evec2[2];
				}
				averageTensor.evec0[0] /= countPointType;
				averageTensor.evec0[1] /= countPointType;
				averageTensor.evec0[2] /= countPointType;
				averageTensor.evec1[0] /= countPointType;
				averageTensor.evec1[1] /= countPointType;
				averageTensor.evec1[2] /= countPointType;
				averageTensor.evec2[0] /= countPointType;
				averageTensor.evec2[1] /= countPointType;
				averageTensor.evec2[2] /= countPointType;

				for (size_t j = 0; j < pointIdxRadiusSearch.size(); j++)
				{
					probfeatVal[pointIdxRadiusSearch[j]].csclcp[0] = sum_cs / countPointType;
					probfeatVal[pointIdxRadiusSearch[j]].csclcp[1] = sum_cl / countPointType;
					probfeatVal[pointIdxRadiusSearch[j]].csclcp[2] = sum_cp / countPointType;

					(*aggregated_tensor)[pointIdxRadiusSearch[j]] = averageTensor;
				}
			}
			pointIdxRadiusSearch.clear();
			pointRadiusSquaredDistance.clear();
		}
	}
}