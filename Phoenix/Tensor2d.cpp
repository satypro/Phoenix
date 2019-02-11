#include "Tensor2d.h"
#include <pcl/common/common.h>
#include <pcl/common/eigen.h>
#include <pcl/common/centroid.h>
#include <pcl/point_types.h>
#include <pcl/features/normal_3d.h>
#include <math.h>
#include "eig3.h"
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */


std::vector<int> neighbours;
int center_point = 925;

/***************************************************************************************/
void Tensor2d::setInputCloud(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud)
{
	_inCloud = cloud;
}
/*****************************************************************************************

*****************************************************************************************/
void Tensor2d::setOctree(OctreeXYZ<PointT, LeafContainerT, BranchContainerT > * octree)
{
	_octree = octree;
}
/*****************************************************************************************

*****************************************************************************************/

void Tensor2d::setParams(float rmin, float rmax, float rmaxpt, float epsi, float scale)
{
	_rmin = rmin;
	_rmax = rmax;
	_rmaxpt = rmaxpt;
	_epsi = epsi;
	_scale = scale;

	_sigma = 1.0;
}

/*****************************************************************************************
*****************************************************************************************/
bool Tensor2d::setAccumzero()
{
	if (_inCloud->points.size() == 0)
		return false;

	_accum.resize(_inCloud->points.size());
	averaged_tensors.resize(_inCloud->points.size());

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
		}
	}
	return true;
}
/*****************************************************************************************
*****************************************************************************************/
void Tensor2d::getnormalizeTensor(size_t idx, float weight)
{
	if (weight == 0.0)
		return;

	for (int j = 0; j < 3; j++)
	{
		_accum[idx].evec0[j] = _accum[idx].evec0[j] / weight;
		_accum[idx].evec1[j] = _accum[idx].evec1[j] / weight;
		_accum[idx].evec2[j] = _accum[idx].evec2[j] / weight;
	}

}
/*****************************************************************************************/
void Tensor2d::calculatePartialDerivative(float radius)
{
	std::vector<int> pointIdxRadiusSearch;
	std::vector<float> pointRadiusSquaredDistance;
	pcl::PointXYZ searchPoint;

	for (size_t i = 0; i < _inCloud->points.size(); i++)
	{
		searchPoint = _inCloud->points[i];
		float fx0 = 0.0;
		float fy0 = 0.0;

		//        for(size_t p =0; p < _inCloud->points.size(); p++)
		//        {
		//            if(p!=i)
		//            {
		//                pcl::PointXYZ tempPoint = _inCloud->points[p];
		//                if((tempPoint.x-searchPoint.x)*(tempPoint.x-searchPoint.x)
		//                        + (tempPoint.y-searchPoint.y)*(tempPoint.y-searchPoint.y)
		//                        < radius*radius)
		//                {
		//                    pointIdxRadiusSearch.push_back(p);
		//                }
		//            }
		//        }

		//        if(pointIdxRadiusSearch.size()>0)
		if (_octree->radiusSearch(searchPoint, radius, pointIdxRadiusSearch, pointRadiusSquaredDistance) > 0)
		{
			int N = pointIdxRadiusSearch.size();
			float rxx = 0.0;
			float rxy = 0.0;
			float ryy = 0.0;
			std::vector<float> Wx(N, 0);
			std::vector<float> Wy(N, 0);


			for (size_t h = 0; h < pointIdxRadiusSearch.size(); h++)
			{
				rxx += double(_inCloud->points[pointIdxRadiusSearch[h]].x - searchPoint.x) * double(_inCloud->points[pointIdxRadiusSearch[h]].x - searchPoint.x);
				rxy += double(searchPoint.x - _inCloud->points[pointIdxRadiusSearch[h]].x) * double(searchPoint.y - _inCloud->points[pointIdxRadiusSearch[h]].y);
				ryy += double(_inCloud->points[pointIdxRadiusSearch[h]].y - searchPoint.y) * double(_inCloud->points[pointIdxRadiusSearch[h]].y - searchPoint.y);
			}
			float D = rxx * ryy - rxy * rxy;
			//            std::cout << "rxx = " << rxx << " ryy = " << ryy << " rxy = " << rxy << std::endl;

			for (size_t h = 0; h < pointIdxRadiusSearch.size(); h++)
			{
				Wx[h] = double(_inCloud->points[pointIdxRadiusSearch[h]].x - searchPoint.x) * ryy - double(_inCloud->points[pointIdxRadiusSearch[h]].y - searchPoint.y) * rxy;
				Wy[h] = double(_inCloud->points[pointIdxRadiusSearch[h]].y - searchPoint.y) * rxx - double(_inCloud->points[pointIdxRadiusSearch[h]].x - searchPoint.x) * rxy;
				if (D > 0.00001) {
					Wx[h] /= D;
					Wy[h] /= D;
				}
				//                Wx[h] /= D;
				//                Wy[h] /= D;
				//                std::cout << "Wx[h] = " << Wx[h] << " Wx[h]/D = " << Wx[h]/D << " Wy[h] = " << Wy[h] << " Wy[h]/D = " << Wy[h]/D << std::endl;

				fx0 += Wx[h] * double(_inCloud->points[pointIdxRadiusSearch[h]].z - searchPoint.z);
				fy0 += Wy[h] * double(_inCloud->points[pointIdxRadiusSearch[h]].z - searchPoint.z);
			}

			pointIdxRadiusSearch.clear();
			pointRadiusSquaredDistance.clear();
		}
		Ix.push_back(fx0);
		Iy.push_back(fy0);
	}
}

/*****************************************************************************************/
bool Tensor2d::structureTensor2D(float radius)
{
	calculatePartialDerivative(radius);
	std::vector<int> pointIdxRadiusSearch;
	std::vector<float> pointRadiusSquaredDistance;
	pcl::PointXYZ searchPoint;
	float  weight;


	for (size_t i = 0; i < _inCloud->points.size(); i++)
	{
		searchPoint = _inCloud->points[i];

		//        for(size_t p =0; p < _inCloud->points.size(); p++)
		//        {
		//            if(p!=i)
		//            {
		//                pcl::PointXYZ tempPoint = _inCloud->points[p];
		//                if((tempPoint.x-searchPoint.x)*(tempPoint.x-searchPoint.x)
		//                        + (tempPoint.y-searchPoint.y)*(tempPoint.y-searchPoint.y)
		//                        < radius*radius)
		//                {
		//                    pointIdxRadiusSearch.push_back(p);
		//                }
		//            }
		//        }

		//        if(pointIdxRadiusSearch.size()>0)
		if (_octree->radiusSearch(searchPoint, radius, pointIdxRadiusSearch, pointRadiusSquaredDistance) > 0)
		{
			weight = 0.0;

			for (size_t h = 0; h < pointIdxRadiusSearch.size(); h++)
			{
				ColumnVector Vect1;
				Vect1(0, 0) = double(_inCloud->points[pointIdxRadiusSearch[h]].x - searchPoint.x);
				Vect1(1, 0) = double(_inCloud->points[pointIdxRadiusSearch[h]].y - searchPoint.y);
				Vect1(2, 0) = double(_inCloud->points[pointIdxRadiusSearch[h]].z - searchPoint.z);

				double  s = Vect1.norm(); // if s is zero

//                double coeff = radius - s;

				double coeff = 1 / (2 * 3.141592653*radius*radius)*exp(-1 * s*s / (2 * radius*radius));

				weight += coeff;


				_accum[i].evec0[0] += coeff * Ix[pointIdxRadiusSearch[h]] * Ix[pointIdxRadiusSearch[h]];
				_accum[i].evec0[1] += coeff * Ix[pointIdxRadiusSearch[h]] * Iy[pointIdxRadiusSearch[h]];
				_accum[i].evec0[2] = 0;

				_accum[i].evec1[0] += coeff * Ix[pointIdxRadiusSearch[h]] * Iy[pointIdxRadiusSearch[h]];
				_accum[i].evec1[1] += coeff * Iy[pointIdxRadiusSearch[h]] * Iy[pointIdxRadiusSearch[h]];
				_accum[i].evec1[2] = 0;

				_accum[i].evec2[0] = 0;
				_accum[i].evec2[1] = 0;
				_accum[i].evec2[2] = 0;
			}
			getnormalizeTensor(i, weight);
			pointIdxRadiusSearch.clear();
			pointRadiusSquaredDistance.clear();
		}
	}
	return true;
}
/*****************************************************************************************
Function Name		:	PointClassificaton::probMeasure
Calls			:	PointClassificaton::CurvatureEstimation
Output Params		:	T *curve, T *disc, T *spherical
Return			:	int
*****************************************************************************************/
void Tensor2d::prob_measure(vector <probfeatnode> *probval)
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

	//radius = 0.5*(_rmin+_rmax);
//    _scale = 1.0;
	//radius = 0.009;

	while (radius <= _rmax)
	{
		cout << "radius (new)  " << radius << " start" << endl;
		_sigma = radius;
		setAccumzero();
		structureTensor2D(radius);
		probMeasure(probfeatVal, radius);
		radius += dDeltaRadius;
		cout << "raidus " << radius - dDeltaRadius << " done" << endl;
	}

	// Single tensor
//     for(size_t i = 0; i < _inCloud->points.size(); i++)
//     {
//         probfeatVal[i].csclcp[0] = 0;
//         probfeatVal[i].csclcp[1] = 0;
//         probfeatVal[i].csclcp[2] = 0;
//         for(int j=0;j<3;j++)
//         {
//             averaged_tensors[i].evec0[j] /= _scale;
//             averaged_tensors[i].evec1[j] /= _scale;
//             averaged_tensors[i].evec2[j] /= _scale;
//         }
//     }
//     _scale = 1.0;
	_accum.swap(averaged_tensors);
	//     probMeasure(probfeatVal,radius);

	writeGlyphVars(radius - dDeltaRadius);

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
		probfeatVal[i].sphericity /= _scale;
		probfeatVal[i].anisotropy /= _scale;
		probfeatVal[i].sum_eigen /= _scale;

		//        if(probfeatVal[i].csclcp[2]>probfeatVal[i].csclcp[0] && probfeatVal[i].csclcp[2]>probfeatVal[i].csclcp[1]){
		//            probfeatVal[i].csclcp[1] = 1;
		//            probfeatVal[i].csclcp[0] = 0;
		//            probfeatVal[i].csclcp[2] = 0;
		//        }else{
		//            probfeatVal[i].csclcp[1] = 0;
		//            probfeatVal[i].csclcp[0] = 0;
		//            probfeatVal[i].csclcp[2] = 0;
		//        }

	}



	probval->resize(_inCloud->points.size());
	for (size_t i = 0; i < _inCloud->points.size(); i++)
		(*probval)[i] = probfeatVal[i];

	probfeatVal.clear();
	return;
}


/*****************************************************************************************
Function Name		:	PointClassificaton::probMeasure
Calls			:	PointClassificaton::CurvatureEstimation
Output Params		:	T *curve, T *disc, T *spherical
Return			:	int
Vote analysis according to modrhobahi paper and computing probability according to kreylo's paper.
*****************************************************************************************/
bool Tensor2d::probMeasure(vector<probfeatnode>&probfeatVal, float radius)
{
	if (probfeatVal.size() != _inCloud->points.size())
		return false;

	_glyphData.clear();
	_glyphData.resize(_inCloud->points.size());

	average_eigenvalues();
	for (size_t i = 0; i < _inCloud->points.size(); i++)
	{
		eigendecomposition(_glyphData[i], i);
		computeSaliencyVals(_glyphData[i], i, radius);
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

		//        probfeatVal[i].csclcp[0] = _glyphData[i].csclcp[0];
		//        probfeatVal[i].csclcp[1] = _glyphData[i].csclcp[1];
		//        probfeatVal[i].csclcp[2] = _glyphData[i].csclcp[2];
		probfeatVal[i].sum_eigen = _glyphData[i].evals[0] + _glyphData[i].evals[1] + _glyphData[i].evals[2];
		if (_glyphData[i].evals[0] != 0)
		{
			float len = glyphtemp.evals[1] + glyphtemp.evals[0];
			float lamda0 = glyphtemp.evals[0] / len;
			float lamda1 = glyphtemp.evals[1] / len;

			probfeatVal[i].planarity = 0;
			probfeatVal[i].anisotropy += (lamda0 - lamda1) / lamda0;
			probfeatVal[i].sphericity += lamda1 / lamda0;
		}
		else
		{
			probfeatVal[i].planarity += 0;
			probfeatVal[i].anisotropy += 0;
			probfeatVal[i].sphericity += 0;
		}
	}

	return true;
}

/*****************************************************************************************
*****************************************************************************************/
void Tensor2d::writeGlyphVars(float radius)
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

/*****************************************************************************************
*****************************************************************************************/
void Tensor2d::computeSaliencyVals(glyphVars &glyphtemp, int idx, float radius)
{
	float len = glyphtemp.evals[2] + glyphtemp.evals[1] + glyphtemp.evals[0];

	float cl = 0.0, cp = 0.0, cs = 0.0;

	if (len != 0.0)
	{
		// lamda0>lambda1>lambda2
		float lamda0 = glyphtemp.evals[0] / len;
		float lamda1 = glyphtemp.evals[1] / len;
		float lamda2 = glyphtemp.evals[2] / len;
		Eigen::MatrixXf e0(3, 1);
		e0 << glyphtemp.evecs[0], glyphtemp.evecs[1], glyphtemp.evecs[2];
		e0.normalize();
		Eigen::MatrixXf e1(3, 1);
		e1 << glyphtemp.evecs[3], glyphtemp.evecs[4], glyphtemp.evecs[5];
		e1.normalize();
		Eigen::MatrixXf e2(3, 1);
		e2 << glyphtemp.evecs[6], glyphtemp.evecs[7], glyphtemp.evecs[8];
		e2.normalize();
		Eigen::Matrix3f T;
		T << 0, 0, 0, 0, 0, 0, 0, 0, 0;

		//        if(lamda0>0)
		T += exp(-lamda0 / 0.16) * e0 * e0.transpose();
		//        if(lamda1>0)
		T += exp(-lamda1 / 0.16) * e1 * e1.transpose();
		//        if(lamda2>0)
		T += exp(-lamda2 / 0.16) * e2 * e2.transpose();

		averaged_tensors[idx].evec0[0] += T(0, 0);
		averaged_tensors[idx].evec0[1] += T(0, 1);
		averaged_tensors[idx].evec0[2] += T(0, 2);

		averaged_tensors[idx].evec0[3] += T(1, 0);
		averaged_tensors[idx].evec0[4] += T(1, 1);
		averaged_tensors[idx].evec0[5] += T(1, 2);

		averaged_tensors[idx].evec0[6] += T(2, 0);
		averaged_tensors[idx].evec0[7] += T(2, 1);
		averaged_tensors[idx].evec0[8] += T(2, 2);

		if (lamda0 > 0) {
			cp = (lamda0 - lamda1) / lamda0;
			cl = 1 - cp;
		}

		//        std::cout << lamda0 << " " << lamda1 << " " << lamda2 << endl;

		//        cout << T << endl << endl;
	}

	//    float T[2][2] = {{_accum[idx].evec0[0],_accum[idx].evec0[1]},
	//                     {_accum[idx].evec1[0],_accum[idx].evec1[1]}};

	//    float mu1 = 0.5 * (T[0][0] + T[1][1] + sqrt((T[0][0] - T[1][1])*(T[0][0] - T[1][1]) + 4 * T[0][1]*T[0][1]));
	//    float mu2 = 0.5 * (T[0][0] + T[1][1] - sqrt((T[0][0] - T[1][1])*(T[0][0] - T[1][1]) + 4 * T[0][1]*T[0][1]));

	//    float phi = 0.5 * atan2 (2*T[0][1], (T[0][0] - T[1][1]));
	//    float normal[2] = {cos(phi),sin(phi)};

	//    mu1 = -1/log(mu1);
	//    mu2 = -1/log(mu2);

	////    mu1 /= average_mu1;
	////    mu2 /= average_mu2;

	//    std::cout << mu1 << " ," << mu2 << " radius = " << radius << std::endl;

	////    if(mu1<0.06 || mu2<0.06){
	////            mu1 = 0; mu2 = 0;
	////    }

	//    if(mu1> 0.0){
	//        cl = (mu1-mu2)/(mu1+mu2);
	//        cp = 1 - cl;
	//    }

	//    glyphtemp.evals[0] = mu1;
	//    glyphtemp.evals[1] = mu2;

	////    std::cout << mu1 << " ," << mu2 << " radius = " << radius << std::endl;
	////    std::cout << " cl = " << cl << std::endl;

	////    float len = glyphtemp.evals[2] + glyphtemp.evals[1] + glyphtemp.evals[0];
	////    if(len!= 0.0)
	////    {
	////        cl = (glyphtemp.evals[0] - glyphtemp.evals[1])/len; //ev0>ev1>ev2
	////        cp = (2*(glyphtemp.evals[1] - glyphtemp.evals[2]))/len ;//(2.0 * (eigen_values[1] - eigen_values[0]));
	////        cs = 1 - (cl+cp); //1.0 - cl - cp;
	////    }

	////    if(mu1+mu2 != 0)
	////    {
	////        cl = (mu1-mu2)/(mu1+mu2);
	////        cp = 1 - cl;
	////        cs = 0;
	////    }
	////    else
	////      if(cl+cp < 0.0001){
	////        cs = 1.0;
	////        cp = 0.0;
	////        cl = 0.0;
	////      }
	//    // Calculate cp from local maximas
	//    std::vector<int> pointIdxRadiusSearch;
	//    std::vector<float> pointRadiusSquaredDistance;
	//    pcl::PointXYZ searchPoint;

	//    searchPoint = _inCloud->points[idx];

	//    bool is_junction = true;
	//    if(mu1+mu2 == 0)
	//        is_junction = false;
	//    else{
	//        if(_octree->radiusSearch(searchPoint, radius, pointIdxRadiusSearch, pointRadiusSquaredDistance) > 0)
	//        {
	//            for(size_t h =0; h < pointIdxRadiusSearch.size(); h++)
	//            {
	//                float T_neighbour[2][2];

	//                for(int i = 0; i < 2; i++)
	//                {
	//                    T_neighbour[i][0] = _accum[pointIdxRadiusSearch[h]].evec0[i];
	//                    T_neighbour[i][1] = _accum[pointIdxRadiusSearch[h]].evec1[i];
	//                }

	//                float mu2_neighbour = 0.5 * (T_neighbour[0][0] + T_neighbour[1][1] - sqrt((T_neighbour[0][0] - T_neighbour[1][1])*(T_neighbour[0][0] - T_neighbour[1][1]) + 4 * T_neighbour[0][1]*T_neighbour[0][1]));
	//    //            std::cout << mu2_neighbour << " " << mu2  << std::endl;
	//                if(mu2_neighbour > mu2)
	//                {
	//                    is_junction = false;
	//                    break;
	//                }
	//            }
	//            pointIdxRadiusSearch.clear();
	//            pointRadiusSquaredDistance.clear();
	//        }
	//    }

	//    if(is_junction)
	//    {
	//        cp = 1.0;
	//    }
	//    else
	//    {
	//        if(mu1+mu2 != 0)
	//        {
	////            float mu1_new = mu1;//sqrt(mu1-mu2);
	////            cl = (mu1_new-mu2)/(mu1_new+mu2);
	////            cp = 1 - cl;
	////            cs = 0.0;
	//        }
	//        else
	//            cs = 0.0;
	//    }

	glyphtemp.csclcp[1] = cl;
	glyphtemp.csclcp[0] = cs;
	glyphtemp.csclcp[2] = cp;

}

/*****************************************************************************************
*****************************************************************************************/

void Tensor2d::average_eigenvalues()
{
	double sum_mu1 = 0.0;
	double sum_mu2 = 0.0;
	double mxmu1 = -1000;
	double mxmu2 = -1000;
	double mnmu1 = 1000;
	double mnmu2 = 1000;

	for (int idx = 0; idx < _accum.size(); idx++)
	{
		float T[2][2] = { {_accum[idx].evec0[0],_accum[idx].evec0[1]},
						 {_accum[idx].evec1[0],_accum[idx].evec1[1]} };

		double mu1 = 0.5 * (T[0][0] + T[1][1] + sqrt((T[0][0] - T[1][1])*(T[0][0] - T[1][1]) + 4 * T[0][1] * T[0][1]));
		double mu2 = 0.5 * (T[0][0] + T[1][1] - sqrt((T[0][0] - T[1][1])*(T[0][0] - T[1][1]) + 4 * T[0][1] * T[0][1]));

		sum_mu1 += mu1;
		sum_mu2 += mu2;

		if (mu1 > mxmu1)
			mxmu1 = mu1;
		if (mu2 > mxmu2)
			mxmu2 = mu2;
		if (mu1 < mnmu1)
			mnmu1 = mu1;
		if (mu2 < mnmu2)
			mnmu2 = mu2;
	}
	average_mu1 = sum_mu1 / _accum.size();
	average_mu2 = sum_mu2 / _accum.size();
	max_mu1 = mxmu1;
	max_mu2 = mxmu2;
	min_mu1 = mnmu1;
	min_mu2 = mnmu2;
}


/*****************************************************************************************
*****************************************************************************************/

void Tensor2d::glyphAnalysis(glyphVars &glyphtemp)
{

	double eps = 1e-4;
	double evals[3], uv[2], abc[3];

	double norm = sqrt(glyphtemp.evals[2] * glyphtemp.evals[2] + glyphtemp.evals[1] * glyphtemp.evals[1] + glyphtemp.evals[0] * glyphtemp.evals[0]);

	if (norm != 0)
	{
		glyphtemp.evals[2] = glyphtemp.evals[2] / norm;  //normalized the eigenvalues for superquadric glyph such that sqrt(lamda^2 + lambad^1 +lambda^0) = 1
		glyphtemp.evals[1] = glyphtemp.evals[1] / norm;
		glyphtemp.evals[0] = glyphtemp.evals[0] / norm;

	}

	evals[0] = glyphtemp.evals[2];   //ev0>ev1>ev2
	evals[1] = glyphtemp.evals[1];
	evals[2] = glyphtemp.evals[0];

	//tenGlyphBqdUvEval(uv, evals);
	//tenGlyphBqdAbcUv(abc, uv, 3.0);

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

/*****************************************************************************************
*****************************************************************************************/
bool Tensor2d::eigendecomposition(glyphVars &glyphtemp, size_t idx)
{

	float A[3][3], V[3][3], d[3];

	for (int i = 0; i < 3; i++)
	{
		A[i][0] = _accum[idx].evec0[i];
		A[i][1] = _accum[idx].evec1[i];
		A[i][2] = _accum[idx].evec2[i];
	}

	eigen_decomposition(A, V, d);

	//    if(d[0]!=0)
	//        cout << d[2] << " " << d[1] << " " << d[0] << endl;
	//    if(fabs(d[0])<1.0e-15)
	//        d[0]=0;
	//    if(fabs(d[1])<1.0e-15)
	//        d[1]=0;
	//    if(fabs(d[2])<1.0e-15)
	//        d[2]=0;

	glyphtemp.evals[2] = d[0];   //d[2] > d[1] > d[0]
	glyphtemp.evals[1] = d[1];
	glyphtemp.evals[0] = d[2];

	glyphtemp.evecs[0] = V[0][2];
	glyphtemp.evecs[1] = V[1][2];
	glyphtemp.evecs[2] = V[2][2];

	glyphtemp.evecs[3] = V[0][1];
	glyphtemp.evecs[4] = V[1][1];
	glyphtemp.evecs[5] = V[2][1];

	glyphtemp.evecs[6] = V[0][0];
	glyphtemp.evecs[7] = V[1][0];
	glyphtemp.evecs[8] = V[2][0];

	return true;
}