#include "Processing.h"
#include "FileProcess.h"
#include <iostream>
#include <fstream>
#include "ClassifiersFactory.h"

#define LABEL_FILE "data/label/d4-alg2l.txt"
#define OPTIMAL_FILE "data/label/optimal.txt";
#define CLOUD_FILE "data/3d_4.txt";

bool comparede(double i, double j) { return i < j; }

Processing::Processing() : _inCloud(new cl)
{
}

Processing::~Processing()
{
}

bool Processing::fileRead(string cloudFile, string labelFile, string optimalFile)
{
	bool nErrorCode;
	_filePath = cloudFile;

	if (_filePath == "NULL")
	{
		std::cout << "Invalid File Path" << std::endl;
		return false;
	}
	else
	{
		FileProcess *readObject = new FileProcess;
		nErrorCode = readObject->read(_filePath, _inCloud, _intensity);
		_FileExtension = readObject->getFileType(_filePath);

		string labelFile = labelFile;
		readObject->readLabels(labelFile, labels);

		string optimalFile = optimalFile;
		readObject->readLabels(optimalFile, optimalScale);

		delete readObject;
	}

	if (_FileExtension == "las" || _FileExtension == "ply" || _FileExtension == "off" || _FileExtension == "xyz")
		//scale();

	//buildOctree();

	cout << "total size " << _inCloud->points.size() << endl;

	return nErrorCode;
}

void Processing::scale()
{
	if (_inCloud->points.size() == 0)
		return;

	float minX, minY, minZ, maxX, maxY, maxZ; //,  minIntensity, maxIntensiy,
	float DelX, DelY, DelZ;

	float minxr, maxxr, minyr, maxyr, minzr, maxzr, delxr, delyr, delzr;
	float scale_fac = 2.0;

	vector<double> tempData;

	tempData.resize(_inCloud->points.size());

	for (size_t i = 0; i < _inCloud->points.size(); i++)
		tempData[i] = _inCloud->points[i].x;

	minX = *std::min_element(tempData.begin(), tempData.end(), comparede);
	maxX = *std::max_element(tempData.begin(), tempData.end(), comparede);

	for (size_t i = 0; i < _inCloud->points.size(); i++)
		tempData[i] = _inCloud->points[i].y;

	minY = *std::min_element(tempData.begin(), tempData.end(), comparede);
	maxY = *std::max_element(tempData.begin(), tempData.end(), comparede);

	for (size_t i = 0; i < _inCloud->points.size(); i++)
		tempData[i] = _inCloud->points[i].z;

	minZ = *std::min_element(tempData.begin(), tempData.end(), comparede);
	maxZ = *std::max_element(tempData.begin(), tempData.end(), comparede);

	tempData.clear();

	DelX = maxX - minX;
	DelY = maxY - minY;
	DelZ = maxZ - minZ;

	if (DelX >= DelY && DelX >= DelZ)
	{
		minxr = 0.0;
		maxxr = 1.0 * scale_fac;

		minyr = 0;
		maxyr = DelY / DelX * scale_fac;

		minzr = 0.0;
		maxzr = (DelZ / DelX) * scale_fac;
	}
	else if (DelY >= DelX && DelY >= DelZ)
	{
		minxr = 0;
		maxxr = DelX / DelY * scale_fac;

		minyr = 0.0;
		maxyr = 1.0 * scale_fac;

		minzr = 0.0;
		maxzr = (DelZ / DelY) * scale_fac;
	}
	else if (DelZ >= DelX && DelZ >= DelY)
	{
		//cout<<"notify it "<<endl;
		minxr = 0.0;
		maxxr = DelX / DelZ * scale_fac;

		minyr = 0.0;
		maxyr = DelY / DelZ * scale_fac;

		minzr = 0.0;
		maxzr = 1.0 * scale_fac;
	}

	delxr = maxxr - minxr;
	delyr = maxyr - minyr;
	delzr = maxzr - minzr;

	cout << "DelX " << DelX << endl;
	cout << "DelY " << DelY << endl;
	cout << "DelZ " << DelZ << endl;

	for (size_t i = 0; i < _inCloud->points.size(); i++)
	{
		_inCloud->points[i].x = delxr * (_inCloud->points[i].x - minX) / DelX;
		_inCloud->points[i].y = delyr * (_inCloud->points[i].y - minY) / DelY;
		_inCloud->points[i].z = delzr * (_inCloud->points[i].z - minZ) / DelZ;
	}
}

void Processing::buildStructure(ProcessRequest request)
{
	// get all the file paths and pass it to the fileRead
	fileRead(request.cloudFileName, request.labelFileName, request.optimalFileName);

	Configuration* config = new Configuration();
	config->SetValue("Resolution", "125");
	SearchStructure structure = static_cast<SearchStructure>(request.dataStructure); // KD-Tree or Oct-Tree
	search = SearchNeighbourFactory::GetNeighbourSearchDataStructure(structure);
	search->build(_inCloud, config);
}

void Processing::processCloud(ProcessRequest request)
{
	ClassifierTypes classifierType = static_cast<ClassifierTypes>(request.classifierType); // Descriptor type
	Classifier* classifier = ClassifiersFactory::GetClassifier(classifierType);
	classifier->setSearch(search);
	classifier->setInputCloud(_inCloud);

	vector <probfeatnode> _probval;
	vector <tensorType> _accum; // aggregated tensor

	//for loop each will be passed with the configuration required...
	classifier->process(request, &_probval, &_accum);
}

pcl::PointCloud<pcl::PointXYZ>::Ptr Processing::getCloud()
{
	return _inCloud;
}

// This is the request to process against each point
// It can take the selection here for the discriptor type for each point
// But the Search Structure whether Kdtree or Octree will be fixed though.
void Processing::processPoint(ProcessRequest request, int index)
{
	ClassifierTypes classifierType = static_cast<ClassifierTypes>(request.classifierType); // Descriptor type
	Classifier* classifier = ClassifiersFactory::GetClassifier(classifierType);
	classifier->setSearch(search);
	classifier->setInputCloud(_inCloud);

	probfeatnode _probval;
	tensorType _accum; // aggregated tensor

	classifier->processPoint(request, &_probval, &_accum, index);
}