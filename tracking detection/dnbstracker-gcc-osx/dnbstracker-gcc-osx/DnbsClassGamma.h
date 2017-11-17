#pragma once
#include "dnbsclassbeta.h"

class DnbsClassGamma :
	public DnbsClassBeta
{
public:
	DnbsClassGamma(void);
	~DnbsClassGamma(void);

public:
	double inprod_limit;
	vector<int> indexOfCentBases;
	vector< vector<int> > memberOfGroup;

public:
	// Main interface of the D-NBS class
	// Third Extensions

	void computeDNBS(
		int MaxNumBases,
		int MaxNumForegrounds,
		int MaxNumBackgrounds,
		IplImage * frame,
		CvRect rect,
		int margin_threshold_x,
		int margin_threshold_y,
		int distance_threshold,
		double lambda,
		int method = 0,
		int samp_width = 2,
		int samp_height = 2,
		double prune_ratio = 0.5,
		CvMat ** previous_foregrounds = NULL,
		double coherence = -1);

	bool genHaarBases_inProd(double minInProd);

	void computeDNBS_v1_branchbound(
		int MaxNumBases,
		double lambda,
		vector<CvMat*> foreground,
		vector<CvMat*> background,
		double inprod_limit);

	void computeDNBS_v1_cluster(
		int MaxNumBases,
		double lambda,
		vector<CvMat*> foreground,
		vector<CvMat*> background,
		double ratio);

	bool genHaarBases_innProdFast(double minInnProd);

	void computeDNBS_v1_cluster_v2(
		int MaxNumBases,
		double lambda,
		vector<CvMat*> foreground,
		vector<CvMat*> background,
		double ratio);

	void computeDNBS_v2_cluster(
		int MaxNumBases,
		double lambda,
		vector<CvMat*> foreground,
		vector<CvMat*> background,
		double ratio);

	void computeDNBS_v1_branchbound_latest(
		int MaxNumBases,
		double lambda,
		vector<CvMat*> foreground,
		vector<CvMat*> background);

	vector<int> computeDNBS_v1_customFeatSet(
		int MaxNumBases,
		double lambda,
		const vector<CvMat*> & foreground,
		const vector<CvMat*> & background,
		const vector<int> & featIndex);

	vector<int> genDownScaleFeatures(int sampwidth, int sampheight);
	vector<int> genDiffuseFeatures(const vector<int>& index, double inprod);


	void computeNBS_v1_cluster_v2(int MaxNumBases, double ratio);

	void sample_backgrounds_supress(CvRect rect,
		int MaxNumBackgrounds,
		int distance_threshold,
		int x_min, int x_max,
		int y_min, int y_max,
		vector< CvRect > & pnt,
		vector< pair<double, int> > & arr,
		double ** weight,
		vector < CvRect > & samples,
		double nonoverlap_ratio = 0.7);
};
