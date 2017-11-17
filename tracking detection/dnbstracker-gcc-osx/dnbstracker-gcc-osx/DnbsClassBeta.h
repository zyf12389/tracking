// DnbsClassBeta.h
// This class is inherited from NaiveNbsClass,
// the original nonorthogonal binary subspace.
// It declares most of the methods and interfaces
// that are used in discriminative nonorthogonal
// binary subspace tracking.
//
// Reference:
// Ang Li, Feng Tang, Yanwen Guo and Hai Tao.
// Discriminative Nonorthogonal Binary Subspace
// Tracking. In: Proc. ECCV 2010.
//
// Ang Li (mailto:nju.angli@gmail.com)
// Dept. Computer Science and Technology
// Nanjing University
// Nov.8, 2010.


#pragma once
#include "headers.h"
#include "NaiveNbsClass.h"

class DnbsClassBeta
    : public NaiveNbsClass
{
public:
    DnbsClassBeta(void);
    ~DnbsClassBeta(void);

public:
    NaiveNbsClass naivenbs;
    vector<CvRect> bgsample;
	vector<CvMat*> fgrdSample;
	vector<CvMat*> bgrdSample;

public:
	// Main interface of the D-NBS class
    void computeNBS(
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
        CvMat ** previous_foregrounds = NULL,
        double coherence = -1);

    void updateNBS(
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
        CvMat ** previous_foregrounds = NULL,
        double coherence = -1);

	// Select background samples
    void sample_backgrounds(CvRect rect,
        int MaxNumBackgrounds,
        int distance_threshold,
        int x_min, int x_max,
        int y_min, int y_max,
        vector< CvRect > & pnt,
        vector< pair<double, int> > & arr,
        double ** weight,
        vector < CvRect > & samples);

	// Version 1
	void computeDNBS_v1(
		int MaxNumBases,
		double lambda,
		vector<CvMat*> foreground,
		vector<CvMat*> background);

	// Version 2: In-developing faster edition
    void computeDNBS_v2(
        int MaxNumBases,
        IplImage * frame,
        double lambda,
        vector< CvRect > & samples);

	// Calculate the reconstruction error
	double reconerr(void);
	double reconerr(IplImage * img);

	// Calculate the coefficient matrix of D-NBS
    void computeCoeffMatrix(int K = 0);


	// Previous editions for coherence pruning
	void computeNBS_discriminant(
		int MaxNumBases,
		IplImage * frame,
		CvRect rect,
		double lambda,
		vector< CvRect > & samples);

	void computeNBS_discriminant(
		int MaxNumBases,
		IplImage * frame,
		CvRect rect,
		double lambda,
		vector< CvRect > & samples,
		double coherence);

public:
    /// Extensions in the beta edition of dnbs class

	void computeDNBS_v1_pruning(
		int MaxNumBases,
		double lambda,
		vector<CvMat*> foreground,
		vector<CvMat*> background,
		int samp_width = 2,
		int samp_height = 2);

	void computeDNBS_v1_downsample(
		int MaxNumBases,
		double lambda,
		vector<CvMat*> foreground,
		vector<CvMat*> background,
		int samp_width = 2,
		int samp_height = 2);

	void computeDNBS_v1_pruning_2nd(
		int MaxNumBases,
		double lambda,
		vector<CvMat*> foreground,
		vector<CvMat*> background,
		int samp_width = 2,
		int samp_height = 2);

	double reconstruct_error(CvMat*);
	double score(double lambda);

	void genHaarBases_coherence(double mu_coherence);
};
