// NaiveNbsClass.h
// C++ Header for Nonorthogonal Binary Subspace Decomposition
// Ang LI (mailto:nju.angli@gmail.com)
// Nov.15, 2009
//
// Attention: 
//   Codes are designed upon Intel OpenCV 2.0, and 
//   have been tested successfully on Microsoft Visual 
//   Studio 2008 under Windows XP, Nov.15 2009.
//
// References:
// [1] 'Non-orthogonal Binary Subspace and its Applications in Computer Vision',
//     Hai Tao, Ryan Crabb and Feng Tang
// [2] 'Optimised Orthogonal Matching Pursuit Approach', 
//     Laura Rebollo-Neira and David Lowe

#pragma once
#include "HaarBase.h"
#include "headers.h"

const double sqeps = 1e-18;

class NaiveNbsClass
{
public:
    NaiveNbsClass(void);
    ~NaiveNbsClass(void);

public:
    // Size of the image patch (imagW, imagH)
    int imagW, imagH;
    // Image 'patch'
    IplImage * patch;
    // List of candidate binary bases
    vector<HaarBase> listOfBases;
    
// NBS Representation
    // Indeces of selected bases
    vector<HaarBase> indexOfBases;
    // Corresponding coefficients of selected bases
    vector<NBS_MAT_TYPE> coeffOfBases;
    // Projection Coeff Matrix
    vector<CvMat*> coeffMat;
    CvMat * invPhiTPhi;
    CvMat * target;

    // Output Directory
    string workSpace;
    FILE * outfile;

    int first;

public: // Public Methods

    CvRect updateBoxSSD(IplImage * frame, CvRect rect, int margin_x, int margin_y);

    // Importing a image patch
    IplImage * importPatch(const char * path);
    IplImage * importPatch(IplImage * patch);
    CvMat * importTarget(CvMat * target);

    // Generate Haar-like Bases
    bool genHaarBases();

    // Computing the Nonorthogonal Binary Subspace
    void computeNBS(int MaxNumBases);
    void computeNBS(int MaxNumBases, double coherence);
	void computeNBS(int MaxNumBases, int samp_width, int samp_height);
    void computeNBS_efficient(int MaxNumBases);
    CvMat * calcPhi_efficient(int n, vector<CvMat*>& Phi, vector<CvMat*>& intPhi, vector<double> & denomphi);

    // Reconstruct the image from the selected bases
    IplImage * reConstruct();

    // Make a comparison between 'patch' and 'repImage'
    double compare(IplImage* repImage);

    // Pre-computing the coeff matrix PHI*(PHI^T*PHI)^-1
    void coefMatrix(int K = 0);

    vector<NBS_MAT_TYPE> compCoeff(CvMat * intImage, int K = 0);
    CvMat* compProjection(CvMat * intImage, int K = 0);

    NBS_MAT_TYPE innerProduct(CvMat * intFrame, int left, int right, int top, int bottom);

    void clear();
    void setWorkspace(const char * path);
    IplImage * showSelection();
    IplImage * showReconstruction();
    IplImage * showReconstruction(vector<double>& coeff);
    void updateCoeff(IplImage * image);

    double sum() const;
    double l2norm() const;
    double l2normScale(double hscale, double vscale) const;

    CvMat * image2mat(IplImage * image);
    void parseImage1(IplImage * bigImage, CvMat * littleMat, int tlr, int tlc);
    void parseImage2(IplImage * bigImage, CvMat * littleMat, int tlr, int tlc);
    void parseImage(IplImage * bigImage, IplImage * littleImage, int tlr, int tlc);

    void demo(void);

protected: // Private Inner Methods

    // Computing the integral image of 'imag'
    CvMat * calcIntegral(CvMat * imag);

    // Computing inner product by brute force
    NBS_MAT_TYPE innerProd(CvMat * matA, CvMat * matB);
    NBS_MAT_TYPE innerProd(CvMat * matA, CvMat * matB, int x, int y);

    // Computing Phi_K (as Gram Schmidt technique)
    // i.e. Phi_K* = Base_LK - SIGMA Phi_i <Phi_i, Base_LK> 
    //      Phi_K = Phi_K* / ||Phi_K*||
    CvMat * calcPhi(int n, vector<CvMat*>& Phi, vector<CvMat*>& intPhi);
    
    // Computing mat / x
    CvMat * calcMatDiv(CvMat * mat, NBS_MAT_TYPE x);
    
    // Output the matrix 'mat'
    void printMatrix(CvMat * mat);

    IplImage * mat2image(CvMat * mat);


protected: // Private Inner Mathematical Methods

    // Computing the square of norm of the matrix 'imag'
    NBS_MAT_TYPE SqrNorm(CvMat * imag);

public: // Members added on Oct.26, 2010
	int ind_base[200][200];

	int calc_index_haarbase(int htl, int wtl, int h, int w, int imagH, int imagW);

public:
	vector< pair<double, CvPoint> > bksimilar;
	CvRect updateBoxSSD_simp(
		IplImage * frame, CvRect rect, int left, int right, int top, int bottom, 
		CvMat * intFrame, CvMat * int2Frame, bool bSortBkSimilarity);
};
