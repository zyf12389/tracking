// HaarBase.cpp
// Implementation for a C++ Class 'HaarBase'
//   for storing Haar-like Binary Bases and
//   related computation (i.e. innerProduct)
//
// Ang LI (mailto:nju.angli@gmail.com)
// Dept. Computer Science and Technology
// Nanjing University
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
#include "headers.h"

const double eps = 1e-9;
const double inf = 1e100;
const double sigma = 1e-6;
#define NBS_CV_MAT_TYPE CV_64FC1
#define NBS_MAT_TYPE double
inline double Sqr(double x) { return x * x; }

class HaarBase
{
public:
    HaarBase(void);
    ~HaarBase(void);

public:
    int nBox;       // Number of boxes = 1 or 2
    int htl, wtl;   // Position of top-left box-pixel
    int h, w;       // Height and width of box/boxes
    NBS_MAT_TYPE norm;  // Norm of the base vector

public:
    // Constructor
    HaarBase(int nBox, int htl, int wtl, int h, int w);

    // Efficiently Computing the inner product of
    // this Haar-like base vector and an image which
    // is represented already as a integral image
    // (i.e. 'intImage'), within at most 8 add operations
    NBS_MAT_TYPE innerProd(CvMat * intImage);
	NBS_MAT_TYPE innerProd_fast(CvMat * intImage);
	NBS_MAT_TYPE innerProd(HaarBase another); // simpler edition
    NBS_MAT_TYPE innerProd(HaarBase another, int imagW);
    NBS_MAT_TYPE innerProd(CvMat * intImage, int left, int right, int top, int bottom);
    NBS_MAT_TYPE innerProdScale(CvMat * intImage, int left, int right, int top, int bottom, double hscale, double vscale);

    double int_innerProd(CvMat * intImage);
    double int_innerProd(CvMat * intImage, int left, int right, int top, int bottom);

    // Transform the base vector into a Matrix style
    CvMat * toMatrix(int row, int col) const;

    // Get the element/pixel at position (r, c),
    // imagW indicating the width of the image
    NBS_MAT_TYPE element(int r, int c, int imagW) const;
    NBS_MAT_TYPE elementdbl(double r, double c, int imagW) const;

    NBS_MAT_TYPE intersectRegion(
        int htla, int wtla, int ha, int wa,
        int htlb, int wtlb, int hb, int wb);

public: // for DEBUG
    // Print out the text description of the base
    void print() const;
};
