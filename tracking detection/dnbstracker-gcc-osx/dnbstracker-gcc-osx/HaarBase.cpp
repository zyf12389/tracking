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

#include "HaarBase.h"

using std::max;
using std::min;

HaarBase::HaarBase(void) {
    this->nBox = 0;
}

HaarBase::~HaarBase(void){}

HaarBase::HaarBase(int nBox, int htl, int wtl, int h, int w) {
    this->nBox = nBox;
    this->htl = htl;
    this->wtl = wtl;
    this->h = h;
    this->w = w;
    this->norm = (nBox == 1 ? sqrt((NBS_MAT_TYPE)h * w) : sqrt((NBS_MAT_TYPE)h * w * 2));
}

NBS_MAT_TYPE HaarBase::innerProd(CvMat * intImage) {
    if (nBox == 1)
        return (
        CV_MAT_ELEM(*intImage, NBS_MAT_TYPE, htl + h, wtl + w)
        - CV_MAT_ELEM(*intImage, NBS_MAT_TYPE, htl, wtl + w)
        - CV_MAT_ELEM(*intImage, NBS_MAT_TYPE, htl + h, wtl)
        + CV_MAT_ELEM(*intImage, NBS_MAT_TYPE, htl, wtl)
        ) / norm;
    /* for two-box Haar features
    else if (nBox == 2)
        return (
        CV_MAT_ELEM(*intImage, NBS_MAT_TYPE, htl + h - 1, wtl + w - 1)
        - (htl ? CV_MAT_ELEM(*intImage, NBS_MAT_TYPE, htl - 1, wtl + w - 1) : 0)
        - (wtl ? CV_MAT_ELEM(*intImage, NBS_MAT_TYPE, htl + h - 1, wtl - 1) : 0)
        + (htl && wtl ? CV_MAT_ELEM(*intImage, NBS_MAT_TYPE, htl - 1, wtl - 1) : 0)
        + CV_MAT_ELEM(*intImage, NBS_MAT_TYPE, htl + h - 1, intImage->width - wtl - 1)
        - (htl ? CV_MAT_ELEM(*intImage, NBS_MAT_TYPE, htl - 1, intImage->width - wtl - 1) : 0)
        - CV_MAT_ELEM(*intImage, NBS_MAT_TYPE, htl + h - 1, intImage->width - wtl - w - 1)
        + (htl ? CV_MAT_ELEM(*intImage, NBS_MAT_TYPE, htl - 1, intImage->width - wtl - w - 1) : 0)
        ) / norm;
    fprintf(stderr, "ERROR(HAARBASE::INNERPROD) : nBox > 2\n");
    */
    return -1;
}

NBS_MAT_TYPE HaarBase::innerProd_fast(CvMat * intImage) {
    return (
        CV_MAT_ELEM(*intImage, NBS_MAT_TYPE, htl + h, wtl + w)
        - CV_MAT_ELEM(*intImage, NBS_MAT_TYPE, htl, wtl + w)
        - CV_MAT_ELEM(*intImage, NBS_MAT_TYPE, htl + h, wtl)
        + CV_MAT_ELEM(*intImage, NBS_MAT_TYPE, htl, wtl)
        ) / norm;
}


double HaarBase::int_innerProd(CvMat * intImage) {
    if (nBox == 1)
        return (
        CV_MAT_ELEM(*intImage, int, htl + h, wtl + w)
        - CV_MAT_ELEM(*intImage, int, htl, wtl + w )
        - CV_MAT_ELEM(*intImage, int, htl + h, wtl)
        + CV_MAT_ELEM(*intImage, int, htl, wtl)
        ) / norm;
    /* for 2-box Haar features
    else if (nBox == 2)
        return (
        CV_MAT_ELEM(*intImage, int, htl + h - 1, wtl + w - 1)
        - (htl ? CV_MAT_ELEM(*intImage, int, htl - 1, wtl + w - 1) : 0)
        - (wtl ? CV_MAT_ELEM(*intImage, int, htl + h - 1, wtl - 1) : 0)
        + (htl && wtl ? CV_MAT_ELEM(*intImage, int, htl - 1, wtl - 1) : 0)
        + CV_MAT_ELEM(*intImage, int, htl + h - 1, intImage->width - wtl - 1)
        - (htl ? CV_MAT_ELEM(*intImage, int, htl - 1, intImage->width - wtl - 1) : 0)
        - CV_MAT_ELEM(*intImage, int, htl + h - 1, intImage->width - wtl - w - 1)
        + (htl ? CV_MAT_ELEM(*intImage, int, htl - 1, intImage->width - wtl - w - 1) : 0)
        ) / norm;
    fprintf(stderr, "ERROR(HAARBASE::INNERPROD) : nBox > 2\n");
    */
    return -1;
}

NBS_MAT_TYPE HaarBase::innerProd(CvMat * intImage, int left, int right, int top, int bottom) {
    if (nBox == 1)
        return (
        CV_MAT_ELEM(*intImage, NBS_MAT_TYPE, top + htl + h, left + wtl + w)
        - CV_MAT_ELEM(*intImage, NBS_MAT_TYPE, top + htl, left + wtl + w)
        - CV_MAT_ELEM(*intImage, NBS_MAT_TYPE, top + htl + h, left + wtl)
        + CV_MAT_ELEM(*intImage, NBS_MAT_TYPE, top + htl, left + wtl)
        ) / norm;
    /* For 2-box Haar features
    else if (nBox == 2)
        return (
        CV_MAT_ELEM(*intImage, NBS_MAT_TYPE, top + htl + h - 1, left + wtl + w - 1)
        - (top + htl ? CV_MAT_ELEM(*intImage, NBS_MAT_TYPE, top + htl - 1, left + wtl + w - 1) : 0)
        - (left + wtl ? CV_MAT_ELEM(*intImage, NBS_MAT_TYPE, top + htl + h - 1, left + wtl - 1) : 0)
        + (top + htl && left + wtl ? CV_MAT_ELEM(*intImage, NBS_MAT_TYPE, top + htl - 1, left + wtl - 1) : 0)
        + CV_MAT_ELEM(*intImage, NBS_MAT_TYPE, top + htl + h - 1, right - wtl - 1)
        - (top + htl ? CV_MAT_ELEM(*intImage, NBS_MAT_TYPE, top + htl - 1, right - wtl - 1) : 0)
        - CV_MAT_ELEM(*intImage, NBS_MAT_TYPE, top + htl + h - 1, right - wtl - w - 1)
        + (top + htl ? CV_MAT_ELEM(*intImage, NBS_MAT_TYPE, top + htl - 1, right - wtl - w - 1) : 0)
        ) / norm;
    fprintf(stderr, "ERROR(HAARBASE::INNERPROD) : nBox > 2\n");
    */
    return -1;
}

NBS_MAT_TYPE HaarBase::innerProdScale(CvMat * intImage,
                                      int left, int right,
                                      int top, int bottom,
                                      double hscale, double vscale) {
    int newhtl = (int)(this->htl * vscale);
    int newwtl = (int)(this->wtl * hscale);
    int neww = (int)((this->w + this->wtl) * hscale + 0.5) - newwtl;
    int newh = (int)((this->h + this->htl) * vscale + 0.5) - newhtl;
    if (top + newhtl + newh - 1 >= bottom)
        printf ("bottom %d >= %d, %d\n", top + newhtl + newh - 1, bottom,
        (int)((this->h + this->htl - 1) * vscale));
    if (left + newwtl + neww - 1 >= right)
        printf ("right %d >= %d, %d\n", left + newwtl + neww - 1, right,
        (int)((this->w + this->wtl - 1) * hscale));
    return (
    CV_MAT_ELEM(*intImage, NBS_MAT_TYPE, top + newhtl + newh, left + newwtl + neww)
    - CV_MAT_ELEM(*intImage, NBS_MAT_TYPE, top + newhtl, left + newwtl + neww)
    - CV_MAT_ELEM(*intImage, NBS_MAT_TYPE, top + newhtl + newh, left + newwtl)
    + CV_MAT_ELEM(*intImage, NBS_MAT_TYPE, top + newhtl, left + newwtl)
    ) / norm;
}

double HaarBase::int_innerProd(CvMat * intImage, int left, int right, int top, int bottom) {
    if (nBox == 1)
        return (
        CV_MAT_ELEM(*intImage, int, top + htl + h, left + wtl + w)
        - CV_MAT_ELEM(*intImage, int, top + htl, left + wtl + w)
        - CV_MAT_ELEM(*intImage, int, top + htl + h, left + wtl)
        + CV_MAT_ELEM(*intImage, int, top + htl, left + wtl)
        ) / norm;
    /* For 2-box Haar features
    else if (nBox == 2)
        return (
        CV_MAT_ELEM(*intImage, int, top + htl + h - 1, left + wtl + w - 1)
        - (top + htl ? CV_MAT_ELEM(*intImage, int, top + htl - 1, left + wtl + w - 1) : 0)
        - (left + wtl ? CV_MAT_ELEM(*intImage, int, top + htl + h - 1, left + wtl - 1) : 0)
        + (top + htl && left + wtl ? CV_MAT_ELEM(*intImage, int, top + htl - 1, left + wtl - 1) : 0)
        + CV_MAT_ELEM(*intImage, int, top + htl + h - 1, right - wtl - 1)
        - (top + htl ? CV_MAT_ELEM(*intImage, int, top + htl - 1, right - wtl - 1) : 0)
        - CV_MAT_ELEM(*intImage, int, top + htl + h - 1, right - wtl - w - 1)
        + (top + htl ? CV_MAT_ELEM(*intImage, int, top + htl - 1, right - wtl - w - 1) : 0)
        ) / norm;
    fprintf(stderr, "ERROR(HAARBASE::INNERPROD) : nBox > 2\n");
    */
    return -1;
}

CvMat * HaarBase::toMatrix(int row, int col) const {
	// Transform to matrix
    CvMat * ret = cvCreateMat(row, col, NBS_CV_MAT_TYPE);
    for(int i = 0; i < ret->rows; ++ i)
        for(int j = 0; j < ret->cols; ++ j)
            CV_MAT_ELEM(*ret, NBS_MAT_TYPE, i, j) = 0;
    if (this->nBox == 1)
        for(int i = htl; i < htl + h; ++ i)
            for(int j = wtl; j < wtl + w; ++ j)
                CV_MAT_ELEM(*ret, NBS_MAT_TYPE, i, j) = 1. / norm;
    else if (this->nBox == 2)
        for(int i = htl; i < htl + h; ++ i)
            for(int j = wtl; j < wtl + w; ++ j) {
                CV_MAT_ELEM(*ret, NBS_MAT_TYPE, i, j)
                = CV_MAT_ELEM(*ret, NBS_MAT_TYPE, i, col - j - 1)
                = 1. / norm;
            }
    else fprintf(stderr, "ERROR(HAARBASE::TOMATRIX): NBOX > 2\n");
    return ret;
}

NBS_MAT_TYPE HaarBase::element(int r, int c, int imagW) const {
	// Get the pixel in the feature
    if (nBox == 1)
        return (r >= htl && r < htl + h && c >= wtl && c < wtl + w ? 1./norm : 0);
    else return (r >= htl && r < htl + h && (c >= wtl && c < wtl + w || c >= imagW - wtl - w && c < imagW - wtl) ? 1./norm : 0);
}

NBS_MAT_TYPE HaarBase::elementdbl(double r, double c, int imagW) const {
    if (nBox == 1)
        return (r >= htl && r < htl + h && c >= wtl && c < wtl + w ? 1./norm : 0);
    else return (r >= htl && r < htl + h && (c >= wtl && c < wtl + w || c >= imagW - wtl - w && c < imagW - wtl) ? 1./norm : 0);
}

void HaarBase::print() const{
    printf("nBox = %d, htl = %d, wtl = %d, h = %d, w = %d\n", nBox, htl, wtl, h, w);
}


NBS_MAT_TYPE HaarBase::intersectRegion(
    int htla, int wtla, int ha, int wa,
    int htlb, int wtlb, int hb, int wb) {
		// Return the area of intersection of two Haar features
    int left = max(wtla, wtlb);
    int right = min(wtla + wa, wtlb + wb);
    int top = max(htla, htlb);
    int bottom = min(htla + ha, htlb + hb);
    if (right > left && bottom > top)
        return (right - left) * (bottom - top);
    return 0;
}

NBS_MAT_TYPE HaarBase::innerProd(HaarBase another, int imagW) {
	// Calculate the inner product between two Haar features
    if (nBox == 1)
        if (another.nBox == 1)
        return intersectRegion(
                htl, wtl, h, w,
                another.htl, another.wtl, another.h, another.w
            ) / norm / another.norm;
        else return (
            intersectRegion(
                htl, wtl, h, w,
                another.htl, another.wtl, another.h, another.w)
            +
            intersectRegion(
                htl, wtl, h, w,
                another.htl, imagW - another.wtl - another.w, another.h, another.w)
            ) / norm / another.norm;
    else if (another.nBox == 1)
        return (
            intersectRegion(
                another.htl, another.wtl, another.h, another.w,
                htl, wtl, h, w)
            +
            intersectRegion(
                another.htl, another.wtl, another.h, another.w,
                htl, imagW - wtl - w, h, w)
            ) / norm / another.norm;
    else return intersectRegion(
                htl, wtl, h, w,
                another.htl, another.wtl, another.h, another.w)
                * 2 / norm / another.norm;
}


NBS_MAT_TYPE HaarBase::innerProd(HaarBase another) {
	return intersectRegion(
		htl, wtl, h, w,
		another.htl, another.wtl, another.h, another.w
		) / norm / another.norm;
}