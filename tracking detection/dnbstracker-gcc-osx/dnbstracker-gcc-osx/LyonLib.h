
// LyonLib.h includes some useful manipulation function
// relevant to integral images and drawing rectangles under
// the OpenCV Library.
// 
// Ang Li (mailto:nju.angli@gmail.com)
// Dept. Computer Science and Technology
// Nanjing University
// Nov.8, 2010.

#pragma once
#include "headers.h"

namespace lyon {
    /* Non-orthogonal Binary Subspace */
#define NBS_CV_MAT_TYPE CV_64FC1
#define NBS_MAT_TYPE double

    CvMat * image2mat(IplImage * image);
    IplImage * mat2image(CvMat * mat);
    CvMat * calcIntegral(CvMat * imag);
    CvMat * calcIntegral(IplImage * imag, int left, int right, int top, int bottom);
    CvMat * calcSqrIntegral(IplImage * imag, int left, int right, int top, int bottom);
    double calcIntSum(CvMat * intImage, int left, int right, int top, int bottom);

    /* Standard Utilities */
    const double inf = 1e100;
    template<typename T> inline T sqr(T x) { return x * x; }
    bool createFolder(const char * m_strFolderPath);
    void undersize(int w, int h, int upperW, int upperH, int & newW, int & newH);
    
    IplImage * create_gray_from_rgb(const IplImage * p_image);
    IplImage * selectedPatch(IplImage * pImage, CvRect rect);
    void drawRect(IplImage * frame, CvRect rect);
    void drawRect(IplImage * frame, CvRect rect, CvScalar c);
    void plot_boxes(IplImage * frame, IplImage * showframe, std::vector<CvRect> rect);

    CvMat * template_normalize(CvMat * patch);
    double ave_template(CvMat * patch);
    double norm_template(CvMat * patch);

    //void drawIpl(CDC* pDC, IplImage* pImage, int startx, int starty, int width, int height);

	void settime(double t);
	double gettime(void);
}
