// Implementation of LyonLib.h
// LyonLib.h includes some useful manipulation function
// relevant to integral images and drawing rectangles under
// the OpenCV Library.
// 
// Ang Li (mailto:nju.angli@gmail.com)
// Dept. Computer Science and Technology
// Nanjing University
// Nov.8, 2010.

#include "LyonLib.h"

namespace lyon {

	static double global_timecost;
	void settime(double t) {global_timecost = t;}
	double gettime(void) {return global_timecost;}

    CvMat * image2mat(IplImage * image) {
        CvMat * ret = cvCreateMat(image->height, image->width, NBS_CV_MAT_TYPE);
        for(int i = 0; i < ret->rows; ++ i)
            for(int j = 0; j < ret->cols; ++ j)
                CV_MAT_ELEM(*ret, NBS_MAT_TYPE, i, j)
                = CV_IMAGE_ELEM(image, uchar, i, j * 3);
        return ret;
    }

    IplImage * mat2image(CvMat * mat) {
        IplImage * ret = cvCreateImage(cvSize(mat->width, mat->height), IPL_DEPTH_8U, 3);
        for(int i = 0; i < ret->height; ++ i)
            for(int j = 0; j < ret->width; ++ j)
                CV_IMAGE_ELEM(ret, uchar, i, j * 3)
                = CV_IMAGE_ELEM(ret, uchar, i, j * 3 + 1)
                = CV_IMAGE_ELEM(ret, uchar, i, j * 3 + 2)
                = (uchar)(CV_MAT_ELEM(*mat, NBS_MAT_TYPE, i, j)
                < 256 ? CV_MAT_ELEM(*mat, NBS_MAT_TYPE, i, j)  > 0 ?
                CV_MAT_ELEM(*mat, NBS_MAT_TYPE, i, j) + 0.5
                : 0 : 255);
        return ret;
    }

    CvMat * calcIntegral(CvMat * imag) {
        CvMat * ret;
        ret = cvCreateMat(imag->rows + 1, imag->cols + 1, NBS_CV_MAT_TYPE);
        cvZero(ret);
        for(int i = 0; i < imag->rows; ++ i)
            for(int j = 0; j < imag->cols; ++ j)
                CV_MAT_ELEM(*ret, NBS_MAT_TYPE, i + 1, j + 1)
                = CV_MAT_ELEM(*ret, NBS_MAT_TYPE, i, j + 1)
                - CV_MAT_ELEM(*ret, NBS_MAT_TYPE, i, j)
                + CV_MAT_ELEM(*ret, NBS_MAT_TYPE, i + 1, j)
                + CV_MAT_ELEM(*imag, NBS_MAT_TYPE, i, j);
        return ret;
    }

    CvMat * calcIntegral(IplImage * imag, int left, int right, int top, int bottom) {
        CvMat * ret;
        ret = cvCreateMat(bottom - top + 1, right - left + 1, NBS_CV_MAT_TYPE);
        cvZero(ret);
        for(int i = top; i < bottom; ++ i)
            for(int j = left; j < right; ++ j)
                CV_MAT_ELEM(*ret, NBS_MAT_TYPE, i - top + 1, j - left + 1)
                = CV_MAT_ELEM(*ret, NBS_MAT_TYPE, i - top, j - left + 1)
                - CV_MAT_ELEM(*ret, NBS_MAT_TYPE, i - top, j - left)
                + CV_MAT_ELEM(*ret, NBS_MAT_TYPE, i - top + 1, j - left)
                + (imag->nChannels == 3
				? static_cast<NBS_MAT_TYPE>(CV_IMAGE_ELEM(imag, uchar, i, j * 3))
				: static_cast<NBS_MAT_TYPE>(CV_IMAGE_ELEM(imag, uchar, i, j)));
        return ret;
    }

    CvMat * calcSqrIntegral(IplImage * imag, int left, int right, int top, int bottom) {
        CvMat * ret;
        ret = cvCreateMat(bottom - top + 1, right - left + 1, NBS_CV_MAT_TYPE);
        cvZero(ret);
        for(int i = top; i < bottom; ++ i)
            for(int j = left; j < right; ++ j)
                CV_MAT_ELEM(*ret, NBS_MAT_TYPE, i - top + 1, j - left + 1)
                = CV_MAT_ELEM(*ret, NBS_MAT_TYPE, i - top, j - left + 1)
                - CV_MAT_ELEM(*ret, NBS_MAT_TYPE, i - top, j - left)
                + CV_MAT_ELEM(*ret, NBS_MAT_TYPE, i - top + 1, j - left)
                + (imag->nChannels == 3
				? sqr(static_cast<NBS_MAT_TYPE>(CV_IMAGE_ELEM(imag, uchar, i, j * 3)))
				: sqr(static_cast<NBS_MAT_TYPE>(CV_IMAGE_ELEM(imag, uchar, i, j))));
        return ret;
    }

    double calcIntSum(CvMat * intImage, int left, int right, int top, int bottom) {
        return CV_MAT_ELEM(*intImage, NBS_MAT_TYPE, bottom, right)
            - CV_MAT_ELEM(*intImage, NBS_MAT_TYPE, top, right)
            - CV_MAT_ELEM(*intImage, NBS_MAT_TYPE, bottom, left)
            + CV_MAT_ELEM(*intImage, NBS_MAT_TYPE, top, left);
    }

    void undersize(int w, int h, int upperW, int upperH, int & newW, int & newH) {
        double scale = max(
            (double)h / upperH,
            (double)w / upperW);
        newH = static_cast<int>(h / scale);
        newW = static_cast<int>(w / scale);
    }

    IplImage * create_gray_from_rgb(const IplImage * p_image) {
        IplImage * ret
            = cvCreateImage(cvSize(p_image->width, p_image->height), p_image->depth, 3);

        for (int i = 0; i < p_image->height; ++ i)
            for (int j = 0; j < p_image->width; ++ j) {
                uchar temp
                    = static_cast<uchar>(
                    0.212671 * CV_IMAGE_ELEM(p_image, uchar, i, j * 3)
                    + 0.715160 * CV_IMAGE_ELEM(p_image, uchar, i, j * 3 + 1)
                    + 0.072169 * CV_IMAGE_ELEM(p_image, uchar, i, j * 3 + 2)
                    );
                CV_IMAGE_ELEM(ret, uchar, i, j * 3)
                    = CV_IMAGE_ELEM(ret, uchar, i, j * 3 + 1)
                    = CV_IMAGE_ELEM(ret, uchar, i, j * 3 + 2)
                    = temp;
            }

            return ret;
    }

    IplImage * selectedPatch(IplImage * pImage, CvRect rect) {
        IplImage * ret = cvCreateImage(cvSize(rect.width, rect.height), IPL_DEPTH_8U, 3);
        for(int i = 0; i < rect.width; ++ i)
            for(int j = 0; j < rect.height; ++ j) {
                CV_IMAGE_ELEM(ret, uchar, j, i * 3)
                = CV_IMAGE_ELEM(ret, uchar, j, i * 3 + 1)
                = CV_IMAGE_ELEM(ret, uchar, j, i * 3 + 2)
                = static_cast<uchar>(
                0.212671 * CV_IMAGE_ELEM(pImage, uchar, rect.y + j, (rect.x + i) * 3)
                + 0.715160 * CV_IMAGE_ELEM(pImage, uchar, rect.y + j, (rect.x + i) * 3 + 1)
                + 0.072169 * CV_IMAGE_ELEM(pImage, uchar, rect.y + j, (rect.x + i) * 3 + 2)
                );
            }
        return ret;
    }

    void drawRect(IplImage * frame, CvRect rect) {
        cvRectangle(
            frame,
            cvPoint(rect.x, rect.y),
            cvPoint(rect.x + rect.width, rect.y + rect.height),
            cvScalar(0x00, 0xff, 0x00));
    }

    void drawRect(IplImage * frame, CvRect rect, CvScalar c) {
        cvRectangle(
            frame,
            cvPoint(rect.x, rect.y),
            cvPoint(rect.x + rect.width, rect.y + rect.height),
            c, 2);
    }

    void plot_boxes(IplImage * frame, IplImage * showframe, vector<CvRect> rect) {
        for (int i = 0; i < (int)rect.size(); ++ i)
            drawRect(showframe, rect[i], cvScalar(0, 0, 0xff));
    }

    double ave_template(CvMat * patch) {
        double ave = 0;
        for (int x = 0; x < patch->width; ++ x)
            for (int y = 0; y < patch->height; ++ y)
                ave += CV_MAT_ELEM(*patch, NBS_MAT_TYPE, y, x);
        ave /= patch->width * patch->height;
        return ave;
    }

    double norm_template(CvMat * patch) {
        double l2norm = 0;
        for (int x = 0; x < patch->width; ++ x)
            for (int y = 0; y < patch->height; ++ y)
                l2norm += lyon::sqr(CV_MAT_ELEM(*patch, NBS_MAT_TYPE, y, x));
        l2norm = sqrt(l2norm);
        return l2norm;
    }

    CvMat * template_normalize(CvMat * patch) {
        CvMat * ret
            = cvCloneMat(patch);
        double ave = 0;
        for (int x = 0; x < patch->width; ++ x)
            for (int y = 0; y < patch->height; ++ y)
                ave += CV_MAT_ELEM(*ret, NBS_MAT_TYPE, y, x);
        ave /= patch->width * patch->height;
        for (int x = 0; x < patch->width; ++ x)
            for (int y = 0; y < patch->height; ++ y)
                CV_MAT_ELEM(*ret, NBS_MAT_TYPE, y, x) -= ave;
        double l2norm = 0;
        for (int x = 0; x < patch->width; ++ x)
            for (int y = 0; y < patch->height; ++ y)
                l2norm += lyon::sqr(CV_MAT_ELEM(*ret, NBS_MAT_TYPE, y, x));
        l2norm = sqrt(l2norm);
        for (int x = 0; x < patch->width; ++ x)
            for (int y = 0; y < patch->height; ++ y)
                CV_MAT_ELEM(*ret, NBS_MAT_TYPE, y, x) /= l2norm;
        return ret;
    }
}

/*
void drawIpl(CDC* pDC, IplImage* pImage, int startx, int starty, int width, int height)
{
    BITMAPINFOHEADER ptempinfoheader;
    memset(&ptempinfoheader, 0, sizeof(BITMAPINFOHEADER));

    ptempinfoheader.biHeight	=	pImage->height;
    ptempinfoheader.biWidth		=	pImage->width;
    ptempinfoheader.biSize		=	sizeof(BITMAPINFOHEADER);
    ptempinfoheader.biPlanes	=	1;
    ptempinfoheader.biBitCount	=	pImage->depth * pImage->nChannels;
    ptempinfoheader.biSizeImage	=	pImage->imageSize;

    uchar buffer[sizeof(BITMAPINFOHEADER) + 1024];
    BITMAPINFO * bmi = (BITMAPINFO*) buffer;
    int bmp_w = pImage->width, bmp_h = pImage->height;
    bmi->bmiHeader = ptempinfoheader;

    SetStretchBltMode(pDC->m_hDC, COLORONCOLOR);

    //FillBitmapInfo( bmi, bmp_w, bmp_h, 24, pImage->origin );
    RGBQUAD* palette = bmi->bmiColors;
    for(int i = 0; i < 256; i++) {
        palette[i].rgbBlue = palette[i].rgbGreen = palette[i].rgbRed = (BYTE)i;
        palette[i].rgbReserved = 0;
    }

    if (pImage->origin == IPL_ORIGIN_TL){
        pImage->origin = IPL_ORIGIN_BL;
        cvConvertImage(pImage, pImage, CV_CVTIMG_FLIP);
    }

    ::StretchDIBits(
        pDC->m_hDC,
        startx, starty, width, height,
        0, 0, pImage->width, pImage->height,
        pImage->imageData, bmi, DIB_RGB_COLORS, SRCCOPY
        );
}
*/
