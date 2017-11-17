// NaiveNbsClass.cpp
// Implementation for a C++ Class 'NaiveNbsClass'
//   focused on Nonorthogonal Binary Subspace Decomposition
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

#include "NaiveNbsClass.h"
#include <time.h>
#include <stdlib.h>
#include "LyonLib.h"

using std::max;
using std::min;

const int ImageSizeThreshold = 100;

NaiveNbsClass::NaiveNbsClass(void) {
    this->invPhiTPhi = NULL;
    this->patch = NULL;
    this->outfile = NULL;
    this->workSpace = "";
    this->target = NULL;
}

NaiveNbsClass::~NaiveNbsClass(void) {
    for(int i = 0; i < (int)coeffMat.size(); ++ i)
        if (coeffMat[i] != NULL)
            cvReleaseMat(&coeffMat[i]);
    if (this->patch)
        cvReleaseImage(&patch);
    if (this->invPhiTPhi)
        cvReleaseMat(&invPhiTPhi);
    if (this->outfile)
        fclose(outfile);
    if (this->target != NULL)
        cvReleaseMat(&target);
}

void NaiveNbsClass::clear() {
    this->imagH = 0;
    this->imagW = 0;
    this->indexOfBases.clear();
    this->listOfBases.clear();
    this->coeffOfBases.clear();

    for(int i = 0; i < (int)coeffMat.size(); ++ i)
        if (coeffMat[i] != NULL)
            cvReleaseMat(&coeffMat[i]);
    coeffMat.clear();

    if (invPhiTPhi != NULL)
        cvReleaseMat(&invPhiTPhi);
    if (outfile != NULL)
        fclose(outfile);
    if (patch != NULL)
        cvReleaseImage(&patch);
}

IplImage * NaiveNbsClass::importPatch(const char * path) {
    // Importing the image patch
    if (patch != NULL)
        cvReleaseImage(&patch);
    patch = cvLoadImage(path);
    this->imagH = patch->height;
    this->imagW = patch->width;
    if (this->imagH > ImageSizeThreshold
        || this->imagW > ImageSizeThreshold) {
            printf("ERROR: Image is too large\n");
            return NULL;
    }
    if (target != NULL)
        cvReleaseMat(&target);
    target = image2mat(patch);
    return cvCloneImage(patch);
}

IplImage * NaiveNbsClass::importPatch(IplImage * patch) {
    if (this->patch != NULL)
        cvReleaseImage(&this->patch);
    this->patch = cvCloneImage(patch);
    this->imagH = patch->height;
    this->imagW = patch->width;
    if (this->imagH > ImageSizeThreshold
        || this->imagW > ImageSizeThreshold) {
            printf("ERROR: Image is too large\n");
            return NULL;
    }
    if (target != NULL)
        cvReleaseMat(&target);
    target = image2mat(patch);
    return patch;
}

CvMat * NaiveNbsClass::calcIntegral(CvMat * imag) {
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

NBS_MAT_TYPE NaiveNbsClass::innerProd(CvMat * matA, CvMat * matB) {
    NBS_MAT_TYPE ret = 0;
    for(int i = 0; i < matA->rows; ++ i)
        for(int j = 0; j < matA->cols; ++ j)
            ret += CV_MAT_ELEM(*matA, NBS_MAT_TYPE, i, j)
                * CV_MAT_ELEM(*matB, NBS_MAT_TYPE, i, j);
    return ret;
}

NBS_MAT_TYPE NaiveNbsClass::innerProd(CvMat * matA, CvMat * matB, int x, int y) {
    NBS_MAT_TYPE ret = 0;
    for(int i = 0; i < matA->rows; ++ i)
        for(int j = 0; j < matA->cols; ++ j)
            ret += CV_MAT_ELEM(*matA, NBS_MAT_TYPE, i, j)
            * CV_MAT_ELEM(*matB, NBS_MAT_TYPE, y + i, x + j);
    return ret;
}

CvMat * NaiveNbsClass::calcPhi(int n, vector<CvMat*>& Phi, vector<CvMat*>& intPhi) {
//    printf("n = %d\n", n);
    CvMat * ret = listOfBases[n].toMatrix(imagH, imagW);
    for(int i = 0; i < (int)Phi.size(); ++ i) {
        NBS_MAT_TYPE temp = listOfBases[n].innerProd(intPhi[i]);
        for(int r = 0; r < ret->rows; ++ r)
            for(int c = 0; c < ret->cols; ++ c)
                CV_MAT_ELEM(*ret, NBS_MAT_TYPE, r, c)
                -= CV_MAT_ELEM(*Phi[i], NBS_MAT_TYPE, r, c) * temp;
    }
    return ret;
}

CvMat * NaiveNbsClass::calcPhi_efficient(int n, vector<CvMat*>& Phi, vector<CvMat*>& intPhi, vector<double> & denomphi) {
    //    printf("n = %d\n", n);
    CvMat * ret = listOfBases[n].toMatrix(imagH, imagW);
    for(int i = 0; i < (int)Phi.size(); ++ i) {
        NBS_MAT_TYPE temp = listOfBases[n].innerProd(intPhi[i]);
        for(int r = 0; r < ret->rows; ++ r)
            for(int c = 0; c < ret->cols; ++ c)
                CV_MAT_ELEM(*ret, NBS_MAT_TYPE, r, c)
                -= CV_MAT_ELEM(*Phi[i], NBS_MAT_TYPE, r, c) * temp / denomphi[i];
    }
    return ret;
}


CvMat * NaiveNbsClass::calcMatDiv(CvMat * mat, NBS_MAT_TYPE x) {
    CvMat * ret = cvCloneMat(mat);
    for(int i = 0; i < ret->rows; ++ i)
        for(int j = 0; j < ret->cols; ++ j)
            CV_MAT_ELEM(*ret, NBS_MAT_TYPE, i, j) /= x;
    return ret;
}

void NaiveNbsClass::printMatrix(CvMat * mat) {
    for(int i = 0; i < mat->rows; ++ i) {
        for(int j = 0; j < mat->cols; ++ j)
            printf("%.2lf ", CV_MAT_ELEM(*mat, NBS_MAT_TYPE, i, j));
        putchar('\n');
    }
}

void NaiveNbsClass::computeNBS(int MaxNumBases) {

// Construction of Nonorthogonal Binary Space
//      by Optimized Orthogonal Matching Pursuit(OOMP)
// Forward Biorthogonalization Technique
//      [Revised Edition] with Complexity O(KWH(WH+K))
// Ang LI (mailto:nju.angli@gmail.com)
// Nov.15, 2009
//
// References:
// [1] 'Non-orthogonal Binary Subspace and its Applications in Computer Vision',
//     Hai Tao, Ryan Crabb and Feng Tang
// [2] 'Optimised Orthogonal Matching Pursuit Approach',
//     Laura Rebollo-Neira and David Lowe

// Generate all the Haar-like bases, stored in 'listOfBases'
//    genHaarBases();
//    printf("HaarBases Been Generated.\n");
//    cvWaitKey();

// Definitions
    // Matrix of the original image
    CvMat * matPatch;
    // Integral image of the original image
    CvMat * intPatch;
    // K is the number of selected bases
    int K;
    // N is the total number of candidate bases
    int N = listOfBases.size();

    // B and D are two arrays of parameters in deciding the
    // optimized index of the next to-be-selected binary base.
    // i.e. selecting the one that maximising E_n(=|B_n|^2/D_n)
    NBS_MAT_TYPE * B = new NBS_MAT_TYPE[N];
    NBS_MAT_TYPE * D = new NBS_MAT_TYPE[N];

    // Indicating whether the base is selected to form the subspace
    bool * selected = new bool[N];

    // Index of each selected binary bases
    vector<int> L;
    // Projection Coefficients
    vector<NBS_MAT_TYPE> Coeff;

    // Residue[K] stands for the residue when subspace is
    // spanned by the first K binary bases.
    // i.e. Residue[K] = f - Proj_V[K](f)
    vector<NBS_MAT_TYPE> Residue;

    // Phi's and their integral image
    vector<CvMat*> Phi;
    vector<CvMat*> intPhi;

    // Inner product of each Phi[K] and image F(i.e. matPatch)
    vector<NBS_MAT_TYPE> innerPhiF;

    // Beta's and their integral image
    vector<CvMat*> Beta;
    vector<CvMat*> intBeta;

    // Temporal memories
    double curMaxE;
    int tempL;

// Pre-processing, converting image 'patch' to type of matrix
    matPatch = cvCloneMat(target);
/*
    matPatch = cvCreateMat(patch->height, patch->width, NBS_CV_MAT_TYPE);
    for(int i = 0; i < matPatch->rows; ++ i)
        for(int j = 0; j < matPatch->cols; ++ j)
        // GRAY = 0.212671*R + 0.715160*G + 0.072169*B
            CV_MAT_ELEM(*matPatch, NBS_MAT_TYPE, i, j)
            = 0.212671 * CV_IMAGE_ELEM(patch, uchar, i, j * 3)
            + 0.715160 * CV_IMAGE_ELEM(patch, uchar, i, j * 3 + 1)
            + 0.072169 * CV_IMAGE_ELEM(patch, uchar, i, j * 3 + 2);
*/
// Initialization
    //printf("Initialization Begin!\n");
    double dtime = clock();
    // Reset memories
    memset(selected, false, sizeof(bool) * N);
    // Computing the integral image of the original image
    intPatch = calcIntegral(matPatch);
    // Process of initializing parameters and get the first choice
    tempL = 0;
    for(int i = 0; i < N; ++ i) {
        B[i] = listOfBases[i].innerProd(intPatch);
        D[i] = 1;
        if (fabs(B[i]) > fabs(B[tempL])) tempL = i;
    }
    L.push_back(tempL);
    selected[tempL] = true;
    //printf ("Standard L[0] = %d\n", L[0]);

    // Computing parameters
    Phi.push_back(listOfBases[L[0]].toMatrix(imagH, imagW));
    intPhi.push_back(calcIntegral(Phi[0]));
    Beta.push_back(cvCloneMat(Phi[0]));
    intBeta.push_back(calcIntegral(Beta[0]));
    Coeff.push_back(B[L[0]]);
    //printf("norm = %.3lf\n", cvNorm(matPatch));
    Residue.push_back(SqrNorm(matPatch) - Sqr(Coeff[0]));
    //printf("Rep.%d Residue[%d] = %.10lf, Coeff[%d] = %.3lf\n", 0, 0, sqrt(Residue[0]), 0, Coeff[0]);
    //listOfBases[L[0]].print();
    innerPhiF.push_back(innerProd(Phi[0], matPatch));
    //outfile = fopen("nbslogs.txt", "w+");
    //fprintf(outfile, "numbase\terror\ttime\n");
    //fprintf(outfile, "%d\t%.3lf\t%.3lf\n", 1, sqrt(Residue[0] / imagH / imagW), (clock() - dtime) / (double)CLOCKS_PER_SEC);
    K = 1;

//    fprintf(outfile, "Size of Image Template = %d width, %d height in pixel\n", imagW, imagH);
//    fprintf(outfile, "Number of bases to be selected = %d\n", MaxNumBases);
//    fprintf(outfile, "Time for pre-processing = %.3lf (seconds)\n", (clock() - dtime) / CLOCKS_PER_SEC);

    //int trackbar = 0;
    //cvNamedWindow("Computing NBS");
    //cvCreateTrackbar("nbs", "Computing NBS", &trackbar, MaxNumBases, NULL);
// Main Selecting Procedures :
    //printf("Main Procedures Begin!\n");
    // Total Time Complexity is O(KN + K^2WH) = O(KWH(WH+K))

    double mu_coherent = 1;
    dtime = clock();
    int * next = new int[N + 1];
    for (int i = 0; i < N; ++ i)
        next[i] = i + 1;
    next[N] = 0;

    while (K < MaxNumBases && Residue[K - 1] > sigma) {
        //cvSetTrackbarPos("nbs", "Computing NBS", K);
        // O(N + KWH) for each repetition
        // Selecting the best fit Base in O(N)
        curMaxE = -inf; tempL = -1;
        NBS_MAT_TYPE tempPhiBase, tempE;
        for(int i = next[N]; i < N; i = next[i])
            if (!selected[i]) {
                tempPhiBase = listOfBases[i].innerProd(intPhi[K - 1]);
                B[i] -= innerPhiF[K - 1] * tempPhiBase;
                D[i] -= Sqr(tempPhiBase);
                tempE = fabs(B[i]) < eps ? 0 : (Sqr(B[i]) / D[i]);
                if (tempE > curMaxE) {
                    curMaxE = tempE;
                    tempL = i;
                }
            }
        L.push_back(tempL);
        selected[tempL] = true;

        if (mu_coherent < 1) {
            int pre_i = N;

            for (int i = next[N]; i < N; i = next[i]) {
                if (i == tempL
                    || listOfBases[i].innerProd(listOfBases[tempL], imagW)
                    >= mu_coherent) {
                        next[pre_i] = next[i];
                } else pre_i = next[i];
            }
        }

        // Computing parameters
        Residue.push_back(Residue[K - 1] - curMaxE);        // O(1)
        CvMat * tempMat = calcPhi(L[K], Phi, intPhi);       // O(KWH)
        Phi.push_back(calcMatDiv(tempMat, sqrt(D[L[K]])));  // O(WH)
        Beta.push_back(calcMatDiv(tempMat, D[L[K]]));       // O(WH)
        intBeta.push_back(calcIntegral(Beta[K]));           // O(WH)
        cvReleaseMat(&tempMat);
        intPhi.push_back(calcIntegral(Phi[K]));             // O(WH)
        innerPhiF.push_back(innerProd(Phi[K], matPatch));   // O(WH)
        Coeff.push_back(B[L[K]] / D[L[K]]);                 // O(1)
        //printf("Rep.%d Residue[%d] = %.10lf, Coeff[%d] = %.3lf\n", K, K, sqrt(Residue[K]), K, Coeff[K]);

        // Post-processing at the end of each selection
        NBS_MAT_TYPE temp;
        for(int i = 0; i < K; ++ i) {   // O(KWH)
            temp = listOfBases[L[K]].innerProd(intBeta[i]);
            for(int r = 0; r < Beta[i]->rows; ++ r)
                for(int c = 0; c < Beta[i]->cols; ++ c)
                    CV_MAT_ELEM(*Beta[i], NBS_MAT_TYPE, r, c)
                    -= CV_MAT_ELEM(*Beta[K], NBS_MAT_TYPE, r, c) * temp;
            // 'temp' is real, so its conjugate is itself;
            Coeff[i] -= temp * Coeff[K];
            cvReleaseMat(&intBeta[i]);
            intBeta[i] = calcIntegral(Beta[i]);
        }

        //fprintf(outfile, "%d\t%.3lf\t%.3lf\n", K + 1, sqrt(Residue[K] / imagW / imagH), (clock() - dtime) / (double)CLOCKS_PER_SEC);
        // Now, get the next binary base
        ++ K;
    }

    delete [] next;
//    printf("It's been done! time = %.2lf:D\n", (clock() - dtime) / (double)CLOCKS_PER_SEC);
    /*
    char out[100];
    sprintf(out, "%.2lfs", (clock() - dtime) / CLOCKS_PER_SEC);
    AfxMessageBox(out);
    */

//    fprintf(outfile, "Time for selecting bases = %.3lf (seconds)\n", (clock() - dtime) / (double)CLOCKS_PER_SEC);
    // Copy Answers
    this->indexOfBases.clear();
    for(int i = 0; i < (int)L.size(); ++ i)
        this->indexOfBases.push_back(listOfBases[L[i]]);
    this->coeffOfBases = Coeff;

    //fclose(outfile);

    // Post-Processing for efficient computation of projection
//    this->coefMatrix();

    // Output each selection
 //   fprintf(outfile, "The L2 distance between the original and reconstructed image vector = %.3lf\n", sqrt(Residue[K - 1]));
    /*
    for(int i = 0; i < K; ++ i) {
        printf("Rep.%d Residue[%d] = %.10lf, Coeff[%d] = %.3lf\n", K, K, sqrt(Residue[K]), K, Coeff[K]);
    */

    // FOR DEBUG: Checking Projection

//    printf("calc project\n");
//    CvMat * proj = this->compProjection(intPatch);
//    printf("norm = %.3lf\n", cvNorm(proj, matPatch));
//    IplImage * imagproj = mat2image(proj);
//    cvNamedWindow("hehehe");
//    cvShowImage("hehehe", imagproj);
//    vector<NBS_MAT_TYPE> co = this->compCoeff(intPatch);
/*
    for(int i = 0; i < (int)co.size(); ++ i)
        if (fabs(co[i] - this->coeffOfBases[i]) > 4)
            printf("%.3lf\n", co[i] - this->coeffOfBases[i]); */
//    compare(imagproj);
//    cvWaitKey();
//    cvDestroyWindow("hehehe");
//    cvReleaseMat(&proj);
//    cvReleaseImage(&imagproj);
    //reConstruct();

    // Release Memories
    cvReleaseMat(&matPatch);
    cvReleaseMat(&intPatch);
    delete [] B;
    delete [] D;
    delete [] selected;
    for(int i = 0; i < (int)Phi.size(); ++ i) {
        cvReleaseMat(&Phi[i]);
        cvReleaseMat(&intPhi[i]);
        cvReleaseMat(&Beta[i]);
        cvReleaseMat(&intBeta[i]);
    }
}

void NaiveNbsClass::computeNBS(int MaxNumBases, double coherence) {

    // Construction of Nonorthogonal Binary Space
    //      by Optimized Orthogonal Matching Pursuit(OOMP)
    // Forward Biorthogonalization Technique
    //      [Revised Edition] with Complexity O(KWH(WH+K))
    // Ang LI (mailto:nju.angli@gmail.com)
    // Nov.15, 2009
    //
    // References:
    // [1] 'Non-orthogonal Binary Subspace and its Applications in Computer Vision',
    //     Hai Tao, Ryan Crabb and Feng Tang
    // [2] 'Optimised Orthogonal Matching Pursuit Approach',
    //     Laura Rebollo-Neira and David Lowe

    // Generate all the Haar-like bases, stored in 'listOfBases'
    //    genHaarBases();
    //    printf("HaarBases Been Generated.\n");
    //    cvWaitKey();

    // Definitions
    // Matrix of the original image
    CvMat * matPatch;
    // Integral image of the original image
    CvMat * intPatch;
    // K is the number of selected bases
    int K;
    // N is the total number of candidate bases
    int N = listOfBases.size();

    // B and D are two arrays of parameters in deciding the
    // optimized index of the next to-be-selected binary base.
    // i.e. selecting the one that maximising E_n(=|B_n|^2/D_n)
    NBS_MAT_TYPE * B = new NBS_MAT_TYPE[N];
    NBS_MAT_TYPE * D = new NBS_MAT_TYPE[N];

    // Indicating whether the base is selected to form the subspace
    bool * selected = new bool[N];

    // Index of each selected binary bases
    vector<int> L;
    // Projection Coefficients
    vector<NBS_MAT_TYPE> Coeff;

    // Residue[K] stands for the residue when subspace is
    // spanned by the first K binary bases.
    // i.e. Residue[K] = f - Proj_V[K](f)
    vector<NBS_MAT_TYPE> Residue;

    // Phi's and their integral image
    vector<CvMat*> Phi;
    vector<CvMat*> intPhi;

    // Inner product of each Phi[K] and image F(i.e. matPatch)
    vector<NBS_MAT_TYPE> innerPhiF;

    // Beta's and their integral image
    vector<CvMat*> Beta;
    vector<CvMat*> intBeta;

    // Temporal memories
    double curMaxE;
    int tempL;

    // Pre-processing, converting image 'patch' to type of matrix
    matPatch = cvCloneMat(target);
    /*
    matPatch = cvCreateMat(patch->height, patch->width, NBS_CV_MAT_TYPE);
    for(int i = 0; i < matPatch->rows; ++ i)
    for(int j = 0; j < matPatch->cols; ++ j)
    // GRAY = 0.212671*R + 0.715160*G + 0.072169*B
    CV_MAT_ELEM(*matPatch, NBS_MAT_TYPE, i, j)
    = 0.212671 * CV_IMAGE_ELEM(patch, uchar, i, j * 3)
    + 0.715160 * CV_IMAGE_ELEM(patch, uchar, i, j * 3 + 1)
    + 0.072169 * CV_IMAGE_ELEM(patch, uchar, i, j * 3 + 2);
    */
    // Initialization
    //printf("Initialization Begin!\n");
    double dtime = clock();
    // Reset memories
    memset(selected, false, sizeof(bool) * N);
    // Computing the integral image of the original image
    intPatch = calcIntegral(matPatch);
    // Process of initializing parameters and get the first choice
    tempL = 0;
    for(int i = 0; i < N; ++ i) {
        B[i] = listOfBases[i].innerProd(intPatch);
        D[i] = 1;
        if (fabs(B[i]) > fabs(B[tempL])) tempL = i;
    }
    L.push_back(tempL);
    selected[tempL] = true;
    //printf ("Standard L[0] = %d\n", L[0]);

    // Computing parameters
    Phi.push_back(listOfBases[L[0]].toMatrix(imagH, imagW));
    intPhi.push_back(calcIntegral(Phi[0]));
    Beta.push_back(cvCloneMat(Phi[0]));
    intBeta.push_back(calcIntegral(Beta[0]));
    Coeff.push_back(B[L[0]]);
    Residue.push_back(SqrNorm(matPatch) - Sqr(Coeff[0]));
    innerPhiF.push_back(innerProd(Phi[0], matPatch));
    K = 1;

    //int trackbar = 0;
    //double mu_coherent = 1;
    dtime = clock();
    int * next = new int[N + 1];
    int prev = N;
    for (int i = 0; i < N; ++ i) {
        double temp = listOfBases[i].innerProd(listOfBases[L[0]], imagW);
        if (temp <= coherence + eps) {
            next[prev] = i;
            prev = i;
        }
    }
    next[prev] = -1;

    while (K < MaxNumBases && Residue[K - 1] > sigma) {
        //cvSetTrackbarPos("nbs", "Computing NBS", K);
        // O(N + KWH) for each repetition
        // Selecting the best fit Base in O(N)
        curMaxE = -inf; tempL = -1;
        NBS_MAT_TYPE tempPhiBase, tempE;
        for(int i = next[N]; i != -1; i = next[i]) {
            tempPhiBase = listOfBases[i].innerProd(intPhi[K - 1]);
            B[i] -= innerPhiF[K - 1] * tempPhiBase;
            D[i] -= Sqr(tempPhiBase);
            tempE = fabs(B[i]) < eps ? 0 : (Sqr(B[i]) / D[i]);
            if (tempE > curMaxE) {
                curMaxE = tempE;
                tempL = i;
            }
        }

//         for(int i = 0; i < N; ++ i)
//             if (!selected[i]) {
//                 if (listOfBases[i].innerProd(listOfBases[L[K-1]], imagW) > coherence) {
//                     selected[i] = true;
//                     continue;
//                 }
//                 tempPhiBase = listOfBases[i].innerProd(intPhi[K - 1]);
//                 B[i] -= innerPhiF[K - 1] * tempPhiBase;
//                 D[i] -= Sqr(tempPhiBase);
//                 tempE = fabs(B[i]) < eps ? 0 : (Sqr(B[i]) / D[i]);
//                 if (tempE > curMaxE) {
//                     curMaxE = tempE;
//                     tempL = i;
//                 }
//             }

        L.push_back(tempL);
        selected[tempL] = true;

        int prev = N;
        for (int i = next[N]; i != -1; i = next[i]) {
            double temp = listOfBases[i].innerProd(listOfBases[tempL], imagW);
            if (temp <= coherence + eps) {
                next[prev] = i;
                prev = i;
            }
        }
        next[prev] = -1;

        // Computing parameters
        Residue.push_back(Residue[K - 1] - curMaxE);        // O(1)
        CvMat * tempMat = calcPhi(L[K], Phi, intPhi);       // O(KWH)
        Phi.push_back(calcMatDiv(tempMat, sqrt(D[L[K]])));  // O(WH)
        Beta.push_back(calcMatDiv(tempMat, D[L[K]]));       // O(WH)
        intBeta.push_back(calcIntegral(Beta[K]));           // O(WH)
        cvReleaseMat(&tempMat);
        intPhi.push_back(calcIntegral(Phi[K]));             // O(WH)
        innerPhiF.push_back(innerProd(Phi[K], matPatch));   // O(WH)
        Coeff.push_back(B[L[K]] / D[L[K]]);                 // O(1)
        //printf("Rep.%d Residue[%d] = %.10lf, Coeff[%d] = %.3lf\n", K, K, sqrt(Residue[K]), K, Coeff[K]);

        // Post-processing at the end of each selection
        NBS_MAT_TYPE temp;
        for(int i = 0; i < K; ++ i) {   // O(KWH)
            temp = listOfBases[L[K]].innerProd(intBeta[i]);
            for(int r = 0; r < Beta[i]->rows; ++ r)
                for(int c = 0; c < Beta[i]->cols; ++ c)
                    CV_MAT_ELEM(*Beta[i], NBS_MAT_TYPE, r, c)
                    -= CV_MAT_ELEM(*Beta[K], NBS_MAT_TYPE, r, c) * temp;
            // 'temp' is real, so its conjugate is itself;
            Coeff[i] -= temp * Coeff[K];
            cvReleaseMat(&intBeta[i]);
            intBeta[i] = calcIntegral(Beta[i]);
        }

        //fprintf(outfile, "%d\t%.3lf\t%.3lf\n", K + 1, sqrt(Residue[K] / imagW / imagH), (clock() - dtime) / (double)CLOCKS_PER_SEC);
        // Now, get the next binary base
        ++ K;
    }

    delete [] next;

    //    fprintf(outfile, "Time for selecting bases = %.3lf (seconds)\n", (clock() - dtime) / (double)CLOCKS_PER_SEC);
    // Copy Answers
    this->indexOfBases.clear();
    for(int i = 0; i < (int)L.size(); ++ i)
        this->indexOfBases.push_back(listOfBases[L[i]]);
    this->coeffOfBases = Coeff;

    // Release Memories
    cvReleaseMat(&matPatch);
    cvReleaseMat(&intPatch);
    delete [] B;
    delete [] D;
    delete [] selected;
    for(int i = 0; i < (int)Phi.size(); ++ i) {
        cvReleaseMat(&Phi[i]);
        cvReleaseMat(&intPhi[i]);
        cvReleaseMat(&Beta[i]);
        cvReleaseMat(&intBeta[i]);
    }
}

void NaiveNbsClass::computeNBS_efficient(int MaxNumBases) {

    // Construction of Nonorthogonal Binary Space
    //      by Optimized Orthogonal Matching Pursuit(OOMP)
    // Ang LI (mailto:nju.angli@gmail.com)
    // Mar.14, 2010

    // Definitions
    // Matrix of the original image
    CvMat * matPatch;
    // Integral image of the original image
    CvMat * intPatch;
    // K is the number of selected bases
    int K;
    // N is the total number of candidate bases
    int N = listOfBases.size();

    // D[i] = ||gamma_i||^2 is the denominator item for i-th base
    double * D = new double[N];
    // Index of each selected binary bases
    vector<int> L;
    // Phi's and their integral image
    vector<CvMat*> Phi;
    vector<CvMat*> intPhi;

    double max_item;
    int opt_base;
    double inner_base_residue;
    double inner_base_phi;
    double curr_item;
    vector<double> denomphi;
    vector<double> nomphi;
    CvMat * fcurr = cvCloneMat(target);
    CvMat * int_fcurr = calcIntegral(fcurr);

    // Pre-processing, converting image 'patch' to type of matrix
    matPatch = cvCloneMat(target);
    // Computing the integral image of the original image
    intPatch = calcIntegral(matPatch);
    // Process of initializing parameters and get the first choice
    opt_base = 0;
    max_item = -inf;
    double nom_item;
    for(int i = 0; i < N; ++ i) {
        inner_base_residue
            = listOfBases[i].innerProd(int_fcurr);
        D[i] = 1;
        if (fabs(inner_base_residue) > max_item) {
            max_item = fabs(inner_base_residue);
            nom_item = inner_base_residue;
            opt_base = i;
        }
    }
    L.push_back(opt_base);

    // Computing parameters
    denomphi.push_back(D[opt_base]);
    nomphi.push_back(nom_item);
    Phi.push_back(listOfBases[L[0]].toMatrix(imagH, imagW));
    intPhi.push_back(calcIntegral(Phi[0]));

    K = 1;
    double mu_coherent = 1;
    int * next = new int[N + 1];
    for (int i = 0; i < N; ++ i)
        next[i] = i + 1;
    next[N] = 0;

    double temp_coef;

    double timeInSecs
        = clock();
    //double totaltime = 0;

    while (K < MaxNumBases) {

        temp_coef = nomphi[K - 1] / denomphi[K - 1];
        for (int r = 0; r < fcurr->rows; ++ r) {
            for (int c = 0; c < fcurr->cols; ++ c) {
                CV_MAT_ELEM(*fcurr, double, r, c)
                    -= CV_MAT_ELEM(*Phi[K - 1], double, r, c)
                    * temp_coef;
            }
        }
        cvReleaseMat(&int_fcurr);
        int_fcurr = calcIntegral(fcurr);

        max_item = -inf;
        opt_base = -1;
        double nom_item;
        for(int i = next[N]; i < N; i = next[i]) {
            inner_base_residue
                = listOfBases[i].innerProd(int_fcurr);
            inner_base_phi
                = listOfBases[i].innerProd(intPhi[K - 1]);
            D[i] -= Sqr(inner_base_phi) / denomphi[K - 1];
            curr_item = fabs(inner_base_residue) < eps
                ? 0 : (Sqr(inner_base_residue) / D[i]);
            if (curr_item > max_item) {
                nom_item = inner_base_residue;
                max_item = curr_item;
                opt_base = i;
            }
        }
        //printf ("opt-base = %d, max_item = %lf\n", opt_base, max_item);
        L.push_back(opt_base);
        denomphi.push_back(D[opt_base]);
        nomphi.push_back(nom_item);
        /*
        for (int i = next[N]; next[i] != opt_base; i = next[i]) {
            if (next[i] == opt_base)
                next[i] = next[opt_base];
        }*/

        if (mu_coherent < 1) {
            for (int pre_i = N, i = next[N]; i < N; i = next[i]) {
                if (i == opt_base
                    || listOfBases[i].innerProd(listOfBases[opt_base], imagW)
                    >= mu_coherent) {
                        next[pre_i] = next[i];
                } else pre_i = i;
            }
        }

        // Computing parameters
        Phi.push_back(calcPhi_efficient(L[K], Phi, intPhi, denomphi));          // O(KWH)
        intPhi.push_back(calcIntegral(Phi[K]));             // O(WH)

        ++ K;
    }


    timeInSecs
        = (clock() - timeInSecs) / CLOCKS_PER_SEC;
    printf ("time cost = %.10lf seconds\n", timeInSecs);

    // Copy Answers
    this->indexOfBases.clear();
    for(int i = 0; i < (int)L.size(); ++ i)
        this->indexOfBases.push_back(listOfBases[L[i]]);
    //this->coefMatrix();
    //this->coeffOfBases = this->compCoeff(intPatch);

    // Release Memories
    cvReleaseMat(&matPatch);
    cvReleaseMat(&intPatch);
    cvReleaseMat(&fcurr);
    cvReleaseMat(&int_fcurr);
    delete [] next;
    delete [] D;
    for(int i = 0; i < (int)Phi.size(); ++ i) {
        cvReleaseMat(&Phi[i]);
        cvReleaseMat(&intPhi[i]);
    }
}


NBS_MAT_TYPE NaiveNbsClass::SqrNorm(CvMat * imag) {
    NBS_MAT_TYPE ret = 0;
    for(int i = 0; i < imag->height; ++ i)
        for(int j = 0; j < imag->width; ++ j)
            ret += Sqr(CV_MAT_ELEM(*imag, NBS_MAT_TYPE, i, j));
    return ret;
}

int NaiveNbsClass::calc_index_haarbase(int htl, int wtl, int h, int w, int imagH, int imagW) {
	return ind_base[w][h] + wtl * (imagH - h + 1) + htl;
}

bool NaiveNbsClass::genHaarBases() {
    // Generate all the possible single-box and
    // vertical symmetric two-box Haar-like bases

    // Initialization
    imagW = target->width;
    imagH = target->height;
    int & w = imagW;
    int & h = imagH;
    listOfBases.clear();

	memset(ind_base, 0, sizeof(ind_base));
	// Single-box Haar-like functions
	for(int wf = 1; wf <= w; ++ wf)
		for(int hf = 1; hf <= h; ++ hf) {
			ind_base[wf][hf] = listOfBases.size();
			// the box is of width 'wf' and height 'hf'
			for(int wtl = 0; wtl < w - wf + 1; ++ wtl)
				for(int htl = 0; htl < h - hf + 1; ++ htl) {
					// top-left pixel is positioned at (htl, wtl)
					listOfBases.push_back(
						HaarBase(1, htl, wtl, hf, wf)
						);
					if (htl == 0 && wtl == 0 && hf == h && wf == w)
						first = listOfBases.size() - 1;
				}
		}

    // Vertically symmetric two-box Haar-like functions
        /*
    for(int wf = 1; wf <= w / 2; ++ wf)
        for(int hf = 1; hf <= h; ++ hf) {
            // the left box is of width 'wf' and height 'hf'
            // it must be on the left half of the whole image
            for(int wtl = 0; wtl < w / 2 - wf + 1; ++ wtl)
                for(int htl = 0; htl < h - hf + 1; ++ htl) {
                    // top-left pixel is positioned at (htl, wtl)
                    listOfBases.push_back(
                        HaarBase(2, htl, wtl, hf, wf)
                    );
                }
        }
*/
    // All bases are stored in vector 'listOfBases'
//    printf("num of bases = %d\n", listOfBases.size());
    return true;
}

IplImage * NaiveNbsClass::reConstruct() {
    CvMat * temp = cvCreateMat(imagH, imagW, NBS_CV_MAT_TYPE);
    for(int i = 0; i < temp->rows; ++ i)
        for(int j = 0; j < temp->cols; ++ j)
            CV_MAT_ELEM(*temp, NBS_MAT_TYPE, i, j) = 0;
    for(int i = 0; i < (int)indexOfBases.size(); ++ i) {
        for(int r = 0; r < temp->rows; ++ r)
            for(int c = 0; c < temp->cols; ++ c)
                CV_MAT_ELEM(*temp, NBS_MAT_TYPE, r, c)
                += coeffOfBases[i] * indexOfBases[i].element(r, c, imagW);
//        printf("Coeff[%d] = %.3lf\n", i, coeffOfBases[i]);
//        listOfBases[indexOfBases[i]].print();
    }
//    printMatrix(temp);
    IplImage * ret = mat2image(temp);
    /*cvCreateImage(cvGetSize(this->patch), patch->depth, patch->nChannels);
    for(int i = 0; i < ret->height; ++ i)
        for(int j = 0; j < ret->width; ++ j)
            CV_IMAGE_ELEM(ret, uchar, i, j * 3)
            = CV_IMAGE_ELEM(ret, uchar, i, j * 3 + 1)
            = CV_IMAGE_ELEM(ret, uchar, i, j * 3 + 2)
            = (uchar)(CV_MAT_ELEM(*temp, NBS_MAT_TYPE, i, j) + 0.5);
            */
    //compare(ret);
    cvReleaseMat(&temp);
    return ret;
}

CvMat * NaiveNbsClass::image2mat(IplImage * image) {
    CvMat * ret = cvCreateMat(image->height, image->width, NBS_CV_MAT_TYPE);
    for(int i = 0; i < ret->rows; ++ i)
        for(int j = 0; j < ret->cols; ++ j)
            CV_MAT_ELEM(*ret, NBS_MAT_TYPE, i, j)
            = CV_IMAGE_ELEM(image, uchar, i, j * 3);
    return ret;
}

double NaiveNbsClass::compare(IplImage * repImage) {
    IplImage * subImage = cvCreateImage(cvGetSize(patch), patch->depth, patch->nChannels);
    cvSub(repImage, patch, subImage);
    cvNamedWindow("subImage");
    cvShowImage("subImage", subImage);
    cvWaitKey();
    cvDestroyWindow("subImage");
    cvReleaseImage(&subImage);
    double dist = 0;
//    double N1 = cvNorm(patch);
//    double N2 = cvNorm(repImage);
//    printf("N1 = %.3lf, N2 = %.3lf\n", N1, N2);
    for(int i = 0; i < repImage->width; ++ i)
        for(int j = 0; j < repImage->height; ++ j)
            dist += Sqr((double)CV_IMAGE_ELEM(patch, uchar, j, i * 3)
            - (double)CV_IMAGE_ELEM(repImage, uchar, j, i * 3));
/*
    FILE * fout;
    fout = fopen("D:/research/NBS&BPCA/repimage.txt", "w+");
    for(int j = 0; j < repImage->height; ++ j){
        for(int i = 0; i < repImage->width; ++ i)
            fprintf(fout, "%d\t", CV_IMAGE_ELEM(repImage, uchar, j, i * 3));
        fprintf(fout, "\n");
    }
    fclose(fout);
*/
    //dist = sqrt(dist);
    printf("distinction between 2 images = %.3lf\n", sqrt(dist));
//    printf("distinction between 2 images = %.3lf\n", cvNorm(image2mat(patch), image2mat(repImage)));
    return dist;
}

void NaiveNbsClass::coefMatrix(int K) {
    if (K == 0)
        K = this->indexOfBases.size();
    CvMat * PhiTPhi = cvCreateMat(K, K, NBS_CV_MAT_TYPE);

    for(int i = 0; i < K; ++ i)
        for(int j = 0; j < K; ++ j)
            CV_MAT_ELEM(*PhiTPhi, NBS_MAT_TYPE, i, j)
            = indexOfBases[i].innerProd(indexOfBases[j], imagW);

    if (invPhiTPhi != NULL) {
        cvReleaseMat(&invPhiTPhi);
        invPhiTPhi = NULL;
    }
    invPhiTPhi = cvCreateMat(K, K, NBS_CV_MAT_TYPE);

    // Compute (Phi^T*Phi)^-1
    cvInvert(PhiTPhi, invPhiTPhi);

    cvReleaseMat(&PhiTPhi);
    for(int i = 0; i < (int)this->coeffMat.size(); ++ i)
        cvReleaseMat(&this->coeffMat[i]);
    this->coeffMat.clear();
    for(int i = 0; i < K; ++ i) {
        CvMat * temp
            = cvCreateMat(imagH, imagW, NBS_CV_MAT_TYPE);
        cvZero(temp);
        for(int j = 0; j < K; ++ j) {
            NBS_MAT_TYPE scale
                = CV_MAT_ELEM(*invPhiTPhi, NBS_MAT_TYPE, j, i);
            for(int r = 0; r < temp->rows; ++ r)
                for(int c = 0; c < temp->cols; ++ c)
                    CV_MAT_ELEM(*temp, NBS_MAT_TYPE, r, c)
                    += indexOfBases[j].element(r, c, temp->cols) * scale;
        }
        this->coeffMat.push_back(temp);
    }
}

vector<NBS_MAT_TYPE> NaiveNbsClass::compCoeff(CvMat *intImage, int K) {
    if (K == 0)
        K = indexOfBases.size();
    vector<NBS_MAT_TYPE> PhiX;
    for(int i = 0; i < K; ++ i)
        PhiX.push_back(indexOfBases[i].innerProd(intImage));
    vector<NBS_MAT_TYPE> ret;
    for(int i = 0; i < K; ++ i) {
        NBS_MAT_TYPE temp = 0;
        for(int j = 0; j < K; ++ j)
            temp += PhiX[j] * CV_MAT_ELEM(*invPhiTPhi, NBS_MAT_TYPE, i, j);
        ret.push_back(temp);
    }
    return ret;
}

void NaiveNbsClass::updateCoeff(IplImage * image) {
    if (coeffMat.size() == 0)
        coefMatrix();
    CvMat * imag = image2mat(image);
    CvMat * intImage = calcIntegral(imag);
    coeffOfBases = compCoeff(intImage);
    cvReleaseMat(&imag);
    cvReleaseMat(&intImage);
}

CvMat * NaiveNbsClass::compProjection(CvMat * intImage, int K) {
    if (K == 0)
        K = indexOfBases.size();
    vector<NBS_MAT_TYPE> PhiX;
    for(int i = 0; i < K; ++ i)
        PhiX.push_back(indexOfBases[i].innerProd(intImage));
    CvMat * ret = cvCreateMat(imagH, imagW, NBS_CV_MAT_TYPE);
    cvZero(ret);
    for(int i = 0; i < K; ++ i) {
        for(int r = 0; r < imagH; ++ r)
            for(int c = 0; c < imagW; ++ c)
                CV_MAT_ELEM(*ret, NBS_MAT_TYPE, r, c)
                += CV_MAT_ELEM(*coeffMat[i], NBS_MAT_TYPE, r, c)
                * PhiX[i];
    }
    return ret;
}

IplImage * NaiveNbsClass::mat2image(CvMat * mat) {
    IplImage * ret = cvCreateImage(cvSize(mat->width, mat->height), IPL_DEPTH_8U, 3);
    for(int i = 0; i < ret->height; ++ i)
        for(int j = 0; j < ret->width; ++ j)
            CV_IMAGE_ELEM(ret, uchar, i, j * 3)
            = CV_IMAGE_ELEM(ret, uchar, i, j * 3 + 1)
            = CV_IMAGE_ELEM(ret, uchar, i, j * 3 + 2)
            = (uchar)(CV_MAT_ELEM(*mat, NBS_MAT_TYPE, i, j)
            < 255 ? CV_MAT_ELEM(*mat, NBS_MAT_TYPE, i, j)  > 0 ?
            CV_MAT_ELEM(*mat, NBS_MAT_TYPE, i, j) + 0.5
            : 0 : 255);
    return ret;
}

void NaiveNbsClass::parseImage1(IplImage * bigImage, CvMat * littleMat, int tlr, int tlc) {
    for(int r = 0; r < littleMat->rows; ++ r)
        for(int c = 0; c < littleMat->cols; ++ c)
            CV_IMAGE_ELEM(bigImage, uchar, tlr + r, (tlc + c) * 3)
            = CV_IMAGE_ELEM(bigImage, uchar, tlr + r, (tlc + c) * 3 + 1)
            = CV_IMAGE_ELEM(bigImage, uchar, tlr + r, (tlc + c) * 3 + 2)
            = CV_MAT_ELEM(*littleMat, NBS_MAT_TYPE, r, c) > eps ? 255 : 0;
}

void NaiveNbsClass::parseImage2(IplImage * bigImage, CvMat * littleMat, int tlr, int tlc) {
    for(int r = 0; r < littleMat->rows; ++ r)
        for(int c = 0; c < littleMat->cols; ++ c)
            CV_IMAGE_ELEM(bigImage, uchar, tlr + r, (tlc + c) * 3)
            = CV_IMAGE_ELEM(bigImage, uchar, tlr + r, (tlc + c) * 3 + 1)
            = CV_IMAGE_ELEM(bigImage, uchar, tlr + r, (tlc + c) * 3 + 2)
            = (uchar)CV_MAT_ELEM(*littleMat, NBS_MAT_TYPE, r, c);
}

void NaiveNbsClass::parseImage(IplImage * bigImage, IplImage * littleImage, int tlr, int tlc) {
    for(int r = 0; r < littleImage->height; ++ r)
        for(int c = 0; c < littleImage->width; ++ c)
            CV_IMAGE_ELEM(bigImage, uchar, tlr + r, (tlc + c) * 3)
            = CV_IMAGE_ELEM(bigImage, uchar, tlr + r, (tlc + c) * 3 + 1)
            = CV_IMAGE_ELEM(bigImage, uchar, tlr + r, (tlc + c) * 3 + 2)
            = (uchar)CV_IMAGE_ELEM(littleImage, uchar, r, c * 3);
}

IplImage * NaiveNbsClass::showSelection() {
    int numrow = 10;
    int margin = 2;
    int rows = (indexOfBases.size() + numrow - 1) / numrow + 1;
    int cols = numrow;
    IplImage * svimage = cvCreateImage(cvSize((imagW + margin) * cols - margin, (imagH + margin) * rows - margin), IPL_DEPTH_8U, 3);

    for(int i = 0; i < svimage->height; ++ i)
        for(int j = 0; j < svimage->width; ++ j)
            CV_IMAGE_ELEM(svimage, uchar, i, j * 3)
            = CV_IMAGE_ELEM(svimage, uchar, i, j * 3 + 1)
            = CV_IMAGE_ELEM(svimage, uchar, i, j * 3 + 2)
            = 0;

    for (int i = 1; i < cols; ++ i)
        for(int j = 0; j < svimage->height; ++ j)
            for(int k = 0; k < margin; ++ k) {
                CV_IMAGE_ELEM(svimage, uchar, j, (i * (imagW + margin) - margin + k) * 3)
                = 200;
                CV_IMAGE_ELEM(svimage, uchar, j, (i * (imagW + margin) - margin + k) * 3 + 1)
                = 255;
                CV_IMAGE_ELEM(svimage, uchar, j, (i * (imagW + margin) - margin + k) * 3 + 2)
                = 200;
            }

    for (int i = 1; i < rows; ++ i)
        for(int j = 0; j < svimage->width; ++ j)
            for(int k = 0; k < margin; ++ k) {
                CV_IMAGE_ELEM(svimage, uchar, i * (imagH + margin) - margin + k, j * 3)
                    = 200;
                CV_IMAGE_ELEM(svimage, uchar, i * (imagH + margin) - margin + k, j * 3 + 1)
                    = 255;
                CV_IMAGE_ELEM(svimage, uchar, i * (imagH + margin) - margin + k, j * 3 + 2)
                    = 200;
            }

    parseImage(svimage, this->patch, 0, 0);

    for(int i = 0; i < (int)indexOfBases.size(); ++ i) {
        int tlc = (i % 10) * (imagW + margin);
        int tlr = (i / 10 + 1) * (imagH + margin);

        CvMat * mat = indexOfBases[i].toMatrix(imagH, imagW);
        parseImage1(svimage, mat, tlr, tlc);
    }

    return svimage;
}

IplImage * NaiveNbsClass::showReconstruction() {
    int numrow = 10;
    int margin = 2;
    int rows = (indexOfBases.size() + numrow - 1) / numrow + 1;
    int cols = numrow;
    IplImage * svimage = cvCreateImage(cvSize((imagW + margin) * cols - margin, (imagH + margin) * rows - margin), IPL_DEPTH_8U, 3);

    for(int i = 0; i < svimage->height; ++ i)
        for(int j = 0; j < svimage->width; ++ j)
            CV_IMAGE_ELEM(svimage, uchar, i, j * 3)
            = CV_IMAGE_ELEM(svimage, uchar, i, j * 3 + 1)
            = CV_IMAGE_ELEM(svimage, uchar, i, j * 3 + 2)
            = 0;

    for (int i = 1; i < cols; ++ i)
        for(int j = 0; j < svimage->height; ++ j)
            for(int k = 0; k < margin; ++ k) {
                CV_IMAGE_ELEM(svimage, uchar, j, (i * (imagW + margin) - margin + k) * 3)
                = 200;
                CV_IMAGE_ELEM(svimage, uchar, j, (i * (imagW + margin) - margin + k) * 3 + 1)
                = 255;
                CV_IMAGE_ELEM(svimage, uchar, j, (i * (imagW + margin) - margin + k) * 3 + 2)
                = 200;
            }

    for (int i = 1; i < rows; ++ i)
        for(int j = 0; j < svimage->width; ++ j)
            for(int k = 0; k < margin; ++ k) {
                CV_IMAGE_ELEM(svimage, uchar, i * (imagH + margin) - margin + k, j * 3)
                    = 200;
                CV_IMAGE_ELEM(svimage, uchar, i * (imagH + margin) - margin + k, j * 3 + 1)
                    = 255;
                CV_IMAGE_ELEM(svimage, uchar, i * (imagH + margin) - margin + k, j * 3 + 2)
                    = 200;
            }

    parseImage(svimage, this->patch, 0, 0);

    CvMat * mat = cvCloneMat(target);
    CvMat * intImage = this->calcIntegral(mat);
    cvReleaseMat(&mat);
    for(int i = 0; i < (int)indexOfBases.size(); ++ i) {
        int tlc = (i % 10) * (imagW + margin);
        int tlr = (i / 10 + 1) * (imagH + margin);

        this->coefMatrix(i + 1);
        CvMat *temp = this->compProjection(intImage, i + 1);

        this->parseImage(svimage, mat2image(temp), tlr, tlc);
        cvReleaseMat(&temp);
    }

    cvReleaseMat(&intImage);
    return svimage;
}

// Setting the workspace of current nbs-class
void NaiveNbsClass::setWorkspace(const char * path) {
    system((string("mkdir \"") + path + "\"").c_str());
    this->workSpace = path;
    outfile = fopen((workSpace + "/info.txt").c_str(), "w+");
}

IplImage * NaiveNbsClass::showReconstruction(vector<double>& coeff) {
    int numrow = 10;
    int margin = 2;
    int rows = (coeff.size() + numrow - 1) / numrow + 1;
    int cols = numrow;
    IplImage * svimage = cvCreateImage(cvSize((imagW + margin) * cols - margin, (imagH + margin) * rows - margin), IPL_DEPTH_8U, 3);

    for(int i = 0; i < svimage->height; ++ i)
        for(int j = 0; j < svimage->width; ++ j)
            CV_IMAGE_ELEM(svimage, uchar, i, j * 3)
            = CV_IMAGE_ELEM(svimage, uchar, i, j * 3 + 1)
            = CV_IMAGE_ELEM(svimage, uchar, i, j * 3 + 2)
            = 0;

    for (int i = 1; i < cols; ++ i)
        for(int j = 0; j < svimage->height; ++ j)
            for(int k = 0; k < margin; ++ k) {
                CV_IMAGE_ELEM(svimage, uchar, j, (i * (imagW + margin) - margin + k) * 3)
                = 200;
                CV_IMAGE_ELEM(svimage, uchar, j, (i * (imagW + margin) - margin + k) * 3 + 1)
                = 255;
                CV_IMAGE_ELEM(svimage, uchar, j, (i * (imagW + margin) - margin + k) * 3 + 2)
                = 200;
            }

    for (int i = 1; i < rows; ++ i)
        for(int j = 0; j < svimage->width; ++ j)
            for(int k = 0; k < margin; ++ k) {
                CV_IMAGE_ELEM(svimage, uchar, i * (imagH + margin) - margin + k, j * 3)
                    = 200;
                CV_IMAGE_ELEM(svimage, uchar, i * (imagH + margin) - margin + k, j * 3 + 1)
                    = 255;
                CV_IMAGE_ELEM(svimage, uchar, i * (imagH + margin) - margin + k, j * 3 + 2)
                    = 200;
            }

    parseImage(svimage, this->patch, 0, 0);

    CvMat * temp = cvCreateMat(imagH, imagW, NBS_CV_MAT_TYPE);
    for(int i = 0; i < temp->rows; ++ i)
        for(int j = 0; j < temp->cols; ++ j)
            CV_MAT_ELEM(*temp, NBS_MAT_TYPE, i, j) = 0;

    CvMat * mat = this->image2mat(patch);
    CvMat * intImage = this->calcIntegral(mat);
    cvReleaseMat(&mat);
    for(int i = 0; i < (int)coeff.size(); ++ i) {

        for(int r = 0; r < temp->rows; ++ r)
            for(int c = 0; c < temp->cols; ++ c)
                CV_MAT_ELEM(*temp, NBS_MAT_TYPE, r, c)
                += coeff[i] * indexOfBases[i].element(r, c, imagW);

        int tlc = (i % 10) * (imagW + margin);
        int tlr = (i / 10 + 1) * (imagH + margin);

        this->parseImage(svimage, mat2image(temp), tlr, tlc);
    }

    cvReleaseMat(&intImage);
    return svimage;
}

// NaiveNbsClass::sum
// Calculate the sum of pixel-levels of the nbs-approximated image
double NaiveNbsClass::sum() const {
    double ret = 0;
    for (int i = 0; i < (int)indexOfBases.size(); ++ i)
        ret += coeffOfBases[i] * indexOfBases[i].norm;
    return ret;
}

// NaiveNbsClass::l2norm
// Calculate ||I||^2, where I is the approximated image
double NaiveNbsClass::l2norm() const {
    double ret = 0;
    for (int i = 0; i < imagH; ++ i)
        for (int j = 0; j < imagW; ++ j) {
            double temp = 0;
            for (int k = 0; k < (int)indexOfBases.size(); ++ k)
                temp += coeffOfBases[k] * indexOfBases[k].element(i, j, imagW);
            ret += temp * temp;
        }
    return ret;
}

// Compute L2-Norm with scale changes
double NaiveNbsClass::l2normScale(double hscale, double vscale) const {
    double ret = 0;
    int new_imagH = (int)(imagH * vscale + 0.5);
    int new_imagW = (int)(imagW * hscale + 0.5);
    for (int i = 0; i < new_imagH; ++ i) {
        for (int j = 0; j < new_imagW; ++ j) {
            double temp = 0;
            for (int k = 0; k < (int)indexOfBases.size(); ++ k)
                temp += coeffOfBases[k] * indexOfBases[k].elementdbl(i / hscale, j / vscale, imagW);
            ret += temp * temp;
        }
    }
    return ret;
}


NBS_MAT_TYPE NaiveNbsClass::innerProduct(CvMat * intFrame, int left, int right, int top, int bottom) {
    NBS_MAT_TYPE ret = 0;
    for(int i = 0; i < (int)coeffOfBases.size(); ++ i)
        ret += coeffOfBases[i] * indexOfBases[i].innerProd(intFrame, left, right, top, bottom);
    return ret;
}

CvMat * NaiveNbsClass::importTarget(CvMat * target) {
    imagH = target->height;
    imagW = target->width;
    if (this->target != NULL)
        cvReleaseMat(&this->target);
    this->target = cvCloneMat(target);
    if (patch != NULL)
        cvReleaseImage(&patch);
    patch = mat2image(target);
    return this->target;
}

void NaiveNbsClass::demo(void) {
    IplImage * reconstruction
        = this->showReconstruction();
    IplImage * selection
        = this->showSelection();
    cvNamedWindow("DisNbsClass::Demo::Reconstruction");
    cvNamedWindow("DisNbsClass::Demo::Selection");
    cvShowImage("DisNbsClass::Demo::Reconstruction", reconstruction);
    cvShowImage("DisNbsClass::Demo::Selection", selection);
    cvWaitKey();
    cvDestroyWindow("DisNbsClass::Demo::Reconstruction");
    cvDestroyWindow("DisNbsClass::Demo::Selection");
}

CvRect NaiveNbsClass::updateBoxSSD(IplImage * frame, CvRect rect, int margin_x, int margin_y) {
    int left = max(0, rect.x - margin_x);
    int right = min(frame->width - rect.width, rect.x + margin_x);
    int top = max(0, rect.y - margin_y);
    int bottom = min(frame->height - rect.height, rect.y + margin_y);

    CvMat * intFrame
        = lyon::calcIntegral(
        frame, left, right + rect.width,
        top, bottom + rect.height);

    CvMat * int2Frame
        = lyon::calcSqrIntegral(
        frame, left, right + rect.width,
        top, bottom + rect.height);

    double globalssd = lyon::inf;
    CvRect ret;

    for(int r = top; r < bottom; ++ r)
        for(int c = left; c < right; ++ c) {
            double ssd = 0;
            ssd += lyon::calcIntSum(int2Frame, c - left, c - left + rect.width, r - top, r - top + rect.height);
            for(int i = 0; i < (int)coeffOfBases.size(); ++ i)
                ssd -= 2. * coeffOfBases[i] * indexOfBases[i].innerProd(intFrame, c - left, c - left + rect.width, r - top, r - top + rect.height);

            if (ssd < globalssd) {
                globalssd = ssd;
                ret = cvRect(c, r, rect.width, rect.height);
            }
        }

        cvReleaseMat(&intFrame);
        cvReleaseMat(&int2Frame);
        return ret;
}

bool compareBkSimilar(const pair<double, CvPoint> & A, const pair<double, CvPoint> & B) {
	return A.first < B.first;
}

CvRect NaiveNbsClass::updateBoxSSD_simp(
	IplImage * frame, CvRect rect, int left, int right, int top, int bottom, 
	CvMat * intFrame, CvMat * int2Frame, bool bSortBkSimilarity) 
{
	double globalssd = lyon::inf;
	CvRect ret;

	if (bSortBkSimilarity)
		bksimilar.clear();
	for(int r = top; r < bottom; ++ r)
		for(int c = left; c < right; ++ c) {
			double ssd = 0;
			ssd += lyon::calcIntSum(int2Frame, c - left, c - left + rect.width, r - top, r - top + rect.height);
			for(int i = 0; i < (int)coeffOfBases.size(); ++ i)
				ssd -= 2. * coeffOfBases[i] * indexOfBases[i].innerProd(intFrame, c - left, c - left + rect.width, r - top, r - top + rect.height);
			if (bSortBkSimilarity)
				bksimilar.push_back(make_pair(ssd, cvPoint(c, r)));
			if (ssd < globalssd) {
				globalssd = ssd;
				ret = cvRect(c, r, rect.width, rect.height);
			}
		}

	if (bSortBkSimilarity)
		sort(bksimilar.begin(), bksimilar.end(), compareBkSimilar);
	printf("Search SSD = %lf\n", (globalssd + l2norm()) / target->width / target->height);
	return ret;
}

void NaiveNbsClass::computeNBS(int MaxNumBases, int samp_width, int samp_height) {

	// Ang Li @ Nanjing, China. Jan. 26, 2011.

	// Definitions
	// Matrix of the original image
	CvMat * matPatch;
	// Integral image of the original image
	CvMat * intPatch;
	// K is the number of selected bases
	int K;
	// N is the total number of candidate bases
	int N = listOfBases.size();

	// B and D are two arrays of parameters in deciding the
	// optimized index of the next to-be-selected binary base.
	// i.e. selecting the one that maximising E_n(=|B_n|^2/D_n)
	NBS_MAT_TYPE * B = new NBS_MAT_TYPE[N];
	NBS_MAT_TYPE * D = new NBS_MAT_TYPE[N];

	// Indicating whether the base is selected to form the subspace
	bool * selected = new bool[N];

	// Index of each selected binary bases
	vector<int> L;
	// Projection Coefficients
	vector<NBS_MAT_TYPE> Coeff;

	// Residue[K] stands for the residue when subspace is
	// spanned by the first K binary bases.
	// i.e. Residue[K] = f - Proj_V[K](f)
	vector<NBS_MAT_TYPE> Residue;

	// Phi's and their integral image
	vector<CvMat*> Phi;
	vector<CvMat*> intPhi;

	// Inner product of each Phi[K] and image F(i.e. matPatch)
	vector<NBS_MAT_TYPE> innerPhiF;

	// Beta's and their integral image
	vector<CvMat*> Beta;
	vector<CvMat*> intBeta;

	// Temporal memories
	double curMaxE;
	int tempL;

	// Pre-processing, converting image 'patch' to type of matrix
	matPatch = cvCloneMat(target);

	/// Select representative features from the dictionary
	vector<int> sampIndex;
	for (int i = 0; i < N; ++ i) {
		if (listOfBases[i].htl % samp_height == 0
			&& listOfBases[i].wtl % samp_width == 0
			&& listOfBases[i].w % samp_width == 0
			&& listOfBases[i].h % samp_height == 0)
			sampIndex.push_back(i);
	}
	int sampN = sampIndex.size();

	int * progress = new int[N];

	// Initialization
	double dtime = clock();
	// Reset memories
	memset(selected, false, sizeof(bool) * N);
	// Computing the integral image of the original image
	intPatch = calcIntegral(matPatch);
	// Process of initializing parameters and get the first choice
	tempL = 0;
	for(int i = 0; i < N; ++ i) {
		progress[i] = 0;
		B[i] = listOfBases[i].innerProd(intPatch);
		D[i] = 1;
		if (fabs(B[i]) > fabs(B[tempL])) tempL = i;
	}
	L.push_back(tempL);
	selected[tempL] = true;

	// Computing parameters
	Phi.push_back(listOfBases[L[0]].toMatrix(imagH, imagW));
	intPhi.push_back(calcIntegral(Phi[0]));
	Beta.push_back(cvCloneMat(Phi[0]));
	intBeta.push_back(calcIntegral(Beta[0]));
	Coeff.push_back(B[L[0]]);
	Residue.push_back(SqrNorm(matPatch) - Sqr(Coeff[0]));
	innerPhiF.push_back(innerProd(Phi[0], matPatch));

	K = 1;

	dtime = clock();

	while (K < MaxNumBases && Residue[K - 1] > sigma) {
		// O(N + KWH) for each repetition
		// Selecting the best fit Base in O(N)
		curMaxE = -inf; tempL = -1;
		NBS_MAT_TYPE tempPhiBase, tempE;
		for(int i = 0; i < sampN; ++ i) {
			int ii = sampIndex[i];
			if (!selected[ii]) {
				progress[ii] = K;
				tempPhiBase = listOfBases[ii].innerProd(intPhi[K - 1]);
				B[ii] -= innerPhiF[K - 1] * tempPhiBase;
				D[ii] -= Sqr(tempPhiBase);
				tempE = fabs(B[ii]) < eps ? 0 : (Sqr(B[ii]) / D[ii]);
				if (tempE > curMaxE) {
					curMaxE = tempE;
					tempL = ii;
				}
			}
		}

		if (tempL == -1) break;

		int optIndex = tempL;
		HaarBase & optbase = listOfBases[optIndex];

		int dX = samp_width;
		int dY = samp_height;

		for (int htl = max(0, -dY + optbase.htl); htl < +dY + optbase.htl+1; ++ htl)
			for (int wtl = max(0, -dX + optbase.wtl); wtl < +dX + optbase.wtl+1; ++ wtl)
				for (int h = max(1, -dY + optbase.h + optbase.htl - htl); h < min(imagH - htl, +dY + optbase.h + optbase.htl - htl) + 1; ++ h)
					for (int w = max(1, -dX + optbase.w + optbase.wtl - wtl); w < min(imagW - wtl, +dX + optbase.w + optbase.wtl - wtl) + 1; ++ w) {
						if (htl % dY == 0 && wtl % dX == 0 && h % dY == 0 && w % dX == 0)
							continue;

						int ii = calc_index_haarbase(htl, wtl, h, w, imagH, imagW);

						for (int k = progress[ii]; k < K; ++ k) {
							tempPhiBase = listOfBases[ii].innerProd(intPhi[k]);
							B[ii] -= innerPhiF[k] * tempPhiBase;
							D[ii] -= tempPhiBase * tempPhiBase;
						}

						progress[ii] = K;

						tempE = fabs(B[ii]) < eps ? 0 : (Sqr(B[ii]) / D[ii]);

						if (tempE > curMaxE) {
							curMaxE = tempE;
							tempL = ii;
						}
					}

		L.push_back(tempL);
		selected[tempL] = true;

		// Computing parameters
		Residue.push_back(Residue[K - 1] - curMaxE);        // O(1)
		CvMat * tempMat = calcPhi(L[K], Phi, intPhi);       // O(KWH)
		Phi.push_back(calcMatDiv(tempMat, sqrt(D[L[K]])));  // O(WH)
		Beta.push_back(calcMatDiv(tempMat, D[L[K]]));       // O(WH)
		intBeta.push_back(calcIntegral(Beta[K]));           // O(WH)
		cvReleaseMat(&tempMat);
		intPhi.push_back(calcIntegral(Phi[K]));             // O(WH)
		innerPhiF.push_back(innerProd(Phi[K], matPatch));   // O(WH)
		Coeff.push_back(B[L[K]] / D[L[K]]);                 // O(1)

		// Post-processing at the end of each selection
		NBS_MAT_TYPE temp;
		for(int i = 0; i < K; ++ i) {   // O(KWH)
			temp = listOfBases[L[K]].innerProd(intBeta[i]);
			for(int r = 0; r < Beta[i]->rows; ++ r)
				for(int c = 0; c < Beta[i]->cols; ++ c)
					CV_MAT_ELEM(*Beta[i], NBS_MAT_TYPE, r, c)
					-= CV_MAT_ELEM(*Beta[K], NBS_MAT_TYPE, r, c) * temp;
			// 'temp' is real, so its conjugate is itself;
			Coeff[i] -= temp * Coeff[K];
			cvReleaseMat(&intBeta[i]);
			intBeta[i] = calcIntegral(Beta[i]);
		}

		// Now, get the next binary base
		++ K;
	}

	// Copy Answers
	this->indexOfBases.clear();
	for(int i = 0; i < (int)L.size(); ++ i)
		this->indexOfBases.push_back(listOfBases[L[i]]);
	this->coeffOfBases = Coeff;

	// Release Memories
	delete [] progress;
	cvReleaseMat(&matPatch);
	cvReleaseMat(&intPatch);
	delete [] B;
	delete [] D;
	delete [] selected;
	for(int i = 0; i < (int)Phi.size(); ++ i) {
		cvReleaseMat(&Phi[i]);
		cvReleaseMat(&intPhi[i]);
		cvReleaseMat(&Beta[i]);
		cvReleaseMat(&intBeta[i]);
	}
}