// DnbsClassBeta.cpp
// Implementation of 'DnbsClassBeta' class.
//
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

#include "DnbsClassBeta.h"
#include "LyonLib.h"
#include "headers.h"

DnbsClassBeta::DnbsClassBeta(void) {
}

DnbsClassBeta::~DnbsClassBeta(void) {
}

void DnbsClassBeta::computeNBS(int MaxNumBases,
                               int MaxNumForegrounds,
                               int MaxNumBackgrounds,
                               IplImage * frame,
                               CvRect rect,
                               int margin_threshold_x,
                               int margin_threshold_y,
                               int distance_threshold,
                               double lambda,
                               int method,
							   int samp_width,
							   int samp_height,
                               CvMat ** previous_foregrounds,
                               double coherence) {
    // Compute the original NBS of foreground
    // and then sample the background templates
    // and re-compute the NBS distinguish the
    // foregrounds and the backgrounds
    // Ang Li, Mar.4, 2010.

    // First compute the original NBS representation
    naivenbs.importTarget(target);
    naivenbs.listOfBases = listOfBases;
    if (method == 2 || method == 3)
		naivenbs.computeNBS(30, samp_width, samp_height);
	else naivenbs.computeNBS(30);

	/*
    if (MaxNumForegrounds == 1 && MaxNumBackgrounds == 0) {
        indexOfBases = naivenbs.indexOfBases;
        coefMatrix();
        return ;
    }
	*/

    // Second update the NBS by sampling background templates

    updateNBS(
        MaxNumBases,
        MaxNumForegrounds,
        MaxNumBackgrounds,
        frame,
        rect,
        margin_threshold_x,
        margin_threshold_y,
        distance_threshold,
        lambda,
        method,
		samp_width,
		samp_height,
        previous_foregrounds,
        coherence);
}

void DnbsClassBeta::sample_backgrounds(CvRect rect,
                                       int MaxNumBackgrounds,
                                       int distance_threshold,
                                       int x_min, int x_max,
                                       int y_min, int y_max,
                                       vector< CvRect > & pnt,
                                       vector< pair<double, int> > & arr,
                                       double ** weight,
                                       vector < CvRect > & samples) {
    // Iteratively sample 'MaxNumBackgrounds' backgrounds
    // (1) by avoiding the similarity (SSD/NCC) of any pair
    // less than a given threshold.
    // or (2) by avoiding the distance between any pair
    // of samples' locations less than a give threshold.
    // (Here we choose (2))

    samples.clear();

    // Clear the region around the foreground
    for (int x = -distance_threshold; x < distance_threshold; ++ x)
        for (int y = -distance_threshold; y < distance_threshold; ++ y)
            if (x + rect.x >= x_min && x + rect.x < x_max
                    && y + rect.y >= y_min && y + rect.y < y_max)
                weight[x + rect.x - x_min][y + rect.y - y_min] = -1;

    for (int i = 0; i < (int)arr.size() && (int)samples.size() < MaxNumBackgrounds; ++ i) {
        if (weight[pnt[arr[i].second].x - x_min][pnt[arr[i].second].y - y_min] != -1)
            samples.push_back(pnt[arr[i].second]);
        //else continue;
        int dest_x = pnt[arr[i].second].x;
        int dest_y = pnt[arr[i].second].y;
        for (int x = -distance_threshold; x < distance_threshold; ++ x)
            for (int y = -distance_threshold; y < distance_threshold; ++ y)
                if (x + dest_x >= x_min && x + dest_x < x_max
                        && y + dest_y >= y_min && y + dest_y < y_max)
                    weight[x + dest_x - x_min][y + dest_y - y_min] = -1;
    }
}

void DnbsClassBeta::updateNBS(int MaxNumBases,
                              int MaxNumForegrounds,
                              int MaxNumBackgrounds,
                              IplImage * frame,
                              CvRect rect,
                              int margin_threshold_x,
                              int margin_threshold_y,
                              int distance_threshold,
                              double lambda,
                              int method,
							  int samp_width,
							  int samp_height,
                              CvMat ** previous_foregrounds,
                              double coherence) {
    // Sampling background templates according to
    // the previously constructed NBS of the foreground
    // Ang Li, Mar.4, 2010.

    // Get the window of neighborhood
    int margin_x = margin_threshold_x;//max(margin_threshold_x, rect.width);
    int margin_y = margin_threshold_y;//max(margin_threshold_y, rect.height);

    // Construct the SSD distribution of the neighborhood
    double fnorm = naivenbs.l2norm();
    //double denom_beta = 256. * 256. * rect.width * rect.height;
    CvMat * intframe = lyon::calcIntegral(frame, 0, frame->width, 0, frame->height);
    CvMat * int2frame = lyon::calcSqrIntegral(frame, 0, frame->width, 0, frame->height);

    int x_min = max(rect.x - margin_x, 0);
    int x_max = min(rect.x + margin_x, frame->width - imagW);
    int y_min = max(rect.y - margin_y, 0);
    int y_max = min(rect.y + margin_y, frame->height - imagH);
    CvSize neighborSize = cvSize(x_max - x_min, y_max - y_min);
    double ** weight = (double**)malloc(sizeof(double*) * neighborSize.width);
    for (int i = 0; i < neighborSize.width; ++ i) {
        weight[i] = (double*)malloc(sizeof(double) * neighborSize.height);
        memset(weight[i], 0, sizeof(double) * neighborSize.height);
    }

    vector< pair<double, int> > arr;
    vector< CvRect > pnt;

    // Calculate SSD between neighbor patches and the target
    // approximated by previous constructed NBS
    // weight[x - x_min][y - y_min] = SSD(x,y)
    for (int x = x_min; x < x_max; ++ x) {
        for (int y = y_min; y < y_max; ++ y) {
            double gnorm
            = lyon::calcIntSum(int2frame, x, x + imagW, y, y + imagH);
            double iproduct
            = naivenbs.innerProduct(intframe, x, x + imagW, y, y + imagH);
            double weight_of_patch = (gnorm - 2. * iproduct + fnorm);
            weight[x - x_min][y - y_min] = weight_of_patch;
            pnt.push_back(cvRect(x, y, target->width, target->height));
            arr.push_back(make_pair(weight_of_patch, pnt.size() - 1));
        }
	}

    // Sample background templates in light of the SSD distribution
    // Sort 'arr' in ascending order of weights
    sort(arr.begin(), arr.end());

    // Set sample vectors and their weights
    vector < CvRect > samples;
    sample_backgrounds(
        rect, MaxNumBackgrounds, distance_threshold,
        x_min, x_max, y_min, y_max, pnt, arr, weight, samples);

    printf ("Number of sampled backgrounds = %d\n", samples.size());

    // bgsample records the sampled backgrounds
    bgsample.clear();
    for (int i = 0; i < (int)samples.size(); ++ i)
        bgsample.push_back(samples[i]);

    // Collect positive and negative samples
    vector<CvMat*> foreground;
    vector<CvMat*> background;
    foreground.push_back(target);
    for (int i = 0; i < MaxNumForegrounds - 1; ++ i)
        if (previous_foregrounds[i] != NULL)
            foreground.push_back(previous_foregrounds[i]);
    for (int i = 0; i < (int)samples.size(); ++ i) {
        IplImage * patch
        = lyon::selectedPatch(frame, samples[i]);
        background.push_back(image2mat(patch));
        cvReleaseImage(&patch);
    }

    // Re-compute the NBS of foreground that discriminates the backgrounds
    switch (method) {
        case 0:
            computeDNBS_v1(MaxNumBases, lambda, foreground, background);
            break;
        case 1:
            computeDNBS_v2(MaxNumBases, frame, lambda, samples);
            break;
        case 2:
            computeDNBS_v1_pruning(
				MaxNumBases, lambda, foreground, background, samp_width, samp_height);
            break;
		case 3:
			computeDNBS_v1_downsample(
				MaxNumBases, lambda, foreground, background, samp_width, samp_height);
        default:
            break;
    }

	// Back-up
	for (int i = 0; i < (int)fgrdSample.size(); ++ i)
		cvReleaseMat(&fgrdSample[i]);
	for (int i = 0; i < (int)bgrdSample.size(); ++ i)
		cvReleaseMat(&bgrdSample[i]);
	fgrdSample.clear();
	bgrdSample.clear();
	for (int i = 0; i < (int)foreground.size(); ++ i)
		fgrdSample.push_back(cvCloneMat(foreground[i]));
	for (int i = 0; i < (int)background.size(); ++ i)
		bgrdSample.push_back(cvCloneMat(background[i]));

    // Release arrays and images
    for (int i = 0; i < (int)background.size(); ++ i)
        cvReleaseMat(&background[i]);
    for (int i = 0; i < neighborSize.width; ++ i)
        free(weight[i]);
    free(weight);
    cvReleaseMat(&intframe);
    cvReleaseMat(&int2frame);
}

void DnbsClassBeta::computeDNBS_v1(
    int MaxNumBases,
    double lambda,
    vector<CvMat*> foreground,
    vector<CvMat*> background) {
    // Discriminative Non-orthogonal Binary Subspace
    // Formulation.
    // \arg\min_Phi{1/Nf\sum||fi-R_Phi{fi}||^2
    //              -lambda/Nb\sum_i||bi-R_Phi{bi}||^2}
    // where, lambda is a trade-off.
    // Ang Li, Nanjing University, Oct.7 2010.

    // Supposing frame is gray-level image with 3 channels

    // Record the width and height
    int width = foreground[0]->width;
    int height = foreground[0]->height;

    int num_foregrounds = foreground.size();
    int num_backgrounds = background.size();

    // Definitions
    // K is the number of selected bases
    int K;
    // N is the total number of candidate bases
    int N = listOfBases.size();

    // B and D are two arrays of parameters in deciding the
    // optimized index of the next to-be-selected binary base.
    // i.e. selecting the one that maximising E_n(=|BF_n||BG_n|/D_n)
    vector<NBS_MAT_TYPE*> BF;
    for (int i = 0; i < num_foregrounds; ++ i) {
        NBS_MAT_TYPE * tempBF
        = new NBS_MAT_TYPE[N];
        BF.push_back(tempBF);
    }

    vector<NBS_MAT_TYPE*> BG;
    for (int i = 0; i < num_backgrounds; ++ i) {
        NBS_MAT_TYPE * tempBG
        = new NBS_MAT_TYPE[N];
        BG.push_back(tempBG);
    }
    NBS_MAT_TYPE * D = new NBS_MAT_TYPE[N];

    // Indicating whether the base is selected to form the subspace
    bool * selected = new bool[N];

    // Index of each selected binary bases
    vector<int> L;
    // Projection Coefficients
    vector<NBS_MAT_TYPE> Coeff;

    // Phi's and their integral image
    vector<CvMat*> Phi;
    vector<CvMat*> intPhi;

    // Inner product of each Phi[K] and image F(i.e. matPatch)
    vector< vector<NBS_MAT_TYPE> > innerPhiF;
    for (int i = 0; i < num_foregrounds; ++ i)
        innerPhiF.push_back(vector<NBS_MAT_TYPE>());
    // Inner product of each Phi[K] and image G
    vector< vector<NBS_MAT_TYPE> > innerPhiG;
    for (int i = 0; i < num_backgrounds; ++ i)
        innerPhiG.push_back(vector<NBS_MAT_TYPE>());

    // Beta's and their integral image
    vector<CvMat*> Beta;
    vector<CvMat*> intBeta;

    // Initialization
    // Reset memories
    memset(selected, 0, sizeof(bool) * N);
    // Computing the integral image of the original image
    vector<CvMat*> integralF;
    for (int i = 0; i < num_foregrounds; ++ i)
        integralF.push_back(calcIntegral(foreground[i]));
    vector<CvMat*> integralG;
    for (int i = 0; i < num_backgrounds; ++ i)
        integralG.push_back(calcIntegral(background[i]));

    // Process of initializing parameters and get the first choice
    int tempL = 0;
    double maxFunc = -inf;
    double curFunc;
    for(int i = 0; i < N; ++ i) {
        for (int j = 0; j < num_foregrounds; ++ j)
            BF[j][i] = listOfBases[i].innerProd(integralF[j]);
        for (int j = 0; j < num_backgrounds; ++ j)
            BG[j][i] = listOfBases[i].innerProd(integralG[j]);
        D[i] = 1;

        curFunc = 0;
        for (int j = 0; j < num_foregrounds; ++ j)
            curFunc += BF[j][i] * BF[j][i] / num_foregrounds;
        for (int j = 0; j < num_backgrounds; ++ j)
            curFunc -= BG[j][i] * BG[j][i] / num_backgrounds * lambda;

        if (curFunc > maxFunc) {
            tempL = i;
            maxFunc = curFunc;
        }
    }
    L.push_back(tempL);
    selected[tempL] = true;

    // Computing parameters
    Phi.push_back(listOfBases[L[0]].toMatrix(height, width));
    intPhi.push_back(calcIntegral(Phi[0]));
    Beta.push_back(cvCloneMat(Phi[0]));
    intBeta.push_back(calcIntegral(Beta[0]));
    //Coeff.push_back(BF[L[0]]);
    for (int i = 0; i < num_foregrounds; ++ i)
        innerPhiF[i].push_back(innerProd(Phi[0], foreground[i]));
    for (int i = 0; i < num_backgrounds; ++ i)
        innerPhiG[i].push_back(innerProd(Phi[0], background[i]));

    //FILE * fout = fopen("res.txt", "a+");

    int curtime = clock();

    K = 1;
    // Total Time Complexity is O(KN + K^2WH) = O(KWH(WH+K))
    while (K < MaxNumBases) {
        // Selecting the best fit Base
        maxFunc = -inf;
        tempL = -1;
        double tempPhiBase;
        for(int i = 0; i < N; ++ i) {
            //if (!selected[i]) {
            tempPhiBase = listOfBases[i].innerProd(intPhi[K - 1]);
            for (int j = 0; j < num_foregrounds; ++ j)
                BF[j][i] -= innerPhiF[j][K - 1] * tempPhiBase;
            for (int j = 0; j < num_backgrounds; ++ j)
                BG[j][i] -= innerPhiG[j][K - 1] * tempPhiBase;
            D[i] -= tempPhiBase * tempPhiBase;

            curFunc = 0;
            for (int j = 0; j < num_foregrounds; ++ j)
                curFunc += BF[j][i] * BF[j][i] / num_foregrounds;
            for (int j = 0; j < num_backgrounds; ++ j)
                curFunc -= BG[j][i] * BG[j][i] / num_backgrounds * lambda;

            curFunc = fabs(curFunc) < sqeps ? 0 : (curFunc / D[i]);

            if (curFunc > maxFunc) {
                maxFunc = curFunc;
                tempL = i;
            }
        }

        if (tempL == -1) break;
        L.push_back(tempL);
        selected[tempL] = true;

        // Computing parameters
        CvMat * tempMat = calcPhi(L[K], Phi, intPhi);       // O(KWH)
        Phi.push_back(calcMatDiv(tempMat, sqrt(D[L[K]])));  // O(WH)
        Beta.push_back(calcMatDiv(tempMat, D[L[K]]));       // O(WH)
        intBeta.push_back(calcIntegral(Beta[K]));           // O(WH)
        cvReleaseMat(&tempMat);
        intPhi.push_back(calcIntegral(Phi[K]));             // O(WH)
        for (int i = 0; i < num_foregrounds; ++ i)
            innerPhiF[i].push_back(innerProd(Phi[K], foreground[i]));
        for (int i = 0; i < num_backgrounds; ++ i)
            innerPhiG[i].push_back(innerProd(Phi[K], background[i]));

        ++ K;
    }

	lyon::settime((clock() - curtime) / (double)CLOCKS_PER_SEC);
    // Copy Answers
    this->indexOfBases.clear();
    for(int i = 0; i < (int)L.size(); ++ i)
        this->indexOfBases.push_back(listOfBases[L[i]]);
    CvMat * intTarget = calcIntegral(target);
    this->coefMatrix();
    this->coeffOfBases = compCoeff(intTarget);


    //demo();

    // Release Memories
    delete [] D;
    delete [] selected;
    for (int i = 0; i < num_foregrounds; ++ i) {
        delete [] BF[i];
    }
    for (int i = 0; i < num_backgrounds; ++ i) {
        delete [] BG[i];
    }
    for(int i = 0; i < (int)Phi.size(); ++ i) {
        cvReleaseMat(&Phi[i]);
        cvReleaseMat(&intPhi[i]);
        cvReleaseMat(&Beta[i]);
        cvReleaseMat(&intBeta[i]);
    }
}



void DnbsClassBeta::computeDNBS_v2(
    int MaxNumBases,
    IplImage * frame,
    double lambda,
    vector< CvRect > & samples) {
    // Discriminative Nonorthogonal Binary Subspace
    //
    // Description. The expected subspace is specially designed to
    // discriminate the possible foreground objects and background
    // templates with nonorthogonal Haar-like box bases.
    //
    // Formulation. The objective function of D-NBS is as follows,
    // \arg\min_Phi {1/N_f||F-R_Phi{F}||_F^2 - lambda/N_b||B-R_Phi{B}||_F^2}
    // where, ||*||_F denotes the Frobenius norm. Phi is the anticipated D-NBS.
    // F is the set of foreground patches and B is the set of background templates.
    // N_f is the size of set F (i.e. the number of selected foregrounds) and
    // N_b is the size of background set B, lambda is the trade-off.
    //
    // Notification. This is the ever fastest edition. :)
    // And it is still in developing.
    //
    // Ang Li
    // Department of Computer Science and Technology
    // Nanjing University
    // mailto:nju.angli@gmail.com
    // May 13, 2010.

    // frame is assumed gray-level image with 3 channels

    // Recording the width and height
    int width = target->width;
    int height = target->height;

    // Recording the foreground matrix
    CvMat * fmat = cvCloneMat(target);

    // Recording the background matrices
    int num_backgrounds = samples.size();
    vector<CvMat*> bmat;
    for (int i = 0; i < num_backgrounds; ++ i) {
        CvMat * tempmat
        = cvCreateMat(height, width, NBS_CV_MAT_TYPE);
        for (int r = 0; r < height; ++ r)
            for (int c = 0; c < width; ++ c)
                CV_MAT_ELEM(*tempmat, NBS_MAT_TYPE, r, c)
                = CV_IMAGE_ELEM(
                      frame, uchar, samples[i].y + r,
                      (samples[i].x + c) * 3);
        bmat.push_back(tempmat);
    }

    // Definitions
    // K is the number of selected bases
    int K;
    // N is the total number of candidate bases
    int N = listOfBases.size();
    // Lvalue is an array which is dedicated to the selection of bases.
    double * Lvalue = new double[N];
    // D[i] = ||phi_i - R_Phi(phi_i)||^2, appearing in the denominators.
    double * D = new double[N];

    // Indicating whether the base is selected to form the subspace
    bool * selected = new bool[N];
    memset(selected, 0, sizeof(bool) * N);

    // Index of each selected binary base
    vector<int> L;

    // Phi's and their integral image
    vector<CvMat*> Phi;
    vector<CvMat*> intPhi;

    // Computing the integral image of the original image
    CvMat * integralF
    = calcIntegral(fmat);
    vector<CvMat*> integralB;
    for (int i = 0; i < num_backgrounds; ++ i)
        integralB.push_back(calcIntegral(bmat[i]));
    CvMat * epsilon_f
    = cvCloneMat(fmat);
    vector<CvMat*> epsilon_b;
    for (int i = 0; i < num_backgrounds; ++ i)
        epsilon_b.push_back(cvCloneMat(bmat[i]));
    CvMat * int_epsilon_f
    = cvCloneMat(integralF);
    vector<CvMat*> int_epsilon_b;
    for (int i = 0; i < num_backgrounds; ++ i)
        int_epsilon_b.push_back(cvCloneMat(integralB[i]));

    // Process of initializing parameters and get the first choice
    // Complexity: O(N*N_s)
    int tempL = 0;
    double maxFunc = -inf;
    double curFunc;
    for(int i = 0; i < N; ++ i) {
        curFunc = 0;
        for (int j = 0; j < num_backgrounds; ++ j)
            curFunc += lyon::sqr(listOfBases[i].innerProd(integralB[j]));
        Lvalue[i] = lyon::sqr(listOfBases[i].innerProd(integralF))
                    - lambda * curFunc / num_backgrounds;
        if (Lvalue[i] > maxFunc) {
            tempL = i;
            maxFunc = Lvalue[i];
        }
        D[i] = 1;
    }
    L.push_back(tempL);
    selected[tempL] = true;

    double alpha_f;
    double * alpha_b = new double[num_backgrounds];
    double Sk;
    CvMat * Ik = cvCreateMat(height, width, NBS_CV_MAT_TYPE);

    // Computing parameters
    K = 1;
    Phi.push_back(listOfBases[L[0]].toMatrix(imagH, imagW));
    intPhi.push_back(calcIntegral(Phi[0]));

    //FILE * fout = fopen("res.txt", "a+");

    int curtime = clock();

    // Theoretical Time Complexity is O(KN + KWHN_s + K^2WH) = O(KWH(WH+N_s+K))
    while (K < MaxNumBases) {
        // Notations:
        //  * base_{k-1} = listOfBases[L[K-1]]

        // Compute alpha
        // alpha(x) = <phi_{k-1}, x> = <base_{k-1}, epsilon_{k-2}(x)>
        alpha_f = listOfBases[L[K-1]].innerProd(int_epsilon_f);
        for (int i = 0; i < num_backgrounds; ++ i)
            alpha_b[i] = listOfBases[L[K-1]].innerProd(int_epsilon_b[i]);

        // u[k-1] = D[L[k-1]] = ||phi_{k-1}||^2
        double u_ksub1 = D[L[K-1]];
        double sqrt_u_ksub1 = sqrt(u_ksub1);

        // Compute S_k
        // S_k = 1/u[k-1](1/N_f\sum{alpha^2(f)}-lambda/N_b\sum{alpha^2(b)})
        Sk = alpha_f * alpha_f / u_ksub1;
        for (int i = 0; i < num_backgrounds; ++ i)
            Sk -= lambda / num_backgrounds * alpha_b[i] * alpha_b[i] / u_ksub1;

        // Compute I_k
        // ita_k(x) = alpha(x)epsilon_{k-2}(x)
        // I_k = -2/sqrt(u[k-1])(1/N_f\sum{ita_k(f)}-lambda/N_b\sum{ita_k(b)})
        for (int r = 0; r < Ik->rows; ++ r)
            for (int c = 0; c < Ik->cols; ++ c) {
                double curr
                = alpha_f * CV_MAT_ELEM(*epsilon_f, double, r, c);
                for (int i = 0; i < num_backgrounds; ++ i)
                    curr -= lambda / num_backgrounds * alpha_b[i]
                            * CV_MAT_ELEM(*epsilon_b[i], double, r, c);
                CV_MAT_ELEM(*Ik, double, r, c) = -2 * curr / sqrt_u_ksub1;
            }
        // The integral image of I_k
        CvMat * intIk = calcIntegral(Ik);

        // Update epsilon_{k-2}(samples)
        // epsilon_k(x) = epsilon_{k-1}(x) - phi_k<phi_k,x>/||phi_k||^2
        for (int r = 0; r < epsilon_f->rows; ++ r)
            for (int c = 0; c < epsilon_f->cols; ++ c)
                CV_MAT_ELEM(*epsilon_f, double, r, c)
                -= CV_MAT_ELEM(*Phi[K-1], double, r, c) * alpha_f / sqrt_u_ksub1;
        for (int i = 0; i < num_backgrounds; ++ i)
            for (int r = 0; r < epsilon_f->rows; ++ r)
                for (int c = 0; c < epsilon_f->cols; ++ c)
                    CV_MAT_ELEM(*epsilon_b[i], double, r, c)
                    -= CV_MAT_ELEM(*Phi[K-1], double, r, c) * alpha_b[i] / sqrt_u_ksub1;
        // Calculate the integral images of epsilon_f's and epsilon_b's.
        cvReleaseMat(&int_epsilon_f);
        int_epsilon_f = calcIntegral(epsilon_f);
        for (int i = 0; i < num_backgrounds; ++ i) {
            cvReleaseMat(&int_epsilon_b[i]);
            int_epsilon_b[i] = calcIntegral(epsilon_b[i]);
        }

        // Selecting the best binary base
        maxFunc = -inf;
        tempL = -1;
        double inner_basei_phi;
        double inner_basei_Ik;
        double temp;
        for(int i = 0; i < N; ++ i)
            if (!selected[i]) {
                inner_basei_phi = listOfBases[i].innerProd(intPhi[K - 1]);
                inner_basei_Ik = listOfBases[i].innerProd(intIk);
                Lvalue[i] += (inner_basei_Ik + inner_basei_phi * Sk) * inner_basei_phi;

                // Update D[i] = ||phi_i - R_Phi(phi_i)||^2
                D[i] -= inner_basei_phi * inner_basei_phi;

                temp = fabs(D[i]) < 1e-10 ? 0 : (Lvalue[i] / D[i]);

                if (temp > maxFunc) {
                    maxFunc = temp;
                    tempL = i;
                }
            }

        if (tempL == -1) break;
        L.push_back(tempL);
        selected[tempL] = true;
        //fprintf (fout, "%d, tempL = %d, Lvalue = %.20e, D = %.20e\n", K, tempL, Lvalue[tempL], D[tempL]);

        // Computing parameters
        CvMat * tempMat = calcPhi(L[K], Phi, intPhi);  // O(KWH)
        Phi.push_back(calcMatDiv(tempMat, sqrt(D[L[K]])));  // O(WH)
        intPhi.push_back(calcIntegral(Phi[K]));        // O(WH)
        cvReleaseMat(&tempMat);

        ++ K;
    }

	lyon::settime((clock() - curtime) / (double)CLOCKS_PER_SEC);
    printf("cost time: %.6lf\n", (clock() - curtime) / (double)CLOCKS_PER_SEC);
    //fclose(fout);

    // Copy Answers
    this->indexOfBases.clear();
    for(int i = 0; i < (int)L.size(); ++ i)
        this->indexOfBases.push_back(listOfBases[L[i]]);
    computeCoeffMatrix();
    this->coeffOfBases = compCoeff(integralF);

    //FILE * fres = fopen("error.txt", "a+");
    //fprintf(fres, "%d\t%.3lf\t%.3lf\n", num_backgrounds, compErrorReconstruction(), (clock() - curtime) / (double)CLOCKS_PER_SEC);
    //fclose(fres);
    //demo();

    // Release Memories
    cvReleaseMat(&fmat);
    cvReleaseMat(&integralF);
    cvReleaseMat(&epsilon_f);
    cvReleaseMat(&int_epsilon_f);
    for (int i = 0; i < num_backgrounds; ++ i) {
        cvReleaseMat(&bmat[i]);
        cvReleaseMat(&integralB[i]);
        cvReleaseMat(&epsilon_b[i]);
        cvReleaseMat(&int_epsilon_b[i]);
    }
    delete [] D;
    delete [] selected;
    for(int i = 0; i < (int)Phi.size(); ++ i) {
        cvReleaseMat(&Phi[i]);
        cvReleaseMat(&intPhi[i]);
    }

    delete [] Lvalue;
    delete [] alpha_b;
    cvReleaseMat(&Ik);
}

void DnbsClassBeta::computeCoeffMatrix(int K) {
    // Computing the coefficient matrix of the D-NBS

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
        = cvCreateMat(target->rows, target->cols, NBS_CV_MAT_TYPE);
        cvZero(temp);
        for(int j = 0; j < K; ++ j) {
            double scale
            = CV_MAT_ELEM(*invPhiTPhi, NBS_MAT_TYPE, j, i);
            for(int r = 0; r < temp->rows; ++ r)
                for(int c = 0; c < temp->cols; ++ c)
                    CV_MAT_ELEM(*temp, NBS_MAT_TYPE, r, c)
                    += indexOfBases[j].element(r, c, temp->cols) * scale;
        }
        this->coeffMat.push_back(temp);
    }
}


double DnbsClassBeta::reconerr(void) {
    // Calculate the reconstruction error of the foreground template
    // through the discriminative nonorthogonal binary subspace.

    CvMat * intTarget
        = lyon::calcIntegral(target);

    double ssd = l2norm();
    for (int r = 0; r < target->rows; ++ r)
        for (int c = 0; c < target->cols; ++ c)
            ssd += lyon::sqr(CV_MAT_ELEM(*target, double, r, c));
    for(int i = 0; i < (int)coeffOfBases.size(); ++ i)
        ssd -= 2. * coeffOfBases[i]
               * indexOfBases[i].innerProd(
                   intTarget, 0, target->width, 0, target->height);
    cvReleaseMat(&intTarget);
    return sqrt(ssd / target->width / target->height);
}



void DnbsClassBeta::computeNBS_discriminant(
    int MaxNumBases,
    IplImage * frame,
    CvRect rect,
    double lambda,
    vector< CvRect > & samples) {
    // Computation of Non-orthogonal Binary Subspace
    // that discriminate the foregrounds and backgrounds
    // by optimizing the equation
    // \arg\min_Phi{(1-lambda)||f-R_Phi{f}||^2
    //              -lambda\sum_i{w_i||bi-R_Phi{bi}||^2}}
    // where, f is the foreground and bi is the i-th
    // background template (top-left is samples[i].second)
    // w_i is the weight for i-th background, i.e.
    // w_i = samples[i].first. lambda is a trade-off.
    // Ang Li, Mar.4, 2010.

    // Supposing frame is gray-level image with 3 channels

    // Record the width and height
    int width = target->width;
    int height = target->height;

    // Record the foreground matrix
    CvMat * fmat = cvCloneMat(target);

    // Record the background matrices
    int num_backgrounds = samples.size();
    vector<CvMat*> gmat;
    for (int i = 0; i < num_backgrounds; ++ i) {
        CvMat * tempmat
        = cvCreateMat(height, width, NBS_CV_MAT_TYPE);
        for (int r = 0; r < height; ++ r)
            for (int c = 0; c < width; ++ c)
                CV_MAT_ELEM(*tempmat, NBS_MAT_TYPE, r, c)
                = CV_IMAGE_ELEM(
                      frame, uchar, samples[i].y + r,
                      (samples[i].x + c) * 3);
        gmat.push_back(tempmat);
    }

    // Definitions
    // K is the number of selected bases
    int K;
    // N is the total number of candidate bases
    int N = listOfBases.size();

    // B and D are two arrays of parameters in deciding the
    // optimized index of the next to-be-selected binary base.
    // i.e. selecting the one that maximising E_n(=|BF_n||BG_n|/D_n)
    NBS_MAT_TYPE * BF = new NBS_MAT_TYPE[N];
    vector<NBS_MAT_TYPE*> BG;
    for (int i = 0; i < num_backgrounds; ++ i) {
        NBS_MAT_TYPE * tempBG
        = new NBS_MAT_TYPE[N];
        BG.push_back(tempBG);
    }
    NBS_MAT_TYPE * D = new NBS_MAT_TYPE[N];

    // Indicating whether the base is selected to form the subspace
    bool * selected = new bool[N];

    // Index of each selected binary bases
    vector<int> L;
    // Projection Coefficients
    vector<NBS_MAT_TYPE> Coeff;

    // Phi's and their integral image
    vector<CvMat*> Phi;
    vector<CvMat*> intPhi;

    // Inner product of each Phi[K] and image F(i.e. matPatch)
    vector<NBS_MAT_TYPE> innerPhiF;
    // Inner product of each Phi[K] and image G
    vector< vector<NBS_MAT_TYPE> > innerPhiG;
    for (int i = 0; i < num_backgrounds; ++ i)
        innerPhiG.push_back(vector<NBS_MAT_TYPE>());

    // Beta's and their integral image
    vector<CvMat*> Beta;
    vector<CvMat*> intBeta;

    // Initialization
    // Reset memories
    memset(selected, 0, sizeof(bool) * N);
    // Computing the integral image of the original image
    CvMat * integralF
    = calcIntegral(fmat);
    vector<CvMat*> integralG;
    for (int i = 0; i < num_backgrounds; ++ i)
        integralG.push_back(calcIntegral(gmat[i]));

    // Process of initializing parameters and get the first choice
    int tempL = 0;
    double maxFunc = -inf;
    double curFunc;
    for(int i = 0; i < N; ++ i) {
        BF[i] = listOfBases[i].innerProd(integralF);
        for (int j = 0; j < num_backgrounds; ++ j)
            BG[j][i] = listOfBases[i].innerProd(integralG[j]);
        D[i] = 1;
        curFunc = 0;
        for (int j = 0; j < num_backgrounds; ++ j)
            curFunc += BG[j][i] * BG[j][i];
        curFunc = - curFunc * lambda / num_backgrounds + BF[i] * BF[i];
        if (curFunc > maxFunc) {
            tempL = i;
            maxFunc = curFunc;
        }
    }
    L.push_back(tempL);
    selected[tempL] = true;

    // Computing parameters
    Phi.push_back(listOfBases[L[0]].toMatrix(imagH, imagW));
    intPhi.push_back(calcIntegral(Phi[0]));
    Beta.push_back(cvCloneMat(Phi[0]));
    intBeta.push_back(calcIntegral(Beta[0]));
    Coeff.push_back(BF[L[0]]);
    innerPhiF.push_back(innerProd(Phi[0], fmat));
    for (int i = 0; i < num_backgrounds; ++ i)
        innerPhiG[i].push_back(innerProd(Phi[0], gmat[i]));

    //FILE * fout = fopen("res.txt", "a+");

    /*    int curtime = clock();*/

    K = 1;
    // Total Time Complexity is O(KN + K^2WH) = O(KWH(WH+K))
    while (K < MaxNumBases) {
        // Selecting the best fit Base
        maxFunc = -inf;
        tempL = -1;
        double tempPhiBase;
        //double lvalue;
        //double di;
        for(int i = 0; i < N; ++ i) {
            //if (!selected[i]) {
            tempPhiBase = listOfBases[i].innerProd(intPhi[K - 1]);
            BF[i] -= innerPhiF[K - 1] * tempPhiBase;
            for (int j = 0; j < num_backgrounds; ++ j)
                BG[j][i] -= innerPhiG[j][K - 1] * tempPhiBase;
            D[i] -= tempPhiBase * tempPhiBase;
            curFunc = 0;
            for (int j = 0; j < num_backgrounds; ++ j)
                curFunc += BG[j][i] * BG[j][i];
            curFunc = - curFunc / num_backgrounds * lambda + BF[i] * BF[i];

            curFunc = fabs(curFunc) < sqeps ? 0 : (curFunc / D[i]);

            if (curFunc > maxFunc) {
                //lvalue = curFunc * D[i];
                //di = D[i];
                maxFunc = curFunc;
                tempL = i;
            }
        }

        if (tempL == -1) break;
        L.push_back(tempL);
        selected[tempL] = true;

        //dprintf ("%d, Lvalue = %.10e, D = %.10e\n", K, lvalue, di);

        // Computing parameters
        CvMat * tempMat = calcPhi(L[K], Phi, intPhi);       // O(KWH)
        Phi.push_back(calcMatDiv(tempMat, sqrt(D[L[K]])));  // O(WH)
        Beta.push_back(calcMatDiv(tempMat, D[L[K]]));       // O(WH)
        intBeta.push_back(calcIntegral(Beta[K]));           // O(WH)
        cvReleaseMat(&tempMat);
        intPhi.push_back(calcIntegral(Phi[K]));             // O(WH)
        innerPhiF.push_back(innerProd(Phi[K], fmat));       // O(WH)
        for (int i = 0; i < num_backgrounds; ++ i)
            innerPhiG[i].push_back(innerProd(Phi[K], gmat[i]));
        Coeff.push_back(BF[L[K]] / D[L[K]]);                // O(1)

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

        ++ K;
    }

    // Copy Answers
    this->indexOfBases.clear();
    for(int i = 0; i < (int)L.size(); ++ i)
        this->indexOfBases.push_back(listOfBases[L[i]]);
    this->coeffOfBases = Coeff;

    //demo();

    // Release Memories
    cvReleaseMat(&fmat);
    delete [] BF;
    delete [] D;
    delete [] selected;
    for (int i = 0; i < num_backgrounds; ++ i) {
        cvReleaseMat(&gmat[i]);
        delete [] BG[i];
    }
    for(int i = 0; i < (int)Phi.size(); ++ i) {
        cvReleaseMat(&Phi[i]);
        cvReleaseMat(&intPhi[i]);
        cvReleaseMat(&Beta[i]);
        cvReleaseMat(&intBeta[i]);
    }
}


void DnbsClassBeta::computeNBS_discriminant(
    int MaxNumBases,
    IplImage * frame,
    CvRect rect,
    double lambda,
    vector< CvRect > & samples,
    double coherence) {
    // Computation of Non-orthogonal Binary Subspace
    // that discriminate the foregrounds and backgrounds
    // by optimizing the equation
    // \arg\min_Phi{(1-lambda)||f-R_Phi{f}||^2
    //              -lambda\sum_i{w_i||bi-R_Phi{bi}||^2}}
    // where, f is the foreground and bi is the i-th
    // background template (top-left is samples[i].second)
    // w_i is the weight for i-th background, i.e.
    // w_i = samples[i].first. lambda is a trade-off.
    // Ang Li, Mar.4, 2010.

    // Supposing frame is gray-level image with 3 channels

    // Record the width and height
    int width = target->width;
    int height = target->height;

    // Record the foreground matrix
    CvMat * fmat = cvCloneMat(target);

    // Record the background matrices
    int num_backgrounds = samples.size();
    vector<CvMat*> gmat;
    for (int i = 0; i < num_backgrounds; ++ i) {
        CvMat * tempmat
        = cvCreateMat(height, width, NBS_CV_MAT_TYPE);
        for (int r = 0; r < height; ++ r)
            for (int c = 0; c < width; ++ c)
                CV_MAT_ELEM(*tempmat, NBS_MAT_TYPE, r, c)
                = CV_IMAGE_ELEM(
                      frame, uchar, samples[i].y + r,
                      (samples[i].x + c) * 3);
        gmat.push_back(tempmat);
    }

    // Definitions
    // K is the number of selected bases
    int K;
    // N is the total number of candidate bases
    int N = listOfBases.size();

    // B and D are two arrays of parameters in deciding the
    // optimized index of the next to-be-selected binary base.
    // i.e. selecting the one that maximising E_n(=|BF_n||BG_n|/D_n)
    NBS_MAT_TYPE * BF = new NBS_MAT_TYPE[N];
    vector<NBS_MAT_TYPE*> BG;
    for (int i = 0; i < num_backgrounds; ++ i) {
        NBS_MAT_TYPE * tempBG
        = new NBS_MAT_TYPE[N];
        BG.push_back(tempBG);
    }
    NBS_MAT_TYPE * D = new NBS_MAT_TYPE[N];

    // Indicating whether the base is selected to form the subspace
    bool * selected = new bool[N];

    // Index of each selected binary bases
    vector<int> L;
    // Projection Coefficients
    vector<NBS_MAT_TYPE> Coeff;

    // Phi's and their integral image
    vector<CvMat*> Phi;
    vector<CvMat*> intPhi;

    // Inner product of each Phi[K] and image F(i.e. matPatch)
    vector<NBS_MAT_TYPE> innerPhiF;
    // Inner product of each Phi[K] and image G
    vector< vector<NBS_MAT_TYPE> > innerPhiG;
    for (int i = 0; i < num_backgrounds; ++ i)
        innerPhiG.push_back(vector<NBS_MAT_TYPE>());

    // Beta's and their integral image
    vector<CvMat*> Beta;
    vector<CvMat*> intBeta;

    // Initialization
    // Reset memories
    memset(selected, 0, sizeof(bool) * N);
    // Computing the integral image of the original image
    CvMat * integralF
    = calcIntegral(fmat);
    vector<CvMat*> integralG;
    for (int i = 0; i < num_backgrounds; ++ i)
        integralG.push_back(calcIntegral(gmat[i]));

    // Process of initializing parameters and get the first choice
    int tempL = 0;
    double maxFunc = -inf;
    double curFunc;
    for(int i = 0; i < N; ++ i) {
        BF[i] = listOfBases[i].innerProd(integralF);
        for (int j = 0; j < num_backgrounds; ++ j)
            BG[j][i] = listOfBases[i].innerProd(integralG[j]);
        D[i] = 1;
        curFunc = 0;
        for (int j = 0; j < num_backgrounds; ++ j)
            curFunc += BG[j][i] * BG[j][i];
        curFunc = - curFunc * lambda + BF[i] * BF[i];
        if (curFunc > maxFunc) {
            tempL = i;
            maxFunc = curFunc;
        }
    }
    L.push_back(tempL);
    selected[tempL] = true;

    // Computing parameters
    Phi.push_back(listOfBases[L[0]].toMatrix(imagH, imagW));
    intPhi.push_back(calcIntegral(Phi[0]));
    Beta.push_back(cvCloneMat(Phi[0]));
    intBeta.push_back(calcIntegral(Beta[0]));
    Coeff.push_back(BF[L[0]]);
    innerPhiF.push_back(innerProd(Phi[0], fmat));
    for (int i = 0; i < num_backgrounds; ++ i)
        innerPhiG[i].push_back(innerProd(Phi[0], gmat[i]));

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

    K = 1;
    // Total Time Complexity is O(KN + K^2WH) = O(KWH(WH+K))
    while (K < MaxNumBases) {
        // Selecting the best fit Base
        maxFunc = -inf;
        tempL = -1;
        double tempPhiBase;
        for(int i = next[N]; i != -1; i = next[i]) {
            tempPhiBase = listOfBases[i].innerProd(intPhi[K - 1]);
            BF[i] -= innerPhiF[K - 1] * tempPhiBase;
            for (int j = 0; j < num_backgrounds; ++ j)
                BG[j][i] -= innerPhiG[j][K - 1] * tempPhiBase;
            D[i] -= tempPhiBase * tempPhiBase;
            curFunc = 0;
            for (int j = 0; j < num_backgrounds; ++ j)
                curFunc += BG[j][i] * BG[j][i];
            curFunc = - curFunc / num_backgrounds * lambda + BF[i] * BF[i];
            curFunc = fabs(curFunc) < sqeps ? 0 : (curFunc / D[i]);

            if (curFunc > maxFunc) {
                maxFunc = curFunc;
                tempL = i;
            }
        }

        if (tempL == -1) break;
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
        CvMat * tempMat = calcPhi(L[K], Phi, intPhi);       // O(KWH)
        Phi.push_back(calcMatDiv(tempMat, sqrt(D[L[K]])));  // O(WH)
        Beta.push_back(calcMatDiv(tempMat, D[L[K]]));       // O(WH)
        intBeta.push_back(calcIntegral(Beta[K]));           // O(WH)
        cvReleaseMat(&tempMat);
        intPhi.push_back(calcIntegral(Phi[K]));             // O(WH)
        innerPhiF.push_back(innerProd(Phi[K], fmat));       // O(WH)
        for (int i = 0; i < num_backgrounds; ++ i)
            innerPhiG[i].push_back(innerProd(Phi[K], gmat[i]));
        Coeff.push_back(BF[L[K]] / D[L[K]]);                // O(1)

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

        ++ K;
    }

    // Copy Answers
    this->indexOfBases.clear();
    for(int i = 0; i < (int)L.size(); ++ i)
        this->indexOfBases.push_back(listOfBases[L[i]]);
    this->coeffOfBases = Coeff;

    //demo();

    // Release Memories
    cvReleaseMat(&fmat);
    delete [] BF;
    delete [] D;
    delete [] selected;
    delete [] next;
    for (int i = 0; i < num_backgrounds; ++ i) {
        cvReleaseMat(&gmat[i]);
        delete [] BG[i];
    }
    for(int i = 0; i < (int)Phi.size(); ++ i) {
        cvReleaseMat(&Phi[i]);
        cvReleaseMat(&intPhi[i]);
        cvReleaseMat(&Beta[i]);
        cvReleaseMat(&intBeta[i]);
    }
}

double DnbsClassBeta::reconerr(IplImage * img) {

    CvMat * intTarget
    = lyon::calcIntegral(img, 0, img->width, 0, img->height);

    vector<double> coeff = compCoeff(intTarget);

    double ssd = l2norm();
    for (int r = 0; r < img->height; ++ r)
        for (int c = 0; c < img->width; ++ c)
            ssd += img->nChannels == 3 
			? lyon::sqr((double)CV_IMAGE_ELEM(img, uchar, r, c * 3))
			: lyon::sqr((double)CV_IMAGE_ELEM(img, uchar, r, c));
    for(int i = 0; i < (int)coeff.size(); ++ i)
        ssd -= 2. * coeff[i]
               * indexOfBases[i].innerProd(
                   intTarget, 0, img->width, 0, img->height);
    cvReleaseMat(&intTarget);
    return sqrt(ssd / img->width / img->height);
}

double DnbsClassBeta::reconstruct_error(CvMat* tmplt) {
	CvMat * intTarget
		= lyon::calcIntegral(tmplt);
	vector<double> coeff
		= compCoeff(intTarget);
	double ans_l2norm = 0;
	for (int r = 0; r < tmplt->height; ++ r)
		for (int c = 0; c < tmplt->width; ++ c) {
			double curr = CV_MAT_ELEM(*tmplt, double, r, c);
			for (int i = 0; i < (int)indexOfBases.size(); ++ i) 
				curr -= coeff[i] * indexOfBases[i].element(r, c, tmplt->width);
			ans_l2norm += curr * curr;
		}
	return ans_l2norm;
}

double DnbsClassBeta::score(double lambda) {
	double ans_score = 0;
	for (int i = 0; i < (int)fgrdSample.size(); ++ i)
		ans_score += reconstruct_error(fgrdSample[i]) / (int)fgrdSample.size();
	for (int i = 0; i < (int)bgrdSample.size(); ++ i)
		ans_score -= lambda * reconstruct_error(bgrdSample[i]) / (int)bgrdSample.size();
	return ans_score;
}


char* double2string(double num, char* buf, int size) {
	sprintf(buf, "%.6lf", num);
	return buf;
}

// char* itoa(int num, char* buf, int size) {
// 	memset(buf, 0, sizeof(char) * size);
// 	sprintf(buf, "%d", num);
// 	return buf;
// }

/// Below is the extensions to the original discriminative nonorthogonal binary
/// subspace optimization techniques.

void DnbsClassBeta::computeDNBS_v1_pruning(
    int MaxNumBases,
    double lambda,
    vector<CvMat*> foreground,
    vector<CvMat*> background,
	int samp_width,
	int samp_height) {

    // Discriminative Non-orthogonal Binary Subspace
    // Formulation.
    // \arg\min_Phi{1/Nf\sum||fi-R_Phi{fi}||^2
    //              -lambda/Nb\sum_i||bi-R_Phi{bi}||^2}
    // where, lambda is a trade-off.
    // Ang Li, Nanjing University, Oct.7 2010.
    //
    // This edition integrates dense sampling and pruning
    // techniques in order to perform fast optimization of
    // discriminative nonorthogonal binary subspace decomposition.
    // Ang Li, Nanjing University, Nov. 10, 2010.

    // Supposing frame is gray-level image with 3 channels

    // Record the width and height
    int width = foreground[0]->width;
    int height = foreground[0]->height;

    int num_foregrounds = foreground.size();
    int num_backgrounds = background.size();

    // Definitions
    // K is the number of selected bases
    int K;
    // N is the total number of candidate bases
    int N = listOfBases.size();

    // B and D are two arrays of parameters in deciding the
    // optimized index of the next to-be-selected binary base.
    // i.e. selecting the one that maximising E_n(=|BF_n||BG_n|/D_n)
    vector<NBS_MAT_TYPE*> BF;
    for (int i = 0; i < num_foregrounds; ++ i) {
        NBS_MAT_TYPE * tempBF
        = new NBS_MAT_TYPE[N];
        BF.push_back(tempBF);
    }

    vector<NBS_MAT_TYPE*> BG;
    for (int i = 0; i < num_backgrounds; ++ i) {
        NBS_MAT_TYPE * tempBG
        = new NBS_MAT_TYPE[N];
        BG.push_back(tempBG);
    }
    NBS_MAT_TYPE * D = new NBS_MAT_TYPE[N];

    // Indicating whether the base is selected to form the subspace
    bool * selected = new bool[N];

    // Index of each selected binary bases
    vector<int> L;
    // Projection Coefficients
    vector<NBS_MAT_TYPE> Coeff;

    // Phi's and their integral image
    vector<CvMat*> Phi;
    vector<CvMat*> intPhi;

    // Inner product of each Phi[K] and image F(i.e. matPatch)
    vector< vector<NBS_MAT_TYPE> > innerPhiF;
    for (int i = 0; i < num_foregrounds; ++ i)
        innerPhiF.push_back(vector<NBS_MAT_TYPE>());
    // Inner product of each Phi[K] and image G
    vector< vector<NBS_MAT_TYPE> > innerPhiG;
    for (int i = 0; i < num_backgrounds; ++ i)
        innerPhiG.push_back(vector<NBS_MAT_TYPE>());

    // Beta's and their integral image
    vector<CvMat*> Beta;
    vector<CvMat*> intBeta;

    // Initialization
    // Reset memories
    memset(selected, 0, sizeof(bool) * N);
    // Computing the integral image of the original image
    vector<CvMat*> integralF;
    for (int i = 0; i < num_foregrounds; ++ i)
        integralF.push_back(calcIntegral(foreground[i]));
    vector<CvMat*> integralG;
    for (int i = 0; i < num_backgrounds; ++ i)
        integralG.push_back(calcIntegral(background[i]));

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

    // Process of initializing parameters and get the first choice
    int tempL = 0;
    double maxFunc = -inf;
    double curFunc;
    for(int i = 0; i < N; ++ i) {
        progress[i] = 0;
        for (int j = 0; j < num_foregrounds; ++ j)
            BF[j][i] = listOfBases[i].innerProd(integralF[j]);
        for (int j = 0; j < num_backgrounds; ++ j)
            BG[j][i] = listOfBases[i].innerProd(integralG[j]);
        D[i] = 1;

        curFunc = 0;
        for (int j = 0; j < num_foregrounds; ++ j)
            curFunc += BF[j][i] * BF[j][i] / num_foregrounds;
        for (int j = 0; j < num_backgrounds; ++ j)
            curFunc -= BG[j][i] * BG[j][i] / num_backgrounds * lambda;

        if (curFunc > maxFunc) {
            tempL = i;
            maxFunc = curFunc;
        }
    }
    L.push_back(tempL);
    selected[tempL] = true;

    // Computing parameters
    Phi.push_back(listOfBases[L[0]].toMatrix(height, width));
    intPhi.push_back(calcIntegral(Phi[0]));
    Beta.push_back(cvCloneMat(Phi[0]));
    intBeta.push_back(calcIntegral(Beta[0]));
    //Coeff.push_back(BF[L[0]]);
    for (int i = 0; i < num_foregrounds; ++ i)
        innerPhiF[i].push_back(innerProd(Phi[0], foreground[i]));
    for (int i = 0; i < num_backgrounds; ++ i)
        innerPhiG[i].push_back(innerProd(Phi[0], background[i]));

    int curtime = clock();

    K = 1;
    // Total Time Complexity is O(KN + K^2WH) = O(KWH(WH+K))
    while (K < MaxNumBases) {
        /**
        xmlDocPtr doc = xmlNewDoc(BAD_CAST"1.0");
        xmlNodePtr root_node = xmlNewNode(NULL, BAD_CAST"staterror");
        xmlDocSetRootElement(doc, root_node);
        static char buf[255];*/

        // Selecting the best fit Base
        maxFunc = -inf;
        tempL = -1;
        double tempPhiBase;
        for(int i = 0; i < sampN; ++ i) {
            //if (!selected[i]) {
            int ii = sampIndex[i];
            progress[ii] = K;
            tempPhiBase = listOfBases[ii].innerProd(intPhi[K - 1]);
            for (int j = 0; j < num_foregrounds; ++ j)
                BF[j][ii] -= innerPhiF[j][K - 1] * tempPhiBase;
            for (int j = 0; j < num_backgrounds; ++ j)
                BG[j][ii] -= innerPhiG[j][K - 1] * tempPhiBase;
            D[ii] -= tempPhiBase * tempPhiBase;

            curFunc = 0;
            for (int j = 0; j < num_foregrounds; ++ j)
                curFunc += BF[j][ii] * BF[j][ii] / num_foregrounds;
            for (int j = 0; j < num_backgrounds; ++ j)
                curFunc -= BG[j][ii] * BG[j][ii] / num_backgrounds * lambda;

            curFunc = fabs(curFunc) < sqeps ? 0 : (curFunc / D[ii]);

/**
            xmlNodePtr node = xmlNewNode(NULL, BAD_CAST "base");
            xmlAddChild(root_node, node);
            xmlNewProp(node, BAD_CAST "index", BAD_CAST itoa(i, buf, 255));
            xmlNewProp(node, BAD_CAST "func", BAD_CAST double2string(curFunc, buf, 255)); */

            if (curFunc > maxFunc) {
                maxFunc = curFunc;
                tempL = ii;
            }
        }

/**
        /// Save XML documents
        //Dumping document to stdio or file
        xmlSaveFormatFileEnc((string("exp/") + itoa(K, buf, 255) + string(".xml")).c_str(), doc, "UTF-8", 1);
        xmlFreeDoc(doc);
        xmlCleanupParser();
        xmlMemoryDump();//debug memory for regression tests */

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
                            for (int j = 0; j < num_foregrounds; ++ j)
                                BF[j][ii] -= innerPhiF[j][k] * tempPhiBase;
                            for (int j = 0; j < num_backgrounds; ++ j)
                                BG[j][ii] -= innerPhiG[j][k] * tempPhiBase;
                            D[ii] -= tempPhiBase * tempPhiBase;
                        }

                        progress[ii] = K;

                        curFunc = 0;
                        for (int j = 0; j < num_foregrounds; ++ j)
                            curFunc += BF[j][ii] * BF[j][ii] / num_foregrounds;
                        for (int j = 0; j < num_backgrounds; ++ j)
                            curFunc -= BG[j][ii] * BG[j][ii] / num_backgrounds * lambda;

                        curFunc = fabs(curFunc) < sqeps ? 0 : (curFunc / D[ii]);

                        if (curFunc > maxFunc) {
                            maxFunc = curFunc;
                            tempL = ii;
                        }
                    }

        L.push_back(tempL);
        selected[tempL] = true;

        // Computing parameters
        CvMat * tempMat = calcPhi(L[K], Phi, intPhi);       // O(KWH)
        Phi.push_back(calcMatDiv(tempMat, sqrt(D[L[K]])));  // O(WH)
        Beta.push_back(calcMatDiv(tempMat, D[L[K]]));       // O(WH)
        intBeta.push_back(calcIntegral(Beta[K]));           // O(WH)
        cvReleaseMat(&tempMat);
        intPhi.push_back(calcIntegral(Phi[K]));             // O(WH)
        for (int i = 0; i < num_foregrounds; ++ i)
            innerPhiF[i].push_back(innerProd(Phi[K], foreground[i]));
        for (int i = 0; i < num_backgrounds; ++ i)
            innerPhiG[i].push_back(innerProd(Phi[K], background[i]));

        ++ K;
    }

    lyon::settime((clock() - curtime) / (double)CLOCKS_PER_SEC);

    // Copy Answers
    this->indexOfBases.clear();
    for(int i = 0; i < (int)L.size(); ++ i)
        this->indexOfBases.push_back(listOfBases[L[i]]);
    CvMat * intTarget = calcIntegral(target);
    this->coefMatrix();
    this->coeffOfBases = compCoeff(intTarget);

    //demo();

    // Release Memories
    delete [] progress;
    delete [] D;
    delete [] selected;
    for (int i = 0; i < num_foregrounds; ++ i) {
        delete [] BF[i];
    }
    for (int i = 0; i < num_backgrounds; ++ i) {
        delete [] BG[i];
    }
    for(int i = 0; i < (int)Phi.size(); ++ i) {
        cvReleaseMat(&Phi[i]);
        cvReleaseMat(&intPhi[i]);
        cvReleaseMat(&Beta[i]);
        cvReleaseMat(&intBeta[i]);
    }
}

void DnbsClassBeta::computeDNBS_v1_downsample(
	int MaxNumBases,
	double lambda,
	vector<CvMat*> foreground,
	vector<CvMat*> background,
	int samp_width,
	int samp_height) {

		// Discriminative Non-orthogonal Binary Subspace
		// Formulation.
		// \arg\min_Phi{1/Nf\sum||fi-R_Phi{fi}||^2
		//              -lambda/Nb\sum_i||bi-R_Phi{bi}||^2}
		// where, lambda is a trade-off.
		// Ang Li, Nanjing University, Oct.7 2010.
		//
		// This edition integrates dense sampling and pruning
		// techniques in order to perform fast optimization of
		// discriminative nonorthogonal binary subspace decomposition.
		// Ang Li, Nanjing University, Nov. 10, 2010.

		// Supposing frame is gray-level image with 3 channels

		// Record the width and height
		int width = foreground[0]->width;
		int height = foreground[0]->height;

		int num_foregrounds = foreground.size();
		int num_backgrounds = background.size();

		// Definitions
		// K is the number of selected bases
		int K;
		// N is the total number of candidate bases
		int N = listOfBases.size();

		// B and D are two arrays of parameters in deciding the
		// optimized index of the next to-be-selected binary base.
		// i.e. selecting the one that maximising E_n(=|BF_n||BG_n|/D_n)
		vector<NBS_MAT_TYPE*> BF;
		for (int i = 0; i < num_foregrounds; ++ i) {
			NBS_MAT_TYPE * tempBF
				= new NBS_MAT_TYPE[N];
			BF.push_back(tempBF);
		}

		vector<NBS_MAT_TYPE*> BG;
		for (int i = 0; i < num_backgrounds; ++ i) {
			NBS_MAT_TYPE * tempBG
				= new NBS_MAT_TYPE[N];
			BG.push_back(tempBG);
		}
		NBS_MAT_TYPE * D = new NBS_MAT_TYPE[N];

		// Indicating whether the base is selected to form the subspace
		bool * selected = new bool[N];

		// Index of each selected binary bases
		vector<int> L;
		// Projection Coefficients
		vector<NBS_MAT_TYPE> Coeff;

		// Phi's and their integral image
		vector<CvMat*> Phi;
		vector<CvMat*> intPhi;

		// Inner product of each Phi[K] and image F(i.e. matPatch)
		vector< vector<NBS_MAT_TYPE> > innerPhiF;
		for (int i = 0; i < num_foregrounds; ++ i)
			innerPhiF.push_back(vector<NBS_MAT_TYPE>());
		// Inner product of each Phi[K] and image G
		vector< vector<NBS_MAT_TYPE> > innerPhiG;
		for (int i = 0; i < num_backgrounds; ++ i)
			innerPhiG.push_back(vector<NBS_MAT_TYPE>());

		// Beta's and their integral image
		vector<CvMat*> Beta;
		vector<CvMat*> intBeta;

		// Initialization
		// Reset memories
		memset(selected, 0, sizeof(bool) * N);
		// Computing the integral image of the original image
		vector<CvMat*> integralF;
		for (int i = 0; i < num_foregrounds; ++ i)
			integralF.push_back(calcIntegral(foreground[i]));
		vector<CvMat*> integralG;
		for (int i = 0; i < num_backgrounds; ++ i)
			integralG.push_back(calcIntegral(background[i]));

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
		memset(progress, -1, N * sizeof(int));

		// Process of initializing parameters and get the first choice
		int tempL = 0;
		double maxFunc = -inf;
		double curFunc;
		for(int i = 0; i < N; ++ i) {
			progress[i] = 0;
			for (int j = 0; j < num_foregrounds; ++ j)
				BF[j][i] = listOfBases[i].innerProd(integralF[j]);
			for (int j = 0; j < num_backgrounds; ++ j)
				BG[j][i] = listOfBases[i].innerProd(integralG[j]);
			D[i] = 1;

			curFunc = 0;
			for (int j = 0; j < num_foregrounds; ++ j)
				curFunc += BF[j][i] * BF[j][i] / num_foregrounds;
			for (int j = 0; j < num_backgrounds; ++ j)
				curFunc -= BG[j][i] * BG[j][i] / num_backgrounds * lambda;

			if (curFunc > maxFunc) {
				tempL = i;
				maxFunc = curFunc;
			}
		}
		L.push_back(tempL);
		selected[tempL] = true;

		// Computing parameters
		Phi.push_back(listOfBases[L[0]].toMatrix(height, width));
		intPhi.push_back(calcIntegral(Phi[0]));
		Beta.push_back(cvCloneMat(Phi[0]));
		intBeta.push_back(calcIntegral(Beta[0]));
		//Coeff.push_back(BF[L[0]]);
		for (int i = 0; i < num_foregrounds; ++ i)
			innerPhiF[i].push_back(innerProd(Phi[0], foreground[i]));
		for (int i = 0; i < num_backgrounds; ++ i)
			innerPhiG[i].push_back(innerProd(Phi[0], background[i]));

		int curtime = clock();

		K = 1;
		// Total Time Complexity is O(KN + K^2WH) = O(KWH(WH+K))
		while (K < MaxNumBases) {
			/**
			xmlDocPtr doc = xmlNewDoc(BAD_CAST"1.0");
			xmlNodePtr root_node = xmlNewNode(NULL, BAD_CAST"staterror");
			xmlDocSetRootElement(doc, root_node);
			static char buf[255];*/

			// Selecting the best fit Base
			maxFunc = -inf;
			tempL = -1;
			double tempPhiBase;
			for(int i = 0; i < sampN; ++ i) {
				//if (!selected[i]) {
				int ii = sampIndex[i];
				progress[ii] = K;
				tempPhiBase = listOfBases[ii].innerProd(intPhi[K - 1]);
				for (int j = 0; j < num_foregrounds; ++ j)
					BF[j][ii] -= innerPhiF[j][K - 1] * tempPhiBase;
				for (int j = 0; j < num_backgrounds; ++ j)
					BG[j][ii] -= innerPhiG[j][K - 1] * tempPhiBase;
				D[ii] -= tempPhiBase * tempPhiBase;

				curFunc = 0;
				for (int j = 0; j < num_foregrounds; ++ j)
					curFunc += BF[j][ii] * BF[j][ii] / num_foregrounds;
				for (int j = 0; j < num_backgrounds; ++ j)
					curFunc -= BG[j][ii] * BG[j][ii] / num_backgrounds * lambda;

				curFunc = fabs(curFunc) < sqeps ? 0 : (curFunc / D[ii]);

				/**
				xmlNodePtr node = xmlNewNode(NULL, BAD_CAST "base");
				xmlAddChild(root_node, node);
				xmlNewProp(node, BAD_CAST "index", BAD_CAST itoa(i, buf, 255));
				xmlNewProp(node, BAD_CAST "func", BAD_CAST double2string(curFunc, buf, 255)); */

				if (curFunc > maxFunc) {
					maxFunc = curFunc;
					tempL = ii;
				}
			}

			if (tempL == -1) break;

			L.push_back(tempL);
			selected[tempL] = true;

			// Computing parameters
			CvMat * tempMat = calcPhi(L[K], Phi, intPhi);       // O(KWH)
			Phi.push_back(calcMatDiv(tempMat, sqrt(D[L[K]])));  // O(WH)
			Beta.push_back(calcMatDiv(tempMat, D[L[K]]));       // O(WH)
			intBeta.push_back(calcIntegral(Beta[K]));           // O(WH)
			cvReleaseMat(&tempMat);
			intPhi.push_back(calcIntegral(Phi[K]));             // O(WH)
			for (int i = 0; i < num_foregrounds; ++ i)
				innerPhiF[i].push_back(innerProd(Phi[K], foreground[i]));
			for (int i = 0; i < num_backgrounds; ++ i)
				innerPhiG[i].push_back(innerProd(Phi[K], background[i]));

			++ K;
		}

		lyon::settime((clock() - curtime) / (double)CLOCKS_PER_SEC);

		// Copy Answers
		this->indexOfBases.clear();
		for(int i = 0; i < (int)L.size(); ++ i)
			this->indexOfBases.push_back(listOfBases[L[i]]);
		CvMat * intTarget = calcIntegral(target);
		this->coefMatrix();
		this->coeffOfBases = compCoeff(intTarget);

		//demo();

		// Release Memories
		delete [] progress;
		delete [] D;
		delete [] selected;
		for (int i = 0; i < num_foregrounds; ++ i) {
			delete [] BF[i];
		}
		for (int i = 0; i < num_backgrounds; ++ i) {
			delete [] BG[i];
		}
		for(int i = 0; i < (int)Phi.size(); ++ i) {
			cvReleaseMat(&Phi[i]);
			cvReleaseMat(&intPhi[i]);
			cvReleaseMat(&Beta[i]);
			cvReleaseMat(&intBeta[i]);
		}
}

void DnbsClassBeta::genHaarBases_coherence(double mu_coherence) {
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
}

void DnbsClassBeta::computeDNBS_v1_pruning_2nd(
	int MaxNumBases,
	double lambda,
	vector<CvMat*> foreground,
	vector<CvMat*> background,
	int samp_width,
	int samp_height) {

		// Supposing frame is gray-level image with 3 channels

		// Record the width and height
		int width = foreground[0]->width;
		int height = foreground[0]->height;

		int num_foregrounds = foreground.size();
		int num_backgrounds = background.size();

		// Definitions
		// K is the number of selected bases
		int K;
		// N is the total number of candidate bases
		int N = listOfBases.size();

		// B and D are two arrays of parameters in deciding the
		// optimized index of the next to-be-selected binary base.
		// i.e. selecting the one that maximising E_n(=|BF_n||BG_n|/D_n)
		vector<NBS_MAT_TYPE*> BF;
		for (int i = 0; i < num_foregrounds; ++ i) {
			NBS_MAT_TYPE * tempBF
				= new NBS_MAT_TYPE[N];
			BF.push_back(tempBF);
		}

		vector<NBS_MAT_TYPE*> BG;
		for (int i = 0; i < num_backgrounds; ++ i) {
			NBS_MAT_TYPE * tempBG
				= new NBS_MAT_TYPE[N];
			BG.push_back(tempBG);
		}
		NBS_MAT_TYPE * D = new NBS_MAT_TYPE[N];

		// Indicating whether the base is selected to form the subspace
		bool * selected = new bool[N];

		// Index of each selected binary bases
		vector<int> L;
		// Projection Coefficients
		vector<NBS_MAT_TYPE> Coeff;

		// Phi's and their integral image
		vector<CvMat*> Phi;
		vector<CvMat*> intPhi;

		// Inner product of each Phi[K] and image F(i.e. matPatch)
		vector< vector<NBS_MAT_TYPE> > innerPhiF;
		for (int i = 0; i < num_foregrounds; ++ i)
			innerPhiF.push_back(vector<NBS_MAT_TYPE>());
		// Inner product of each Phi[K] and image G
		vector< vector<NBS_MAT_TYPE> > innerPhiG;
		for (int i = 0; i < num_backgrounds; ++ i)
			innerPhiG.push_back(vector<NBS_MAT_TYPE>());

		// Beta's and their integral image
		vector<CvMat*> Beta;
		vector<CvMat*> intBeta;

		// Initialization
		// Reset memories
		memset(selected, 0, sizeof(bool) * N);
		// Computing the integral image of the original image
		vector<CvMat*> integralF;
		for (int i = 0; i < num_foregrounds; ++ i)
			integralF.push_back(calcIntegral(foreground[i]));
		vector<CvMat*> integralG;
		for (int i = 0; i < num_backgrounds; ++ i)
			integralG.push_back(calcIntegral(background[i]));

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

		// Process of initializing parameters and get the first choice
		int tempL = 0;
		double maxFunc = -inf;
		double curFunc;
		for(int i = 0; i < N; ++ i) {
			progress[i] = 0;
			for (int j = 0; j < num_foregrounds; ++ j)
				BF[j][i] = listOfBases[i].innerProd(integralF[j]);
			for (int j = 0; j < num_backgrounds; ++ j)
				BG[j][i] = listOfBases[i].innerProd(integralG[j]);
			D[i] = 1;

			curFunc = 0;
			for (int j = 0; j < num_foregrounds; ++ j)
				curFunc += BF[j][i] * BF[j][i] / num_foregrounds;
			for (int j = 0; j < num_backgrounds; ++ j)
				curFunc -= BG[j][i] * BG[j][i] / num_backgrounds * lambda;

			if (curFunc > maxFunc) {
				tempL = i;
				maxFunc = curFunc;
			}
		}
		L.push_back(tempL);
		selected[tempL] = true;

		// Computing parameters
		Phi.push_back(listOfBases[L[0]].toMatrix(height, width));
		intPhi.push_back(calcIntegral(Phi[0]));
		Beta.push_back(cvCloneMat(Phi[0]));
		intBeta.push_back(calcIntegral(Beta[0]));
		//Coeff.push_back(BF[L[0]]);
		for (int i = 0; i < num_foregrounds; ++ i)
			innerPhiF[i].push_back(innerProd(Phi[0], foreground[i]));
		for (int i = 0; i < num_backgrounds; ++ i)
			innerPhiG[i].push_back(innerProd(Phi[0], background[i]));

		int curtime = clock();

		double * sampScore = new double[sampN];

		K = 1;
		// Total Time Complexity is O(KN + K^2WH) = O(KWH(WH+K))
		while (K < MaxNumBases) {
			/**
			xmlDocPtr doc = xmlNewDoc(BAD_CAST"1.0");
			xmlNodePtr root_node = xmlNewNode(NULL, BAD_CAST"staterror");
			xmlDocSetRootElement(doc, root_node);
			static char buf[255];*/

			// Selecting the best fit Base
			maxFunc = -inf;
			tempL = -1;
			double tempPhiBase;
			for(int i = 0; i < sampN; ++ i) {
				//if (!selected[i]) {
				int ii = sampIndex[i];
				progress[ii] = K;
				tempPhiBase = listOfBases[ii].innerProd(intPhi[K - 1]);
				for (int j = 0; j < num_foregrounds; ++ j)
					BF[j][ii] -= innerPhiF[j][K - 1] * tempPhiBase;
				for (int j = 0; j < num_backgrounds; ++ j)
					BG[j][ii] -= innerPhiG[j][K - 1] * tempPhiBase;
				D[ii] -= tempPhiBase * tempPhiBase;

				curFunc = 0;
				for (int j = 0; j < num_foregrounds; ++ j)
					curFunc += BF[j][ii] * BF[j][ii] / num_foregrounds;
				for (int j = 0; j < num_backgrounds; ++ j)
					curFunc -= BG[j][ii] * BG[j][ii] / num_backgrounds * lambda;

				curFunc = fabs(curFunc) < sqeps ? 0 : (curFunc / D[ii]);

				sampScore[i] = curFunc;

				/**
				xmlNodePtr node = xmlNewNode(NULL, BAD_CAST "base");
				xmlAddChild(root_node, node);
				xmlNewProp(node, BAD_CAST "index", BAD_CAST itoa(i, buf, 255));
				xmlNewProp(node, BAD_CAST "func", BAD_CAST double2string(curFunc, buf, 255)); */

				if (curFunc > maxFunc) {
					maxFunc = curFunc;
					tempL = ii;
				}
			}

			/**
			/// Save XML documents
			//Dumping document to stdio or file
			xmlSaveFormatFileEnc((string("exp/") + itoa(K, buf, 255) + string(".xml")).c_str(), doc, "UTF-8", 1);
			xmlFreeDoc(doc);
			xmlCleanupParser();
			xmlMemoryDump();//debug memory for regression tests */

			if (tempL == -1) break;

			static const double maxFuncScale = 0.5;
			// check features again
			for (int i = 0; i < sampN; ++ i)
				if (sampScore[i] > maxFunc - fabs(maxFunc) * maxFuncScale) {

					int optIndex = sampIndex[i];
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
										for (int j = 0; j < num_foregrounds; ++ j)
											BF[j][ii] -= innerPhiF[j][k] * tempPhiBase;
										for (int j = 0; j < num_backgrounds; ++ j)
											BG[j][ii] -= innerPhiG[j][k] * tempPhiBase;
										D[ii] -= tempPhiBase * tempPhiBase;
									}

									progress[ii] = K;

									curFunc = 0;
									for (int j = 0; j < num_foregrounds; ++ j)
										curFunc += BF[j][ii] * BF[j][ii] / num_foregrounds;
									for (int j = 0; j < num_backgrounds; ++ j)
										curFunc -= BG[j][ii] * BG[j][ii] / num_backgrounds * lambda;

									curFunc = fabs(curFunc) < sqeps ? 0 : (curFunc / D[ii]);

									if (curFunc > maxFunc) {
										maxFunc = curFunc;
										tempL = ii;
									}
								}
				}

				L.push_back(tempL);
				selected[tempL] = true;

				// Computing parameters
				CvMat * tempMat = calcPhi(L[K], Phi, intPhi);       // O(KWH)
				Phi.push_back(calcMatDiv(tempMat, sqrt(D[L[K]])));  // O(WH)
				Beta.push_back(calcMatDiv(tempMat, D[L[K]]));       // O(WH)
				intBeta.push_back(calcIntegral(Beta[K]));           // O(WH)
				cvReleaseMat(&tempMat);
				intPhi.push_back(calcIntegral(Phi[K]));             // O(WH)
				for (int i = 0; i < num_foregrounds; ++ i)
					innerPhiF[i].push_back(innerProd(Phi[K], foreground[i]));
				for (int i = 0; i < num_backgrounds; ++ i)
					innerPhiG[i].push_back(innerProd(Phi[K], background[i]));

				++ K;
		}

		lyon::settime((clock() - curtime) / (double)CLOCKS_PER_SEC);

		// Copy Answers
		this->indexOfBases.clear();
		for(int i = 0; i < (int)L.size(); ++ i)
			this->indexOfBases.push_back(listOfBases[L[i]]);
		CvMat * intTarget = calcIntegral(target);
		this->coefMatrix();
		this->coeffOfBases = compCoeff(intTarget);

		//demo();

		// Release Memories
		delete [] progress;
		delete [] D;
		delete [] selected;
		for (int i = 0; i < num_foregrounds; ++ i) {
			delete [] BF[i];
		}
		for (int i = 0; i < num_backgrounds; ++ i) {
			delete [] BG[i];
		}
		for(int i = 0; i < (int)Phi.size(); ++ i) {
			cvReleaseMat(&Phi[i]);
			cvReleaseMat(&intPhi[i]);
			cvReleaseMat(&Beta[i]);
			cvReleaseMat(&intBeta[i]);
		}
}