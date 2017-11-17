// DnbsTracker.cpp 
// This source file contains only one function:
//		int startDnbsTracker(
//			VideoStruct & videoSetting,
//			DnbsParamStruct & dnbsParam);
// which starts the D-NBS tracking process.
//
// Ang Li (mailto:nju.angli@gmail.com)
// Dept. Computer Science and Technology
// Nanjing University
// Nov.8, 2010.

#include "headers.h"
#include "DnbsTracker.h"
#include "LyonLib.h"
#include "DnbsClassAlpha.h"
#include "DnbsClassGamma.h"
#include "commonlib.h"
#include "utilities.h"

struct compactStruct {
	VideoStruct * pVideoSetting;
	DnbsParamStruct * pDnbsParam;
};

//UINT dnbsTrackerThread(LPVOID pParam);
int coreDnbsTracker(VideoStruct & videoSetting,
					DnbsParamStruct & dnbsParam);

//static HANDLE hMutex = NULL;
static IplImage * g_showframe = NULL;
static bool g_esc_thread = false;

int startDnbsTracker(VideoStruct &videoSetting, DnbsParamStruct &dnbsParam) {
    
    coreDnbsTracker(videoSetting, dnbsParam);
    
    /*
	compactStruct *thCompactStruct = new compactStruct;
	thCompactStruct->pVideoSetting = &videoSetting;
	thCompactStruct->pDnbsParam = &dnbsParam;

	hMutex = CreateMutex(NULL, FALSE, NULL);

	CWinThread * pThread 
		= AfxBeginThread(dnbsTrackerThread, (LPVOID)thCompactStruct);
	
	string win_title = videoSetting.title;

	while (!g_esc_thread) {
		WaitForSingleObject(hMutex, INFINITE);
		cvShowImage(win_title.c_str(), g_showframe);
		ReleaseMutex(hMutex);
		char c = cvWaitKey(10);
		if (c == 27) break;
	}

	g_esc_thread = true;

	cvWaitKey();
	cvDestroyWindow(win_title.c_str());

	delete thCompactStruct;

	WaitForSingleObject(pThread->m_hThread,INFINITE);*/
	return 0;
}

//UINT dnbsTrackerThread(LPVOID pParam) {
//	compactStruct* pCompactStruct
//		= (compactStruct*)pParam;
//	coreDnbsTracker(
//		*(pCompactStruct->pVideoSetting),
//		*(pCompactStruct->pDnbsParam));
//
//	return 0;
//}

xmlNodePtr createFrameXMLNode(xmlNodePtr& father, int framenumber, CvRect boundbox) {
	xmlNodePtr newNode
		= xmlNewNode(NULL, BAD_CAST"frame");
	xmlAddChild(father, newNode);
	char buf[255];
	xmlNewProp(newNode, BAD_CAST"num", BAD_CAST itoa(framenumber, buf, 255));
	xmlNodePtr boxNode
		= xmlNewNode(NULL, BAD_CAST"bound");
	xmlAddChild(newNode, boxNode);
	xmlNewProp(boxNode, BAD_CAST "x", BAD_CAST itoa(boundbox.x, buf, 255));
	xmlNewProp(boxNode, BAD_CAST "y", BAD_CAST itoa(boundbox.y, buf, 255));
	xmlNewProp(boxNode, BAD_CAST "width", BAD_CAST itoa(boundbox.width, buf, 255));
	xmlNewProp(boxNode, BAD_CAST "height", BAD_CAST itoa(boundbox.height, buf, 255));

	return newNode;
}

//#include <Shlwapi.h>

int coreDnbsTracker(VideoStruct & videoSetting,
                     DnbsParamStruct & dnbsParam) {

    DnbsClassGamma nbsobj;

	// Record tracking results to XML
	xmlDocPtr xmldoc = NULL;
	xmlNodePtr rootNode = NULL;
	if (videoSetting.save_resxml) {
		xmldoc = xmlNewDoc(BAD_CAST"1.0");
		rootNode = xmlNewNode(NULL, BAD_CAST"trackresult");
		xmlNewProp(rootNode, BAD_CAST"title", BAD_CAST videoSetting.title.c_str());
		xmlDocSetRootElement(xmldoc, rootNode);
	}

	CvCapture * pVideoCapture = NULL;
	int current_frame = 0;
	IplImage * qframe = NULL;
	IplImage * frame = NULL;
	IplImage * showframe = NULL;
	char imagepath[255];
	int frame_width;
	int frame_height;

	if (videoSetting.srctype == 0) { 
		// take in video
		pVideoCapture = cvCreateFileCapture(videoSetting.videoPath.c_str());
		while (current_frame <= videoSetting.frameCreate) {
			int err = cvGrabFrame(pVideoCapture);
			if (err == 0) {
				printf("ERROR: Unable to Grab Frames\r\n");
				return 0;
			}
			++ current_frame;
		}
		-- current_frame;
		qframe = cvRetrieveFrame(pVideoCapture);
		frame = lyon::create_gray_from_rgb(qframe);
		showframe = cvCloneImage(frame);
	} else { 
		// take in images
		current_frame = videoSetting.frameCreate;
		sprintf(imagepath, videoSetting.imagePattern.c_str(), current_frame);
		printf("path: %s\n", imagepath);
		qframe = cvLoadImage(imagepath);
		frame = lyon::create_gray_from_rgb(qframe);
		showframe = cvCloneImage(frame);
		cvReleaseImage(&qframe);
	}

	g_showframe = cvCloneImage(showframe);

	frame_width = frame->width;
	frame_height = frame->height;
	CvRect boundbox = videoSetting.boundbox;
	int margin_width = videoSetting.marginWidth;
	int margin_height = videoSetting.marginHeight;
	double lambda
		= dnbsParam.lambda;
	int num_of_foregrounds
		= dnbsParam.numForegrounds;
	int num_of_backgrounds
		= dnbsParam.numBackgrounds;
	int num_of_base
		= dnbsParam.numBases;
	int version
		= dnbsParam.version;
	int distance_threshold = 5;
	double gamma = 0.5;
	int samp_width = dnbsParam.samp_width;
	int samp_height = dnbsParam.samp_height;
	double min_innprod = dnbsParam.min_innprod;
	double searchratio = dnbsParam.searchratio;

    IplImage * patch
        = lyon::selectedPatch(frame, boundbox);
    CvMat * target
        = lyon::image2mat(patch);
	CvMat * intInitial
		= lyon::calcIntegral(target);
    cvReleaseImage(&patch);

	CvVideoWriter * pVideoWriter = NULL;
	string pathSave = videoSetting.path_save_resavi;
	enum {NOSAVE, SAVE2VIDEO, SAVE2IMAGES} savetype = NOSAVE;
	if (videoSetting.save_resavi) {
		if (dir_exists(pathSave)) {
			savetype = SAVE2IMAGES;
		} else savetype = SAVE2VIDEO;
	} else savetype = NOSAVE;
	if (savetype == SAVE2VIDEO) {
		printf("Saving Result AVI\n");
		pVideoWriter 
			= cvCreateVideoWriter(
				videoSetting.path_save_resavi.c_str(),
				CV_FOURCC('X','V','I','D'),
				20,
				cvSize(frame_width, frame_height),
				1);
	}

    lyon::drawRect(showframe, boundbox, CV_RGB(0xff, 0, 0));

	CvMat ** previous_foreground
		= (CvMat**)malloc(sizeof(CvMat*) * (num_of_foregrounds));
	for (int i = 0; i < num_of_foregrounds; ++ i)
		previous_foreground[i] = NULL;
	int float_pointer = 0;

	// Time statistics
	double timestat_preproc = 0;
	double timestat_train = 0;
	int number_of_train = 0;
	double timestat_track = 0;
	int number_of_track = 0;
	int currtime;

	nbsobj.importTarget(target);

	currtime = clock();
	nbsobj.genHaarBases();
	if (version == 5 || version == 6 || version == 7 || version == 9 || version == 10)
		nbsobj.genHaarBases_innProdFast(min_innprod);
	timestat_preproc = (clock() - currtime) / (double)CLOCKS_PER_SEC;

	if (num_of_foregrounds > 0) {
		if (previous_foreground[float_pointer] != NULL)
			cvReleaseMat(&previous_foreground[float_pointer]);
		previous_foreground[float_pointer] = target;
		float_pointer = (float_pointer + 1) % (num_of_foregrounds);
	}

	//nbsobj.NaiveNbsClass::computeNBS(50);
	//vector<HaarBase> baseInitial = nbsobj.listOfBases;
	currtime = clock();
    nbsobj.computeDNBS(
        num_of_base, num_of_foregrounds, num_of_backgrounds, frame, 
		boundbox, margin_width + 20, margin_height + 20,
        distance_threshold, lambda, version, samp_width, samp_height, 
		searchratio, previous_foreground);
	timestat_train += (clock() - currtime) / (double)CLOCKS_PER_SEC;
	++ number_of_train;
    printf("time cost: %.6lf\n", (clock() - currtime) / (double)CLOCKS_PER_SEC);

    nbsobj.coefMatrix();
    for (int i = 0; i < (int)nbsobj.bgsample.size(); ++ i)
		lyon::drawRect(showframe, nbsobj.bgsample[i], cvScalar(0xff, 0, 0));

//	WaitForSingleObject(hMutex, INFINITE);
	cvCopyImage(showframe, g_showframe);
    cvShowImage(videoSetting.title.c_str(), g_showframe);
    cvWaitKey(33);
//	ReleaseMutex(hMutex);

// 	cvShowImage("result", showframe);
// 	cvWaitKey();

	if (xmldoc != NULL) {
		createFrameXMLNode(rootNode, current_frame, boundbox);
	}

    if (pVideoWriter != NULL)
        cvWriteFrame(pVideoWriter, showframe);
	if (savetype == SAVE2IMAGES) {
		char filename[255];
        sprintf(filename, "%s/%04d.jpg", pathSave.c_str(), current_frame);
		cvSaveImage(filename, showframe);
	}

    int step = 0;
    CvMat * first_template = cvCloneMat(target);

	FILE * fstatSSD = fopen("fstatSSD.txt", "w+");

	bool newly_updated = true;
	double markSSD;

	int left = max(0, boundbox.x - margin_width);
	int right = min(frame->width - boundbox.width, boundbox.x + margin_width);
	int top = max(0, boundbox.y - margin_height);
	int bottom = min(frame->height - boundbox.height, boundbox.y + margin_height);
	struct FrameLog {
		CvRect boundbox;
		int left, right, top, bottom;
		IplImage* original;
		CvMat* integral;
		CvMat* integral2;
		FrameLog(): original(0), integral(0), integral2(0) {}
		~FrameLog() {
//			if (original) cvReleaseImage(&original);
//			if (integral) cvReleaseMat(&integral);
//			if (integral2) cvReleaseMat(&integral2);
		}
	};
	static const int maxFrameCache = 5;
	queue< FrameLog > frameCache;
	FrameLog framelog;
	framelog.left = left; framelog.right = right;
	framelog.top = top; framelog.bottom = bottom;
	framelog.boundbox = boundbox;
	framelog.original = cvCloneImage(frame);
	framelog.integral
		= lyon::calcIntegral(
		frame, left, right + boundbox.width,
		top, bottom + boundbox.height);
	framelog.integral2
		= lyon::calcSqrIntegral(
		frame, left, right + boundbox.width,
		top, bottom + boundbox.height);
	frameCache.push(framelog);

	bool losetrack = false;
	int confident_cnt = 0;

	double timeframe = 0;
	int avenumframe = 0;
	double avetimeframe = 0;

    while (true) {
		timeframe = clock();
        ++ current_frame;
        if (videoSetting.frameFinish > 0 && current_frame > videoSetting.frameFinish)
            break;

		if (videoSetting.srctype == 0) {
			qframe = cvQueryFrame(pVideoCapture);
			if (qframe == NULL) break;
			frame = lyon::create_gray_from_rgb(qframe);
			showframe = cvCloneImage(frame);
		} else {
			sprintf(imagepath, videoSetting.imagePattern.c_str(), current_frame);
			qframe = cvLoadImage(imagepath);
			if (qframe == NULL) break;
			frame = lyon::create_gray_from_rgb(qframe);
			showframe = cvCloneImage(frame);
			cvReleaseImage(&qframe);
		}

		currtime = clock();

		CvMat * intFrame
			= lyon::calcIntegral(
			frame, left, right + boundbox.width,
			top, bottom + boundbox.height);

		CvMat * int2Frame
			= lyon::calcSqrIntegral(
			frame, left, right + boundbox.width,
			top, bottom + boundbox.height);

		CvRect prevBoundBox = boundbox;
        boundbox = nbsobj.updateBoxSSD_simp(
			frame, boundbox, left, right, top, bottom, intFrame, int2Frame, true);
		// Origin: boundbox = nbsobj.updateBoxSSD(frame, boundbox, margin_width, margin_height);
		timestat_track += (clock() - currtime) / (double)CLOCKS_PER_SEC;
		++ number_of_track;

        IplImage * foreground
			= lyon::selectedPatch(frame, boundbox);
		CvMat * reference_template
			= lyon::image2mat(foreground);
		CvMat * target
			= cvCloneMat(reference_template);

// 		if (losetrack) {
// 			//boundbox = prevBoundBox;
// 			static const int lost_margin_width = 200;
// 			static const int lost_margin_height = 200;
// 
// 			left = max(0, boundbox.x - lost_margin_width);
// 			right = min(frame->width - boundbox.width, boundbox.x + lost_margin_width);
// 			top = max(0, boundbox.y - lost_margin_height);
// 			bottom = min(frame->height - boundbox.height, boundbox.y + lost_margin_height);
// 		} else {
			left = max(0, boundbox.x - margin_width);
			right = min(frame->width - boundbox.width, boundbox.x + margin_width);
			top = max(0, boundbox.y - margin_height);
			bottom = min(frame->height - boundbox.height, boundbox.y + margin_height);
//		}

		// Check if losing track
		CvMat * int_foreground
			= lyon::calcIntegral(target);
		//nbsobj.importTarget(reference_template);
		//vector<HaarBase> prevBases = nbsobj.listOfBases;
		vector<double> prevCoefficients = nbsobj.coeffOfBases;
		nbsobj.coeffOfBases = nbsobj.compCoeff(int_foreground);
		FrameLog & backLog = frameCache.back();
		CvRect back_boundbox = nbsobj.updateBoxSSD_simp(
			backLog.original, backLog.boundbox, 
			backLog.left, backLog.right, backLog.top, backLog.bottom, 
			backLog.integral, backLog.integral2, true);
		double back_error 
			= sqrt(Sqr(back_boundbox.x - backLog.boundbox.x) 
			+ Sqr(back_boundbox.y - backLog.boundbox.y));
		printf("Backward error = %lf\n", back_error);
		static const double thresh_backward_error = .5 * min(boundbox.width, boundbox.height);
		static const int max_confident_cnt = 3;
		if (back_error > thresh_backward_error) {
			losetrack = true;
			//nbsobj.coeffOfBases = prevCoefficients;
			//nbsobj.listOfBases = baseInitial;
			//nbsobj.coefMatrix();
// 			nbsobj.importPatch(frameCache.back().original);
// 			nbsobj.computeDNBS(
// 				num_of_base, num_of_foregrounds, num_of_backgrounds, frame,
// 				boundbox, margin_width + 20, margin_height + 20,
// 				distance_threshold, lambda, version, samp_width, samp_height, 
// 				searchratio, previous_foreground);
// 			nbsobj.coefMatrix();
			//nbsobj.coeffOfBases = nbsobj.compCoeff(intInitial);
			//nbsobj.coeffOfBases = nbsobj.compCoeff(frameCache.front().integral);
			nbsobj.coeffOfBases = prevCoefficients;
			confident_cnt = max_confident_cnt;
			goto LABEL_NOUPDATE;
		} else if (losetrack) {
			-- confident_cnt;
			if (confident_cnt == 0)
				losetrack = false;
			else {
				nbsobj.coeffOfBases = prevCoefficients;
				//nbsobj.listOfBases = baseInitial;
				//nbsobj.coefMatrix();
				//nbsobj.coeffOfBases = nbsobj.compCoeff(intInitial);
				//nbsobj.coeffOfBases = nbsobj.compCoeff(frameCache.front().integral);
				goto LABEL_NOUPDATE;
			}
		}

        {
		// Log frame info
		framelog.left = left; framelog.right = right;
		framelog.top = top; framelog.bottom = bottom;
		framelog.boundbox = boundbox;
		framelog.original = cvCloneImage(frame);
		framelog.integral
			= lyon::calcIntegral(
			frame, left, right + boundbox.width,
			top, bottom + boundbox.height);
		framelog.integral2
			= lyon::calcSqrIntegral(
			frame, left, right + boundbox.width,
			top, bottom + boundbox.height);
		frameCache.push(framelog);
		while (frameCache.size() > maxFrameCache) frameCache.pop();

#define _TEMPLATE_UPDATE_ 0
#if (_TEMPLATE_UPDATE_ == 0)
		for (int r = 0; r < reference_template->rows; ++ r)
			for (int c = 0; c < reference_template->cols; ++ c)
				CV_MAT_ELEM(*reference_template, double, r, c)
				= (1-gamma) * CV_MAT_ELEM(*reference_template, double, r, c)
				+ gamma * CV_MAT_ELEM(*first_template, double, r, c);
#elif (_TEMPLATE_UPDATE == 1)
		static const double beta = 0.3;
		for (int r = 0; r < reference_template->rows; ++ r)
			for (int c = 0; c < reference_template->cols; ++ c)
				CV_MAT_ELEM(*reference_template, double, r, c)
				= beta * CV_MAT_ELEM(*reference_template, double, r, c)
				+ (1 - beta) * CV_MAT_ELEM(*nbsobj.target, double, r, c);
#endif

        CvMat * intReference
            = lyon::calcIntegral(reference_template);
        nbsobj.importTarget(reference_template);
        nbsobj.coeffOfBases = nbsobj.compCoeff(intReference);

		if (num_of_foregrounds > 0) {
			if (previous_foreground[float_pointer] != NULL)
				cvReleaseMat(&previous_foreground[float_pointer]);
			previous_foreground[float_pointer] = target;
			float_pointer = (float_pointer + 1) % (num_of_foregrounds);
		}

        ++ step;
#define _SUBSPACE_UPDATE_ 1
#if (_SUBSPACE_UPDATE_ == 1)
		cvReleaseMat(&intFrame);
		cvReleaseMat(&int2Frame);
        if (step % 5 == 0 ) {
#elif (_SUBSPACE_UPDATE_ == 2)
		// Reconstruction Error
		double reconstruct_error 
			= sqrt(nbsobj.reconstruct_error(target) / target->width / target->height);
		printf("New Reconstruct. Error = %lf\n", reconstruct_error);
		// Background Distracting Power
		// 'Nd' samples with lowest SSD values and 
		// overlap area less than 'maxOverlapRatio' with the target.
		int Nd = 5;
		double maxOverlapRatio = 0.3;
		double nbsNorm = nbsobj.l2norm();
		double globalMinSSD = lyon::inf;
		static vector<CvPoint> sampleDistractors;
		sampleDistractors.clear();
		sampleDistractors.push_back(cvPoint(boundbox.x, boundbox.y));
		double selfSSD = 0;
		selfSSD += lyon::calcIntSum(
			int2Frame, boundbox.x - left, boundbox.x - left + boundbox.width,
			boundbox.y - top, boundbox.y - top + boundbox.height);
		for(int j = 0; j < (int)nbsobj.coeffOfBases.size(); ++ j)
			selfSSD -= 2. * nbsobj.coeffOfBases[j] 
					* nbsobj.indexOfBases[j].innerProd(
						intFrame, boundbox.x - left, boundbox.x - left + boundbox.width, 
						boundbox.y - top, boundbox.y - top + boundbox.height);
		selfSSD += nbsNorm;
		selfSSD /= boundbox.width * boundbox.height;
		for (int i = 0; i < (int)nbsobj.bksimilar.size() && Nd; ++ i) {
			CvPoint & topleft = nbsobj.bksimilar[i].second;
			bool nonSupressed = true;
			for (int j = 0; j < sampleDistractors.size(); ++ j) {
				int width = boundbox.width - abs(sampleDistractors[j].x - topleft.x);
				int height = boundbox.height - abs(sampleDistractors[j].y - topleft.y);
				if (width * height > maxOverlapRatio * boundbox.width * boundbox.height) {
					nonSupressed = false;
					break;
				}
			}
			if (nonSupressed) {
				lyon::drawRect(showframe, 
					cvRect(topleft.x, topleft.y, boundbox.width, boundbox.height), 
					CV_RGB(0, 0xff, 0));
				// Overlap area less than 'maxOverlapRatio' w.r.t the target.
				int r = topleft.y;
				int c = topleft.x;
				double ssd = 0;
				ssd += lyon::calcIntSum(
					int2Frame, c - left, c - left + boundbox.width,
					r - top, r - top + boundbox.height);
				for(int j = 0; j < (int)nbsobj.coeffOfBases.size(); ++ j)
					ssd -= 2. * nbsobj.coeffOfBases[j] 
					* nbsobj.indexOfBases[j].innerProd(
						intFrame, c - left, c - left + boundbox.width, 
						r - top, r - top + boundbox.height);
				ssd += nbsNorm;
				if (ssd < globalMinSSD) {
					globalMinSSD = ssd;
				}
				Nd = Nd - 1;
				sampleDistractors.push_back(topleft);
			}
		}
		cvReleaseMat(&intFrame);
		cvReleaseMat(&int2Frame);
		static const double threshDistractSSD = 0; //2000;
		globalMinSSD = globalMinSSD / boundbox.width / boundbox.height;
		printf("selfSSD = %lf, globalMinSSD = %lf\n", selfSSD, globalMinSSD);
		fprintf(fstatSSD, "%d\t%lf\n", current_frame, globalMinSSD);
		if (newly_updated) {
			newly_updated = false;
			markSSD = globalMinSSD;
		}
		static const double urgentRatioSSD = 0.5;
		static const double constrastRatioSSD = 0.8;
		if (/*reconstruct_error > 40
			||*/ globalMinSSD < urgentRatioSSD * markSSD
			/*|| selfSSD > constrastRatioSSD * globalMinSSD*/) {
#endif
			printf("%d, %d\n", samp_width, samp_height);
            currtime = clock();
            nbsobj.computeDNBS(
                num_of_base, num_of_foregrounds, num_of_backgrounds, frame,
                boundbox, margin_width + 20, margin_height + 20,
                distance_threshold, lambda, version, samp_width, samp_height, 
				searchratio, previous_foreground);
			timestat_train += (clock() - currtime) / (double)CLOCKS_PER_SEC;
			++ number_of_train;
            printf("time cost: %.6lf\r\n", (clock() - currtime)/(double)CLOCKS_PER_SEC);
            nbsobj.coefMatrix();
            for (int i = 0; i < (int)nbsobj.bgsample.size(); ++ i)
                lyon::drawRect(showframe, nbsobj.bgsample[i], cvScalar(0xff, 0, 0));
            if (first_template != NULL)
                cvReleaseMat(&first_template);
            first_template = cvCloneMat(nbsobj.target);
			newly_updated = true;
			printf("Reconstruct. Error = %lf\n", nbsobj.reconerr());
        }

//  		cvShowImage("result", showframe);
//  		if (cvWaitKey() == 27) break;

// 		if (num_of_foregrounds > 1) {
// 			if (previous_foreground[float_pointer] != NULL)
// 				cvReleaseMat(&previous_foreground[float_pointer]);
// 			previous_foreground[float_pointer] = target;
// 			float_pointer = (float_pointer + 1) % (num_of_foregrounds - 1);
// 		}

        cvReleaseImage(&foreground);
        cvReleaseMat(&reference_template);
        cvReleaseMat(&intReference);
        }

LABEL_NOUPDATE:
		if (g_esc_thread) break;
		if (!losetrack)
			lyon::drawRect(showframe, boundbox, CV_RGB(0xff, 0, 0));
		else lyon::drawRect(showframe, boundbox, CV_RGB(0, 0xff, 0xff));
//		WaitForSingleObject(hMutex, INFINITE);
            cvCopyImage(showframe, g_showframe);
            cvShowImage(videoSetting.title.c_str(), g_showframe);
            cvWaitKey(33);
//		ReleaseMutex(hMutex);

// 		if (losetrack) {
// 			//boundbox = prevBoundBox;
// 			static const int lost_margin_width = 200;
// 			static const int lost_margin_height = 200;
// 
// 			left = max(0, boundbox.x - lost_margin_width);
// 			right = min(frame->width - boundbox.width, boundbox.x + lost_margin_width);
// 			top = max(0, boundbox.y - lost_margin_height);
// 			bottom = min(frame->height - boundbox.height, boundbox.y + lost_margin_height);
// 
// // 			left = 0;
// // 			right = max(0, frame->width - boundbox.width);
// // 			top = 0;
// // 			bottom = max(0, frame->height - boundbox.height);
// 			//CvMat* recon = nbsobj.compProjection(intInitial);
// 			//cv::Mat matrecon = cv::Mat(recon);
// 			//cv::imshow("dnbs", matrecon);
// 			//cvShowImage("dnbs", recon);
// 			//cvWaitKey();
// 			//cvReleaseMat(&recon);;
// 		}

		if (xmldoc != NULL) {
			createFrameXMLNode(rootNode, current_frame, boundbox);
		}

        printf("Frame %d\r\n", current_frame);

        if (pVideoWriter != NULL)
            cvWriteFrame(pVideoWriter, showframe);
		if (savetype == SAVE2IMAGES) {
			char filename[255];
			sprintf(filename, "%s/%04d.jpg", pathSave.c_str(), current_frame);
			cvSaveImage(filename, showframe);
		}

        cvReleaseImage(&showframe);
        cvReleaseImage(&frame);

		timeframe = (clock() - timeframe) / CLOCKS_PER_SEC;
		avetimeframe += timeframe;
		++ avenumframe;

    }

	avetimeframe /= avenumframe;
	fprintf(stderr, "fps: %.2lf\n", 1 / avetimeframe);
	fclose(fstatSSD);

	timestat_train = timestat_train / number_of_train;
	timestat_track = timestat_track / number_of_track;

	if (xmldoc != NULL) {
		char buf[255];
		xmlNewProp(rootNode, BAD_CAST"frame_create", BAD_CAST itoa(videoSetting.frameCreate, buf, 255));
		xmlNewProp(rootNode, BAD_CAST"frame_finish", BAD_CAST itoa(current_frame - 1, buf, 255));
		xmlNewProp(rootNode, BAD_CAST"avetime_train", BAD_CAST double2string(timestat_train, buf));
		xmlNewProp(rootNode, BAD_CAST"avetime_track", BAD_CAST double2string(timestat_track, buf));
		xmlNewProp(rootNode, BAD_CAST"time_preproc", BAD_CAST double2string(timestat_preproc, buf));
		xmlSaveFormatFileEnc(
			videoSetting.path_save_resxml.c_str(), 
			xmldoc, "UTF-8", 1);
		/*free the document */
		xmlFreeDoc(xmldoc);
		xmlCleanupParser();
		xmlMemoryDump();//debug memory for regression tests 
	}

    if (pVideoWriter != NULL)
        cvReleaseVideoWriter(&pVideoWriter);

	if (pVideoCapture != NULL)
		cvReleaseCapture(&pVideoCapture);

	for (int i = 0; i < num_of_foregrounds - 1; ++ i)
		if (previous_foreground[i] != NULL)
			cvReleaseMat(&previous_foreground[i]);
	free(previous_foreground);

	printf("Finish.\r\n");

	g_esc_thread = true;

    return 1;
}
