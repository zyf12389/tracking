#pragma once
#include <cv.h>
#include "cxcore.h"
#include "highgui.h"
#include <libxml/parser.h>

#include <iostream>
#include <vector>
#include <string>
#include <queue>
using namespace std;

// VideoStruct stores the information for sequences
struct VideoStruct {
	int srctype; // 0: video; 1: image sequences
	string title;
	CvRect boundbox;
	int frameCreate;
	int frameFinish;
	int frameCurrent;
	int marginWidth;
	int marginHeight;

	string videoPath;
	string imagePattern;

	bool save_resxml; // save result xml
	string path_save_resxml; // path for result xml
	bool save_resavi; // save result avi
	string path_save_resavi; // path for result avi

	VideoStruct() {
		srctype = 0;
		title = "Unknown";
		videoPath = "";
		imagePattern = "";
		frameCreate = -1;
		frameFinish = -1;
		frameCurrent = -1;
		save_resxml = false;
		save_resavi = false;
		marginHeight = 20;
		marginWidth = 40;
	}
};

// DnbsParamStruct stores the information for D-NBS parameters
struct DnbsParamStruct {
	double lambda;
	int numBases;
	int numForegrounds;
	int numBackgrounds;
	int version;
	int samp_width;
	int samp_height;
	double min_innprod;
	double searchratio;

	DnbsParamStruct() {
		lambda = 0.25;
		numBases = 30;
		numForegrounds = 1;
		numBackgrounds = 3;
		version = 0;
		samp_width = 2;
		samp_height = 2;
		min_innprod = 0;
		searchratio = 1;
	}
};
