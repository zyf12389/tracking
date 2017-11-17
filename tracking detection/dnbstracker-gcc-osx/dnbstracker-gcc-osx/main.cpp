// dnbstracker_win32.cpp 
// 
// This source file is the main entrance of this application.
//
// It implements a robust visual object tracker based on
// the discriminative non-orthogonal binary subspace (D-NBS).
//
// This is a win32 console application, and
// there are 3 ways to import the data and parameters
// and to start the D-NBS tracker:
// [1] Passing arguments through the command line.
//		Guidance will be released soon.
// [2] Starting a GUI to input paths and parameters.
//		%% -gui
// [3] (Preferred) Passing arguments by formatted XML.
//		%% -xml [PathOfXMLfile]
//
// Environment:
// This application is developped under Microsoft Visual
// Studio 2008, using MFC, Intel OpenCV 2.0, and libxml.
// To successfully execute the application, following are
// necessary requirements:
// [1] Windows Series System
// [2] Visual Studio 2008 / Redistributable Package or above
// [3] Intel OpenCV 2.0 (dll provided)
// [4] libxml (dll provided)
// 
// Tested successful under Windows 7 Ultimate Ed., Nov.8, 2010.
//
// Reference:
// Ang Li, Feng Tang, Yanwen Guo, and Hai Tao. Discriminative
// Nonorthogonal Binary Subspace Tracking. In: Proceedings of
// the 11th European Conference on Computer Vision (ECCV), 2010.
//
// Ang Li, (mailto:nju.angli@gmail.com)
// Dept. Computer Science and Technology
// Nanjing University
// Nov. 8, 2010.

//using namespace std;


#include "headers.h"
#include "commonlib.h"
#include "DnbsTracker.h"

typedef char TCHAR;

struct VideoStruct;
struct DnbsParamStruct;

extern int startDnbsTracker(VideoStruct & videoStruct,
							DnbsParamStruct & dnbsStruct);

extern int saveInitialXML(const string& path,
						  const VideoStruct & videoSetting, 
						  const DnbsParamStruct & dnbsParam);

extern int readInitialXML(const string& path,
						  VideoStruct & videoSetting, 
						  DnbsParamStruct & dnbsParam);

int AnalyseCommandLine(int argc, TCHAR* argv[]);

vector<string> split(const string& s, char c);

int printSettings(const VideoStruct & videoSetting, 
				  const DnbsParamStruct & dnbsParam);

void showGuidance();

int main(int argc, TCHAR* argv[])
{
	int nRetCode = 0;
    
    if (argc == 1) {
        showGuidance();
        nRetCode = -1;
    } else
        nRetCode = AnalyseCommandLine(argc, argv);
    
	return nRetCode;
}

int AnalyseCommandLine(int argc, TCHAR* argv[]) {
	int loadtype = 0; // 0: cmd; 1: xml; 2: gui
	string path_inixml; // path for init xml file
	bool save_inixml = false; // save init xml
	string path_save_inixml; // path for init xml to save

	string strParam;
	VideoStruct videoSetting;
	DnbsParamStruct dnbsParam;

	for (int p = 1; p < argc; ++ p) {
		if (strlen(argv[p]) > 0 && argv[p][0] == '-') {
			// new argument
			if (!strcmp(argv[p], "-xml")) {
				loadtype = 1;
				path_inixml = argv[++p];
			} else if (!strcmp(argv[p], "-gui")) {
				loadtype = 2;
			} else if (!strcmp(argv[p], "-video")) {
				videoSetting.srctype = 0;
				videoSetting.videoPath = argv[++p];
			} else if (!strcmp(argv[p], "-image")) {
				videoSetting.srctype = 1;
				videoSetting.imagePattern = argv[++p];
			} else if (!strcmp(argv[p], "-from")) {
				videoSetting.frameCreate = atoi(argv[++p]);
			} else if (!strcmp(argv[p], "-to")) {
				videoSetting.frameFinish = atoi(argv[++p]);
			} else if (!strcmp(argv[p], "-title")) {
				videoSetting.title = argv[++p];
			} else if (!strcmp(argv[p], "-oinixml")) {
				save_inixml = true;
				path_save_inixml = argv[++p];
			} else if (!strcmp(argv[p], "-oxml")) {
				videoSetting.save_resxml = true;
				videoSetting.path_save_resxml = argv[++p];
			} else if (!strcmp(argv[p], "-oavi")) {
				videoSetting.save_resavi = true;
				videoSetting.path_save_resavi = argv[++p];
			} else if (!strcmp(argv[p], "-set")) {
				strParam = argv[++p];
			} else if (!strcmp(argv[p], "-box")) {
				string boxParam = argv[++p];
				vector<string> vBoxParams
					= split(boxParam, ';');
				for (int i = 0; i < (int)vBoxParams.size(); ++ i) {
					vector<string> locParams
						= split(vBoxParams[i], '=');
					if (locParams.size() != 2)
						perror("Parameter setting errors.");
					if (locParams[0] == "x") {
						videoSetting.boundbox.x = atoi(locParams[1].c_str());
					} else if (locParams[0] == "y") {
						videoSetting.boundbox.y = atoi(locParams[1].c_str());
					} else if (locParams[0] == "width") {
						videoSetting.boundbox.width = atoi(locParams[1].c_str());
					} else if (locParams[0] == "height") {
						videoSetting.boundbox.height = atoi(locParams[1].c_str());
					} else perror("Unknown parameter(s).");
				}
			} else if (!strcmp(argv[p], "-margin")) {
				string marParam = argv[++p];
				vector<string> vMarParams
					= split(marParam, ';');
				for (int i = 0; i < (int)vMarParams.size(); ++ i) {
					vector<string> locParams
						= split(vMarParams[i], '=');
					if (locParams.size() != 2)
						perror("Parameter setting errors.");
					if (locParams[0] == "width") {
						videoSetting.marginWidth = atoi(locParams[1].c_str());
					} else if (locParams[0] == "height") {
						videoSetting.marginHeight = atoi(locParams[1].c_str());
					} else perror("Unknown parameter(s).");
				}
			}
		}
	}
	// Arguments resolved.
	//uisetting dlg_uiset(&videoSetting, &dnbsParam);

	switch (loadtype) {
		case 2: 
			// load from graphical user interfaces
			//dlg_uiset.DoModal();
			break;

		case 1:
			// load from xml initial file
			readInitialXML(path_inixml, videoSetting, dnbsParam);

		default:
			// load from command arguments
			vector<string> vParams
				= split(strParam, ';');
			for (int i = 0; i < (int)vParams.size(); ++ i) {
				vector<string> locParams
					= split(vParams[i], '=');
				if (locParams.size() != 2)
					perror("Parameter setting errors.");
				if (locParams[0] == "nbase") {
					dnbsParam.numBases = atoi(locParams[1].c_str());
				} else if (locParams[0] == "nfgrd") {
					dnbsParam.numForegrounds = atoi(locParams[1].c_str());
				} else if (locParams[0] == "nbgrd") {
					dnbsParam.numBackgrounds = atoi(locParams[1].c_str());
				} else if (locParams[0] == "vers") {
					dnbsParam.version = atoi(locParams[1].c_str());
				} else if (locParams[0] == "lambda") {
					dnbsParam.lambda = atof(locParams[1].c_str());
				} else if (locParams[0] == "samp_width") {
					dnbsParam.samp_width = atoi(locParams[1].c_str());
				} else if (locParams[0] == "samp_height") {
					dnbsParam.samp_height = atoi(locParams[1].c_str());
				} else if (locParams[0] == "iprod") {
					dnbsParam.min_innprod = atof(locParams[1].c_str());
				} else if (locParams[0] == "ratio") {
					dnbsParam.searchratio = atof(locParams[1].c_str());
				} else perror("Unknown parameter(s).");
			}
	}

	if (save_inixml) {
		saveInitialXML(path_save_inixml, videoSetting, dnbsParam);
	}
	printSettings(videoSetting, dnbsParam);

    /*
	if (dlg_uiset.cancel)
		goto RETURN;*/

	// Start the D-NBS Tracker
	startDnbsTracker(videoSetting, dnbsParam);

//RETURN:
	return 0;
}
vector<string> split(const string& s, char c) {
	vector<string> vRet;
	string::size_type i = 0;
	string::size_type j = s.find(c);
	while (j != string::npos) {
		vRet.push_back(s.substr(i, j-i));
		i = ++j;
		j = s.find(c, j);
		if (j == string::npos)
			vRet.push_back(s.substr(i, s.length()));
	}
	return vRet;
}

int printSettings(const VideoStruct & videoSetting, const DnbsParamStruct & dnbsParam) {
	// Show all the settings on the console screen.
	printf("[TRACKER: VIDEO SETTING]\n");
	printf("Title: %s\n", videoSetting.title.c_str());
	printf("Bounding Box: x=%d,y=%d,width=%d,height=%d\n", 
		videoSetting.boundbox.x, videoSetting.boundbox.y,
		videoSetting.boundbox.width, videoSetting.boundbox.height);
	printf("FrameCreate: %d\nFrameFinish: %d\n", 
		videoSetting.frameCreate, videoSetting.frameFinish);
	printf("Margin: width=%d,height=%d\n",
		videoSetting.marginWidth, videoSetting.marginHeight);
	printf("Type: %s\n", videoSetting.srctype ? "Image Sequence" : "Video File");
	if (videoSetting.srctype == 0)
		printf("Video Path: %s\n", videoSetting.videoPath.c_str());
	else printf("Image Pattern: %s\n", videoSetting.imagePattern.c_str());
	if (videoSetting.save_resxml)
		printf("Result XML path: %s\n", videoSetting.path_save_resxml.c_str());
	if (videoSetting.save_resavi)
		printf("Result AVI path: %s\n", videoSetting.path_save_resavi.c_str());

	// DNBS Parameters
	printf("\n[TRACKER: DNBS PARAMETER]\n");
	printf("lambda = %.3lf\n", dnbsParam.lambda);
	printf("nbase = %d\n", dnbsParam.numBases);
	printf("nfgrd = %d\n", dnbsParam.numForegrounds);
	printf("nbgrd = %d\n", dnbsParam.numBackgrounds);
	printf("vers = %d\n", dnbsParam.version);

	return 0;
}

#include <iomanip>
void showGuidance() {
	int width = 10;
	cout << "Discriminative Nonorthogonal Binary Subspace Tracker" << endl
		<< endl 
		<< "COMMAND LINE GUIDANCE" << endl
		<< "Parameters" << endl
		<< setw(width) << "-xml"	<< "\t" << "Initial the tracker by XML file" << endl
		<< setw(width) << " "		<< "\t"	<< "format:  -xml <xmlpath>" << endl<< endl
		<< setw(width) << "-gui"	<< "\t" << "Initial the tracker through GUI" << endl
		<< setw(width) << " "		<< "\t"	<< "format:  -gui" << endl<< endl
		<< setw(width) << "-video"	<< "\t" << "Make tracker take in video file" << endl
		<< setw(width) << " "		<< "\t"	<< "format:  -video <videopath>" << endl<< endl
		<< setw(width) << "-image"	<< "\t" << "Make tracker take in image sequence" << endl
		<< setw(width) << " "		<< "\t"	<< "format:  -image <imageformat>" << endl<< endl
		<< setw(width) << "-title"	<< "\t" << "Set the title for this sequence" << endl
		<< setw(width) << " "		<< "\t"	<< "format:  -title <titlename=\"Unknown\">" << endl<< endl
		<< setw(width) << "-from"	<< "\t" << "Set the frame it creates from" << endl
		<< setw(width) << " "		<< "\t"	<< "format:  -from <int_from=0>" << endl<< endl
		<< setw(width) << "-to"		<< "\t" << "Set the frame it finishes at" << endl
		<< setw(width) << " "		<< "\t"	<< "format:  -to <int_to=-1>" << endl<< endl
		<< setw(width) << "-oinixml"<< "\t" << "Save all the parameters to XML initial file " << endl
		<< setw(width) << " "		<< "\t"	<< "for future use - quickly restarting the tracker." << endl
		<< setw(width) << " "		<< "\t"	<< "format:  -oinixml <oini_xmlpath>" << endl<< endl
		<< setw(width) << "-oxml"	<< "\t" << "Save the tracking results to an XML" << endl
		<< setw(width) << " "		<< "\t"	<< "format:  -oxml <o_xmlpath>" << endl<< endl
		<< setw(width) << "-oavi"	<< "\t" << "Save the tracking result video to an AVI" << endl
		<< setw(width) << " "		<< "\t"	<< "format:  -oavi <avipath>" << endl<< endl
		<< setw(width) << "-set"	<< "\t" << "Set parameters for D-NBS" << endl
		//<< setw(width) << " "		<< "\t" << "   Arguments:" << endl
		<< setw(width) << " "		<< "\t" << "   1. nbase - number of bases (default: 30)" << endl
		<< setw(width) << " "		<< "\t" << "   2. nfgrd - number of foreground samples (default: 1)" << endl
		<< setw(width) << " "		<< "\t" << "   3. nbgrd - number of background samples (default: 3)" << endl
		<< setw(width) << " "		<< "\t" << "   4. lambda - parameter in D-NBS formulation (default: 0.25)" << endl
		<< setw(width) << " "		<< "\t" << "   5. vers  - 0 or 1, version of implementation (default: 0)" << endl
		<< setw(width) << " "		<< "\t"	<< "format:  -set <param1>=<num1>;...;<paramN>=<numN>" << endl<< endl
		<< setw(width) << "-margin"	<< "\t" << "Set width and height of searching margin" << endl
		<< setw(width) << " "		<< "\t"	<< "format:  -margin width=<num1>;height=<num2>" << endl<< endl
		<< setw(width) << "-box"	<< "\t" << "Initial the bounding box at the first frame" << endl
		<< setw(width) << " "		<< "\t"	<< "format:  -box x=<num1>;y=<num2>;width=<num3>;height=<num4>" << endl;

	system("pause");
}