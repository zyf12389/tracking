// Useful tools, including:
//	(1) string-number transformation;
//	(2) Reading and writing XML initial files
//		for D-NBS Trackers.
//
// Ang Li (mailto:nju.angli@gmail.com)
// Nov.8, 2010.

#pragma once
#include "headers.h"

char* itoa(int num, char* buf, int size);
char* double2string(double num, char* buf);
int readInitialXML(const string& path, VideoStruct & videoSetting, DnbsParamStruct & dnbsParam);
int saveInitialXML(const string& path, const VideoStruct & videoSetting, const DnbsParamStruct & dnbsParam);
