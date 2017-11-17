// Useful tools, including:
//	(1) string-number transformation;
//	(2) Reading and writing XML initial files
//		for D-NBS Trackers.
//
// Ang Li (mailto:nju.angli@gmail.com)
// Nov.8, 2010.

#include "commonlib.h"

char* itoa(int num, char* buf, int size) {
	memset(buf, 0, sizeof(char) * size);
	sprintf(buf, "%d", num);
	return buf;
}

char* double2string(double num, char* buf) {
	sprintf(buf, "%.6lf", num);
	return buf;
}

int readInitialXML(const string& path, VideoStruct & videoSetting, DnbsParamStruct & dnbsParam) {

	xmlDocPtr doc
		= xmlReadFile(path.c_str(), "UTF-8", XML_PARSE_RECOVER);
	if (NULL == doc) {
		printf("Unable to open the XML file.\n");
		return -1;
	}

	xmlNodePtr rootNode = xmlDocGetRootElement(doc);
	if (NULL == rootNode) {
		printf("The XML file is empty.\n");
		xmlFreeDoc(doc);
		return -1;
	}

	if (xmlStrcmp(rootNode->name, BAD_CAST "init")) {
		printf("XML format error.\n");
		return -1;
	}

	xmlNodePtr node = NULL;

	// Looking for 'sequence'
	node = rootNode->xmlChildrenNode; 
	while (node != NULL && xmlStrcmp(node->name, BAD_CAST"sequence"))
		node = node->next;
	if (node == NULL)
		perror("XML format error");

	videoSetting.title = (const char*)xmlGetProp(node, BAD_CAST"title");
	string srctype = (char*)xmlGetProp(node, BAD_CAST"type");
	if (srctype == "video")
		videoSetting.srctype = 0;
	else videoSetting.srctype = 1;
	videoSetting.frameCreate = atoi((char*)xmlGetProp(node, BAD_CAST"frame_create"));
	videoSetting.frameFinish = atoi((char*)xmlGetProp(node, BAD_CAST"frame_finish"));

	node = node->xmlChildrenNode;
	while (node != NULL) {
		if (!xmlStrcmp(node->name, BAD_CAST"videopath")) {
			videoSetting.videoPath = (char*)xmlNodeGetContent(node);
		} else if (!xmlStrcmp(node->name, BAD_CAST"imagepattern")) {
			videoSetting.imagePattern = (char*)xmlNodeGetContent(node);
		} else if (!xmlStrcmp(node->name, BAD_CAST"margin")) {
			videoSetting.marginWidth = atoi((char*)xmlGetProp(node, BAD_CAST"width"));
			videoSetting.marginHeight = atoi((char*)xmlGetProp(node, BAD_CAST"height"));
		} else if (!xmlStrcmp(node->name, BAD_CAST"boundbox")) {
			videoSetting.boundbox.x = atoi((char*)xmlGetProp(node, BAD_CAST"x"));
			videoSetting.boundbox.y = atoi((char*)xmlGetProp(node, BAD_CAST"y"));
			videoSetting.boundbox.width = atoi((char*)xmlGetProp(node, BAD_CAST"width"));
			videoSetting.boundbox.height = atoi((char*)xmlGetProp(node, BAD_CAST"height"));
		}
		node = node->next;
	}

	node = rootNode->xmlChildrenNode;
	while (node != NULL && xmlStrcmp(node->name, BAD_CAST"dnbsparam"))
		node = node->next;
	if (node == NULL)
		perror("XML format error");

	dnbsParam.lambda = atof((char*)xmlGetProp(node, BAD_CAST"lambda"));
	dnbsParam.numBases = atoi((char*)xmlGetProp(node, BAD_CAST"nbase"));
	dnbsParam.numForegrounds = atoi((char*)xmlGetProp(node, BAD_CAST"nfgrd"));
	dnbsParam.numBackgrounds = atoi((char*)xmlGetProp(node, BAD_CAST"nbgrd"));
	dnbsParam.version = atoi((char*)xmlGetProp(node, BAD_CAST"vers"));
	dnbsParam.samp_height = atoi((char*)xmlGetProp(node, BAD_CAST"samp_height"));
	dnbsParam.samp_width = atoi((char*)xmlGetProp(node, BAD_CAST"samp_width"));
	if (xmlHasProp(node, BAD_CAST"iprod") != NULL)
	dnbsParam.min_innprod = atof((char*)xmlGetProp(node, BAD_CAST"iprod"));
	if (xmlHasProp(node, BAD_CAST"ratio") != NULL)
	dnbsParam.searchratio = atof((char*)xmlGetProp(node, BAD_CAST"ratio"));

	return 0;
}


int saveInitialXML(const string& path, const VideoStruct & videoSetting, const DnbsParamStruct & dnbsParam) {
	// Temporary char buffer
	char buf[255];

	// Creates a new document, a node and set it as a root node
	xmlDocPtr doc 
		= xmlNewDoc(BAD_CAST "1.0");
	xmlNodePtr rootNode
		= xmlNewNode(NULL, BAD_CAST "init");
	xmlDocSetRootElement(doc, rootNode);

	//Create XML node for video settings
	xmlNodePtr sequenceNode
		= xmlNewNode(NULL, BAD_CAST "sequence");
	xmlAddChild(rootNode, sequenceNode);
	xmlNewProp(sequenceNode, BAD_CAST "title", BAD_CAST videoSetting.title.c_str());
	xmlNewProp(sequenceNode, BAD_CAST "type", videoSetting.srctype == 0 ? BAD_CAST "video" : BAD_CAST "image");
	if (videoSetting.srctype == 0)
		xmlNewChild(sequenceNode, NULL, BAD_CAST "videopath", BAD_CAST videoSetting.videoPath.c_str());
	else xmlNewChild(sequenceNode, NULL, BAD_CAST "imagepattern", BAD_CAST videoSetting.imagePattern.c_str());
	xmlNewProp(sequenceNode, BAD_CAST "frame_create", BAD_CAST itoa(videoSetting.frameCreate, buf, 255));
	xmlNewProp(sequenceNode, BAD_CAST "frame_finish", BAD_CAST itoa(videoSetting.frameFinish, buf, 255));

	xmlNodePtr marginNode
		= xmlNewNode(NULL, BAD_CAST "margin");
	xmlNewProp(marginNode, BAD_CAST "width", BAD_CAST itoa(videoSetting.marginWidth, buf, 255));
	xmlNewProp(marginNode, BAD_CAST "height", BAD_CAST itoa(videoSetting.marginHeight, buf, 255));
	xmlAddChild(sequenceNode, marginNode);

	xmlNodePtr boxNode
		= xmlNewNode(NULL, BAD_CAST "boundbox");
	xmlNewProp(boxNode, BAD_CAST "x", BAD_CAST itoa(videoSetting.boundbox.x, buf, 255));
	xmlNewProp(boxNode, BAD_CAST "y", BAD_CAST itoa(videoSetting.boundbox.y, buf, 255));
	xmlNewProp(boxNode, BAD_CAST "width", BAD_CAST itoa(videoSetting.boundbox.width, buf, 255));
	xmlNewProp(boxNode, BAD_CAST "height", BAD_CAST itoa(videoSetting.boundbox.height, buf, 255));
	xmlAddChild(sequenceNode, boxNode);

	// Create XML Node for DNBS Parameters
	xmlNodePtr dnbsNode
		= xmlNewNode(NULL, BAD_CAST "dnbsparam");
	xmlAddChild(rootNode, dnbsNode);
	xmlNewProp(dnbsNode, BAD_CAST "lambda", BAD_CAST double2string(dnbsParam.lambda, buf));
	xmlNewProp(dnbsNode, BAD_CAST "nbase", BAD_CAST itoa(dnbsParam.numBases, buf, 255));
	xmlNewProp(dnbsNode, BAD_CAST "nfgrd", BAD_CAST itoa(dnbsParam.numForegrounds, buf, 255));
	xmlNewProp(dnbsNode, BAD_CAST "nbgrd", BAD_CAST itoa(dnbsParam.numBackgrounds, buf, 255));
	xmlNewProp(dnbsNode, BAD_CAST "vers", BAD_CAST itoa(dnbsParam.version, buf, 255));
	xmlNewProp(dnbsNode, BAD_CAST"samp_width", BAD_CAST itoa(dnbsParam.samp_width, buf, 255));
	xmlNewProp(dnbsNode, BAD_CAST"samp_height", BAD_CAST itoa(dnbsParam.samp_height, buf, 255));
	xmlNewProp(dnbsNode, BAD_CAST"iprod", BAD_CAST double2string(dnbsParam.min_innprod, buf));
	xmlNewProp(dnbsNode, BAD_CAST"ratio", BAD_CAST double2string(dnbsParam.searchratio, buf));

	//Dumping document to stdio or file
	xmlSaveFormatFileEnc(path.c_str(), doc, "UTF-8", 1);
	/*free the document */
	xmlFreeDoc(doc);
	xmlCleanupParser();
	xmlMemoryDump();//debug memory for regression tests 

	return 0;
}

