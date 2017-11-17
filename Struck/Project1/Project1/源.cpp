#include<iostream>
#include<fstream>
#include<opencv2/opencv.hpp>
#include<opencv2/highgui/highgui.hpp>
using namespace std;
using namespace cv;
int main()
{
	string imgpath = "C://Users//ww//Desktop//visual tracking//KCF//tracker_release2//tracker_release2//sequence//dog1//img";
	string filepath = "C://Users//ww//Desktop//visual tracking//muster//MUSTer_code_v1.1//MUSTer_code_v1.1//Results//dog1.txt";
	string imgformat = imgpath + "/img%05d.jpg";
	string fileframe = "C://Users//ww//Desktop//visual tracking//KCF//tracker_release2//tracker_release2//sequence//dog1//dog1_frames.txt";
	string framesline;
	string pos;
	ifstream in(fileframe);
	ifstream inFile(filepath);
	getline(in, framesline);
	int lx, ly, width, height;
	getline(inFile,pos);
	int startframe, endframe;
	sscanf(framesline.c_str(), "%d,%d", &startframe, &endframe);
	sscanf(pos.c_str(),"%d,%d,%d,%d",&lx,&ly,&width,&height);
	char filename[256];
	
	VideoCapture capture("VideoTest.avi");
	VideoWriter writer("VideoTest.avi", CV_FOURCC('X', 'V', 'I', 'D'), 25.0, Size(640, 480));


	for (int i = startframe; i <= endframe; i++)
	{
		sprintf(filename, imgformat.c_str(), i);
		cout << filename << endl;
		Mat img = imread(filename);
		if (!img.data)
			cout << "no image\n";
		Mat gray = Mat::zeros(img.rows, img.cols, CV_8UC1);
		cvtColor(img, gray, CV_RGB2GRAY);
	
		rectangle(gray,Point(lx,ly),Point(lx+width-1,ly+height-1),Scalar(255,255,255));
		imshow("video", gray);
		waitKey(27);
		writer << gray;
		
		getline(inFile,pos);
		sscanf(pos.c_str(), "%d,%d,%d,%d", &lx, &ly, &width, &height);
	}
	Mat frame;
	while (1)
	{
		capture >> frame;
		imshow("ds", frame);
		waitKey(0);
	}
	waitKey(-1);
	in.close();
	inFile.close();
	system("pause");
	return 0;
}