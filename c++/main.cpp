#include "opencv2/core.hpp"
#include "opencv2/imgproc.hpp"
#include "opencv2/highgui.hpp"
#include "opencv2/videoio.hpp"
#include "LineDetect.h"

#include <iostream>
#include <vector>

using namespace cv;
using namespace std;

void DrawPoints(vector<int> &corners, cv::Mat& img, string title)
{
	cv::Mat dest = img.clone();

	for (size_t i = 0; i < corners.size() / 2; i++)
	{
		CvScalar color = cvScalar(0, 0, 255);
		circle(dest, cv::Point(corners[2 * i + 0], corners[2 * i + 1]), 4, color, 3);
	}

	cv::namedWindow(title, CV_WINDOW_NORMAL);
	cv::imshow(title, dest);
	cv::waitKey(0);
}

int main(int argc, char **argv)
{
    //string imgpath = "./2.jpg";
    string imgpath = argv[1];

	//read
	Mat img = imread(imgpath, 0);

    int corners[8];
    int result = getVertex(img.data, img.cols, img.rows, img.cols, corners);
	if(result < 0)
		cout << "detect vertex error" << endl;

	//show the result
	Mat img_rgb = imread(imgpath, 1);
	vector<int> vertexs(corners, corners + 8);
	DrawPoints(vertexs, img_rgb, "final result");	
}
