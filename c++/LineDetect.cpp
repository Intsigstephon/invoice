#include "LineDetect.h"
#include "intsigImage.h"
#include <algorithm>
#include <iostream>

using itcv::LINE;
using namespace std;

#define DEBUG 0
#if DEBUG
    #include "opencv2/core.hpp"
    #include "opencv2/imgproc.hpp"
    #include "opencv2/highgui.hpp"
    #include "opencv2/videoio.hpp"
    using namespace cv;
#endif

//assume invoice has been resized
const int min_hori_len = 500;    //min horizontal line length
const int min_vert_len = 300;    //min vertical line length
const int valley_thresh = 30;    //valley thresh

const int max_line_gap =  10;    //parallel line dist
const int max_merge_cost = 1000; //max merge cost
 
const int hori_gap_vec[4] = {238, 440, 85, 220};  //hori gap vec
const int vert_gap_vec[4] = {522, 1474};          //vert gap vec

const double vert_ratio = 0.53;     //filter vert line
const double match_thresh = 0.35;
const double dist_ratio = 2.11;     //vert line dist / hori line dist
const double vert_dist_diff = 0.1;  //vert line dist diff

static double GetDistance(itcv::Point &pnt1, itcv::Point &pnt2)
{
	int dx = pnt1.x - pnt2.x;
	int dy = pnt1.y - pnt2.y;
	return sqrt(double(dx)*dx + dy*dy);
}

static double GetDistance(itcv::Point &pnt, itcv::Point &StPnt, itcv::Point &EdPnt)
{
	double a, b;
	if (EdPnt.x == StPnt.x)
		return (int)fabs(double(pnt.x) - StPnt.x);
	if (EdPnt.y == StPnt.y)
		return (int)fabs(double(pnt.y) - StPnt.y);
	a = (EdPnt.y - StPnt.y)*1.0 / (EdPnt.x - StPnt.x);
	b = StPnt.y - a*StPnt.x;
	itcv::Point   pnt2;
	double  x2, y2;
	x2 = (pnt.x + a*pnt.y - a*b)*1.0 / (a*a + 1);
	y2 = a*x2 + b;
	pnt2.x = (long)x2;
	pnt2.y = (long)y2;
	return GetDistance(pnt, pnt2);
}

double Parallel_Line_dist(LINE &a, LINE &b)
{
	double dist1 = GetDistance(a.StPnt, b.StPnt, b.EdPnt);
	double dist2 = GetDistance(a.EdPnt, b.StPnt, b.EdPnt);
	//return (dist1 + dist2) / 2;
	return min(dist1, dist2);
}

static bool	SameLine(LINE &a, LINE &b, int thresh)
{
	double dist1 = GetDistance(a.StPnt, b.StPnt, b.EdPnt);
	double dist2 = GetDistance(a.EdPnt, b.StPnt, b.EdPnt);
	//if (dist1 + dist2 < 2 * thresh)
	if (min(dist1, dist2) < 2 * thresh)
	{
		return true;
	}
	return false;
}

static double DotProduct(double x1, double y1, double x2, double y2)
{
	return x1 * x2 + y1 * y2;
}

static bool SameDir(LINE &a, LINE &b)
{
	if (DotProduct(a.EdPnt.x - a.StPnt.x, a.EdPnt.y - a.StPnt.y, b.EdPnt.x - b.StPnt.x, b.EdPnt.y - b.StPnt.y) >= 0.99 * a.distance() * b.distance())
		return 1;
	return 0;
}

static bool OppoSiteDir(LINE &a, LINE &b)
{
	if (DotProduct(a.EdPnt.x - a.StPnt.x, a.EdPnt.y - a.StPnt.y, b.EdPnt.x - b.StPnt.x, b.EdPnt.y - b.StPnt.y) <= -0.99 * a.distance() * b.distance())
		return 1;
	return 0;
}

static bool EntireOverlap(LINE &a, LINE &b)
{
	if ((a.StPnt.x - b.StPnt.x) * (a.EdPnt.x - b.EdPnt.x) + (a.StPnt.y - b.StPnt.y) * (a.EdPnt.y - b.EdPnt.y) < 0)
		return 1;
	return 0;
}

static bool	lessthan_y(const LINE &a, const LINE& b)
{
	if (a.StPnt.y + a.EdPnt.y < b.StPnt.y + b.EdPnt.y)
	{
		return true;
	}
	return false;
}

static bool	morethan_y(const LINE &a, const LINE& b)
{
	if (a.StPnt.y + a.EdPnt.y > b.StPnt.y + b.EdPnt.y)
	{
		return true;
	}
	return false;
}

static bool	lessthan_x(const LINE &a, const LINE& b)
{
	if (a.StPnt.x + a.EdPnt.x < b.StPnt.x + b.EdPnt.x)
	{
		return true;
	}
	return false;
}

static bool	morethan_x(const LINE &a, const LINE& b)
{
	if (a.StPnt.x + a.EdPnt.x > b.StPnt.x + b.EdPnt.x)
	{
		return true;
	}
	return false;
}

static bool	greatthan_dist(LINE &a, LINE& b)
{
	if (a.distance() > b.distance())
	{
		return true;
	}
	return false;
}

static int GetCrossPoint(double x1, double y1, double x2, double y2,
	double _x1, double _y1, double _x2, double _y2, itcv::Point &CrossPnt)
{
	double dx, _dx, dy, _dy, x, y;
	dx = x2 - x1; dy = y2 - y1; _dx = _x2 - _x1; _dy = _y2 - _y1;
	if (_dx*dy == dx*_dy)	return -1;
	if (dx == 0)
	{
		x = x1;
		y = _y1 + (x1 - _x1)*_dy / _dx;
	}
	else
	{
		x = (dx*_dx*(_y1 - y1) - dx*_dy*_x1 + _dx*dy*x1) / (_dx*dy - dx*_dy);
		y = (dy*(x - x1) + y1*dx) / dx;
	}
	CrossPnt.x = (int)(x + 0.5); CrossPnt.y = (int)(y + 0.5);
	return 0;
}

void MergeLine(vector<LINE> &vec_line, int gap_thresh, bool Ishorizontal, double merge_thresh)
{
	vector<LINE> diff_line;

	//sort line
	if (Ishorizontal)
	{
		sort(vec_line.begin(), vec_line.end(), lessthan_y);
	}
	else
	{
		sort(vec_line.begin(), vec_line.end(), lessthan_x);
	}

	int oldsize = 0;
	do {
		oldsize = diff_line.size();

		for (size_t i = 0; i < vec_line.size(); i++)
		{
			if (!vec_line[i].flag)
			{
				continue;
			}

			for (size_t j = i + 1; j < vec_line.size(); j++)
			{
				if (!vec_line[j].flag)
				{
					continue;
				}

				double para_dist = Parallel_Line_dist(vec_line[i], vec_line[j]);
				if (para_dist > gap_thresh)
				{
					continue;
				}

				double merge_cost = vec_line[i].merge_cost(vec_line[j]);
				if (merge_cost > merge_thresh)
				{
					continue;
				}

				//judge if overlap
				if (EntireOverlap(vec_line[i], vec_line[j]))
				{
					continue;
				}

				vec_line[i].merge(vec_line[j]);
				vec_line[j].flag = false;
			}
		}

		diff_line.resize(0);
		for (size_t i = 0; i < vec_line.size(); i++)
		{
			if (vec_line[i].flag)
			{
				diff_line.push_back(vec_line[i]);
			}
		}
	} while (diff_line.size() != oldsize);

	vec_line.resize(diff_line.size());
	copy(diff_line.begin(), diff_line.end(), vec_line.begin());

}

void FilterLineBasedLen(vector<LINE> &vec_line, double ratio)
{
	double max_distance = 0;

	for (size_t i = 0; i < vec_line.size(); i++)
	{
		double dist = vec_line[i].distance();
		if (dist > max_distance)
		{
			max_distance = dist;
		}
	}

	vector<LINE> temp;
	for (size_t i = 0; i < vec_line.size(); i++)
	{
		double dist = vec_line[i].distance();
		if (dist > max_distance * ratio)
		{
			temp.push_back(vec_line[i]);
		}
	}

	vec_line.resize(temp.size());
	copy(temp.begin(), temp.end(), vec_line.begin());
}

void UniqueLine(vector<LINE> &vec_line, int gap_thresh)
{
	vector<LINE> temp;
	for (size_t i = 0; i < vec_line.size(); i++)
	{
		if (!vec_line[i].flag) continue;
		temp.push_back(vec_line[i]);
		for (size_t j = i + 1; j < vec_line.size(); j++)
		{
			if (vec_line[j].flag && SameLine(vec_line[i], vec_line[j], gap_thresh))
			{
				vec_line[j].flag = false;
			}
		}
	}

	vec_line.resize(temp.size());
	copy(temp.begin(), temp.end(), vec_line.begin());
}

double DecideBorderLine(vector<LINE> &vec_line, vector<double> &vec_gap, vector<LINE>& result)
{
    double min_dist = 100000;
	double min_error = 0;
	int up = 0;
	int down = 0;

	//Get ratio
	vector<double> ratio;
	for (size_t i = 0; i < vec_gap.size() - 1; i++)
	{
		ratio.push_back((double)vec_gap[i + 1] / (double)vec_gap[i]);
	}

	int num_line = ratio.size() + 2;
	for (size_t i = 0; i <= vec_line.size() - num_line; i++)
	{
		vector<LINE> part(num_line);
		copy(vec_line.begin() + i, vec_line.begin() + num_line + i, part.begin());

		vector<double> part_dist;
		vector<double> part_ratio;
		for (size_t j = 0; j < part.size() - 1; j++)
		{
			double dist1 = GetDistance(part[j + 1].StPnt, part[j].StPnt, part[j].EdPnt);
			double dist2 = GetDistance(part[j + 1].EdPnt, part[j].StPnt, part[j].EdPnt);
			part_dist.push_back(dist1 / 2 + dist2 / 2);
		}

		for (size_t j = 0; j < part_dist.size() - 1; j++)
		{
			part_ratio.push_back((double)part_dist[j + 1] / (double)part_dist[j]);
		}

		double dist = 0;
		double max_error = 0;
		for (size_t j = 0; j < ratio.size(); j++)
		{
			double error = ratio[j] - part_ratio[j];
			dist += error * error;
			max_error = max(max_error, fabs(error) / ratio[j]);
		}

		if (dist < min_dist)
		{
			min_dist = dist;
			min_error = max_error;

			//the line index, also each line index is sure 
			up = i;
			down = i + num_line - 1;
		}
	}

	//result.push_back(vec_line[up]);
	//result.push_back(vec_line[down]);
	for (size_t i = 0; i < num_line; i++)
	{
		result.push_back(vec_line[up + i]);
	}

	return min_error;
}

bool GlobalMatch(vector<LINE> &vec_line, bool IsHorizontal, double match_thresh, vector<double> &vec_gap, vector<LINE> &result)
{
	int num_line = vec_gap.size() + 1;

	if (vec_line.size() < num_line)
	{
        if(DEBUG)
		    std::cout << "error: less line detected" << std::endl;
		return false;
	}

	if (IsHorizontal)
	{
		sort(vec_line.begin(), vec_line.end(), lessthan_y);
	}
	else
	{
		sort(vec_line.begin(), vec_line.end(), lessthan_x);
	}

	double min_dist = DecideBorderLine(vec_line, vec_gap, result);

    #if DEBUG
		std::cout << "detected min dist is" << min_dist << endl;
	#endif

	if (min_dist > match_thresh)
	{
        #if DEBUG
		    std::cout << "error, min dist is too large: " << min_dist << endl;
		#endif

		return false;
	}

	return true;  
}

#if DEBUG
void ShowResult(vector<LINE> vec_line, cv::Mat& img, string title)
{
	cv::Mat dest = img.clone();

	for (size_t i = 0; i < vec_line.size(); i++)
	{
		CvScalar color = cvScalar(0, 0, 255);
		line(dest, cv::Point(vec_line[i].StPnt.x, vec_line[i].StPnt.y), cv::Point(vec_line[i].EdPnt.x, vec_line[i].EdPnt.y), color, 3);
	}
	cv::namedWindow(title, CV_WINDOW_NORMAL);
	cv::imshow(title, dest);
	cv::waitKey(0);
}
#endif

void GetFourCorners(LINE &up, LINE &down, LINE &left, LINE &right, int expand, int *corner)
{
	itcv::Point temp;

	//judge if upside_down or leftside_right
	int upside_down = 1;
	if (up.StPnt.y + up.EdPnt.y > down.StPnt.y + down.EdPnt.y)
	{
		upside_down = -1;
	}

	int leftside_right = 1;
	if (left.StPnt.x + left.EdPnt.x > right.StPnt.x + right.EdPnt.x)
	{
		leftside_right = -1;
	}

	GetCrossPoint(up.StPnt.x, up.StPnt.y, up.EdPnt.x, up.EdPnt.y,
		left.StPnt.x, left.StPnt.y, left.EdPnt.x, left.EdPnt.y, temp);		

	corner[0] = temp.x - expand * leftside_right;
	corner[1] = temp.y - expand * upside_down;

	GetCrossPoint(up.StPnt.x, up.StPnt.y, up.EdPnt.x, up.EdPnt.y,
		right.StPnt.x, right.StPnt.y, right.EdPnt.x, right.EdPnt.y, temp);  

	corner[2] = temp.x + expand * leftside_right;
	corner[3] = temp.y - expand * upside_down;

    GetCrossPoint(down.StPnt.x, down.StPnt.y, down.EdPnt.x, down.EdPnt.y,
        right.StPnt.x, right.StPnt.y, right.EdPnt.x, right.EdPnt.y, temp);  

    corner[4] = temp.x + expand * leftside_right;
    corner[5] = temp.y + expand * upside_down;

	GetCrossPoint(down.StPnt.x, down.StPnt.y, down.EdPnt.x, down.EdPnt.y,
		left.StPnt.x, left.StPnt.y, left.EdPnt.x, left.EdPnt.y, temp);		

	corner[6] = temp.x - expand * leftside_right;
	corner[7] = temp.y + expand * upside_down;
}

bool ChooseVertBorderBasedLen(vector<LINE> &vec_line, double dist_base, double diff_ratio, vector<LINE> &result)
{
	#if DEBUG
		cout << "choose vert border by method B" << endl;
	#endif

	sort(vec_line.begin(), vec_line.end(), lessthan_x);
	if(vec_line.size() < 2)
		return false;

	result.push_back(vec_line[0]);
	result.push_back(vec_line[vec_line.size() - 1]);

	double dist = Parallel_Line_dist(result[0], result[1]);
	if(abs(dist - dist_base) > dist_base * diff_ratio)
	{
		#if DEBUG
			cout << "method B failed" << endl;
		#endif
		return false;
	}
	
	#if DEBUG
		cout << "method B succeed" << endl;
	#endif

	return  true;
}

/************************************************************************/
/* Name        : getInvoiceVertex
* Function     : input an source image, output an deskew form image.
* Parameters   : src_data - source gray image data
*                width    - source image width
*                height   - source image height
*                widhtStep- source image width step
*				  vertex   - vertex of the table: x1,y1,x2,y2,x3,y3,x4,y4
* Return       : -1(error), 0 (OK)                                      */
/************************************************************************/
EXPORT_API int getVertex(unsigned char* src_data, int width, int height, int widthstep, int *vertex)
{
	//draw the line
    #if DEBUG
        Mat im(height, width, CV_8UC1, src_data);
        Mat im_rgb;
        cvtColor(im, im_rgb, CV_GRAY2BGR);
    #endif

    int line[100000];
    int line_num = dsccline(src_data, width, height, width, true, min_hori_len, valley_thresh, line);  //based real 
    
	if(line_num < 1)
		return -1;

    #if DEBUG
        cout << line_num << endl;
    #endif

    //change result to vector<LINE>	
	vector<LINE>  lines;
	for(int i =0; i < line_num; i++)
	{
		LINE temp(line[4 * i + 0], line[4 * i + 1], line[4 * i + 2], line[4 * i + 3], 1);
		lines.push_back(temp);
	}

    #if(DEBUG)
	    ShowResult(lines, im_rgb, "hori line");
    #endif

	//do line merge
	MergeLine(lines, max_line_gap, true, max_merge_cost);
    #if(DEBUG)
        cout << "after merge: lines remain: " << lines.size() << endl;
        ShowResult(lines, im_rgb, "after merge");
    #endif

	//unique line
	UniqueLine(lines, max_line_gap);
	if(lines.size() < 5)
		return -1;
    #if(DEBUG)
	    cout << "after unique: lines remain: " << lines.size() << endl;
	    ShowResult(lines, im_rgb, "after unique");
    #endif

	//get hori border
	vector<LINE> result;
    vector<double> hori_gap(hori_gap_vec, hori_gap_vec + 4);

	bool match_reslt = GlobalMatch(lines, true, match_thresh, hori_gap, result);
	if(!match_reslt)
		return -1;

	#if(DEBUG)
        ShowResult(result, im_rgb, "all hori lines");
    #endif

    //cal hori line distance
    if(result.size() < 5)
		return -1;    			//output error
    double vert_base = Parallel_Line_dist(result[0], result[4]);
    #if(DEBUG)
        cout << "distance between hor lines is: " << vert_base << endl;
    #endif

    int vert_line[100000];
    int vert_line_num = dsccline(src_data, width, height, width, false, vert_base * vert_ratio, 30, vert_line);  //based real 
	if(vert_line_num < 1)
		return -1;
    #if(DEBUG)
        cout << "detected vert line num is:" << vert_line_num << endl;
    #endif

    //change result to vector<LINE>	
	vector<LINE>  vert_lines;
	for(int i =0; i < vert_line_num; i++)
	{
		LINE temp(vert_line[4 * i + 0], vert_line[4 * i + 1], vert_line[4 * i + 2], vert_line[4 * i + 3], 1);
		vert_lines.push_back(temp);
	}
    #if(DEBUG)
	    ShowResult(vert_lines, im_rgb, "detect vert line");
    #endif

    //get vert border
	vector<LINE> vert_rslt;
    vector<double> vert_gap(vert_gap_vec, vert_gap_vec + 2);
	bool vert_match_rslt = GlobalMatch(vert_lines, false, match_thresh, vert_gap, vert_rslt);
	if(!vert_match_rslt)
		return -1;
	
	if(vert_rslt.size() < 3)
		return -1;

	//verify the distance
	double vert_dist = Parallel_Line_dist(vert_rslt[0], vert_rslt[2]);
	#if(DEBUG)
        cout << "distance between vert lines is:  " << vert_dist << endl;
    #endif
	double vert_dist_ref = vert_base * dist_ratio;
	double diff = (vert_dist - vert_dist_ref) / vert_dist_ref;
	#if(DEBUG)
        cout << "vert lines diff ratio is:  " << diff << endl;
    #endif

	if(abs(diff) > vert_dist_diff)
	{   
		//another method
		vert_rslt.resize(0);
		if(!ChooseVertBorderBasedLen(vert_lines, vert_dist_ref, vert_dist_diff, vert_rslt))
			return -1;
		vert_rslt.push_back(vert_rslt[1]);
	}


    #if(DEBUG)
        ShowResult(vert_rslt, im_rgb, "all vert line");
    #endif

    //get four borders
    vector<LINE> borders;
    borders.push_back(result[0]);
    borders.push_back(result[4]);
    borders.push_back(vert_rslt[0]);
    borders.push_back(vert_rslt[2]);
    #if(DEBUG)
        ShowResult(borders, im_rgb, "all borders");
    #endif

    //get points
    GetFourCorners(borders[0], borders[1], borders[2], borders[3], 0, vertex);
	return 0;
}
