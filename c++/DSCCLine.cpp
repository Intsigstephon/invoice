#include "DirLine.h"
#include "LineDetect.h"
#include <stdio.h>

/*************************************************************************************
/*	Name:		GetDistance
/*	Function:	Get distance of two points, interface for 'doulbe'
/*	Parameter:	x1 -- X ordinate of point 1
/*				y1 -- Y ordinate of point 1
/*				x2 -- X ordinate of point 2
/*				y2 -- Y ordinate of point 2
/*	Return:		Distance of two points
/*
/**************************************************************************************/
double GetDistance(double x1, double y1, double x2, double y2)
{
	double dx = x1 - x2;
	double dy = y1 - y2;
	return sqrt(dx*dx + dy*dy);
}

/*************************************************************************************
/*	Name:		DSCCLine
/*	Function:	Use DSCC to detect line
/*	Parameter:	img_info     -- image infomation
/*				Ishorizontal -- true : detect horizontal line
								false: detect vertical line
/*				vec_line     -- store line result
/*	Return:		nums of line detected
/**************************************************************************************/
int DSCCLine(itcv::Mat& img_info, bool Ishorizontal, int Min_linelen, int Valley_gray, vector<itcv::LINE> &vec_line)
{
	// Detect the line
	CDirLine dirline;

	dirline.SetDefaultDetectParams(Valley_gray);
	if (Ishorizontal)
	{
		dirline.m_Param.MinHorLineLength = Min_linelen;
	}
	else
	{
		dirline.m_Param.MinVerLineLength = Min_linelen;
	}

	//Get DSCC and filter
	dirline.DSCCFiltering(img_info, Ishorizontal);

	//Chain to chains
	dirline.MergeChains();

	//Chains to lines
	dirline.ChainsToFORMLINE();

	//store to vector
	int num_lines0 = dirline.m_nLine;
	FORMLINE *ptr0 = dirline.m_pLine;

	int num = 0;
	for (int i = 0; i < num_lines0; i++)
	{
		if (GetDistance(ptr0[i].StPnt.x, ptr0[i].StPnt.y, ptr0[i].EdPnt.x, ptr0[i].EdPnt.y) > Min_linelen)
		{
			LINE line = ptr0[i].convertToLine();
			vec_line.push_back(line);
			num++;
		}
	}
	return num;
}

///*************************************************************************************
///*	Name:		DSCCLine
///*	Function:	Use DSCC to detect line
///*	Parameter:	fnImg -- X ordinate of point 1
///*				Ishorizontal -- true : detect horizontal line
//				false: detect vertical line
///*				vec_line     -- store line result
///*	Return:		nums of line detected
///**************************************************************************************/
//int DSCCLine(char* fnImg, bool Ishorizontal, int Min_linelen, int Valley_gray, vector<int> &vec_line)
//{
//	// Get image data and width \ height
//	ImgInfo img_info;
//	if (false == GetRawData(fnImg, img_info))
//	{
//		printf("Read image: %s failed!\n", fnImg);
//		return -1;
//	}
//
//	return DSCCLine(img_info, Ishorizontal, Min_linelen, Valley_gray,vec_line);
//}

/****************************************************************************************************
/*	Name:		DSCCFiltering
/*	Function:	Filter DSCC, only preserve those created by horizontal line segments
/*	Parameter:	img_info		-- Original image
img_filtered    -- Filtered image
/*	Return:		0  -- Succeed
/*				-1 -- Failed
/***************************************************************************************************/
int DSCCFiltering(itcv::Mat& img_info, bool bHorLine, int Valley_gray, itcv::Mat& img_filtered)
{
	// Detect the line
	CDirLine dirline;

	dirline.SetDefaultDetectParams(Valley_gray);

	return dirline.DSCCFiltering(img_info, bHorLine, img_filtered);
}

/*************************************************************************************
/*	Name:		dsccline
/*	Function:	Use DSCC to detect line
/*	Parameter:	unsigned char* src_data : point or image data
				int width: image width
				int height: image height
				int widthstep: usually image width
/*				Ishorizontal -- true : detect horizontal line
				false: detect vertical line
				int Min_linelen: distance thresh
				int Valley_gray: valley_gray, by default:30
/*	Return:		nums of line detected
				line     -- store line detect result: sx,sy,ex,ey
/**************************************************************************************/
EXPORT_API int  dsccline(unsigned char* src_data, int width, int height, int widthstep,
	bool Ishorizontal, int Min_linelen, int Valley_gray,    
	int *line)
{
	itcv::Mat img(height, width, itcv::DEPTH_8U, 1, src_data);

	vector<itcv::LINE> result;
	int num_line = DSCCLine(img, Ishorizontal, Min_linelen, Valley_gray, result);

	//printf("DSCCLine finished, num of line is %d", num_line);
	assert(num_line == result.size());

	//copy line from vector to int array
	for (size_t i = 0; i < result.size(); i++)
	{
		line[4 * i + 0] = result[i].StPnt.x;
		line[4 * i + 1] = result[i].StPnt.y;
		line[4 * i + 2] = result[i].EdPnt.x;
		line[4 * i + 3] = result[i].EdPnt.y;
	}

	return num_line;
}
