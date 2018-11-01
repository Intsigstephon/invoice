#ifndef _DSCC_LINE_H
#define _DSCC_LINE_H
#include <vector>
#include "intsigImage.h"
#include <math.h>
using namespace std;

/*************************************************************************************
/*	Name:		DSCCLine
/*	Function:	Use DSCC to detect line
/*	Parameter:	img_info     -- image  infomation
/*				Ishorizontal -- true : detect horizontal line
				false: detect vertical line
/*				vec_line     -- store line result
/*	Return:		nums of line detected
/**************************************************************************************/
int  DSCCLine(itcv::Mat& img_info, bool Ishorizontal, int Min_linelen, int Valley_gray, vector<itcv::LINE> &vec_line);

/**************************************************************************************
/*	Name:		DSCCFiltering
/*	Function:	Filter DSCC, only preserve those created by horizontal line segments
/*	Parameter:	img_info		-- Original image
				img_filtered    -- Filtered image
/*	Return:		0  -- Succeed
/*				-1 -- Failed
/**************************************************************************************/
int  DSCCFiltering(itcv::Mat& img_info, bool bHorLine, int Valley_gray, itcv::Mat& img_filtered);

#endif
