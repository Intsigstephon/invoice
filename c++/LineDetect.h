#ifndef __LINE_DETECT_H__
#define __LINE_DETECT_H__

#include <string.h>
#include <vector>
using namespace std;

#ifdef __cplusplus
extern "C"
{
#endif

#ifdef WIN32
#define CREATEDLL
#endif

#ifdef CREATEDLL
#define EXPORT_API __declspec(dllexport)
#endif

#ifdef CALLDLL
#define EXPORT_API __declspec(dllimport)
#endif

#ifndef EXPORT_API
#define EXPORT_API
#define __stdcall
#endif

	/*************************************************dsccline/getInvoiceVertex**************************************************/
	/*************************************************************************************
	/*	Name:		dsccline
	/*	Function:	Use DSCC to detect line
	/*	Parameter:	unsigned char* src_data : point to image data
					int width: image width
					int height: image height
					int widthstep: usually image width
	/*				Ishorizontal -- true : detect horizontal line
									false: detect vertical line
					int Min_linelen: distance thresh
					int Valley_gray: valley_gray, by default:30
	/*	Return:		nums of lines detected
					line     -- store line detect result: format: sx,sy,ex,ey....
	/**************************************************************************************/
	EXPORT_API int  dsccline(unsigned char* src_data, int width, int height, int widthstep,
		bool Ishorizontal, int Min_linelen, int Valley_gray,
		int *line);


	/************************************************************************/
	/* Name         : getFormImage
	 * Function     : input an source image, output an deskew form image.
	 * Parameters   : src_data - source gray image data
	 *                width    - source image width
	 *                height   - source image height
	 *                widhtStep- source image width step
	 *				  vertex   - vertex of the table: x1,y1,x2,y2,x3,y3,x4,y4
	 * Return       : -1(error), 0 (OK)                                      */
	 /************************************************************************/
	EXPORT_API int getVertex(unsigned char* src_data, int width, int height, int widthstep, int *vertex);


#ifdef __cplusplus
}

#endif
#endif
