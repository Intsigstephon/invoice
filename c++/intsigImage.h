#ifndef INTSIG_IMAGE_H
#define INTSIG_IMAGE_H

#include <limits.h>

#include <vector>
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <string>
#include <cmath>

namespace itcv
{
	/// The meaning of this value is the following:
	/// SS -- 00 - 0  (sx1, sy1) <--> (sx2, sy2)
	/// SE -- 01 - 1  (sx1, sy1) <--> (ex2, ey2)
	/// ES -- 10 - 2  (ex1, ey1) <--> (sx2, sy2)
	/// EE -- 11 - 3  (ex1, ey1) <--> (ex2, ey2)
	/// 
	#define		SS	0
	#define		SE	1
	#define		ES	2
	#define		EE	3

	enum DEPTH_TYPE {
		DEPTH_8U = 0,
		DEPTH_8S,
		DEPTH_16U,
		DEPTH_16S,
		DEPTH_32S,
		DEPTH_32F,
		DEPTH_64F,
	};

	//! interpolation algorithm
	enum
	{
		INTER_NEAREST = 0,  //!< nearest neighbor interpolation
		INTER_LINEAR,		//!< bilinear interpolation
		INTER_CUBIC,	    //!< bicubic interpolation
		INTER_AREA,			//!< area-based (or super) interpolation
	};

	// #define DEPTH_8U 0
	// #define DEPTH_8S 1
	// #define DEPTH_16U 2
	// #define DEPTH_16S 3
	// #define DEPTH_32S 4
	// #define DEPTH_32F 5
	// #define DEPTH_64F 6

	const int DepthSize[] = {
		sizeof(uint8_t),
		sizeof(int8_t),
		sizeof(uint16_t),
		sizeof(int16_t),
		sizeof(int32_t),
		sizeof(float),
		sizeof(double),
	};

	inline int getDepthSize(int depthType)
	{
		assert(depthType >= 0 && depthType <= 6);
		return DepthSize[depthType];
	}
	
	#ifndef  NULL
	#define  NULL  0 
	#endif

	#ifndef MIN
	#define MIN(a,b)  ((a) > (b) ? (b) : (a))
	#endif

	#ifndef MAX
	#define MAX(a,b)  ((a) < (b) ? (b) : (a))
	#endif

	#ifndef ABS
	#define ABS(x) ((x) < 0 ? (-(x)) : (x))
	#endif

	#ifndef PI
	#define PI 3.1415926535897932384626433832795f
	#endif

	#ifndef Round
	#define Round(a) ((int)((a) + ((a) >= 0 ? 0.5 : -0.5)))
	#endif

	#ifndef Floor
	#define Floor(a) ((int)(a) - ((int)(a) > (a)))
	#endif

	#ifndef Ceil
	#define Ceil(a) ((int)(a) + ((int)(a) < (a)))
	#endif

	#ifndef CAST_TO_8U
	#define CAST_TO_8U(t) (uint8_t)(!((t) & ~255) ? (t) : (t) > 0 ? 255 : 0)
	#endif

	#ifndef COLOR_FORMAT_888 
	#define COLOR_FORMAT_888  0
	#endif

	#ifndef COLOR_FORMAT_565
	#define COLOR_FORMAT_565  1
	#endif

	#ifndef COLOR_FORMAT_RGBA 
	#define COLOR_FORMAT_RGBA 2
	#endif

	#define  COLOR_TO_GRAY   1
	#define	 GRAY_TO_COLOR   2
	#define  RGB_TO_RGBALPHA 3
	#define  RGBALPHA_TO_RGB 4
	#define  RGB888_TO_565   5
	#define  RGB565_TO_888   6

	struct Size
	{
		int width, height;

		Size() :width(0), height(0) {}
		Size(int w, int h) :width(w), height(h) {}
		Size(const Size& sz) :width(sz.width), height(sz.height) {}
		Size& operator=(const Size& sz)
		{
			width = sz.width;
			height = sz.height;
			return *this;
		}
		int area() const 
		{
			return width * height;
		}
	};

	struct Point {
		int x, y;

		Point() :x(0), y(0) {}
		Point(int _x, int _y) :x(_x), y(_y) {}
		Point(const Point &p) : x(p.x), y(p.y) {}
		Point operator+(const Point&p)
		{
			Point retP;
			retP.x = x + p.x;
			retP.y = y + p.y;
			return retP;
		}

		Point& operator=(const Point& p)
		{
			x = p.x; y = p.y;
			return *this;
		}

		bool operator==(const Point & p)
		{
			return (x == p.x && y == p.y);
		}

		Point operator*=(double scale)
		{
			Point tmp;
			tmp.x = x * scale;
			tmp.y = y * scale;
			return tmp;
		}
	};

	struct LINE {
		Point	StPnt;
		Point	EdPnt;
		int		flag;

		LINE()
		{
			StPnt.x = 0; StPnt.y = 0;
			EdPnt.x = 0; EdPnt.y = 0;
			flag = 0;
		}

		LINE(int x1, int y1, int x2, int y2, int f = 0)
		{
			StPnt.x = x1; StPnt.y = y1;
			EdPnt.x = x2; EdPnt.y = y2;
			flag = f;
		}

		LINE(itcv::Point p1, itcv::Point p2)
		{
			StPnt = p1;
			EdPnt = p2;
			flag = 0;
		}

		LINE& operator=(const LINE& p)
		{
			StPnt = p.StPnt;
			EdPnt = p.EdPnt;
			flag = p.flag;
			return *this;
		}

		double distance() const
		{
			double dx = EdPnt.x - StPnt.x;
			double dy = EdPnt.y - StPnt.y;
			return sqrt(dx * dx + dy * dy);
		}

		double get_angle() const
		{
			double angle = 0.;
			if (EdPnt.x != StPnt.x)
			{
				double atg = atan(((double)(StPnt.y - EdPnt.y)) / ((double)(EdPnt.x - StPnt.x)));
				if (EdPnt.x > StPnt.x)
					angle = atg;
				else
				{
					if (EdPnt.y < StPnt.y)
						angle = atg + PI;
					else
						angle = atg - PI;
				}
			}
			else
			{
				if (EdPnt.y < StPnt.y)
					angle = PI / 2;
				else if (EdPnt.y > StPnt.y)
					angle = -PI / 2;
				else
					angle = 0;
			}

			return angle;
		}

		double merge_cost(LINE &b) const
		{
			double dx = StPnt.x - b.StPnt.x;
			double dy = StPnt.y - b.StPnt.y;
			double d = sqrt(dx*dx + dy*dy);
			double min = d;

			dx = StPnt.x - b.EdPnt.x;
			dy = StPnt.y - b.EdPnt.y;
			d = sqrt(dx*dx + dy*dy);
			if (d < min) { min = d; }

			dx = EdPnt.x - b.StPnt.x;
			dy = EdPnt.y - b.StPnt.y;
			d = sqrt(dx*dx + dy*dy);
			if (d < min) { min = d; }

			dx = EdPnt.x - b.EdPnt.x;
			dy = EdPnt.y - b.EdPnt.y;
			d = sqrt(dx*dx + dy*dy);
			if (d < min) { min = d; }

			return min;
		}

		double max_dist(LINE &b, int *pwhich) const
		{
			double dx = StPnt.x - b.StPnt.x;
			double dy = StPnt.y - b.StPnt.y;
			double d = sqrt(dx*dx + dy*dy);
			double max = d;
			int which = SS;

			dx = StPnt.x - b.EdPnt.x;
			dy = StPnt.y - b.EdPnt.y;
			d = sqrt(dx*dx + dy*dy);
			if (d > max) { max = d; which = SE; }

			dx = EdPnt.x - b.StPnt.x;
			dy = EdPnt.y - b.StPnt.y;
			d = sqrt(dx*dx + dy*dy);
			if (d > max) { max = d; which = ES; }

			dx = EdPnt.x - b.EdPnt.x;
			dy = EdPnt.y - b.EdPnt.y;
			d = sqrt(dx*dx + dy*dy);
			if (d > max) { max = d; which = EE; }

			if (pwhich) *pwhich = which;
			return max;
		}

		void  merge(LINE &b)
		{
			int which;
			double dist = max_dist(b, &which);

			switch (which)
			{
			case SS:
				EdPnt.x = b.StPnt.x;
				EdPnt.y = b.StPnt.y;
				break;

			case SE:
				EdPnt.x = b.EdPnt.x;
				EdPnt.y = b.EdPnt.y;
				break;

			case ES:
				StPnt.x = b.StPnt.x;
				StPnt.y = b.StPnt.y;
				break;

			case EE:
				StPnt.x = b.EdPnt.x;
				StPnt.y = b.EdPnt.y;
				break;

			default:
				break;
			}

			flag = 1;
		}
	};

	struct LINEARRAY
	{
		std::vector<LINE> array;
		int index;
	};

	struct Rect
	{
		int x, y, width, height;

		Rect() :x(0), y(0), width(0), height(0) {}
		Rect(int _x, int _y, int _w, int _h) :x(_x), y(_y), width(_w), height(_h) {}
		Rect(const Rect& r) :x(r.x), y(r.y), width(r.width), height(r.height) {}
		Rect(const Point &a, const Point &b)
		{
			x = MIN(a.x, b.x);
			y = MIN(a.y, b.y);
			width = MAX(a.x, b.x) - x;
			height = MAX(a.y, b.y) - y;
		}

		Rect& operator=(const Rect& r)
		{
			x = r.x;
			y = r.y;
			width = r.width;
			height = r.height;
			return *this;
		}
		int area() const 
		{
			return width * height;
		}

		bool contains(int pt_x, int pt_y) const
		{
			if (pt_x > x && pt_x < x + width - 1 && pt_y > y && pt_y < y + height - 1)
			{
				return true;
			}
			else
			{
				return false;
			}
		}
	};

	struct BoundBox
	{
		int left, top, right, bottom;

		BoundBox() :left(0), top(0), right(0), bottom(0) {}
		BoundBox(int _l, int _t, int _r, int _b) :left(_l), top(_t), right(_r), bottom(_b) {}
		BoundBox(const BoundBox& box) :left(box.left), top(box.top), right(box.right), bottom(box.bottom) {}
		BoundBox& operator=(const BoundBox& box)
		{
			left = box.left;
			top = box.top;
			right = box.right;
			bottom = box.bottom;
			return *this;
		}
	};

	/*!
	The 2D range class

	This is the class used to specify a continuous subsequence, i.e. part of a contour, or a column span in a matrix.
	*/
	class  Range
	{
	public:
		Range() : start(0), end(0) {}
		Range(int _start, int _end) : start(_start), end(_end) {}
		int size() const { return end - start; }
		bool empty() const { return start == end; }
		static Range all() { return Range(INT_MIN, INT_MAX); };

		int start, end;
	};

	class Mat 
	{
	public:
		Mat()
			:rows(0), cols(0), depth(0), channels(0), data(0)
		{}

		Mat(int _rows, int _cols, int _depth, int _channels)
			:rows(0), cols(0), depth(0), channels(0), data(0)
		{
			create(_rows, _cols, _depth, _channels);
		}

		Mat(Size size, int _depth, int _channels)
			:rows(0), cols(0), depth(0), channels(0), data(0)
		{
			create(size, _depth, _channels);
		}

		Mat(int _rows, int _cols, int _depth, int _channels, void *_data)
			:rows(_rows), cols(_cols), depth(_depth), channels(_channels), data((uint8_t *)_data)
		{
			releaseFlag = false;
			setWidthStep();
			initROI();
		}

		Mat(const Mat& m):rows(0), cols(0), depth(0), channels(0), data(0)
		{
			create(m.getSize(), m.depth, m.channels);
			memcpy(data, m.data, m.rows * m.widthStep);
		}

		Mat(const Mat& m, Rect rect):rows(0), cols(0), depth(0), channels(0), data(0)
		{
			rect.x = MAX(0, rect.x);
			rect.y = MAX(0, rect.y);
			rect.width = MIN(rect.width, m.cols - rect.x);
			rect.height = MIN(rect.height, m.rows - rect.y);

			create(rect.height, rect.width, m.getDepth(), m.channels);

			int sstep = m.widthStep;
			int dstep = widthStep;

			uint8_t *psrc = (uint8_t *)m.data + rect.y * sstep + rect.x * m.getDepthSize() * m.channels;
			uint8_t *pdst = data;

			Size size = getSize();

			for (int i = 0; i < size.height; i++, psrc += sstep, pdst += dstep)
			{
				memcpy(pdst, psrc, dstep);
			}
		}

		Mat row(int y) const { return Mat(*this, Rect(0, y, cols - 1, y + 1)); }
		Mat col(int x) const { return Mat(*this, Rect(x, 0, x + 1, rows - 1)); }

		Mat rowRange(int startrow, int endrow) const
		{
			return Mat(*this, Rect(0, startrow, cols - 1, endrow));
		}

		Mat rowRange(const Range& r) const
		{
			return Mat(*this, Rect(0, r.start, cols - 1, r.end));
		}

		Mat colRange(int startcol, int endcol) const
		{
			return Mat(*this, Rect(startcol, 0, endcol, rows - 1));
		}

		Mat colRange(const Range& r) const
		{
			return Mat(*this, Rect(r.start, 0, r.end, rows - 1));
		}

		~Mat()
		{
			release();
		}

		void create(int _rows, int _cols, int _depth, int _channels)
		{
			create(Size(_cols, _rows), _depth, _channels);
		}

		void create(Size size, int _depth, int _channels)
		{
			if (rows == size.height && cols == size.width && depth == _depth && channels == _channels)
			{
				return;
			}

			release();

			int bufferSize = size.area() * itcv::getDepthSize(_depth) * _channels;
			data = (uint8_t *)malloc(bufferSize);
			if (data)
			{
				rows = size.height;
				cols = size.width;
				depth = _depth;
				channels = _channels;
				releaseFlag = true;
				setWidthStep();
				initROI();
			}
		}

		void release()
		{
			if (data && releaseFlag)
			{
  				free(data);
			}
			data = 0;
			releaseFlag = false;
		}

		void setROI(const Rect& _roi)
		{
			roi.x = MAX(_roi.x, 0);
			roi.y = MAX(_roi.y, 0);
			roi.width = MIN(_roi.width, cols - roi.x);
			roi.height = MIN(_roi.height, rows - roi.y);
		}

		void initROI()
		{
			roi.x = 0;
			roi.y = 0;
			roi.width = cols;
			roi.height = rows;
		}

		Size getSize() const
		{
			Size size(roi.width, roi.height);
			return size;
		}

		void resetROI()
		{
			initROI();
		}

		void getROIImage(Rect rect, Mat& roiImage)
		{
			rect.x = MAX(0, rect.x);
			rect.y = MAX(0, rect.y);
			rect.width = MIN(rect.width, cols - rect.x);
			rect.height = MIN(rect.height, cols - rect.y);

			roiImage.create(rect.height, rect.width, depth, channels);
		}

		uint8_t *getDataOrigin() const 
		{
			return (uint8_t *)data + roi.y * widthStep + roi.x * getDepthSize() * channels;
		}

		int getDepth() const
		{
			return depth;
		}

		int getDepthSize() const 
		{
			return itcv::getDepthSize(depth);
		}

		int getElementSize() const
		{
			return getDepthSize() * channels;
		}

		void setWidthStep()
		{
			widthStep = cols * getDepthSize() * channels;
		}

		int empty() const 
		{
			return (data == 0 || total() == 0);
		}

		size_t total() const
		{
			return (size_t)cols * rows;
		}

		int getOpenCvMatType()
		{
			int type = 0;

	#ifdef DEBUG_WITH_OPENCV
			switch (depth)
			{
			case DEPTH_8U:
				type = CV_8UC(channels);
				break;
			case DEPTH_8S:
				type = CV_8SC(channels);
				break;
			case DEPTH_16U:
				type = CV_16UC(channels);
				break;
			case DEPTH_16S:
				type = CV_16SC(channels);
				break;
			case DEPTH_32S:
				type = CV_32SC(channels);
				break;
			case DEPTH_32F:
				type = CV_32FC(channels);
				break;
			case DEPTH_64F:
				type = CV_64FC(channels);
				break;
			}
	#endif
			return type;
		}

		int rows, cols;
		int depth;
		int widthStep;
		int channels;
		int releaseFlag;
		Rect roi;
		uint8_t *data;
	};

	typedef union itcv32suf
	{
		int i;
		unsigned u;
		float f;
	}
	itcv32suf;
}

#endif

