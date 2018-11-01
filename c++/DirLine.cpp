#include "DirLine.h"
#include <algorithm>
#include <string.h>
#include <math.h>
using namespace std;
#define GAPBETWEENLINE  5

/*************************************************************************************
/*	Name:		GetAngle
/*	Function:	Calculate the angle of the line formed by two points
/*	Parameter:	Start -- Starting point
/*				End -- Ending point
/*	Return:		Angle of the line
/*
/**************************************************************************************/
static double GetAngle(itcv::Point Start, itcv::Point End)
{
	if (End.x != Start.x)
	{
		double atg = atan(((double)(Start.y - End.y)) / ((double)(End.x - Start.x)));
		if (End.x > Start.x)
			return atg;
		else
		{
			if (End.y < Start.y)
				return atg + PI;
			else
				return atg - PI;
		}
	}
	else
	{
		if (End.y < Start.y)
			return PI / 2;
		else if (End.y > Start.y)
			return -PI / 2;
		else
			return 0;
	}
}

/*************************************************************************************
/*	Name:		GetDistance
/*	Function:	Get distance of two points
/*	Parameter:	pnt1 -- Point 1
/*				pnt2 -- Point 2
/*	Return:		Distance of two points
/*
/**************************************************************************************/
static double GetDistance(itcv::Point &pnt1, itcv::Point &pnt2)
{
	int dx = pnt1.x - pnt2.x;
	int dy = pnt1.y - pnt2.y;
	return sqrt(double(dx)*dx + dy*dy);
}

/*************************************************************************************
/*	Name:		GetDistance
/*	Function:	Get the distance from a point to a line formed by two points
/*	Parameter:	pnt	-- Point
/*				StPnt -- Starting point of the line
/*				EdPnt -- Ending point of the line
/*	Return:		Distance of a point to a line
/*
/**************************************************************************************/
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

/***********************************************************************
/*	Name:		SmoothProject
/*	Fucntion:	Smooth projection
/*	Parameter:	Proj -- Projecton
/*				nProj -- Projection number
/*	Return:		0  -- Correct
/*				-1 -- Error
/***********************************************************************/
static int	SmoothProject(double *Proj, int nProj)
{
	double	*Proj2 = (double*)malloc(sizeof(double)*nProj);
	memcpy(Proj2, Proj, sizeof(double)*nProj);
	for (int i = 1; i<nProj - 1; i++)
		Proj[i] = floor((Proj2[i - 1] + Proj2[i] + Proj2[i + 1]) / 3);
	free(Proj2);
	return 0;
}

/*************************************************************************************
/*	Name:		GetMidValue
/*	Function:	Get the median value of an array, interface for 'double' type
/*	Parameter:	Value -- Data
/*				nTotal -- Number of data
/*	Return:		Middel value
/**************************************************************************************/
static double GetMidValue(double *Value, int nTotal)
{
	int i, j;
	int *Greater = new int[nTotal];

	for (i = 0; i<nTotal; i++)
		Greater[i] = 0;

	for (i = 0; i<nTotal - 1; i++)
		for (j = i + 1; j<nTotal; j++)
			if (Value[i] != Value[j])
			{
				if (Value[i] > Value[j])
					Greater[i]++;
				else
					Greater[j]++;
			}

	double Mid = Value[0];
	for (i = 1; i<nTotal; i++)
		if (Value[i] < Mid)
			Mid = Value[i];
	for (i = 0; i<nTotal; i++)
		if (Value[i]>Mid && Greater[i] <= nTotal / 2)
			Mid = Value[i];
	delete Greater;
	return Mid;
}

/*************************************************************************************
/*	Name:		GetCrossPoint
/*	Function:	Get cross point of two lines
/*	Parameter:	x1	-- X ordinate of starting point of line 1
/*				y1	-- Y ordinate of starting point of line 1
/*				x2	-- X ordinate of ending point of line 1
/*				y2	-- Y ordinate of ending point of line 1
/*				_x1	-- X ordinate of starting point of line 2
/*				_y1	-- Y ordinate of starting point of line 2
/*				_x2	-- X ordinate of ending point of line 2
/*				_y2	-- Y ordinate of ending point of line 2
/*				CrossPnt -- Cross point
/*	Return:		0  -- Succeed
/*				-1 -- Error
/*
/**************************************************************************************/
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

/*************************************************************************************
/*	Name:		SortData
/*	Function:	Sort data ascendly or decendly, interface for 'double' type
/*	Parameter:	Data -- Data
/*				Num  -- Number of data
/*				bIncrease -- TRUE  ascendly
/*							 FALSE descendly
/*	Return:		void
/*
/**************************************************************************************/
static void SortData(double  Data[], int Num, bool bIncrease)
{
	int i, j;
	double tmp;
	if (!bIncrease)
	{
		int Max;
		for (i = 0; i<Num; i++)
		{
			Max = i;
			for (j = i + 1; j<Num; j++)
				if (Data[j] > Data[Max])
					Max = j;
			if (Max != i)
			{
				tmp = Data[i];
				Data[i] = Data[Max];
				Data[Max] = tmp;
			}
		}
	}
	else
	{
		int Min;
		for (i = 0; i<Num; i++)
		{
			Min = i;
			for (j = i + 1; j<Num; j++)
				if (Data[j] < Data[Min])
					Min = j;
			if (Min != i)
			{
				tmp = Data[i];
				Data[i] = Data[Min];
				Data[Min] = tmp;
			}
		}
	}
	return;
}

///////////////////////////////////////////////////////////////////////////////////
//	Methods to detect form line
///////////////////////////////////////////////////////////////////////////////////
CDirLine::CDirLine()
{
    m_pChain = (CHAIN *)NULL ;
    m_nChain = 0 ;
	m_nOldChain = 0;
    m_pChains = (CHAINS *)NULL ;
    m_nChains = 0 ;
    m_pLine = (FORMLINE *)NULL ;
    m_nLine = 0 ;
		
	m_nCurTree = 0;
	m_pTree = NULL;
	memset(m_nChainStart, 0, sizeof(int)*MAXSTRIP);
	m_nStrip = 1;
	m_Param.MaxGap = 15;
	m_bParamsSet = false;
}

CDirLine::~CDirLine()
{
    FreeMem();
}

/****************************************************************************
/*	Name:		FreeMem
/*	Function:	Free memory
/*	Parameter:	void
/*	return:		0 -- Succeed
/*				1 -- Error
/***************************************************************************/
int CDirLine::FreeMem()
{
	if( m_pTree != NULL)
		for(int i = 0; i <= m_nCurTree; i++) 
			delete m_pTree[i]; 
	free(m_pTree);
	m_pTree = NULL;
	m_nCurTree = 0;
	if ( m_nLine>0 || m_pLine!=(FORMLINE *)NULL )
	{
			free ( m_pLine ) ;
			m_pLine = (FORMLINE *)NULL ;
			m_nLine = 0 ;
    }
    if ( m_nChains>0 || m_pChains!=(CHAINS *)NULL )
    {
        free ( m_pChains ) ;
        m_pChains = (CHAINS *)NULL ;
        m_nChains = 0 ;
    }
    if ( m_nChain>0 || m_pChain!=(CHAIN *)NULL )
    {
        free ( m_pChain ) ;
        m_pChain = (CHAIN *)NULL ;
		m_nOldChain = 0;
        m_nChain = 0 ;
    }
    return 0 ;
}

/***************************************************************************
/*	Name:		ChainStatics
/*	Function:	Calculate some statistical properties of a CHAIN
/*	Parameter:  Chain -- CHAIN
/*	Return:		0  -- Succeed
/*				-1 -- Error
/*
/***************************************************************************/
int  _x[25000];
int CDirLine::ChainStatics(CHAIN &Chain)
{
	double sum_x, sum_ys, sum_ye, sum_xx, sum_xys, sum_xye;
	double sum_dxdy, sum_dxdx, sum_dydy, ax, ay, dx, dy /*sum_err, err*/;
	int    *x, *ys, *ye, *w, *cntWider;
	int Len, pNode, j, num, numpixel;
	double  averWidth;

	Len = Chain.Len;
	if (Len <= 5000)
		x = _x;
	else
	{
		x = (int *)malloc(sizeof(int) * 5 * Len);
		if (x == NULL)		return -1;
	}
	ys = x + Len;
	ye = ys + Len;
	w = ye + Len;
	cntWider = w + Len;
	averWidth = 0;
	num = 0;

	//Initialize the bounding box of the CHAIN	
	Chain.rcBound.top = 10000;
	Chain.rcBound.left = 10000;
	Chain.rcBound.bottom = -1;
	Chain.rcBound.right = -1;


	pNode = Chain.pHead;
	for (j = 0; j<Len; j++)
	{
		x[j] = m_pTree[m_nCurTree]->m_pNode[pNode].v.x & 0x1fffffff;
		ys[j] = m_pTree[m_nCurTree]->m_pNode[pNode].v.yvs;
		ye[j] = m_pTree[m_nCurTree]->m_pNode[pNode].v.yve;
		w[j] = ye[j] - ys[j] + 1;
		if (w[j]<m_Param.MaxLineWidth)
		{
			averWidth += w[j];
			num++;
		}
		cntWider[j] = 0;
		pNode = m_pTree[m_nCurTree]->m_pNode[pNode].pRight;

		//Update the bounding box of the CHAIN
		Chain.rcBound.left = min(Chain.rcBound.left, x[j]);
		Chain.rcBound.right = max(Chain.rcBound.right, x[j]);
		Chain.rcBound.top = min(Chain.rcBound.top, ys[j]);
		Chain.rcBound.bottom = max(Chain.rcBound.bottom, ye[j]);
	}
	if (num != 0)	averWidth = averWidth / num;
	else		averWidth = 3;

	sum_x = sum_ys = sum_ye = sum_xx = sum_xys = sum_xye = 0;
	num = numpixel = 0;
	Chain.Width = 0;

	double sx, sxx, sy, syy;
	sx = sy = sxx = syy = 0;
	for (j = 0; j<Len; j++)
	{
		if (w[j] <= 2 * averWidth || Len < 20)
		{
			Chain.Width += ye[j] - ys[j] + 1;
			sum_x += x[j];
			sum_xx += x[j] * x[j];
			sum_ys += ys[j];
			sum_ye += ye[j];
			sum_xys += x[j] * ys[j];
			sum_xye += x[j] * ye[j];
			num++;

			//Calculate new variables 09/03/2004
			int nn = ye[j] - ys[j] + 1;
			sx += 1.0*x[j] * nn;
			sxx += 1.0*x[j] * x[j] * nn;
			sy += 1.0*(ys[j] + ye[j])*nn / 2;
			numpixel += nn;
			for (int k = ys[j]; k <= ye[j]; k++)
				syy += 1.0*k*k;
		}
	}
	Chain.NumPixel = numpixel;
	Chain.SX = sx;
	Chain.SY = sy;
	Chain.SXX = sxx;
	Chain.SYY = syy;

	Chain.Num = num;
	Chain.SumX = sum_x;
	Chain.SumXX = sum_xx;
	Chain.SumY = (sum_ys + sum_ye) / 2;
	Chain.SumXY = (sum_xys + sum_xye) / 2;
	Chain.xs = x[0];
	Chain.xe = x[Len - 1];
	Chain.fYs = (int)fYofChain(Chain, Chain.xs);
	Chain.fYe = (int)fYofChain(Chain, Chain.xe);

	if (num == 0)
		ax = ay = 0;
	else
	{
		ax = (Chain.SumX + Chain.SumY) / num;
		ay = (Chain.SumY - Chain.SumX) / num;
		Chain.Width /= num;
		averWidth = Chain.Width;
	}
	sum_dxdx = sum_dydy = sum_dxdy = 0;
	m_pTree[m_nCurTree]->m_pNode[Chain.pHead].v.x |= 0xc0000000;
	m_pTree[m_nCurTree]->m_pNode[Chain.pTail].v.x |= 0xc0000000;
	pNode = Chain.pHead;
	pNode = m_pTree[m_nCurTree]->m_pNode[pNode].pRight; //Add 2000/4/10
	Chain.Width = 0;
	num = 0;
	double Err = 0;
	for (j = 1; j<Len - 1; j++)//The start and end run-length may be odd, do not add them  2000/3/23
	{
		if (w[j] <= 1.5*averWidth)//2000-3-15
		{
			dx = x[j] + (ys[j] + ye[j]) / 2 - ax; dy = (ys[j] + ye[j]) / 2 - x[j] - ay;
			sum_dxdx += dx*dx; sum_dydy += dy*dy; sum_dxdy += dx*dy;
			Chain.Width += w[j];
			num++;
			if (Len < 50)
				Err += fabs((ys[j] + ye[j]) / 2 - fYofChain(Chain, x[j]));
		}
		else
			m_pTree[m_nCurTree]->m_pNode[pNode].v.x |= 0xc0000000;//毛刺
		pNode = m_pTree[m_nCurTree]->m_pNode[pNode].pRight;
	}
	if (Err > 1.5*num)	//2000-7-1
	{
		Chain.pLeft = -2;
		Chain.pRight = -2;
	}

	if (num == 0)	Chain.Width = averWidth;
	else			Chain.Width = Chain.Width / num;
	if (Len <= 8 || sum_dxdx<1e-7 || sum_dydy<1e-7)
		Chain.r = 1;
	else
		Chain.r = fabs(sum_dxdy / sqrt(sum_dxdx*sum_dydy));

	if (Len>8 && Chain.Width >= 8)//2000-3-6
	{
		double Angle = ::GetAngle(itcv::Point(Chain.xs, Chain.fYs), itcv::Point(Chain.xe, Chain.fYe));
		Chain.Width = Chain.Width*cos(Angle);
	}
	if (Len > 5000) free(x);
	return 0;
}


/**********************************************************************************
/*	Name:		InWhichChains
/*	Function:	Get the CHAINS which the CHAIN belong to
/*	Parameter:	nChain -- CHAIN NO
/*	return:		-1 -- Error
/*				Others -- The CHAINS which the CHAIN belong to
/*
/*********************************************************************************/
int CDirLine::InWhichChains(int nChain)
{
	if (m_pChain[nChain].pLeft == -1 && m_pChain[nChain].pRight == -1)	return -1;
	int nLeft = nChain;
	while (m_pChain[nLeft].pLeft >= 0)
		nLeft = m_pChain[nLeft].pLeft;
	int nRight = nChain;
	while (m_pChain[nRight].pRight >= 0)
		nRight = m_pChain[nRight].pRight;
	for (int i = 0; i<m_nChains; i++)
		if (m_pChains[i].pHead == nLeft && m_pChains[i].pTail == nRight)
			return i;
	return -1;
}

/*******************************************************************
*Name:		fYofChain
*Function:	Get y ordinate using CHAIN statistical property
*
********************************************************************/
double CDirLine::fYofChain(CHAIN &Chain, double x)
{
	double tmp, a, b;
	tmp = Chain.SumX*Chain.SumX - Chain.Num*Chain.SumXX;
	if (fabs(tmp) > 1e-8)
	{
		a = Chain.SumXY*Chain.SumX - Chain.SumY*Chain.SumXX;
		b = Chain.SumX*Chain.SumY - Chain.Num*Chain.SumXY;
		return (a + b*x) / tmp;
	}
	else if (Chain.Num != 0)
		return Chain.SumY / Chain.Num;
	else
		return Chain.SumY;
}

/*******************************************************************
*Name:		fYofLine
*Function:	Get y ordinate using Form Frame Line statistical property
*
********************************************************************/
double CDirLine::fYofLine(FORMLINE &Line, double x)
{
	if (Line.nIndex >= 0)
		return fYofChains(m_pChains[Line.nIndex], x);
	else
	{
		if (m_bIsHorLine)
		{
			if (Line.StPnt.x == Line.EdPnt.x)
				return 1.0*Line.StPnt.y;
			else
				return Line.StPnt.y + 1.0*(Line.EdPnt.y - Line.StPnt.y)*(x - Line.StPnt.x) / (Line.EdPnt.x - Line.StPnt.x);
		}
		else
		{
			if (Line.StPnt.y == Line.EdPnt.y)
				return 1.0*Line.StPnt.x;
			else
				return Line.StPnt.x + 1.0*(Line.EdPnt.x - Line.StPnt.x)*(x - Line.StPnt.y) / (Line.EdPnt.y - Line.StPnt.y);
		}
	}
}

/***************************************************************************
/*	Name:		SetDetectParams
/*	Function:	Set parameters for building Directional Single Connected
/*				Chain from runlengths
/*	Parameter:	nIsHorLine	-- true  Detect horizontal lines
/*							   false Detect vertical lines
/*				Param -- Parameters
/*	Return:		0  -- Succeed
/*				-1 -- Error
/*
/***************************************************************************/
int CDirLine::SetDetectParams(bool bIsHorLine, PARAMETER &Param)
{
	m_bIsHorLine = bIsHorLine;
	m_bParamsSet = true;
	m_Param = Param;
	return 0;
}

/***************************************************************************
/*	Name:		SetDefaultDetectParams
/*	Function:	Set default parameters for building Directional Single Connected
/*				Chain from runlengths
/*	Parameter:
/*	Return:		0  -- Succeed
/*				-1 -- Error
/***************************************************************************/
int CDirLine::SetDefaultDetectParams(int gray)
{
	m_bParamsSet = true;
	m_Param.ValleyGrayDepth = gray;
	m_Param.MinVerLineLength = 30;
	m_Param.MinHorLineLength = 50;
	m_Param.MaxLineWidth = 15;
    //m_Param.MaxLineWidth = 5;
	m_Param.FilterSmallDSCC = true; //we make dscc, if we want to filter the dsccs that is too short
	m_Param.RLS = true;			//run length smoothing
	m_Param.bHasSlantLine = false;  //has slant line
	m_Param.MaxGap = 15;
	return 0;
}

/***************************************************************************
/*	Name:		AquireHorLineData
/*	Function:	Read a vertical image line data to build horizontal DSCC
*	Parameter:	data -- Pointer to image
/*				widthstep   -- Width of the whole image in unsigned chars
/*				n_col -- X offset
/*				start_row -- Start line
/*				end_row -- End line
/*				buf -- Buffer for read data
/*	Return:		0  -- Succeed
/*				-1 -- Error
/*
/***************************************************************************/
int CDirLine::AquireHorLineData(unsigned char *data, int widthstep, int n_col, int start_row, int end_row, int *buf)
{
	int i;
	data += widthstep * start_row + n_col;
	for (i = start_row; i <= end_row; i++, data += widthstep)
		*buf++ = *data;
	return 0;
}

/***************************************************************************
/*	Name:		AquireVerLineData
/*	Function:	Read a horizontal image line data to build vertical DSCC
/*	Parameter:	data -- Pointer to image
/*				widthstep -- Width of the whole image in unsigned chars
/*				n_row --  Y Offset
/*				start_col -- Start X ordinate
/*				end_col   -- End X ordinate
/*				buf -- Buffer for read data
/*	Return:		0  -- Succeed
/*				-1 -- Error
/*
/***************************************************************************/
int CDirLine::AquireVerLineData(unsigned char *data, int widthstep, int n_row, int start_col, int end_col, int *buf)
{
	data += widthstep*n_row + start_col;
	for (int i = start_col; i <= end_col; i++)
		*buf++ = *data++;
	return 0;
}

/***************************************************************************
/*	Name:		ValleyDetect
/*	Function:	Detect valley on a image line for binary/gray image
/*	Parameter:
/*	Return:		0  -- Succeed
/*				-1 -- Error
/***************************************************************************/
int CDirLine::ValleyDetect(int *p, int start, int end, LineValley *valley, int *s, int Depth, int MaxLen)
{
	//0 --- black             255 --- white
	int i, j, nadir_pos, vsp, vep, vs, nadir, mid, t, l;
	int Valleys = 0;
	int sp = 0;

	nadir = vs = l = *p;      //start value   
	nadir_pos = vsp = start;  //start pos

	bool dire = false;
	bool InValley = false;

	for (i = start; i <= end; i++, p++)  // DIB: upside down, which means the data is upside_down
	{
		t = *p;
		if (InValley)
			s[sp++] = t;

		if (dire)     // downward
		{
			if (t > l)  
			{
				if (l < nadir)
				{
					nadir = l;
					nadir_pos = i - 1;
				}
				if (!InValley)
				{
					if (vs - l > Depth)
					{
						InValley = true;
						for (j = i - vsp; j >= 0; j--)
                        {
                            s[sp++] = *(p - j);
                        }
					}
				}
				dire = false;
			}
			else
			{
				if (!InValley)
                {
					if (i - vsp > MaxLen)
                    {
						if (vs - t < Depth)
						{
							vsp = i - 1;
							vs = l;
							nadir_pos = i;
							nadir = t;
						}
                    }
                }
			}
		}
		else    // upward
		{
			if (t < l)  // value is small than l;
			{
				if (!InValley)
				{
					mid = (nadir + vs * 2) / 3;
					if (l >= mid)
					{
						vsp = i - 1;
						vs = l;

						nadir = t;
						nadir_pos = i;
					}
				}
				else
				{
					mid = (nadir + vs * 2) / 3;
					if (l >= mid)
					{
						InValley = false;
						valley->ys = vsp;
						valley->ye = i - 1;
						int count = 0;
						mid = max(l, vs) - Depth;
						for (j = sp - 2; j >= 0; j--)
                        {
							if (s[j] <= mid)
							{
								vep = vsp + j;
								break;
							}
							else
							{
								count++;
						       if (count == MaxLen)
						        {
									valley->ye -= MaxLen;
									count = 0;
								}
							}
                        }
                      //oox  valley->ye -= count;

						mid = vs - Depth;
						for (j = 0; j<sp; j++)
                        {
							if (s[j] <= mid)
							{
								vsp += j;
								break;
							}
                        }
						valley->yvs = vsp;
						valley->yve = vep;

						if( abs(valley->yve - valley->yvs) < MaxLen )
						{
							valley->EdgeGray = max(l, vs);
							valley->gray = nadir;
							valley++;
							Valleys++;
						}

                        //valley->EdgeGray = max(l, vs);
                        //valley->gray = nadir;
                        //valley++;
                        //Valleys++;
						
						vsp = i - 1;
						vs = l;
						nadir = t;
						nadir_pos = i;
						sp = 0;
					}
                    //else
                    //{
                    //    InValley = false;
                    //    valley->ys = vsp;
                    //    valley->ye = i - 1;
                    //    int count = 0;
                    //    mid = max(l, vs) - Depth;
                    //    for (j = sp - 2; j >= 0; j--)
                    //    {
                    //        if (s[j] <= mid)
                    //        {
                    //            vep = vsp + j;
                    //            break;
                    //        }
                    //        else
                    //        {
                    //            count++;
                    //            if (count == MaxLen)
                    //            {
                    //                valley->ye -= MaxLen;
                    //                count = 0;
                    //            }
                    //        }
                    //    }
                    //    valley->ye -= count;

                    //    mid = vs - Depth;
                    //    for (j = 0; j<sp; j++)
                    //    {
                    //        if (s[j] <= mid)
                    //        {
                    //            vsp += j;
                    //            break;
                    //        }
                    //    }
                    //    valley->yvs = vsp;
                    //    valley->yve = vep;

                    //    valley->EdgeGray = max(l, vs);
                    //    valley->gray = nadir;
                    //    valley++;
                    //    Valleys++;

                    //    vsp = i - 1;
                    //    vs = l;
                    //    nadir = t;
                    //    nadir_pos = i;
                    //    sp = 0;


                    //    //InValley = false;
                    //    //vsp = i - 1;
                    //    //vs = l;
                    //    //nadir = t;
                    //    //nadir_pos = i;
                    //    //sp = 0;
                    //}
				}
				dire = true;
			}
		}
		l = t;
	}

	if (InValley)
	{
		valley->ys = vsp;
		valley->ye = i - 1;
		mid = max(l, vs) - Depth;
		for (j = sp - 2; j >= 0; j--)
        {
			if (s[j] <= mid)
			{
				vep = vsp + j;
				break;
			}
        }
		mid = vs - Depth;
		for (j = 0; j < sp; j++)
        {
			if (s[j] <= mid)
			{
				vsp += j;
				break;
			}
        }
		valley->yvs = vsp;
		valley->yve = vep;
        valley->EdgeGray = max(l, vs);
        valley->gray = nadir;
        Valleys++;
	}

	return Valleys;
}


int	CDirLine::GetColumnRunLength ( uint8_t *p, int w, int h, int column, int start, int end, LineValley *valley )
{
    int i, Valleys=0 ;
    bool bInValley ;

    p += w*start + column;

    for ( i=start, bInValley=false ; i<=end ; i++ )
    {
        if ( !*p)
        {
            if ( !bInValley )
            {
                valley->ys = valley->yvs = i ;
                valley->EdgeGray = 255 ;
                valley->gray = 0 ;
                bInValley = true ;
            }
        }
        else
        {
            if ( bInValley )
            {
                valley->ye = valley->yve = i-1 ;
                valley++ ;
                Valleys++ ;
                bInValley = false ;
            }
        }
        if( i == end )	break;
        p+=w;
    }
    if ( bInValley )
    {
        valley->ye = valley->yve = i-1 ;
        Valleys++ ;
    }
    return Valleys ;
}
int	CDirLine::GetRowRunLength ( uint8_t *p, int w, int h, int row, int start, int end, LineValley *valley)
{
    int i, j, Valleys=0 ;
    bool bInValley ;

    uint8_t tmp ;
    p += w*row + start;

    for ( i=start, j=start, bInValley=false, tmp=*p++ ; i<=end ; i++ )
    {
        if ( tmp == 0)
        {
            if ( !bInValley )
            {
                valley->ys = valley->yvs = i ;
                valley->EdgeGray = 255 ;
                valley->gray = 0 ;
                bInValley = true ;
            }
        }
        else
        {
            if ( bInValley )
            {
                valley->ye = valley->yve = i-1 ;
                valley++ ;
                Valleys++ ;
                bInValley = false ;
            }
        }
        tmp= *p++;    
    }
    if ( bInValley )
    {
        valley->ye = valley->yve = i-1 ;
        Valleys++ ;
    }
    return Valleys ;
}

int	CDirLine::ColRunLenSmooth(uint8_t *p, int w, int h, int col, LineValley *valley, int& Valleys)
{
    int T1 = 3; 
    int T2 = 6;
    int MinLen = 20;

    int i=0, j;
    while(i<Valleys-1)
    {
        if(valley[i+1].ys-valley[i].ye <= T1 || 
            (valley[i+1].ys-valley[i].ye<=T2 && valley[i].ye-valley[i].ys > MinLen && valley[i+1].ye-valley[i+1].ys > MinLen) )
        {
            valley[i].ye = valley[i+1].ye;
            valley[i].yve = valley[i+1].yve;
            for(j=i+1; j<Valleys-1; j++)
                valley[j] = valley[j+1];
            Valleys--;
            continue;
        }
        i++;
    }
    return 0;

}
int CDirLine::RowRunLenSmooth(uint8_t *p, int wb, int h, int row, LineValley *valley, int& Valleys)
{
    int T1 = 3; 
    int T2 = 6;
    int MinLen = 20;

    p+=wb*(h-1-row);
    int i=0, j;
    while(i<Valleys-1)
    {
        if(valley[i+1].ys-valley[i].ye <= T1 || 
            (valley[i+1].ys-valley[i].ye<=T2 && valley[i].ye-valley[i].ys > MinLen && valley[i+1].ye-valley[i+1].ys > MinLen) )
        {
            valley[i].ye = valley[i+1].ye;
            valley[i].yve = valley[i+1].yve;
            for(j=i+1; j<Valleys-1; j++)
                valley[j] = valley[j+1];
            Valleys--;
            continue;
        }
        i++;
    }
    return 0;
}



/***************************************************************************
/*	Name:		BuildConnTree
/*	Function:	Extract runlengths
/*	Parameter:	unsigned char *data  -- point to data
int widthstep        -- widthstep every line
int width,
int height
/8				rRange -- Image region to analysed
/*	Return:		0  -- Succeed
/*				-1 -- Error
/*
/***************************************************************************/
int CDirLine::BuildConnTree(unsigned char *data, int widthstep, int width, int height, itcv::BoundBox rRange, bool bwFlag)
{
	if (m_bParamsSet == false)
		SetDefaultDetectParams();

	int w = width;
	int h = height;
	int wb = widthstep;

	int j, k, st, ed, span, valleys;

	if (m_bIsHorLine)
	{
		st = rRange.left;
		ed = rRange.right;
		span = rRange.bottom - rRange.top + 1;
	}
	else
	{
		st = rRange.top;
		ed = rRange.bottom;
		span = rRange.right - rRange.left + 1;
	}

	if (!m_pTree[m_nCurTree]->Initialize(rRange))
		return -1;

	int *t = (int *)malloc(sizeof(int)*span + 16);
	if (t == (int *)NULL) return -1;
	int *s = (int *)malloc(sizeof(int)*span + 16);
	if (s == (int *)NULL) { free(t); return -1; }
	LineValley *v = (LineValley *)malloc(sizeof(LineValley)*span * 2 + 16);
	if (v == (LineValley *)NULL) { free(s); free(t); return -1; }

	unsigned char *p = data;
	if (p == (unsigned char *)NULL || wb <= 0)
	{
		free(v); free(s); free(t); return -1;
	}

	for (j = st; j <= ed; j++)
	{
        if(bwFlag)
        {
            if ( m_bIsHorLine )
            {   
                valleys = GetColumnRunLength ( p, wb, h, j, rRange.top, rRange.bottom, v) ;
                //if(m_Param.RLS)
                //   ColRunLenSmooth(p, wb, h, j, v, valleys);
            }
            else
            {   
                valleys = GetRowRunLength ( p, wb, h, j, rRange.left, rRange.right, v) ;
                //if(m_Param.RLS)	
                //   RowRunLenSmooth(p, wb, h, j, v, valleys);
            }
        }
        else
        {
            if (m_bIsHorLine)
            {
                AquireHorLineData(p, wb, j, rRange.top, rRange.bottom, s);
                valleys = ValleyDetect(s, rRange.top, rRange.bottom, v, t, m_Param.ValleyGrayDepth, m_Param.MaxLineWidth);
            }
            else
            {
                AquireVerLineData(p, wb, j, rRange.left, rRange.right, s);
                valleys = ValleyDetect(s, rRange.left, rRange.right, v, t, m_Param.ValleyGrayDepth, m_Param.MaxLineWidth);
            }
        }


		for (k = 0; k<valleys; k++)
		{
			v[k].x = j;
            if(bwFlag)
                continue;
#if 0
            int vx = j;
            if( m_bIsHorLine )
            {
                if( v[k].ye == rRange.bottom && v[k].ye < h-1 )
                {
                    while(v[k].ye < h-1 )
                    {
                        uint8_t	DownLineuint8_t = *( p+(h-1-v[k].ye-1)*wb+vx);
                        if (DownLineuint8_t != 0)
                            break;
                        v[k].ye++;
                    }
                }
                if( v[k].ys == rRange.top && v[k].ys > 0 )
                {
                    while(v[k].ys > 0 )
                    {
                        uint8_t	UpLineuint8_t = *( p+(h-1-v[k].ys+1)*wb+vx);
                        if(UpLineuint8_t != 0 )	
                            break;
                        v[k].ys --;
                    }
                }
            }
            else
            { 
                if( v[k].ye == rRange.right && v[k].ye < w-1 )
                {
                    
                    while( v[k].ye < w-1 )
                    {
                        uint8_t	RightPixeluint8_t = *( p+(h-1-vx)*wb+(v[k].ye+1));
                        if (RightPixeluint8_t != 0)	break;
                        v[k].ye++;
                    }
                }
                if( v[k].ys == rRange.left && v[k].ys > 0 )
                {
                    
                    while(v[k].ys > 0 )
                    {
                        uint8_t	LeftPixeluint8_t = *( p+(h-1-vx)*wb+(v[k].ys-1));
                        if(LeftPixeluint8_t != 0)		break;
                        v[k].ys--;
                    }
                }
            }
#endif
		}
		m_pTree[m_nCurTree]->AddNewCol(v, valleys, j);
	}

	free(v);
	free(s);
	free(t);
	return 0;
}

/***************************************************************************
/*	Name:		CalTree
/*	Function:	Build Directional Single Connected Chain from runlengths
/*	Parameter:	void
/*	Return:		0  -- Succeed
/*				-1 -- Error
/*
/***************************************************************************/
int CDirLine::CalTree()
{
	int i, Len, pNode, pNext, pThis;
	static int MaxCnt = 5000;
	if (m_pTree[m_nCurTree]->m_pEmptHead == 0)
		return 0;
	unsigned char* nFlag = (unsigned char*)malloc(m_pTree[m_nCurTree]->m_pEmptHead);
	memset(nFlag, 0, m_pTree[m_nCurTree]->m_pEmptHead);

	if (m_pChain == NULL)
	{
		MaxCnt = 5000;
		m_pChain = (CHAIN *)malloc(sizeof(CHAIN)*MaxCnt);
	}

	if (m_pChain == NULL)
		return -1;
	itcv::BoundBox	CurrentRect = m_pTree[m_nCurTree]->m_rcRange;

	for (i = 0; i<m_pTree[m_nCurTree]->m_nDepth; i++)
	{
		pNode = m_pTree[m_nCurTree]->m_pColHead[i];
		while (pNode >= 0)
		{
			if (nFlag[pNode] == 0)
			{
				Len = 1;
				pThis = pNode;
				while (1)
				{
					nFlag[pThis] = 1;
					if (m_pTree[m_nCurTree]->m_pNode[pThis].nRtTotal == 1)
					{
						pNext = m_pTree[m_nCurTree]->m_pNode[pThis].pRight;
						LineValley v = m_pTree[m_nCurTree]->m_pNode[pNext].v;
						if (v.yve - v.yvs >= 100)
						{
							if (v.yvs == CurrentRect.left || v.yve == CurrentRect.right)
								break;
						}
						if (m_pTree[m_nCurTree]->m_pNode[pNext].nLtTotal == 1 && m_pTree[m_nCurTree]->m_pNode[pNext].pLeft == pThis)
						{
							Len++;
							pThis = pNext;
						}
						else
							break;
					}
					else
						break;
				}

#if 1
				if (m_Param.FilterSmallDSCC)
				{
					if (Len <= 2)	continue;
					if (Len <= 4)
					{
						if (m_pTree[m_nCurTree]->m_pNode[pNode].nLtTotal == 0 && m_pTree[m_nCurTree]->m_pNode[pThis].nRtTotal == 0)
							continue;
					}
				}
#endif // 1

				m_pChain[m_nChain].pHead = pNode;
				m_pChain[m_nChain].pTail = pThis;
				m_pChain[m_nChain].pLeft = m_pChain[m_nChain].pRight = -1;
				m_pChain[m_nChain++].Len = Len;
				if (m_nChain == MaxCnt)
				{
					MaxCnt += 3000;
					m_pChain = (CHAIN *)realloc(m_pChain, sizeof(CHAIN)*MaxCnt);
					if (m_pChain == NULL)
						return -1;
				}
			}
			pNode = m_pTree[m_nCurTree]->m_pNode[pNode].pUnder;
		}
	}

	free(nFlag);
	for (i = m_nOldChain; i<m_nChain; i++)
		ChainStatics(m_pChain[i]);
	return 0;
}

/*******************************************************************
*Name:		fYofChains
*Function:	Get y ordinate using CHAINS statistical property
*
********************************************************************/
double CDirLine::fYofChains(CHAINS &Chains, double x)
{
	double tmp, a, b;
	tmp = Chains.SumX*Chains.SumX - Chains.Num*Chains.SumXX;
	if (fabs(tmp) > 1e-8)
	{
		a = Chains.SumXY*Chains.SumX - Chains.SumY*Chains.SumXX;
		b = Chains.SumX*Chains.SumY - Chains.Num*Chains.SumXY;
		return (a + b*x) / tmp;
	}
	else if (Chains.Num != 0)
		return Chains.SumY / Chains.Num;
	else
		return Chains.SumY;
}

/***************************************************************************
/*	Name:		PixelsBetween
/*	Function:	Number of black pixels between a CHAIN and a CHAINS
/*	Parameter:	Chains -- CHAINS
/*				Chain  -- CHAIN
/*				MaxWidth -- Maximum white gap, changed after calling
/*	Return:		0  -- Succeed
/*				-1 -- Error
/*
/***************************************************************************/
int CDirLine::PixelsBetween(CHAINS &Chains, CHAIN &Chain, int &MaxWidth)
{
	int xs, xe, gap, pix = 0;
	if (Chains.xe < Chain.xs)
	{
		gap = Chain.xs - Chains.xe - 1;
		xs = Chains.xe + 1;
		xe = Chain.xs - 1;
	}
	else if (Chains.xs > Chain.xe)
	{
		gap = Chains.xs - Chain.xe - 1;
		xs = Chain.xe + 1;
		xe = Chains.xs - 1;
	}
	else
		return 0;

	int i, pNode, pHead;
	double w, ys, ye;

	pHead = xs - m_pTree[m_nCurTree]->m_nLeft;
	w = min(4, (int)Chains.Width);
	MaxWidth = -1;
	for (i = xs; i <= xe; i++)
	{
		ys = fYofChains(Chains, i) - w;
		ye = ys + 2 * w;
		if ((pNode = m_pTree[m_nCurTree]->m_pColHead[pHead++]) >= 0)
			while (pNode >= 0)
			{
				LineValley	v = m_pTree[m_nCurTree]->m_pNode[pNode].v;
				if (v.yve >= ys)
				{
					if (v.yvs <= ye)
					{
						pix++;
						MaxWidth = max(MaxWidth, v.yve - v.yvs + 1);
					}
					else
						break;
				}
				pNode = m_pTree[m_nCurTree]->m_pNode[pNode].pUnder;
			}
	}
	return	pix;
}

/***************************************************************************
/*	Name:		RightMerge
/*	Function:	Merge DSCC on the right side
/*	Parameter:
/*	Return:		0  -- Succeed
/*				-1 -- Error
/***************************************************************************/
int CDirLine::RightMerge(INTCHAIN *pHeadXTbl, int *pHeadXIndex, int &SeedChains)
{
	int i, j, k, m, p, xl, xe;
	CHAINS Chains = m_pChains[SeedChains];
	xe = Chains.xe;
	xl = m_pTree[m_nCurTree]->GetLeftMostX();
	int T = min(2 * PARAM_CHARWIDTH, max(2 * (Chains.xe - Chains.xs + 1), PARAM_CHARWIDTH));

	double Dist, DistQueue[QUEUE_DEPTH];
	int BestChain[QUEUE_DEPTH];//QUEUE_DEPTH=3
	for (i = 0; i<QUEUE_DEPTH; i++)
	{
		DistQueue[i] = 1e100;
		BestChain[i] = -1;
	}

	double y = Chains.SumY / Chains.Num;
	double dy = 15;
	for (i = xe + 1 - xl; i<min(xe + T - xl, m_pTree[m_nCurTree]->m_nDepth); i++)
	{
		if (i + xl - xe > DistQueue[QUEUE_DEPTH - 1])
			break;
		p = pHeadXIndex[i];//find the best suitable chain ready for merge
		while (p >= 0)
		{
			j = pHeadXTbl[p].n;
			if (m_pChain[j].pLeft == -1 && m_pChain[j].pRight == -1)//2000-3-15
			{
				if (fabs(m_pChain[j].fYs - fYofChains(Chains, m_pChain[j].xs))<dy)//??
				{
					Dist = ChainDistance(SeedChains, j, DistQueue[QUEUE_DEPTH - 1]);
					for (k = 0; k<QUEUE_DEPTH; k++)
					{	//Insert the Dist into DistQueue array
						if (Dist < DistQueue[k])
						{
							for (m = k + 1; m<QUEUE_DEPTH; m++)
								DistQueue[m] = DistQueue[m - 1], BestChain[m] = BestChain[m - 1];
							DistQueue[k] = Dist;
							BestChain[k] = j;
							break;
						}
					}
				}
			}
			p = pHeadXTbl[p].pNext;//Get next chain
		}
	}
	for (k = 0; k<QUEUE_DEPTH; k++)
	{
		if (DistQueue[k] > 1e99)
			break;
		if (m_pChain[BestChain[k]].Width > 8 && m_pChain[BestChain[k]].Width > 2 * Chains.Width)	continue;
		int gap = m_pChain[BestChain[k]].xs - Chains.xe - 1;
		double MErr = sqrt((DistQueue[k] - gap) / 2);
		double MWidth = Chains.Width;
		if (m_pChain[BestChain[k]].Len < 20 && MErr >= max((double)3, MWidth))		continue;
		if (m_pChain[BestChain[k]].Len >= 20 && MErr >= 1.5*MWidth)	continue;

		int width, gap2;
		int pix = PixelsBetween(Chains, m_pChain[BestChain[k]], width);
		gap2 = gap - pix;
		if (gap2 <= 1)
			return AddChain(SeedChains, BestChain[k]);
		if (pix <= 1 && gap<m_Param.MaxGap)
			return	AddChain(SeedChains, BestChain[k]);
		if ((width >= 0 && width<2 * Chains.Width && gap2 < m_Param.MaxGap) || (width >= 2 * Chains.Width && gap<m_Param.MaxGap / 2))
			return AddChain(SeedChains, BestChain[k]);
	}

	return -1;
}

/***************************************************************************
/*	Name:		LeftMerge
/*	Function:	Merge DSCC on the left side
/*	Parameter:
/*	Return:		0  -- Succeed
/*				-1 -- Error
/***************************************************************************/
int CDirLine::LeftMerge(INTCHAIN *pTailXTbl, int *pTailXIndex, int &SeedChains)
{
	int i, j, k, m, p, xl, xs;
	CHAINS Chains = m_pChains[SeedChains];
	xs = Chains.xs;
	xl = m_pTree[m_nCurTree]->GetLeftMostX();

	int T;
	T = min(2 * PARAM_CHARWIDTH, max(2 * (Chains.xe - Chains.xs + 1), PARAM_CHARWIDTH));

	double Dist, DistQueue[QUEUE_DEPTH];
	int BestChain[QUEUE_DEPTH];
	for (i = 0; i<QUEUE_DEPTH; i++)
	{
		DistQueue[i] = 1e100;
		BestChain[i] = -1;
	}
	double y = Chains.SumY / Chains.Num;
	double dy = 15;
	for (i = xs - 1 - xl; i >= max(xs - T - xl, 0); i--)
	{
		if (xs - (i + xl) > DistQueue[QUEUE_DEPTH - 1])
			break;
		p = pTailXIndex[i];
		while (p >= 0)
		{
			j = pTailXTbl[p].n;
			if (m_pChain[j].pRight == -1 && m_pChain[j].pLeft == -1)//2000-3-15
			{
				if (fabs(m_pChain[j].fYe - fYofChains(Chains, m_pChain[j].xe))<dy)//??
				{
					Dist = ChainDistance(SeedChains, j, DistQueue[QUEUE_DEPTH - 1]);
					for (k = 0; k<QUEUE_DEPTH; k++)
					{
						if (Dist < DistQueue[k])
						{
							for (m = k + 1; m<QUEUE_DEPTH; m++)
								DistQueue[m] = DistQueue[m - 1], BestChain[m] = BestChain[m - 1];
							DistQueue[k] = Dist;
							BestChain[k] = j;
							break;
						}
					}
				}
			}
			p = pTailXTbl[p].pNext;
		}
	}
	for (k = 0; k<QUEUE_DEPTH; k++)
	{
		if (DistQueue[k] > 1e99)
			break;
		if (m_pChain[BestChain[k]].Width > 8 && m_pChain[BestChain[k]].Width > 2 * Chains.Width)	continue;
		int gap = Chains.xs - m_pChain[BestChain[k]].xe - 1;
		double MErr = sqrt(DistQueue[k] - gap);
		double MWidth = Chains.Width;
		if (m_pChain[BestChain[k]].Len < 20 && MErr >= max((double)3, MWidth))		continue;
		if (m_pChain[BestChain[k]].Len >= 20 && MErr >= 1.5*MWidth)	continue;

		int width, gap2;
		int pix = PixelsBetween(Chains, m_pChain[BestChain[k]], width);
		gap2 = gap - pix;
		if (gap2 <= 1)
			return AddChain(SeedChains, BestChain[k]);
		if (pix <= 1 && gap<m_Param.MaxGap)
			return	AddChain(SeedChains, BestChain[k]);
		if ((width >= 0 && width<2 * Chains.Width && gap2 < m_Param.MaxGap) || (width >= 2 * Chains.Width && gap<m_Param.MaxGap / 2))
			return AddChain(SeedChains, BestChain[k]);
	}
	return -1;
}

/***************************************************************************
/*	Name:		MergeChains
/*	Function:	Merge DSCC to get CHAINS
/*	Parameter:	void
/*	Return:		0  -- Succeed
/*				-1 -- Error
/***************************************************************************/
int CDirLine::MergeChains()
{
	if (m_nChain == m_nOldChain)
		return -1;
	INTCHAIN *pHeadXTbl = (INTCHAIN *)malloc(sizeof(INTCHAIN) * (m_nChain - m_nOldChain));
	if (pHeadXTbl == NULL)	return -1;
	int *pHeadXIndex = (int *)malloc(sizeof(int) * m_pTree[m_nCurTree]->m_nDepth);
	if (pHeadXIndex == NULL)	return -1;
	SortChainHead(pHeadXTbl, pHeadXIndex);//Speed up the search process

	INTCHAIN *pTailXTbl = (INTCHAIN *)malloc(sizeof(INTCHAIN) * (m_nChain - m_nOldChain));
	if (pTailXTbl == NULL)	return -1;
	int *pTailXIndex = (int *)malloc(sizeof(int) * m_pTree[m_nCurTree]->m_nDepth);
	if (pTailXIndex == NULL)	return -1;
	SortChainTail(pTailXTbl, pTailXIndex);//Speed up the search process

	INTCHAIN *pLenTbl = (INTCHAIN *)malloc(sizeof(INTCHAIN) * (m_nChain - m_nOldChain));
	if (pLenTbl == NULL)	return -1;
	int *pLenIndex = (int *)malloc(sizeof(int) * (m_pTree[m_nCurTree]->m_nDepth + 1));
	if (pLenIndex == NULL)	return -1;
	SortChainLen(m_pTree[m_nCurTree]->m_nDepth, pLenTbl, pLenIndex);//Sort Seed Chain

	int i, j;
	int nOldChains = m_nChains;
	static int MaxCnt = 1000;
	if (m_pChains == NULL)
		m_pChains = (CHAINS *)malloc(sizeof(CHAINS)*MaxCnt);

	if (m_pChains == NULL)
		return -1;

	for (i = m_pTree[m_nCurTree]->m_nDepth; i >= 0; i--)
	{
		int p = pLenIndex[i];
		while (p >= 0)
		{
			j = pLenTbl[p].n;
			CHAIN Chain = m_pChain[j];
			if (Chain.pRight == -1 && Chain.pLeft == -1 && Chain.Num >= Chain.Len / 2)
			{
				InitChains(m_pChains[m_nChains], Chain, j);
				int SeedChains = m_nChains++;
				if (m_nChains == MaxCnt)
				{
					MaxCnt += 1000;
					m_pChains = (CHAINS *)realloc(m_pChains, sizeof(CHAINS)*MaxCnt);
					if (m_pChains == NULL)  	return -1;
				}
				while (RightMerge(pHeadXTbl, pHeadXIndex, SeedChains) == 0);
				while (LeftMerge(pTailXTbl, pTailXIndex, SeedChains) == 0);

				m_pChains[SeedChains].Angle = GetAngle(itcv::Point(m_pChains[SeedChains].xs, m_pChains[SeedChains].fYs),
					itcv::Point(m_pChains[SeedChains].xe, m_pChains[SeedChains].fYe));

				if (m_pChain[j].pLeft<0 && m_pChain[j].pRight<0)
                {
					if (m_pChain[j].xe - m_pChain[j].xs < 15)
                    {
                        m_nChains--;
                    }
					else
					{
						m_pChain[j].pLeft = -2;			m_pChain[j].pRight = -2;
					}
                }
			}
			p = pLenTbl[p].pNext;
		}
	}

	itcv::BoundBox   rcRange = m_pTree[m_nCurTree]->m_rcRange;
	vector<double> tons(12, 0);

	for (i = nOldChains; i<m_nChains; i++)//Delete chains in overlaped area which are extracted twice
	{
		if (m_nCurTree != 0)
		{
			if (m_bIsHorLine && m_nCurTree> 0 && m_pChains[i].fYs < rcRange.top + 10 && m_pChains[i].fYe < rcRange.top + 10)
			{
				DeleteChains(i);
				i--;
			}
			else if (!m_bIsHorLine && m_nCurTree> 0 && m_pChains[i].fYs < rcRange.left + 10 && m_pChains[i].fYe < rcRange.left + 10)
			{
				DeleteChains(i);
				i--;
			}
		}

		if (m_nCurTree != m_nStrip - 1)
		{
			if (m_bIsHorLine && rcRange.bottom < m_rcBoundRange.bottom - 20 && m_pChains[i].fYs > rcRange.bottom - 10 && m_pChains[i].fYe > rcRange.bottom - 10)
			{
				DeleteChains(i);
				i--;
			}
			else if (!m_bIsHorLine && rcRange.right < m_rcBoundRange.right - 20 && m_pChains[i].fYs  > rcRange.right - 10 && m_pChains[i].fYe > rcRange.right - 10)
			{
				DeleteChains(i);
				i--;
			}
		}

		/* delete Chains based Q */
		m_pChains[i].Q = ChainsQuality(m_pChains[i]);


		//if (m_pChains[i].Q < 0.9)
		//{
		//	DeleteChains(i);
		//	i--;
		//}
	}


	free(pLenIndex);
	free(pLenTbl);
	free(pTailXIndex);
	free(pTailXTbl);
	free(pHeadXIndex);
	free(pHeadXTbl);
	return 0;
}

/***************************************************************************
/*	Name:		MergeChains
/*	Function:	Merge two CHAINSs
/*	Parameter:	Chains1 -- CHAINS 1, changed after calling
/*				Chains2 -- CHAINS 2
/*	Return:		0  -- Succeed
/*				-1 -- Error
/***************************************************************************/
int CDirLine::MergeChains(CHAINS &Chains1, CHAINS &Chains2)
{
	if (Chains1.xe < Chains2.xe)
	{
		if (m_pChain[Chains1.pTail].pLeft == Chains2.pHead ||
			m_pChain[Chains2.pHead].pRight == Chains1.pTail)
			return -1;//2000-3-14
		m_pChain[Chains1.pTail].pRight = Chains2.pHead;
		m_pChain[Chains2.pHead].pLeft = Chains1.pTail;
		Chains1.pTail = Chains2.pTail;
		Chains1.xe = Chains2.xe;
	}
	else if (Chains1.xs > Chains2.xs)
	{
		if (m_pChain[Chains1.pHead].pRight == Chains2.pTail ||
			m_pChain[Chains2.pTail].pLeft == Chains1.pHead)
			return -1;//2000-3-14
		m_pChain[Chains1.pHead].pLeft = Chains2.pTail;
		m_pChain[Chains2.pTail].pRight = Chains1.pHead;
		Chains1.pHead = Chains2.pHead;
		Chains1.xs = Chains2.xs;
	}
	else
		return -1;

	Chains1.SumX += Chains2.SumX;
	Chains1.SumY += Chains2.SumY;
	Chains1.SumXX += Chains2.SumXX;
	Chains1.SumXY += Chains2.SumXY;
	Chains1.Width = (Chains1.Width * Chains1.Num + Chains2.Width*Chains2.Num) / (Chains1.Num + Chains2.Num);
	Chains1.Num += Chains2.Num;
	Chains1.NumPixel += Chains2.NumPixel;
	Chains1.fYs = (int)fYofChains(Chains1, Chains1.xs);
	Chains1.fYe = (int)fYofChains(Chains1, Chains1.xe);
	return 0;
}

/***************************************************************************
/*	Name:		SortChainHead
/*	Function:	Sort DSCC according to head to accelerate the merging process
/*	Parameter:
/*	Return:		0  -- Succeed
/*				-1 -- Error
/*
/***************************************************************************/
int CDirLine::SortChainHead(INTCHAIN *pTbl, int *pIndex)
{
	int i, n, xl, x;

	n = m_pTree[m_nCurTree]->m_nDepth;
	for (i = 0; i<n; i++)
		pIndex[i] = -1;
	xl = m_pTree[m_nCurTree]->GetLeftMostX();
	for (i = 0; i<m_nChain - m_nOldChain; i++)
	{
		x = m_pTree[m_nCurTree]->m_pNode[m_pChain[i + m_nOldChain].pHead].v.x & 0x1fffffff;
		n = x - xl;
		pTbl[i].n = i + m_nOldChain;
        pTbl[i].pNext = pIndex[n];
		pIndex[n] = i;
	}
	return 0;
}

/***************************************************************************
/*	Name:		SortChainTail
/*	Function:	Sort DSCC according to tail to accelerate the merging process
/*	Parameter:
/*	Return:		0  -- Succeed
/*				-1 -- Error
/*
/***************************************************************************/
int CDirLine::SortChainTail(INTCHAIN *pTbl, int *pIndex)
{
	int i, n, xl, x;

	n = m_pTree[m_nCurTree]->m_nDepth;
	for (i = 0; i<n; i++)
		pIndex[i] = -1;
	xl = m_pTree[m_nCurTree]->GetLeftMostX();
	for (i = 0; i<m_nChain - m_nOldChain; i++)
	{
		x = m_pTree[m_nCurTree]->m_pNode[m_pChain[i + m_nOldChain].pTail].v.x & 0x1fffffff;
		n = x - xl;
		pTbl[i].n = i + m_nOldChain;
		pTbl[i].pNext = pIndex[n];
		pIndex[n] = i;
	}
	return 0;
}

/***************************************************************************
/*	Name:		SortChainLen
/*	Function:	Sort DSCC according to length to accelerate the merging process
/*	Parameter:
/*	Return:		0  -- Succeed
/*				-1 -- Error
/*
/***************************************************************************/
int CDirLine::SortChainLen(int MaxLen, INTCHAIN *pTbl, int *pIndex)
{
	int i, n;

	for (i = 0; i <= MaxLen; i++)
		pIndex[i] = -1;
	for (i = 0; i<m_nChain - m_nOldChain; i++)
		if ((n = m_pChain[i + m_nOldChain].xe - m_pChain[i + m_nOldChain].xs + 1) <= MaxLen)
		{
			pTbl[i].n = i + m_nOldChain;
			pTbl[i].pNext = pIndex[n];
			pIndex[n] = i;
		}
	return 0;
}

/***************************************************************************
/*	Name:		AddChain
/*	Function:	Add a CHAIN to a CHAINS
/*	Parameter:	nChains -- CHAINS NO.
/*				nChain -- CHAIN NO.
/*	Return:		0  -- Succeed
/*				-1 -- Error
/*
/***************************************************************************/
int CDirLine::AddChain(int &nChains, int nChain)
{
	if (nChains < 0 || nChains >= m_nChains || nChain < 0 || nChain >= m_nChain)
		return -1;
	if (m_pChain[nChain].pLeft<0 && m_pChain[nChain].pRight<0)
	{
		CHAINS tmpChains;
		InitChains(tmpChains, m_pChain[nChain], nChain);
		MergeChains(m_pChains[nChains], tmpChains);
		return 0;
	}
	else if (m_pChain[nChain].pLeft<0 || m_pChain[nChain].pRight<0)
	{
		int nDelChains = InWhichChains(nChain);
		if (nDelChains >= 0)
		{
			MergeChains(m_pChains[nChains], m_pChains[nDelChains]);
			return 0;
		}
		else
			return -2;
	}
	else
		return -1;
}

/***************************************************************************
/*	Name:		DeleteChains
/*	Function:	Delete a CHAINS
/*	Parameter:	nDelChains -- CHAINS NO to delete
/*	Return:		0  -- Succeed
/*				-1 -- Error
/*
/***************************************************************************/
int CDirLine::DeleteChains(int nDelChains)
{
	if (nDelChains >= 0 && nDelChains<m_nChains)
	{
		for (int i = nDelChains; i<m_nChains - 1; i++)
			m_pChains[i] = m_pChains[i + 1];
		m_nChains--;
		return 0;
	}
	else
		return -1;
}

/***************************************************************************
/*	Name:		ChainDistance
/*	Function:	Distance of two CHAINS
/*	Parameter:	Chains1 -- CHAINS 1
/*				Chains2 -- CHAINS 2
/*				Max -- ???
/*	Return:		Distance of two CHAINS
/*
/***************************************************************************/
double CDirLine::ChainDistance(CHAINS &Chains1, CHAINS &Chains2, double Max)
{
	int Len1 = Chains1.xe - Chains1.xs + 1;
	int Len2 = Chains2.xe - Chains2.xs + 1;
	LineValley  v;

	CHAINS comChains;
	comChains = Chains1;
	int num, MaxNum, pNode, nChain, gap;
	double x, Distance, dy, MaxDy;

	MaxNum = Chains2.xe - Chains2.xs + 1;
	Distance = 0;
	num = 0;
	Max *= MaxNum;
	MaxDy = PARAM_CHARWIDTH;

	if (Chains2.xe < Chains1.xs)   // at the left side
	{
		x = Chains2.xe;
		gap = Chains1.xs - Chains2.xe - 1;
		nChain = Chains2.pTail;
		while (Chains2.xe - x<MaxNum && nChain >= 0)
		{
			pNode = m_pChain[nChain].pTail;
			while (Chains2.xe - x<MaxNum && pNode >= 0)
			{
				v = m_pTree[m_nCurTree]->m_pNode[pNode].v;
				if (!(v.x & 0xc0000000))
				{
					x = v.x & 0x1fffffff;
					dy = (v.yvs + v.yve) / 2 - fYofChains(comChains, x);
					Distance += dy*dy;
					if (Distance > Max)
						return 1.7e308;
					num++;
				}
				if (pNode != m_pChain[nChain].pHead)
					pNode = m_pTree[m_nCurTree]->m_pNode[pNode].pLeft;
				else
					break;
			}
			if (nChain != Chains2.pHead)
				nChain = m_pChain[nChain].pLeft;
			else
				break;
		}
	}
	else if (Chains2.xs > Chains1.xe)  // at the right side
	{
		x = Chains2.xs;
		gap = Chains2.xs - Chains1.xe - 1;
		nChain = Chains2.pHead;
		while (x - Chains2.xs<MaxNum && nChain >= 0)
		{
			pNode = m_pChain[nChain].pHead;
			while (x - Chains2.xs<MaxNum && pNode >= 0)
			{
				v = m_pTree[m_nCurTree]->m_pNode[pNode].v;
				if (!(v.x & 0xc0000000))
				{
					x = v.x & 0x1fffffff;
					dy = (v.yvs + v.yve) / 2 - fYofChains(comChains, x);
					Distance += dy*dy;
					if (Distance > Max)
						return 1.7e308;
					num++;
				}
				if (pNode != m_pChain[nChain].pTail)
					pNode = m_pTree[m_nCurTree]->m_pNode[pNode].pRight;
				else
					break;
			}
			if (nChain != Chains2.pTail)
				nChain = m_pChain[nChain].pRight;
			else
				break;
		}
	}
	else
		return 1.7e308;

	if (num > 0)
		return Distance / num + gap;
	else
		return 1.7e308;
}

/***************************************************************************
/*	Name:		ChainDistance
/*	Function:	Distance of a CHAIN to a CHAINS
/*	Parameter:	nChains -- CHAINS NO
/*				nChain -- CHAIN  NO
/*				Max -- ???
/*	Return:		Distance of a CHAIN to a CHAINS
/*
/***************************************************************************/
double CDirLine::ChainDistance(int nChains, int nChain, double Max)
{
	if (m_pChain[nChain].pLeft<0 && m_pChain[nChain].pRight<0)
	{
		CHAINS tmpChains;
		InitChains(tmpChains, m_pChain[nChain], nChain);
		return ChainDistance(m_pChains[nChains], tmpChains, Max);
	}
	else
		return ChainDistance(m_pChains[nChains], m_pChains[InWhichChains(nChain)], Max);
}

/***************************************************************************
/*	Name:		ChainsQuality
/*	Function:	Calculate quality of a CHAINS
/*	Parameter:	Chains -- CHAINS
/*	Return:		Quality of a CHAINS
/*
/***************************************************************************/
double CDirLine::ChainsQuality(CHAINS &Chains)
{
	double gap, gap_paddings, sum_num, len, weight;
	double sum_r;
	int segs, nChain, nNext;

	gap = gap_paddings = sum_num = 0;
	sum_r = 0;
	len = Chains.xe - Chains.xs + 1;
	weight = segs = 0;
	nChain = Chains.pHead;
	while (nChain >= 0)
	{
		CHAIN Chain = m_pChain[nChain];
		sum_num += Chain.Num + 2;//The start and end Run-Lengths are removed, so add 2 to detect dotted lines. 2000/4/22
		sum_r += Chain.r*Chain.Len;
		weight += Chain.Len;
		nNext = Chain.pRight;
		int width = 0;
		int gap2 = 0;
		if (nNext >= 0)
		{
			if (IsConnected2(nChain, nNext, width, gap2) == 0)
			{
				if (width < 2 * Chains.Width)
				{
					gap += gap2;
					sum_num += m_pChain[nNext].xs - Chain.xe - 1 - gap2;
				}
				else	gap += m_pChain[nNext].xs - Chain.xe - 1;
				segs++;
			}
			else if (width < 2 * Chains.Width)
				sum_num += m_pChain[nNext].xs - Chain.xe - 1;
		}
		if (nChain == Chains.pTail)		break;
		else
			nChain = nNext;
	}

	double fgap;
	if (segs == 0)
		fgap = 1;
	else
		fgap = sqrt(1 - gap / segs / len);

	double	sum_dxdy, sum_dxdx, sum_dydy, dx, dy;
	sum_dxdy = sum_dxdx = sum_dydy = 0;
	double  ax = (Chains.SumX + Chains.SumY) / Chains.Num;
	double  ay = (Chains.SumY - Chains.SumX) / Chains.Num;
	int		x, y, nTree;
	LineValley  v;

	int j = Chains.pHead;
	nTree = 0;
	while (j >= m_nChainStart[nTree] && nTree<MAXSTRIP)	nTree++;
	while (j >= 0)
	{
		int pNode = m_pChain[j].pHead;
		while (pNode >= 0)
		{
			v = m_pTree[nTree]->m_pNode[pNode].v;
			x = v.x & 0x1fffffff;
			if (!(v.x & 0xc0000000))
			{
				y = (v.yvs + v.yve) / 2;
				dx = x + y - ax; dy = y - x - ay;
				sum_dxdx += dx*dx; sum_dydy += dy*dy; sum_dxdy += dx*dy;
			}
			if (pNode != m_pChain[j].pTail)
				pNode = m_pTree[nTree]->m_pNode[pNode].pRight;
			else
				break;
		}
		j = m_pChain[j].pRight;
	}
	double tmp = sqrt(sum_dxdx*sum_dydy);
	if (tmp<1e-8)
		Chains.r = 1.0;
	else
		Chains.r = fabs(sum_dxdy / tmp);

	Chains.Angle = GetAngle(itcv::Point(Chains.xs, Chains.fYs), itcv::Point(Chains.xe, Chains.fYe));
	if (fabs(Chains.Angle) < PI * 30 / 180)
		Chains.Q = sum_r / weight * Chains.r * Chains.r * fgap;
	else
		Chains.Q = fgap;

	if (segs == 0)//连续性好的直线
		Chains.Q = Chains.Q * pow(sum_num / (len - gap), 0.25);
	else if (segs <= 2)
		Chains.Q = Chains.Q * pow(sum_num / (len - gap), 0.333);
	else//断裂的直线
		Chains.Q = Chains.Q * pow(sum_num / (len - gap), 0.5);

	return Chains.Q;
}

/***************************************************************************
/*	Name:		IsConnected
/*	Function:	Test if two CHAINs are connected
/*	Parameter:	c1 -- CHAIN 1
/*				c2 -- CHAIN 2
/*				Max -- ???
/*	Return:		TRUE -- Connected
/*				FALSE -- Not connected
/*
/***************************************************************************/
bool CDirLine::IsConnected(int c1, int c2)
{
	int nTree = 0;
	while (c1 >= m_nChainStart[nTree])		nTree++;
	bool  bConnected = false;
	if (m_pChain[c1].xe < m_pChain[c2].xs)
		bConnected = m_pTree[nTree]->IsConnected(m_pChain[c1].pTail, m_pChain[c2].pHead);
	else
		bConnected = m_pTree[nTree]->IsConnected(m_pChain[c1].pHead, m_pChain[c2].pTail);
	return bConnected;
}

/***************************************************************************
/*	Name:		IsConnected2
/*	Function:	Test if two CHAINs are connected within a tolerance
/*	Parameter:	c1 -- CHAIN 1
/*				c2 -- CHAIN 2
/*				width -- ???
/*				gap -- ???
/*	Return:		TRUE -- Connected
/*				FALSE -- Not connected
/*
/***************************************************************************/
int CDirLine::IsConnected2(int c1, int c2, int &width, int &gap)
{
	int nTree = 0;
	while (c1 >= m_nChainStart[nTree])		nTree++;
	if (m_pChain[c1].xe < m_pChain[c2].xs)
		return m_pTree[nTree]->IsConnected2(m_pChain[c1].pTail, m_pChain[c2].pHead, width, gap);
	else
		return m_pTree[nTree]->IsConnected2(m_pChain[c1].pHead, m_pChain[c2].pTail, width, gap);
}

/***************************************************************************
/*	Name:		InitChains
/*	Function:	Initialze a CHAINS structure with a CHAIN
/*	Parameter:	Chains -- CHAINS to initialzie
/*				Chain -- CHAIN to be copied to CHAINS
/*				nChain -- CHAIN NO
/*	Return:		0
/*
/***************************************************************************/
int CDirLine::InitChains(CHAINS &Chains, CHAIN &Chain, int nChain)
{
	Chains.SumX = Chain.SumX;
	Chains.SumY = Chain.SumY;
	Chains.SumXX = Chain.SumXX;
	Chains.SumXY = Chain.SumXY;
	Chains.Num = Chain.Num;
	Chains.NumPixel = Chain.NumPixel;
	Chains.r = Chain.r;
	Chains.xs = Chain.xs;
	Chains.xe = Chain.xe;
	Chains.fYs = Chain.fYs;
	Chains.fYe = Chain.fYe;

	Chains.pHead = Chains.pTail = nChain;
	Chains.Width = Chain.Width;
	Chains.Q = 1;
	return 0;
}

/****************************************************************************
/*	Name:		FilterChain
/*	Function:	Filter chains not like to be a line segment
/*	Parameter:	void
/*	return:		0 -- Succeed
/*				1 -- Error
/***************************************************************************/
int CDirLine::FilterChain()
{
	int		i, j, k;
	int		w, h;
	double	Height;
	double	MaxHeight = 5;
	for (i = 0, k = 0; i<m_nChain; i++)
	{
		bool	bRemove = false;
		w = m_pChain[i].rcBound.right - m_pChain[i].rcBound.left + 1;
		h = m_pChain[i].rcBound.bottom - m_pChain[i].rcBound.top + 1;
		double	vw = m_pChain[i].SXX / m_pChain[i].NumPixel - pow(m_pChain[i].SX / m_pChain[i].NumPixel, 2);
		vw = sqrt(vw);
		double	vh = m_pChain[i].SYY / m_pChain[i].NumPixel - pow(m_pChain[i].SY / m_pChain[i].NumPixel, 2);
		vh = max(1.0, sqrt(vh));
		if (vw >= 3 * vh)
			bRemove = false;
		else
			if (w >= 3 * h)
				bRemove = false;
			else
				if (h >= 2 * MaxHeight && h >= w)
					bRemove = true;
				else
				{
					Height = 0;
					int		pNode = m_pChain[i].pHead;
					itcv::Point	StPnt(m_pChain[i].xs, m_pChain[i].fYs);
					itcv::Point	EdPnt(m_pChain[i].xe, m_pChain[i].fYe);
					LineValley v;

					for (j = 0; j<m_pChain[i].Len; j++)
					{
						v = m_pTree[m_nCurTree]->m_pNode[pNode].v;
						int x = v.x & 0x1fffffff;
						itcv::Point pys(x, v.yvs);
						itcv::Point pye(x, v.yve);
						Height = max(Height, GetDistance(pys, StPnt, EdPnt));
						Height = max(Height, GetDistance(pye, StPnt, EdPnt));

						pNode = m_pTree[m_nCurTree]->m_pNode[pNode].pRight;
						if (pNode < 0)	break;
					}

					if (Height >= MaxHeight)
						bRemove = true;
				}
		if (bRemove == false)
		{
			m_pChain[k] = m_pChain[i];
			k++;
		}
	}
	m_nChain = k;
	return 0;
}

/****************************************************************************************************
/*	Name:		DSCCFiltering
/*	Function:	Filter DSCC, only preserve those created by horizontal line segments
/*	Parameter:	img_info		-- Original image
/*	Return:		0  -- Succeed
/*				-1 -- Failed
/***************************************************************************************************/
int CDirLine::DSCCFiltering(itcv::Mat& img_info, bool bHorLine, bool bwFlag)
{
	FreeMem();
	m_rcBoundRange.left = 0;
	m_rcBoundRange.right = img_info.cols - 1;
	m_rcBoundRange.top = 0;
	m_rcBoundRange.bottom = img_info.rows - 1;

	m_nStrip = 1;
	m_pTree = (CConnTree**)malloc(sizeof(CConnTree*));

	m_pTree[0] = new CConnTree();
	m_pTree[0]->m_IsHorConn = bHorLine;
	m_bIsHorLine = bHorLine;
	m_nCurTree = 0;
	if (BuildConnTree(img_info.data, img_info.widthStep, img_info.cols, img_info.rows, m_rcBoundRange, bwFlag) != 0)
		return -1;

	if (CalTree() != 0) return -1;  // Get all the chain;
	FilterChain();					// Filter chains not like to be a line segment
	m_nChainStart[0] = m_nChain;

#if 0
    itcv::Mat img_filtered(img_info);
    memset(img_filtered.data, 255, img_info.cols * img_info.rows);

    //Save the filtered image data
    int  nTree = 0;
    int  pNode = 0;
    LineValley v;
    for (int i = 0; i < m_nChain; i++)
    {
        if (i >= m_nChainStart[nTree]) nTree++;
        pNode = m_pChain[i].pHead;
        for (int j = 0; j < m_pChain[i].Len; j++)
        {
            v = m_pTree[nTree]->m_pNode[pNode].v;
            int x = v.x & 0x1fffffff;
            if (m_bIsHorLine)
            {
                unsigned char *p = img_filtered.data + img_info.widthStep * v.yvs + x;
                for (int k = v.yvs; k <= v.yve + 1; k++)
                {
                    *p = 0;  //black;
                    p += img_info.widthStep;
                }
            }
            else
            {
                unsigned char *p = img_filtered.data + img_info.widthStep * x + v.yvs;
                for (int k = v.yvs; k <= v.yve + 1; k++)
                {
                    *p = 0;  //black;
                    p++;
                }
            }
            pNode = m_pTree[nTree]->m_pNode[pNode].pRight;
            if (pNode < 0)  break;
        }
    }

    img_filtered.show("filter");
#endif

	return 0;
}

/****************************************************************************************************
/*	Name:		DSCCFiltering
/*	Function:	Filter DSCC, only preserve those created by horizontal line segments
/*	Parameter:	img_info		-- Original image
				img_filtered    -- Filtered image
/*	Return:		0  -- Succeed
/*				-1 -- Failed
/***************************************************************************************************/
int CDirLine::DSCCFiltering(itcv::Mat& img_info, bool bHorLine, itcv::Mat& img_filtered)
{
	FreeMem();
	m_rcBoundRange.left = 0;
	m_rcBoundRange.right = img_info.cols - 1;
	m_rcBoundRange.top = 0;
	m_rcBoundRange.bottom = img_info.rows - 1;

	m_nStrip = 1;
	m_pTree = (CConnTree**)malloc(sizeof(CConnTree*));

	m_pTree[0] = new CConnTree();
	m_pTree[0]->m_IsHorConn = bHorLine;
	m_bIsHorLine = bHorLine;
	m_nCurTree = 0;
	if (BuildConnTree(img_info.data, img_info.widthStep, img_info.cols, img_info.rows, m_rcBoundRange) != 0)
		return -1;

	if (CalTree() != 0) return -1;  // Get all the chain;
	FilterChain();					// Filter chains not like to be a line segment
	m_nChainStart[0] = m_nChain;

	memset(img_filtered.data, 255, img_info.cols * img_info.rows);
	
	//Save the filtered image data
	int  nTree = 0;
	int  pNode = 0;
	LineValley v;
	for (int i = 0; i < m_nChain; i++)
	{
		if (i >= m_nChainStart[nTree]) nTree++;
		pNode = m_pChain[i].pHead;
		for (int j = 0; j < m_pChain[i].Len; j++)
		{
			v = m_pTree[nTree]->m_pNode[pNode].v;
			int x = v.x & 0x1fffffff;
			if (m_bIsHorLine)
			{
				unsigned char *p = img_filtered.data + img_info.widthStep * v.yvs + x;
				for (int k = v.yvs; k <= v.yve + 1; k++)
				{
					*p = 0;  //black;
					p += img_info.widthStep;
				}
			}
			else
			{
				unsigned char *p = img_filtered.data + img_info.widthStep * x + v.yvs;
				for (int k = v.yvs; k <= v.yve + 1; k++)
				{
					*p = 0;  //black;
					p++;
				}
			}
			pNode = m_pTree[nTree]->m_pNode[pNode].pRight;
			if (pNode < 0)  break;
		}
	}

	return 0;
}

/****************************************************************************************************
/*	Name:		EstimateSkew
/*	Function:	Skew angle estimation based on lines
/*	Parameter:	Image	   -- Image
/*				nSkewAngle -- The estimated skew angle
/*	Return:		0  -- Succeed
/*				-1 -- Failed
/***************************************************************************************************/
int CDirLine::EstimateSkew(itcv::Mat &img_info, double &nSkewAngle)    //基于表格线来计算倾斜角，如果是正的，那么表格线只有直线和竖线
{
	FreeMem();
	m_rcBoundRange.left = 0;
	m_rcBoundRange.right = img_info.cols - 1;
	m_rcBoundRange.top = 0;
	m_rcBoundRange.bottom = img_info.rows - 1;

	m_nStrip = 1;
	m_pTree = (CConnTree**)malloc(sizeof(CConnTree*));

	m_pTree[0] = new CConnTree();
	m_pTree[0]->m_IsHorConn = true;
	m_nCurTree = 0;
	if (BuildConnTree(img_info.data, img_info.widthStep, img_info.cols, img_info.rows, m_rcBoundRange) != 0)
		return -1;

	if (CalTree() != 0) return -1;

	//Filter chain
	FilterChain();	 

	//Merge chains
	MergeChains();   
	m_nChainStart[0] = m_nChain;

	int		i;
	int		MinLen = 100;
	double	nAngleHist[200];
	memset(nAngleHist, 0, 200 * sizeof(double));

	double	AngleDelta = PI / 360;
	for (i = 0; i < m_nChains; i++)
	{
		if (m_pChains[i].Num < MinLen)			continue;
		if (fabs(m_pChains[i].Angle) > PI / 4)	continue;

		int	nIndex = (int)(m_pChains[i].Angle / AngleDelta) + 100;
		if (nIndex < 0 || nIndex > 199)
			continue;
		nAngleHist[nIndex] += m_pChains[i].Num;
	}

	SmoothProject(nAngleHist, 200);

	int		nMaxIndex = 0;
	for (i = 0; i<200; i++)
		if (nAngleHist[i] > nAngleHist[nMaxIndex])
			nMaxIndex = i;

	double	nAngle[10000];
	int		Num = 0;
	nSkewAngle = (nMaxIndex - 100) * AngleDelta;

	for (i = 0; i < m_nChains; i++)
	{
		if (m_pChains[i].Num < MinLen)			continue;
		if (fabs(m_pChains[i].Angle) > PI / 4)	continue;
		if (fabs(m_pChains[i].Angle - nSkewAngle) <= AngleDelta && Num < 1000)
		{
			nAngle[Num] = m_pChains[i].Angle;
			Num++;
		}
	}

	if (Num <= 6)//To few samples for estimation
		; //nSkewAngle = 0;
	else
		nSkewAngle = GetMidValue(nAngle, Num);
	return 0;
}

/****************************************************************************
/*	Name:		ChainsToFORMLINE
/*	Function:	Get a frame line from a CHAINS
/*	Parameter:	Line -- A new frame line
/*				Chains -- CHAINS
/*	return:		0 -- Succeed
/*				1 -- Error
/*
/***************************************************************************/
int CDirLine::ChainsToFORMLINE(FORMLINE &Line, CHAINS &Chains)
{
	int xs, xe, ys, ye;

	if (m_bIsHorLine)
	{
		xs = Chains.xs;
		xe = Chains.xe;
		ys = Chains.fYs;
		ye = Chains.fYe;
	}
	else
	{
		ys = Chains.xs;
		ye = Chains.xe;
		xs = Chains.fYs;
		xe = Chains.fYe;
	}
	Line.StPnt = itcv::Point(xs, ys);
	Line.EdPnt = itcv::Point(xe, ye);
	Line.Angle = GetAngle(Line.StPnt, Line.EdPnt);
	Line.Width = Chains.Width;
	Line.Q = Chains.Q;
	Line.nUseType = 0;
	Line.nStyle = 0;
	Line.bSlant = false;
	return 0;
}

/****************************************************************************
/*	Name:		ChainsToFORMLINE
/*	Function:	Get all frame line from all CHAINS
/*	Parameter:	Line -- A new frame line
/*				Chains -- CHAINS
/*	return:		0 -- Succeed
/*				1 -- Error
/*
/***************************************************************************/
int CDirLine::ChainsToFORMLINE()
{
	int num_chains = m_nChains;
	m_nLine = num_chains;
	m_pLine = (FORMLINE *)malloc(num_chains * sizeof(FORMLINE));

	int result = 0;
	for (int i = 0; i < num_chains; i++)
	{
		result += ChainsToFORMLINE(m_pLine[i], m_pChains[i]);
	}

	//Judge
	if (result)
	{
		return 1;
	}
	else
		return 0;
}
