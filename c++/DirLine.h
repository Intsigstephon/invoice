#ifndef __DIRLINE_H__
#define __DIRLINE_H__
#include "ConnTree.h"
#include "DSCCLine.h"
#include "intsigImage.h"
using namespace itcv;

#ifndef	PI
#define PI  3.1415926535
#endif

#define	QUEUE_DEPTH			3
#define	PARAM_CHARWIDTH		40
#define	MAX_LINE_GAP	    1000
#define MAXSTRIP	        100	    // Maximum number of strips to split form image in order to speed up the algorithm

struct  PARAMETER
{
	int        ValleyGrayDepth;		//Valley gray depth, for gray form image recognition
	int        MinVerLineLength;	//MinVerLineLength=0 MinHorLineLength=0  Estimate char width automatically
	int        MinHorLineLength;	//Esle   use this two parameters
	int		   MaxLineWidth;		//Maximum line width
	int		   MaxGap;

	bool	   FilterSmallDSCC;		//TRUE -- Filter small DSCC to speed up the algorithm
									//FALSE -- Do not filter to increase robustness to line broken

	bool	   RLS;					//TRUE -- Runlength smoothing  

	bool	   bHasSlantLine;		//TRUE -- Has slant lines on the form
									//FALSE -- No slant lines, delete all slant lines detected

	bool	   bSceneImg;			//Scene image or document image
};

struct CHAIN
{
	int pHead;  // Head node (run-length), etc left most node
	int pTail;  // Tail node (run-length), etc right most node

	//Before merging pLeft and pRight = -1;  
	//After merging if a chain form a chains indepently, pLeft=pRight=-2; 

	int pLeft;  // Left chain
	int pRight; // Right chain
	int Len;	//Chain length
	int xs;	//Start position
	int xe;	//End position
	int fYs;//fYs = fYofChain(Chain, xs)   2000/3/23
	int fYe;//fYe = fYofChain(Chain, xe)   2000/3/23
	int Num;	//Number of run length
	double r;//Straightness
	double Width;//Average width
				 //Variable to MSE straight line approximation
	double SumX;//Sum X
	double SumY;//Sum Y
	double SumXX;//Sum X*X
	double SumXY;	//Sum X*Y

	itcv::BoundBox   rcBound; //09/03/2002
						  //Variable to calculate the long and short axis of the chain
						  //They are different with SumX, SumY, SumXX, and SumXY, which are run length based.
						  //These new variables are pixel based
						  //Added on 09/03/2004
	int		NumPixel;	//Number of black pixels
	double	SX;
	double	SY;
	double	SYY;
	double	SXX;
};

struct CHAINS
{
	int pHead;  // Head chain, etc left most chain
	int pTail;  // Tail chain, etc right most chain
	int xs;//Start position
	int xe;//End position
	int fYs;//fYs = fYofChains(Chains, xs);        2000/3/23
	int fYe;//fYe = fYofChains(Chains, xe);       2000/3/23
	int Num;//Number of run lengths
	int NumPixel; //Number of black pixels
	double r;//Straightness
	double Q;//Quality
	double	Angle;//Angle
	double Width;//Average width
				 //Variable to MSE straight line approximation
	double SumX;//Sum X
	double SumY;//Sum Y
	double SumXX;//Sum X*X
	double SumXY;//Sum X*Y
};

struct FORMLINE
{
	int nIndex; //nIndex pointer to the Chains
	int nStyle; // style of line, solid(0), dash(1), dot(2) or virtual(3)
	itcv::Point StPnt;//Start point
	itcv::Point EdPnt;//End point
	double Angle;//Angle
	double Width;//Average width
	double Q;//Quality
	bool   bSlant;//Slant or not
				  //Added Zheng Y.F.   1999-9-2
	int	   nUseType;//0 -- unused   1 -- rectangle cell    2 -- other cell type

	LINE convertToLine()
	{
		LINE line;
		line.StPnt = StPnt;
		line.EdPnt = EdPnt;
		line.flag = true;
		return line;
	}
};

//Structure used to sort CHAIN
struct INTCHAIN
{
	int n;//CHAIN NO
	int pNext;//Next CHAIN
};

class CDirLine
{
public :
	double		m_nSkewAngle;

    bool		m_bIsHorLine;   //Horizontal or Veritical line
	itcv::BoundBox	m_rcBoundRange; //Bound region of the image analysed

    int			 m_nLine;   //Line number
    FORMLINE	*m_pLine;   //Pointer to Line info

    bool		m_bParamsSet;
	PARAMETER   m_Param;   //Parameters

	//protected :  //it is inner data, so it is protected;
    int		m_nChain;   //Chain number
    CHAIN	*m_pChain;  //Pointer to Chain info
    int		m_nChains;  //Chains number
    CHAINS	*m_pChains; //Pointer to Chains info

    //Variables used to split a image into several strips to speed up the algorithm
	int		m_nStrip;                //Strip number of the image devided
	int		m_nCurTree;              //Current strip under analysis
    CConnTree** m_pTree;             //Each strip analysis by a CConnTree instance
	int		m_nOldChain;             //Number of CHAINs have been extracted before the current strip analysed
	int		m_nChainStart[MAXSTRIP]; //First CHAIN of each strip

public :
    CDirLine();
    ~CDirLine();
    int		FreeMem(); //Free memory

///////////////////////////////////////////////////////////////////////////////////
//	Methods to create Directional Single Connected Chain
///////////////////////////////////////////////////////////////////////////////////

	//Set parameters
    int		SetDetectParams ( bool bIsHorLine, PARAMETER &Param );
	//Set default parameters
	int		SetDefaultDetectParams(int gray= 30);
    //Get pixels on a horizontal image row
	int     AquireHorLineData(unsigned char *data, int widthstep, int n_col, int start_row, int end_row, int *buf);
	//Get pixels on a vertical image column
	int     AquireVerLineData(unsigned char *data, int widthstep, int n_row, int start_col, int end_col, int *buf);
    //Extract runlength on gray image

    int		GetColumnRunLength ( uint8_t *p, int w, int h, int column, int start, int end, LineValley *valley ) ;
    int		GetRowRunLength ( uint8_t *p, int w, int h, int row, int start, int end, LineValley *valley) ;
    int		ColRunLenSmooth(uint8_t *p, int w, int h, int col, LineValley *valley, int& Valleys);
    int		RowRunLenSmooth(uint8_t *p, int w, int h, int row, LineValley *valley, int& Valleys);	


	int		ValleyDetect ( int *p, int start, int end, LineValley *valley, int *s, int Depth, int MaxLen );
	//Extract runlengths
	int		BuildConnTree(unsigned char *data, int widthstep, int width, int height, itcv::BoundBox rRange, bool bwFlag=false);
	//Create DSCC
    int		CalTree() ;
	//Filter CHAINs unlikely to be a line segment
	int		FilterChain();
	//Filter CHAINs unlikely to be a line segment(mostly word)
	int     DSCCFiltering(itcv::Mat& img_info, bool bHorLine, bool bwFlag = false);
	//Filter CHAINs unlikely to be a line segment, save result
	int     DSCCFiltering(itcv::Mat& img_info, bool bHorLine, itcv::Mat& img_filtered);

///////////////////////////////////////////////////////////////////////////////////
//	Methods to extract CHAINS from CHAIN
///////////////////////////////////////////////////////////////////////////////////

	//Calculate black pixels between a CHAIN to a CHAINS
	int		PixelsBetween ( CHAINS &Chains, CHAIN &Chain, int &MaxWidth );
	//Merge CHAINs on the left side
    int		LeftMerge ( INTCHAIN *pTailXTbl, int *pTailXIndex, int &SeedChains ) ;
    //Merge CHAINs on the right side
	int		RightMerge ( INTCHAIN *pHeadXTbl, int *pHeadXIndex, int &SeedChains ) ;
	//Merge CHAINs to get CHAINSes
    int		MergeChains() ;
	//Merge two CHAINSes
	int		MergeChains ( CHAINS &Chains1, CHAINS &Chains2 ) ;
	//Sort CHAIN head, speedup searching
    int		SortChainHead ( INTCHAIN *pTbl, int *pIndex ) ;
    //Sort CHAIN tail, speedup searching
	int		SortChainTail ( INTCHAIN *pTbl, int *pIndex ) ;
    //Sort CHAIN by length, speedup searching
	int		SortChainLen ( int MaxLen, INTCHAIN *pTbl, int *pIndex ) ;
    //Add a CHAIN to a CHAINS
	int		AddChain ( int &nChains, int nChain ) ;
	//Delete redudent CHAINS
    int		DeleteChains ( int nDelChains ) ;
    //Calculate the distance of two CHAINs
	double	ChainDistance ( int nChains, int nChain, double Max ) ;
    //Calculate the distance of two CHAINs
	double	ChainDistance ( CHAINS &Chains1, CHAINS &Chains2, double Max ) ;
    //Initialize a CHAINS
	int		InitChains ( CHAINS &Chains, CHAIN &Chain, int nChain ) ;
	//Test if two CHAINs are connected
    bool	IsConnected ( int c1, int c2 ) ;
	//Test if two CHAINs are connected within a gap
	int		IsConnected2(int c1, int c2, int&width, int&gap);


///////////////////////////////////////////////////////////////////////////////////
//	Methods to get statistical property of DSCC and lines
///////////////////////////////////////////////////////////////////////////////////

	//Find a CHAIN in which strip
	int		InWhichTree(int nChain);
	//Find a CHAIN in which CHAINS
    int		InWhichChains ( int nChain ) ;
	//Get y ordinate using CHAIN statistical property
	double	fYofChain(CHAIN &Chain, double x);
	//Get y ordinate using CHAINS statistical property
	double	fYofChains(CHAINS &Chains, double x);
	//Get y ordinate using Line statistical property
    double	fYofLine ( FORMLINE &Line, double x ) ;
    //Calculate CHAIN quality
	int		ChainStatics( CHAIN &Chain ) ;
    //Calculate CHAINS quality
	double	ChainsQuality ( CHAINS &Chains ) ;

////////////////////////////////////////////////////////////////////////////////////
//	Methods to get line and skew angle
////////////////////////////////////////////////////////////////////////////////////

	//Get Form Line from CHAINS 
	int		ChainsToFORMLINE(FORMLINE &Line, CHAINS &Chains);
	//Get All Form Line from All CHAINS 
	int		ChainsToFORMLINE();
	//Estimate skew angle
	int		EstimateSkew(itcv::Mat& img_info, double &nSkewAngle);
};
#endif