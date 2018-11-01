#ifndef __CONNTREE_H__
#define __CONNTREE_H__
#include "intsigImage.h"

struct LineValley
{
	int x;		//X coordinate
	int ys;	//Start position
	int ye;	//End position

			//For binary yvs=ys, yve=ye
	int yvs;	//Start position of valley
	int yve;	//End position of valley

	unsigned char EdgeGray;				//Gray level of valley edge
	unsigned char gray;					//gray level of bottom of valley?
};

struct ConnNode
{
	int nLtTotal;//Total left connected runlength
	int pLeft;   //Left runlength
	int nRtTotal;//Total right connected runlength
	int pRight;  //Right runlength
	int pUnder;  //Under ConnNode
	int pAbove;  //Above ConnNode
	LineValley v;//Node information
};

struct ConnComp
{
	int nInitialNode;      //Initial runlength
	int nPixelsNum;		   //Number of black pixels
	itcv::BoundBox rcBound;   //Bounding box of connected component
};

class CConnTree
{
public:
	int			m_IsHorConn;    // Indicate horizontal or vertical runlengths
	itcv::BoundBox	m_rcRange;	    // The region to analysed

    //Variable used to store runlengths
	int		m_nDepth;  
    int	   *m_pColHead;			// Column head pointer.
    int		m_pEmptHead;		// Empty column head
    int		m_nMaxNodes;		// Maximum nodes has been allocated. Variable used to allocate memory dynamically.
	int		m_nLeftMostX;
	int		m_nLeft;
	ConnNode *m_pNode;			// Runlengths extracted

    //Variable used to extract connected components
	unsigned char  *m_pFlag;	// Indicate the runlength used or not
	int		*m_pStack;			// Stack used to store runlength ready for processing
	ConnComp *m_pConnComp;		// Connected components
    int		m_nTotalConnComps;	// Toal connected components extracted

public :
    CConnTree() ;
    ~CConnTree() ;

    bool	Initialize(itcv::BoundBox rcRange);
    int		FreeMem();		// Free memory used by CConnTree

    //Methods for testing nodes connectionship
	bool	IsConnected ( int pNode1, int pNode2) ;//Test if two nodes connected
	int		IsConnected2( int pNode1, int pNode2, int &width, int &gap ) ;//Test if two nodes neighboring within the specific region
	bool	IsRightConnected ( int pNode1, int pNode2 ) ;//If right connected
	int		IsRightConnected2 ( int pNode1, int pNode2, int&width, int&gap ) ;//If right connected
	bool	IsLeftConnected ( int pNode1, int pNode2 ) ;//If left connected
	int		IsLeftConnected2 ( int pNode1, int pNode2, int&width, int&gap );//If left connected

    //Methods for Creating and storing runlengths
	bool	AllocNewNodes(int nodeNums=0) ;//Allocate new nodes
    int		CopyValley ( LineValley *v, int valleys ) ;//Create a new runlengths
    int		AddNewCol ( LineValley *v, int valleys, int nColNum ) ;//Add a new runlength
    int		MakeRightConn ( int nColNum ) ;//Store the new runlegth and adjust its right connectionship with other runlengths
    int		MakeLeftConn ( int nColNum ) ;//Store the new runlegth and adjust its left connectionship with other runlengths
    int		GetLeftMostX() ;//Return the most left ordinate m_nLeftMostX
	
    //Connected Components extraction methods
    int		GetAllConnComps() ;//Get all connected components
    int		GetConnComp ( ConnComp &cc, int nInitial) ;//Get one connected component.
	int		GetConnComp ( ConnComp &cc, int nInitial, unsigned char *pFlag);//Get one connected componet.
    int		*GetConnComp ( ConnComp &cc, int *_nTotalNodes ) ;//Get one connected component.
} ;
#endif // __CONNTREE_H__
