/*

Copyright (C) 1996-1999
Silicon Graphics Computer Systems, Inc.

Permission to use, copy, modify, distribute and sell this software and
its documentation for any purpose is hereby granted without fee, provided
that the above copyright notice appears in all copies and that both that
copyright notice and this permission notice appear in supporting
documentation. Silicon Graphics makes no representations about the
suitability of this software for any purpose. It is provided "as is"
without express or implied warranty.


Copyright (C) 1994
Hewlett-Packard Company

Permission to use, copy, modify, distribute and sell this software and
its documentation for any purpose is hereby granted without fee, provided
that the above copyright notice appears in all copies and that both that
copyright notice and this permission notice appear in supporting
documentation. Hewlett-Packard Company makes no representations about the
suitability of this software for any purpose. It is provided "as is"
without express or implied warranty.

*/

/*
*******************************************************************
*****                                                         *****
*****  Author: Jonathan R. Wells                              *****
*****  Date  : 21/Jan/2002                                    *****
*****                                                         *****
*****  Description:-                                          *****
*****                                                         *****
*****    Binary heap implementation that allows for dynamic   *****
*****    memory allocations.                                  *****
*****                                                         *****
*****    Based on the STL binary heap codes.                  *****
*****                                                         *****
*****  NOTE: There is no check to see if the number of nodes  *****
*****        have been exceeded. The only check for memory    *****
*****        allocation failure is by the 'assert' macro.     *****
*****                                                         *****
*****  History:-                                              *****
*****                                                         *****
*****    23/Jan/2002 - Change the key to a function call.     *****
*****                - Reverse the inequality for key         *****
*****                  comparsion.                            *****
*****                - Remove the KEY_TYPE parameter from     *****
*****                  bh_insert function.                    *****
*****                                                         *****
*****                                                         *****
*****    24/Jan/2002 - Remove the reduntant test from the     *****
*****                  '_createNextBlock' function            *****
*****                - Fixed bug for boundry detection        *****
*****                - Clean up and re-arrange the macros     *****
*****                - Change 'lastAllocatedBlock' to         *****
*****                  'currentBlock' : name is more          *****
*****                  appropriate                            *****
*****                - Add check for empty heap for 'delete'  *****
*****                  and 'delete_min'                       *****
*****                                                         *****
*******************************************************************
*/

/***  Comment out of the following line if you want to have  ***/
/***  the 'assert' macro activated.                          ***/

//#define NDEBUG
#ifndef _HEAP_CODE
#define _HEAP_CODE

#include <stdlib.h>
#include <assert.h>
#include <limits.h>
#include <memory.h>
//#include <sys/types.h>
//#include <stdint.h>

#ifdef _MSC_VER
// if the compiler does not recognise this type, change it to another int type 8 bytes long
// like long long int
typedef  __int64 ULINTH; //this type myst be 8 bytes long
#else
typedef unsigned long long int ULINTH; //this type myst be 8 bytes long
//#define ULINTH unsigned long long int ULINT //this type myst be 8 bytes long
#endif

struct datum {
	unsigned int Idx1;
	unsigned int Idx2;
	float val;
};

struct datum BADDATUM = { 0xFFFFFFFF , 0xFFFFFFFF , 10e11 };


/***  User define data type  ***/
//#define DATA_TYPE ULINTH  // was UINT 4 bytes - this is for MPI to pack the key
#define DATA_TYPE datum  // GB change 10/10/2020

//#define DATA_TYPE void*

/***  Type for the index in binary heap.  ***/

#if ULONG_MAX == 0xFFFFFFFFUL

#define INDEX_TYPE unsigned  long int //long

#elif UINT_MAX == 0xFFFFFFFFU

#define INDEX_TYPE unsigned long int

#else

#error HELP!!! Unknown integer size. Need an integer size of 4 bytes.

#endif

/*** Type for the key. User define.  **/
#define KEY_TYPE float
#define BADINDEX BADDATUM

/*
*********************************************************************
*****                                                           *****
*****  Structure: node_t                                        *****
*****                                                           *****
*****  Description:-                                            *****
*****                                                           *****
*****    The individual node for the binary heap tree.          *****
*****                                                           *****
*****  Fields:-                                                 *****
*****                                                           *****
*****    index : The location within the heap tree. This field  *****
*****            is used for deleting a node from the binary    *****
*****            tree.                                          *****
*****    key   : The value on which the node is sorted on.      *****
*****    data  : The user define data.                          *****
*****                                                           *****
*********************************************************************
*/

typedef struct node_s
{
  DATA_TYPE  data;
} node_t;



/*
**********************************************************************
*****                                                            *****
*****  Structure: heap_block_t                                   *****
*****                                                            *****
*****  Description:-                                             *****
*****                                                            *****
*****    Contains the array of node pointers sorted by the KEY.  *****
*****                                                            *****
*****  Fields:-                                                  *****
*****                                                            *****
*****    node : node pointer array                               *****
*****                                                            *****
**********************************************************************
*/

typedef struct heap_block_s
{
  node_t* node;   //GB *
} heap_block_t;


/*
************************************************************************
*****                                                              *****
*****  Structure: bheap_t                                          *****
*****                                                              *****
*****  Description:-                                               *****
*****                                                              *****
*****    The main structure that tie all the other structures      *****
*****    together.                                                 *****
*****                                                              *****
*****  Fields:-                                                    *****
*****                                                              *****
*****    block              : Array of blocks which each block     *****
*****                         is part of binary heap.              *****
*****    nodeCount          : Total node(s) in the heap.           *****
*****    emptyBlock         : Number of empty blocks in the heap.  *****
*****    lastAllocatedBlock : Currently used block                 *****
*****                                                              *****
************************************************************************
*/

typedef struct bheap_s
{
  heap_block_t* block;

  INDEX_TYPE nodeCount, emptyBlocks, currentBlock;
} bheap_t;


/*
****************************************************************************
*****                                                                  *****
*****  The relationship between the structures and the memory usage.   *****
*****  -------------------------------------------------------------   *****
*****                                                                  *****
*****                                                                  *****
*****    node_t size      : user defined                               *****
*****    heap_block_t size: 4 bytes                                    *****
*****    bheap_t size     : 12 bytes                                   *****
*****                                                                  *****
*****                                                                  *****
*****    /----------------------------------------------------------\  *****
*****    |         Block         |         Index         |  Memory  |  *****
*****    |-------+---------------+-------+---------------|          |  *****
*****    | Bytes |    Entries    | Bytes |    Entries    | Overhead |  *****
*****    |-------+---------------+-------+---------------+----------|  *****
*****    |    0  |             1 |    4  | 4,294,967,296 |     4B   |  *****
*****    |    1  |           256 |    3  |    16,777,216 |     1KB  |  *****
*****    |    2  |        65,536 |    2  |        65,536 |   256KB  |  *****
*****    |    3  |    16,777,216 |    1  |           256 |    64MB  |  *****
*****    |    4  | 4,294,967,296 |    0  |             1 |    16GB  |  *****
*****    \----------------------------------------------------------/  *****
*****                                                                  *****
*****                                                                  *****
*****  Format of 'index' for node_t.                                   *****
*****                                                                  *****
*****    This variable is an unsigned long integer that is broken up   *****
*****    into two parts. The first part will be the block number and   *****
*****    the second part will be index within that block. Each part    *****
*****    of the 'index' is kept on the byte boundry for speed of       *****
*****    calculation.                                                  *****
*****                                                                  *****
*****    ie. Used logical 'shift' and 'and' operator instead of using  *****
*****    multiplication and division operator because they are         *****
*****    expensive in CPU time.                                        *****
*****                                                                  *****
*****    Therefore, the block must be kept on the power of '2'         *****
*****    boundry. The following formula will be used instead for       *****
*****    calculation the correct block and index.                      *****
*****                                                                  *****
*****      block = value >> blockPart                                  *****
*****      index = value &  indexPart                                  *****
*****                                                                  *****
*****    instead of using                                              *****
*****                                                                  *****
*****      block = value / blockPart                                   *****
*****      index = value % indexPart                                   *****
*****                                                                  *****
****************************************************************************
*/

/*

*************************************
***                               ***
***     Define the boundries      ***
***     --------------------      ***
***                               ***
***  /-------------------------\  ***
***  | Block | Index | IDX_SHF |  ***
***  |-------+-------+---------|  ***
***  |   0   |   4   |    32   |  ***
***  |   1   |   3   |    24   |  ***
***  |   2   |   2   |    16   |  ***
***  |   3   |   1   |     8   |  ***
***  |   4   |   0   |     0   |  ***
***  \-------------------------/  ***
***                               ***
***  Block and Index are the      ***
***  numbers of byte(s) to use.   ***
***  IDX_SHF is the value that    ***
***  need to be set for that row  ***
***                               ***
***  NOTE: Do not use 32 or 0     ***
***        for IDX_SHF as the     ***
***        code is not setup for  ***
***        those values.          ***
***                               ***
*************************************
*/

#define IDX_SHF 16
#define IDX_MASK ((1 << (IDX_SHF)) - 1)

#define BLK_SHF (32 - (IDX_SHF))

/***  Macros for calculating the correct location of the node  ***/
#define BLOCK(A) ((A) >> (IDX_SHF))
#define INDEX(A) ((A) & (IDX_MASK))

/***  The upper limits - 4GB  ***/
#define MAX_NODES UINT_MAX

/***  Define the ranges  ***/
#define MAX_INDEXES (1 << (IDX_SHF))
#define MAX_BLOCKS (((MAX_NODES) / (MAX_INDEXES) + 1))


/*
***************************
***  Private functions  ***
***************************
*/

/*
**************************************************************
*****                                                    *****
*****  Function: _createNextBlock                        *****
*****                                                    *****
*****  Description:-                                     *****
*****                                                    *****
*****    Allocate the next block for the heap providing  *****
*****    there is no empty blocks and memory is          *****
*****    available.                                      *****
*****                                                    *****
*****  Parameters:-                                      *****
*****                                                    *****
*****    bheap_t* : The heap to allocate the new block.  *****
*****                                                    *****
*****  Return: none                                      *****
*****                                                    *****
**************************************************************
*/

inline void _createNextBlock(bheap_t* theHeap)
{
  theHeap->currentBlock++;

  if(theHeap->emptyBlocks == 0)
  {
/* JRW - 24/Jan/2002 : This test is not needed!!!!
                       Not sure why I did the code this way!
    if(theHeap->currentBlock == 0)
    {
      theHeap->block[theHeap->currentBlock].node =
        (node_t** ) calloc(MAX_INDEXES + 1, sizeof(node_t* ));
      assert(theHeap->block[theHeap->currentBlock].node != NULL);
    }
    else
    {
      theHeap->block[theHeap->currentBlock].node =
        (node_t** ) calloc(MAX_INDEXES, sizeof(node_t* ));
      assert(theHeap->block[theHeap->currentBlock].node != NULL);
    }
*/

	  //GB *
    theHeap->block[theHeap->currentBlock].node =
      (node_t* ) calloc(MAX_INDEXES, sizeof(node_t ));
    assert(theHeap->block[theHeap->currentBlock].node != NULL);
  }
  else
    theHeap->emptyBlocks--;
}



/*
**************************************************************
*****                                                    *****
*****  Function: _getKey                                 *****
*****                                                    *****
*****  Description:-                                     *****
*****                                                    *****
*****    Retrieve the key value.                         *****
*****                                                    *****
*****  Parameters:-                                      *****
*****                                                    *****
*****    node_t*    : Node to get the key from.          *****
*****                                                    *****
*****  Return: KEY_TYPE                                  *****
*****                                                    *****
**************************************************************
*/

//GB external functions
 KEY_TYPE _getKey(node_t theNode);

 INDEX_TYPE _getIndex(node_t theNode);

 void _setIndex(node_t theNode, INDEX_TYPE I);

/*
**************************************************************
*****                                                    *****
*****  Function: _insert                                 *****
*****                                                    *****
*****  Description:-                                     *****
*****                                                    *****
*****    The main 'insert' function. Insert the new      *****
*****    node in the correct location in the heap.       *****
*****                                                    *****
*****  Parameters:-                                      *****
*****                                                    *****
*****    bheap_t*   : The heap to store the new node.    *****
*****    INDEX_TYPE : insert point                       *****
*****    node_t*    : Node to be inserted.               *****
*****                                                    *****
*****  Return: none                                      *****
*****                                                    *****
**************************************************************
*/

inline void _insert(bheap_t* theHeap, INDEX_TYPE _H, node_t theNode)
{
  INDEX_TYPE _I;

  for(_I = (_H - 1) / 2;
      ((0 < _H) &&

/* JRW - 23/Jan/2002
       (theHeap->block[BLOCK(_I)].node[INDEX(_I)]->key >= theNode->key));
*/
       (_getKey(theHeap->block[BLOCK(_I)].node[INDEX(_I)]) < _getKey(theNode)));

      _I = (_H - 1) / 2)
  {
    theHeap->block[BLOCK(_H)].node[INDEX(_H)] =
      theHeap->block[BLOCK(_I)].node[INDEX(_I)];
	_setIndex(theHeap->block[BLOCK(_H)].node[INDEX(_H)] , _H);
//    theHeap->block[BLOCK(_H)].node[INDEX(_H)]->index = _H;

    _H = _I;
  }

  theHeap->block[BLOCK(_H)].node[INDEX(_H)] = theNode;

  _setIndex(theHeap->block[BLOCK(_H)].node[INDEX(_H)] , _H);

//  theHeap->block[BLOCK(_H)].node[INDEX(_H)]->index = _H;
}



/*
**************************************************************
*****                                                    *****
*****  Function: _delete                                 *****
*****                                                    *****
*****  Description:-                                     *****
*****                                                    *****
*****    The main 'delete' function. Delete the node     *****
*****    given by delete point form the heap.            *****
*****                                                    *****
*****  Parameters:-                                      *****
*****                                                    *****
*****    bheap_t*   : The heap to delete node from.      *****
*****    INDEX_TYPE : delete point                       *****
*****                                                    *****
*****  Return: none                                      *****
*****                                                    *****
**************************************************************
*/

inline void _delete(bheap_t* theHeap, INDEX_TYPE _H)
{
  INDEX_TYPE _K = 2 * _H + 2;

  theHeap->nodeCount--;

  for(; _K < theHeap->nodeCount; _K = 2 * _K + 2)
  {
/* JRW - 23/Jan/2002
    if(theHeap->block[BLOCK(_K)].node[INDEX(_K)]->key >=
       theHeap->block[BLOCK(_K - 1)].node[INDEX(_K - 1)]->key)
*/
    if(_getKey(theHeap->block[BLOCK(_K)].node[INDEX(_K)]) <
       _getKey(theHeap->block[BLOCK(_K - 1)].node[INDEX(_K - 1)]))

      _K--;

    theHeap->block[BLOCK(_H)].node[INDEX(_H)] =
      theHeap->block[BLOCK(_K)].node[INDEX(_K)];

    _setIndex(theHeap->block[BLOCK(_H)].node[INDEX(_H)] , _H);

    //theHeap->block[BLOCK(_H)].node[INDEX(_H)]->index = _H;

    _H = _K;
  }

  if(_K == theHeap->nodeCount)
  {
    theHeap->block[BLOCK(_H)].node[INDEX(_H)] =
      theHeap->block[BLOCK(_K - 1)].node[INDEX(_K - 1)];

    _setIndex(theHeap->block[BLOCK(_H)].node[INDEX(_H)] , _H);

    //theHeap->block[BLOCK(_H)].node[INDEX(_H)]->index = _H;

    _H = _K - 1;
  }

  _insert(theHeap, _H,
          theHeap->block[BLOCK(theHeap->nodeCount)].node[INDEX(theHeap->nodeCount)]);
}



/*
**************************
***  Public Functions  ***
**************************
*/


/*
**************************************************************
*****                                                    *****
*****  Function: bh_alloc                                *****
*****                                                    *****
*****  Description:-                                     *****
*****                                                    *****
*****    Initialise a new binary heap.                   *****
*****                                                    *****
*****  Parameters:-                                      *****
*****                                                    *****
*****    none                                            *****
*****                                                    *****
*****  Return: bheap_t* - newly created binary heap.     *****
*****                                                    *****
**************************************************************
*/

inline bheap_t* bh_alloc(void)
{
  bheap_t* theHeap;

  theHeap = (bheap_t* ) malloc(sizeof(bheap_t));
  assert(theHeap != NULL);

  theHeap->block =
    (heap_block_t* ) calloc(MAX_BLOCKS, sizeof(heap_block_t));
  assert(theHeap->block != NULL);

  theHeap->nodeCount = theHeap->emptyBlocks = 0;

/***  The starting point is 0 but set to -1 because the
      _createNextBlock will increament the value
      before using it
***/
  theHeap->currentBlock = (INDEX_TYPE) -1;
  _createNextBlock(theHeap);

  return theHeap;
}



/*
**************************************************************
*****                                                    *****
*****  Function: bh_free                                 *****
*****                                                    *****
*****  Description:-                                     *****
*****                                                    *****
*****    Destroy the binary heap.                        *****
*****                                                    *****
*****  Parameters:-                                      *****
*****                                                    *****
*****    bheap_t* : the heap to be destroy.              *****
*****                                                    *****
*****  Return: none                                      *****
*****                                                    *****
**************************************************************
*/

inline void bh_free(bheap_t* theHeap)
{
  INDEX_TYPE loop;

/*********************************************************/
/*****                                               *****/
/*****  NOTE: unsigned long int is used; therefore,  *****/
/*****        cannot count down to -1!!!!            *****/
/*****                                               *****/
/*********************************************************/

  for(loop = theHeap->currentBlock + theHeap->emptyBlocks;
      loop > 0; loop--)
    free(theHeap->block[loop].node);


/***  To free the first block!  ***/
  free(theHeap->block[0].node);

  free(theHeap->block);
  free(theHeap);
}



/*
**************************************************************
*****                                                    *****
*****  Function: bh_delete                               *****
*****                                                    *****
*****  Description:-                                     *****
*****                                                    *****
*****    Delete the node from the heap.                  *****
*****                                                    *****
*****  Parameters:-                                      *****
*****                                                    *****
*****    bheap_t*   : The heap to delete the new node    *****
*****                 from.                              *****
*****    node_t*    : Node to be deleted.                *****
*****                                                    *****
*****  Return: none                                      *****
*****                                                    *****
**************************************************************
*/

inline void bh_delete(bheap_t* theHeap, INDEX_TYPE theIndex)
{
  _delete(theHeap, theIndex);

  if((theHeap->nodeCount & IDX_MASK) == IDX_MASK)
  {
    theHeap->currentBlock--;
    theHeap->emptyBlocks++;
  }

  //GB
// ????  free(theNode);
}



/*
**************************************************************
*****                                                    *****
*****  Function: bh_delete_min                           *****
*****                                                    *****
*****  Description:-                                     *****
*****                                                    *****
*****    Delete the node that have the minimum key from  *****
*****    the heap.                                       *****
*****                                                    *****
*****  Parameters:-                                      *****
*****                                                    *****
*****    bheap_t*   : The heap to delete the new node    *****
*****                 from.                              *****
*****                                                    *****
*****  Return: DATA_TYPE - the user define data of the   *****
*****                      deleted node.                 *****
*****                      If heap is empty the return   *****
*****                      NULL.                         *****
*****                                                    *****
**************************************************************
*/

inline DATA_TYPE bh_delete_min(bheap_t* theHeap)
{
  if(theHeap->nodeCount > 0)
  {
    node_t theNode = theHeap->block[0].node[0];

    DATA_TYPE data = theNode.data;

    bh_delete(theHeap, _getIndex(theNode) );

    return data;
  }
  else
    return BADINDEX;
}



/*
**************************************************************
*****                                                    *****
*****  Function: bh_insert                               *****
*****                                                    *****
*****  Description:-                                     *****
*****                                                    *****
*****    Build a new node and then add it to the heap.   *****
*****                                                    *****
*****  Parameters:-                                      *****
*****                                                    *****
*****    bheap_t*   : The heap to delete the new node    *****
*****                 from.                              *****
*****    DATA_TYPE  : the user define data               *****
*****    KEY_TYPE   : the key to be sorted by.           *****
*****                 (Now removed.)                     *****
*****                                                    *****
*****  Return: node_t* - the newly created node          *****
*****                                                    *****
**************************************************************
*/

/* JRW - 23/Jan/2002
inline node_t* bh_insert(bheap_t* theHeap, DATA_TYPE data, KEY_TYPE key)
*/
inline INDEX_TYPE bh_insert(bheap_t* theHeap, DATA_TYPE data)
{
//GB
	node_t newNode;

//  newNode = (node_t* ) malloc(sizeof(node_t));
//  assert(newNode != NULL);

  newNode.data = data;

/* JRW - 23/Jan/2002
  newNode->key = key;
*/

  _insert(theHeap, theHeap->nodeCount, newNode);

  theHeap->nodeCount++;

/* JRW - 24/Jan/2002 : This is wrong. Should be as what now given.
  if((theHeap->nodeCount & IDX_POS_MASK) == (MAX_INDEXES - 1))
*/
  if((theHeap->nodeCount & IDX_MASK) == 0)
    _createNextBlock(theHeap);

  return _getIndex(newNode);
}



/*
**************************************************************
*****                                                    *****
*****  Function: bh_return_min                           *****
*****                                                    *****
*****  Description:-                                     *****
*****                                                    *****
*****    Return the user define data of the node that    *****
*****    have the minimum key.                           *****
*****                                                    *****
*****  Parameters:-                                      *****
*****                                                    *****
*****    bheap_t*   : The heap to get the data from.     *****
*****                                                    *****
*****  Return: DATA_TYPE - the user define data          *****
*****                      if heap is empty then return  *****
*****                      NULL                          *****
*****                                                    *****
**************************************************************
*/

inline DATA_TYPE bh_return_min(bheap_t* theHeap)
{
  if(theHeap->nodeCount > 0)
    return theHeap->block[0].node[0].data;
  else
    return BADINDEX;
}


inline DATA_TYPE bh_return_node(bheap_t* theHeap, INDEX_TYPE theIndex)
{
  if(theHeap->nodeCount > theIndex)
    return     theHeap->block[BLOCK(theIndex)].node[INDEX(theIndex)].data;
  else
    return BADINDEX;
}

#endif
