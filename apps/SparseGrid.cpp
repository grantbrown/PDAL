#include <SparseGrid.hpp>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <boost/cstdint.hpp>
#include <cstdarg>
#include <stack>
#include <algorithm>
#include <math.h>


//Sparse Grid Implementation


SparseGrid::SparseGrid(int _xmin, int _ymin,
                       int _xmax, int _ymax,
                       boost::uint64_t _nPoints, int _dim_width) 
{
    xmin = _xmin;
    xmax = _xmax;
    ymin = _ymin;
    ymax = _ymax;
    
    numPoints = _nPoints;
    dim_width = _dim_width;
    
    set_bins();
    initializeGrid();
}

int SparseGrid::set_bins()
{
    xdensity = (static_cast<double>(numPoints))/(static_cast<double>(xmax-xmin));
    ydensity = (static_cast<double>(numPoints))/(static_cast<double>(ymax-ymin));
    
    std::cout << "X,YDensity: " << xdensity << ", " << ydensity <<  std::endl;

    xinterval = (numPoints/dim_width)/xdensity;
    yinterval = (numPoints/dim_width)/ydensity;

    std::cout << "X,Y Interval: " << xinterval << ", " << yinterval << std::endl;

    xbins = (xmax - xmin)/(xinterval);    
    ybins = (ymax - ymin)/(yinterval);    

    std::cout << "X,Y Bins: " << xbins << ", " << ybins << std::endl;
    return(0);
}

int SparseGrid::initializeGrid()
{
    grid = new std::vector<SparseGridNode*>;
    int gridsize = xbins*ybins;
    grid -> resize(gridsize);
    // Initialize Grid
    std::cout << "Initializing Grid" << std::endl;
    for (int i =0; i < gridsize; i++)
    {
        //*i = new SparseGridNode(0);
        (*grid)[i] = new SparseGridNode(0);
    }
    std::cout << "Grid Initialized" << std::endl;
    return(0);
}

int SparseGrid::subset_and_regrid(int _dim_width)
{
    boost::uint64_t* newsize = new boost::uint64_t;
    std::stack<SparseGridNode*>* good_nodes = getValidPointRefs(newsize);
    //for (int i = 0; i < xbins*ybins; i++)
    //{
    //    delete (*grid)[i];
    //}
    numPoints = *newsize;
    dim_width = _dim_width;
    set_bins();
    
    initializeGrid();
    PointData* pt_dat;
    SparseGridNode* grd_node;
    boost::int32_t x;
    boost::int32_t y;
    boost::int64_t pt_idx;

    int itrs = 0;
    //Reinserting points...
    std::cout << "Reinserting Points" << std::endl;
    while (!(good_nodes -> empty()))
    {
        grd_node = (good_nodes -> top());
        //std::cout << "Empty:" << ((grd_node -> point_stack) -> empty()) << std::endl;
        //std::cout << "Not Empty:" << !((grd_node -> point_stack) -> empty()) << std::endl;
        while (!((grd_node -> point_stack) -> empty()))
        {
            //std::cout << itrs << std::endl;
            pt_dat = ((grd_node -> point_stack) -> top());
            x = pt_dat -> x;
            y = pt_dat -> y;
            pt_idx = pt_dat -> pt_idx;
            insertPoint(x,y,pt_idx);
            ((grd_node -> point_stack) ->  pop());
        }
        itrs += 1;
        (good_nodes -> pop());
    }
    std::cout << "Points Inserted" << std::endl;
    return(0);
}

int SparseGrid::getIndex(int xidx, int yidx)
{
    return(xidx + yidx*xbins);
}
bool SparseGrid::isValid(int xidx, int yidx)
{
    for (int x = -1; x <= 1; x += 1)
    {
        for (int y = -1; y <= 1; y += 1)
        {
            if (x == 0 && y == 0)
            {
                continue;
            }
            if ((xidx + x < 0) || (xidx + x >= xbins) ||
                (yidx + y < 0) || (yidx + y >= ybins))
            {
                return(true);
            }
            if ((((*grid)[getIndex(xidx + x, yidx + y)]) -> count) == 0)
            {
                return(true);
            }
        }
    }
    return(false);
}

int SparseGrid::insertPoint(boost::int32_t X, boost::int32_t Y, boost::int64_t point_idx)
{
    int _xidx = static_cast<boost::uint32_t>(floor((X-xmin)/xinterval));
    _xidx = _xidx < xbins ? _xidx : xbins-1;
    int _yidx = static_cast<boost::uint32_t>(floor((Y-ymin)/yinterval));
    _yidx = _yidx < ybins ? _yidx : ybins-1;
    int idx = getIndex(_xidx, _yidx);
    SparseGridNode* gridnode = (*grid)[idx];
    if (isValid(_xidx, _yidx))
    {
        PointData* newpt = new PointData;
        newpt -> x = X; newpt -> y = Y; newpt -> pt_idx = point_idx;
        gridnode -> count ++;
        gridnode -> point_stack -> push(newpt);
    }
    return(0);
}

std::stack<SparseGridNode*>* SparseGrid::getValidPointRefs(boost::uint64_t* count_ref)
{
    *count_ref = 0;
    std::stack<SparseGridNode*>* outstack = new std::stack<SparseGridNode*>;
    for (int x = 0; x < xbins; x++)
    {
    for (int y = 0; y < ybins; y++)
    {
            if (isValid(x,y))
            {
                SparseGridNode* gridnode = (*grid)[getIndex(x,y)];
                outstack -> push(gridnode);
                (*count_ref) += ((gridnode -> point_stack) -> size());
            }

        }
    }
    return(outstack);
}

std::stack<boost::uint64_t>* SparseGrid::getValidPointIdx()
{
    PointData* idx;
    std::stack<boost::uint64_t>* outstack = new std::stack<boost::uint64_t>;
    for (int x = 0; x < xbins; x++)
    {
    for (int y = 0; y < ybins; y++)
    {
            if (isValid(x,y))
            {
                SparseGridNode* gridnode = (*grid)[getIndex(x,y)];
                while (!((gridnode -> point_stack) -> empty()))
                {
                    idx = (gridnode -> point_stack) -> top();
                    outstack -> push(idx -> pt_idx);
                    (gridnode -> point_stack) -> pop();
                }
            }

        }
    }
    return(outstack);
}

SparseGrid::~SparseGrid()
{
    for (int i = 0; i < grid -> size(); i ++)
    {
        delete (*grid)[i];
    }
}

//SparseGridNode Implementation

SparseGridNode::SparseGridNode(int _count)
{
    count = _count;
    point_stack = new std::stack<PointData*>;
}

SparseGridNode::~SparseGridNode()
{
    delete point_stack;
}
