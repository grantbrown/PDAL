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


SparseGrid::SparseGrid(int _xmin, int _ymin, int _zmin,
                       int _xmax, int _ymax, int _zmax,
                       int _nPoints) 
{
    xmin = _xmin;
    xmax = _xmax;
    ymin = _ymin;
    ymax = _ymax;
    zmin = _zmin;
    zmax = _zmax;
    
    numPoints = _nPoints;

    xdensity = (static_cast<double>(numPoints))/(static_cast<double>(xmax-xmin));
    ydensity = (static_cast<double>(numPoints))/(static_cast<double>(ymax-ymin));
    zdensity = (static_cast<double>(numPoints))/(static_cast<double>(zmax-zmin));
    
    std::cout << "X,Y,Z Density: " << xdensity << ", " << ydensity << ", " << zdensity << std::endl;

    xinterval = (numPoints/100)/xdensity;
    yinterval = (numPoints/100)/ydensity;
    zinterval = (numPoints/100)/zdensity;

    std::cout << "X,Y,Z Interval: " << xinterval << ", " << yinterval << ", " << zinterval << std::endl;


    xbins = (xmax - xmin)/(xinterval);    
    ybins = (ymax - ymin)/(yinterval);    
    zbins = (zmax - zmin)/(zinterval);    

    std::cout << "X,Y,Z Bins: " << xbins << ", " << ybins << ", " << zbins << std::endl;


    initializeGrid();

}

int SparseGrid::initializeGrid()
{
    // This is not very efficient, but the grid size is much smaller than 
    // the point size.
    grid = new std::vector<SparseGridNode*>;
    int gridsize = xbins*ybins*zbins;
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

int SparseGrid::getIndex(int xidx, int yidx, int zidx)
{
    return(xidx + yidx*xbins + zidx*xbins*ybins);
}
bool SparseGrid::isValid(int xidx, int yidx, int zidx)
{
    for (int x = -1; x <= 1; x += 1)
    {
        for (int y = -1; y <= 1; y += 1)
        {
            for (int z = -1; z <= 1; z += 1)
            {
                if (x == 0 && y == 0 && z == 1)
                {
                    continue;
                }
                if ((xidx + x < 0) || (xidx + x >= xbins) ||
                    (yidx + y < 0) || (yidx + y >= ybins) ||
                    (zidx + z < 0) || (zidx + z >= zbins))
                {
                    return(true);
                }
                if ((((*grid)[getIndex(xidx + x, yidx + y, zidx + z)]) -> count) < 8)
                {
                    return(true);
                }
            }
        }
    }
    return(false);
}

int SparseGrid::insertPoint(boost::int32_t X, boost::int32_t Y, boost::int32_t Z, boost::int64_t point_idx)
{
    int _xidx = static_cast<boost::uint32_t>(floor((X-xmin)/xinterval));
    _xidx = _xidx < xbins ? _xidx : xbins-1;
    int _yidx = static_cast<boost::uint32_t>(floor((Y-ymin)/yinterval));
    _yidx = _yidx < ybins ? _yidx : ybins-1;
    int _zidx = static_cast<boost::uint32_t>(floor((Z-zmin)/zinterval));
    _zidx = _zidx < zbins ? _zidx : zbins-1;
    int idx = getIndex(_xidx, _yidx, _zidx);
    SparseGridNode* gridnode = (*grid)[idx];
    if (isValid(_xidx, _yidx, _zidx))
    {
        gridnode -> count ++;
        gridnode -> point_stack -> push(point_idx);
    }
    return(0);
}

int SparseGrid::getValidPoints(int* goodPointIndices)
{
    return(1);
}

//SparseGridNode Implementation

SparseGridNode::SparseGridNode(int _count)
{
    count = _count;
    point_stack = new std::stack<boost::int64_t>;
}
