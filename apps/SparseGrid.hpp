#include <stdio.h>
#include <iostream>
#include <stack>
#include <vector>
#include <boost/cstdint.hpp>




struct SparseGridNode
{
    SparseGridNode(int _count);
    int count;
    std::stack<boost::int64_t>* point_stack;
};

class SparseGrid 
{
    public:
        SparseGrid(int _xmin, int _ymin, int _zmin,
                int _xmax, int _ymax, int _zmax,
                int _nPoints);

        int insertPoint(boost::int32_t X, boost::int32_t Y, boost::int32_t Z, boost::int64_t point_idx);
        int getIndex(int xidx, int yidx, int zidx);
        bool isValid(int xidx, int yidx, int zidx);
        int getValidPoints(int* goodPointIndices);

        std::vector<SparseGridNode*>* grid; 

        int xmin; int ymin; int zmin; 
        int xmax; int ymax; int zmax;
        int xbins; int ybins; int zbins;
        int numPoints;
        double xdensity; double ydensity; double zdensity;
        double xinterval; double yinterval; double zinterval;
    private:
        int initializeGrid();
    
};
