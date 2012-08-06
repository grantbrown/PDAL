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
        SparseGrid(int _xmin, int _ymin,
                int _xmax, int _ymax,
                int _nPoints, int _dim_width = 50);

        int insertPoint(boost::int32_t X, boost::int32_t Y, boost::int64_t point_idx);
        int getIndex(int xidx, int yidx);
        bool isValid(int xidx, int yidx);
        std::stack<boost::uint64_t>* getValidPoints();

        std::vector<SparseGridNode*>* grid; 

        int xmin; int ymin; 
        int xmax; int ymax; 
        int xbins; int ybins; 
        int numPoints; int dim_width;
        double xdensity; double ydensity;
        double xinterval; double yinterval;
    private:
        int initializeGrid();
    
};
