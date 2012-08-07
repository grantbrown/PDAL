#include <stdio.h>
#include <iostream>
#include <stack>
#include <vector>
#include <boost/cstdint.hpp>


struct PointData
{
    boost::int32_t x;
    boost::int32_t y;
    boost::int64_t pt_idx;
};

struct SparseGridNode
{
    SparseGridNode(int _count);
    ~SparseGridNode();
    int count;
    std::stack<PointData*>* point_stack;
};

class SparseGrid 
{
    public:
        SparseGrid(int _xmin, int _ymin,
                int _xmax, int _ymax,
                boost::uint64_t _nPoints, int _dim_width = 50);
        ~SparseGrid();

        int insertPoint(boost::int32_t X, boost::int32_t Y, boost::int64_t point_idx);
        int subset_and_regrid(int _dim_width);
        int set_bins();

        int getIndex(int xidx, int yidx);
        bool isValid(int xidx, int yidx);
        std::stack<SparseGridNode*>* getValidPointRefs(boost::uint64_t* count_ref);
        std::stack<boost::uint64_t>* getValidPointIdx();

        std::vector<SparseGridNode*>* grid; 

        int xmin; int ymin; 
        int xmax; int ymax; 
        int xbins; int ybins; 
        boost::uint64_t numPoints; int dim_width;
        double xdensity; double ydensity;
        double xinterval; double yinterval;
    private:
        int initializeGrid();
    
};
