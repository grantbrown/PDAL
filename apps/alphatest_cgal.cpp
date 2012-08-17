/******************************************************************************
* Copyright (c) 2012, Howard Butler (hobu.inc@gmail.com)
*
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following
* conditions are met:
*
*     * Redistributions of source code must retain the above copyright
*       notice, this list of conditions and the following disclaimer.
*     * Redistributions in binary form must reproduce the above copyright
*       notice, this list of conditions and the following disclaimer in
*       the documentation and/or other materials provided
*       with the distribution.
*     * Neither the name of Hobu, Inc. or Flaxen Geo Consulting nor the
*       names of its contributors may be used to endorse or promote
*       products derived from this software without specific prior
*       written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
* "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
* LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
* FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
* COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
* INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
* BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
* OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
* AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
* OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
* OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY
* OF SUCH DAMAGE.
****************************************************************************/

#include <iostream>
#include <fstream>
#include <boost/scoped_ptr.hpp>

#include <pdal/Stage.hpp>
#include <pdal/StageIterator.hpp>
#include <pdal/FileUtils.hpp>
#include <pdal/PointBuffer.hpp>
#include <pdal/filters/Index.hpp>
#include <pdal/filters/Stats.hpp>

#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/json_parser.hpp>

#include <boost/lexical_cast.hpp>

#include "AppSupport.hpp"
#include "Application.hpp"
#include <cstdarg>
#include <math.h>
#include <tr1/unordered_map>
#include "SparseGrid.hpp"

#ifdef PDAL_HAVE_GEOS
#include <geos_c.h>

#include <iomanip>


#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/convex_hull_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;


//Rip out GEOS stuff
namespace alphatest_cgal
{
    static void _GEOSErrorHandler(const char *fmt, ...)
    {
        va_list args;

        va_start(args, fmt);
        char buf[1024];  

        vsnprintf( buf, sizeof( buf), fmt, args);
        std::cout << "GEOS Error: " << buf << std::endl;

        va_end(args);
    }

    static void _GEOSWarningHandler(const char *fmt, ...)
    {
        va_list args;

        char buf[1024];  
        vsnprintf( buf, sizeof( buf), fmt, args);
        std::cout << "GEOS warning: " << buf << std::endl;

        va_end(args);
    }

} // geos

#endif

using namespace pdal;


class AlphaShapeQuery : public Application
{
public:
    AlphaShapeQuery(int argc, char* argv[]);
    int execute(); // overrride

private:
    void addSwitches(); // overrride
    void validateSwitches(); // overrride

    boost::uint64_t numPoints;

    std::string m_inputFile;
    std::string m_outputFile;
    std::string m_wkt;
    float m_Alpha;
    int m_subset;

#ifdef PDAL_HAVE_GEOS
	GEOSContextHandle_t m_geosEnvironment;
#endif
	
    pdal::Options m_options;
};


AlphaShapeQuery::AlphaShapeQuery(int argc, char* argv[])
    : Application(argc, argv, "alphatest")
    , m_inputFile("")
{
    return;
}


void AlphaShapeQuery::validateSwitches()
{
    
    if (m_inputFile == "")
    {
        throw app_usage_error("--input/-i required");
    }
    if (m_outputFile == "")
    {
        throw app_usage_error("--output/-o required");
    }
    return;
}

void AlphaShapeQuery::addSwitches()
{
    namespace po = boost::program_options;

    po::options_description* file_options = new po::options_description("file options");
    
    file_options->add_options()
        ("input,i", po::value<std::string>(&m_inputFile)->default_value(""), "input file name")
        ("output,o", po::value<std::string>(&m_outputFile)->default_value(""), "output file name")
        ("point", po::value< std::vector<float> >()->multitoken(), "A 2d or 3d point to use for querying")
        ("wkt", po::value<std::string>(&m_wkt)->default_value(""), "WKT object to use for querying")
        ("alpha,a", po::value<float>(&m_Alpha)->default_value(0.5), "Alpha Parameter")
        ("nsub,n", po::value<int>(&m_subset)->default_value(8), "How many subset operations?");

    addSwitchSet(file_options);
    po::options_description* processing_options = new po::options_description("processing options");
    processing_options->add_options();
    addSwitchSet(processing_options);
    addPositionalSwitch("input", 1);
    addPositionalSwitch("output", 2);


    return;
}

int AlphaShapeQuery::execute()
{

    using boost::lexical_cast;
    Options readerOptions;
    {
        if (m_usestdin)
            m_inputFile = "STDIN";
        readerOptions.add<std::string>("filename", m_inputFile);
        readerOptions.add<bool>("debug", isDebug());
        readerOptions.add<boost::uint32_t>("verbose", getVerboseLevel());
    }

    Stage* stage = AppSupport::makeReader(readerOptions);


    

    stage -> initialize();
    int pointbuffersize = 1000;

    pdal::Options options = m_options + readerOptions;
    
    boost::uint64_t numPoints = stage->getNumPoints();
    const Schema& stage_schema = stage->getSchema();
    PointBuffer* data = new PointBuffer(stage_schema, pointbuffersize);
    boost::scoped_ptr<StageRandomIterator>* iter = new boost::scoped_ptr<StageRandomIterator>(stage->createRandomIterator(*data));
    const Schema& buffer_schema = data->getSchema();

    int dim = 3;
    Dimension const & dimx = buffer_schema.getDimension("X");
    Dimension const & dimy = buffer_schema.getDimension("Y");
    Dimension const & dimz = buffer_schema.getDimension("Z");

    int xmin;
    int ymin;
    int zmin;
    int xmax;
    int ymax;
    int zmax;
    int _x;
    int _y;
    int _z;    
    bool first = true;
    std::cout << "Reading Data and Calculating Extremes" << std::endl;
    boost::uint64_t itrs = 0;
    while (itrs < numPoints)
    {
        (**iter).read(*data);

        for (boost::uint32_t i = 0; i < (data -> getNumPoints()); i ++)
        {
            itrs += 1;
            _x = (data -> getField<boost::int32_t>(dimx,i));
            _y = (data -> getField<boost::int32_t>(dimy,i));

            if (!first)
            {
                if (_x > xmax){xmax = _x;}
                if (_x < xmin){xmin = _x;}
                if (_y > ymax){ymax = _y;}
                if (_y < ymin){ymin = _y;}
            }
        }
        if (first)
        {
            first = false;
            xmin = _x; 
            xmax = xmin; 
            ymin = _y; 
            ymax = ymin;
        }
    }
    std::cout << "Extremes:" << std::endl;
    std::cout << "X: " << xmin << ", " << xmax << std::endl;
    std::cout << "Y: " << ymin << ", " << ymax << std::endl;

    std::cout << "Building Grid:" << std::endl;
    SparseGrid* grid = new SparseGrid(xmin, ymin, 
                                     xmax, ymax,
                                     numPoints, 10000);

    (**iter).seek(0); 

    itrs = 0;
    while (itrs < numPoints)
    {
        (**iter).read(*data);

        for (boost::uint32_t i = 0; i < (data -> getNumPoints()); i ++)
        {
            grid -> insertPoint(
            (data -> getField<boost::int32_t>(dimx,i)),
            (data -> getField<boost::int32_t>(dimy,i)), itrs);
            itrs += 1;
        }
    }
    std::cout << "Grid Built, subsetting." << std::endl;
    
    boost::uint64_t cells = 0;
    for (int itr = 0; itr < m_subset; itr ++)
    {
        cells = (150 + 100*itr)*(150 + 100*itr);
        grid -> subset_and_regrid(cells); 
    }


    //Find alpha heuristically    
    int xrange = xmax-xmin;
    int yrange = ymax-ymin;
    int num = (xrange > yrange ? xrange : yrange);
    int denom = (xrange <= yrange ? xrange : yrange);

    double ratio = static_cast<double>(num)/denom;
    boost::uint64_t area = (xrange)*(yrange); 
    double density = static_cast<double>(numPoints)/area;
    
    std::cout << "num = " << num << std::endl;
    std::cout << "denom = " << denom << std::endl;
    std::cout << "ratio = " << ratio << std::endl;
    std::cout << "density = " << density << std::endl;

    double  CalcAlpha = (70000 + 0.04838 * num - (0.0497*denom) - (17490 * ratio));

    std::cout << "Alpha Calculated to be: " << CalcAlpha << std::endl;
    if (CalcAlpha < 0)
    {
        std::cout << "WARNING, negative alpha calculated. Using dummy value." << std::endl;
        CalcAlpha = 50000;
    }

    std::cout << "Xrange: " << (xmax - xmin) << std::endl;
    std::cout << "Yrange: " << (ymax - ymin) << std::endl;
    std::cout << "Num Points: " << numPoints << std::endl;

    std::cout << "Getting New Valid Points" << std::endl; 
    std::stack<boost::uint64_t>* goodpoints = grid -> getValidPointIdx();

    PointBuffer* outdata = new PointBuffer(buffer_schema, 1);
    const Schema& buffer_schema_2 = outdata->getSchema();

    Dimension const & dimx2 = buffer_schema_2.getDimension("X");
    Dimension const & dimy2 = buffer_schema_2.getDimension("Y");
    Dimension const & dimz2 = buffer_schema_2.getDimension("Z");
    
    boost::uint64_t good_point = 0;
    int itr = 0;


    while (!(goodpoints -> empty()))
    {
        good_point = goodpoints -> top(); 

        (**iter).seek(good_point); 
        (**iter).read(*outdata);
        _x = (outdata -> getField<boost::int32_t>(dimx2,0));
        _y = (outdata -> getField<boost::int32_t>(dimy2,0));

        // Dump stuff into CGAL here. 


        itr += 1;
        goodpoints -> pop();
    }



    // Do magic here.

    std::cout << "Writing to: " << m_outputFile << std::endl;
    std::ofstream outfile;
    char* filepath = new char[m_outputFile.size() + 1];
    filepath[m_outputFile.size() + 1] = 0;
    memcpy(filepath, m_outputFile.c_str(), m_outputFile.size());
    outfile.open(filepath);
    outfile <<  std::fixed << std::setprecision(0);

    // Write WKT output here. 

    outfile.close();

    delete filepath;
    delete data;
    delete outdata;
    delete grid;
    delete stage;

    return 0;
}


int main(int argc, char* argv[])
{
    AlphaShapeQuery app(argc, argv);
    return(app.run());
}

