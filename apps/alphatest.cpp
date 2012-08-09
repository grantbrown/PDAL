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
#include "SparseGrid.hpp"

#ifdef PDAL_HAVE_GEOS
#include <geos_c.h>


#include <pcl/ModelCoefficients.h>
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/sample_consensus/method_types.h>
#include <pcl/sample_consensus/model_types.h>
#include <pcl/filters/passthrough.h>
#include <pcl/segmentation/sac_segmentation.h>
#include <pcl/surface/concave_hull.h>









namespace alphatest
{
    static void _GEOSErrorHandler(const char *fmt, ...)
    {
        va_list args;

        va_start(args, fmt);
        char buf[1024];  

        vsnprintf( buf, sizeof( buf), fmt, args);
        std::cerr << "GEOS Error: " << buf << std::endl;

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
        ;

    addSwitchSet(file_options);

    po::options_description* processing_options = new po::options_description("processing options");
    
    processing_options->add_options()

        ;
    
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
    /*
    Options writerOptions;
    {
        writerOptions.add<std::string>("filename", m_outputFile);
        writerOptions.add<bool>("debug", isDebug());
        writerOptions.add<boost::uint32_t>("verbose", getVerboseLevel());
    }
    */

    Stage* stage = AppSupport::makeReader(readerOptions);


    
    //writer -> initialize();
    stage -> initialize();
    //Writer* writer = AppSupport::makeWriter(readerOptions, *stage);
    int pointbuffersize = 1000;




    pdal::Options options = m_options + readerOptions;
    
    boost::uint64_t numPoints = stage->getNumPoints();
    const Schema& schema = stage->getSchema();
    PointBuffer* data = new PointBuffer(schema, pointbuffersize);
    boost::scoped_ptr<StageRandomIterator>* iter = new boost::scoped_ptr<StageRandomIterator>(stage->createRandomIterator(*data));

    int dim = 3;
    Dimension const & dimx = schema.getDimension("X");
    Dimension const & dimy = schema.getDimension("Y");
    Dimension const & dimz = schema.getDimension("Z");

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
            //_x = (data -> getField<boost::int32_t>(dimx,i));
            //_y = (data -> getField<boost::int32_t>(dimy,i));

            grid -> insertPoint(
            (data -> getField<boost::int32_t>(dimx,i)),
            (data -> getField<boost::int32_t>(dimy,i)), itrs);

            //grid -> insertPoint(memcpy(_x), memcpy(_y), itrs);
            itrs += 1;
        }
    }
    std::cout << "Grid Built, subsetting." << std::endl;
    

    grid -> subset_and_regrid(200*200);
    grid -> subset_and_regrid(500*500);
    grid -> subset_and_regrid(1000*1000);



    std::cout << "Getting New Valid Points" << std::endl; 
    std::stack<boost::uint64_t>* goodpoints = grid -> getValidPointIdx();

    PointBuffer* outdata = new PointBuffer(schema, 1);

    //test code:
    
    boost::uint64_t good_point = 0;
    std::ofstream outfile;
    outfile.open("keep_indices.txt");
    int itr = 0;
    pcl::PointCloud<pcl::PointXYZ> cloud;
    cloud.points.resize(goodpoints -> size());
    while (!(goodpoints -> empty()))
    {

        //std::cout << "Writing Point: " << itr << std::endl;
        good_point = goodpoints -> top(); 
        outfile << good_point << "\n";
        (**iter).seek(good_point); 
        (**iter).read(*outdata);
        _x = (data -> getField<boost::int32_t>(dimx,itr));
        _y = (data -> getField<boost::int32_t>(dimy,itr));
        _z = (data -> getField<boost::int32_t>(dimz,itr)); 
        cloud.points[itr].x = _x;
        cloud.points[itr].y = _y;
        cloud.points[itr].z = _z;
        itr += 1;
        //writer -> write();
        goodpoints -> pop();
    }


    
    
    //boost::property_tree::ptree stats_tree = static_cast<pdal::filters::iterators::sequential::Stats*>(iter->get())->toPTree();
    //boost::property_tree::ptree tree;
    //tree.add_child("stats", stats_tree);
    //write_xml(std::cout, tree);
    //std::cout << stats_tree.get<double>("X.minimum") << std::endl;    

    
    //delete[] &xyz;
    std::cout << std::endl;

    delete iter;
    delete stage;
    return 0;
}


int main(int argc, char* argv[])
{
    AlphaShapeQuery app(argc, argv);
    return app.run();
}

