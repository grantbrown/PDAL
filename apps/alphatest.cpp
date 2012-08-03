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




#ifdef PDAL_HAVE_FLANN
#include <flann/flann.hpp>
#endif

#ifdef PDAL_HAVE_GEOS
#include <geos_c.h>


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
    std::string m_wkt;

#ifdef PDAL_HAVE_GEOS
	GEOSContextHandle_t m_geosEnvironment;
#endif
	
    pdal::Options m_options;
};


AlphaShapeQuery::AlphaShapeQuery(int argc, char* argv[])
    : Application(argc, argv, "pcquery")
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
    return;
}

void AlphaShapeQuery::addSwitches()
{
    namespace po = boost::program_options;

    po::options_description* file_options = new po::options_description("file options");
    

    file_options->add_options()
        ("input,i", po::value<std::string>(&m_inputFile)->default_value(""), "input file name")
        ("point", po::value< std::vector<float> >()->multitoken(), "A 2d or 3d point to use for querying")
        ("wkt", po::value<std::string>(&m_wkt)->default_value(""), "WKT object to use for querying")
        ;

    addSwitchSet(file_options);

    po::options_description* processing_options = new po::options_description("processing options");
    
    processing_options->add_options()

        ;
    
    addSwitchSet(processing_options);

    addPositionalSwitch("input", 1);

    return;
}



int AlphaShapeQuery::execute()
{
    using namespace flann;
    using boost::lexical_cast;
    Options readerOptions;
    {
        if (m_usestdin)
            m_inputFile = "STDIN";
        readerOptions.add<std::string>("filename", m_inputFile);
        readerOptions.add<bool>("debug", isDebug());
        readerOptions.add<boost::uint32_t>("verbose", getVerboseLevel());
    }

#ifdef PDAL_HAVE_FLANN

    Stage* stage = AppSupport::makeReader(readerOptions);
    stage -> initialize();
#endif



    pdal::Options options = m_options + readerOptions;
    
    //pdal::filters::Stats* filter = new pdal::filters::Stats(*stage, options);

    //filter->initialize();

    boost::uint64_t numPoints = stage->getNumPoints();
    std::cout << numPoints << std::endl;
    const Schema& schema = stage->getSchema();
    PointBuffer* data = new PointBuffer(schema, 1);
    //PointBuffer data(schema, 0);
    boost::scoped_ptr<StageSequentialIterator>* iter = new boost::scoped_ptr<StageSequentialIterator>(stage->createSequentialIterator(*data));
    //boost::scoped_ptr<StageSequentialIterator>* iter = new boost::scoped_ptr<StageSequentialIterator>(filter->createSequentialIterator(*data));


    int dim = 3;
    std::vector<boost::int32_t> xyz;
    xyz.resize(numPoints*dim);
   
    Dimension const & dimx = schema.getDimension("X");
    Dimension const & dimy = schema.getDimension("Y");
    Dimension const & dimz = schema.getDimension("Z");
    boost::uint32_t itr = 0;

    int xmin;
    int ymin;
    int zmin;
    int xmax;
    int ymax;
    int zmax;
    int _x;
    int _y;
    int _z;
    while (!((**iter).atEnd()))
    {
        (**iter).read(*data);
        xyz[itr] =   (data -> getField<boost::int32_t>(dimx,0));
        xyz[itr+1] = (data -> getField<boost::int32_t>(dimy,0));
        xyz[itr+2] = (data -> getField<boost::int32_t>(dimz,0));
        if (itr == 0)
        {
            xmin = xyz[itr];
            xmax = xmin; 
            ymin = xyz[itr+1];
            ymax = ymin;
            zmin = xyz[itr + 2];
            zmax = zmin;
        }
        else
        {
            _x = xyz[itr];
            _y = xyz[itr + 1];
            _z = xyz[itr + 2];

            if (_x > xmax){xmax = _x;}
            else if (_x < xmin){xmin = _x;}
            if (_y > ymax){ymax = _y;}
            else if (_y < ymin){ymin = _y;}
            if (_z > zmax){zmax = _z;}
            else if (_z < zmin){zmin = _z;} 
        }
        itr ++;
    }
    std::cout << "Data Read" << std::endl;

    std::cout << "X: " << xmin << ", " << xmax << std::endl;
    std::cout << "Y: " << ymin << ", " << ymax << std::endl;
    std::cout << "Z: " << zmin << ", " << zmax << std::endl;

    //boost::property_tree::ptree stats_tree = static_cast<pdal::filters::iterators::sequential::Stats*>(iter->get())->toPTree();
    //boost::property_tree::ptree tree;
    //tree.add_child("stats", stats_tree);
    //write_xml(std::cout, tree);
    //std::cout << stats_tree.get<double>("X.minimum") << std::endl;    

    
    //delete[] &xyz;
    std::cout << std::endl;
    

    delete iter;
    //delete filter;
    delete stage;

    
    return 0;
}


int main(int argc, char* argv[])
{
    AlphaShapeQuery app(argc, argv);
    return app.run();
}

