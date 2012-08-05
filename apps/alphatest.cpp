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
#include <math.h>




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
    int pointbuffersize = 1000;
#endif



    pdal::Options options = m_options + readerOptions;
    
    //pdal::filters::Stats* filter = new pdal::filters::Stats(*stage, options);

    //filter->initialize();

    boost::uint64_t numPoints = stage->getNumPoints();
    const Schema& schema = stage->getSchema();
    PointBuffer* data = new PointBuffer(schema, pointbuffersize);
    //PointBuffer data(schema, 0);
    boost::scoped_ptr<StageSequentialIterator>* iter = new boost::scoped_ptr<StageSequentialIterator>(stage->createSequentialIterator(*data));
    //boost::scoped_ptr<StageSequentialIterator>* iter = new boost::scoped_ptr<StageSequentialIterator>(filter->createSequentialIterator(*data));



    int dim = 3;
    std::vector<boost::int32_t> xyz;
    xyz.resize(numPoints*dim);
    for (boost::uint32_t i = 0; i < numPoints*dim; i ++)
    {
        xyz[i] = 9999;
    }
   
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
    
    bool first = true;
    std::cout << "Reading Data and Calculating Extremes" << std::endl;
    int k;
    int idx1 = 0;
    int idx2;
    int idx3;
    while (!((**iter).atEnd()))
    {
        (**iter).read(*data);

        for (boost::uint32_t i = 0; i < (data -> getNumPoints()); i ++)
        {
            _x = (data -> getField<boost::int32_t>(dimx,i));
            _y = (data -> getField<boost::int32_t>(dimy,i));
            _z = (data -> getField<boost::int32_t>(dimz,i));

            idx2 = idx1 + 1;
            idx3 = idx2 + 1;

            xyz[idx1] = _x; 
            xyz[idx2] = _y;
            xyz[idx3] = _z;

            if ((xyz[idx1] == 9999) ||  (xyz[idx2] == 9999) || (xyz[idx3] == 9999))
            {
                std::cout << "Bad Point: " << std::endl;
                std::cout << "Index: " << idx1 << ", " << idx2 << ", " << idx3 << std::endl;
                std::cout << xyz[idx1] << ", " << xyz[idx2] << ", " << xyz[idx3] << std::endl;
                std::cout << _x << ", " << _y << ", " <<  _z << std::endl; 
            }
            if (!first)
            {
                _x = xyz[idx1];
                _y = xyz[idx2];
                _z = xyz[idx3]; 
                if (_x > xmax){xmax = _x;}
                else if (_x < xmin){xmin = _x;}
                if (_y > ymax){ymax = _y;}
                else if (_y < ymin){ymin = _y;}
                if (_z > zmax){zmax = _z;}
                else if (_z < zmin){zmin = _z;} 
            }
            idx1 += 3;
        }

        if (first)
        {
            first = false;
            xmin = xyz[0];
            xmax = xmin; 
            ymin = xyz[1];
            ymax = ymin;
            zmin = xyz[2];
            zmax = zmin;

        }
    }
    std::cout << "Total Iterations: " << idx1 << std::endl;
    std::cout << "Extremes:" << std::endl;
    std::cout << "X: " << xmin << ", " << xmax << std::endl;
    std::cout << "Y: " << ymin << ", " << ymax << std::endl;
    std::cout << "Z: " << zmin << ", " << zmax << std::endl;

    std::cout << "Building Bins" << std::endl;
    double xdensity = (static_cast<double>(numPoints))/(static_cast<double>(xmax-xmin));
    double ydensity = (static_cast<double>(numPoints))/(static_cast<double>(ymax-ymin));
    double zdensity = (static_cast<double>(numPoints))/(static_cast<double>(zmax-zmin));

    std::cout << "X Density: " << xdensity << std::endl;
    std::cout << "Y Density: " << ydensity << std::endl;
    std::cout << "Z Density: " << zdensity << std::endl;

    double xinterval = 10000/xdensity;
    double yinterval = 10000/ydensity;
    double zinterval = 10000/zdensity;

    std::cout << "...x interval: " << xinterval << std::endl;
    std::cout << "...y interval: " << yinterval << std::endl;
    std::cout << "...z interval: " << zinterval << std::endl;
 
    int xbins = (xmax-xmin)/(xinterval);
    int ybins = (xmax-xmin)/(yinterval);
    int zbins = (xmax-xmin)/(zinterval);

    std::cout << "...x bins: " << xbins << std::endl;
    std::cout << "...y bins: " << ybins << std::endl;
    std::cout << "...z bins: " << zbins << std::endl;

    
    std::vector<boost::uint32_t> point_bin_membership;
    point_bin_membership.resize(numPoints);
    
    std::vector<boost::uint32_t> point_bins;
    point_bins.resize(xbins * ybins * zbins);
    std::cout << "...Initializing Bins" << std::endl;
    for (int i = 0; i < xbins*ybins*zbins; i ++)
    {
        point_bins[i] = 0;
    }
    std::cout << "...Bins Initialized" << std::endl;

    boost::uint32_t _xidx;
    boost::uint32_t _yidx;
    boost::uint32_t _zidx;
    boost::uint32_t idx;

    std::cout << "Filling Bins, numPoints = " << numPoints << std::endl;
    for (boost::uint32_t i = 0; i<numPoints; i++)
    {
        _x = xyz[i*3];
        _y = xyz[i*3 + 1];
        _z = xyz[i*3 + 2];
        _xidx = static_cast<boost::uint32_t>(floor((_x - xmin)/xinterval)); 
        _yidx = static_cast<boost::uint32_t>(floor((_y - ymin)/yinterval)); 
        _zidx = static_cast<boost::uint32_t>(floor((_z - zmin)/zinterval)); 
        idx = _xidx + _yidx*xbins + _zidx*xbins*ybins;
        if (idx > xbins*ybins*zbins)
        {
            std::cout << "Bad Index: " << idx << std::endl;
            std::cout << "...X Index: " << _xidx << std::endl;
            std::cout << "...Y Index: " << _yidx << std::endl;
            std::cout << "...Z Index: " << _zidx << std::endl;
            std::cout << "...X: " << _x << std::endl;
            std::cout << "...Y: " << _y << std::endl;
            std::cout << "...Z: " << _x << std::endl;
            std::cout << "...Point Number: " << i << std::endl;


        }
        point_bin_membership[i] = idx;
        point_bins[idx] += 1;

    }
    std::cout << "Bins Populated." << std::endl;






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

