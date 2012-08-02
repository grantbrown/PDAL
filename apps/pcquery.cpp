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

#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/json_parser.hpp>

#include <boost/lexical_cast.hpp>

#include "AppSupport.hpp"
#include "Application.hpp"
#include <cstdarg>

#ifdef PDAL_HAVE_FLANN
#include <flann/flann.h>
#endif

#ifdef PDAL_HAVE_GEOS
#include <geos_c.h>


extern "C" 
{
        struct FLANNParameters DEFAULT_FLANN_PARAMETERS;
}

namespace pcquery
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


class PcQuery : public Application
{
public:
    PcQuery(int argc, char* argv[]);
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


PcQuery::PcQuery(int argc, char* argv[])
    : Application(argc, argv, "pcquery")
    , m_inputFile("")
{
    return;
}


void PcQuery::validateSwitches()
{
    
    if (m_inputFile == "")
    {
        throw app_usage_error("--input/-i required");
    }
    if (m_wkt.size())
    {
#ifdef PDAL_HAVE_GEOS
        // read the WKT using GEOS
        m_geosEnvironment = initGEOS_r(pcquery::_GEOSWarningHandler, pcquery::_GEOSErrorHandler);
        GEOSWKTReader* reader = GEOSWKTReader_create_r(m_geosEnvironment);
        GEOSGeometry* geom = GEOSWKTReader_read_r(m_geosEnvironment, reader, m_wkt.c_str());
        if (!geom)
            throw app_runtime_error("unable to ingest given WKT string");
#endif
            
    }
    return;
}


void PcQuery::addSwitches()
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



int PcQuery::execute()
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



    /*    
    flann::Matrix<float> X = new float[npts];
    flann::Matrix<float> Y = new float[npts];
    flann::Matrix<float> Z = new float[npts];
    */
#endif



    pdal::Options options = m_options + readerOptions;
    
    pdal::filters::Index* filter = new pdal::filters::Index(*stage, options);

    filter->initialize();

    boost::uint64_t numPoints = stage->getNumPoints();
    std::cout << numPoints << std::endl;
    const Schema& schema = stage->getSchema();
    PointBuffer* data = new PointBuffer(schema, 1);
    //PointBuffer data(schema, 0);
    boost::scoped_ptr<StageSequentialIterator>* iter = new boost::scoped_ptr<StageSequentialIterator>(stage->createSequentialIterator(*data));


    /*GET INFORMATION FROM WKT HERE*/
    int t_numPoints = 1;

     float* t_xyz = new float[3];
     t_xyz[0] = 0.0;
     t_xyz[1] = 0.0;
     t_xyz[2] = 0.0;




    /*END TEST WKT CODE*/
    
    int dim = 3;
    float* xyz = new float[numPoints*dim];
   
    //boost::property_tree::ptree data_container;

    boost::uint32_t itr = 0;
    while (!((**iter).atEnd()))
    {
        (**iter).read(*data);
        xyz[itr] = lexical_cast<float>(data -> getField<boost::uint32_t>(schema.getDimension("X"),0));
        xyz[itr+1] = lexical_cast<float>(data -> getField<boost::uint32_t>(schema.getDimension("Y"),0));
        xyz[itr+2] = lexical_cast<float>(data -> getField<boost::uint32_t>(schema.getDimension("Z"),0));
        itr ++;
        if (itr > numPoints)
        {
            std::cout << "Iterate Past End!" << std::endl;
        }
    }
    std::cout << "Data Read" << std::endl;
    int nn = 1;
    int* result = new int[dim*numPoints];
    float * dists = new float[numPoints*nn];
    float speedup;


    //Index<L2<float> > index(xyz, flann::KDTreeIndexParams(4));
    //index.buildIndex();
    //

    struct FLANNParameters p;
    flann_index_t index_id;

    p = DEFAULT_FLANN_PARAMETERS;
   

    //set Parameters
    //struct FLANNParameters p = DEFAULT_FLANN_PARAMETERS; 
    //p.algorithm=KDTREE;
    //p.target_precision = 0.9;

    //flann::flann_find_nearest_neighbors(xyz,numPoints, dim, t_xyz, t_numPoints, result, dists, nn, &p);
    
    delete[] xyz;
    delete[] result;
    delete[] dists;
    
    
    
    std::cout << std::endl;
    
    delete data;
    delete iter;
    delete filter;
    delete stage;






    
    return 0;
}


int main(int argc, char* argv[])
{
    PcQuery app(argc, argv);
    return app.run();
}

