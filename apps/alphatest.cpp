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
#include <pcl/ModelCoefficients.h>
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/sample_consensus/method_types.h>
#include <pcl/sample_consensus/model_types.h>
#include <pcl/filters/passthrough.h>
#include <pcl/filters/project_inliers.h>
#include <pcl/segmentation/sac_segmentation.h>
#include <pcl/surface/concave_hull.h>
#include <pcl/kdtree/kdtree_flann.h>

/*
#include <sensor_msgs/PointCloud2.h>


#include <pcl/visualization/common/common.h>
#include <pcl/ros/conversions.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <vtkTextActor.h>
#include <vtkTextProperty.h>
#include <vtkCellData.h>
#include <vtkWorldPointPicker.h>
#include <vtkPropPicker.h>
#include <vtkPlatonicSolidSource.h>
#include <vtkLoopSubdivisionFilter.h>
#include <vtkTriangle.h>
#include <vtkTransform.h>
#include <vtkVisibleCellSelector.h>
#include <vtkSelection.h>
#include <vtkPointPicker.h>
#include <vtkSmartPointer.h>
*/

namespace alphatest
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
        ("nsub,n", po::value<int>(&m_subset)->default_value(8), "How many subset operations?")


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
            //_x = (data -> getField<boost::int32_t>(dimx,i));
            //_y = (data -> getField<boost::int32_t>(dimy,i));

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
    //double density_alpha = (((xmax - xmin)*(ymax-ymin))/cells);
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
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_ptr (new pcl::PointCloud<pcl::PointXYZ>);

    (*cloud_ptr).points.resize(goodpoints -> size());
    

    //place centered points into cloud
    while (!(goodpoints -> empty()))
    {
        good_point = goodpoints -> top(); 

        (**iter).seek(good_point); 
        (**iter).read(*outdata);
        _x = (outdata -> getField<boost::int32_t>(dimx2,0));
        _y = (outdata -> getField<boost::int32_t>(dimy2,0));
        //_z = (outdata -> getField<boost::int32_t>(dimz2,0)); 
        (*cloud_ptr).points[itr].x = (_x - (xmin + (xmax-xmin)/2));
        (*cloud_ptr).points[itr].y = (_y - (ymin + (ymax-ymin)/2));
        //(*cloud_ptr).points[itr].z = _z;
        itr += 1;
        goodpoints -> pop();
    }
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_projected (new pcl::PointCloud<pcl::PointXYZ>);
    pcl::ModelCoefficients::Ptr coefficients(new pcl::ModelCoefficients);
    pcl::PointIndices::Ptr inliers (new pcl::PointIndices);
    pcl::SACSegmentation<pcl::PointXYZ> seg;
    seg.setOptimizeCoefficients(true);
    seg.setModelType(pcl::SACMODEL_PLANE);
    seg.setMethodType (pcl::SAC_RANSAC);
    seg.setDistanceThreshold(10);
    seg.setInputCloud(cloud_ptr -> makeShared());
    seg.segment(*inliers, *coefficients);
    std::cout << "After Segmentation: " << inliers -> indices.size() << std::endl;
    std::cout << "Model coefficients: " << coefficients->values[0] << " " 
                                        << coefficients->values[1] << " "
                                        << coefficients->values[2] << " " 
                                        << coefficients->values[3] << std::endl;
    pcl::ProjectInliers<pcl::PointXYZ> proj;
    proj.setModelType(pcl::SACMODEL_PLANE);
    proj.setInputCloud(cloud_ptr);
    proj.setModelCoefficients(coefficients);
    proj.filter(*cloud_projected);

    std::cout << "After Projection: " << cloud_projected  -> points.size() << std::endl;
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_hull(new pcl::PointCloud<pcl::PointXYZ>);

    pcl::ConcaveHull<pcl::PointXYZ> chull;
    chull.setInputCloud(cloud_projected);

    chull.setDimension(2);
    std::vector<pcl::Vertices>* polygons = new std::vector<pcl::Vertices>;


    int csize = 100000;
    double alpha = m_Alpha;

    // Try to avoid doing nearest neighbors search with too many points.
    while (csize > 10000)
    {
        std::cout << "Using Alpha = " << alpha << " Cloud size is: ";
        chull.setAlpha(alpha);
        chull.reconstruct(*cloud_hull, *polygons);
        std::cout << cloud_hull -> points.size() << std::endl;
        csize = cloud_hull -> points.size();
        chull.setInputCloud(cloud_hull);
        alpha *= 1.1;
    }

    /*
    while (csize == 0)
    {
        chull.setAlpha(2*mult);
        std::cout << "Trying Alpha: " << 2*mult << std::endl;
        chull.reconstruct(*cloud_hull, *polygons);
       
        csize = cloud_hull -> points.size();
        chull.setInputCloud(cloud_projected);
        mult += 10;
    }
    int calc_alpha = (2*mult)*10;
    std::cout << "Decided on Alpha = " << calc_alpha << std::endl;
    chull.setInputCloud(cloud_projected);
    for (double _a = 1; _a < 1.5; _a += 0.2)
    {
        std::cout << "Creating shape, Alpha = " << pow(calc_alpha,_a) << std::endl;
        chull.setAlpha(pow(calc_alpha,_a));
        chull.reconstruct(*cloud_hull, *polygons);
        chull.setInputCloud(cloud_hull);
        std::cout << "New Cloud Size: " << cloud_hull -> points.size() << std::endl;
    }
    */

    /*
    std::cout << "Alpha: " << m_Alpha << std::endl;
    chull.setAlpha(m_Alpha);
    chull.reconstruct(*cloud_hull, *polygons);
    if (cloud_hull -> points.size() == 0)
    {
        std::cout << "Warning, looks like pcl failed to make an alpha shape. Trying iterative approach." << std::endl;
        for (int _alpha = 10000; _alpha >= 1000; _alpha -= 1000)
        {
            std::cout << "Trying alpha = " << _alpha << std::endl;
            chull.setAlpha(_alpha);
            chull.reconstruct(*cloud_hull, *polygons);
            std::cout << "Concave hull has: " << cloud_hull -> points.size() << std::endl; 
            chull.setInputCloud(cloud_hull);
        }
        //chull.setAlpha(m_Alpha);
        //chull.reconstruct(*cloud_hull, *polygons);
    }
    */
    std::cout << "Concave hull has: " << cloud_hull -> points.size() << std::endl; 
    std::cout << "Polygon vec has size: " << polygons -> size() << std::endl;
    //for (int i = 0; i < polygons -> size(); i++)
    //{
    //    std::cout << "Vertex " << i << ", size: " << (*polygons)[i].vertices.size() << std::endl;
    //}
    std::cout << "Chull dimension: " << chull.getDimension() << std::endl;

    pcl::PointCloud<pcl::PointXYZ>::Ptr output_cloud_hull(new pcl::PointCloud<pcl::PointXYZ>);

    boost::int32_t lastx;
    boost::int32_t lasty;
    boost::int32_t lastz;

    _x = (*cloud_hull).points[0].x;
    _y = (*cloud_hull).points[0].y;
    //_z = (*cloud_hull).points[0].z;
    std::stack<boost::int32_t>* outpts = new std::stack<boost::int32_t>;
    outpts -> push(0);
    lastx = _x;
    lasty = _y;
    //lastz = _z;
    for (int i = 1; i < cloud_hull -> points.size(); i ++ )
    {
        _x = (*cloud_hull).points[i].x;
        _y = (*cloud_hull).points[i].y;
        //_z = (*cloud_hull).points[i].z;

        //if ((_x != lastx) || (_y != lasty) || (_z != lastz)){
        if ((_x != lastx) || (_y != lasty)){
            outpts -> push(i);
            lastx = _x;
            lasty = _y; 
            //lastz = _z;
        }
    }
    int outsize = outpts -> size();
    output_cloud_hull -> resize(outpts -> size());
    boost::uint64_t pt_idx;
    int i = 0;
    while (! (outpts -> empty()))
    {
        pt_idx = outpts -> top();
        _x = (*cloud_hull).points[pt_idx].x;
        _y = (*cloud_hull).points[pt_idx].y;
        //_z = (*cloud_hull).points[pt_idx].z;

        output_cloud_hull -> points[i].x = _x; 
        output_cloud_hull -> points[i].y = _y;         
        //output_cloud_hull -> points[i].z = _z; 
        outpts -> pop();
        i ++;
    } 


    std::cout << "Writing to: " << m_outputFile << std::endl;
    ofstream outfile;
    char* filepath = new char[m_outputFile.size() + 1];
    filepath[m_outputFile.size() + 1] = 0;
    memcpy(filepath, m_outputFile.c_str(), m_outputFile.size());
    outfile.open(filepath);
    outfile <<  fixed << setprecision(0);
    outfile << "MULTIPOLYGON((";


    for (int poly = 0; poly < polygons -> size(); poly ++ )
    { 
        outfile << "(";
        for (int vtx = 0; vtx < (*polygons)[poly].vertices.size(); vtx ++)
        {
            outfile << (output_cloud_hull -> points)[(*polygons)[poly].vertices[vtx]].x + xmin + (xmax-xmin)/2  << " " <<
                       (output_cloud_hull -> points)[(*polygons)[poly].vertices[vtx]].y + ymin + (ymax-ymin)/2; 

            if (vtx == ((*polygons)[poly].vertices.size() - 1))
            {
                //outfile <<", " <<(output_cloud_hull -> points)[(*polygons)[poly].vertices[0]].x + xmin + (xmax-xmin)/2  << " " <<
                //           (output_cloud_hull -> points)[(*polygons)[poly].vertices[0]].y + ymin + (ymax-ymin)/2; 
           
                outfile << ")";
            }
            else 
            {
                outfile << ",\n";
            }
        }
        if (poly == (polygons -> size() - 1))
        {
            outfile << "))";
        }
        else
        {
            outfile << ",\n";
        }
    }
    outfile.close();

    /*


    sensor_msgs::PointCloud2 msg_alpha;
    pcl::toROSMsg(*cloud_hull, msg_alpha);
    pcl::PolygonMesh mesh_alpha;
    mesh_alpha.cloud = msg_alpha;
    mesh_alpha.polygons = *polygons;
    

    //Code yanked from addPolylineFromPolygonMesh
    // Can probably get rid of toROSMsg calls, and keep everything in pcl
    vtkSmartPointer<vtkPoints> poly_points = vtkSmartPointer<vtkPoints>::New();
    pcl::PointCloud<pcl::PointXYZ> point_cloud;
    pcl::fromROSMsg (mesh_alpha.cloud, point_cloud);
    poly_points -> SetNumberOfPoints(point_cloud.points.size());
    size_t i;
    for (i = 0; i < point_cloud.points.size(); ++i)
    {
        poly_points -> InsertPoint(i, point_cloud.points[i].x, point_cloud.points[i].y, point_cloud.points[i].z);
    }
    vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkPolyData> polyData;
    for (i =0; i < mesh_alpha.polygons.size(); i++)
    {
        vtkSmartPointer<vtkPolyLine> polyLine = vtkSmartPointer<vtkPolyLine>::New();
        polyLine -> GetPointIds() -> SetNumberOfIds(mesh_alpha.polygons[i].vertices.size());
        for (unsigned int k = 0; k < mesh_alpha.polygons[i].vertices.size(); k++)
        {
            polyLine -> GetPointIds() -> SetId(k, mesh_alpha.polygons[i].vertices[k]);
        }
        cells -> InsertNextCell(polyLine);
    }
    std::cout << poly_points -> size() << std::endl;
    polyData -> SetPoints(poly_points);
    polyData -> SetLines(cells);
    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper -> SetInput(polyData);

    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor -> SetMapper(mapper);
    */
    //pcl::PCDWriter pclwriter;
    //pclwriter.write(m_outputFile, *output_cloud_hull, false);
    //pclwriter.write(m_outputFile, *cloud_hull, false);

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

