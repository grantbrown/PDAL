/******************************************************************************
* Copyright (c) 2011, Michael P. Gerlek (mpg@flaxen.com)
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

#ifndef INCLUDED_DRIVERS_LAS_SUPPORT_HPP
#define INCLUDED_DRIVERS_LAS_SUPPORT_HPP

#include <pdal/pdal_export.hpp>

#include <pdal/Schema.hpp>
#include <pdal/Stage.hpp>

#include <iostream>

#ifdef PDAL_HAVE_LASZIP
#include <laszip/laszip.hpp>
#include <laszip/lasunzipper.hpp>
#include <laszip/laszipper.hpp>

#else
typedef struct LASzipH *LASzip;
typedef struct LASunzipperH *LASunzipper;
#endif

namespace pdal
{
namespace drivers
{
namespace las
{

class SummaryData;
class ZipPoint;

enum PointFormat
{
    PointFormat0 = 0,         // base
    PointFormat1 = 1,         // base + time
    PointFormat2 = 2,         // base + color
    PointFormat3 = 3,         // base + time + color
    PointFormat4 = 4,         // base + time + wave
    PointFormat5 = 5,         // base + time + color + wave  (NOT SUPPORTED)
    PointFormatUnknown = 99
};

class PDAL_DLL PointDimensions
{
public:
    PointDimensions(const Schema& schema, std::string const& ns);

    pdal::Dimension const* X;
    pdal::Dimension const* Y;
    pdal::Dimension const* Z;

    pdal::Dimension const* Intensity;
    pdal::Dimension const* ReturnNumber;
    pdal::Dimension const* NumberOfReturns;
    pdal::Dimension const* ScanDirectionFlag;
    pdal::Dimension const* EdgeOfFlightLine;
    pdal::Dimension const* Classification;
    pdal::Dimension const* ScanAngleRank;
    pdal::Dimension const* UserData;
    pdal::Dimension const* PointSourceId;

    pdal::Dimension const* Time;

    pdal::Dimension const* Red;
    pdal::Dimension const* Green;
    pdal::Dimension const* Blue;
};

class PDAL_DLL Support
{
public:
    static void registerFields(Stage& stage, Schema& schema, PointFormat pointFormat);
    static void setScaling(Schema& schema, double scaleX, double scaleY, double scaleZ, double offsetX, double offsetY, double offsetZ);

    static bool hasTime(PointFormat);
    static bool hasColor(PointFormat);
    static bool hasWave(PointFormat);
    static boost::uint16_t getPointDataSize(PointFormat pointFormat);

    // assumes the stream position is pointing to the first byte of the header
    // this function updates the header's min/max xyz fields, and the point return counts fields
    static void rewriteHeader(std::ostream& stream, const SummaryData& data);
};


}
}
} // namespace

#endif
