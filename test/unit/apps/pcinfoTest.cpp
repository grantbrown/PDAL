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

#include <boost/test/unit_test.hpp>
#include <pdal/FileUtils.hpp>
#include "Support.hpp"

#include <iostream>
#include <sstream>
#include <string>


BOOST_AUTO_TEST_SUITE(pcinfoTest)


static std::string appName()
{
    const std::string app = Support::binpath(Support::exename("pcinfo"));
    BOOST_CHECK(pdal::FileUtils::fileExists(app));
    return app;
}


#ifdef PDAL_COMPILER_MSVC
BOOST_AUTO_TEST_CASE(pcinfoTest_no_input)
{
    const std::string cmd = appName();

    std::string output;
    int stat = Support::run_command(cmd, output);
    BOOST_CHECK_EQUAL(stat, 1);

    const std::string expected = "Usage error: no action option specified";
    BOOST_CHECK_EQUAL(output.substr(0, expected.length()), expected);

    return;
}
#endif


BOOST_AUTO_TEST_CASE(pcinfo_test_common_opts)
{
    const std::string cmd = appName();

    std::string output;
    int stat = Support::run_command(cmd + " -h", output);
    BOOST_CHECK_EQUAL(stat, 0);

    stat = Support::run_command(cmd + " --version", output);
    BOOST_CHECK_EQUAL(stat, 0);

    return;
}


BOOST_AUTO_TEST_CASE(pcinfo_test_switches)
{
#if 1
    const std::string cmd = appName();

    std::string inputLas = Support::datapath("apps/simple.las");
    std::string inputLaz = Support::datapath("apps/simple.laz");

    std::string output;
    std::string expected;

    int stat = 0;

    // does the default work?
    stat = Support::run_command(cmd + " " + inputLas, output);
    BOOST_CHECK_EQUAL(stat, 1);
    expected = "Usage error: no action option specified";
    BOOST_CHECK_EQUAL(output.substr(0, expected.length()), expected);

    // does --input work?
    stat = Support::run_command(cmd + " --input=" + inputLas, output);
    BOOST_CHECK_EQUAL(stat, 1);
    expected = "Usage error: no action option specified";
    BOOST_CHECK_EQUAL(output.substr(0, expected.length()), expected);

    // does -i work?
    stat = Support::run_command(cmd + " -i " + inputLas, output);
    BOOST_CHECK_EQUAL(stat, 1);
    expected = "Usage error: no action option specified";
    BOOST_CHECK_EQUAL(output.substr(0, expected.length()), expected);

#ifdef PDAL_HAVE_LASZIP
    // does it work for .laz?
    stat = Support::run_command(cmd + " " + inputLaz, output);
    BOOST_CHECK_EQUAL(stat, 1);
    expected = "Usage error: no action option specified";
    BOOST_CHECK_EQUAL(output.substr(0, expected.length()), expected);
#endif

#ifdef PDAL_HAVE_LIBLAS
    // does --liblas work?
    stat = Support::run_command(cmd + " --liblas " + inputLas, output);
    BOOST_CHECK_EQUAL(stat, 1);
    expected = "Usage error: no action option specified";
    BOOST_CHECK_EQUAL(output.substr(0, expected.length()), expected);
#endif

#ifdef PDAL_HAVE_LIBLAS
#ifdef PDAL_HAVE_LASZIP
    // does --liblas work for .laz too?
    stat = Support::run_command(cmd + " --liblas " + inputLaz, output);
    BOOST_CHECK_EQUAL(stat, 1);
    expected = "Usage error: no action option specified";
    BOOST_CHECK_EQUAL(output.substr(0, expected.length()), expected);
#endif
#endif
#endif
    return;
}


BOOST_AUTO_TEST_CASE(pcinfo_test_dumps)
{
    const std::string cmd = appName();

    const std::string inputLas = Support::datapath("apps/simple.las");
    const std::string inputLaz = Support::datapath("apps/simple.laz");
    bool were_equal(false);
    std::string output;

    int stat = 0;

    std::ostringstream command;

    // dump a single point to json

    std::string pt_test = Support::temppath("pcinfo_point.txt");
    command << cmd + " --point=1 " + inputLas + " > " + pt_test;
    stat = Support::run_command(command.str(), output);
    BOOST_CHECK_EQUAL(stat, 0);
    were_equal = Support::compare_text_files(pt_test, Support::datapath("apps/pcinfo_point.txt"));
    BOOST_CHECK(were_equal);
    if (were_equal)
        pdal::FileUtils::deleteFile(pt_test);
    else
        std::cout << command.str() << std::endl;

    // dump summary of all points to json
    command.str("");

    std::string stats_test = Support::temppath("pcinfo_stats.txt");
    command << cmd + " --stats " + inputLas + " --seed 1234" +" > " + stats_test; 
    stat = Support::run_command(command.str(), output);
    BOOST_CHECK_EQUAL(stat, 0);
#if defined(PDAL_PLATFORM_WIN32)
    were_equal = Support::compare_text_files(stats_test, Support::datapath("apps/pcinfo_stats-win32.txt"));
#else
    were_equal = Support::compare_text_files(stats_test, Support::datapath("apps/pcinfo_stats.txt"));
#endif
    BOOST_CHECK(were_equal);
    if (were_equal)
        pdal::FileUtils::deleteFile(stats_test);
    else
        std::cout << command.str() << std::endl;

    // dump schema to json
    command.str("");

    std::string schema_test = Support::temppath("pcinfo_schema.txt");
    command << cmd + " --schema " + inputLas +" > " + schema_test;
    stat = Support::run_command(command.str(), output);
    BOOST_CHECK_EQUAL(stat, 0);
    were_equal = Support::compare_text_files(schema_test, Support::datapath("apps/pcinfo_schema.txt"));
    BOOST_CHECK(were_equal);
    if (were_equal)
        pdal::FileUtils::deleteFile(schema_test);
    else
        std::cout << command.str() << std::endl;

    // dump stage info to json
    command.str("");

    std::string stage_test = Support::temppath("pcinfo_stage.txt");
    command << cmd + " --stage " + inputLas +" > " + stage_test;
    stat = Support::run_command(command.str(), output);
    BOOST_CHECK_EQUAL(stat, 0);

#ifdef PDAL_HAVE_GDAL
    unsigned int check = Support::diff_text_files(stage_test, Support::datapath("apps/pcinfo_stage.txt"), 15);
    BOOST_CHECK_EQUAL(check, 0u);
#else
    unsigned int check = Support::diff_text_files(stage_test, Support::datapath("apps/pcinfo_stage_nosrs.txt"), 15);
    BOOST_CHECK_EQUAL(check, 0u);
#endif
    if (check == 0u)
        pdal::FileUtils::deleteFile(stage_test);
    else
        std::cout << command.str() << std::endl;

    return;
}


BOOST_AUTO_TEST_SUITE_END()
