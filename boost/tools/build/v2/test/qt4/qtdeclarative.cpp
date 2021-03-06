// (c) Copyright Juergen Hunold 2011
// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#define BOOST_TEST_MODULE QtMultimedia

#include <QApplication>
#include <QDeclarativeView>

#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE( defines)
{
    BOOST_CHECK_EQUAL(BOOST_IS_DEFINED(QT_CORE_LIB), true);
    BOOST_CHECK_EQUAL(BOOST_IS_DEFINED(QT_GUI_LIB), true);
    BOOST_CHECK_EQUAL(BOOST_IS_DEFINED(QT_XML_LIB), true);
    BOOST_CHECK_EQUAL(BOOST_IS_DEFINED(QT_DECLARATIVE_LIB), true);
}


BOOST_AUTO_TEST_CASE( declarative )
{
    QApplication app(pdalboost::unit_test::framework::master_test_suite().argc,
                     pdalboost::unit_test::framework::master_test_suite().argv);
    QDeclarativeView view;
}
