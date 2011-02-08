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

#include <cassert>

#include "libpc/FauxReader.hpp"
#include "libpc/Utils.hpp"

using std::vector;
using std::string;


FauxReader::FauxReader(const Bounds& bounds, int numPoints)
  : Reader()
{
  Header& header = getHeader();
  PointLayout& layout = header.getLayout();
  
  header.setNumPoints(numPoints);
  header.setBounds(bounds);

  int index;
  
  index = layout.addField(Field(Field::XPos, Field::F32, true));
  index = layout.addField(Field(Field::YPos, Field::F32, true));
  index = layout.addField(Field(Field::ZPos, Field::F32, true));
  index = layout.addField(Field(Field::Time, Field::F64, true));

  header.dump();

  return;
}


void FauxReader::readNextPoints(PointData& data)
{
  // make up some data and put it into the buffer

  int numPoints = data.getNumPoints();
  assert(m_lastPointRead + numPoints <= getHeader().getNumPoints());

  const PointLayout& layout = data.getLayout();
  Header& header = getHeader();

  int fieldIndexT = layout.findFieldIndex(Field::Time);
  assert(fieldIndexT != -1);

  float v = (float)m_lastPointRead;

  const Bounds& bounds = header.getBounds();

  for (int pointIndex=0; pointIndex<numPoints; pointIndex++)
  {
    const float x = (float)Utils::random(bounds.m_minX,bounds.m_maxX);
    const float y = (float)Utils::random(bounds.m_minY,bounds.m_maxY);
    const float z = (float)Utils::random(bounds.m_minZ,bounds.m_maxZ);

    data.setValid(pointIndex);

    data.setX(pointIndex, x);
    data.setY(pointIndex, y);
    data.setZ(pointIndex, z);
    
    data.setField_F64(pointIndex, fieldIndexT, v * 0.1);

    ++v;
  }

  m_lastPointRead += numPoints;

  return;
}