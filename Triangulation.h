/*
 Copyright (c) 2015 Simon Geilfus

 This code uses a class extracted and adapted from the OpenCV library
 for use with the Cinder C++ library, http://libcinder.org and all the 
 relevant portion of code is tied to OpenCV license that can be found 
 in Triangulation.cpp file
 
 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:
 
 The above copyright notice and this permission notice shall be included in all
 copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 SOFTWARE.
 */

#pragma once

#include "cinder/Rect.h"
#include "cinder/Vector.h"

namespace Delaunay {
	
	//! triangulates the list of points and returns a list of triangles indices
	std::vector<uint32_t> getTriangleIndices( const ci::Rectf &rect, const std::vector<ci::vec2> &points = std::vector<ci::vec2>() );
};