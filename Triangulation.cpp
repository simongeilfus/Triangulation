/*
 IMPORTANT: READ BEFORE DOWNLOADING, COPYING, INSTALLING OR USING.
 
 By downloading, copying, installing or using the software you agree to this license.
 If you do not agree to this license, do not download, install,
 copy or use the software.
 
 
 Intel License Agreement
 For Open Source Computer Vision Library
 
 Copyright (C) 2000, Intel Corporation, all rights reserved.
 Third party copyrights are property of their respective owners.
 
 Redistribution and use in source and binary forms, with or without modification,
 are permitted provided that the following conditions are met:
 
 * Redistribution's of source code must retain the above copyright notice,
 this list of conditions and the following disclaimer.
 
 * Redistribution's in binary form must reproduce the above copyright notice,
 this list of conditions and the following disclaimer in the documentation
 and/or other materials provided with the distribution.
 
 * The name of Intel Corporation may not be used to endorse or promote products
 derived from this software without specific prior written permission.
 
 This software is provided by the copyright holders and contributors "as is" and
 any express or implied warranties, including, but not limited to, the implied
 warranties of merchantability and fitness for a particular purpose are disclaimed.
 In no event shall the Intel Corporation or contributors be liable for any direct,
 indirect, incidental, special, exemplary, or consequential damages
 (including, but not limited to, procurement of substitute goods or services;
 loss of use, data, or profits; or business interruption) however caused
 and on any theory of liability, whether in contract, strict liability,
 or tort (including negligence or otherwise) arising in any way out of
 the use of this software, even if advised of the possibility of such damage.
 */

#include "Triangulation.h"

#include "cinder/CinderAssert.h"
#include "cinder/Log.h"

using namespace ci;
using namespace std;

namespace {
	
	class Subdiv2D {
	public:
		
		Subdiv2D();
		Subdiv2D( ci::Rectf rect, const std::vector<ci::vec2> &points = std::vector<ci::vec2>() );
		void initDelaunay( ci::Rectf rect, const std::vector<ci::vec2> &points = std::vector<ci::vec2>() );
		
		enum {
			PTLOC_ERROR        = -2,
			PTLOC_OUTSIDE_RECT = -1,
			PTLOC_INSIDE       = 0,
			PTLOC_VERTEX       = 1,
			PTLOC_ON_EDGE      = 2
		};
		
		enum {
			NEXT_AROUND_ORG   = 0x00,
			NEXT_AROUND_DST   = 0x22,
			PREV_AROUND_ORG   = 0x11,
			PREV_AROUND_DST   = 0x33,
			NEXT_AROUND_LEFT  = 0x13,
			NEXT_AROUND_RIGHT = 0x31,
			PREV_AROUND_LEFT  = 0x20,
			PREV_AROUND_RIGHT = 0x02
		};
		
		struct Triangle {
			Triangle( const ci::vec2 &a, const ci::vec2 &b, const ci::vec2 &c ) : mA( a ), mB( b ), mC( c ) {}
			ci::vec2 mA, mB, mC;
		};
		
		int insert(ci::vec2 pt, int indice = -1 );
		void insert(const std::vector<ci::vec2>& ptvec);
		int locate(ci::vec2 pt, int& edge, int& vertex);
		
		int findNearest(ci::vec2 pt, ci::vec2* nearestPt = 0);
		void getEdgeList(std::vector<ci::vec4>& edgeList) const;
		void getTriangleList(std::vector<Triangle>& triangleList) const;
		std::vector<uint32_t> getTriangleIndices();
		void getVoronoiFacetList(const std::vector<int>& idx, std::vector<std::vector<ci::vec2> >& facetList,
								 std::vector<ci::vec2>& facetCenters);
		
		ci::vec2 getVertex(int vertex, int* firstEdge = 0) const;
		
		int getEdge( int edge, int nextEdgeType ) const;
		int nextEdge(int edge) const;
		int rotateEdge(int edge, int rotate) const;
		int symEdge(int edge) const;
		int edgeOrg(int edge, ci::vec2* orgpt = 0) const;
		int edgeDst(int edge, ci::vec2* dstpt = 0) const;
		
	protected:
		int newEdge();
		void deleteEdge(int edge);
		int newPoint(ci::vec2 pt, bool isvirtual, int firstEdge = 0, int indice = -1 );
		void deletePoint(int vtx);
		void setEdgePoints( int edge, int orgPt, int dstPt );
		void splice( int edgeA, int edgeB );
		int connectEdges( int edgeA, int edgeB );
		void swapEdges( int edge );
		int isRightOf(ci::vec2 pt, int edge) const;
		void calcVoronoi();
		void clearVoronoi();
		void checkSubdiv() const;
		
		struct Vertex
		{
			Vertex();
			Vertex(ci::vec2 pt, bool _isvirtual, int _firstEdge=0, int indice = -1 );
			bool isvirtual() const;
			bool isfree() const;
			
			int firstEdge;
			int type;
			ci::vec2 pt;
			int originalIndice;
		};
		
		struct QuadEdge
		{
			QuadEdge();
			QuadEdge(int edgeidx);
			bool isfree() const;
			
			int next[4];
			int pt[4];
		};
		
		std::vector<Vertex> vtx;
		std::vector<QuadEdge> qedges;
		int freeQEdge;
		int freePoint;
		bool validGeometry;
		
		int recentEdge;
		ci::vec2 topLeft;
		ci::vec2 bottomRight;
	};
	
	int Subdiv2D::nextEdge(int edge) const
	{
		CI_ASSERT((size_t)(edge >> 2) < qedges.size());
		return qedges[edge >> 2].next[edge & 3];
	}
	
	int Subdiv2D::rotateEdge(int edge, int rotate) const
	{
		return (edge & ~3) + ((edge + rotate) & 3);
	}
	
	int Subdiv2D::symEdge(int edge) const
	{
		return edge ^ 2;
	}
	
	int Subdiv2D::getEdge(int edge, int nextEdgeType) const
	{
		CI_ASSERT((size_t)(edge >> 2) < qedges.size());
		edge = qedges[edge >> 2].next[(edge + nextEdgeType) & 3];
		return (edge & ~3) + ((edge + (nextEdgeType >> 4)) & 3);
	}
	
	int Subdiv2D::edgeOrg(int edge, vec2* orgpt) const
	{
		CI_ASSERT((size_t)(edge >> 2) < qedges.size());
		int vidx = qedges[edge >> 2].pt[edge & 3];
		if( orgpt )
		{
			CI_ASSERT((size_t)vidx < vtx.size());
			*orgpt = vtx[vidx].pt;
		}
		return vidx;
	}
	
	int Subdiv2D::edgeDst(int edge, vec2* dstpt) const
	{
		CI_ASSERT((size_t)(edge >> 2) < qedges.size());
		int vidx = qedges[edge >> 2].pt[(edge + 2) & 3];
		if( dstpt )
		{
			CI_ASSERT((size_t)vidx < vtx.size());
			*dstpt = vtx[vidx].pt;
		}
		return vidx;
	}
	
	
	vec2 Subdiv2D::getVertex(int vertex, int* firstEdge) const
	{
		CI_ASSERT((size_t)vertex < vtx.size());
		if( firstEdge )
			*firstEdge = vtx[vertex].firstEdge;
		return vtx[vertex].pt;
	}
	
	
	Subdiv2D::Subdiv2D()
	{
		validGeometry = false;
		freeQEdge = 0;
		freePoint = 0;
		recentEdge = 0;
	}
	
	Subdiv2D::Subdiv2D( Rectf rect, const std::vector<ci::vec2> &points )
	{
		validGeometry = false;
		freeQEdge = 0;
		freePoint = 0;
		recentEdge = 0;
		
		initDelaunay( rect, points );
	}
	
	
	Subdiv2D::QuadEdge::QuadEdge()
	{
		next[0] = next[1] = next[2] = next[3] = 0;
		pt[0] = pt[1] = pt[2] = pt[3] = 0;
	}
	
	Subdiv2D::QuadEdge::QuadEdge(int edgeidx)
	{
		CI_ASSERT((edgeidx & 3) == 0);
		next[0] = edgeidx;
		next[1] = edgeidx+3;
		next[2] = edgeidx+2;
		next[3] = edgeidx+1;
		
		pt[0] = pt[1] = pt[2] = pt[3] = 0;
	}
	
	bool Subdiv2D::QuadEdge::isfree() const
	{
		return next[0] <= 0;
	}
	
	Subdiv2D::Vertex::Vertex()
	{
		firstEdge = 0;
		type = -1;
	}
	
	Subdiv2D::Vertex::Vertex(vec2 _pt, bool _isvirtual, int _firstEdge, int indice )
	{
		firstEdge = _firstEdge;
		type = (int)_isvirtual;
		pt = _pt;
		originalIndice = indice;
	}
	
	bool Subdiv2D::Vertex::isvirtual() const
	{
		return type > 0;
	}
	
	bool Subdiv2D::Vertex::isfree() const
	{
		return type < 0;
	}
	
	void Subdiv2D::splice( int edgeA, int edgeB )
	{
		int& a_next = qedges[edgeA >> 2].next[edgeA & 3];
		int& b_next = qedges[edgeB >> 2].next[edgeB & 3];
		int a_rot = rotateEdge(a_next, 1);
		int b_rot = rotateEdge(b_next, 1);
		int& a_rot_next = qedges[a_rot >> 2].next[a_rot & 3];
		int& b_rot_next = qedges[b_rot >> 2].next[b_rot & 3];
		std::swap(a_next, b_next);
		std::swap(a_rot_next, b_rot_next);
	}
	
	void Subdiv2D::setEdgePoints(int edge, int orgPt, int dstPt)
	{
		qedges[edge >> 2].pt[edge & 3] = orgPt;
		qedges[edge >> 2].pt[(edge + 2) & 3] = dstPt;
		vtx[orgPt].firstEdge = edge;
		vtx[dstPt].firstEdge = edge ^ 2;
	}
	
	int Subdiv2D::connectEdges( int edgeA, int edgeB )
	{
		int edge = newEdge();
		
		splice(edge, getEdge(edgeA, NEXT_AROUND_LEFT));
		splice(symEdge(edge), edgeB);
		
		setEdgePoints(edge, edgeDst(edgeA), edgeOrg(edgeB));
		return edge;
	}
	
	void Subdiv2D::swapEdges( int edge )
	{
		int sedge = symEdge(edge);
		int a = getEdge(edge, PREV_AROUND_ORG);
		int b = getEdge(sedge, PREV_AROUND_ORG);
		
		splice(edge, a);
		splice(sedge, b);
		
		setEdgePoints(edge, edgeDst(a), edgeDst(b));
		
		splice(edge, getEdge(a, NEXT_AROUND_LEFT));
		splice(sedge, getEdge(b, NEXT_AROUND_LEFT));
	}
	
	static double triangleArea( vec2 a, vec2 b, vec2 c )
	{
		return ((double)b.x - a.x) * ((double)c.y - a.y) - ((double)b.y - a.y) * ((double)c.x - a.x);
	}
	
	int Subdiv2D::isRightOf(vec2 pt, int edge) const
	{
		vec2 org, dst;
		edgeOrg(edge, &org);
		edgeDst(edge, &dst);
		double cw_area = triangleArea( pt, dst, org );
		
		return (cw_area > 0) - (cw_area < 0);
	}
	
	int Subdiv2D::newEdge()
	{
		if( freeQEdge <= 0 )
		{
			qedges.push_back(QuadEdge());
			freeQEdge = (int)(qedges.size()-1);
		}
		int edge = freeQEdge*4;
		freeQEdge = qedges[edge >> 2].next[1];
		qedges[edge >> 2] = QuadEdge(edge);
		return edge;
	}
	
	void Subdiv2D::deleteEdge(int edge)
	{
		CI_ASSERT((size_t)(edge >> 2) < (size_t)qedges.size());
		splice( edge, getEdge(edge, PREV_AROUND_ORG) );
		int sedge = symEdge(edge);
		splice(sedge, getEdge(sedge, PREV_AROUND_ORG) );
		
		edge >>= 2;
		qedges[edge].next[0] = 0;
		qedges[edge].next[1] = freeQEdge;
		freeQEdge = edge;
	}
	
	int Subdiv2D::newPoint(vec2 pt, bool isvirtual, int firstEdge, int indice )
	{
		if( freePoint == 0 )
		{
			vtx.push_back(Vertex());
			freePoint = (int)(vtx.size()-1);
		}
		int vidx = freePoint;
		freePoint = vtx[vidx].firstEdge;
		vtx[vidx] = Vertex(pt, isvirtual, firstEdge, indice);
		
		return vidx;
	}
	
	void Subdiv2D::deletePoint(int vidx)
	{
		CI_ASSERT( (size_t)vidx < vtx.size() );
		vtx[vidx].firstEdge = freePoint;
		vtx[vidx].type = -1;
		freePoint = vidx;
	}
	
	int Subdiv2D::locate(vec2 pt, int& _edge, int& _vertex)
	{
		int vertex = 0;
		
		int i, maxEdges = (int)(qedges.size() * 4);
		
		if( qedges.size() < (size_t)4 )
			CI_LOG_E( "Subdiv2D Error: Subdivision is empty" );
		
		if( pt.x < topLeft.x || pt.y < topLeft.y || pt.x >= bottomRight.x || pt.y >= bottomRight.y )
			CI_LOG_E( "Subdiv2D OutOfRange" );
		
		int edge = recentEdge;
		CI_ASSERT(edge > 0);
		
		int location = PTLOC_ERROR;
		
		int right_of_curr = isRightOf(pt, edge);
		if( right_of_curr > 0 )
		{
			edge = symEdge(edge);
			right_of_curr = -right_of_curr;
		}
		
		for( i = 0; i < maxEdges; i++ )
		{
			int onext_edge = nextEdge( edge );
			int dprev_edge = getEdge( edge, PREV_AROUND_DST );
			
			int right_of_onext = isRightOf( pt, onext_edge );
			int right_of_dprev = isRightOf( pt, dprev_edge );
			
			if( right_of_dprev > 0 )
			{
				if( right_of_onext > 0 || (right_of_onext == 0 && right_of_curr == 0) )
				{
					location = PTLOC_INSIDE;
					break;
				}
				else
				{
					right_of_curr = right_of_onext;
					edge = onext_edge;
				}
			}
			else
			{
				if( right_of_onext > 0 )
				{
					if( right_of_dprev == 0 && right_of_curr == 0 )
					{
						location = PTLOC_INSIDE;
						break;
					}
					else
					{
						right_of_curr = right_of_dprev;
						edge = dprev_edge;
					}
				}
				else if( right_of_curr == 0 &&
						isRightOf( vtx[edgeDst(onext_edge)].pt, edge ) >= 0 )
				{
					edge = symEdge( edge );
				}
				else
				{
					right_of_curr = right_of_onext;
					edge = onext_edge;
				}
			}
		}
		
		recentEdge = edge;
		
		if( location == PTLOC_INSIDE )
		{
			vec2 org_pt, dst_pt;
			edgeOrg(edge, &org_pt);
			edgeDst(edge, &dst_pt);
			
			double t1 = fabs( pt.x - org_pt.x );
			t1 += fabs( pt.y - org_pt.y );
			double t2 = fabs( pt.x - dst_pt.x );
			t2 += fabs( pt.y - dst_pt.y );
			double t3 = fabs( org_pt.x - dst_pt.x );
			t3 += fabs( org_pt.y - dst_pt.y );
			
			if( t1 < FLT_EPSILON )
			{
				location = PTLOC_VERTEX;
				vertex = edgeOrg( edge );
				edge = 0;
			}
			else if( t2 < FLT_EPSILON )
			{
				location = PTLOC_VERTEX;
				vertex = edgeDst( edge );
				edge = 0;
			}
			else if( (t1 < t3 || t2 < t3) &&
					fabs( triangleArea( pt, org_pt, dst_pt )) < FLT_EPSILON )
			{
				location = PTLOC_ON_EDGE;
				vertex = 0;
			}
		}
		
		if( location == PTLOC_ERROR )
		{
			edge = 0;
			vertex = 0;
		}
		
		_edge = edge;
		_vertex = vertex;
		
		return location;
	}
	
	
	inline int
	isPtInCircle3( vec2 pt, vec2 a, vec2 b, vec2 c)
	{
		const double eps = FLT_EPSILON*0.125;
		double val = ((double)a.x * a.x + (double)a.y * a.y) * triangleArea( b, c, pt );
		val -= ((double)b.x * b.x + (double)b.y * b.y) * triangleArea( a, c, pt );
		val += ((double)c.x * c.x + (double)c.y * c.y) * triangleArea( a, b, pt );
		val -= ((double)pt.x * pt.x + (double)pt.y * pt.y) * triangleArea( a, b, c );
		
		return val > eps ? 1 : val < -eps ? -1 : 0;
	}
	
	
	int Subdiv2D::insert(vec2 pt, int indice)
	{
		int curr_point = 0, curr_edge = 0, deleted_edge = 0;
		int location = locate( pt, curr_edge, curr_point );
		
		if( location == PTLOC_ERROR )
			CI_LOG_E( "Subdiv2D BadSize" );
		
		if( location == PTLOC_OUTSIDE_RECT )
			CI_LOG_E( "Subdiv2D OutOfRange" );
		
		if( location == PTLOC_VERTEX )
			return curr_point;
		
		if( location == PTLOC_ON_EDGE )
		{
			deleted_edge = curr_edge;
			recentEdge = curr_edge = getEdge( curr_edge, PREV_AROUND_ORG );
			deleteEdge(deleted_edge);
		}
		else if( location == PTLOC_INSIDE )
			;
		else
			CI_LOG_E( "Subdiv2D::locate returned invalid location " << location );
		
		assert( curr_edge != 0 );
		validGeometry = false;
		
		curr_point = newPoint(pt, false,0,indice);
		int base_edge = newEdge();
		int first_point = edgeOrg(curr_edge);
		setEdgePoints(base_edge, first_point, curr_point);
		splice(base_edge, curr_edge);
		
		do
		{
			base_edge = connectEdges( curr_edge, symEdge(base_edge) );
			curr_edge = getEdge(base_edge, PREV_AROUND_ORG);
		}
		while( edgeDst(curr_edge) != first_point );
		
		curr_edge = getEdge( base_edge, PREV_AROUND_ORG );
		
		int i, max_edges = (int)(qedges.size()*4);
		
		for( i = 0; i < max_edges; i++ )
		{
			int temp_dst = 0, curr_org = 0, curr_dst = 0;
			int temp_edge = getEdge( curr_edge, PREV_AROUND_ORG );
			
			temp_dst = edgeDst( temp_edge );
			curr_org = edgeOrg( curr_edge );
			curr_dst = edgeDst( curr_edge );
			
			if( isRightOf( vtx[temp_dst].pt, curr_edge ) > 0 &&
			   isPtInCircle3( vtx[curr_org].pt, vtx[temp_dst].pt,
							 vtx[curr_dst].pt, vtx[curr_point].pt ) < 0 )
			{
				swapEdges( curr_edge );
				curr_edge = getEdge( curr_edge, PREV_AROUND_ORG );
			}
			else if( curr_org == first_point )
				break;
			else
				curr_edge = getEdge( nextEdge( curr_edge ), PREV_AROUND_LEFT );
		}
		
		return curr_point;
	}
	
	void Subdiv2D::insert(const std::vector<vec2>& ptvec)
	{
		for( size_t i = 0; i < ptvec.size(); i++ )
			insert(ptvec[i],i);
	}
	
	void Subdiv2D::initDelaunay( Rectf rect, const std::vector<ci::vec2> &points )
	{
		float big_coord = 3.f * max( rect.getWidth(), rect.getHeight() );
		float rx = (float) rect.getX1();
		float ry = (float) rect.getY1();
		
		vtx.clear();
		qedges.clear();
		
		recentEdge = 0;
		validGeometry = false;
		
		topLeft = vec2( rx, ry );
		bottomRight = vec2( rx + rect.getWidth(), ry + rect.getHeight() );
		
		vec2 ppA( rx + big_coord, ry );
		vec2 ppB( rx, ry + big_coord );
		vec2 ppC( rx - big_coord, ry - big_coord );
		
		vtx.push_back(Vertex());
		qedges.push_back(QuadEdge());
		
		freeQEdge = 0;
		freePoint = 0;
		
		int pA = newPoint(ppA, false);
		int pB = newPoint(ppB, false);
		int pC = newPoint(ppC, false);
		
		int edge_AB = newEdge();
		int edge_BC = newEdge();
		int edge_CA = newEdge();
		
		setEdgePoints( edge_AB, pA, pB );
		setEdgePoints( edge_BC, pB, pC );
		setEdgePoints( edge_CA, pC, pA );
		
		splice( edge_AB, symEdge( edge_CA ));
		splice( edge_BC, symEdge( edge_AB ));
		splice( edge_CA, symEdge( edge_BC ));
		
		recentEdge = edge_AB;
		
		insert( points );
	}
	
	
	void Subdiv2D::clearVoronoi()
	{
		size_t i, total = qedges.size();
		
		for( i = 0; i < total; i++ )
			qedges[i].pt[1] = qedges[i].pt[3] = 0;
		
		total = vtx.size();
		for( i = 0; i < total; i++ )
		{
			if( vtx[i].isvirtual() )
				deletePoint((int)i);
		}
		
		validGeometry = false;
	}
	
	
	static vec2 computeVoronoiPoint(vec2 org0, vec2 dst0, vec2 org1, vec2 dst1)
	{
		double a0 = dst0.x - org0.x;
		double b0 = dst0.y - org0.y;
		double c0 = -0.5*(a0 * (dst0.x + org0.x) + b0 * (dst0.y + org0.y));
		
		double a1 = dst1.x - org1.x;
		double b1 = dst1.y - org1.y;
		double c1 = -0.5*(a1 * (dst1.x + org1.x) + b1 * (dst1.y + org1.y));
		
		double det = a0 * b1 - a1 * b0;
		
		if( det != 0 )
		{
			det = 1. / det;
			return vec2((float) ((b0 * c1 - b1 * c0) * det),
						(float) ((a1 * c0 - a0 * c1) * det));
		}
		
		return vec2(FLT_MAX, FLT_MAX);
	}
	
	
	void Subdiv2D::calcVoronoi()
	{
		// check if it is already calculated
		if( validGeometry )
			return;
		
		clearVoronoi();
		int i, total = (int)qedges.size();
		
		// loop through all quad-edges, except for the first 3 (#1, #2, #3 - 0 is reserved for "NULL" pointer)
		for( i = 4; i < total; i++ )
		{
			QuadEdge& quadedge = qedges[i];
			
			if( quadedge.isfree() )
				continue;
			
			int edge0 = (int)(i*4);
			vec2 org0, dst0, org1, dst1;
			
			if( !quadedge.pt[3] )
			{
				int edge1 = getEdge( edge0, NEXT_AROUND_LEFT );
				int edge2 = getEdge( edge1, NEXT_AROUND_LEFT );
				
				edgeOrg(edge0, &org0);
				edgeDst(edge0, &dst0);
				edgeOrg(edge1, &org1);
				edgeDst(edge1, &dst1);
				
				vec2 virt_point = computeVoronoiPoint(org0, dst0, org1, dst1);
				
				if( fabs( virt_point.x ) < FLT_MAX * 0.5 &&
				   fabs( virt_point.y ) < FLT_MAX * 0.5 )
				{
					quadedge.pt[3] = qedges[edge1 >> 2].pt[3 - (edge1 & 2)] =
					qedges[edge2 >> 2].pt[3 - (edge2 & 2)] = newPoint(virt_point, true);
				}
			}
			
			if( !quadedge.pt[1] )
			{
				int edge1 = getEdge( edge0, NEXT_AROUND_RIGHT );
				int edge2 = getEdge( edge1, NEXT_AROUND_RIGHT );
				
				edgeOrg(edge0, &org0);
				edgeDst(edge0, &dst0);
				edgeOrg(edge1, &org1);
				edgeDst(edge1, &dst1);
				
				vec2 virt_point = computeVoronoiPoint(org0, dst0, org1, dst1);
				
				if( fabs( virt_point.x ) < FLT_MAX * 0.5 &&
				   fabs( virt_point.y ) < FLT_MAX * 0.5 )
				{
					quadedge.pt[1] = qedges[edge1 >> 2].pt[1 + (edge1 & 2)] =
					qedges[edge2 >> 2].pt[1 + (edge2 & 2)] = newPoint(virt_point, true);
				}
			}
		}
		
		validGeometry = true;
	}
	
	
	static int
	isRightOf2( const vec2& pt, const vec2& org, const vec2& diff )
	{
		double cw_area = ((double)org.x - pt.x)*diff.y - ((double)org.y - pt.y)*diff.x;
		return (cw_area > 0) - (cw_area < 0);
	}
	
	
	int Subdiv2D::findNearest(vec2 pt, vec2* nearestPt)
	{
		if( !validGeometry )
			calcVoronoi();
		
		int vertex = 0, edge = 0;
		int loc = locate( pt, edge, vertex );
		
		if( loc != PTLOC_ON_EDGE && loc != PTLOC_INSIDE )
			return vertex;
		
		vertex = 0;
		
		vec2 start;
		edgeOrg(edge, &start);
		vec2 diff = pt - start;
		
		edge = rotateEdge(edge, 1);
		
		int i, total = (int)vtx.size();
		
		for( i = 0; i < total; i++ )
		{
			vec2 t;
			
			for(;;)
			{
				CI_ASSERT( edgeDst(edge, &t) > 0 );
				if( isRightOf2( t, start, diff ) >= 0 )
					break;
				
				edge = getEdge( edge, NEXT_AROUND_LEFT );
			}
			
			for(;;)
			{
				CI_ASSERT( edgeOrg( edge, &t ) > 0 );
				
				if( isRightOf2( t, start, diff ) < 0 )
					break;
				
				edge = getEdge( edge, PREV_AROUND_LEFT );
			}
			
			vec2 tempDiff;
			edgeDst(edge, &tempDiff);
			edgeOrg(edge, &t);
			tempDiff -= t;
			
			if( isRightOf2( pt, t, tempDiff ) >= 0 )
			{
				vertex = edgeOrg(rotateEdge( edge, 3 ));
				break;
			}
			
			edge = symEdge( edge );
		}
		
		if( nearestPt && vertex > 0 )
			*nearestPt = vtx[vertex].pt;
		
		return vertex;
	}
	
	void Subdiv2D::getEdgeList(std::vector<ci::vec4>& edgeList) const
	{
		edgeList.clear();
		
		for( size_t i = 4; i < qedges.size(); i++ )
		{
			if( qedges[i].isfree() )
				continue;
			if( qedges[i].pt[0] > 0 && qedges[i].pt[2] > 0 )
			{
				vec2 org = vtx[qedges[i].pt[0]].pt;
				vec2 dst = vtx[qedges[i].pt[2]].pt;
				edgeList.push_back(ci::vec4(org.x, org.y, dst.x, dst.y));
			}
		}
	}
	
	void Subdiv2D::getTriangleList(std::vector<Triangle>& triangleList) const
	{
		triangleList.clear();
		int i, total = (int)(qedges.size()*4);
		std::vector<bool> edgemask(total, false);
		
		for( i = 4; i < total; i += 2 )
		{
			if( edgemask[i] )
				continue;
			vec2 a, b, c;
			int edge = i;
			edgeOrg(edge, &a);
			edgemask[edge] = true;
			edge = getEdge(edge, NEXT_AROUND_LEFT);
			edgeOrg(edge, &b);
			edgemask[edge] = true;
			edge = getEdge(edge, NEXT_AROUND_LEFT);
			edgeOrg(edge, &c);
			edgemask[edge] = true;
			triangleList.push_back( Triangle( a, b, c ) );
		}
	}
	std::vector<uint32_t> Subdiv2D::getTriangleIndices()
	{
		std::vector<uint32_t> triangleIndices;
		int i, total = (int)(qedges.size()*4);
		std::vector<bool> edgemask(total, false);
		
		for( i = 4; i < total; i += 2 )
		{
			if( edgemask[i] )
				continue;
			vec2 a, b, c;
			ivec3 triangle;
			int edge = i;
			triangle.x = edgeOrg(edge, &a);
			edgemask[edge] = true;
			edge = getEdge(edge, NEXT_AROUND_LEFT);
			triangle.y = edgeOrg(edge, &b);
			edgemask[edge] = true;
			edge = getEdge(edge, NEXT_AROUND_LEFT);
			triangle.z = edgeOrg(edge, &c);
			edgemask[edge] = true;
			
			triangle.x = vtx[triangle.x].originalIndice;
			triangle.y = vtx[triangle.y].originalIndice;
			triangle.z = vtx[triangle.z].originalIndice;
			
			
			if( triangle.x != -1 && triangle.y != -1 && triangle.z != -1 ) {
				triangleIndices.push_back( triangle.z );
				triangleIndices.push_back( triangle.y );
				triangleIndices.push_back( triangle.x );
			}
		}
		
		return triangleIndices;
	}
	
	void Subdiv2D::getVoronoiFacetList(const std::vector<int>& idx,
									   std::vector<std::vector<vec2> >& facetList,
									   std::vector<vec2>& facetCenters)
	{
		calcVoronoi();
		facetList.clear();
		facetCenters.clear();
		
		std::vector<vec2> buf;
		
		size_t i, total;
		if( idx.empty() )
			i = 4, total = vtx.size();
		else
			i = 0, total = idx.size();
		
		for( ; i < total; i++ )
		{
			int k = idx.empty() ? (int)i : idx[i];
			
			if( vtx[k].isfree() || vtx[k].isvirtual() )
				continue;
			int edge = rotateEdge(vtx[k].firstEdge, 1), t = edge;
			
			// gather points
			buf.clear();
			do
			{
				buf.push_back(vtx[edgeOrg(t)].pt);
				t = getEdge( t, NEXT_AROUND_LEFT );
			}
			while( t != edge );
			
			facetList.push_back(buf);
			facetCenters.push_back(vtx[k].pt);
		}
	}
	
	
	void Subdiv2D::checkSubdiv() const
	{
		int i, j, total = (int)qedges.size();
		
		for( i = 0; i < total; i++ )
		{
			const QuadEdge& qe = qedges[i];
			
			if( qe.isfree() )
				continue;
			
			for( j = 0; j < 4; j++ )
			{
				int e = (int)(i*4 + j);
				int o_next = nextEdge(e);
				int o_prev = getEdge(e, PREV_AROUND_ORG );
				int d_prev = getEdge(e, PREV_AROUND_DST );
				int d_next = getEdge(e, NEXT_AROUND_DST );
				
				// check points
				CI_ASSERT( edgeOrg(e) == edgeOrg(o_next));
				CI_ASSERT( edgeOrg(e) == edgeOrg(o_prev));
				CI_ASSERT( edgeDst(e) == edgeDst(d_next));
				CI_ASSERT( edgeDst(e) == edgeDst(d_prev));
				
				if( j % 2 == 0 )
				{
					CI_ASSERT( edgeDst(o_next) == edgeOrg(d_prev));
					CI_ASSERT( edgeDst(o_prev) == edgeOrg(d_next));
					CI_ASSERT( getEdge(getEdge(getEdge(e,NEXT_AROUND_LEFT),NEXT_AROUND_LEFT),NEXT_AROUND_LEFT) == e );
					CI_ASSERT( getEdge(getEdge(getEdge(e,NEXT_AROUND_RIGHT),NEXT_AROUND_RIGHT),NEXT_AROUND_RIGHT) == e);
				}
			}
		}
	}
} // anonymous namespace

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

namespace Delaunay {
	
	std::vector<uint32_t> getTriangleIndices( const ci::Rectf &rect, const std::vector<ci::vec2> &points )
	{
		Subdiv2D delaunay( rect, points );
		return delaunay.getTriangleIndices();
	}
};