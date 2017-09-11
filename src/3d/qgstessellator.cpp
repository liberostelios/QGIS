#include "qgstessellator.h"

#include "qgscurve.h"
#include "qgspoint.h"
#include "qgspolygon.h"

#include "poly2tri/poly2tri.h"

#include <gsl/gsl_multifit.h>
#include <gsl/gsl_blas.h>

#include <QtDebug>

#include <QVector3D>

#include <iostream>

static void make_quad( float x0, float y0, float x1, float y1, float zLow, float zHigh, QVector<float> &data, bool addNormals )
{
  float dx = x1 - x0;
  float dy = -( y1 - y0 );

  // perpendicular vector in plane to [x,y] is [-y,x]
  QVector3D vn( -dy, 0, dx );
  vn.normalize();

  // triangle 1
  data << x0 << zHigh << -y0;
  if ( addNormals )
    data << vn.x() << vn.y() << vn.z();
  data << x1 << zHigh << -y1;
  if ( addNormals )
    data << vn.x() << vn.y() << vn.z();
  data << x0 << zLow  << -y0;
  if ( addNormals )
    data << vn.x() << vn.y() << vn.z();

  // triangle 2
  data << x0 << zLow  << -y0;
  if ( addNormals )
    data << vn.x() << vn.y() << vn.z();
  data << x1 << zHigh << -y1;
  if ( addNormals )
    data << vn.x() << vn.y() << vn.z();
  data << x1 << zLow  << -y1;
  if ( addNormals )
    data << vn.x() << vn.y() << vn.z();
}


QgsTessellator::QgsTessellator( double originX, double originY, bool addNormals )
  : originX( originX )
  , originY( originY )
  , addNormals( addNormals )
{
  stride = 3 * sizeof( float );
  if ( addNormals )
    stride += 3 * sizeof( float );
}


static bool _isRingCounterClockWise( const QgsCurve &ring )
{
  double a = 0;
  int count = ring.numPoints();
  QgsVertexId::VertexType vt;
  QgsPoint pt, ptPrev;
  ring.pointAt( 0, ptPrev, vt );
  for ( int i = 1; i < count + 1; ++i )
  {
    ring.pointAt( i % count, pt, vt );
    a += ptPrev.x() * pt.y() - ptPrev.y() * pt.x();
    ptPrev = pt;
  }
  return a > 0; // clockwise if a is negative
}

static void _makeWalls( const QgsCurve &ring, bool ccw, float extrusionHeight, QVector<float> &data, bool addNormals, double originX, double originY )
{
  // we need to find out orientation of the ring so that the triangles we generate
  // face the right direction
  // (for exterior we want clockwise order, for holes we want counter-clockwise order)
  bool is_counter_clockwise = _isRingCounterClockWise( ring );

  QgsVertexId::VertexType vt;
  QgsPoint pt;

  QgsPoint ptPrev;
  ring.pointAt( is_counter_clockwise == ccw ? 0 : ring.numPoints() - 1, ptPrev, vt );
  for ( int i = 1; i < ring.numPoints(); ++i )
  {
    ring.pointAt( is_counter_clockwise == ccw ? i : ring.numPoints() - i - 1, pt, vt );
    float x0 = ptPrev.x() - originX, y0 = ptPrev.y() - originY;
    float x1 = pt.x() - originX, y1 = pt.y() - originY;
    float height = pt.z();
    // make a quad
    make_quad( x0, y0, x1, y1, height, height + extrusionHeight, data, addNormals );
    ptPrev = pt;
  }
}

static QgsPoint cross_product( const QgsPoint u, const QgsPoint v)
{
  QgsPoint result;

  result.setX(u.y() * v.z() - u.z() * v.y());
  result.setY(u.z() * v.x() - u.x() * v.z());
  result.addZValue(u.x() * v.y() - u.y() * v.x());

  return result;
}

static double dot_product( const QgsPoint u, const QgsPoint v)
{
    return u.x() * v.x() + u.y() * v.y() + u.z() * v.z();
}

static QgsPoint multiply( const QgsPoint p, const double factor )
{
    return QgsPoint(p.x() * factor, p.y() * factor, p.z() * factor);
}

static QgsPoint add( const QgsPoint p1, const QgsPoint p2)
{
    return QgsPoint(p1.x() + p2.x(), p1.y() + p2.y(), p1.z() + p2.z());
}

void QgsTessellator::addPolygon( const QgsPolygonV2 &polygon, float extrusionHeight )
{
  const QgsCurve *exterior = polygon.exteriorRing();

  QList< std::vector<p2t::Point *> > polylinesToDelete;
  QHash<p2t::Point *, double> z;

  std::vector<p2t::Point *> polyline;
  polyline.reserve( exterior->numPoints() );

  QgsVertexId::VertexType vt;
  QgsPoint pt, ptPrev;

  // Compute the best fitting plane
  gsl_matrix *X, *cov;
  gsl_vector *y, *c;
  double chisq;

  X = gsl_matrix_alloc(exterior->numPoints() - 1, 3);
  y = gsl_vector_alloc(exterior->numPoints() - 1);

  c = gsl_vector_alloc(3);
  cov = gsl_matrix_alloc(3, 3);

  for ( int i = 0; i < exterior->numPoints() - 1; ++i )
  {
    exterior->pointAt( i, pt, vt );
    gsl_matrix_set(X, i, 0, pt.x() - originX);
    gsl_matrix_set(X, i, 1, pt.y() - originY);
    gsl_matrix_set(X, i, 2, 1);

    gsl_vector_set(y, i, qIsNaN( pt.z() ) ? 0 : pt.z());
  }

  gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc(exterior->numPoints() - 1, 3);
  gsl_multifit_linear(X, y, c, cov, &chisq, work);
  gsl_multifit_linear_free(work);

  double pa, pb, pc, pd;

#define C(i) (gsl_vector_get(c, (i)))

  pa = C(0);
  pb = C(1);
  pc = -1;
  pd = C(2);

  QgsPoint pNormal(pa, -pb, pc), pOrigin(0, 0, pd), pXVector(1, 0, pa);
  QgsPoint pYVector = cross_product(pNormal, pXVector);
  pYVector.setY(-pYVector.y());
  pNormal = multiply(pNormal, 1.0f / pNormal.distance3D(0, 0, 0));
  pXVector = multiply(pXVector, 1.0f / pXVector.distance3D(0, 0, 0));
  pYVector = multiply(pYVector, 1.0f / pYVector.distance3D(0, 0, 0));
  p2t::Point *ptFirst = new p2t::Point();

  std::cout << "Doing cdt for: ";

  for ( int i = 0; i < exterior->numPoints() - 1; ++i )
  {
    exterior->pointAt( i, pt, vt );
    if ( i == 0 || pt != ptPrev )
    {
      QgsPoint tempPt( pt.x() - originX - pOrigin.x(), pt.y() - originY - pOrigin.y(), (qIsNaN( pt.z() ) ? 0 : pt.z()) - pOrigin.z() );
      double x = dot_product(tempPt, pXVector);
      double y = dot_product(tempPt, pYVector);

      p2t::Point *pt2 = new p2t::Point( x, y );
      bool found = false;
      for (std::vector<p2t::Point *>::iterator it = polyline.begin(); it != polyline.end(); it++)
      {
          if (*pt2 == **it)
          {
              found = true;
              std::cout << "XXXX - FOUND a duplicate!";
              return;
          }
      }
      if (pt2->x != ptFirst->x && pt2->y != ptFirst->y)
      {
        polyline.push_back( pt2 );
      }

      if ( i == 0 )
      {
        ptFirst = pt2;
      }

      float zPt = qIsNaN( pt.z() ) ? 0 : pt.z();
      z[pt2] = zPt;

      std::cout << tempPt << " | ";
    }
    ptPrev = pt;
  }
  polylinesToDelete << polyline;
  std::cout << std::endl << std::flush;

  p2t::CDT *cdt = new p2t::CDT( polyline );

//  // polygon holes
//  for ( int i = 0; i < polygon.numInteriorRings(); ++i )
//  {
//    std::vector<p2t::Point *> holePolyline;
//    holePolyline.reserve( exterior->numPoints() );
//    const QgsCurve *hole = polygon.interiorRing( i );
//    for ( int j = 0; j < hole->numPoints() - 1; ++j )
//    {
//      hole->pointAt( j, pt, vt );
//      if ( j == 0 || pt != ptPrev )
//      {
//        p2t::Point *pt2 = new p2t::Point( pt.x() - originX, pt.y() - originY );
//        holePolyline.push_back( pt2 );
//        float zPt = qIsNaN( pt.z() ) ? 0 : pt.z();
//        z[pt2] = zPt;
//      }
//      ptPrev = pt;
//    }
//    cdt->AddHole( holePolyline );
//    polylinesToDelete << holePolyline;
//  }

  // TODO: robustness (no nearly duplicate points, invalid geometries ...)

  cdt->Triangulate();

  std::vector<p2t::Triangle *> triangles = cdt->GetTriangles();

  for ( size_t i = 0; i < triangles.size(); ++i )
  {
    p2t::Triangle *t = triangles[i];
    for ( int j = 0; j < 3; ++j )
    {
      p2t::Point *p = t->GetPoint( j );
      double zPt = z[p];
      QgsPoint nPoint = add(add(pOrigin, multiply(pXVector, p->x)), multiply(pYVector, p->y));
      double fx = nPoint.x();
      double fy = nPoint.y();
      double fz = extrusionHeight + (qIsNaN(zPt) ? 0 : zPt);
      data << fx << fz << -fy;
      if ( addNormals )
        data << 0.f << 1.f << 0.f;
    }
  }

  delete cdt;
  for ( int i = 0; i < polylinesToDelete.count(); ++i )
    qDeleteAll( polylinesToDelete[i] );

  // add walls if extrusion is enabled
  if ( extrusionHeight != 0 )
  {
    _makeWalls( *exterior, false, extrusionHeight, data, addNormals, originX, originY );

    for ( int i = 0; i < polygon.numInteriorRings(); ++i )
      _makeWalls( *polygon.interiorRing( i ), true, extrusionHeight, data, addNormals, originX, originY );
  }
}
