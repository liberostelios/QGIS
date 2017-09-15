#include "qgstessellator.h"

#include "qgscurve.h"
#include "qgspoint.h"
#include "qgspolygon.h"

#include "poly2tri/poly2tri.h"

#include <gsl/gsl_multifit.h>
#include <gsl/gsl_blas.h>

#include <QtDebug>

#include <QVector3D>

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

void QgsTessellator::addPolygon( const QgsPolygonV2 &polygon, float extrusionHeight )
{
  const QgsCurve *exterior = polygon.exteriorRing();

  QList< std::vector<p2t::Point *> > polylinesToDelete;
  QHash<p2t::Point *, double> z;

  std::vector<p2t::Point *> polyline;
  polyline.reserve( exterior->numPoints() );

  QgsVertexId::VertexType vt;
  QgsPoint pt, ptPrev, ptFirst;

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
    gsl_matrix_set(X, i, 0, pt.x());
    gsl_matrix_set(X, i, 1, pt.y());
    gsl_matrix_set(X, i, 2, 1);

    gsl_vector_set(y, i, qIsNaN( pt.z() ) ? 0 : pt.z());

    if (i == 0)
    {
      ptFirst = pt;
    }
  }

  gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc(exterior->numPoints() - 1, 3);
  gsl_multifit_linear(X, y, c, cov, &chisq, work);

  bool switchPlane = chisq > 0.001;

  if (switchPlane)
  {
      for ( int i = 0; i < exterior->numPoints() - 1; ++i )
      {
        exterior->pointAt( i, pt, vt );
        gsl_matrix_set(X, i, 0, pt.x());
        gsl_matrix_set(X, i, 1, qIsNaN( pt.z() ) ? 0 : pt.z());
        gsl_matrix_set(X, i, 2, 1);

        gsl_vector_set(y, i, pt.y());
      }

      gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc(exterior->numPoints() - 1, 3);
      gsl_multifit_linear(X, y, c, cov, &chisq, work);
  }
  gsl_multifit_linear_free(work);

  double pa, pb, pc, pd;

#define C(i) (gsl_vector_get(c, (i)))

  if (switchPlane)
  {
      pa = C(0);
      pb = -1;
      pc = C(1);
      pd = C(2);
  }
  else
  {
      pa = C(0);
      pb = C(1);
      pc = -1;
      pd = C(2);
  }

  QVector3D pNormal(pa, pb, pc), pOrigin(ptFirst.x(), ptFirst.y(), ptFirst.z()), pXVector;
  if (pc > 0.001 || pc < -0.001)
  {
      pXVector = QVector3D(1, 0, -pa/pc);
  }
  else
  {
      pXVector = QVector3D(1, -pa/pb, 0);
  }
  QVector3D pYVector = QVector3D::normal(pNormal, pXVector);
  pYVector.setY(pYVector.y());
  pNormal.normalize();
  pXVector.normalize();

  for ( int i = 0; i < exterior->numPoints() - 1; ++i )
  {
    exterior->pointAt( i, pt, vt );
    if ( i == 0 || pt != ptPrev )
    {
      QVector3D tempPt( pt.x(), pt.y(), (qIsNaN( pt.z() ) ? 0 : pt.z()) );
      float x = QVector3D::dotProduct(tempPt - pOrigin, pXVector);
      float y = QVector3D::dotProduct(tempPt - pOrigin, pYVector);

      p2t::Point *pt2 = new p2t::Point( x, y );
      bool found = false;
      for (std::vector<p2t::Point *>::iterator it = polyline.begin(); it != polyline.end(); it++)
      {
          if (*pt2 == **it)
          {
              found = true;
          }
      }

      if (found)
      {
          continue;
      }

      polyline.push_back(pt2);

      float zPt = qIsNaN( pt.z() ) ? 0 : pt.z();
      z[pt2] = zPt;

    }
    ptPrev = pt;
  }
  polylinesToDelete << polyline;

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

  if (polyline.size() == 3)
  {
      for (std::vector<p2t::Point*>::iterator it = polyline.begin(); it != polyline.end(); it++)
      {
          p2t::Point *p = *it;
          double zPt = z[p];
          QVector3D nPoint = pOrigin + pXVector * p->x + pYVector * p->y;
          double fx = nPoint.x() - originX;
          double fy = nPoint.y() - originY;
          double fz = extrusionHeight + (qIsNaN(zPt) ? 0 : zPt);
          data << fx << fz << -fy;
          if ( addNormals )
            data << pNormal.x() << pNormal.z() << - pNormal.y();
      }

      return;
  }
  else if (polyline.size() < 3)
  {
      return;
  }

  cdt->Triangulate();

  std::vector<p2t::Triangle *> triangles = cdt->GetTriangles();

  for ( size_t i = 0; i < triangles.size(); ++i )
  {
    p2t::Triangle *t = triangles[i];
    for ( int j = 0; j < 3; ++j )
    {
      p2t::Point *p = t->GetPoint( j );
      double zPt = z[p];
      QVector3D nPoint = pOrigin + pXVector * p->x + pYVector * p->y;
      double fx = nPoint.x() - originX;
      double fy = nPoint.y() - originY;
      double fz = extrusionHeight + (qIsNaN(zPt) ? 0 : zPt);
      data << fx << fz << -fy;
      if ( addNormals )
        data << -pNormal.x() << -pNormal.z() << pNormal.y();
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
