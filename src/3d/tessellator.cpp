#include "tessellator.h"

#include "qgscurve.h"
#include "qgspoint.h"
#include "qgspolygon.h"

#include "poly2tri/poly2tri.h"

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


Tessellator::Tessellator( double originX, double originY, bool addNormals )
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

void Tessellator::addPolygon( const QgsPolygonV2 &polygon, float extrusionHeight )
{
  const QgsCurve *exterior = polygon.exteriorRing();

  QList< std::vector<p2t::Point *> > polylinesToDelete;
//  QHash<p2t::Point *, float> z;

  std::vector<p2t::Point *> polyline;
  polyline.reserve( exterior->numPoints() );

  QgsVertexId::VertexType vt;
  QgsPoint pt;

  // Find the best fitting plane
  std::list<Kernel::Point_3> pointsInPolygon;
  for ( int i = 0; i < exterior->numPoints() - 1; ++i )
  {
    exterior->pointAt( i, pt, vt );

    pointsInPolygon.push_back(Kernel::Point_3(pt.x() - originX, pt.y() - originY, pt.z()));
  }

  for ( int i = 0; i < polygon.numInteriorRings(); ++i )
  {
    std::vector<p2t::Point *> holePolyline;
    holePolyline.reserve( exterior->numPoints() );
    const QgsCurve *hole = polygon.interiorRing( i );
    for ( int j = 0; j < hole->numPoints() - 1; ++j )
    {
      hole->pointAt( j, pt, vt );
      pointsInPolygon.push_back(Kernel::Point_3(pt.x() - originX, pt.y() - originY, pt.z()));
    }
  }

  Kernel::Plane_3 bestPlane;
  linear_least_squares_fitting_3(pointsInPolygon.begin(), pointsInPolygon.end(), bestPlane, CGAL::Dimension_tag<0>());

  // Triangulate the projection of the edges to the plane
  Triangulation triangulation;
  QgsPointSequence listOfPoints;
  exterior->points(listOfPoints);
  QList<QgsPoint>::iterator currentPoint = listOfPoints.begin();
  Triangulation::Vertex_handle currentVertex = triangulation.insert(bestPlane.to_2d(Kernel::Point_3(currentPoint->x() - originX, currentPoint->y() - originY, currentPoint->z())));
  ++currentPoint;
  Triangulation::Vertex_handle previousVertex;
  while (currentPoint != listOfPoints.end()) {
    previousVertex = currentVertex;
    currentVertex = triangulation.insert(bestPlane.to_2d(Kernel::Point_3(currentPoint->x() - originX, currentPoint->y() - originY, currentPoint->z())));
    if (previousVertex != currentVertex) triangulation.insert_constraint(previousVertex, currentVertex);
    ++currentPoint;
  }

  for ( int i = 0; i < polygon.numInteriorRings(); ++i )
  {
    const QgsCurve *hole = polygon.interiorRing( i );
    hole->points(listOfPoints);
    if (listOfPoints.size() < 4) {
      std::cout << "\tRing with < 4 points! Skipping..." << std::endl;
      continue;
    }

    currentPoint = listOfPoints.begin();
    currentVertex = triangulation.insert(bestPlane.to_2d(Kernel::Point_3(currentPoint->x() - originX, currentPoint->y() - originY, currentPoint->z())));
    while (currentPoint != listOfPoints.end()) {
      previousVertex = currentVertex;
      currentVertex = triangulation.insert(bestPlane.to_2d(Kernel::Point_3(currentPoint->x() - originX, currentPoint->y() - originY, currentPoint->z())));
      if (previousVertex != currentVertex) triangulation.insert_constraint(previousVertex, currentVertex);
      ++currentPoint;
    }
  }

  // Label the triangles to find out interior/exterior
  if (triangulation.number_of_faces() == 0) {
//      std::cout << "Degenerate face produced no triangles. Skipping..." << std::endl;
    return;
  }

  for (Triangulation::All_faces_iterator currentFace = triangulation.all_faces_begin(); currentFace != triangulation.all_faces_end(); ++currentFace) {
    currentFace->info() = std::pair<bool, bool>(false, false);
  }

  std::list<Triangulation::Face_handle> toCheck;
  triangulation.infinite_face()->info() = std::pair<bool, bool>(true, false);
  CGAL_assertion(triangulation.infinite_face()->info().first == true);
  CGAL_assertion(triangulation.infinite_face()->info().second == false);
  toCheck.push_back(triangulation.infinite_face());

  while (!toCheck.empty()) {
    CGAL_assertion(toCheck.front()->info().first);
    for (int neighbour = 0; neighbour < 3; ++neighbour) {
      if (toCheck.front()->neighbor(neighbour)->info().first) {
        // Note: validation code. But here we assume that some triangulations will be invalid anyway.
//          if (triangulation.is_constrained(Triangulation::Edge(toCheck.front(), neighbour))) CGAL_assertion(toCheck.front()->neighbor(neighbour)->info().second != toCheck.front()->info().second);
//          else CGAL_assertion(toCheck.front()->neighbor(neighbour)->info().second == toCheck.front()->info().second);
      } else {
        toCheck.front()->neighbor(neighbour)->info().first = true;
        CGAL_assertion(toCheck.front()->neighbor(neighbour)->info().first == true);
        if (triangulation.is_constrained(Triangulation::Edge(toCheck.front(), neighbour))) {
          toCheck.front()->neighbor(neighbour)->info().second = !toCheck.front()->info().second;
          toCheck.push_back(toCheck.front()->neighbor(neighbour));
        } else {
          toCheck.front()->neighbor(neighbour)->info().second = toCheck.front()->info().second;
          toCheck.push_back(toCheck.front()->neighbor(neighbour));
        }
      }
    }

    toCheck.pop_front();
  }

  // Project the triangles back to 3D and add
  for (Triangulation::Finite_faces_iterator currentFace = triangulation.finite_faces_begin(); currentFace != triangulation.finite_faces_end(); ++currentFace) {
    if (currentFace->info().second) {
      for (unsigned int currentVertexIndex = 0; currentVertexIndex < 3; ++currentVertexIndex) {
        Kernel::Point_3 point3 = bestPlane.to_3d(currentFace->vertex(currentVertexIndex)->point());
        data << point3.x() << point3.z() + 0.01 << -point3.y();
        if ( addNormals ) {
          data << bestPlane.orthogonal_vector().x() << bestPlane.orthogonal_vector().z() << -bestPlane.orthogonal_vector().y();
        }
      }
    }
  }

  // add walls if extrusion is enabled
  if ( extrusionHeight != 0 )
  {
    _makeWalls( *exterior, false, extrusionHeight, data, addNormals, originX, originY );

    for ( int i = 0; i < polygon.numInteriorRings(); ++i )
      _makeWalls( *polygon.interiorRing( i ), true, extrusionHeight, data, addNormals, originX, originY );
  }
}
