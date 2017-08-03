#ifndef TESSELLATOR_H
#define TESSELLATOR_H

class QgsPolygonV2;

#include <QVector>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_3.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel cgalK;
typedef CGAL::Triangulation_3<cgalK>      cgalTriangulation;
typedef cgalTriangulation::Finite_facets_iterator cgalFinite_facets_iterator;
typedef cgalTriangulation::Locate_type    cgalLocate_type;
typedef cgalTriangulation::Point          cgalPoint;
typedef cgalTriangulation::Triangle       cgalTriangle;

class Tessellator
{
  public:
    Tessellator( double originX, double originY, bool addNormals );

    void addPolygon( const QgsPolygonV2 &polygon, float extrusionHeight );

    // input:
    // - origin X/Y
    // - whether to add walls
    // - stream of geometries
    // output:
    // - vertex buffer data (+ index buffer data ?)

    double originX, originY;
    bool addNormals;
    //QByteArray data;
    QVector<float> data;
    int stride;  //!< Size of one vertex entry in bytes
};

#endif // TESSELLATOR_H
