#ifndef TESSELLATOR_H
#define TESSELLATOR_H

class QgsPolygonV2;

#include <QVector>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_3.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/linear_least_squares_fitting_3.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Exact_predicates_tag Tag;
typedef CGAL::Triangulation_vertex_base_2<Kernel> VertexBase;
typedef CGAL::Constrained_triangulation_face_base_2<Kernel> FaceBase;
typedef CGAL::Triangulation_face_base_with_info_2<std::pair<bool, bool>, Kernel, FaceBase> FaceBaseWithInfo;
typedef CGAL::Triangulation_data_structure_2<VertexBase, FaceBaseWithInfo> TriangulationDataStructure;
typedef CGAL::Constrained_Delaunay_triangulation_2<Kernel, TriangulationDataStructure, Tag> Triangulation;

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
