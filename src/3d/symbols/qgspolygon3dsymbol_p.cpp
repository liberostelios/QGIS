#include "qgspolygon3dsymbol_p.h"

#include "qgspolygon3dsymbol.h"
#include "qgstessellatedpolygongeometry.h"
#include "qgs3dmapsettings.h"
#include "terraingenerator.h"
#include "qgs3dutils.h"

#include <Qt3DCore/QTransform>

#include "qgsvectorlayer.h"
#include "qgsmultipolygon.h"



QgsPolygon3DSymbolEntity::QgsPolygon3DSymbolEntity( const Qgs3DMapSettings &map, QgsVectorLayer *layer, const QgsPolygon3DSymbol &symbol, Qt3DCore::QNode *parent )
  : Qt3DCore::QEntity( parent )
{
  addEntityForSelectedPolygons( map, layer, symbol );
  addEntityForNotSelectedPolygons( map, layer, symbol );
}

void QgsPolygon3DSymbolEntity::addEntityForSelectedPolygons( const Qgs3DMapSettings &map, QgsVectorLayer *layer, const QgsPolygon3DSymbol &symbol )
{
  // build the default material
  Qt3DExtras::QPhongMaterial *mat = material( symbol );

  // update the material with selection colors
  mat->setDiffuse( map.selectionColor() );
  mat->setAmbient( map.selectionColor().darker() );

  // build a transform function
  Qt3DCore::QTransform *tform = new Qt3DCore::QTransform;
  tform->setTranslation( QVector3D( 0, 0, 0 ) );

  // build the feature request to select features
  QgsFeatureRequest req;
  req.setDestinationCrs( map.crs );
  req.setFilterFids( layer->selectedFeatureIds() );

  // build the entity
  QgsPolygon3DSymbolEntityNode *entity = new QgsPolygon3DSymbolEntityNode( map, layer, symbol, req );
  entity->addComponent( mat );
  entity->addComponent( tform );
  entity->setParent( this );
}

void QgsPolygon3DSymbolEntity::addEntityForNotSelectedPolygons( const Qgs3DMapSettings &map, QgsVectorLayer *layer, const QgsPolygon3DSymbol &symbol )
{
  // build the default material
  Qt3DExtras::QPhongMaterial *mat = material( symbol );

  // build a transform function
  Qt3DCore::QTransform *tform = new Qt3DCore::QTransform;
  tform->setTranslation( QVector3D( 0, 0, 0 ) );

  // build the feature request to select features
  QgsFeatureRequest req;
  req.setDestinationCrs( map.crs );

  QgsFeatureIds notSelected = layer->allFeatureIds();
  notSelected.subtract( layer->selectedFeatureIds() );
  req.setFilterFids( notSelected );

  // build the entity
  QgsPolygon3DSymbolEntityNode *entity = new QgsPolygon3DSymbolEntityNode( map, layer, symbol, req );
  entity->addComponent( mat );
  entity->addComponent( tform );
  entity->setParent( this );
}

Qt3DExtras::QPhongMaterial *QgsPolygon3DSymbolEntity::material( const QgsPolygon3DSymbol &symbol ) const
{
  Qt3DExtras::QPhongMaterial *material = new Qt3DExtras::QPhongMaterial;
  material->setAmbient( symbol.material().ambient() );
  material->setDiffuse( symbol.material().diffuse() );
  material->setSpecular( symbol.material().specular() );
  material->setShininess( symbol.material().shininess() );
  return material;
}

QgsPolygon3DSymbolEntityNode::QgsPolygon3DSymbolEntityNode( const Qgs3DMapSettings &map, QgsVectorLayer *layer, const QgsPolygon3DSymbol &symbol, const QgsFeatureRequest &req, Qt3DCore::QNode *parent )
  : Qt3DCore::QEntity( parent )
{
  addComponent( renderer( map, symbol, layer, req ) );
}

Qt3DRender::QGeometryRenderer *QgsPolygon3DSymbolEntityNode::renderer( const Qgs3DMapSettings &map, const QgsPolygon3DSymbol &symbol, const QgsVectorLayer *layer, const QgsFeatureRequest &request )
{
  QgsPointXY origin( map.originX, map.originY );
  QList<QgsPolygonV2 *> polygons;
  QgsFeature f;
  QgsFeatureIterator fi = layer->getFeatures( request );
  while ( fi.nextFeature( f ) )
  {
    if ( f.geometry().isNull() )
      continue;

    QgsGeometry geom = f.geometry();

    // segmentize curved geometries if necessary
    if ( QgsWkbTypes::isCurvedType( geom.geometry()->wkbType() ) )
      geom = QgsGeometry( geom.geometry()->segmentize() );

//    if ( !geom.isGeosValid() )
//    {
//      // invalid geometries break tessellation
//      qDebug() << "skipping invalid geometry" << f.id();
//      continue;
//    }

    QgsAbstractGeometry *g = geom.geometry();

    if ( QgsWkbTypes::flatType( g->wkbType() ) == QgsWkbTypes::Polygon )
    {
      QgsPolygonV2 *poly = static_cast<QgsPolygonV2 *>( g );
      QgsPolygonV2 *polyClone = poly->clone();
      Qgs3DUtils::clampAltitudes( polyClone, symbol.altitudeClamping(), symbol.altitudeBinding(), symbol.height(), map );
      polygons.append( polyClone );
    }
    else if ( QgsWkbTypes::flatType( g->wkbType() ) == QgsWkbTypes::MultiPolygon )
    {
      QgsMultiPolygonV2 *mpoly = static_cast<QgsMultiPolygonV2 *>( g );
      for ( int i = 0; i < mpoly->numGeometries(); ++i )
      {
        QgsAbstractGeometry *g2 = mpoly->geometryN( i );
        Q_ASSERT( QgsWkbTypes::flatType( g2->wkbType() ) == QgsWkbTypes::Polygon );
        QgsPolygonV2 *polyClone = static_cast<QgsPolygonV2 *>( g2 )->clone();
        Qgs3DUtils::clampAltitudes( polyClone, symbol.altitudeClamping(), symbol.altitudeBinding(), symbol.height(), map );
        polygons.append( polyClone );
      }
    }
    else
      qDebug() << "not a polygon";
  }

  mGeometry = new QgsTessellatedPolygonGeometry;
  mGeometry->setPolygons( polygons, origin, symbol.extrusionHeight() );

  Qt3DRender::QGeometryRenderer *renderer = new Qt3DRender::QGeometryRenderer;
  renderer->setGeometry( mGeometry );

  return renderer;
}
