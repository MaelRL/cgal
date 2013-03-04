#ifndef SCENE_STARSET3_ITEM_H
#define SCENE_STARSET3_ITEM_H

#include "Scene_starset3_item_config.h"
#include "StarSet_type.h"
#include <CGAL_demo/Scene_item_with_display_list.h>

#include <QVector>
#include <QColor>
#include <set>

#include <QtCore/qglobal.h>
#include <CGAL/gl.h>
#include <QGLViewer/manipulatedFrame.h>
#include <QGLViewer/qglviewer.h>

struct Scene_starset3_item_priv;

class SCENE_STARSET3_ITEM_EXPORT Scene_starset3_item
  : public Scene_item_with_display_list
{
  Q_OBJECT
public:
  typedef qglviewer::ManipulatedFrame ManipulatedFrame;

  Scene_starset3_item(const Criteria& criteria,
                      const Metric& metric,
                      const Constrain_surface* const surface,
                      const int nb_initial_points);
  ~Scene_starset3_item();

  const Surface_star_set& star_set() const;
  Surface_star_set& star_set();
  
  bool manipulatable() const 
  {
    return true;
  }

  ManipulatedFrame* manipulatedFrame() 
  {
    return frame;
  }

  void setPosition(float x, float y, float z) 
  {
    frame->setPosition(x, y, z);
  }

  void setNormal(float x, float y, float z) 
  {
    frame->setOrientation(x, y, z, 0.f);
  }

  Kernel::Plane_3 plane() const;

  bool isFinite() const { return true; }
  bool isEmpty() const  { return star_set().empty(); }
  int nbStars() const { return (int)star_set().size(); }

  Bbox bbox() const;

  Scene_starset3_item* clone() const { return 0; }

  QString toolTip() const;
//  virtual QPixmap graphicalToolTip() const;

  // Indicate if rendering mode is supported
  bool supportsRenderingMode(RenderingMode m) const {
    return (m != Gouraud); // CHECK THIS!
  }

  void direct_draw() const;
  
  //io
  bool save(std::ofstream& out) const;

  // data item
  inline const Scene_item* data_item() const;
  inline void set_data_item(const Scene_item* data_item);
  
  // rebuild histogram
//  inline void update_histogram();
  
  // Call this if starset has been modified
  void starset_changed();

public slots:
  inline void data_item_destroyed();
  virtual void setColor(QColor c);
  
//private:
//  void build_histogram();
  void compute_color_map(const QColor& c);
//  QColor get_histogram_color(const double v) const;

public:
  int& draw_star()                 { return m_draw_star_id; }
  const int& draw_star() const     { return m_draw_star_id; }
  bool& draw_dual()                { return m_draw_dual; }
  bool& draw_poles()        { return m_draw_poles; }
  bool& draw_initial_points()        { return m_draw_initial_points; }
  bool& draw_surface_delaunay_balls() { return m_draw_surface_delaunay_balls; }
  bool& draw_inconsistent_facets() { return m_draw_inconsistent_facets; }
  bool& draw_metric_field()        { return m_draw_metric_field; }
  double draw_metric_eps() const { return m_draw_metric_eps; }
  bool& draw_mesh_3()              { return m_draw_mesh_3; }

protected:
  Scene_starset3_item_priv* d;

  qglviewer::ManipulatedFrame* frame;
  
private:
  //QPixmap histogram_;
  const Scene_item* data_item_;
  
  typedef std::set<int> Indices;
  Indices indices_;

  bool m_draw_dual;
  bool m_draw_poles;
  bool m_draw_initial_points;
  bool m_draw_surface_delaunay_balls;
  int m_draw_star_id;
  bool m_draw_inconsistent_facets;
  bool m_draw_metric_field;
  double m_draw_metric_eps;
  bool m_draw_mesh_3;
};

inline
const Scene_item*
Scene_starset3_item::data_item() const
{
  return data_item_;
}

inline
void
Scene_starset3_item::set_data_item(const Scene_item* data_item)
{
  data_item_ = data_item;
  
  if ( NULL != data_item )
  {
    connect(data_item, SIGNAL(aboutToBeDestroyed()),
            this, SLOT(data_item_destroyed()));
  }
}

//inline
//void
//Scene_starset3_item::update_histogram()
//{
//  build_histogram();
//}

inline
void
Scene_starset3_item::data_item_destroyed()
{
  set_data_item(NULL);
}

#endif // SCENE_STARSET3_ITEM_H
