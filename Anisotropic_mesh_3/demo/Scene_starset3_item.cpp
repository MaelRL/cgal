#include "Scene_starset3_item.h"

#include <QVector>
#include <QColor>
#include <QPixmap>
#include <QPainter>

#include <map>
#include <vector>
#include <CGAL/gl.h>
#include <CGAL/Mesh_3/dihedral_angle_3.h>

#include <CGAL_demo/Scene_item_with_display_list.h>
#include <CGAL_demo/Scene_interface.h>
#include <QtCore/qglobal.h>
#include <CGAL/gl.h>
#include <QGLViewer/manipulatedFrame.h>
#include <QGLViewer/qglviewer.h>

namespace {
  void CGALglcolor(QColor c, int dv = 0)
  {
    if ( 0 != dv )
    {
// workaround for Qt-4.2.
#if QT_VERSION < 0x040300
#  define darker dark
#endif
      c = c.darker(dv);
#undef darker
    }
    
    ::glColor4f(c.red()/255.0, c.green()/255.0, c.blue()/255.0, c.alpha()/255.0);
  }
}

//template<typename C3t3>
//std::vector<int>
//create_histogram(const C3t3& c3t3, double& min_value, double& max_value);


double complex_diag(const Scene_item* item) 
{
  const Scene_item::Bbox& bbox = item->bbox();
  const double& xdelta = bbox.xmax-bbox.xmin;
  const double& ydelta = bbox.ymax-bbox.ymin;
  const double& zdelta = bbox.zmax-bbox.zmin;
  const double diag = std::sqrt(xdelta*xdelta +
                                ydelta*ydelta +
                                zdelta*zdelta);
  return diag * 0.7;
}


struct Scene_starset3_item_priv 
{
  Scene_starset3_item_priv(const Criteria& criteria,
                           const Metric& metric,
                           const Constrain_surface* const surface,
                           const int nb_initial_points)
    : star_set(criteria, metric, surface, nb_initial_points)
  { }

  Surface_star_set star_set;
  QVector<QColor> colors;
};


Scene_starset3_item::
Scene_starset3_item(const Criteria& criteria,
                    const Metric& metric,
                    const Constrain_surface* const surface,
                    const int nb_initial_points)
  : d(new Scene_starset3_item_priv(criteria, metric, surface, nb_initial_points)),
  frame(new ManipulatedFrame()),
  //histogram_(),
  data_item_(NULL),
  indices_(),
  m_draw_surface_star_set(true),
  m_draw_cell(false),
  m_draw_dual(false),
  m_draw_poles(false),
  m_draw_initial_points(false),
  m_draw_surface_delaunay_balls(false),
  m_draw_star_id(-1),
  m_draw_inconsistent_facets(false),
  m_draw_metric_field(false),
  m_draw_metric_eps(metric.epsilon),
  m_draw_mesh_3(false)
{
  connect(frame, SIGNAL(modified()), this, SLOT(changed()));
  starset_changed();
}

Scene_starset3_item::~Scene_starset3_item()
{
  delete frame;
  delete d;
}

const Surface_star_set& 
Scene_starset3_item::star_set() const {
  return d->star_set;
}

Surface_star_set& 
Scene_starset3_item::star_set()
{
  return d->star_set;
}

Kernel::Plane_3 
Scene_starset3_item::plane() const 
{
  const qglviewer::Vec& pos = frame->position();
  const qglviewer::Vec& n =
    frame->inverseTransformOf(qglviewer::Vec(0.f, 0.f, 1.f));
  return Kernel::Plane_3(n[0], n[1],  n[2], - n * pos);
}

Scene_item::Bbox 
Scene_starset3_item::bbox() const
{
  if(star_set().empty())
    return Bbox();
  else 
  {
    bool initialized = false;
    CGAL::Bbox_3 result;
    std::size_t N = star_set().m_stars.size();
    for(unsigned int i = 0; i < N; i++)
    {
      if(star_set().m_stars[i]->is_surface_star())
      {
        if(!initialized){
          result = star_set().m_stars[i]->center_point().bbox();
          initialized = true;
        }
        else
          result = result + star_set().m_stars[i]->center_point().bbox();
      }
    }
    return Bbox(result.xmin(), result.ymin(), result.zmin(),
                result.xmax(), result.ymax(), result.zmax());
  }
}

bool
Scene_starset3_item::save(std::ofstream& out) const
{
  (const_cast< Scene_starset3_item* >(this))->star_set().output(out, true);
  return !out.bad();
}

QString 
Scene_starset3_item::toolTip() const 
{
  return tr("<p>3D cell star set: <br />"
            "<b>%5</b></p>"
            "<p>Number of surface vertices : %1<br />"
            "Total number of vertices (memory): %2<br />"
            "Number of surface facets: %3<br />"
            "Number of volume tetrahedra: %4</p>")
  .arg(star_set().number_of_surface_stars())
  .arg(star_set().total_number_of_vertices())
  .arg(star_set().count_restricted_facets())
  .arg(star_set().number_of_tets_in_star_set())
  .arg(this->name());
}

void
Scene_starset3_item::direct_draw() const 
{
  ::glPushMatrix();
  ::glMultMatrixd(frame->matrix());
  QGLViewer::drawGrid((float)complex_diag(this));
  ::glPopMatrix();

  if(isEmpty())
    return;

  CGALglcolor(QColor(0,0,0));

  GLboolean lighting = ::glIsEnabled(GL_LIGHTING);
  GLboolean two_side;
  ::glGetBooleanv(GL_LIGHT_MODEL_TWO_SIDE, &two_side);
  ::glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
  if(!lighting)
    ::glEnable(GL_LIGHTING);

  const Kernel::Plane_3& plane = this->plane();

  if(m_draw_surface_star_set)
    star_set().gl_draw(plane, true/*draw_edges*/, m_draw_star_id);
  if(m_draw_cell)
    star_set().gl_draw_cell(plane, m_draw_star_id);
  if(m_draw_poles)
    star_set().gl_draw_poles(plane);
  if(m_draw_initial_points)
    star_set().gl_draw_initial_points(plane);
  if(m_draw_dual)
    star_set().gl_draw_dual(m_draw_star_id);
  if(m_draw_surface_delaunay_balls)
    star_set().gl_draw_surface_delaunay_balls(plane, m_draw_star_id);
  if(m_draw_inconsistent_facets)
    star_set().gl_draw_inconsistent_facets(m_draw_star_id);
  if(m_draw_metric_field)
  {
    Scene_item::Bbox mf_bbox = this->bbox();
    double bbox_min = (std::min)(mf_bbox.depth(), mf_bbox.height());
    bbox_min = (std::min)(bbox_min, mf_bbox.width());
    star_set().gl_draw_metric(plane, bbox_min, draw_metric_eps(), m_draw_star_id);
  }
  if(m_draw_mesh_3)
    star_set().constrain_surface()->gl_draw_intermediate_mesh_3(plane);
  
  if(!two_side)
    ::glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);
  if(lighting)
    ::glEnable(GL_LIGHTING);
  else
    ::glDisable(GL_LIGHTING);
}


/*
QPixmap
Scene_starset3_item::graphicalToolTip() const
{
  if ( ! histogram_.isNull() )
  {
    return histogram_;
  }
  else
  {
    const_cast<Scene_c3t3_item&>(*this).build_histogram();
    return histogram_;
  }
}*/


void
Scene_starset3_item::setColor(QColor c)
{
  color_ = c;
  compute_color_map(c);
}

void
Scene_starset3_item::setSelectedPoint(double x, double y, double z)
{
  std::cerr << "Scene_starset3_item " << qPrintable(this->name())
            << " received point(" << x << ", " << y << ", " << z << ")\n";

  typedef Surface_star_set::Kd_traits TreeTraits;
  typedef CGAL::Orthogonal_k_neighbor_search<TreeTraits> Neighbor_search;
  typedef Surface_star_set::Kd_tree::Star_pmap Star_pmap;

  Neighbor_search search(this->star_set().m_kd_tree
                         , Kernel::Point_3(x, y, z)
                         , 1 // N
                         , 0 // epsilon
                         , true // search nearest
                         , this->star_set().m_kd_tree.traits().point_property_map()
                         );

  if(search.begin() != search.end()) {
    std::cerr << "Picked start #" << search.begin()->first << std::endl;
  }
}

void
Scene_starset3_item::starset_changed()
{
  // Update colors
  // Fill indices map and get max subdomain value
  indices_.clear();
  
  std::size_t max = d->star_set.size();
/*  for(C3t3::Cells_in_complex_iterator cit = this->c3t3().cells_in_complex_begin(),
      end = this->c3t3().cells_in_complex_end() ; cit != end; ++cit)
  {
    max = (std::max)(max, cit->subdomain_index());
    indices_.insert(cit->subdomain_index());
  }
  */
  d->colors.resize((int)(max+1));
  compute_color_map(color_);
  
  // Rebuild histogram
  //build_histogram();
  this->changed();
}


void
Scene_starset3_item::compute_color_map(const QColor& c)
{
  std::size_t nb_stars = d->star_set.size();
  for(std::size_t i = 0; i < nb_stars; i++)
  {
    double hue = c.hueF() + 1./nb_stars * i;
    if(hue > 1)
      hue -= 1.;
    d->colors[i] = QColor::fromHsvF(hue, c.saturationF(), c.valueF());
  }
}

#include "Scene_starset3_item.moc"
