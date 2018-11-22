#include <QApplication>
#include <QMainWindow>
#include <QTime>
#include <QAction>
#include <QObject>
#include <QDockWidget>

#include "SMesh_type.h"
#include "Scene.h"
#include <CGAL/Three/Scene_group_item.h>
#include "Scene_surface_mesh_item.h"
#include "Scene_polygon_soup_item.h"
#include "Scene_polylines_item.h"
#include "Messages_interface.h"

#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include <CGAL/convex_hull_2.h>

#include <map>
#include <vector>
#include <cstdlib>
#include <ctime>

#include "ui_Surface_mesh_approximation_dockwidget.h"
#include "VSA_wrapper.h"

using namespace CGAL::Three;

class Polyhedron_demo_surface_mesh_approximation_plugin :
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

  typedef boost::property_map<SMesh, CGAL::face_patch_id_t<int> >::type Face_id_map;
  struct Patch_id_pmap : public boost::put_get_helper<std::size_t, Patch_id_pmap>
  {
  public:
    typedef boost::read_write_property_map_tag category;
    typedef std::size_t value_type;
    typedef int& reference;
    typedef face_descriptor key_type;

    Patch_id_pmap(Face_id_map fid_map): m_fid_map(fid_map) {}

    friend void put(Patch_id_pmap pmap, key_type f, value_type id) {
      put(pmap.m_fid_map, f, int(id));
    }

    friend reference get(Patch_id_pmap pmap, key_type f) {
      return get(pmap.m_fid_map, f);
    }

    Face_id_map m_fid_map;
  };

  typedef VSA_wrapper::Indexed_triangle Indexed_triangle;
  typedef std::map<Scene_surface_mesh_item *, VSA_wrapper *> SM_wrapper_map;
  typedef std::pair<Scene_surface_mesh_item *, VSA_wrapper *> SM_wrapper_pair;

public:
  Polyhedron_demo_surface_mesh_approximation_plugin() {
    std::srand(time(0));
  }

  void init(QMainWindow *main_window,
    Scene_interface *scene_interface,
    Messages_interface *message_interface) {
    mw = main_window;
    scene = scene_interface;
    mi = message_interface;

    actionSurfaceMeshApproximation = new QAction("Surface Mesh Approximation", mw);
    actionSurfaceMeshApproximation->setProperty("subMenuName", "Triangulated Surface Mesh Approximation");
    connect(actionSurfaceMeshApproximation, SIGNAL(triggered()), this, SLOT(on_actionSurfaceMeshApproximation_triggered()));

    dock_widget = new QDockWidget("Mesh approximation parameters", mw);
    dock_widget->setVisible(true);
    ui_widget.setupUi(dock_widget);
    mw->addDockWidget(Qt::LeftDockWidgetArea, dock_widget);

    // connect ui actions
    connect(ui_widget.comboMetric, SIGNAL(currentIndexChanged(int)), this, SLOT(on_comboMetric_currentIndexChanged(const int)));
    connect(ui_widget.buttonSeeding, SIGNAL(clicked()), this, SLOT(on_buttonSeeding_clicked()));
    connect(ui_widget.buttonAdd, SIGNAL(clicked()), this, SLOT(on_buttonAdd_clicked()));
    connect(ui_widget.buttonFit, SIGNAL(clicked()), this, SLOT(on_buttonFit_clicked()));
    connect(ui_widget.buttonMeshing, SIGNAL(clicked()), this, SLOT(on_buttonMeshing_clicked()));
    connect(ui_widget.buttonTeleport, SIGNAL(clicked()), this, SLOT(on_buttonTeleport_clicked()));
    connect(ui_widget.buttonSplit, SIGNAL(clicked()), this, SLOT(on_buttonSplit_clicked()));
    // scene item delete action
    if (Scene *scene = dynamic_cast<Scene *>(scene_interface))
      connect(scene, SIGNAL(itemAboutToBeDestroyed(CGAL::Three::Scene_item *)),
        this, SLOT(itemAboutToBeDestroyed(CGAL::Three::Scene_item *)));
  }

  void closure() {
    // dock_widget->hide();
  }

  QList<QAction *> actions() const {
    return QList<QAction*>() << actionSurfaceMeshApproximation;
  }

  bool applicable(QAction *) const {
    return qobject_cast<Scene_surface_mesh_item *>(scene->item(scene->mainSelectionIndex()));
  }

public Q_SLOTS:
  void on_actionSurfaceMeshApproximation_triggered();
  void on_buttonSeeding_clicked();
  void on_buttonFit_clicked();
  void on_buttonAdd_clicked();
  void on_buttonTeleport_clicked();
  void on_buttonSplit_clicked();
  void on_buttonMeshing_clicked();
  void on_comboMetric_currentIndexChanged(const int);
  void itemAboutToBeDestroyed(CGAL::Three::Scene_item *);

private:
  QAction *actionSurfaceMeshApproximation;
  Ui::Surface_mesh_approximation ui_widget;
  QDockWidget *dock_widget;

  QMainWindow *mw;
  Scene_interface *scene;
  Messages_interface *mi;

  SM_wrapper_map m_sm_wrapper_map;
}; // end Polyhedron_demo_surface_mesh_approximation_plugin

void Polyhedron_demo_surface_mesh_approximation_plugin::on_actionSurfaceMeshApproximation_triggered()
{ dock_widget->show(); }

void Polyhedron_demo_surface_mesh_approximation_plugin::on_buttonSeeding_clicked() {
  Scene_surface_mesh_item *sm_item = qobject_cast<Scene_surface_mesh_item *>(
    scene->item(scene->mainSelectionIndex()));
  if (!sm_item) {
    mi->information(QString("No surface mesh item selected."));
    return;
  }

  QApplication::setOverrideCursor(Qt::WaitCursor);

  SMesh *pmesh = sm_item->face_graph();
  SM_wrapper_map::iterator search = m_sm_wrapper_map.find(sm_item);
  if (search == m_sm_wrapper_map.end())
    search = m_sm_wrapper_map.insert(SM_wrapper_pair(sm_item, new VSA_wrapper(*pmesh))).first;
  VSA_wrapper &approx = *search->second;
  approx.set_metric(static_cast<VSA_wrapper::Metric>(
    ui_widget.comboMetric->currentIndex()));

  QTime time;
  time.start();
  approx.initialize_seeds(CGAL::parameters::seeding_method(
    static_cast<VSA::Seeding_method>(ui_widget.comboMethod->currentIndex()))
    .max_number_of_proxies((ui_widget.cb_nb_proxies->isChecked() ?
      boost::optional<std::size_t>(ui_widget.nb_proxies->value()) : boost::none))
    .min_error_drop((ui_widget.cb_error_drop->isChecked() ?
      boost::optional<EPICK::FT>(ui_widget.error_drop->value()) : boost::none))
    .number_of_relaxations(ui_widget.nb_relaxations->value()));
  approx.run(ui_widget.nb_iterations->value());

  Patch_id_pmap pidmap(get(CGAL::face_patch_id_t<int>(), *sm_item->face_graph()));
  approx.output(CGAL::parameters::face_proxy_map(pidmap));

  mi->information(QString("Done, #proxies = %1. (%2 ms)").arg(
    approx.number_of_proxies()).arg(time.elapsed()));

  sm_item->color_vector() = approx.proxy_colors();
  sm_item->setItemIsMulticolor(true);
  sm_item->invalidateOpenGLBuffers();
  scene->itemChanged(scene->item_id(sm_item));

  QApplication::restoreOverrideCursor();
}

void Polyhedron_demo_surface_mesh_approximation_plugin::on_buttonFit_clicked() {
  Scene_surface_mesh_item *sm_item = qobject_cast<Scene_surface_mesh_item *>(
    scene->item(scene->mainSelectionIndex()));
  if (!sm_item) {
    mi->information(QString("No surface mesh item selected."));
    return;
  }
  SM_wrapper_map::iterator search = m_sm_wrapper_map.find(sm_item);
  if (search == m_sm_wrapper_map.end() || !search->second->initialized()) {
    mi->information(QString("Please initialize seeds first."));
    return;
  }
  SMesh *pmesh = search->first->face_graph();
  VSA_wrapper &approx = *search->second;

  QApplication::setOverrideCursor(Qt::WaitCursor);

  approx.run(1);
  Patch_id_pmap pidmap(get(CGAL::face_patch_id_t<int>(), *pmesh));
  approx.output(CGAL::parameters::face_proxy_map(pidmap));

  mi->information(QString("Fit one iteration, #proxies = %1.").arg(approx.number_of_proxies()));

  sm_item->color_vector() = approx.proxy_colors();
  sm_item->setItemIsMulticolor(true);
  sm_item->invalidateOpenGLBuffers();
  scene->itemChanged(scene->item_id(sm_item));

  QApplication::restoreOverrideCursor();
}

void Polyhedron_demo_surface_mesh_approximation_plugin::on_buttonAdd_clicked() {
  Scene_surface_mesh_item *sm_item = qobject_cast<Scene_surface_mesh_item *>(
    scene->item(scene->mainSelectionIndex()));
  if (!sm_item) {
    mi->information(QString("No surface mesh item selected."));
    return;
  }
  SM_wrapper_map::iterator search = m_sm_wrapper_map.find(sm_item);
  if (search == m_sm_wrapper_map.end() || !search->second->initialized()) {
    mi->information(QString("Please initialize seeds first."));
    return;
  }
  SMesh *pmesh = search->first->face_graph();
  VSA_wrapper &approx = *search->second;

  QApplication::setOverrideCursor(Qt::WaitCursor);

  if (approx.add_one_proxy() == 0) {
    mi->information(QString("No proxy added, #proxies = %1.").arg(approx.number_of_proxies()));
    QApplication::restoreOverrideCursor();
    return;
  }
  mi->information(QString("One proxy added, #proxies = %1.").arg(approx.number_of_proxies()));

  Patch_id_pmap pidmap(get(CGAL::face_patch_id_t<int>(), *pmesh));
  approx.output(CGAL::parameters::face_proxy_map(pidmap));

  sm_item->color_vector() = approx.proxy_colors();
  sm_item->setItemIsMulticolor(true);
  sm_item->invalidateOpenGLBuffers();
  scene->itemChanged(scene->item_id(sm_item));

  QApplication::restoreOverrideCursor();
}

void Polyhedron_demo_surface_mesh_approximation_plugin::on_buttonTeleport_clicked() {
  Scene_surface_mesh_item *sm_item = qobject_cast<Scene_surface_mesh_item *>(
    scene->item(scene->mainSelectionIndex()));
  if (!sm_item) {
    mi->information(QString("No surface mesh item selected."));
    return;
  }
  SM_wrapper_map::iterator search = m_sm_wrapper_map.find(sm_item);
  if (search == m_sm_wrapper_map.end() || !search->second->initialized()) {
    mi->information(QString("Please initialize seeds first."));
    return;
  }
  SMesh *pmesh = search->first->face_graph();
  VSA_wrapper &approx = *search->second;

  QApplication::setOverrideCursor(Qt::WaitCursor);

  if (approx.teleport_one_proxy() == 0) {
    mi->information(QString("No proxy teleported, #proxies = %1.").arg(approx.number_of_proxies()));
    QApplication::restoreOverrideCursor();
    return;
  }
  Patch_id_pmap pidmap(get(CGAL::face_patch_id_t<int>(), *pmesh));
  approx.output(CGAL::parameters::face_proxy_map(pidmap));
  mi->information(QString("One proxy teleported, #proxies = %1.").arg(approx.number_of_proxies()));

  sm_item->color_vector() = approx.proxy_colors();
  sm_item->setItemIsMulticolor(true);
  sm_item->invalidateOpenGLBuffers();
  scene->itemChanged(scene->item_id(sm_item));

  QApplication::restoreOverrideCursor();
}

void Polyhedron_demo_surface_mesh_approximation_plugin::on_buttonSplit_clicked() {
  Scene_surface_mesh_item *sm_item = qobject_cast<Scene_surface_mesh_item *>(
    scene->item(scene->mainSelectionIndex()));
  if (!sm_item) {
    mi->information(QString("No surface mesh item selected."));
    return;
  }
  SM_wrapper_map::iterator search = m_sm_wrapper_map.find(sm_item);
  if (search == m_sm_wrapper_map.end() || !search->second->initialized()) {
    mi->information(QString("Please initialize seeds first."));
    return;
  }
  SMesh *pmesh = search->first->face_graph();
  VSA_wrapper &approx = *search->second;

  QApplication::setOverrideCursor(Qt::WaitCursor);

  if (!approx.split(ui_widget.split_proxy_idx->value(),
    ui_widget.split_nb_sections->value(),
    ui_widget.split_nb_relaxations->value())) {
    mi->information(QString("No proxy splitted, #proxies = %1.").arg(approx.number_of_proxies()));
    QApplication::restoreOverrideCursor();
    return;
  }
  mi->information(QString("One proxy splitted, #proxies = %1.").arg(approx.number_of_proxies()));

  Patch_id_pmap pidmap(get(CGAL::face_patch_id_t<int>(), *pmesh));
  approx.output(CGAL::parameters::face_proxy_map(pidmap));

  sm_item->color_vector() = approx.proxy_colors();
  sm_item->setItemIsMulticolor(true);
  sm_item->invalidateOpenGLBuffers();
  scene->itemChanged(scene->item_id(sm_item));

  QApplication::restoreOverrideCursor();
}

void Polyhedron_demo_surface_mesh_approximation_plugin::on_buttonMeshing_clicked() {
  Scene_surface_mesh_item *sm_item = qobject_cast<Scene_surface_mesh_item *>(
    scene->item(scene->mainSelectionIndex()));
  if (!sm_item) {
    mi->information(QString("No surface mesh item selected."));
    return;
  }
  SM_wrapper_map::iterator search = m_sm_wrapper_map.find(sm_item);
  if (search == m_sm_wrapper_map.end() || !search->second->initialized()) {
    mi->information(QString("Please initialize seeds first."));
    return;
  }
  SMesh *pmesh = search->first->face_graph();
  VSA_wrapper &approx = *search->second;

  QApplication::setOverrideCursor(Qt::WaitCursor);

  approx.extract_mesh(CGAL::parameters::subdivision_ratio(ui_widget.chord_error->value())
    .relative_to_chord(ui_widget.is_relative_to_chord->isChecked())
    .with_dihedral_angle(ui_widget.with_dihedral_angle->isChecked())
    .optimize_anchor_location(ui_widget.if_optimize_anchor_location->isChecked())
    .pca_plane(ui_widget.pca_plane->isChecked()));

  std::vector<Point_3> anchor_points;
  std::vector<Indexed_triangle> indexed_triangles;
  std::vector<vertex_descriptor> anchor_vertices;
  std::vector<std::vector<std::size_t> > patch_polygons;

  Patch_id_pmap pidmap(get(CGAL::face_patch_id_t<int>(), *sm_item->face_graph()));
  approx.output(CGAL::parameters::face_proxy_map(pidmap)
    .anchors(std::back_inserter(anchor_points))
    .triangles(std::back_inserter(indexed_triangles)));
  approx.anchor_vertices(std::back_inserter(anchor_vertices));
  approx.indexed_boundary_polygons(std::back_inserter(patch_polygons));

  // create a group
  Scene_group_item *group = new Scene_group_item(tr("Approximation of %1").arg(sm_item->name()));
  scene->addItem(group);

  // add triangle soup item
  Scene_polygon_soup_item *psoup_item = new Scene_polygon_soup_item();
  std::vector<std::vector<std::size_t> > polygons;
  BOOST_FOREACH(const Indexed_triangle &t, indexed_triangles) {
    std::vector<std::size_t> polygon;
    polygon.push_back(t[0]);
    polygon.push_back(t[1]);
    polygon.push_back(t[2]);
    polygons.push_back(polygon);
  }
  psoup_item->load(anchor_points, polygons);
  psoup_item->setName(tr("Triangle soup").arg(sm_item->name()));
  psoup_item->setColor(Qt::lightGray);
  psoup_item->setRenderingMode(FlatPlusEdges);
  scene->addItem(psoup_item);
  scene->changeGroup(psoup_item, group);

  // add patch border item
  Scene_polylines_item *borders_item = new Scene_polylines_item();
  BOOST_FOREACH(const std::vector<std::size_t> &border, patch_polygons) {
    std::vector<Point_3> polyline;
    for (std::size_t i = 0; i <= border.size(); ++i)
      polyline.push_back(anchor_points[border[i % border.size()]]);
    borders_item->polylines.push_back(polyline);
  }
  borders_item->setName(tr("Patch polygons").arg(sm_item->name()));
  borders_item->setColor(Qt::red);
  borders_item->invalidateOpenGLBuffers();
  scene->addItem(borders_item);
  scene->changeGroup(borders_item, group);

  // add anchor vertex and anchor point correspondence item
  Scene_polylines_item *corres_item = new Scene_polylines_item();
  for (std::size_t i = 0; i < anchor_vertices.size(); ++i) {
    std::vector<Point_3> polyline;
    polyline.push_back(pmesh->point(anchor_vertices[i]));
    polyline.push_back(anchor_points[i]);
    corres_item->polylines.push_back(polyline);
  }
  corres_item->setName(tr("Anchors").arg(sm_item->name()));
  corres_item->setColor(Qt::blue);
  corres_item->invalidateOpenGLBuffers();
  scene->addItem(corres_item);
  scene->changeGroup(corres_item, group);

  // add patch convex hull item
  std::vector<std::vector<EPICK::Triangle_3> > patch_triangles(approx.number_of_proxies());
  BOOST_FOREACH(face_descriptor f, faces(*pmesh)) {
    halfedge_descriptor h = halfedge(f, *pmesh);
    patch_triangles[get(pidmap, f)].push_back(EPICK::Triangle_3(
      pmesh->point(source(h, *pmesh)),
      pmesh->point(target(h, *pmesh)),
      pmesh->point(target(next(h, *pmesh), *pmesh))));
  }
  std::vector<EPICK::Plane_3> patch_planes;
  BOOST_FOREACH(const std::vector<EPICK::Triangle_3> &tris, patch_triangles) {
    EPICK::Plane_3 fit_plane;
    CGAL::linear_least_squares_fitting_3(
      tris.begin(), tris.end(), fit_plane, CGAL::Dimension_tag<2>());
    patch_planes.push_back(fit_plane);
  }
  std::vector<std::vector<Point_3> > patch_points(approx.number_of_proxies());
  BOOST_FOREACH(vertex_descriptor v, vertices(*pmesh)) {
    BOOST_FOREACH(halfedge_descriptor h, CGAL::halfedges_around_target(v, *pmesh)) {
      if (!CGAL::is_border(h, *pmesh)) {
        const std::size_t fidx = get(pidmap, face(h, *pmesh));
        patch_points[fidx].push_back(pmesh->point(v));
      }
    }
  }
  std::vector<Point_3> cvx_hull_points;
  std::vector<std::vector<std::size_t> > cvx_hulls;
  for (std::size_t i = 0; i < approx.number_of_proxies(); ++i) {
    const std::vector<Point_3> &pts = patch_points[i];
    const EPICK::Plane_3 plane = patch_planes[i];
    const Point_3 origin = plane.projection(pts.front());

    EPICK::Vector_3 base1 = plane.base1();
    EPICK::Vector_3 base2 = plane.base2();
    base1 = base1 / std::sqrt(base1.squared_length());
    base2 = base2 / std::sqrt(base2.squared_length());

    EPICK::Line_3 base_linex(origin, base1);
    EPICK::Line_3 base_liney(origin, base2);

    std::vector<EPICK::Point_2> pts_2d;
    BOOST_FOREACH(const Point_3 &p, pts) {
      const Point_3 point = plane.projection(p);
      EPICK::Vector_3 vecx(origin, base_linex.projection(point));
      EPICK::Vector_3 vecy(origin, base_liney.projection(point));
      double x = std::sqrt(vecx.squared_length());
      double y = std::sqrt(vecy.squared_length());
      x = vecx * base1 < 0 ? -x : x;
      y = vecy * base2 < 0 ? -y : y;
      pts_2d.push_back(EPICK::Point_2(x, y));
    }

    std::vector<EPICK::Point_2> cvx_hull_2d;
    CGAL::convex_hull_2(pts_2d.begin(), pts_2d.end(), std::back_inserter(cvx_hull_2d));

    cvx_hulls.push_back(std::vector<std::size_t>());
    BOOST_FOREACH(const EPICK::Point_2 &p, cvx_hull_2d) {
      cvx_hulls.back().push_back(cvx_hull_points.size());
      cvx_hull_points.push_back(origin + p.x() * base1 + p.y() * base2);
    }
  }
  std::vector<CGAL::Color> fcolors;
  BOOST_FOREACH(const QColor &c, approx.proxy_colors())
    fcolors.push_back(CGAL::Color(c.red(), c.green(), c.blue()));
  Scene_polygon_soup_item *cvx_hull_item = new Scene_polygon_soup_item();
  cvx_hull_item->load(cvx_hull_points, cvx_hulls, fcolors, std::vector<CGAL::Color>());
  cvx_hull_item->setName(tr("Patch planes").arg(sm_item->name()));
  cvx_hull_item->setColor(Qt::yellow);
  cvx_hull_item->setRenderingMode(FlatPlusEdges);
  scene->addItem(cvx_hull_item);
  scene->changeGroup(cvx_hull_item, group);

  mi->information(QString("Meshing done."));

  QApplication::restoreOverrideCursor();
}

void Polyhedron_demo_surface_mesh_approximation_plugin::on_comboMetric_currentIndexChanged(const int) {
  Scene_surface_mesh_item *sm_item = qobject_cast<Scene_surface_mesh_item *>(
    scene->item(scene->mainSelectionIndex()));
  if (!sm_item) {
    mi->information(QString("No surface mesh item selected."));
    return;
  }
  SM_wrapper_map::iterator search = m_sm_wrapper_map.find(sm_item);
  if (search == m_sm_wrapper_map.end())
    return;
  search->second->initialized() = false;

  sm_item->color_vector().clear();
  sm_item->setItemIsMulticolor(false);
  sm_item->invalidateOpenGLBuffers();
  scene->itemChanged(scene->item_id(sm_item));
}

void Polyhedron_demo_surface_mesh_approximation_plugin::itemAboutToBeDestroyed(CGAL::Three::Scene_item *scene_item)
{
  if (Scene_surface_mesh_item *sm_item = qobject_cast<Scene_surface_mesh_item *>(scene_item)) {
    SM_wrapper_map::iterator search = m_sm_wrapper_map.find(sm_item);
    if (search != m_sm_wrapper_map.end()) {
      delete search->second;
      m_sm_wrapper_map.erase(sm_item);
    }
  }
}

#include "Surface_mesh_approximation_plugin.moc"
