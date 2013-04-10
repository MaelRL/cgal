#include "config.h"

#include <CGAL_demo/Plugin_helper.h>
#include <CGAL_demo/Plugin_interface.h>
#include <CGAL_demo/Messages_interface.h>
#include "ui_Aniso_meshing_dialog.h"
#include "ui_Choose_star_dialog.h"
#include "ui_Choose_pickvalid_point_dialog.h"

#include <QObject>
#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include <QMessageBox>
#include <QInputDialog>
#include <QFileDialog>
#include <QTimer>
#include <QVariant>

#include "Scene_polyhedron_item.h"
#include "Scene_constrained_surface_implicit_item.h"

#include "Polyhedron_type_fwd.h"
#include "Implicit_surface_type_fwd.h"

#include "Scene_starset3_item.h"
#include "Anisotropic_meshing_thread.h"
#include "Anisotropic_mesh_function.h"
#include "anisotropic_meshing_options.h"

#include <iostream>
#include <fstream>
#include <math.h>

// Constants
const QColor default_mesh_color(45,169,70);


// declare the CGAL functions

  //polyhedron functions
Criteria* build_param_and_metric(const CGAL::Anisotropic_mesh_3::Constrain_surface_3_polyhedral<Kernel, Polyhedron>* const p_domain,
                                 Anisotropic_mesh_parameters & param,
                                 CGAL::Anisotropic_mesh_3::Metric_field<Kernel>* & mf,
                                 const double epsilon,
                                 const double approximation,
                                 const double radius_edge_ratio,
                                 const double sliverity,
                                 const double circumradius,
                                 const double distortion,
                                 const double beta,
                                 const double delta,
                                 const std::size_t max_times_to_try_in_picking_region,
                                 const int dim,
                                 const Metric_options& metric,
                                 const double en_factor);

Anisotropic_meshing_thread* cgal_code_anisotropic_mesh_3(const Polyhedron*,
                                 const double epsilon,
                                 const double approximation,
                                 const double radius_edge_ratio,
                                 const double sliverity,
                                 const double circumradius,
                                 const double distortion,
                                 const double beta,
                                 const double delta,
                                 const std::size_t max_times_to_try_in_picking_region,
                                 const int dim,
                                 const int nb_initial_points,
                                 const Metric_options& metric,
                                 const bool pick_valid_causes_stop,
                                 const bool pick_valid_use_cube_probing,
                                 const int pick_valid_max_failures,
                                 const double en_factor);

  //implicit functions
Criteria* build_param_and_metric(const Implicit_surface* p_domain,
                                 Anisotropic_mesh_parameters & param,
                                 CGAL::Anisotropic_mesh_3::Metric_field<Kernel>* & mf,
                                 const double epsilon,
                                 const double approximation,
                                 const double radius_edge_ratio,
                                 const double sliverity,
                                 const double circumradius,
                                 const double distortion,
                                 const double beta,
                                 const double delta,
                                 const std::size_t max_times_to_try_in_picking_region,
                                 const int dim,
                                 const Metric_options& metric,
                                 const double en_factor);


Anisotropic_meshing_thread* cgal_code_anisotropic_mesh_3(const Implicit_surface*,
                                 const double epsilon,
                                 const double approximation,
                                 const double radius_edge_ratio,
                                 const double sliverity,
                                 const double circumradius,
                                 const double distortion,
                                 const double beta,
                                 const double delta,
                                 const std::size_t max_times_to_try_in_picking_region,
                                 const int dim,
                                 const int nb_initial_points,
                                 const Metric_options& metric,
                                 const bool pick_valid_causes_stop,
                                 const bool pick_valid_use_cube_probing,
                                 const int pick_valid_max_failures,
                                 const double en_factor);

double get_approximate(double d, int precision, int& decimals);


class Anisotropic_mesh_3_plugin : 
  public QObject,
  protected Plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(Plugin_interface)
  
  typedef Plugin_helper Base;
public:
  Anisotropic_mesh_3_plugin();
  
  using Base::init;
  virtual void init(QMainWindow* mainWindow, 
                    Scene_interface* scene_interface,
                    Messages_interface* msg_interface)
  {
    this->scene = scene_interface;
    this->mw = mainWindow;
    this->msg = msg_interface;

    actionAnisotropicMeshing = this->getActionFromMainWindow(mw, "actionAnisotropicMeshing");
    if(actionAnisotropicMeshing)
      connect(actionAnisotropicMeshing, SIGNAL(triggered()), this, SLOT(anisotropic_mesh_3()));

    actionResumeMeshing = this->getActionFromMainWindow(mw, "actionResumeMeshing");
    if(actionResumeMeshing)
      connect(actionResumeMeshing, SIGNAL(triggered()), this, SLOT(resume_aniso_mesh_3()));

    actionDraw_surface_star_set = this->getActionFromMainWindow(mw, "actionDraw_surface_star_set");
    if(actionDraw_surface_star_set)
      connect(actionDraw_surface_star_set, SIGNAL(toggled(bool)), this, SLOT(view_surface_star_set(bool)));

    actionDraw_cell = this->getActionFromMainWindow(mw, "actionDraw_cell");
    if(actionDraw_cell)
      connect(actionDraw_cell, SIGNAL(toggled(bool)), this, SLOT(view_cell(bool)));

    actionDraw_dual_edges = this->getActionFromMainWindow(mw, "actionDraw_dual_edges");
    if(actionDraw_dual_edges)
     connect(actionDraw_dual_edges, SIGNAL(toggled(bool)), this, SLOT(view_dual_edges(bool)));

    actionDraw_poles = this->getActionFromMainWindow(mw, "actionDraw_poles");
    if(actionDraw_poles)
     connect(actionDraw_poles, SIGNAL(toggled(bool)), this, SLOT(view_poles(bool)));

    actionDraw_initial_points = this->getActionFromMainWindow(mw, "actionDraw_initial_points");
    if(actionDraw_initial_points)
     connect(actionDraw_initial_points, SIGNAL(toggled(bool)), this, SLOT(view_initial_points(bool)));

    actionDraw_surface_delaunay_balls = this->getActionFromMainWindow(mw, "actionDraw_surface_delaunay_balls");
    if(actionDraw_surface_delaunay_balls)
      connect(actionDraw_surface_delaunay_balls, SIGNAL(toggled(bool)), this, SLOT(view_surface_delaunay_balls(bool)));

    actionDraw_one_star = this->getActionFromMainWindow(mw, "actionDraw_one_star");
    if(actionDraw_one_star)
      connect(actionDraw_one_star, SIGNAL(triggered()), this, SLOT(view_one_star()));

    actionDraw_one_pickvalid_point = this->getActionFromMainWindow(mw, "actionDraw_one_pickvalid_point");
    if(actionDraw_one_pickvalid_point)
      connect(actionDraw_one_pickvalid_point, SIGNAL(triggered()), this, SLOT(view_one_pickvalid_point()));

    actionDraw_inconsistent_facets = this->getActionFromMainWindow(mw, "actionDraw_inconsistent_facets");
    if(actionDraw_inconsistent_facets)
      connect(actionDraw_inconsistent_facets, SIGNAL(toggled(bool)), this, SLOT(view_inconsistent_facets(bool)));

    actionDraw_metric_field = this->getActionFromMainWindow(mw, "actionDraw_metric_field");
    if(actionDraw_metric_field)
      connect(actionDraw_metric_field, SIGNAL(toggled(bool)), this, SLOT(view_metric_field(bool)));

    actionDraw_mesh_3 = this->getActionFromMainWindow(mw, "actionDraw_mesh_3");
    if(actionDraw_mesh_3)
      connect(actionDraw_mesh_3, SIGNAL(toggled(bool)), this, SLOT(view_mesh_3(bool)));

    actionDraw_distortion = this->getActionFromMainWindow(mw, "actionDraw_distortion");
    if(actionDraw_distortion)
      connect(actionDraw_distortion, SIGNAL(toggled(bool)), this, SLOT(view_distortion(bool)));
  }

  virtual QList<QAction*> actions() const
  {
    return QList<QAction*>() << actionAnisotropicMeshing << actionResumeMeshing;
  }
  
  bool applicable() const {
    return true; //checked within aniso functions anyway

//    const Scene_interface::Item_id index = scene->mainSelectionIndex();

//    // -----------------------------------
//    // Check if selected item is meshable
//    // -----------------------------------
//    Scene_polyhedron_item* poly_item =
//      qobject_cast<Scene_polyhedron_item*>(scene->item(index));

//    Scene_constrained_surface_implicit_item* function_item =
//      qobject_cast<Scene_constrained_surface_implicit_item*>(scene->item(index));

//    return poly_item || function_item;
  }

public slots:
  //meshing
  void anisotropic_mesh_3();
  void resume_aniso_mesh_3();
  //view
  void view_surface_star_set(bool);
  void view_cell  (bool);
  void view_dual_edges(bool);
  void view_poles(bool);
  void view_initial_points(bool);
  void view_surface_delaunay_balls(bool);
  void view_one_star();
  void view_all_stars();
  void view_one_pickvalid_point();
  void view_all_pickvalid_point();
  void view_inconsistent_facets(bool);
  void view_metric_field(bool);
  void view_mesh_3(bool);
  void view_distortion(bool);
  //others
  void meshing_done_resume(Anisotropic_meshing_thread* t);
  void meshing_done_nresume(Anisotropic_meshing_thread* t);
  void meshing_done(Anisotropic_meshing_thread* t, const bool resumed);
  void status_report(QString str);
  
private:
  void launch_thread(Anisotropic_meshing_thread* mesh_thread,
                const bool resumed = false);
  void treat_result(Scene_item& source_item,
               Scene_starset3_item& result_item,
               const bool replacing = false) const;

private:
  QAction* actionAnisotropicMeshing;
  QAction* actionResumeMeshing;
  QAction* actionDraw_surface_star_set;
  QAction* actionDraw_cell;
  QAction* actionDraw_dual_edges;
  QAction* actionDraw_poles;
  QAction* actionDraw_initial_points;
  QAction* actionDraw_surface_delaunay_balls;
  QAction* actionDraw_one_star;
  QAction* actionDraw_one_pickvalid_point;
  QAction* actionDraw_inconsistent_facets;
  QAction* actionDraw_metric_field;
  QAction* actionDraw_mesh_3;
  QAction* actionDraw_distortion;

  Messages_interface* msg;
  QMessageBox* message_box_;
  Scene_item* source_item_;
  
}; // end class Anisotropic_mesh_3_plugin

Anisotropic_mesh_3_plugin::
Anisotropic_mesh_3_plugin()
  : actionAnisotropicMeshing(NULL)
  , actionResumeMeshing(NULL)
  , actionDraw_surface_star_set(NULL)
  , actionDraw_cell(NULL)
  , actionDraw_dual_edges(NULL)
  , actionDraw_poles(NULL)
  , actionDraw_initial_points(NULL)
  , actionDraw_surface_delaunay_balls(NULL)
  , actionDraw_one_star(NULL)
  , actionDraw_one_pickvalid_point(NULL)
  , actionDraw_inconsistent_facets(NULL)
  , actionDraw_distortion(NULL)
  , msg(NULL)
  , message_box_(NULL)
  , source_item_(NULL)
{
}

void Anisotropic_mesh_3_plugin::anisotropic_mesh_3()
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();

  // -----------------------------------
  // Check if selected item is meshable
  // -----------------------------------
  Scene_polyhedron_item* poly_item = 
    qobject_cast<Scene_polyhedron_item*>(scene->item(index));

  Scene_constrained_surface_implicit_item* function_item =
    qobject_cast<Scene_constrained_surface_implicit_item*>(scene->item(index));

  // Get item
  Domain_type t;
  Scene_item* item = NULL;
  if( NULL != poly_item )          
  { 
    item = poly_item; 
    t = POLYHEDRAL_SURFACE; 
  }
  else if( NULL != function_item ) 
  { 
    item = function_item; 
    t = IMPLICIT_SURFACE; 
  }
  else if ( NULL == item )
  {
    QMessageBox::warning(mw,tr(""),tr("Selected object can't be meshed")); 
    return;
  }
  
  // -----------------------------------
  // Create Mesh dialog
  // -----------------------------------
  QDialog dialog(mw);
  Ui::Aniso_meshing_dialog ui;
  ui.setupUi(&dialog);
  connect(ui.buttonBox, SIGNAL(accepted()),
          &dialog, SLOT(accept()));
  connect(ui.buttonBox, SIGNAL(rejected()),
          &dialog, SLOT(reject()));

  // Connect checkboxes to spinboxes
  connect(ui.noRadiusEdgeRatio,   SIGNAL(toggled(bool)),
          ui.radius_edge_ratio,   SLOT(setEnabled(bool)));
  connect(ui.noSliverity,         SIGNAL(toggled(bool)),
          ui.sliverity,           SLOT(setEnabled(bool)));
  connect(ui.noCircumradius,      SIGNAL(toggled(bool)),
          ui.circumradius,        SLOT(setEnabled(bool)));
  connect(ui.noDistortion,        SIGNAL(toggled(bool)),
          ui.distortion,          SLOT(setEnabled(bool)));
  connect(ui.noApproximation,     SIGNAL(toggled(bool)),
          ui.approximation,       SLOT(setEnabled(bool)));

   // Set default parameters
  Scene_interface::Bbox bbox = item->bbox();
  ui.objectName->setText(item->name());
  ui.objectNameSize->setText(tr("Object bbox size (w,h,d):  <b>%1</b>,  <b>%2</b>,  <b>%3</b>")
                             .arg(bbox.width(),0,'g',3)
                             .arg(bbox.height(),0,'g',3)
                             .arg(bbox.depth(),0,'g',3) );  
  double max_size = (std::max)((std::max)(bbox.width(),bbox.height()),bbox.depth());
  ui.radius_edge_ratio->setValue(3.0);
  ui.sliverity->setValue(0.3);
  ui.circumradius->setValue(1.0);
  ui.distortion->setValue(1.8);
  ui.approximation->setValue(0.01*max_size);
  ui.epsilon->setValue(0.1);
  ui.beta->setValue(2.5);
  ui.delta->setValue(0.3);
  
  ui.resume_memory->setEnabled(false);

  ui.comboBox_metric->addItem("Euclidean", EUCLIDEAN);
  ui.comboBox_metric->addItem("Hyperbolic shock", HYPERBOLIC_SHOCK);
  if(t == POLYHEDRAL_SURFACE)
  {
    ui.comboBox_metric->addItem("Polyhedron curvature", POLYHEDRON_CURVATURE);
    ui.comboBox_metric->setCurrentIndex(2);//polyhedron curvature
  }
  else //IMPLICIT_SURFACE
  {
    ui.comboBox_metric->addItem("Implicit curvature", IMPLICIT_CURVATURE);
    ui.comboBox_metric->addItem("Torus (1/r, 1/R)", TORUS_NAIVE);
    ui.comboBox_metric->setCurrentIndex(2);//implicit curvature
  }
 
  ui.maxTries->setDecimals(0);
  ui.maxTries->setSingleStep(1);
  ui.maxTries->setValue(60); // default value
  ui.dimension->setCurrentIndex(0);
  ui.nbInitialPoints->setValue(10);
  ui.en_factor->setValue(1.);

  // -----------------------------------
  // Get values
  // -----------------------------------
  int i = dialog.exec();
  if( i == QDialog::Rejected ) { return; }
  
  // 0 means parameter is not considered
  const double radius_edge_ratio = !ui.noRadiusEdgeRatio->isChecked() ? 0 : ui.radius_edge_ratio->value(); 
  const double sliverity = !ui.noSliverity->isChecked() ? 0 : ui.sliverity->value(); 
  const double circumradius = !ui.noCircumradius->isChecked() ? 0 : ui.circumradius->value(); 
  const double distortion = !ui.noDistortion->isChecked() ? 0 : ui.distortion->value();
  const double approximation = !ui.noApproximation->isChecked() ? 0. : ui.approximation->value();
  const double epsilon = ui.epsilon->value();
  const double beta = ui.beta->value();
  const double delta = ui.delta->value();
  const int max_times_to_try_in_picking_region = ui.maxTries->value();
  const bool pick_valid_causes_stop = ui.pick_valid_causes_stop->isChecked();
  const bool pick_valid_use_cube_probing = ui.cube_probing->isChecked();
  const int pick_valid_max_failures = ui.maxPickFailures->value();
  int dim = -1;
  if(ui.dimension->currentText().compare(QString("Surface")) == 0)
    dim = 2;
  else if(ui.dimension->currentText().compare(QString("Volume")) == 0)
    dim = 3;
  const int nb_initial_points = ui.nbInitialPoints->value();
  int metric_index = ui.comboBox_metric->currentIndex();
  Metric_options metric = (Metric_options)ui.comboBox_metric->itemData(metric_index).toInt();
  const double en_factor = ui.en_factor->value();

  // -----------------------------------
  // Dispatch mesh process
  // -----------------------------------
  QApplication::setOverrideCursor(Qt::WaitCursor);
  
  Anisotropic_meshing_thread* thread = NULL;
  
  // Polyhedron
  if ( NULL != poly_item )
  {
    Polyhedron* pPoly = poly_item->polyhedron();
    if( NULL == pPoly )
    {
      QMessageBox::critical(mw,tr(""),tr("ERROR: no data in selected item")); 
      return;
    }
    thread = cgal_code_anisotropic_mesh_3(pPoly, epsilon,
      approximation, radius_edge_ratio, sliverity, circumradius,
      distortion, beta, delta, max_times_to_try_in_picking_region,
      dim, nb_initial_points, metric, pick_valid_causes_stop,
      pick_valid_use_cube_probing, pick_valid_max_failures,
      en_factor);
  }
  //// Function
  else if( NULL != function_item )
  {
    Implicit_surface* pFunction = function_item->function();
    if( NULL == pFunction )
    {
      QMessageBox::critical(mw,tr(""),tr("ERROR: no data in selected item")); 
      return;
    }
    thread = cgal_code_anisotropic_mesh_3(pFunction, epsilon,
      approximation, radius_edge_ratio, sliverity, circumradius,
      distortion, beta, delta, max_times_to_try_in_picking_region,
      dim, nb_initial_points, metric, pick_valid_causes_stop,
      pick_valid_use_cube_probing, pick_valid_max_failures,
      en_factor);
  }

  if ( NULL == thread )
  {
    QMessageBox::critical(mw,tr(""),tr("ERROR: no thread created")); 
    return;
  }
  
  // Launch thread
  source_item_ = item;
  launch_thread(thread);
  QApplication::restoreOverrideCursor();
}

void Anisotropic_mesh_3_plugin::resume_aniso_mesh_3(){
  // ------------------------------------------
  // Check if selected item can be resumed from
  // ------------------------------------------
  Scene_item* item = NULL;
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  Scene_starset3_item* ssetitem = qobject_cast<Scene_starset3_item*>(scene->item(index));

  if( NULL == ssetitem )
  {
    QMessageBox::warning(mw,tr(""),tr("Nothing to resume from"));
    return;
  }
  else{
    item = ssetitem;
  }

  typedef CGAL::Anisotropic_mesh_3::Constrain_surface_3_polyhedral<Kernel, Polyhedron>  Constrain_surface_polyhedral;
  const Constrain_surface_polyhedral* poly_surf =
      dynamic_cast<const Constrain_surface_polyhedral*>(ssetitem->star_set().constrain_surface());

  const Implicit_surface* function_surf =
      dynamic_cast<const Implicit_surface*>(ssetitem->star_set().constrain_surface());

  // Get item
  if( NULL == poly_surf && NULL == function_surf)
  {
    QMessageBox::warning(mw,tr(""),tr("Unknown surface type"));
    return;
  }

  // ------------------
  // Create Mesh dialog
  // ------------------
  QDialog dialog(mw);
  Ui::Aniso_meshing_dialog ui;
  ui.setupUi(&dialog);
  connect(ui.buttonBox, SIGNAL(accepted()),
          &dialog, SLOT(accept()));
  connect(ui.buttonBox, SIGNAL(rejected()),
          &dialog, SLOT(reject()));

  // Connect checkboxes to spinboxes
  connect(ui.noRadiusEdgeRatio,   SIGNAL(toggled(bool)),
          ui.radius_edge_ratio,   SLOT(setEnabled(bool)));
  connect(ui.noSliverity,         SIGNAL(toggled(bool)),
          ui.sliverity,           SLOT(setEnabled(bool)));
  connect(ui.noCircumradius,      SIGNAL(toggled(bool)),
          ui.circumradius,        SLOT(setEnabled(bool)));
  connect(ui.noDistortion,        SIGNAL(toggled(bool)),
          ui.distortion,          SLOT(setEnabled(bool)));
  connect(ui.noApproximation,     SIGNAL(toggled(bool)),
          ui.approximation,       SLOT(setEnabled(bool)));

   // Set default parameters
  Scene_interface::Bbox bbox = item->bbox();
  ui.objectName->setText(item->name());
  ui.objectNameSize->setText(tr("Object bbox size (w,h,d):  <b>%1</b>,  <b>%2</b>,  <b>%3</b>")
                             .arg(bbox.width(),0,'g',3)
                             .arg(bbox.height(),0,'g',3)
                             .arg(bbox.depth(),0,'g',3) );
  double max_size = (std::max)((std::max)(bbox.width(),bbox.height()),bbox.depth());
  ui.radius_edge_ratio->setValue(3.0);
  ui.sliverity->setValue(0.3);
  ui.circumradius->setValue(1.0);
  ui.distortion->setValue(1.8);
  ui.approximation->setValue(0.01*max_size);
  ui.epsilon->setValue(0.1);
  ui.beta->setValue(2.5);
  ui.delta->setValue(0.3);
  ui.en_factor->setValue(1.);

  ui.comboBox_metric->setEnabled(false);
  //disable the metric field when restarting : not yet implemented (todo?)

  ui.maxTries->setDecimals(0);
  ui.maxTries->setSingleStep(1);
  ui.maxTries->setValue(60); // default value
  ui.dimension->setCurrentIndex(0);
  ui.nbInitialPoints->setEnabled(false);

  // -----------------------------------
  // Get values
  // -----------------------------------
  int i = dialog.exec();
  if( i == QDialog::Rejected ) { return; }

  // 0 means parameter is not considered
  const bool keep_same_parameters = ui.resume_memory->isChecked();

  const double radius_edge_ratio = !ui.noRadiusEdgeRatio->isChecked() ? 0 : ui.radius_edge_ratio->value();
  const double sliverity = !ui.noSliverity->isChecked() ? 0 : ui.sliverity->value();
  const double circumradius = !ui.noCircumradius->isChecked() ? 0 : ui.circumradius->value();
  const double distortion = !ui.noDistortion->isChecked() ? 0 : ui.distortion->value();
  const double approximation = !ui.noApproximation->isChecked() ? 0. : ui.approximation->value();
  const double epsilon = ui.epsilon->value();
  const double beta = ui.beta->value();
  const double delta = ui.delta->value();
  const int max_times_to_try_in_picking_region = ui.maxTries->value();
  const bool pick_valid_causes_stop = ui.pick_valid_causes_stop->isChecked();
  const bool pick_valid_use_cube_probing = ui.cube_probing->isChecked();
  const int pick_valid_max_failures = ui.maxPickFailures->value();
  int dim = -1;
  if(ui.dimension->currentText().compare(QString("Surface")) == 0)
    dim = 2;
  else if(ui.dimension->currentText().compare(QString("Volume")) == 0)
    dim = 3;
  int metric_index = ui.comboBox_metric->currentIndex();
  Metric_options metric = (Metric_options)ui.comboBox_metric->itemData(metric_index).toInt();
  const double en_factor = ui.en_factor->value();

  // -----------------------------------
  // Dispatch mesh process
  // -----------------------------------
  QApplication::setOverrideCursor(Qt::WaitCursor);

  Anisotropic_meshing_thread* thread = NULL;

  // POLYHEDRON -------------------------------------------------------------------------------

  if ( NULL != poly_surf )
  {
    Anisotropic_mesh_parameters param;
    Criteria* criteria = NULL;
    typedef CGAL::Anisotropic_mesh_3::Metric_field<Kernel> Metric_field;
    Metric_field* mf = NULL;

    if( !keep_same_parameters )
    {
      criteria = build_param_and_metric(poly_surf, param, mf, epsilon, approximation,
                                        radius_edge_ratio, sliverity, circumradius,
                                        distortion, beta, delta, max_times_to_try_in_picking_region,
                                        dim, metric, en_factor);

      ssetitem->star_set().set_criteria(criteria);
      {// can't handle metric changes for now
        //ssetitem->star_set().set_metric_field(mf);
        delete mf;
      }
      ssetitem->star_set().update_stars_criteria(); //update criteria in stars
      ssetitem->star_set().fill_refinement_queue(); //reset and refill since criteria has been modified
    }
    else
    {
      param.dim = dim;
    }

    typedef Anisotropic_mesh_function<Constrain_surface_polyhedral, Metric_field> AMesh_function;
    AMesh_function* p_mesh_function = new AMesh_function(ssetitem->star_set(), param,
                                         ssetitem->star_set().criteria(),
                                         ssetitem->star_set().metric_field(),
                                         pick_valid_causes_stop, pick_valid_use_cube_probing,
                                         pick_valid_max_failures);

    thread = new Anisotropic_meshing_thread(p_mesh_function, ssetitem);
  }
  //// FUNCTION -------------------------------------------------------------------------------
  else if( NULL != function_surf )
  {
    Anisotropic_mesh_parameters param;
    Criteria* criteria = NULL;
    typedef CGAL::Anisotropic_mesh_3::Metric_field<Kernel> Metric_field;
    Metric_field* mf = NULL;

    if( !keep_same_parameters )
    {
      criteria = build_param_and_metric(function_surf, param, mf, epsilon, approximation,
                                        radius_edge_ratio, sliverity, circumradius,
                                        distortion, beta, delta, max_times_to_try_in_picking_region,
                                        dim, metric, en_factor);

      ssetitem->star_set().set_criteria(criteria);
      {// can't handle metric changes for now
        //ssetitem->star_set().set_metric_field(mf);
        delete mf;
      }
      ssetitem->star_set().update_stars_criteria(); //update criteria in stars
      ssetitem->star_set().fill_refinement_queue(); //reset and refill since criteria has been modified
    }
    else
    {
      param.dim = dim;
    }

    typedef Anisotropic_mesh_function<Implicit_surface, Metric_field> AMesh_function;
    AMesh_function* p_mesh_function = new AMesh_function(ssetitem->star_set(), param,
                                         ssetitem->star_set().criteria(),
                                         ssetitem->star_set().metric_field(),
                                         pick_valid_causes_stop, pick_valid_use_cube_probing,
                                         pick_valid_max_failures);

    thread = new Anisotropic_meshing_thread(p_mesh_function, ssetitem);
  }

  if ( NULL == thread )
  {
    QMessageBox::critical(mw,tr(""),tr("ERROR: no thread created"));
    return;
  }

  // Launch thread
  source_item_ = item;
  launch_thread(thread, true /*resumed*/);
  QApplication::restoreOverrideCursor();
}

void Anisotropic_mesh_3_plugin::view_one_star()
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();

  // Get item
  Scene_starset3_item* ssetitem =
    qobject_cast<Scene_starset3_item*>(scene->item(index));
  if( NULL == ssetitem )
  {
    QMessageBox::warning(mw,tr(""),tr("Selected object is not a star set.")); 
    return;
  }

  // Create dialog
  QDialog dialog(mw);
  Ui::Choose_star_dialog ui;
  ui.setupUi(&dialog);
  connect(ui.buttonBox, SIGNAL(accepted()), &dialog, SLOT(accept()));
  connect(ui.buttonBox, SIGNAL(rejected()), &dialog, SLOT(reject()));
  connect(ui.displayAll, SIGNAL(clicked()), this, SLOT(view_all_stars()));
  connect(ui.displayAll, SIGNAL(clicked()), &dialog, SLOT(reject()));

   // Set default parameters
  int indexMax = ssetitem->nbStars()-1;
  ui.objectQuestion->setText(
    tr("Which star do you want to display? (from 0 to %1)").arg(indexMax));
  ui.starId->setDecimals(0);
  ui.starId->setSingleStep(1);
  ui.starId->setValue((std::max)(1, ssetitem->draw_star()));
  ui.starId->setMinimum(0);
  ui.starId->setMaximum(indexMax);
   
  // Get value and use it
  int i = dialog.exec();  
  if(i != QDialog::Rejected) 
    ssetitem->draw_star() = (int)ui.starId->value();  
  
  ssetitem->starset_changed();
}

void Anisotropic_mesh_3_plugin::view_all_stars()
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  Scene_starset3_item* ssetitem =
    qobject_cast<Scene_starset3_item*>(scene->item(index));
  if(ssetitem != NULL)
  {
    ssetitem->draw_star() = -1;
    ssetitem->starset_changed();
  }
}

void Anisotropic_mesh_3_plugin::view_one_pickvalid_point()
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();

  // Get item
  Scene_starset3_item* ssetitem =
    qobject_cast<Scene_starset3_item*>(scene->item(index));
  if( NULL == ssetitem )
  {
    QMessageBox::warning(mw,tr(""),tr("Selected object is not a star set."));
    return;
  }

  if(!ssetitem->star_set().pick_valid_failed())
  {
    QMessageBox::warning(mw,tr(""),tr("Nothing has failed. Yet."));
    return;
  }

  // Create dialog
  QDialog dialog(mw);
  Ui::Choose_pickvalid_point_dialog ui;
  ui.setupUi(&dialog);
  connect(ui.buttonBox, SIGNAL(accepted()), &dialog, SLOT(accept()));
  connect(ui.buttonBox, SIGNAL(rejected()), &dialog, SLOT(reject()));
  connect(ui.displayAll, SIGNAL(clicked()), this, SLOT(view_all_pickvalid_point()));
  connect(ui.displayAll, SIGNAL(clicked()), &dialog, SLOT(reject()));

   // Set default parameters
  int indexMax = ((int)ssetitem->star_set().pickvalid_problematic_facets_size())-1;
  ui.objectQuestion->setText(
    tr("Which point do you want to display? (from 0 to %1)").arg(indexMax));
  ui.PickvalidPointId->setDecimals(0);
  ui.PickvalidPointId->setSingleStep(1);
  ui.PickvalidPointId->setValue(0);
  ui.PickvalidPointId->setMinimum(0);
  ui.PickvalidPointId->setMaximum(indexMax);

  // Get value and use it
  int i = dialog.exec();
  if(i != QDialog::Rejected)
    ssetitem->draw_pickvalid_point() = (int)ui.PickvalidPointId->value();

  ssetitem->starset_changed();
}

void Anisotropic_mesh_3_plugin::view_all_pickvalid_point()
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  Scene_starset3_item* ssetitem = 
    qobject_cast<Scene_starset3_item*>(scene->item(index));
  if(ssetitem != NULL)
  {
    ssetitem->draw_pickvalid_point() = -1;
    ssetitem->starset_changed();
  }
}

void Anisotropic_mesh_3_plugin::view_surface_star_set(bool b)
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  Scene_starset3_item* ssetitem =
    qobject_cast<Scene_starset3_item*>(scene->item(index));
  if(ssetitem != NULL)
  {
    ssetitem->draw_surface_star_set() = b;
    ssetitem->starset_changed();
  }
}

void Anisotropic_mesh_3_plugin::view_cell(bool b)
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  Scene_starset3_item* ssetitem =
    qobject_cast<Scene_starset3_item*>(scene->item(index));
  if(ssetitem != NULL)
  {
    ssetitem->draw_cell() = b;
    ssetitem->starset_changed();
  }
}

void Anisotropic_mesh_3_plugin::view_dual_edges(bool b)
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  Scene_starset3_item* ssetitem =
    qobject_cast<Scene_starset3_item*>(scene->item(index));
  if(ssetitem != NULL)
  {
    ssetitem->draw_dual() = b;
    ssetitem->starset_changed();
  }
}

void Anisotropic_mesh_3_plugin::view_poles(bool b)
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  Scene_starset3_item* ssetitem =
    qobject_cast<Scene_starset3_item*>(scene->item(index));
  if(ssetitem != NULL)
  {
    ssetitem->draw_poles() = b;
    ssetitem->starset_changed();
  }
}

void Anisotropic_mesh_3_plugin::view_initial_points(bool b)
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  Scene_starset3_item* ssetitem =
    qobject_cast<Scene_starset3_item*>(scene->item(index));
  if(ssetitem != NULL)
  {
    ssetitem->draw_initial_points() = b;
    ssetitem->starset_changed();
  }
}

void Anisotropic_mesh_3_plugin::view_surface_delaunay_balls(bool b)
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  Scene_starset3_item* ssetitem =
    qobject_cast<Scene_starset3_item*>(scene->item(index));
  if(ssetitem != NULL)
  {
    ssetitem->draw_surface_delaunay_balls() = b;
    ssetitem->starset_changed();
  }
}

void Anisotropic_mesh_3_plugin::view_inconsistent_facets(bool b)
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  Scene_starset3_item* ssetitem = 
    qobject_cast<Scene_starset3_item*>(scene->item(index));
  if(ssetitem != NULL)
  {
    ssetitem->draw_inconsistent_facets() = b;
    ssetitem->starset_changed();
  }
}

void Anisotropic_mesh_3_plugin::view_distortion(bool b)
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  Scene_starset3_item* ssetitem = 
    qobject_cast<Scene_starset3_item*>(scene->item(index));
  if(ssetitem != NULL)
  {
    ssetitem->draw_distortion() = b;
    ssetitem->starset_changed();
  }
}

void Anisotropic_mesh_3_plugin::view_metric_field(bool b)
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  Scene_starset3_item* ssetitem = 
    qobject_cast<Scene_starset3_item*>(scene->item(index));
  if(ssetitem != NULL)
  {
    ssetitem->draw_metric_field() = b;
    ssetitem->starset_changed();
  }
}

void Anisotropic_mesh_3_plugin::view_mesh_3(bool b)
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  Scene_starset3_item* ssetitem = 
    qobject_cast<Scene_starset3_item*>(scene->item(index));
  if(ssetitem != NULL)
  {
    ssetitem->draw_mesh_3() = b;
    ssetitem->starset_changed();
  }
}

void
Anisotropic_mesh_3_plugin::
launch_thread(Anisotropic_meshing_thread* mesh_thread,
              const bool resumed /*= false*/)
{
  // -----------------------------------
  // Create message box with stop button
  // -----------------------------------
  message_box_ = new QMessageBox(QMessageBox::NoIcon,
                                 "Anisotropic meshing",
                                 "Anisotropic mesh generation in progress...",
                                 QMessageBox::Cancel,
                                 mw);
  
  message_box_->setDefaultButton(QMessageBox::Cancel);
  QAbstractButton* cancelButton = message_box_->button(QMessageBox::Cancel);
  cancelButton->setText(tr("Stop"));
  
  QObject::connect(cancelButton, SIGNAL(clicked()),
                   mesh_thread,  SLOT(stop()));
  
  message_box_->show();
  
  // --------------------------------------
  // Connect main thread to meshing thread
  // --------------------------------------

  if(resumed)
  {
    QObject::connect(mesh_thread, SIGNAL(done(Anisotropic_meshing_thread*)),
                     this, SLOT(meshing_done_resume(Anisotropic_meshing_thread*)));
  }
  else
  {
    QObject::connect(mesh_thread, SIGNAL(done(Anisotropic_meshing_thread*)),
                     this, SLOT(meshing_done_nresume(Anisotropic_meshing_thread*)));
  }
  QObject::connect(mesh_thread, SIGNAL(status_report(QString)),
                   this, SLOT(status_report(QString)));
  
  // -----------------------------------
  // Launch mesher
  // -----------------------------------
  mesh_thread->start();
}

void
Anisotropic_mesh_3_plugin::
status_report(QString str)
{
  if ( NULL == message_box_ ) { return; }
  
  message_box_->setInformativeText(str);
}

void
Anisotropic_mesh_3_plugin::
meshing_done_resume(Anisotropic_meshing_thread* thread)
{
  meshing_done(thread, true);
}

void
Anisotropic_mesh_3_plugin::
meshing_done_nresume(Anisotropic_meshing_thread* thread)
{
  meshing_done(thread, false);
}

void
Anisotropic_mesh_3_plugin::
meshing_done(Anisotropic_meshing_thread* thread, const bool resumed)
{
  // Print message in console
  QString str = QString("Anisotropic meshing of \"%1\" done in %2s<br>")
    .arg(source_item_->name())
    .arg(thread->time());
  
  Q_FOREACH( QString param, thread->parameters_log() )
  {
    str.append(QString("( %1 )<br>").arg(param));
  }
  
  msg->information(qPrintable(str));
  
  // Treat new starset3 item
  Scene_starset3_item* result_item = thread->item();
  treat_result(*source_item_, *result_item, resumed);
  
  // close message box
  message_box_->close();
  message_box_ = NULL;
  
  // free memory
  // TODO: maybe there is another way to do that
  delete thread;
}


void
Anisotropic_mesh_3_plugin::
treat_result(Scene_item& source_item,
             Scene_starset3_item& result_item,
             const bool replacing /*= false*/) const
{
  if(!replacing)
    result_item.setName(tr("%1 [3D Anisotropic Mesh]").arg(source_item.name()));
  else
    result_item.setName(source_item.name());

  result_item.starset_changed();
  
  const Scene_item::Bbox& bbox = result_item.bbox();
  result_item.setPosition((bbox.xmin + bbox.xmax)/2.f,
                          (bbox.ymin + bbox.ymax)/2.f,
                          (bbox.zmin + bbox.zmax)/2.f);
    
  result_item.setColor(default_mesh_color);
  result_item.setRenderingMode(source_item.renderingMode());
  result_item.set_data_item(&source_item);
  
  source_item.setVisible(false);
  result_item.setVisible(true);
    
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  scene->itemChanged(index);

  if(replacing){
    scene->replaceItem(index, &result_item);
    scene->setSelectedItem(index);
  }
  else{
    Scene_interface::Item_id new_item_id = scene->addItem(&result_item);
    scene->setSelectedItem(new_item_id);
  }
}


double
get_approximate(double d, int precision, int& decimals)
{
  if ( d<0 ) { return 0; }
  
  double i = std::pow(10.,precision-1);
  
  decimals = 0;
  while ( d > i*10 ) { d = d/10.; ++decimals; }
  while ( d < i ) { d = d*10.; --decimals; }
  
  return std::floor(d)*std::pow(10.,decimals);
}


Q_EXPORT_PLUGIN2(Anisotropic_mesh_3_plugin, Anisotropic_mesh_3_plugin)

#include "mesh_3_anisotropic_plugin.moc"
