//
// Author(s)     : Jane Tournois, Stephane Tayeb
//
//******************************************************************************
// File Description : 
//******************************************************************************

#ifndef CGAL_DEMO_ANISOTROPIC_MESH_3_MESH_FUNCTION_H
#define CGAL_DEMO_ANISOTROPIC_MESH_3_MESH_FUNCTION_H

#define CGAL_MESH_3_ANISO_MESHER_STATUS_ACTIVATED 1

#include <QStringList>
#include <QString>

#include <Anisotropic_meshing_thread.h>
#include <StarSet_type.h>

#include <CGAL/Anisotropic_mesher_3.h>
#include <CGAL/Anisotropic_surface_mesher_3.h>
#include <CGAL/Anisotropic_tet_mesher_3.h>

template < typename Domain_, typename Metric_field >
class Anisotropic_mesh_function
  : public Anisotropic_mesh_function_interface
{
  typedef CGAL::Anisotropic_mesh_3::Anisotropic_mesher_3_base            Mesher_base;
  typedef CGAL::Anisotropic_mesh_3::Anisotropic_mesher_3<Kernel>         AMesher;
  typedef CGAL::Anisotropic_mesh_3::Anisotropic_surface_mesher_3<Kernel> ASMesher;
  typedef CGAL::Anisotropic_mesh_3::Anisotropic_tet_mesher_3<Kernel>     ATMesher;

private:
  Starset_with_info& m_starset_with_info;
  bool m_continue;
  Mesher_base* m_amesher;

public:
  Anisotropic_mesh_function(Starset_with_info&);
  ~Anisotropic_mesh_function();

  // Launch
  virtual void launch();

  // Stop
  virtual void stop();

  // Logs
  virtual QStringList parameters_log() const;
  virtual QString status(double time_period) const;
};

// -----------------------------------
// Class Mesh_function
// -----------------------------------
template < typename D_, typename Metric_field>
Anisotropic_mesh_function<D_, Metric_field>::
Anisotropic_mesh_function(Starset_with_info& ss_with_info_)
:
  m_starset_with_info(ss_with_info_)
, m_continue(true)
, m_amesher(NULL)
{ }

template < typename D_, typename Metric_field>
Anisotropic_mesh_function<D_, Metric_field>::
~Anisotropic_mesh_function()
{
  delete m_amesher;
}

template < typename D_, typename Metric_field>
void
Anisotropic_mesh_function<D_, Metric_field>::
launch()
{
  int dimension = m_starset_with_info.criteria()->dimension;
  if(dimension == 2)
    m_amesher = new ASMesher(m_starset_with_info,
                             m_starset_with_info.constrain_surface(),
                             m_starset_with_info.criteria(),
                             m_starset_with_info.metric_field());
  else if(dimension ==3)
    m_amesher = new ATMesher(m_starset_with_info,
                             m_starset_with_info.constrain_surface(),
                             m_starset_with_info.criteria(),
                             m_starset_with_info.metric_field());
  else if(dimension == 6)
    m_amesher = new AMesher(m_starset_with_info,
                            m_starset_with_info.constrain_surface(),
                            m_starset_with_info.criteria(),
                            m_starset_with_info.metric_field());

  m_amesher->initialize();

  while(!m_amesher->is_algorithm_done() && m_continue)
  {
    m_amesher->one_step();
  }

  m_amesher->report();
}


template < typename D_, typename Metric_field>
void
Anisotropic_mesh_function<D_, Metric_field>::
stop()
{
  m_continue = false;
}


template < typename D_, typename Metric_field>
QStringList
Anisotropic_mesh_function<D_, Metric_field>::
parameters_log() const
{
  const Criteria* criteria = m_starset_with_info.criteria();

  return QStringList()
  << QString("general criteria")
  << QString("approximation error : %1").arg(criteria->approximation)
  << QString("maximum distortion (gamma_0): %1").arg(criteria->distortion)
  << QString("checking threshold (beta): %1").arg(criteria->beta)
  << QString("picking region radius (delta): %1").arg(criteria->delta)
  << QString("max. time to try in picking region: %1").arg(criteria->max_times_to_try_in_picking_region)

  << QString("facet criteria :")
  << QString("radius-edge-ratio (rho_0): %1").arg(criteria->facet_radius_edge_ratio)
  << QString("maximum circumradius (r_0): %1").arg(criteria->facet_circumradius)

  << QString("cell criteria : ")
  << QString("sliverity ratio (sigma_0): %1").arg(criteria->sliverity)
  << QString("radius-edge-ratio (rho_0): %1").arg(criteria->cell_radius_edge_ratio)
  << QString("maximum circumradius (r_0): %1").arg(criteria->cell_circumradius);
}

template < typename D_, typename Metric_field>
QString
Anisotropic_mesh_function<D_, Metric_field>::
status(double time_period) const
{
  // If m_amesher is not yet created, it means that either launch() has not
  // been called or that initial points have not been found
  if ( NULL == m_amesher )
  {
    return QString("Initialization in progress...");
  }

  // Get status and return a string corresponding to it
//  typename Mesher::Mesher_status s = m_amesher->status();

  QString result;
  /* = QString("Vertices: %1 \n")
                           "Vertices inserted last %2s: %3 \n\n")
    .arg(s.vertices);
    .arg(time_period)
    .arg(s.vertices - m_last_report.vertices);

  m_last_report = s;
  */
  return result;
}

#endif // CGAL_DEMO_ANISOTROPIC_MESH_3_MESH_FUNCTION_H
