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

#include <CGAL/Anisotropic_surface_mesher_3.h>
//#include <Anisotropic_mesh_3/include/CGAL/Anisotropic_tet_mesher_3.h>

#include <Domain/Constrain_surface_3_sphere.h>

//#include "Anisotropic_mesh_3/include/CGAL/Polyhedral_curvature_metric_field.h"
//#include "Anisotropic_mesh_3/include/CGAL/Metric_field.h"
#include <Metric_field/Arctan_metric_field.h>
#include <CGAL/Euclidean_metric_field.h>

#include <StarSet_type.h>
#include <Anisotropic_meshing_thread.h>


struct Anisotropic_mesh_parameters
{
  double approximation;
  double radius_edge_ratio;
  double sliverity;
  double circumradius;
  double distortion;
  double beta;
  double delta;
  int max_times_to_try_in_picking_region;
  int dim;
  
  inline QStringList log() const;
};


template < typename Domain_, typename Metric_field >
class Anisotropic_mesh_function
  : public Anisotropic_mesh_function_interface
{
  typedef Domain_ Domain;

  typedef Surface_star_set::Criteria Criteria;
//  typedef Surface_star_set::Metric_field Metric_field;
  
public:
  Anisotropic_mesh_function(Surface_star_set& starset,
                            //Domain* domain, 
                            const Anisotropic_mesh_parameters& p,
                            Criteria* criteria,
                            Metric_field* metrix_field);
  // Note: 'this' takes the ownership of 'criteria' and 'metrix_field'.
  
  ~Anisotropic_mesh_function();
  
  // Launch
  virtual void launch();
  
  // Stop
  virtual void stop();
  
  // Logs
  virtual QStringList parameters_log() const;
  virtual QString status(double time_period) const;

  typedef CGAL::Anisotropic_mesh_3::Anisotropic_surface_mesher_3<Kernel> SMesher;
//  typedef CGAL::Anisotropic_mesh_3::Anisotropic_tet_mesher_3<Kernel> TMesher;

private:
  Surface_star_set& starset_;
//  Domain* domain_;
  Anisotropic_mesh_parameters p_;
  bool continue_;
  
  SMesher* smesher_;

  Criteria* criteria_;
  Metric_field* metrix_field_;

//  TMesher* tmesher_;
//  mutable typename Mesher::Mesher_status last_report_;
};



// -----------------------------------
// Class Mesh_parameters
// -----------------------------------
inline
QStringList
Anisotropic_mesh_parameters::
log() const
{
  return QStringList()
  << QString("approximation: %1").arg(approximation)
  << QString("radius edge ratio: %1").arg(radius_edge_ratio)
  << QString("sliverity: %1").arg(sliverity)
  << QString("circumradius: %1").arg(circumradius)
  << QString("distortion: %1").arg(distortion)
  << QString("beta: %1").arg(beta)
  << QString("delta: %1").arg(delta)
  << QString("max tries: %1").arg(max_times_to_try_in_picking_region)
  << QString("dimension: %1").arg(dim);
}


// -----------------------------------
// Class Mesh_function
// -----------------------------------
template < typename D_, typename Metric_field>
Anisotropic_mesh_function<D_, Metric_field>::Anisotropic_mesh_function(
  Surface_star_set& starset, 
//  Domain* domain, 
  const Anisotropic_mesh_parameters& param,
  Criteria* criteria,
  Metric_field* metrix_field)
: starset_(starset) 
, p_(param)
, continue_(true)
, smesher_(NULL)
, criteria_(criteria)
, metrix_field_(metrix_field)
//, domain_(domain)
//, tmesher_(NULL)
//, last_report_(0,0,0)
{
}


template < typename D_, typename Metric_field>
Anisotropic_mesh_function<D_, Metric_field>::
~Anisotropic_mesh_function()
{
  delete metrix_field_;
  delete criteria_;
  delete smesher_;
//  delete domain_;
//  delete tmesher_;
}


template < typename D_, typename Metric_field>
void
Anisotropic_mesh_function<D_, Metric_field>::
launch()
{
  // Build mesher and launch refinement process
  if(p_.dim == 2)
  {
    smesher_ = new SMesher(starset_);
    //smesher_->refine_all();//p_.max_times_to_try_in_picking_region);

//------------------------------------------------------------------------
//The following is a copy of refine_all from surface_star_set_3.h
//to use a continue boolean. If you edit the algorithm here,
//make sure to edit in the above file as well.

    const int max_count = INT_MAX;//p_.max_times_to_try_in_picking_region);

#ifdef ANISO_VERBOSE
    std::cout << "\nRefine all...";
    std::clock_t start_time = clock();
#endif
    smesher_->star_set.fill_refinement_queue();

#ifdef ANISO_VERBOSE
    smesher_->star_set.vertex_with_picking_count = 0;
    smesher_->star_set.vertex_without_picking_count =
            (int)smesher_->star_set.m_stars.size();
    std::cout << "There are ";
    std::cout << smesher_->star_set.count_restricted_facets();
    std::cout << " restricted facets.\n";
#endif
    std::size_t nbv = smesher_->star_set.m_stars.size();
    while(nbv < max_count && continue_)
    {
      if(nbv == 4) //dimension 3 reached
        smesher_->star_set.update_bboxes();

      if(nbv % 100 == 0)
      {
        smesher_->star_set.clean_stars();//remove useless vertices
#ifdef ANISO_VERBOSE
        std::cerr << " " << nbv << " vertices, ";
        std::cerr << smesher_->star_set.duration(start_time) << " sec.,\t";
        smesher_->star_set.m_refine_queue.print();
#endif
#ifdef ANISO_DEBUG
        std::ostringstream oss;
        oss << "out_" << nbv << ".off";
        smesher_->star_set.output(oss.str().c_str(), false/*consistent_only*/);
#endif
       }

       if(!smesher_->star_set.refine())
       {
         smesher_->star_set.clean_stars();
         //debug_show_distortions();
         break;
       }
       nbv = smesher_->star_set.m_stars.size();
    }

#ifdef ANISO_VERBOSE
    double time = smesher_->star_set.duration(start_time);
    std::cout << "\nRefinement done (" << nbv << " vertices in " << time << " seconds)\n";
    if(smesher_->star_set.is_consistent(true/*verbose*/))
      std::cout << "Triangulation is consistent.\n";
    else
      std::cout << "Triangulation is not consistent.\n";

    std::cout << "Vertices via picking: ";
    std::cout << smesher_->star_set.vertex_with_picking_count << std::endl;
    std::cout << "Vertex non-picking: ";
    std::cout << smesher_->star_set.vertex_without_picking_count << std::endl;
    std::cout << "picking rate:       ";
    std::cout << (double)smesher_->star_set.vertex_with_picking_count /
      (double)(smesher_->star_set.vertex_without_picking_count +
               smesher_->star_set.vertex_with_picking_count) << std::endl;
    std::cout << "Approximation error : ";
    std::cout << smesher_->star_set.compute_approximation_error() << std::endl;

    smesher_->star_set.report();
    histogram_vertices_per_star<Surface_star_set>(smesher_->star_set);
#endif
#ifdef USE_ANISO_TIMERS
    report_timers();
#endif

// end of copy ----------------------------------------------------------

  }
  else if(p_.dim == 3)
  {
    std::cout << "Tet mesher : ToDo !! \n";
    //tmesher_ = new TMesher(*domain_, metric_field, criteria);
    //tmesher_->refine();
  }
  else
    std::cerr << "ToDo : Anisotropic meshing for volume + surface.\n";
}


template < typename D_, typename Metric_field>
void
Anisotropic_mesh_function<D_, Metric_field>::
stop()
{
  continue_ = false;
}


template < typename D_, typename Metric_field>
QStringList
Anisotropic_mesh_function<D_, Metric_field>::
parameters_log() const
{
  return p_.log();
}


template < typename D_, typename Metric_field>
QString
Anisotropic_mesh_function<D_, Metric_field>::
status(double time_period) const
{
  // If smesher_ is not yet created, it means that either launch() has not
  // been called or that initial points have not been found
  if ( NULL == smesher_ )
  {
    return QString("Initialization in progress...");
  }
  
  // Get status and return a string corresponding to it
//  typename Mesher::Mesher_status s = smesher_->status();
  
  QString result;
  /* = QString("Vertices: %1 \n")
                           "Vertices inserted last %2s: %3 \n\n")
    .arg(s.vertices);
    .arg(time_period)
    .arg(s.vertices - last_report_.vertices);
  
  last_report_ = s;
  */
  return result;
}

#endif // CGAL_DEMO_ANISOTROPIC_MESH_3_MESH_FUNCTION_H
