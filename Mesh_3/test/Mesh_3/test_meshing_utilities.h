// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
// Author(s)     : Stephane Tayeb
//
//******************************************************************************
// File Description : Test meshing utilities.
//******************************************************************************

#include <CGAL/Cartesian.h>
#include <CGAL/Simple_cartesian.h>

#include "test_utilities.h"

#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/refine_mesh_3.h>
#include <CGAL/optimize_mesh_3.h>

#include <vector>

// IO
#include <fstream>
#include <iostream>


template <typename K>
struct Tester
{
  template<typename C3t3, typename Domain, typename Criteria>
  void verify(C3t3& c3t3,
              const Domain& domain,
              const Criteria& criteria,
              const unsigned int min_vertices_expected = 0,
              const unsigned int max_vertices_expected = (unsigned int)(-1),
              const unsigned int min_facets_expected = 0,
              const unsigned int max_facets_expected = (unsigned int)(-1),
              const unsigned int min_cells_expected = 0,
              const unsigned int max_cells_expected = (unsigned int)(-1) ) const
  {
    typedef typename C3t3::size_type size_type;
    typedef typename C3t3::Triangulation::Geom_traits Gt;
    typedef typename CGAL::Mesh_3::Min_dihedral_angle_criterion<Gt> Sliver_criterion;
    typedef typename CGAL::Mesh_3::Slivers_exuder<C3t3, Sliver_criterion> Exuder;

    // Store mesh properties
    size_type v = c3t3.triangulation().number_of_vertices();
    size_type f = c3t3.number_of_facets();
    size_type c = c3t3.number_of_cells();

    // Verify
    verify_c3t3(c3t3,
                min_vertices_expected,
                max_vertices_expected,
                min_facets_expected,
                max_facets_expected,
                min_cells_expected,
                max_cells_expected);

    // Refine again and verify nothing changed
    std::cerr << "Refining again...\n";
    refine_mesh_3(c3t3,domain,criteria,
                  CGAL::parameters::no_exude(),
                  CGAL::parameters::no_perturb());

    verify_c3t3(c3t3,v,v,f,f,c,c);
    
    // Exude.
    // Vertex number should not change (obvious)
    // Facet number should not change as exuder preserves boundary facets
    // Quality should increase
    C3t3 exude_c3t3(c3t3);
    std::cerr << "Exude...\n";
    CGAL::exude_mesh_3(exude_c3t3);
    verify_c3t3(exude_c3t3,v,v,f,f);
    verify_c3t3_quality(c3t3,exude_c3t3);
    
    // Perturb.
    // Vertex number should not change (obvious)
    // Quality should increase
    C3t3 perturb_c3t3(c3t3);
    std::cerr << "Perturb...\n";
    CGAL::perturb_mesh_3(perturb_c3t3, domain, CGAL::parameters::time_limit=5);
    verify_c3t3(perturb_c3t3,v,v);
    verify_c3t3_quality(c3t3,perturb_c3t3);
    
    // Odt-smoothing
    // Vertex number should not change (obvious)
    C3t3 odt_c3t3(c3t3);
    std::cerr << "Odt...\n";
    CGAL::odt_optimize_mesh_3(odt_c3t3, domain, CGAL::parameters::time_limit=5,
                              CGAL::parameters::convergence=0.001, CGAL::parameters::sliver_bound=0.0005);
    verify_c3t3(odt_c3t3,v,v);
    
    // Lloyd-smoothing
    // Vertex number should not change (obvious)
    C3t3 lloyd_c3t3(c3t3);
    std::cerr << "Lloyd...\n";
    CGAL::lloyd_optimize_mesh_3(lloyd_c3t3, domain, CGAL::parameters::time_limit=5,
                                CGAL::parameters::convergence=0.001, CGAL::parameters::sliver_bound=0.0005);
    verify_c3t3(lloyd_c3t3,v,v);
  }

  template<typename C3t3>
  void verify_c3t3(const C3t3& c3t3,
                   const unsigned int min_vertices_expected = 0,
                   const unsigned int max_vertices_expected = (unsigned int)(-1),
                   const unsigned int min_facets_expected = 0,
                   const unsigned int max_facets_expected = (unsigned int)(-1),
                   const unsigned int min_cells_expected = 0,
                   const unsigned int max_cells_expected = (unsigned int)(-1)) const
  {
    //-------------------------------------------------------
    // Verifications
    //-------------------------------------------------------
    std::cerr << "\tNumber of cells: " << c3t3.number_of_cells() << "\n";
    std::cerr << "\tNumber of facets: " << c3t3.number_of_facets() << "\n";
    std::cerr << "\tNumber of vertices: " << c3t3.number_of_vertices() << "\n";
        
    std::size_t dist_facets ( std::distance(c3t3.facets_begin(), 
                                            c3t3.facets_end()) );
    std::size_t dist_cells ( std::distance(c3t3.cells_begin(), 
                                            c3t3.cells_end()) );

    assert(min_vertices_expected <= c3t3.number_of_vertices());
    assert(max_vertices_expected >= c3t3.number_of_vertices());

    assert(min_facets_expected <= c3t3.number_of_facets());
    assert(max_facets_expected >= c3t3.number_of_facets());
    assert(dist_facets == c3t3.number_of_facets());
    
    assert(min_cells_expected <= c3t3.number_of_cells());
    assert(max_cells_expected >= c3t3.number_of_cells());
    assert(dist_cells == c3t3.number_of_cells());
  }
  
  template<typename C3t3>
  void verify_c3t3_quality(const C3t3& original_c3t3,
                           const C3t3& modified_c3t3) const
  {
    double original = min_value(original_c3t3);
    double modified = min_value(modified_c3t3);
    
    std::cout << "\tQuality before optimization: " << original
              << " - Quality after optimization: " << modified << std::endl;
    
    assert(original <= modified);
  }
  
  template<typename C3t3>
  double min_value(const C3t3& c3t3) const
  {
    typedef typename C3t3::Triangulation::Geom_traits                 Gt;
    typedef typename CGAL::Mesh_3::Min_dihedral_angle_criterion<Gt>   Criterion;
    typedef typename C3t3::Cell_iterator                              Cell_iterator;
    
    Criterion criterion;
    double min_value = Criterion::max_value;
    
    for ( Cell_iterator cit = c3t3.cells_begin(), end = c3t3.cells_end() ;
         cit != end ; ++cit )
    {
      min_value = (std::min)(min_value, criterion(c3t3.triangulation().tetrahedron(cit)));
    }
    
    return min_value;
  }
};

