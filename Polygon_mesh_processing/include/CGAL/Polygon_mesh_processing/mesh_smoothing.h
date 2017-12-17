// Copyright (c) 2017 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
// Author(s)     : Konstantinos Katrioplas (konst.katrioplas@gmail.com)

#ifndef CGAL_POLYGON_MESH_PROCESSING_SMOOTHING_H
#define CGAL_POLYGON_MESH_PROCESSING_SMOOTHING_H

#include <CGAL/Polygon_mesh_processing/internal/Smoothing/mesh_smoothing_impl.h>
#include <CGAL/Polygon_mesh_processing/internal/Smoothing/curvature_flow_impl.h>
#include <CGAL/Polygon_mesh_processing/internal/Isotropic_remeshing/remesh_impl.h>
#include <CGAL/Polygon_mesh_processing/internal/Smoothing/evaluation.h>

#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>

#include <CGAL/Polygon_mesh_processing/distance.h>
#define CGAL_PMP_SMOOTHING_VERBOSE

#ifdef CGAL_PMP_SMOOTHING_VERBOSE
#include <CGAL/Timer.h>
#endif

#if defined(CGAL_LINKED_WITH_TBB)
#define TAG CGAL::Parallel_tag
#else
#define TAG CGAL::Sequential_tag
#endif

namespace CGAL {
namespace Polygon_mesh_processing {

/*!
* \ingroup PMP_meshing_grp
* smoothes a triangulated region of a polygon mesh using angle based criteria.
* This function improves the angles of triangle faces by iteratively moving non
* constrained vertices.
* A weighted scheme is also proposed to favorise the removal of small angles.
* Optionally, the points are reprojected after each iteration.
*
* @todo Make the projection of points optional.
* @tparam PolygonMesh model of `MutableFaceGraph`.
*         If `PolygonMesh` has an internal property map for `CGAL::face_index_t`,
*         and no `face_index_map` is given
*         as a named parameter, then the internal one should be initialized.
* @tparam FaceRange range of `boost::graph_traits<PolygonMesh>::%face_descriptor`,
          model of `Range`. Its iterator type is `ForwardIterator`.
* @tparam NamedParameters a sequence of \ref namedparameters.
*
* @param pmesh a polygon mesh with triangulated surface patches to be smoothed.
* @param faces the range of triangular faces defining one or several surface patches to be smoothed.
* @param np optional sequence of \ref namedparameters among the ones listed below.
*
* \cgalNamedParamsBegin
*  \cgalParamBegin{geom_traits} a geometric traits class instance, model of `Kernel`.
*    Exact constructions kernels are not supported by this function.
*  \cgalParamEnd
*  \cgalParamBegin{vertex_point_map} the property map with the points associated
*    to the vertices of `pmesh`. Instance of a class model of `ReadWritePropertyMap`.
*  \cgalParamEnd
*  \cgalParamBegin{face_index_map} a property map containing the index of each face of `pmesh`.
*  \cgalParamEnd
*  \cgalParamBegin{number_of_iterations} the number of iterations for the
*    sequence of the smoothing iterations performed.
*  \cgalParamEnd
*  \cgalParamBegin{edge_is_constrained_map} a property map containing the
*    constrained-or-not status of each edge of `pmesh`. Vertices that belong to constrained
*    edges are not modified at all during smoothing.
*  \cgalParamEnd
*  \cgalParamBegin{vertex_is_constrained_map} a property map containing the
*    constrained-or-not status of each vertex of `pmesh`. A constrained vertex
*    cannot be modified at all during smoothing.
*  \cgalParamEnd
*  \cgalParamBegin{use_weights} If `true`, small angles carry more weight than larger ones.
*  \cgalParamEnd
* \cgalNamedParamsEnd
*/
template<typename PolygonMesh, typename FaceRange, typename NamedParameters>
void smooth_angles(const FaceRange& faces, PolygonMesh& pmesh, const NamedParameters& np)
{
    using boost::choose_param;
    using boost::get_param;

#ifdef CGAL_PMP_SMOOTHING_VERBOSE
    CGAL::Timer t;
    std::cout << "Smoothing parameters...";
    std::cout.flush();
    t.start();
#endif

    // geom traits
    typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type GeomTraits;

    // vpmap
    typedef typename GetVertexPointMap<PolygonMesh, NamedParameters>::type VertexPointMap;
    VertexPointMap vpmap = choose_param(get_param(np, internal_np::vertex_point),
                                 get_property_map(CGAL::vertex_point, pmesh));

    // fimap
    typedef typename GetFaceIndexMap<PolygonMesh, NamedParameters>::type FIMap;
    FIMap fimap = choose_param(get_param(np, internal_np::face_index),
                             get_property_map(face_index, pmesh));

    // vcmap
    typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;
    typedef typename boost::lookup_named_param_def <
        internal_np::vertex_is_constrained_t,
        NamedParameters,
        internal::No_constraint_pmap<vertex_descriptor>//default
      > ::type VCMap;
    VCMap vcmap = choose_param(get_param(np, internal_np::vertex_is_constrained),
                               internal::No_constraint_pmap<vertex_descriptor>());

    // ecmap
    typedef typename boost::lookup_named_param_def <
          internal_np::edge_is_constrained_t,
          NamedParameters,
          internal::Border_constraint_pmap<PolygonMesh, FaceRange, FIMap>
        > ::type ECMap;
    ECMap ecmap = (boost::is_same<ECMap, internal::Border_constraint_pmap<PolygonMesh, FaceRange, FIMap> >::value)
    ? choose_param(get_param(np, internal_np::edge_is_constrained),
                   internal::Border_constraint_pmap<PolygonMesh, FaceRange, FIMap>(pmesh, faces, fimap))
    : choose_param(get_param(np, internal_np::edge_is_constrained),
                   internal::Border_constraint_pmap<PolygonMesh, FaceRange, FIMap>());

    // nb_iterations
    unsigned int nb_iterations = choose_param(get_param(np, internal_np::number_of_iterations), 1);

    //use weighted angles - should always be true?
    bool use_weights = choose_param(get_param(np, internal_np::use_weights), true);

    // convergence precision
    double precision = choose_param(get_param(np, internal_np::number_of_iterations), 1);

    internal::Compatible_remesher<PolygonMesh, VertexPointMap, VCMap, ECMap, GeomTraits>
            remesher(pmesh, vpmap, vcmap, ecmap);

#ifdef CGAL_PMP_SMOOTHING_VERBOSE
    t.stop();
    std::cout << " done ("<< t.time() <<" sec)." << std::endl;
    std::cout << "Removing degenerate faces..." << std::endl;
    t.reset(); t.start();
#endif
    remesher.remove_degenerate_faces();

#ifdef CGAL_PMP_SMOOTHING_VERBOSE
    t.stop();
    std::cout << " done ("<< t.time() <<" sec)." << std::endl;
    std::cout << "Initializing..." << std::endl;
    t.reset(); t.start();
#endif
    remesher.init_smoothing(faces);

#ifdef CGAL_PMP_SMOOTHING_VERBOSE
    t.stop();
    std::cout << " done ("<< t.time() <<" sec)." << std::endl;
    std::cout << "#iter = " << nb_iterations << std::endl;
    std::cout << "Smoothing ..." << std::endl;
    t.reset(); t.start();
#endif

    for(unsigned int i=0; i<nb_iterations; ++i)
    {

#ifdef CGAL_PMP_SMOOTHING_VERBOSE
        std::cout << " * Iteration " << (i + 1) << " *" << std::endl;
#endif

        remesher.angle_relaxation(use_weights);
        remesher.project_to_surface();
    }


#ifdef CGAL_PMP_SMOOTHING_VERBOSE
    t.stop();
    std::cout << "Smoothing done in ";
    std::cout << t.time() << " sec." << std::endl;
#endif

}

template<typename PolygonMesh, typename NamedParameters>
void smooth_angles(PolygonMesh& pmesh, const NamedParameters& np)
{
    smooth_angles(faces(pmesh), pmesh, np);
}

template<typename PolygonMesh>
void smooth_angles(PolygonMesh& pmesh)
{
    smooth_angles(faces(pmesh), pmesh, parameters::all_default());
}

/*!
* \ingroup PMP_meshing_grp
* smoothes a triangulated region of a polygon mesh using area based criteria.
* This function tries to make the area distribution of triangle as uniform as possible
* by moving non constrained vertices.
* Optionally, the points are reprojected after each iteration.
* @todo make the reprojection optional.
*
* @tparam PolygonMesh model of `MutableFaceGraph`.
*         If `PolygonMesh` has an internal property map for `CGAL::face_index_t`,
*         and no `face_index_map` is given
*         as a named parameter, then the internal one should be initialized.
* @tparam FaceRange range of `boost::graph_traits<PolygonMesh>::%face_descriptor`,
          model of `Range`. Its iterator type is `ForwardIterator`.
* @tparam NamedParameters a sequence of \ref namedparameters
*
* @param pmesh a polygon mesh with triangulated surface patches to be smoothed.
* @param faces the range of triangular faces defining one or several surface patches to be smoothed.
* @param np optional sequence of \ref namedparameters among the ones listed below.
*
* \cgalNamedParamsBegin
*  \cgalParamBegin{geom_traits} a geometric traits class instance, model of `Kernel`.
*    Exact constructions kernels are not supported by this function.
*  \cgalParamEnd
*  \cgalParamBegin{vertex_point_map} the property map with the points associated
*    to the vertices of `pmesh`. Instance of a class model of `ReadWritePropertyMap`.
*  \cgalParamEnd
*  \cgalParamBegin{face_index_map} a property map containing the index of each face of `pmesh`.
*  \cgalParamEnd
*  \cgalParamBegin{number_of_iterations} the number of iterations for the
*    sequence of the smoothing iterations performed.
*  \cgalParamEnd
*  \cgalParamBegin{edge_is_constrained_map} a property map containing the
*    constrained-or-not status of each edge of `pmesh`. Vertices that belong to constrained
*    edges are not modified at all during smoothing.
*  \cgalParamEnd
*  \cgalParamBegin{vertex_is_constrained_map} a property map containing the
*    constrained-or-not status of each vertex of `pmesh`. A constrained vertex
*    cannot be modified at all during smoothing.
*  \cgalParamEnd
*  \cgalParamBegin{gradient_descent_precision} The precision to be met while the energy of each
*    triangle is being minimized. Corresponds to the
*    relative energy of triangles between iterations. Triangle energy is defined based on its area compared to
*    the average area of all triangles adjacent to the vertex that is being moved.  Defaults to 1e-5.
*  \cgalParamEnd
* \cgalNamedParamsEnd
*/
template<typename PolygonMesh, typename FaceRange, typename NamedParameters>
void smooth_areas(const FaceRange& faces, PolygonMesh& pmesh, const NamedParameters& np)
{
    using boost::choose_param;
    using boost::get_param;

#ifdef CGAL_PMP_SMOOTHING_VERBOSE
    CGAL::Timer t;
    std::cout << "Smoothing parameters...";
    std::cout.flush();
    t.start();
#endif

    //geom traits
    typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type GeomTraits;

    //vpmap
    typedef typename GetVertexPointMap<PolygonMesh, NamedParameters>::type VertexPointMap;
    VertexPointMap vpmap = choose_param(get_param(np, internal_np::vertex_point),
                                 get_property_map(CGAL::vertex_point, pmesh));

    //fimap
    typedef typename GetFaceIndexMap<PolygonMesh, NamedParameters>::type FIMap;
    FIMap fimap = choose_param(get_param(np, internal_np::face_index),
                             get_property_map(face_index, pmesh));

    //vcmap
    typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;
    typedef typename boost::lookup_named_param_def <
        internal_np::vertex_is_constrained_t,
        NamedParameters,
        internal::No_constraint_pmap<vertex_descriptor>//default
      > ::type VCMap;
    VCMap vcmap = choose_param(get_param(np, internal_np::vertex_is_constrained),
                               internal::No_constraint_pmap<vertex_descriptor>());

    //ecmap
    typedef typename boost::lookup_named_param_def <
          internal_np::edge_is_constrained_t,
          NamedParameters,
          internal::Border_constraint_pmap<PolygonMesh, FaceRange, FIMap>
        > ::type ECMap;
    ECMap ecmap = (boost::is_same<ECMap, internal::Border_constraint_pmap<PolygonMesh, FaceRange, FIMap> >::value)
    ? choose_param(get_param(np, internal_np::edge_is_constrained),
                   internal::Border_constraint_pmap<PolygonMesh, FaceRange, FIMap>(pmesh, faces, fimap))
    : choose_param(get_param(np, internal_np::edge_is_constrained),
                   internal::Border_constraint_pmap<PolygonMesh, FaceRange, FIMap>());

    //nb_iterations
    unsigned int nb_iterations = choose_param(get_param(np, internal_np::number_of_iterations), 1);

    //gradient descent precision
    double gd_precision = choose_param(get_param(np, internal_np::gradient_descent_precision), 1e-5);

    internal::Compatible_remesher<PolygonMesh, VertexPointMap, VCMap, ECMap, GeomTraits>
            remesher(pmesh, vpmap, vcmap, ecmap);

#ifdef CGAL_PMP_SMOOTHING_VERBOSE
    t.stop();
    std::cout << " done ("<< t.time() <<" sec)." << std::endl;
    std::cout << "Removing degenerate faces..." << std::endl;
    t.reset(); t.start();
#endif
    remesher.remove_degenerate_faces();

#ifdef CGAL_PMP_SMOOTHING_VERBOSE
    t.stop();
    std::cout << " done ("<< t.time() <<" sec)." << std::endl;
    std::cout << "Initializing..." << std::endl;
    t.reset(); t.start();
#endif
    remesher.init_smoothing(faces);

#ifdef CGAL_PMP_SMOOTHING_VERBOSE
    t.stop();
    std::cout << " done ("<< t.time() <<" sec)." << std::endl;
    std::cout << "#iter = " << nb_iterations << std::endl;
    std::cout << "Smoothing ..." << std::endl;
    t.reset(); t.start();
#endif

    for(unsigned int i=0; i<nb_iterations; ++i)
    {

#ifdef CGAL_PMP_SMOOTHING_VERBOSE
        std::cout << " * Iteration " << (i + 1) << " *" << std::endl;
#endif

        remesher.area_relaxation(gd_precision);
        remesher.project_to_surface();
    }

#ifdef CGAL_PMP_SMOOTHING_VERBOSE
    t.stop();
    std::cout << "Smoothing done in ";
    std::cout << t.time() << " sec." << std::endl;
#endif

}

template<typename PolygonMesh, typename NamedParameters>
void smooth_areas(PolygonMesh& pmesh, const NamedParameters& np)
{
    smooth_areas(faces(pmesh), pmesh, np);
}

template<typename PolygonMesh>
void smooth_areas(PolygonMesh& pmesh)
{
    smooth_areas(faces(pmesh), pmesh, parameters::all_default());
}

/*!
* \ingroup PMP_meshing_grp
* smoothes a triangulated region of a polygon mesh using area and angle based criteria.
* Each iteration alternates a call to `smooth_angles()`, `smooth_areas()`, and `smooth_angles()`.
* Convergence is realized based on the Hausdorff distance of the mesh between previous
* and current iteration or until a specified maximum number of iterations.
*
* @tparam PolygonMesh model of `MutableFaceGraph`.
*         If `PolygonMesh` has an internal property map for `CGAL::face_index_t`,
*         and no `face_index_map` is given
*         as a named parameter, then the internal one should be initialized.
* @tparam FaceRange range of `boost::graph_traits<PolygonMesh>::%face_descriptor`,
          model of `Range`. Its iterator type is `ForwardIterator`.
* @tparam NamedParameters a sequence of \ref namedparameters
*
* @param pmesh a polygon mesh with triangulated surface patches to be smoothed.
* @param faces the range of triangular faces defining one or several surface patches to be smoothed.
* @param np optional sequence of \ref namedparameters among the ones listed below.
*
* @todo remove Hausdorff and use ratio between max (or mean?) displacement.
* \cgalNamedParamsBegin
*  \cgalParamBegin{geom_traits} a geometric traits class instance, model of `Kernel`.
*    Exact constructions kernels are not supported by this function.
*  \cgalParamEnd
*  \cgalParamBegin{vertex_point_map} the property map with the points associated
*    to the vertices of `pmesh`. Instance of a class model of `ReadWritePropertyMap`.
*  \cgalParamEnd
*  \cgalParamBegin{face_index_map} a property map containing the index of each face of `pmesh`.
*  \cgalParamEnd
*  \cgalParamBegin{number_of_iterations} the maximum number of iterations for the
*    sequence of smoothing iterations performed.
*  \cgalParamEnd
*  \cgalParamBegin{edge_is_constrained_map} a property map containing the
*    constrained-or-not status of each edge of `pmesh`. Vertices that belong to constrained
*    edges are not modified at all during smoothing.
*  \cgalParamEnd
*  \cgalParamBegin{vertex_is_constrained_map} a property map containing the
*    constrained-or-not status of each vertex of `pmesh`. A constrained vertex
*    cannot be modified at all during smoothing.
*  \cgalParamEnd
*  \cgalParamBegin{use_weights} If `true`, small angles carry more weight than larger ones.
*  \cgalParamEnd
*  \cgalParamBegin{distance_precision} The Hausdorff distance between the current mesh and the mesh at the iteration.
*    Defaults to 0.01.
*  \cgalParamEnd
* \cgalNamedParamsEnd
*/
template<typename PolygonMesh, typename FaceRange, typename NamedParameters>
void compatible_smoothing(const FaceRange& faces, PolygonMesh& pmesh, const NamedParameters& np)
{
    using boost::choose_param;
    using boost::get_param;

#ifdef CGAL_PMP_SMOOTHING_VERBOSE
    CGAL::Timer t;
    std::cout << "Smoothing parameters...";
    std::cout.flush();
    t.start();
#endif

    // geom traits
    typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type GeomTraits;

    // vpmap
    typedef typename GetVertexPointMap<PolygonMesh, NamedParameters>::type VertexPointMap;
    VertexPointMap vpmap = choose_param(get_param(np, internal_np::vertex_point),
                                 get_property_map(CGAL::vertex_point, pmesh));

    // fimap
    typedef typename GetFaceIndexMap<PolygonMesh, NamedParameters>::type FIMap;
    FIMap fimap = choose_param(get_param(np, internal_np::face_index),
                             get_property_map(face_index, pmesh));

    // vcmap
    typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;
    typedef typename boost::lookup_named_param_def <
        internal_np::vertex_is_constrained_t,
        NamedParameters,
        internal::No_constraint_pmap<vertex_descriptor>//default
      > ::type VCMap;
    VCMap vcmap = choose_param(get_param(np, internal_np::vertex_is_constrained),
                               internal::No_constraint_pmap<vertex_descriptor>());

    // ecmap
    typedef typename boost::lookup_named_param_def <
          internal_np::edge_is_constrained_t,
          NamedParameters,
          internal::Border_constraint_pmap<PolygonMesh, FaceRange, FIMap>
        > ::type ECMap;
    ECMap ecmap = (boost::is_same<ECMap, internal::Border_constraint_pmap<PolygonMesh, FaceRange, FIMap> >::value)
    ? choose_param(get_param(np, internal_np::edge_is_constrained),
                   internal::Border_constraint_pmap<PolygonMesh, FaceRange, FIMap>(pmesh, faces, fimap))
    : choose_param(get_param(np, internal_np::edge_is_constrained),
                   internal::Border_constraint_pmap<PolygonMesh, FaceRange, FIMap>());

    // nb_iterations
    unsigned int nb_iterations = choose_param(get_param(np, internal_np::number_of_iterations), 20);

    // use weighted angles
    bool use_weights = choose_param(get_param(np, internal_np::use_weights), true);

    // gradient descent precision
    double gd_precision = choose_param(get_param(np, internal_np::gradient_descent_precision), 0.0001);

    // convergence precision
    double dist_precision = choose_param(get_param(np, internal_np::distance_precision), 0.01);

    internal::Compatible_remesher<PolygonMesh, VertexPointMap, VCMap, ECMap, GeomTraits>
            remesher(pmesh, vpmap, vcmap, ecmap);

#ifdef CGAL_PMP_SMOOTHING_VERBOSE
    t.stop();
    std::cout << " done ("<< t.time() <<" sec)." << std::endl;
    std::cout << "Removing degenerate faces..." << std::endl;
    t.reset(); t.start();
#endif
    remesher.remove_degenerate_faces();

#ifdef CGAL_PMP_SMOOTHING_VERBOSE
    t.stop();
    std::cout << " done ("<< t.time() <<" sec)." << std::endl;
    std::cout << "Initializing..." << std::endl;
    t.reset(); t.start();
#endif
    remesher.init_smoothing(faces);


#ifdef CGAL_PMP_SMOOTHING_VERBOSE
    t.stop();
    std::cout << " done ("<< t.time() <<" sec)." << std::endl;
    std::cout << "#iter = " << nb_iterations << std::endl;
    std::cout << "Smoothing ..." << std::endl;
    t.reset(); t.start();
#endif

    unsigned int count_iterations = 0;
    PolygonMesh previous_mesh;

    while(true)
    {
        previous_mesh = pmesh;
        remesher.angle_relaxation(use_weights);
        remesher.area_relaxation(gd_precision);
        remesher.angle_relaxation(use_weights);
        remesher.project_to_surface();

        double dist = approximate_Hausdorff_distance
            <TAG>(previous_mesh, pmesh, parameters::number_of_points_per_area_unit(1000));

        count_iterations++;

#ifdef CGAL_PMP_SMOOTHING_VERBOSE
        std::cout<<"iteration:"<<count_iterations<<std::endl;
#endif

        if(dist < dist_precision)
        {

#ifdef CGAL_PMP_SMOOTHING_VERBOSE
            t.stop();
            std::cout << "Smoothing done in ";
            std::cout << t.time() << " sec." << std::endl;
            std::cout << "Convergence to relative hausdorff distance has been achieved." << std::endl;
#endif
            break;
        }

        if(count_iterations == nb_iterations)
        {

#ifdef CGAL_PMP_SMOOTHING_VERBOSE
            t.stop();
            std::cout << "Smoothing done in ";
            std::cout << t.time() << " sec." << std::endl;
            std::cout << "Maximum number of iterations has been achieved." << std::endl;
#endif
            break;
        }
    }

}

template<typename PolygonMesh, typename NamedParameters>
void compatible_smoothing(PolygonMesh& pmesh, const NamedParameters& np)
{
    compatible_smoothing(faces(pmesh), pmesh, np);
}

template<typename PolygonMesh>
void compatible_smoothing(PolygonMesh& pmesh)
{
    compatible_smoothing(faces(pmesh), pmesh, parameters::all_default());
}




// not documented
template<typename PolygonMesh, typename GeomTraits>
void angles_evaluation(PolygonMesh& pmesh, const char* filename)
{
    internal::Quality_evaluator<PolygonMesh, GeomTraits> evaluator(pmesh);

    evaluator.gather_angles();
    evaluator.extract_angles(filename);
}

template<typename PolygonMesh, typename GeomTraits>
void areas_evaluation(PolygonMesh& pmesh, const char* filename)
{
    internal::Quality_evaluator<PolygonMesh, GeomTraits> evaluator(pmesh);

    evaluator.measure_areas();
    evaluator.extract_areas(filename);
}

template<typename PolygonMesh, typename GeomTraits>
void aspect_ratio_evaluation(PolygonMesh& pmesh, const char* filename)
{
    internal::Quality_evaluator<PolygonMesh, GeomTraits> evaluator(pmesh);

    evaluator.calc_aspect_ratios();
    evaluator.extract_aspect_ratios(filename);
}





} // namespace Polygon_mesh_processing
} // namespace CGAL











#endif // CGAL_POLYGON_MESH_PROCESSING_SMOOTHING_H
