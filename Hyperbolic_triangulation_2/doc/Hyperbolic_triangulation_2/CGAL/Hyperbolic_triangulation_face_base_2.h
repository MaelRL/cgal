// Copyright (c) 2016-2017 INRIA Nancy - Grand Est (France).
// All rights reserved.

namespace CGAL { 

/*!
\ingroup PkgHyperbolicTriangulation2VertexFaceClasses

The class `Hyperbolic_triangulation_face_base_2` is designed as one of the default models for the
concept `HyperbolicTriangulationFaceBase_2`.

\tparam Gt must be a model of `HyperbolicDelaunayTriangulationTraits_2`.
\tparam Fb must be a model of  `TriangulationDSFaceBase_2`. Defaults to `Triangulation_ds_face_base_2<>`.


\sa Hyperbolic_triangulation_face_base_with_info_2

\cgalModels HyperbolicTriangulationFaceBase_2
*/


template < typename Gt, typename Fb >
class Hyperbolic_triangulation_face_base_2 : public Fb {

public:

  /// \name Creation
  /// @{
    /*!
      Default constructor
    */
    Hyperbolic_triangulation_face_base_2();

    /*!
      Creates a face to which the vertices `v0, v1, v2` are incident.
    */
    Hyperbolic_triangulation_face_base_2(Vertex_handle v0, 
  			                                 Vertex_handle v1, 
  			                                 Vertex_handle v2);
      
    /*!
      Creates a face to which the vertices `v0, v1, v2` are incident, and the faces `n0, n1, n2` are neighbors.
    */
    Hyperbolic_triangulation_face_base_2(Vertex_handle v0, 
  			                                 Vertex_handle v1, 
  			                                 Vertex_handle v2,
  			                                 Face_handle n0, 
  			                                 Face_handle n1, 
  			                                 Face_handle n2);
  /// @}

};


} //namespace CGAL 

