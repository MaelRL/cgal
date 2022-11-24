/*!
\ingroup PkgIsosurfacing3Concepts
\cgalConcept

The concept `IsosurfacingDomainWithGradient` describes the set of requirements to be
fulfilled by any class used as input data for some isosurfacing algorithms.

\cgalHasModel `CGAL::Isosurfacing::Explicit_cartesian_grid_domain`
\cgalHasModel `CGAL::Isosurfacing::Implicit_cartesian_grid_domain`
\cgalHasModel `CGAL::Isosurfacing::Implicit_octree_domain`

*/

class IsosurfacingDomainWithGradient : public IsosurfacingDomain {
public:
    /// \name Types
    /// @{

    /*!
    The vector type.
    */
    typedef unspecified_type Vector;

    /// @}

    /// \name Operations
    /// The following member function must exist.
    /// @{

    /*!
    Returns the gradient at the position `p`
    */
    Vector gradient(const Point& p) const;

    /// @}
};
