#ifndef CGAL_MARCHING_CUBES_3_H
#define CGAL_MARCHING_CUBES_3_H

#include "Isosurfacing_3/internal/Marching_cubes_3_internal.h"
#include "util.h"

#include <array>
#include <mutex>

namespace CGAL {

    template<class Domain_, class PointRange, class PolygonRange>
    void make_triangle_mesh_using_marching_cubes( const Domain_& domain, const typename Domain_::FT iso_value, PointRange& points,
                                                  PolygonRange& polygons ) {
        ScopeTimer timer;

        std::mutex mutex;

        const std::size_t size_k = domain.size_x();
        const std::size_t size_j = domain.size_y();
        const std::size_t size_i = domain.size_z();

        // TODO: look at polygon mesh processing for tbb (also linking)

        #pragma omp parallel for
        for( int k = 0; k < size_k - 1; k++ ) {
            for( int j = 0; j < size_j - 1; j++ ) {
                for( int i = 0; i < size_i - 1; i++ ) {
                    internal::Marching_cubes_3::marching_cubes_cell( i, j, k, domain, iso_value, points, polygons, mutex );
                }
            }
        }
    }

}    // namespace CGAL

#endif    // CGAL_MARCHING_CUBES_3_H
