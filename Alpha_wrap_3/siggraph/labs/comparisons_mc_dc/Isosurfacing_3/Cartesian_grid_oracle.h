#ifndef CGAL_CARTESIAN_GRID_ORACLE_H
#define CGAL_CARTESIAN_GRID_ORACLE_H

#include "Cartesian_grid_3.h"

#include <array>

namespace CGAL {

    template<class GeomTraits>
    class Cartesian_grid_oracle {
      public:
        typedef GeomTraits Geom_traits;
        typedef typename Geom_traits::FT FT;
        typedef typename Geom_traits::Point_3 Point_3;
        typedef typename Geom_traits::Vector_3 Vector_3;

      public:
        Cartesian_grid_oracle( const Cartesian_grid_3<Geom_traits>& grid ) : grid( &grid ) {}

        std::size_t size_x() const { return grid->xdim(); }
        std::size_t size_y() const { return grid->ydim(); }
        std::size_t size_z() const { return grid->zdim(); }

        std::size_t lex_index( const std::size_t& i, const std::size_t& j, const std::size_t& k ) const {
            return k * grid->zdim() * grid->ydim() + j * grid->xdim() + i;
        }

        /// <summary>
        /// compute unique edge global index.
        /// </summary>
        /// <param name="e">local edge index</param>
        /// <param name="i_idx">i-index of cell</param>
        /// <param name="j_idx">j-index of cell</param>
        /// <param name="k_idx">k-index of cell</param>
        /// <returns></returns>
        std::size_t e_glIndex( const std::size_t& e, const std::size_t& i_idx, const std::size_t& j_idx, const std::size_t& k_idx ) const {
            const unsigned long long gei_pattern_ = 670526590282893600ull;
            const size_t i                        = i_idx + (size_t)( ( gei_pattern_ >> 5 * e ) & 1 );            // global_edge_id[eg][0];
            const size_t j                        = j_idx + (size_t)( ( gei_pattern_ >> ( 5 * e + 1 ) ) & 1 );    // global_edge_id[eg][1];
            const size_t k                        = k_idx + (size_t)( ( gei_pattern_ >> ( 5 * e + 2 ) ) & 1 );    // global_edge_id[eg][2];
            const size_t offs                     = (size_t)( ( gei_pattern_ >> ( 5 * e + 3 ) ) & 3 );
            return ( 3 * lex_index( i, j, k ) + offs );
        }

        /// compute the gradient of the scalar field with central difference
        void gradient( std::array<Vector_3, 8>& n, const std::array<float, 8>& s, const std::size_t& i, const std::size_t& j,
                       const std::size_t& k ) const {
            auto index = []( const int dim, const int ii ) { return ( ii < 0 ? 0 : ii >= dim ? ( dim - 1 ) : ii ); };

            const FT& dx = grid->voxel_x();
            const FT& dy = grid->voxel_y();
            const FT& dz = grid->voxel_z();

            // vertex 0
            n[0] = { 0.5f * ( s[1] - value( index( grid->xdim(), i - 1 ), j, k ) ) / dx,
                     0.5f * ( s[2] - value( i, index( grid->ydim(), j - 1 ), k ) ) / dy,
                     0.5f * ( s[4] - value( i, j, index( grid->zdim(), k - 1 ) ) ) / dz };

            // vertex 1
            n[1] = { 0.5f * ( value( index( grid->xdim(), i + 2 ), j, k ) - s[0] ) / dx,
                     0.5f * ( s[3] - value( i + 1, index( grid->ydim(), j - 1 ), k ) ) / dy,
                     0.5f * ( s[5] - value( i + 1, j, index( grid->zdim(), k - 1 ) ) ) / dz };

            // vertex 2
            n[2] = { 0.5f * ( s[3] - value( index( grid->xdim(), i - 1 ), j + 1, k ) ) / dx,
                     0.5f * ( value( i, index( grid->ydim(), j + 2 ), k ) - s[0] ) / dy,
                     0.5f * ( s[6] - value( i, j + 1, index( grid->zdim(), k - 1 ) ) ) / dz };

            // vertex 3
            n[3] = { 0.5f * ( value( index( grid->xdim(), i + 2 ), j + 1, k ) - s[2] ) / dx,
                     0.5f * ( value( i + 1, index( grid->ydim(), j + 2 ), k ) - s[1] ) / dy,
                     0.5f * ( s[7] - value( i + 1, j + 1, index( grid->zdim(), k - 1 ) ) ) / dz };

            // vertex 4
            n[4] = { 0.5f * ( s[5] - value( index( grid->xdim(), i - 1 ), j, k + 1 ) ) / dx,
                     0.5f * ( s[6] - value( i, index( grid->ydim(), j - 1 ), k + 1 ) ) / dy,
                     0.5f * ( value( i, j, index( grid->zdim(), k + 2 ) ) - s[0] ) / dz };

            // vertex 5
            n[5] = { 0.5f * ( value( index( grid->xdim(), i + 2 ), j, k + 1 ) - s[4] ) / dx,
                     0.5f * ( s[7] - value( i + 1, index( grid->ydim(), j - 1 ), k + 1 ) ) / dy,
                     0.5f * ( value( i + 1, j, index( grid->zdim(), k + 2 ) ) - s[1] ) / dz };

            // vertex 6
            n[6] = { 0.5f * ( s[7] - value( index( grid->xdim(), i - 1 ), j + 1, k + 1 ) ) / dx,
                     0.5f * ( value( i, index( grid->ydim(), j + 2 ), k + 1 ) - s[4] ) / dy,
                     0.5f * ( value( i, j + 1, index( grid->zdim(), k + 2 ) ) - s[2] ) / dz };

            // vertex 7
            n[7] = { 0.5f * ( value( index( grid->xdim(), i + 2 ), j + 1, k + 1 ) - s[6] ) / dx,
                     0.5f * ( value( i + 1, index( grid->ydim(), j + 2 ), k + 1 ) - s[5] ) / dy,
                     0.5f * ( value( i + 1, j + 1, index( grid->zdim(), k + 2 ) ) - s[3] ) / dz };
        }

        Point_3 position( const std::size_t x, const std::size_t y, const std::size_t z ) const {
            const FT vx = grid->voxel_x();
            const FT vy = grid->voxel_y();
            const FT vz = grid->voxel_z();

            return Point_3( x * vx + grid->offset_x(), y * vy + grid->offset_y(), z * vz + grid->offset_z() );
        }

        FT value( const std::size_t x, const std::size_t y, const std::size_t z ) const { return grid->value( x, y, z ); }

      private:
        const Cartesian_grid_3<Geom_traits>* grid;
    };

}    // namespace CGAL

#endif    // CGAL_CARTESIAN_GRID_ORACLE_H
