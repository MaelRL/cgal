#ifndef CGAL_MARCHING_CUBES_3_INTERNAL_MARCHING_CUBES_3_H
#define CGAL_MARCHING_CUBES_3_INTERNAL_MARCHING_CUBES_3_H

#include "Tables.h"

#include <mutex>

namespace CGAL {

    namespace internal {
        namespace Marching_cubes_3 {

            template<class Point_3, typename FT>
            Point_3 vertex_interpolation( const Point_3& p0, const Point_3& p1, const FT d0, const FT d1, const FT iso_value ) {
                FT mu = -1.0;

                // don't divide by 0
                if( abs( d1 - d0 ) < 0.000001 ) {
                    mu = 0.5;    // if both points have the same value, assume isolevel is in the middle
                } else {
                    mu = ( iso_value - d0 ) / ( d1 - d0 );
                }

                if( mu < 0.0 || mu > 1.0 ) {
                    printf( "ERROR: isolevel is not between points\n" );    // TODO: error handling
                }

                // linear interpolation
                return Point_3( p1.x() * mu + p0.x() * ( 1 - mu ), p1.y() * mu + p0.y() * ( 1 - mu ), p1.z() * mu + p0.z() * ( 1 - mu ) );
            }

            template<class Domain_, class PointRange, class PolygonRange>
            void marching_cubes_cell( const std::size_t x, const std::size_t y, const std::size_t z, const Domain_& domain,
                                      const typename Domain_::FT iso_value, PointRange& points, PolygonRange& polygons, std::mutex& mutex ) {
                typedef std::array<std::size_t, 3> Idx_3;
                typedef typename Domain_::FT FT;
                typedef typename Domain_::Point_3 Point_3;

                const Idx_3 idx0 = { x + 0, y + 1, z + 0 };
                const Idx_3 idx1 = { x + 1, y + 1, z + 0 };
                const Idx_3 idx2 = { x + 1, y + 0, z + 0 };
                const Idx_3 idx3 = { x + 0, y + 0, z + 0 };
                const Idx_3 idx4 = { x + 0, y + 1, z + 1 };
                const Idx_3 idx5 = { x + 1, y + 1, z + 1 };
                const Idx_3 idx6 = { x + 1, y + 0, z + 1 };
                const Idx_3 idx7 = { x + 0, y + 0, z + 1 };

                const Point_3 pos0 = domain.position( idx0[0], idx0[1], idx0[2] );
                const Point_3 pos1 = domain.position( idx1[0], idx1[1], idx1[2] );
                const Point_3 pos2 = domain.position( idx2[0], idx2[1], idx2[2] );
                const Point_3 pos3 = domain.position( idx3[0], idx3[1], idx3[2] );
                const Point_3 pos4 = domain.position( idx4[0], idx4[1], idx4[2] );
                const Point_3 pos5 = domain.position( idx5[0], idx5[1], idx5[2] );
                const Point_3 pos6 = domain.position( idx6[0], idx6[1], idx6[2] );
                const Point_3 pos7 = domain.position( idx7[0], idx7[1], idx7[2] );

                const FT dist0 = domain.value( idx0[0], idx0[1], idx0[2] );
                const FT dist1 = domain.value( idx1[0], idx1[1], idx1[2] );
                const FT dist2 = domain.value( idx2[0], idx2[1], idx2[2] );
                const FT dist3 = domain.value( idx3[0], idx3[1], idx3[2] );
                const FT dist4 = domain.value( idx4[0], idx4[1], idx4[2] );
                const FT dist5 = domain.value( idx5[0], idx5[1], idx5[2] );
                const FT dist6 = domain.value( idx6[0], idx6[1], idx6[2] );
                const FT dist7 = domain.value( idx7[0], idx7[1], idx7[2] );

                unsigned int cubeindex = 0;
                // set bit if corresponding corner is below iso
                cubeindex |= ( dist0 < iso_value ) << 0;
                cubeindex |= ( dist1 < iso_value ) << 1;
                cubeindex |= ( dist2 < iso_value ) << 2;
                cubeindex |= ( dist3 < iso_value ) << 3;
                cubeindex |= ( dist4 < iso_value ) << 4;
                cubeindex |= ( dist5 < iso_value ) << 5;
                cubeindex |= ( dist6 < iso_value ) << 6;
                cubeindex |= ( dist7 < iso_value ) << 7;

                Point_3 vertlist[12];
                // bitmap of edges that intersect the iso-plane(s)
                const int edges = edgeTable[cubeindex];

                if( edges == 0 ) {
                    return;
                }

                // interpolation of vertices incident to the edge, according to diagram
                if( edges & ( 1 << 0 ) ) {
                    vertlist[0] = vertex_interpolation( pos0, pos1, dist0, dist1, iso_value );
                }
                if( edges & ( 1 << 1 ) ) {
                    vertlist[1] = vertex_interpolation( pos1, pos2, dist1, dist2, iso_value );
                }
                if( edges & ( 1 << 2 ) ) {
                    vertlist[2] = vertex_interpolation( pos2, pos3, dist2, dist3, iso_value );
                }
                if( edges & ( 1 << 3 ) ) {
                    vertlist[3] = vertex_interpolation( pos3, pos0, dist3, dist0, iso_value );
                }
                if( edges & ( 1 << 4 ) ) {
                    vertlist[4] = vertex_interpolation( pos4, pos5, dist4, dist5, iso_value );
                }
                if( edges & ( 1 << 5 ) ) {
                    vertlist[5] = vertex_interpolation( pos5, pos6, dist5, dist6, iso_value );
                }
                if( edges & ( 1 << 6 ) ) {
                    vertlist[6] = vertex_interpolation( pos6, pos7, dist6, dist7, iso_value );
                }
                if( edges & ( 1 << 7 ) ) {
                    vertlist[7] = vertex_interpolation( pos7, pos4, dist7, dist4, iso_value );
                }
                if( edges & ( 1 << 8 ) ) {
                    vertlist[8] = vertex_interpolation( pos0, pos4, dist0, dist4, iso_value );
                }
                if( edges & ( 1 << 9 ) ) {
                    vertlist[9] = vertex_interpolation( pos1, pos5, dist1, dist5, iso_value );
                }
                if( edges & ( 1 << 10 ) ) {
                    vertlist[10] = vertex_interpolation( pos2, pos6, dist2, dist6, iso_value );
                }
                if( edges & ( 1 << 11 ) ) {
                    vertlist[11] = vertex_interpolation( pos3, pos7, dist3, dist7, iso_value );
                }

                std::lock_guard<std::mutex> lock( mutex );

                for( int i = 0; triTable[cubeindex][i] != -1; i += 3 ) {
                    const Point_3& p0 = vertlist[triTable[cubeindex][i + 0]];
                    const Point_3& p1 = vertlist[triTable[cubeindex][i + 1]];
                    const Point_3& p2 = vertlist[triTable[cubeindex][i + 2]];

                    const std::size_t p0_idx = points.size();    // TODO: not allowed

                    points.push_back( p0 );
                    points.push_back( p1 );
                    points.push_back( p2 );

                    polygons.push_back( {} );
                    auto& triangle = polygons.back();

                    triangle.push_back( p0_idx + 0 );
                    triangle.push_back( p0_idx + 1 );
                    triangle.push_back( p0_idx + 2 );
                }
            }

        }    // namespace Marching_cubes_3
    }        // namespace internal

}    // namespace CGAL

#endif    // CGAL_MARCHING_CUBES_3_INTERNAL_MARCHING_CUBES_3_H
