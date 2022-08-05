#ifndef CGAL_DUAL_CONTOURING_3_H
#define CGAL_DUAL_CONTOURING_3_H

#include "Isosurfacing_3/internal/Marching_cubes_3_internal.h"
#include "util.h"

#include <Eigen/SVD>
#include <array>
#include <map>
#include <vector>

namespace CGAL {

    /// <summary>
    /// Compute vertex position for Dual Contouring
    /// </summary>
    /// <typeparam name="Domain_"></typeparam>
    /// <param name="domain"></param>
    /// <param name="iso_value"></param>
    /// <param name="i"></param>
    /// <param name="j"></param>
    /// <param name="k"></param>
    /// <returns> true, if there is a point in the cell</returns>
    template<class Domain_, bool use_bbox = false>
    bool get_vertex_position( const Domain_& domain, const typename Domain_::FT iso_value, const std::size_t& i, const std::size_t& j,
                              const std::size_t& k, typename Domain_::Point_3& point ) {
        typedef typename Domain_::Point_3 Point_3;
        typedef typename Domain_::Vector_3 Vector_3;

        std::array<float, 8> s;
        s[0] = domain.value( i + 0, j + 0, k + 0 );
        s[1] = domain.value( i + 1, j + 0, k + 0 );
        s[2] = domain.value( i + 0, j + 1, k + 0 );
        s[3] = domain.value( i + 1, j + 1, k + 0 );
        s[4] = domain.value( i + 0, j + 0, k + 1 );
        s[5] = domain.value( i + 1, j + 0, k + 1 );
        s[6] = domain.value( i + 0, j + 1, k + 1 );
        s[7] = domain.value( i + 1, j + 1, k + 1 );

        std::array<bool, 8> b;
        b[0] = s[0] <= iso_value;
        b[1] = s[1] <= iso_value;
        b[2] = s[2] <= iso_value;
        b[3] = s[3] <= iso_value;
        b[4] = s[4] <= iso_value;
        b[5] = s[5] <= iso_value;
        b[6] = s[6] <= iso_value;
        b[7] = s[7] <= iso_value;

        unsigned int cubeindex = 0;
        // set bit if corresponding corner is below iso
        cubeindex |= b[0] << 0;
        cubeindex |= b[1] << 1;
        cubeindex |= b[2] << 2;
        cubeindex |= b[3] << 3;
        cubeindex |= b[4] << 4;
        cubeindex |= b[5] << 5;
        cubeindex |= b[6] << 6;
        cubeindex |= b[7] << 7;

        if( cubeindex == 0 || cubeindex == 255 ) {
            return false;
        }

        std::array<Vector_3, 8> pos;
        pos[0] = domain.position( i + 0, j + 0, k + 0 ) - CGAL::ORIGIN;
        pos[1] = domain.position( i + 1, j + 0, k + 0 ) - CGAL::ORIGIN;
        pos[2] = domain.position( i + 0, j + 1, k + 0 ) - CGAL::ORIGIN;
        pos[3] = domain.position( i + 1, j + 1, k + 0 ) - CGAL::ORIGIN;
        pos[4] = domain.position( i + 0, j + 0, k + 1 ) - CGAL::ORIGIN;
        pos[5] = domain.position( i + 1, j + 0, k + 1 ) - CGAL::ORIGIN;
        pos[6] = domain.position( i + 0, j + 1, k + 1 ) - CGAL::ORIGIN;
        pos[7] = domain.position( i + 1, j + 1, k + 1 ) - CGAL::ORIGIN;

        point = CGAL::ORIGIN + ( pos[0] + 0.5 * ( pos[7] - pos[0] ) );    // set point to voxel center

        std::array<Vector_3, 8> normals;
        domain.gradient( normals, s, i, j, k );

        // compute edge intersections
        std::vector<Point_3> edge_intersections;
        std::vector<Vector_3> edge_intersection_normals;

        if( b[0] != b[1] ) {    // e0
            const FT u           = ( s[0] - iso_value ) / ( s[0] - s[1] );
            const Point_3 p_lerp = CGAL::ORIGIN + ( ( 1 - u ) * pos[0] + u * pos[1] );
            edge_intersections.push_back( p_lerp );
            const Vector_3 n_lerp = ( 1 - u ) * normals[0] + u * normals[1];
            edge_intersection_normals.push_back( n_lerp );
        }
        if( b[1] != b[3] ) {    // e1
            const FT u           = ( s[1] - iso_value ) / ( s[1] - s[3] );
            const Point_3 p_lerp = CGAL::ORIGIN + ( ( 1 - u ) * pos[1] + u * pos[3] );
            edge_intersections.push_back( p_lerp );
            const Vector_3 n_lerp = ( 1 - u ) * normals[1] + u * normals[3];
            edge_intersection_normals.push_back( n_lerp );
        }
        if( b[2] != b[3] ) {    // e2
            const FT u           = ( s[2] - iso_value ) / ( s[2] - s[3] );
            const Point_3 p_lerp = CGAL::ORIGIN + ( ( 1 - u ) * pos[2] + u * pos[3] );
            edge_intersections.push_back( p_lerp );
            const Vector_3 n_lerp = ( 1 - u ) * normals[2] + u * normals[3];
            edge_intersection_normals.push_back( n_lerp );
        }
        if( b[0] != b[2] ) {    // e3
            const FT u           = ( s[0] - iso_value ) / ( s[0] - s[2] );
            const Point_3 p_lerp = CGAL::ORIGIN + ( ( 1 - u ) * pos[0] + u * pos[2] );
            edge_intersections.push_back( p_lerp );
            const Vector_3 n_lerp = ( 1 - u ) * normals[0] + u * normals[2];
            edge_intersection_normals.push_back( n_lerp );
        }
        if( b[4] != b[5] ) {    // e4
            const FT u           = ( s[4] - iso_value ) / ( s[4] - s[5] );
            const Point_3 p_lerp = CGAL::ORIGIN + ( ( 1 - u ) * pos[4] + u * pos[5] );
            edge_intersections.push_back( p_lerp );
            const Vector_3 n_lerp = ( 1 - u ) * normals[4] + u * normals[5];
            edge_intersection_normals.push_back( n_lerp );
        }
        if( b[5] != b[7] ) {    // e5
            const FT u           = ( s[5] - iso_value ) / ( s[5] - s[7] );
            const Point_3 p_lerp = CGAL::ORIGIN + ( ( 1 - u ) * pos[5] + u * pos[7] );
            edge_intersections.push_back( p_lerp );
            const Vector_3 n_lerp = ( 1 - u ) * normals[5] + u * normals[7];
            edge_intersection_normals.push_back( n_lerp );
        }
        if( b[6] != b[7] ) {    // e6
            const FT u           = ( s[6] - iso_value ) / ( s[6] - s[7] );
            const Point_3 p_lerp = CGAL::ORIGIN + ( ( 1 - u ) * pos[6] + u * pos[7] );
            edge_intersections.push_back( p_lerp );
            const Vector_3 n_lerp = ( 1 - u ) * normals[6] + u * normals[7];
            edge_intersection_normals.push_back( n_lerp );
        }
        if( b[4] != b[6] ) {    // e7
            const FT u           = ( s[4] - iso_value ) / ( s[4] - s[6] );
            const Point_3 p_lerp = CGAL::ORIGIN + ( ( 1 - u ) * pos[4] + u * pos[6] );
            edge_intersections.push_back( p_lerp );
            const Vector_3 n_lerp = ( 1 - u ) * normals[4] + u * normals[6];
            edge_intersection_normals.push_back( n_lerp );
        }
        if( b[0] != b[4] ) {    // e8
            const FT u           = ( s[0] - iso_value ) / ( s[0] - s[4] );
            const Point_3 p_lerp = CGAL::ORIGIN + ( ( 1 - u ) * pos[0] + u * pos[4] );
            edge_intersections.push_back( p_lerp );
            const Vector_3 n_lerp = ( 1 - u ) * normals[0] + u * normals[4];
            edge_intersection_normals.push_back( n_lerp );
        }
        if( b[1] != b[5] ) {    // e9
            const FT u           = ( s[1] - iso_value ) / ( s[1] - s[5] );
            const Point_3 p_lerp = CGAL::ORIGIN + ( ( 1 - u ) * pos[1] + u * pos[5] );
            edge_intersections.push_back( p_lerp );
            const Vector_3 n_lerp = ( 1 - u ) * normals[1] + u * normals[5];
            edge_intersection_normals.push_back( n_lerp );
        }
        if( b[3] != b[7] ) {    // e10
            const FT u           = ( s[3] - iso_value ) / ( s[3] - s[7] );
            const Point_3 p_lerp = CGAL::ORIGIN + ( ( 1 - u ) * pos[3] + u * pos[7] );
            edge_intersections.push_back( p_lerp );
            const Vector_3 n_lerp = ( 1 - u ) * normals[3] + u * normals[7];
            edge_intersection_normals.push_back( n_lerp );
        }
        if( b[2] != b[6] ) {    // e11
            const FT u           = ( s[2] - iso_value ) / ( s[2] - s[6] );
            const Point_3 p_lerp = CGAL::ORIGIN + ( ( 1 - u ) * pos[2] + u * pos[6] );
            edge_intersections.push_back( p_lerp );
            const Vector_3 n_lerp = ( 1 - u ) * normals[2] + u * normals[6];
            edge_intersection_normals.push_back( n_lerp );
        }

        // SVD QEM
        if( true ) {
            Eigen::Matrix3d A;
            A.setZero();
            Eigen::Vector3d b;
            b.setZero();
            for( int i = 0; i < edge_intersections.size(); ++i ) {
                Eigen::Vector3d n_k = { edge_intersection_normals[i].x(), edge_intersection_normals[i].y(), edge_intersection_normals[i].z() };
                Eigen::Vector3d p_k = { edge_intersections[i].x(), edge_intersections[i].y(), edge_intersections[i].z() };
                double d_k          = n_k.transpose() * p_k;

                Eigen::Matrix3d A_k = n_k * n_k.transpose();
                Eigen::Vector3d b_k = d_k * n_k;
                A += A_k;
                b += b_k;
            }

            Eigen::JacobiSVD<Eigen::MatrixXd> svd( A, Eigen::ComputeThinU | Eigen::ComputeThinV );
            // set threshold as in Peter Lindstrom's paper, "Out-of-Core
            // Simplification of Large Polygonal Models"
            svd.setThreshold( 1e-3 );

            // Init x hat
            Eigen::Vector3d x_hat;
            x_hat << point.x(), point.y(), point.z();

            // Lindstrom formula for QEM new position for singular matrices
            Eigen::VectorXd v_svd = x_hat + svd.solve( b - A * x_hat );

            point = Point_3( v_svd[0], v_svd[1], v_svd[2] );
        }

        // bbox
        if constexpr( use_bbox ) {
            CGAL::Bbox_3 bbox = ( CGAL::ORIGIN + pos[0] ).bbox() + ( CGAL::ORIGIN + pos[7] ).bbox();

            FT x  = std::min<FT>( std::max<FT>( point.x(), bbox.xmin() ), bbox.xmax() );
            FT y  = std::min<FT>( std::max<FT>( point.y(), bbox.ymin() ), bbox.ymax() );
            FT z  = std::min<FT>( std::max<FT>( point.z(), bbox.zmin() ), bbox.zmax() );
            point = Point_3( x, y, z );
        }

        return true;
    }

    template<class Domain_, class PointRange, class PolygonRange>
    void make_quad_mesh_using_dual_contouring( const Domain_& domain, const typename Domain_::FT iso_value, PointRange& points,
                                               PolygonRange& polygons, const bool use_bbox = false ) {
        typedef typename Domain_::FT FT;
        typedef typename Domain_::Point_3 Point_3;

        ScopeTimer timer;

        // compute dc-vertices
        std::map<size_t, size_t> map_voxel_to_point_id;
        std::map<size_t, Point_3> map_voxel_to_point;
        size_t points_counter = 0;
        std::map<size_t, std::array<size_t, 4>> quads;

        const std::size_t size_k = domain.size_x();
        const std::size_t size_j = domain.size_y();
        const std::size_t size_i = domain.size_z();

        // save all points
        if( use_bbox ) {
            for( int k = 0; k < size_k - 1; k++ ) {
                for( int j = 0; j < size_j - 1; j++ ) {
                    for( int i = 0; i < size_i - 1; i++ ) {
                        const size_t voxel_id = domain.lex_index( i, j, k );
                        Point_3 p;
                        if( get_vertex_position<Domain_, true>( domain, iso_value, i, j, k, p ) ) {
                            map_voxel_to_point[voxel_id]    = p;
                            map_voxel_to_point_id[voxel_id] = points_counter++;
                        }
                    }
                }
            }
        } else {
            for( int k = 0; k < size_k - 1; k++ ) {
                for( int j = 0; j < size_j - 1; j++ ) {
                    for( int i = 0; i < size_i - 1; i++ ) {
                        const size_t voxel_id = domain.lex_index( i, j, k );
                        Point_3 p;
                        if( get_vertex_position( domain, iso_value, i, j, k, p ) ) {
                            map_voxel_to_point[voxel_id]    = p;
                            map_voxel_to_point_id[voxel_id] = points_counter++;
                        }
                    }
                }
            }
        }

        // save all quads
        for( int k = 0; k < size_k - 1; k++ ) {
            for( int j = 0; j < size_j - 1; j++ ) {
                for( int i = 0; i < size_i - 1; i++ ) {
                    const double s0 = domain.value( i + 0, j + 0, k + 0 );
                    const double s1 = domain.value( i + 1, j + 0, k + 0 );
                    const double s2 = domain.value( i + 0, j + 1, k + 0 );
                    const double s3 = domain.value( i + 1, j + 1, k + 0 );
                    const double s4 = domain.value( i + 0, j + 0, k + 1 );
                    const double s5 = domain.value( i + 1, j + 0, k + 1 );
                    const double s6 = domain.value( i + 0, j + 1, k + 1 );
                    const double s7 = domain.value( i + 1, j + 1, k + 1 );

                    unsigned int cubeindex = 0;
                    // set bit if corresponding corner is below iso
                    cubeindex |= ( s0 <= iso_value ) << 0;
                    cubeindex |= ( s1 <= iso_value ) << 1;
                    cubeindex |= ( s2 <= iso_value ) << 2;
                    cubeindex |= ( s3 <= iso_value ) << 3;
                    cubeindex |= ( s4 <= iso_value ) << 4;
                    cubeindex |= ( s5 <= iso_value ) << 5;
                    cubeindex |= ( s6 <= iso_value ) << 6;
                    cubeindex |= ( s7 <= iso_value ) << 7;

                    if( cubeindex == 0 || cubeindex == 255 ) {
                        continue;
                    }

                    // e0 x
                    if( j > 0 && k > 0 ) {
                        if( s0 > iso_value && s1 <= iso_value ) {
                            const size_t a     = domain.lex_index( i, j, k );
                            const size_t b     = domain.lex_index( i, j, k - 1 );
                            const size_t c     = domain.lex_index( i, j - 1, k - 1 );
                            const size_t d     = domain.lex_index( i, j - 1, k );
                            const size_t e_idx = domain.e_glIndex( 0, i, j, k );
                            quads[e_idx]       = { a, b, c, d };
                        } else if( s0 <= iso_value && s1 > iso_value ) {
                            const size_t a     = domain.lex_index( i, j, k );
                            const size_t b     = domain.lex_index( i, j - 1, k );
                            const size_t c     = domain.lex_index( i, j - 1, k - 1 );
                            const size_t d     = domain.lex_index( i, j, k - 1 );
                            const size_t e_idx = domain.e_glIndex( 0, i, j, k );
                            quads[e_idx]       = { a, b, c, d };
                        }
                    }
                    // e3 y
                    if( i > 0 && k > 0 ) {
                        if( s0 > iso_value && s2 <= iso_value ) {
                            const size_t a     = domain.lex_index( i, j, k );
                            const size_t b     = domain.lex_index( i - 1, j, k );
                            const size_t c     = domain.lex_index( i - 1, j, k - 1 );
                            const size_t d     = domain.lex_index( i, j, k - 1 );
                            const size_t e_idx = domain.e_glIndex( 3, i, j, k );
                            quads[e_idx]       = { a, b, c, d };
                        } else if( s0 <= iso_value && s2 > iso_value ) {
                            const size_t a     = domain.lex_index( i, j, k );
                            const size_t b     = domain.lex_index( i, j, k - 1 );
                            const size_t c     = domain.lex_index( i - 1, j, k - 1 );
                            const size_t d     = domain.lex_index( i - 1, j, k );
                            const size_t e_idx = domain.e_glIndex( 3, i, j, k );
                            quads[e_idx]       = { a, b, c, d };
                        }
                    }
                    // e8 z
                    if( i > 0 && j > 0 ) {
                        if( s0 > iso_value && s4 <= iso_value ) {
                            const size_t a     = domain.lex_index( i, j, k );
                            const size_t b     = domain.lex_index( i, j - 1, k );
                            const size_t c     = domain.lex_index( i - 1, j - 1, k );
                            const size_t d     = domain.lex_index( i - 1, j, k );
                            const size_t e_idx = domain.e_glIndex( 8, i, j, k );
                            quads[e_idx]       = { a, b, c, d };
                        } else if( s0 <= iso_value && s4 > iso_value ) {
                            const size_t a     = domain.lex_index( i, j, k );
                            const size_t b     = domain.lex_index( i - 1, j, k );
                            const size_t c     = domain.lex_index( i - 1, j - 1, k );
                            const size_t d     = domain.lex_index( i, j - 1, k );
                            const size_t e_idx = domain.e_glIndex( 8, i, j, k );
                            quads[e_idx]       = { a, b, c, d };
                        }
                    }
                }
            }
        }

        // write points and quads in ranges
        points.resize( points_counter );
        for( const auto& vtop: map_voxel_to_point ) {
            points[map_voxel_to_point_id[vtop.first]] = vtop.second;
        }

        polygons.reserve( quads.size() );
        for( const auto& q: quads ) {
            std::vector<std::size_t> vertex_ids;
            for( const auto& v_id: q.second ) {
                vertex_ids.push_back( map_voxel_to_point_id[v_id] );
            }
            polygons.push_back( vertex_ids );
        }
    }

}    // namespace CGAL

#endif    // CGAL_DUAL_CONTOURING_3_H
