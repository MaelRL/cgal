#include <CGAL/Starset.h>
#include <CGAL/Stretched_Delaunay_3.h>

#include <Metric_field/Euclidean_metric_field.h>
#include <Metric_field/Custom_metric_field.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <Eigen/Dense>

#include <omp.h>
#include <boost/array.hpp>

#include <iostream>
#include <vector>

using namespace CGAL::Anisotropic_mesh_3;

typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
typedef typename K::FT                                       FT;
typedef Stretched_Delaunay_3<K>                              Star;
typedef typename Star::Point_3                               Point_3;
typedef typename Star::Metric                                Metric;
typedef typename Star::Traits                                Traits;

typedef boost::array<std::size_t, 3>                         Triangle;
typedef boost::array<std::size_t, 4>                         Tet;
typedef std::set<std::size_t>                                Simplex;

typedef typename Eigen::Matrix<double, 3, 1>                 Vector3d;

typedef typename CGAL::Anisotropic_mesh_3::Euclidean_metric_field<K>* MF;
//typedef typename CGAL::Anisotropic_mesh_3::Custom_metric_field<K>* MF;

#define FILTER_SEEDS_OUTSIDE_GRID
#define USE_RECURSIVE_UPDATES

#define TMP_REFINEMENT_UGLY_HACK

// using a lot of global variables, it's ugly but passing them by function is tedious
std::size_t vertices_nv = 10;

// grid related stuff
Point_3 center(0.,0.,0.);
const FT grid_side = 4.0;
FT offset_x = center.x() - grid_side/2.; // offset is the bottom left point
FT offset_y = center.y() - grid_side/2.;
FT offset_z = center.z() - grid_side/2.;
FT n = 50.; // number of points per side of the grid
FT sq_n = n*n;
FT step = grid_side / n;

// the metric field and the seeds
MF mf;
std::vector<Point_3> seeds;
std::vector<Metric> seeds_m;

// how big of an update the new value needs to be
FT recursive_tolerance = 1e-4;

#ifdef TMP_REFINEMENT_UGLY_HACK
FT best_ref_x = 1e30, best_ref_y = 1e30, best_ref_z = 1e30;
#endif

enum FMM_state
{
  CHANGED = 0, // changed implies known. It is used to not add elements that are already in the queue
  KNOWN,
  TRIAL,
  FAR
};

template<typename Gp>
struct Grid_point_comparer
{
  bool operator()(Gp const * const gp1, Gp const * const gp2)
  {
    return gp1->distance_to_closest_seed < gp2->distance_to_closest_seed;
  }
};

int insert_new_seed(const FT x, const FT y, const FT z)
{
#ifdef FILTER_SEEDS_OUTSIDE_GRID
    if(x < offset_x || x > center.x() + grid_side/2. ||
       y < offset_y || y > center.y() + grid_side/2. ||
       z < offset_z || z > center.z() + grid_side/2.)
    {
      std::cout << "filtered : " << x << " " << y << " " << z << std::endl;
      return seeds.size();
    }
#endif
  seeds.push_back(Point_3(x, y, z));
  seeds_m.push_back(mf->compute_metric(seeds.back()));
  return seeds.size();
}

int build_seeds()
{
  std::ifstream in("bambimboum.mesh");
  std::string word;
  std::size_t useless, nv, dim;
  FT r_x, r_y, r_z;

  in >> word >> useless; //MeshVersionFormatted i
  in >> word >> dim; //Dimension d
  in >> word >> nv;
  std::cout << "file has : " << nv << " pts" << std::endl;
  assert(dim == 3);

  std::size_t min_nv = (std::min)(nv, vertices_nv);

  seeds.reserve(min_nv);
  seeds_m.reserve(min_nv);

  for(std::size_t i=0; i<nv; ++i)
  {
    in >> r_x >> r_y >> r_z >> useless;
//    insert_new_seed(offset_x, offset_y, offset_z);
//    insert_new_seed(1.6, offset_y, offset_z);
//    insert_new_seed(1.6, 1.6, 1.6);
//    insert_new_seed(offset_x, offset_y+grid_side/2., offset_z);
//    insert_new_seed(offset_x, offset_y, offset_z+grid_side/2.);
    insert_new_seed(r_x, r_y, r_z);

    if(seeds.size() >= vertices_nv)
      break;
  }
  std::cout << "seeds: " << seeds.size() << std::endl;
  return seeds.size();
}

struct Grid_point
{
  typedef boost::array<Grid_point*, 6> Neighbors;

  Point_3 point;
  std::size_t index;
  mutable FT distance_to_closest_seed;
  mutable std::size_t closest_seed_id;
  mutable  FMM_state state;
  Neighbors neighbors; // 6 neighbors ORDERED (above, left, back, right, front, below)
                                          // if no neighbor (borders for example) then NULL
  Metric metric;

  const Grid_point* ancestor; // used to debug & draw nice things

  bool compute_closest_seed_1D(const Grid_point* gp)
  {
//    std::cout << "compute closest 1D " << index << " to " << gp->index << std::endl;
//    std::cout << "gp distance to closest is : " << gp->distance_to_closest_seed << std::endl;

    CGAL_assertion(gp->state == KNOWN || gp->state == CHANGED);
    bool changed = false;

    Vector3d v;
    v(0) = gp->point.x() - point.x();
    v(1) = gp->point.y() - point.y();
    v(2) = gp->point.z() - point.z();
    FT neighbor_d = std::sqrt(v.transpose()*metric.get_mat()*v);
    FT dcs_at_gp = gp->distance_to_closest_seed;
    FT d = dcs_at_gp + neighbor_d;
    if(distance_to_closest_seed - d > recursive_tolerance * d)
    {
      changed = true;
      distance_to_closest_seed = d;
      closest_seed_id = gp->closest_seed_id;
      ancestor = gp;
//      std::cout << "new closest at : " << distance_to_closest_seed << std::endl;
    }
    return changed;
  }

  FT eval_at_p(FT p, const FT a, const FT b, const FT c, const FT d, const FT e) const
  {
    CGAL_assertion(c*p*p + d*p + e >= 0);
    return p*a+(1.-p)*b + std::sqrt(c*p*p + d*p + e);
  }

  FT compute_min_distance_2D(const Grid_point* gp, const Grid_point* gq,
                             FT& p) const
  {
    CGAL_assertion((gp->state == KNOWN || gp->state == CHANGED) &&
                   (gq->state == KNOWN || gq->state == CHANGED) );

    // in 2D, we can solve the minimization analytically.
    // We're solving :
    // min_{p\in [0,1]} f(p) = F(x)*p + F(y)*(1-p) + sqrt(v(p)^t * M * v(p))
    // local minimum is given by f'(p_0) = 0 _but_ we have to compare f(p_0)
    // with f(0) and f(1) as we're only interested in solutions that live in [0,1]
    // If there's no solution in [0,1], then there's no (local) minimum in [0,1]
    // and the minimum is obtained in 0 or in 1.

    const FT dcs_at_gp = gp->distance_to_closest_seed;
    const FT dcs_at_gq = gq->distance_to_closest_seed;
    Vector3d vp, vq;
    vp(0) = gp->point.x() - point.x();
    vp(1) = gp->point.y() - point.y();
    vp(2) = gp->point.z() - point.z();
    vq(0) = gq->point.x() - point.x();
    vq(1) = gq->point.y() - point.y();
    vq(2) = gq->point.z() - point.z();

    Eigen::Matrix3d m = metric.get_mat();
    const FT vpmvp = vp.transpose() * m * vp;
    const FT vpmvq = vp.transpose() * m * vq;
    const FT vqmvq = vq.transpose() * m * vq;

    // expressing f(p) as ap+(1-p)b+sqrt(cp^2+dp+e)
    // the tiny trick is to write f'(p) = 0 and then square the expression to get rid of sqrt

    const FT a = dcs_at_gp;
    const FT b = dcs_at_gq;
    const FT c = vpmvp + vqmvq - 2.*vpmvq;
    const FT d = 2.*(vpmvq - vqmvq);
    const FT e = vqmvq;

    // Write f'(p) as Ap^2 + Bp + C
    const FT A = 4*(a-b)*(a-b)*c-4*c*c;
    const FT B = 4*(a-b)*(a-b)*d-4*c*d;
    const FT C = 4*(a-b)*(a-b)*e-d*d;
    FT delta = B*B-4*A*C;

//    std::cout << "a: " << a << std::endl;
//    std::cout << "b: " << b << std::endl;
//    std::cout << "c: " << c << std::endl;
//    std::cout << "d: " << d << std::endl;
//    std::cout << "e: " << e << std::endl;
//    std::cout << "A: " << A << std::endl;
//    std::cout << "B: " << B << std::endl;
//    std::cout << "C: " << C << std::endl;
//    std::cout << "delta: " << delta << std::endl;

    FT min_id = 0.;
    FT min = eval_at_p(0., a, b, c, d, e);

    FT d1 = eval_at_p(1., a, b, c, d, e);
    if(min > d1)
    {
      min_id = 1.;
      min = d1;
    }

    // get solutions
    CGAL_assertion(A != 0. || B != 0.);
    if(A == 0.)
    {
      p = -C/B;
      FT dp = eval_at_p(p, a, b, c, d, e);
      if(min > d && p >= 0 && p <= 1.)
      {
        min_id = p;
        min = dp;
      }
    }
    else if(delta >= 0)
    {
      FT sqrt_delta = std::sqrt(delta);
      FT denom = 1. / (2.*A);
      double p1 = (-B + sqrt_delta) * denom;
      double p2 = (-B - sqrt_delta) * denom;

      FT dp1 = eval_at_p(p1, a, b, c, d, e);
      FT dp2 = eval_at_p(p2, a, b, c, d, e);

//      std::cout << "p1 : " << p1 << " dp1: " << dp1 << std::endl;
//      std::cout << "p2 : " << p2 << " dp2: " << dp2 << std::endl;

      if(min > dp1 && p1 >= 0 && p1 <= 1.)
      {
        min_id = p1;
        min = dp1;
      }
      if(min > dp2 && p2 >= 0 && p2 <= 1.)
      {
        min_id = p2;
        min = dp2;
      }
    }

    p = min_id;
//    std::cout << "final min at : " << p << " d: " << min << std::endl;
    return min;
  }

  bool compute_closest_seed_2D(const Grid_point* gp, const Grid_point* gq)
  {
//    std::cout << "compute closest 2D " << index << " to " << gp->index << " & " << gq->index << std::endl;

    CGAL_assertion((gp->state == KNOWN || gp->state == CHANGED) &&
                   (gq->state == KNOWN || gq->state == CHANGED) );
    bool changed = false;

    FT p = 0.;
    const FT d = compute_min_distance_2D(gp, gq, p);

    // tolerance to ignore unsignificant changes
    if(distance_to_closest_seed - d > recursive_tolerance *d)
    {
      changed = true;
      distance_to_closest_seed = d;
      FT gp_d = gp->distance_to_closest_seed;
      FT gq_d = gq->distance_to_closest_seed;
      const Grid_point* best_p = (gp_d < gq_d) ? gp : gq; // fixme this works but really shouldn't...
      //const Grid_point* best_p = (p > 0.5) ? gp : gq; // this should be the correct one but it doesn't work ??
      closest_seed_id = best_p->closest_seed_id;
      ancestor = best_p; // a bit ugly, create the point p ?
    }
    return changed;
  }

  FT compute_closest_seed_3D(const Grid_point* gp,
                             const Grid_point* gq,
                             const Grid_point* gr)
  {
    std::cout << "compute closest seed 3D" << std::endl;
    CGAL_assertion((gp->state == KNOWN || gp->state == CHANGED) &&
                   (gq->state == KNOWN || gq->state == CHANGED) &&
                   (gr->state == KNOWN || gr->state == CHANGED));
    bool changed = false;

    // Doing it uglily for now : we sample the face, loop on all the points and take the min
    Vector3d v1, v2;
    v1(0) = gq->point.x() - gp->point.x();
    v1(1) = gq->point.y() - gp->point.y();
    v1(2) = gq->point.z() - gp->point.z();
    v2(0) = gr->point.x() - gp->point.x();
    v2(1) = gr->point.y() - gp->point.y();
    v2(2) = gr->point.z() - gp->point.z();

    std::size_t pt_n = 3, tot = pt_n*(pt_n+1);
    FT increment = 1./static_cast<FT>(pt_n);

    // below is a hack with a single loop that is equivalent to the following :

    //for(std::size_t ni=0; ni<=pt_n; ni++)
    //  for(std::size_t nj=0; nj<=pt_n-ni; nj++)

    // that is, we build a grid on a triangle. To parallelize these two loops
    // we create a single loop and to make it less tedious to count, we create
    // a grid of the parallelogram and then ignore the pts that are outside
    // the triangle (when ni+nj > pt_n)

    // actually using the parallel version very much slows the code atm.... :-|

//#pragma omp parallel for
    for(std::size_t nij = 0; nij <= tot; ++nij)
    {
      std::size_t ni = nij / (pt_n+1);
      std::size_t nj = nij % (pt_n+1);

      if(ni+nj > pt_n)
        continue;

      FT i = ni * increment;
      FT j = nj * increment;

//      std::cout << "nij: " << nij << " " << ni << " " << nj << " " << i << " " << j << std::endl;

      Point_3 p(gp->point.x() + i*v1(0) + j*v2(0),
                gp->point.y() + i*v1(1) + j*v2(1),
                gp->point.z() + i*v1(2) + j*v2(2));

      // compute the bary weights
      FT lambda_p, lambda_q, lambda_r;
      Vector3d v3;
      v3(0) = p.x() - gp->point.x();
      v3(1) = p.y() - gp->point.y();
      v3(2) = p.z() - gp->point.z();

      FT d11 = v1.dot(v1);
      FT d12 = v1.dot(v2);
      FT d22 = v2.dot(v2);
      FT d31 = v3.dot(v1);
      FT d32 = v3.dot(v2);
      FT den = 1. / (d11 * d22 - d12 * d12);
      lambda_q = (d22 * d31 - d12 * d32) * den;
      lambda_r = (d11 * d32 - d12 * d31) * den;
      lambda_p = 1. - lambda_q - lambda_r;

      // interpolated value at P
      FT gp_d = gp->distance_to_closest_seed;
      FT gq_d = gq->distance_to_closest_seed;
      FT gr_d = gr->distance_to_closest_seed;
      FT dcs_at_p = lambda_p * gp_d + lambda_q * gq_d + lambda_r * gr_d;

      CGAL_assertion(lambda_p >= -1e-10 && lambda_p <= 1.+1e-10 &&
                     lambda_q >= -1e-10 && lambda_q <= 1.+1e-10 &&
                     lambda_r >= -1e-10 && lambda_r <= 1.+1e-10);

      Vector3d v;
      v(0) = p.x() - point.x();
      v(1) = p.y() - point.y();
      v(2) = p.z() - point.z();
      FT neighbor_d = std::sqrt(v.transpose()*metric.get_mat()*v);
      FT d = dcs_at_p + neighbor_d;

//#pragma omp critical
//{
//      std::cout << "Debug compute_closest_seed_3D" << std::endl;
//      std::cout << "the three points : " << std::endl;
//      std::cout << gp->point.x() << " " << gp->point.x() << " " << gp->point.z() << std::endl;
//      std::cout << gq->point.y() << " " << gq->point.y() << " " << gq->point.z() << std::endl;
//      std::cout << gr->point.z() << " " << gr->point.z() << " " << gr->point.z() << std::endl;
//      std::cout << "ij:" << i << " " << j << " and p : " << p.x() << " " << p.y() << " " << p.z() << std::endl;
//      std::cout << "barywhite: " << lambda_p << " " << lambda_q << " " << lambda_r << std::endl;
//      std::cout << "check: ";
//      std::cout << lambda_p*gp->point.x()+lambda_q*gq->point.x()+lambda_r*gr->point.x() << " ";
//      std::cout << lambda_p*gp->point.y()+lambda_q*gq->point.y()+lambda_r*gr->point.y() << " ";
//      std::cout << lambda_p*gp->point.z()+lambda_q*gq->point.z()+lambda_r*gr->point.z() << std::endl;
//}
      // tolerance to ignore unsignificant changes
      if(distance_to_closest_seed - d > recursive_tolerance * d)
      {
#pragma omp critical
{
        changed = true;
        distance_to_closest_seed = d;
        const Grid_point* best_p;
        // theoretically it should be that
        //          best_p = (lambda_p > lambda_q) ? ((lambda_p > lambda_r) ? gp : gr) :
        //                                           ((lambda_q > lambda_r) ? gq : gr);
        // but again, below is better ?
        best_p = (gp_d > gq_d) ? ((gp_d > gr_d) ? gp : gr) :
                                 ((gq_d > gr_d) ? gq : gr);

        closest_seed_id = best_p->closest_seed_id;
        ancestor = best_p;
}
      }
    }
    std::cout << "end" << std::endl;
    return changed;
  }

  boost::array<const Grid_point*, 4> get_adjacent_neighbors(const Grid_point* gp,
                                                            const Grid_point* gq)
  {
    // Find, in the neighbors of p, the points adjacent to q
    // for example, if q = right, the adjacent neighbors are back, above, front, below

    CGAL_assertion(gp && gq && gp != gq);
    boost::array<const Grid_point*, 4> adj_n;
    for(std::size_t i=0; i<adj_n.size(); ++i)
      adj_n[i] = NULL;

    // an annoying case is the border of the grid: not all neighbors are != NULL
    // so if we want to find the adjacent neighbors, we have to know what is the
    // neighbor opposite of gq (AND THAT OPPOSITE NEIGHBOR VERY WELL COULD BE 'NULL'!!)
    // we therefore use simply the index of the neighbors: we know the opposite couples:
    // 0-5 1-3 2-4

    // first of all, we need to find the index of gq in gp->neighbors
    std::size_t pos_gq = 6; // 6 so that it fails an assert later if gq not found
    for(std::size_t i=0; i<gp->neighbors.size(); ++i)
      if(gp->neighbors[i] == gq)
        pos_gq = i;

    std::size_t pos = 0.;
    for(std::size_t i=0; i<gp->neighbors.size(); ++i)
    {
      const Grid_point* gr = gp->neighbors[i]; // a potential adjacent neighbor
      if(gr == gq) // we don't want to grab gq, only its neighbors
        continue;

      if(pos_gq == 0 && i == 5 || pos_gq == 5 && i == 0 ||
         pos_gq == 2 && i == 4 || pos_gq == 4 && i == 2 ||
         pos_gq == 3 && i == 1 || pos_gq == 1 && i == 3)
        continue; // this is the neighbor opposite of gq in gp, we ignore it

      if(!gr) // everything that is left is a neighbor... as long as it exists
      {
        pos++; // even if there's no neighbor, we must increase the pos
        continue; // to not consider incorrect faces later on
      }

      adj_n[pos++] = gr;
    }
    CGAL_assertion(pos == 4);
    return adj_n;
  }

  bool compute_closest_seed(const Grid_point* ancestor)
  {
//    std::cout << "compute closest high level " << index << std::endl;
    bool changed = false;

    // get the neighbors of 'this' that are adjacent of ancestor
    // then 4 possible 3D faces: this-ancestor-2 from adj_n
    // if a tet fails, we have two possible 2D faces: this-ancestor-1 from tet
    // if no face is available, then compute 1D: this-ancestor
    boost::array<const Grid_point*, 4> adj_n = get_adjacent_neighbors(this, ancestor);

    for(std::size_t i=0; i<adj_n.size(); ++i)
    {
      const Grid_point* gp = adj_n[i];
      const Grid_point* gq = adj_n[(i+1)%(adj_n.size())];

      if(gp && (gp->state == KNOWN || gp->state == CHANGED))
      {
        if(gq && (gq->state == KNOWN || gq->state == CHANGED))
        {
          // the tet 'this - ancestor - gp - gq'
          if(compute_closest_seed_3D(ancestor, gp, gq))
            changed = true;
        }
        else
        { // the face 'this - ancestor - gp'
          if(compute_closest_seed_2D(ancestor, gp))
              changed = true;
        }
      }
      else if(gq && (gq->state == KNOWN || gq->state == CHANGED))
      { // the face 'this - ancestor - gq'
        if(compute_closest_seed_2D(ancestor, gq))
          changed = true;
      }
      else
      {
        if(compute_closest_seed_1D(ancestor))
          changed = true;
      }
    }
    return changed;
  }

  template<typename PQ>
  void update_neighbors_distances(PQ& changed_pq,
                                  PQ& trial_pq) const
  {
//    std::cout << "update neigh" << std::endl;

    // consider all the known neighbors and compute the value to that
    CGAL_assertion(state == KNOWN || state == CHANGED);

    for(std::size_t i=0; i<neighbors.size(); ++i)
    {
      Grid_point* gp = neighbors[i]; // a neighbor of 'this'
      if(!gp)
        continue;

     std::cout << "gp: " << gp->index << " is go, status: " << gp->state << std::endl;

#ifdef USE_RECURSIVE_UPDATES
      if(gp->state == KNOWN || gp->state == CHANGED) // that's the recursive part
      {
        // recompute the distance in the direction gp-'this'
        FT mem = gp->distance_to_closest_seed;
        std::size_t mem_ancestor = (gp->ancestor)?gp->ancestor->index:-1;
        if(gp->compute_closest_seed(this))
        {
          std::cout.precision(17);
          std::cout << "new value at " << gp->index << " : " << gp->distance_to_closest_seed;
          std::cout << " with ancestor: " << gp->ancestor->index;
          std::cout << " (mem: " << mem << " ancestor: " << mem_ancestor << ") ";
          std::cout << " diff: " << gp->distance_to_closest_seed-mem << std::endl;
          if(gp->state == KNOWN)
          {
            // if it's already changed we only need to reorder the changed queue but we don't need to insert anything
            std::cout << "and inserting it" << std::endl;
            gp->state = CHANGED;
            changed_pq.push_back(gp);
            std::push_heap(changed_pq.begin(), changed_pq.end(), Grid_point_comparer<Grid_point>());
          }
        }
      }
      else
#endif
      if(gp->state == TRIAL)
      {
        gp->compute_closest_seed(this); // need to change the order of the element
        // note that we don't insert in trial_pq since it's already in
        // the sort after this function will change its position if needs be
      }
      else if(gp->state == FAR)
      {
        gp->state = TRIAL;
        gp->compute_closest_seed(this); // always returns true here so no need to test it
        trial_pq.push_back(gp);
        std::push_heap(trial_pq.begin(), trial_pq.end(), Grid_point_comparer<Grid_point>());
      }
    }
  }

  bool is_contained(const Point_3& p,
                    const Grid_point* next) const
  {
    // above : n0
    //    n2
    //    |
    // n1-gp-n3
    //    |
    //    n4
    // below : n5

    Grid_point *n0 = neighbors[0], *n1 = neighbors[1], *n2 = neighbors[2],
               *n3 = neighbors[3], *n4 = neighbors[4], *n5 = neighbors[5];

    bool is_p_geq_n0_z = (n1) ? p.z() <= n0->point.z() : false;
    bool is_p_leq_n1_x = (n1) ? p.x() <= n1->point.x() : false;
    bool is_p_geq_n2_y = (n2) ? p.y() >= n2->point.y() : false;
    bool is_p_geq_n3_x = (n3) ? p.x() >= n3->point.x() : false;
    bool is_p_leq_n4_y = (n4) ? p.y() <= n4->point.y() : false;
    bool is_p_leq_n5_z = (n5) ? p.z() <= n5->point.z() : false;

    if(is_p_geq_n0_z) { next = n0; return false; }
    if(is_p_leq_n1_x) { next = n1; return false; }
    if(is_p_geq_n2_y) { next = n2; return false; }
    if(is_p_geq_n3_x) { next = n3; return false; }
    if(is_p_leq_n4_y) { next = n4; return false; }
    if(is_p_leq_n5_z) { next = n5; return false; }
    next = this; // we're inside !
    return true;
  }

  void initialize_from_point(const Point_3& p,
                             const std::size_t seed_id) const
  {
    Vector3d v;
    v(0) = p.x() - point.x();
    v(1) = p.y() - point.y();
    v(2) = p.z() - point.z();
    FT d = std::sqrt(v.transpose()*metric.get_mat()*v);
    closest_seed_id = seed_id;
    distance_to_closest_seed = d;
    state = TRIAL;
  }

  Grid_point(const Point_3& point_, const std::size_t index_)
    : point(point_), index(index_), distance_to_closest_seed(1e30),
      closest_seed_id(-1), state(FAR), neighbors(),
      metric(mf->compute_metric(point)), ancestor(NULL)
  {
    for(std::size_t i=0; i<neighbors.size(); ++i)
      neighbors[i] = NULL;
  }
};

struct Geo_grid
{
  typedef std::vector<Grid_point*> Grid_point_vector;

  std::vector<Grid_point> points;
  Grid_point_vector changed_points;
  Grid_point_vector trial_points;

  const Grid_point* locate_point(const Point_3& p) const
  {
    bool found = false;
    const Grid_point& next = points[0];
    while(!found)
      found = (&next)->is_contained(p, &next);
    return &next;
  }

  void locate_and_initialize(const Point_3& p,
                             const std::size_t seed_id)
  {
    // WARNING: only usable for a full grid since we assume the step size is constant
    int index_x = std::floor((p.x()-offset_x)/step);
    int index_y = std::floor((p.y()-offset_y)/step);
    int index_z = std::floor((p.z()-offset_z)/step);

    Grid_point* gp = &(points[index_z*sq_n + index_y*n + index_x]);
    if(gp->closest_seed_id != -1)
    {
      std::cerr << "WARNING: a new seed is overwriting the closest seed id";
      std::cerr << " of a grid point! previous index is: " << gp->closest_seed_id << std::endl;
    }

    // if gp is found the same for 2 different initial points, then your grid
    // is not dense enough... (or you're refining !)

    gp->initialize_from_point(p, seed_id);
    trial_points.push_back(gp);
    std::push_heap(trial_points.begin(), trial_points.end(),Grid_point_comparer<Grid_point>());

    std::cout << "step reminder: " << step << std::endl;
    std::cout << "looking for p: " << p.x() << " " << p.y() << " " << p.z() << std::endl;
    std::cout << "found gp: " << gp->index << " [" << gp->point.x() << ", " << gp->point.y() << " " << gp->point.z() << "]" << std::endl;
    std::cout << "index: " << index_x << " " << index_y << " " << index_z << std::endl;
    std::cout << "offset: " << offset_x << " " << offset_y << " " << offset_z << std::endl;

    std::sort_heap(trial_points.begin(), trial_points.end(), Grid_point_comparer<Grid_point>());
  }

  void initialize_geo_grid()
  {
    std::cout << "grid ini" << std::endl;

    // create the grid points (either a full grid or a smart grid (todo) )
    for(unsigned int k=0; k<n; ++k)
    {
      for(unsigned int j=0; j<n; ++j)
      {
        for(unsigned int i=0; i<n; ++i)
        {
          Point_3 p(offset_x+i*step, offset_y+j*step, offset_z+k*step); // fill from bot left to top right
          Grid_point gp(p, i + j*n + k*sq_n);
          points.push_back(gp);
        }
      }
    }

    // assign the neighbors (todo: find a better way...)
    for(unsigned int k=0; k<n; ++k)
    {
      for(unsigned int j=0; j<n; ++j)
      {
        for(unsigned int i=0; i<n; ++i)
        {
          std::size_t curr_id = i + j*n + k*sq_n;
          if(k != n-1) // there is a neighbor above
            points[curr_id].neighbors[0] = &(points[i + j*n + (k+1)*sq_n]);
          if(i != 0) // there is a neighbor left
            points[curr_id].neighbors[1] = &(points[i-1 + j*n + k*sq_n]);
          if(j != n-1) // there is a neighbor back
            points[curr_id].neighbors[2] = &(points[i + (j+1)*n + k*sq_n]);
          if(i != n-1) // there is a neighbor right
            points[curr_id].neighbors[3] = &(points[i+1 + j*n + k*sq_n]);
          if(j != 0) // there is a neighbor front
            points[curr_id].neighbors[4] = &(points[i + (j-1)*n + k*sq_n]);
          if(k != 0) // there is a neighbor below
            points[curr_id].neighbors[5] = &(points[i + j*n + (k-1)*sq_n]);
        }
      }
    }
    std::cout << "assigned neighbors" << std::endl;

    // seeds belong to a quad, find it and initialize the pts of the quads
    for(std::size_t i=0; i<seeds.size(); ++i)
    {
      const Point_3& p = seeds[i];
//      const Grid_point* gp = locate_point(p);
//      gp->initialize_from_point(p, i);
      locate_and_initialize(p, i);
    }

    std::cout << "grid initialization end" << std::endl;
  }

  bool print_states() const
  {
    int changed_count=0, known_count=0, trial_count=0, far_count=0;

    for(std::size_t i=0; i<points.size(); ++i)
    {
      const Grid_point& gp = points[i];
      if(gp.state == CHANGED)
        changed_count++;
      if(gp.state == KNOWN)
        known_count++;
      if(gp.state == TRIAL)
        trial_count++;
      if(gp.state == FAR)
        far_count++;
    }
    std::cout << "changed: " << changed_count << " known: " << known_count;
    std::cout << " trial: " << trial_count << " far: " << far_count << std::endl;

    return false;
  }

  void geo_grid_loop()
  {
    std::cout << "main loop" << std::endl;

    bool is_cp_empty = changed_points.empty();
    bool is_t_empty = trial_points.empty();

    Grid_point* gp;
    while(!is_cp_empty || !is_t_empty)
    {
      std::cout << "changed heap: " << std::endl;
      for (std::vector<Grid_point*>::iterator it = changed_points.begin();
           it != changed_points.end(); ++it)
        std::cout << (*it)->index << " " << (*it)->distance_to_closest_seed << std::endl;
      std::cout << std::endl;

      std::cout << "trial heap: " << std::endl;
      for (std::vector<Grid_point*>::iterator it = trial_points.begin();
           it != trial_points.end(); ++it)
        std::cout << (*it)->index << " " << (*it)->distance_to_closest_seed << std::endl;
      std::cout << std::endl;

      if(!is_cp_empty)
      {
        gp = changed_points.front();
        std::pop_heap(changed_points.begin(), changed_points.end(), Grid_point_comparer<Grid_point>());
        changed_points.pop_back();
      }
      else // !is_t_empty
      {
        gp = trial_points.front();
        std::pop_heap(trial_points.begin(), trial_points.end(), Grid_point_comparer<Grid_point>());
        trial_points.pop_back();
      }

      if(print_states())
        break;
      std::cout << "Queue sizes. Trial: " << trial_points.size() << " Changed: " << changed_points.size() << std::endl;
      std::cout << "picked: " << gp->index << " | " << gp->point.x() << " " << gp->point.y() << " " << gp->point.z() << std::endl;

      gp->state = KNOWN;
      gp->update_neighbors_distances(changed_points, trial_points);
      std::sort_heap(trial_points.begin(), trial_points.end(), Grid_point_comparer<Grid_point>());
      std::sort_heap(changed_points.begin(), changed_points.end(), Grid_point_comparer<Grid_point>());

      is_cp_empty = changed_points.empty();
      is_t_empty = trial_points.empty();
    }

    std::cout << "DOUBLE CHECK THAT WE ARE REALLY FINISHED : " << std::endl;
    for(std::size_t i=0; i<points.size(); ++i)
    {
      Grid_point* gp = &(points[i]);
      for(std::size_t i; i<gp->neighbors.size(); ++i)
      {
        const Grid_point* n = gp->neighbors[i];
        if(n)
          CGAL_assertion(!(gp->compute_closest_seed(n)));
      }
    }
  }

  void build_grid()
  {
    initialize_geo_grid();
    geo_grid_loop();
  }

  Point_3 compute_refinement_point()
  {
    // todo
  }

  void refresh_grid_after_new_seed_creation()
  {
    std::size_t seed_id = seeds.size() - 1;
    const Point_3& p = seeds.back();
    locate_and_initialize(p, seed_id);
    geo_grid_loop();
  }

  void refine_grid()
  {
#ifdef TMP_REFINEMENT_UGLY_HACK
    std::cout << "insert: " << best_ref_x << " " << best_ref_y << " " << best_ref_z << std::endl;
    Point_3 p(best_ref_x, best_ref_y, best_ref_z);
    CGAL_assertion(best_ref_x != 1e30 && best_ref_y != 1e30 && best_ref_z != 1e30);
#else
    Point_3 p = compute_refinement_point();
#endif
    vertices_nv = insert_new_seed(p.x(), p.y(), p.z());
    refresh_grid_after_new_seed_creation();

#ifdef TMP_REFINEMENT_UGLY_HACK
    best_ref_x = 1e30; best_ref_y = 1e30; best_ref_z = 1e30;
#endif
  }

  void add_tet_to_simplices(std::set<Tet>& simplices, const Grid_point& gp,
                            const Grid_point* n1, const Grid_point* n2, const Grid_point* n3) const
  {
    Tet t;
    if(n1 && n2 && n3)
    {
      t[0] = gp.index;
      t[1] = n1->index;
      t[2] = n2->index;
      t[3] = n3->index;
      std::sort(t.begin(), t.end());
      simplices.insert(t);
    }
  }

  void add_tet_to_simplices(std::set<Tet>& simplices, const Grid_point& gp,
                            std::size_t i1, std::size_t i2, std::size_t i3) const
  {
    return add_tet_to_simplices(simplices, gp, gp.neighbors[i1], gp.neighbors[i2], gp.neighbors[i3]);
  }

  void output_grid(const std::string str_base) const
  {
    std::ofstream out((str_base + ".mesh").c_str());
    std::ofstream out_bb((str_base + ".bb").c_str());
    out << "MeshVersionFormatted 1" << std::endl;
    out << "Dimension 3" << std::endl;
    out << "Vertices" << std::endl;
    out << points.size() << std::endl;
    out_bb << "3 1 " << points.size() << " 2" << std::endl;

    std::set<Tet> simplices;
    for(std::size_t i=0; i<points.size(); ++i)
    {
      const Grid_point& gp = points[i];
      out << gp.point.x() << " " << gp.point.y() << " " << gp.point.z() << " " << i+1 << std::endl;

//      out_bb << gp.closest_seed_id << std::endl;
      out_bb << gp.distance_to_closest_seed << std::endl;

      // compute the tets... todo
      const boost::array<Grid_point*, 6>& ns = gp.neighbors; // just to get a shorter name...
      for(std::size_t i=0; i<ns.size(); ++i)
      {
        // lazily doing something overlapping for now : the eight corner tetrahedra
#if 0
        add_tet_to_simplices(simplices, gp, 0, 1, 2);
        add_tet_to_simplices(simplices, gp, 0, 2, 3);
        add_tet_to_simplices(simplices, gp, 0, 3, 4);
        add_tet_to_simplices(simplices, gp, 0, 4, 1);
        add_tet_to_simplices(simplices, gp, 5, 1, 2);
        add_tet_to_simplices(simplices, gp, 5, 2, 3);
        add_tet_to_simplices(simplices, gp, 5, 3, 4);
        add_tet_to_simplices(simplices, gp, 5, 4, 1);
#else
        add_tet_to_simplices(simplices, gp, 0, 2, 3);
        add_tet_to_simplices(simplices, gp, 0, 1, 4);
        add_tet_to_simplices(simplices, gp, 3, 4, 5);
        add_tet_to_simplices(simplices, gp, 1, 2, 5);

        // the last one is a bit annoying since we have to grab neighbors of neighbors
        const Grid_point* n4 = gp.neighbors[4];
        const Grid_point* n3 = gp.neighbors[3];
        if(n4 && n3)
          add_tet_to_simplices(simplices, gp, n4->neighbors[3],
                               n4->neighbors[0], n3->neighbors[0]);
#endif
      }
    }

    out << "Tetrahedra" << std::endl;
    out << simplices.size() << std::endl;
    for(typename std::set<Tet>::iterator it = simplices.begin();
                                             it != simplices.end(); ++it)
    {
      const Tet& tet = *it;
      std::set<std::size_t> materials;
      for(std::size_t i=0; i<tet.size(); ++i)
      {
        out << tet[i] + 1 << " ";
        materials.insert(points[tet[i]].closest_seed_id);
      }
      std::size_t mat = (materials.size()==1)?(*(materials.begin())):(seeds.size());
      out << mat << std::endl;
    }

    std::set<Triangle> triangles;
    for(typename std::set<Tet>::iterator it = simplices.begin();
                                         it != simplices.end(); ++it)
    {
      const Tet& tet = *it;
      Triangle tr; // could use a permutation to make it look nicer I suppose...
      tr[0] = tet[0]; tr[1] = tet[1]; tr[2] = tet[2];
      std::sort(tr.begin(), tr.end()); triangles.insert(tr);
      tr[0] = tet[0]; tr[1] = tet[1]; tr[2] = tet[3];
      std::sort(tr.begin(), tr.end()); triangles.insert(tr);
      tr[0] = tet[0]; tr[1] = tet[2]; tr[2] = tet[3];
      std::sort(tr.begin(), tr.end()); triangles.insert(tr);
      tr[0] = tet[1]; tr[1] = tet[2]; tr[2] = tet[3];
      std::sort(tr.begin(), tr.end()); triangles.insert(tr);
    }

    out << "Triangles" << std::endl;
    out << triangles.size() << std::endl;
    for(typename std::set<Triangle>::iterator it = triangles.begin();
                                              it != triangles.end(); ++it)
    {
      const Triangle& tr = *it;
      std::set<std::size_t> materials;
      for(std::size_t i=0; i<tr.size(); ++i)
      {
        out << tr[i] + 1 << " ";
        materials.insert(points[tr[i]].closest_seed_id);
      }
      std::size_t mat = (materials.size()==1)?(*(materials.begin())):(seeds.size());
      out << mat << std::endl;
    }
  }

  std::set<Simplex> compute_dual() const
  {
    // equivalent of a "smart dual" for grid_generator
    std::cout << "dual computations" << std::endl;

#ifdef TMP_REFINEMENT_UGLY_HACK
    FT farthest_dual_distance = 0.;
#endif

    std::set<Simplex> dual_simplices;
    // a dual exists if a simplex has 4 different values
#pragma omp parallel for
    for(std::size_t i=0; i<points.size(); ++i)
    {
      const Grid_point& gp = points[i];
      Simplex dual_simplex;
      dual_simplex.insert(gp.closest_seed_id);

      boost::array<Grid_point*, 6>::const_iterator it = gp.neighbors.begin(),
                                                   iend = gp.neighbors.end();
      for(; it!=iend; ++it)
      {
        const Grid_point* gq = *it;
        if(!gq)
          continue;

        dual_simplex.insert(gq->closest_seed_id);
      }

       // no one cares if dual_simplex.size() == 1
      if(dual_simplex.size() == 2)
      {
        // that'd be an edge in the dual if we cared about it
      }
      else if(dual_simplex.size() == 3)
      {
        // that'd be a triangle in the dual if we cared about it
      }
      else if(dual_simplex.size() >= 4)
      {
        if(dual_simplex.size() == 4)
#pragma omp critical
          dual_simplices.insert(dual_simplex); // a tetrahedron!
        else
          std::cout << "WARNING: COSPHERICITY..." << std::endl;
#ifdef TMP_REFINEMENT_UGLY_HACK
        //todo : this could be gp + its neighbors instead of gp alone
        if(gp.distance_to_closest_seed > farthest_dual_distance)
        {
#pragma omp critical
{
          farthest_dual_distance = gp.distance_to_closest_seed;
          best_ref_x = gp.point.x();
          best_ref_y = gp.point.y();
          best_ref_z = gp.point.z();
}
        }
#endif
      }
    }
#ifdef TMP_REFINEMENT_UGLY_HACK
    if(dual_simplices.empty())
      std::cerr << "WARNING: no simplices captured, so no best_ref_x,y set" << std::endl;
#endif
    return dual_simplices;
  }

  void output_dual() const
  {
    const std::set<Simplex>& dual_simplices = compute_dual();
    std::cout << "captured: " << dual_simplices.size() << " simplices" << std::endl;

    std::ofstream out("geo_grid_dual.mesh");
    out << "MeshVersionFormatted 1" << std::endl;
    out << "Dimension 3" << std::endl;
    out << "Vertices" << std::endl;
    out << seeds.size() << std::endl;
    for(int i=0; i<seeds.size(); ++i)
      out << seeds[i].x() << " " << seeds[i].y() << " " << seeds[i].z() << " " << i+1 << std::endl;

    out << "Tetrahedra" << std::endl;
    out << dual_simplices.size() << std::endl;
    for(typename std::set<Simplex>::iterator it = dual_simplices.begin();
                                             it != dual_simplices.end(); ++it)
    {
      const Simplex& s = *it;
      CGAL_assertion(s.size() == 4);
      typename Simplex::iterator sit = s.begin();
      typename Simplex::iterator siend = s.end();
      for(; sit!=siend; ++sit)
        out << *sit+1 << " ";
      out << "1" << std::endl;
    }

    std::set<Triangle> triangles;
    for(typename std::set<Simplex>::iterator it = dual_simplices.begin();
                                             it != dual_simplices.end(); ++it)
    {
      const Simplex& s = *it;
      Tet tet;
      typename Simplex::iterator sit = s.begin();
      typename Simplex::iterator siend = s.end();
      for(std::size_t i=0; sit!=siend; ++sit, ++i)
        tet[i] = *sit;

      Triangle tr; // could use a permutation to make it look nicer I suppose...
      tr[0] = tet[0]; tr[1] = tet[1]; tr[2] = tet[2];
      std::sort(tr.begin(), tr.end()); triangles.insert(tr);
      tr[0] = tet[0]; tr[1] = tet[1]; tr[2] = tet[3];
      std::sort(tr.begin(), tr.end()); triangles.insert(tr);
      tr[0] = tet[0]; tr[1] = tet[2]; tr[2] = tet[3];
      std::sort(tr.begin(), tr.end()); triangles.insert(tr);
      tr[0] = tet[1]; tr[1] = tet[2]; tr[2] = tet[3];
      std::sort(tr.begin(), tr.end()); triangles.insert(tr);
    }

    out << "Triangles" << std::endl;
    out << triangles.size() << std::endl;
    for(typename std::set<Triangle>::iterator it = triangles.begin();
                                              it != triangles.end(); ++it)
    {
      const Triangle& tr = *it;
      for(std::size_t i=0; i<tr.size(); ++i)
        out << tr[i] + 1 << " ";
      out << "1" << std::endl;
    }

    out << "End" << std::endl;
  }

  Geo_grid() : points(), changed_points(), trial_points()
  {
    std::make_heap(changed_points.begin(), changed_points.end(), Grid_point_comparer<Grid_point>());
    std::make_heap(trial_points.begin(), trial_points.end(), Grid_point_comparer<Grid_point>());
  }
};

void initialize()
{
  mf = new Euclidean_metric_field<K>(1,1,1);
//  mf = new Custom_metric_field<K>();

  vertices_nv = build_seeds();
}

int main(int, char**)
{
  std::clock_t start;
  double duration;
  start = std::clock();

  std::freopen("geo_grid_log.txt", "w", stdout);
  std::srand(0);
  initialize();

  Geo_grid gg;
  gg.build_grid();
  gg.output_grid("geo_grid");

  int n_refine = 0;
  for(int i=0; i<n_refine; ++i)
  {
    gg.refine_grid();
    std::ostringstream out;
    out << "geo_grid_ref_" << i;
    gg.output_grid(out.str());
  }

  gg.output_grid("geo_grid");
  gg.output_dual();

  duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
  std::cerr << "dur: " << duration << std::endl;
}
