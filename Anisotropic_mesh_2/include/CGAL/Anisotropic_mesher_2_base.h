#ifndef CGAL_ANISOTROPIC_MESH_2_MESHER_BASE_H
#define CGAL_ANISOTROPIC_MESH_2_MESHER_BASE_H

namespace CGAL
{
namespace Anisotropic_mesh_2
{

template<typename PointType>
struct No_condition
{
  bool operator()(const PointType& p) const
  {
    return true;
  }
};

class Anisotropic_mesher_2_base
{
  typedef Anisotropic_mesher_2_base Self;

public:
  virtual double refine_mesh(const bool refine_consistency = true) = 0;
  virtual void resume_from_mesh_file(const char* filename) = 0;

  // Step-by-step methods
  virtual void initialize() = 0;
  virtual void one_step() = 0;
  virtual bool is_algorithm_done() = 0;
  virtual void report() = 0;

public:
  Anisotropic_mesher_2_base() { }
  virtual ~Anisotropic_mesher_2_base() { }

private:
  Anisotropic_mesher_2_base(const Self& src);
  Self& operator=(const Self& src);
}; // Anisotropic_mesher_2_base

}  // Anisotropic_mesh_2
}  // CGAL

#endif // CGAL_ANISOTROPIC_MESH_2_MESHER_BASE_H
