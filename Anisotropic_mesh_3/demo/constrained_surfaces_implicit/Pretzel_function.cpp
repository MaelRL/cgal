
#include <QObject>
#include <constrained_surfaces_implicit/Constrained_surface_implicit_interface.h>
#include <Domain/Constrain_surface_3_pretzel.h>
#include <Implicit_surface_type.h>

using CGAL::Anisotropic_mesh_3::Constrain_surface_3_pretzel;

class Pretzel_function :
  public QObject,
  public Constrained_surface_implicit_interface
{
  Q_OBJECT
  Q_INTERFACES(Constrained_surface_implicit_interface)

private:
  Constrain_surface_3_pretzel<Kernel>* pSurface;

public:
  Implicit_surface* surface() const { return pSurface; }

  virtual QString name() const
  {
    std::ostringstream o;
    o << "Pretzel (" << pSurface->get_c() << ", " << pSurface->get_d() << ")";
    return QString(o.str().data());
  }

  virtual std::string parameters() const
  {
    std::ostringstream o;
    o << pSurface->get_c() << ", " << pSurface->get_d();
    return o.str();
  }

  virtual bool set_parameters(const std::string& str)
  { 
    std::vector<double> vec = this->get_values(str);
    if(vec.size() != 2)
      return false;

    pSurface->set_c(vec[0]);
    pSurface->set_d(vec[1]);
    return true; 
  }

  Pretzel_function() : pSurface(new Constrain_surface_3_pretzel<Kernel>()) {}
  ~Pretzel_function() { delete pSurface; }
};


#include <QtPlugin>
Q_EXPORT_PLUGIN2(Pretzel_function, Pretzel_function)
#include "Pretzel_function.moc"
