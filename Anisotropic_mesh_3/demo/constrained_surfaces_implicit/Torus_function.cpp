
#include <QObject>
#include "Constrained_surface_implicit_interface.h"
#include <Anisotropic_mesh_3/include/Domain/Constrain_surface_3_torus.h>
#include <Anisotropic_mesh_3/demo/Implicit_surface_type.h>

class Torus_function :
  public QObject,
  public Constrained_surface_implicit_interface
{
  Q_OBJECT
  Q_INTERFACES(Constrained_surface_implicit_interface)

private:
  Constrain_surface_3_torus<Kernel>* pSurface;

public:
  Implicit_surface* surface() const { return pSurface; }

  virtual QString name() const
  {
    std::ostringstream o;
    o << "Torus (" << pSurface->get_r() << ", " << pSurface->get_R() << ")";
    return QString(o.str().data());
  }

  virtual std::string parameters() const
  {
    std::ostringstream o;
    o << pSurface->get_r() << ", " << pSurface->get_R();
    return o.str();
  }

  virtual bool set_parameters(const std::string& str)
  { 
    std::vector<double> vec = this->get_values(str);
    if(vec.size() != 2)
      return false;

    pSurface->set_r(vec[0]);
    pSurface->set_R(vec[1]);
    return true; 
  }

  Torus_function() : pSurface(new Constrain_surface_3_torus<Kernel>()) {}
  ~Torus_function() { delete pSurface; }

};


#include <QtPlugin>
Q_EXPORT_PLUGIN2(Torus_function, Torus_function)
#include "Torus_function.moc"
