
#include <QObject>
#include "Constrained_surface_implicit_interface.h"
#include <Anisotropic_mesh_3/include/Domain/Constrain_surface_3_sphere.h>
#include <Anisotropic_mesh_3/demo/Implicit_surface_type.h>

class Sphere_function :
  public QObject,
  public Constrained_surface_implicit_interface
{
  Q_OBJECT
  Q_INTERFACES(Constrained_surface_implicit_interface)

private:
  Constrain_surface_3_sphere<Kernel>* pSurface;

public:  
  Implicit_surface* surface() const { return pSurface; }
  
  virtual QString name() const
  {
    return QString("Unit Sphere");
  }

  virtual std::string parameters() const { return std::string(); }

  virtual bool set_parameters(const std::string& str){ return true; } //nothing to do

  Sphere_function() : pSurface(new Constrain_surface_3_sphere<Kernel>()) {}
  ~Sphere_function() { delete pSurface; }

};


#include <QtPlugin>
Q_EXPORT_PLUGIN2(Sphere_function, Sphere_function)
#include "Sphere_function.moc"
