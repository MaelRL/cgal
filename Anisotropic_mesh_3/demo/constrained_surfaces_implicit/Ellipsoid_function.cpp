#include <QObject>
#include "Constrained_surface_implicit_interface.h"
#include <Domain/Constrain_surface_3_ellipse.h>
#include "Implicit_surface_type.h"

class Ellipsoid_function :
  public QObject,
  public Constrained_surface_implicit_interface
{
  Q_OBJECT
  Q_INTERFACES(Constrained_surface_implicit_interface)

private:
  Constrain_surface_3_ellipse<Kernel>* pSurface;

public:
  Implicit_surface* surface() const { return pSurface; }

  virtual QString name() const
  {
    std::ostringstream o;
    o << "Ellipsoid (" << pSurface->get_a() << ", " << pSurface->get_b() << ", " << pSurface->get_c() << ")";
    return QString(o.str().data());
  }

  virtual std::string parameters() const
  {
    std::ostringstream o;
    o << pSurface->get_a() << ", " << pSurface->get_b() << ", " << pSurface->get_c();
    return o.str();
  }

  virtual bool set_parameters(const std::string& str)
  {
    std::vector<double> vec = this->get_values(str);
    if(vec.size() != 3)
      return false;

    pSurface->set_a(vec[0]);
    pSurface->set_b(vec[1]);
    pSurface->set_c(vec[2]);
    return true;
  }

  Ellipsoid_function() : pSurface(new Constrain_surface_3_ellipse<Kernel>()) {}
  ~Ellipsoid_function() { delete pSurface; }

};


#include <QtPlugin>
Q_EXPORT_PLUGIN2(Ellipsoid_function, Ellipsoid_function)
#include "Ellipsoid_function.moc"
