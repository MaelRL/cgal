#include <QObject>
#include "Constrained_surface_implicit_interface.h"
#include <Domain/Constrain_surface_3_double_ellipsoid.h>
#include "Implicit_surface_type.h"

class Double_ellipsoid_function :
  public QObject,
  public Constrained_surface_implicit_interface
{
  Q_OBJECT
  Q_INTERFACES(Constrained_surface_implicit_interface)

private:
  Constrain_surface_3_double_ellipsoid<Kernel>* pSurface;

public:
  Implicit_surface* surface() const { return pSurface; }

  virtual QString name() const
  {
    std::ostringstream o;
    o << "Double ellipsoid ( all the get to fill later )";
    return QString(o.str().data());
  }

  virtual std::string parameters() const
  {
    std::ostringstream o;
    return o.str();
  }

  virtual bool set_parameters(const std::string& str)
  {
    std::vector<double> vec = this->get_values(str);
    return true;
  }

  Double_ellipsoid_function() : pSurface(new Constrain_surface_3_double_ellipsoid<Kernel>()) {}
  ~Double_ellipsoid_function() { delete pSurface; }

};


#include <QtPlugin>
Q_EXPORT_PLUGIN2(Double_ellipsoid_function, Double_ellipsoid_function)
#include "Double_ellipsoid_function.moc"
