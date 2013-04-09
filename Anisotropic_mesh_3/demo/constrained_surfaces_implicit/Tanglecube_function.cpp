
#include <QObject>
#include "Constrained_surface_implicit_interface.h"
#include <Domain/Constrain_surface_3_tanglecube.h>
#include "Implicit_surface_type.h"

class Tanglecube_function :
  public QObject,
  public Constrained_surface_implicit_interface
{
  Q_OBJECT
  Q_INTERFACES(Constrained_surface_implicit_interface)

private:
  Constrain_surface_3_tanglecube<Kernel>* pSurface;

public:
  Implicit_surface* surface() const { return pSurface; }

  virtual QString name() const
  {
    std::ostringstream o;
    o << "Tanglecube (" << pSurface->get_stretch() << ")";
    return QString(o.str().data());
  }

  virtual std::string parameters() const
  {
    std::ostringstream o;
    o << pSurface->get_stretch();
    return o.str();
  }

  virtual bool set_parameters(const std::string& str)
  { 
    std::vector<double> vec = this->get_values(str);
    if(vec.size() != 1)
      return false;

    pSurface->set_stretch(vec[0]);
    return true;
  }

  Tanglecube_function() : pSurface(new Constrain_surface_3_tanglecube<Kernel>()) {}
  ~Tanglecube_function() { delete pSurface; }
};


#include <QtPlugin>
Q_EXPORT_PLUGIN2(Tanglecube_function, Tanglecube_function)
#include "Tanglecube_function.moc"
