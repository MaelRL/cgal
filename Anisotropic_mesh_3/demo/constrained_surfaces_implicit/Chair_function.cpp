
#include <QObject>
#include <constrained_surfaces_implicit/Constrained_surface_implicit_interface.h>
#include <Domain/Constrain_surface_3_chair.h>
#include <Implicit_surface_type.h>

class Chair_function :
  public QObject,
  public Constrained_surface_implicit_interface
{
  Q_OBJECT
  Q_INTERFACES(Constrained_surface_implicit_interface)

private:
  Constrain_surface_3_chair<Kernel>* pSurface;

public:
  Implicit_surface* surface() const { return pSurface; }

  virtual QString name() const
  {
    std::ostringstream o;
    o << "Chair (" << pSurface->get_a() << ", " << pSurface->get_b() << ", " << pSurface->get_k() << ")";
    return QString(o.str().data());
  }

  virtual std::string parameters() const
  {
    std::ostringstream o;
    o << pSurface->get_a() << ", " << pSurface->get_b() << ", " << pSurface->get_k();
    return o.str();
  }

  virtual bool set_parameters(const std::string& str)
  { 
    std::vector<double> vec = this->get_values(str);
    if(vec.size() != 3)
      return false;

    pSurface->set_a(vec[0]);
    pSurface->set_b(vec[1]);
    pSurface->set_k(vec[2]);
    return true; 
  }

  Chair_function() : pSurface(new Constrain_surface_3_chair<Kernel>()) {}
  ~Chair_function() { delete pSurface; }
};


#include <QtPlugin>
Q_EXPORT_PLUGIN2(Chair_function, Chair_function)
#include "Chair_function.moc"
