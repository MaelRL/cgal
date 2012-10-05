
#ifndef CGAL_CONSTRAINED_IMPLICIT_INTERFACE_H
#define CGAL_CONSTRAINED_IMPLICIT_INTERFACE_H

#include <QObject>
#include <QString>
#include <CGAL_demo/Scene_interface.h>
#include <Implicit_surface_type.h>


class Constrained_surface_implicit_interface 
{
public:
  virtual QString name() const = 0;
  virtual Implicit_surface* surface() const = 0;
  virtual std::string parameters() const = 0;
  virtual bool set_parameters(const std::string& str) = 0;

public:
  std::vector<double> get_values(const std::string& str) const
  {
    std::cout << "get_values in " << str << std::endl;
    std::vector<double> values;    
    std::size_t str_end = str.size();
    if(str_end == 0)
      return values;

    std::string tmp_str = str;
    std::size_t index = 0;
    do
    {       
      tmp_str = tmp_str.substr(index, std::string::npos);
      index = tmp_str.find(",");
      std::string str_index = tmp_str.substr(0, index);

      std::istringstream is(str_index);
      double d;
      if(! (is >> d))
        std::cerr << "Warning : wrong parameters ("<< str <<")!\n";
      values.push_back(d);

      if(index == std::string::npos)
        break;
      index++;//to skip the ","
    }
    while(index < str_end);

    return values;
  }
};

Q_DECLARE_INTERFACE(Constrained_surface_implicit_interface,
  "com.geometryfactory.Mesh3Demo.Constrained_surface_implicit_interface/1.0")

#endif // CGAL_CONSTRAINED_IMPLICIT_INTERFACE_H
