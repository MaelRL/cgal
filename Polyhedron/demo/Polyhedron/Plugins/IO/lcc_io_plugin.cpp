#include <CGAL/Three/Three.h>
#include <CGAL/Three/Polyhedron_demo_io_plugin_interface.h>
#include <CGAL/Combinatorial_map_save_load.h>
#include <CGAL/Polyhedron_3_to_lcc.h>

#include "Scene_lcc_item.h"

#include <iostream>
#include <fstream>

class LCC_io_plugin :
    public QObject,
    public CGAL::Three::Polyhedron_demo_io_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_io_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.IOPluginInterface/1.0" FILE "lcc_io_plugin.json")
  
public:
  bool isDefaultLoader(const CGAL::Three::Scene_item *item) const 
  { 
    if(qobject_cast<const Scene_lcc_item*>(item)) 
      return true; 
    return false;
  }
  QString name() const { return "lcc_plugin"; }
  QString nameFilters() const { return 
        "OFF files (*.off);;"
        "3-map files (*.3map)"; }
  
  QString saveNameFilters() const {
    return 
        "3-map files (*.3map)"; 
  }
  
  bool canLoad() const { return true; }
  CGAL::Three::Scene_item* load(QFileInfo fileinfo){
    // Open file
    std::ifstream ifs(fileinfo.filePath().toUtf8());
    if(!ifs) {
      std::cerr << "Error! Cannot open file " << (const char*)fileinfo.filePath().toUtf8() << std::endl;
      return nullptr;
    }
    
    Scene_lcc_item::LCC lcc;
    QString ext = fileinfo.suffix();
    bool res = true;
    if(ext == "off")
      CGAL::import_from_polyhedron_3_flux < Scene_lcc_item::LCC > (lcc, ifs);
    else
    {
      res = CGAL::load_combinatorial_map(ifs, lcc);
    }
    if(!res)
    {
      return nullptr;
    }
    Scene_lcc_item* new_item = new Scene_lcc_item(lcc);
    new_item->setName(fileinfo.fileName());
    new_item->invalidateOpenGLBuffers();
    return new_item;
  }
  
  
  bool canSave(const CGAL::Three::Scene_item*){return true;}
  bool save(const CGAL::Three::Scene_item*, QFileInfo fileinfo){
    qDebug() << fileinfo.fileName() << "Implement me, you fool !";
    return false;
  }
  
};


#include "lcc_io_plugin.moc"
