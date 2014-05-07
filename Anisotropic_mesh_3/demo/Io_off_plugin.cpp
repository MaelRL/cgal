#include "Scene_polyhedron_item.h"
#include "Scene_polygon_soup.h"
#include "Scene_starset3_item.h"
#include "Polyhedron_type.h"

#include <CGAL_demo/Io_plugin_interface.h>
#include <fstream>

class Io_off_plugin :
  public QObject,
  public Io_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(Io_plugin_interface)

public:
  QStringList nameFilters() const;
  bool canLoad() const;
  Scene_item* load(QFileInfo fileinfo);

  bool canSave(const Scene_item*);
  bool save(const Scene_item*, QFileInfo fileinfo);
  virtual QString name() const { return "Io_off_plugin"; }
};

QStringList Io_off_plugin::nameFilters() const {
  return QStringList() << "OFF files (*.off)";
}

bool Io_off_plugin::canLoad() const {
  return true;
}


Scene_item* 
Io_off_plugin::load(QFileInfo fileinfo) 
{
  // Open file
  std::ifstream in(fileinfo.filePath().toUtf8());
  if(!in) {
    std::cerr << "Error! Cannot open file " << (const char*)fileinfo.filePath().toUtf8() << std::endl;
    return NULL;
  }
  else
    std::cout << "File " << (const char*)fileinfo.filePath().toUtf8() << " opened." << std::endl;
    
  // Try to read .off in a polyhedron
  Scene_polyhedron_item* item = new Scene_polyhedron_item();
  item->setName(fileinfo.baseName());
  if(!item->load(in))
  {
    delete item;

    // Try to read .off in a polygon soup
    Scene_polygon_soup* soup_item = new Scene_polygon_soup;
    soup_item->setName(fileinfo.baseName());
    in.close();
    std::ifstream in2(fileinfo.filePath().toUtf8());
    if(!soup_item->load(in2)) 
    {
      delete soup_item;
      return 0;
    }
    return soup_item;
  }

  return item;
}

bool Io_off_plugin::canSave(const Scene_item* item)
{
  // This plugin supports polyhedrons and polygon soups
  return qobject_cast<const Scene_polyhedron_item*>(item) ||
         qobject_cast<const Scene_polygon_soup*>(item) ||
         qobject_cast<const Scene_starset3_item*>(item);
}

bool Io_off_plugin::save(const Scene_item* item, QFileInfo fileinfo)
{
  // This plugin supports polyhedrons and polygon soups
  const Scene_polyhedron_item* poly_item = 
    qobject_cast<const Scene_polyhedron_item*>(item);
  const Scene_polygon_soup* soup_item = 
    qobject_cast<const Scene_polygon_soup*>(item);
  const Scene_starset3_item* starset_item =
    qobject_cast<const Scene_starset3_item*>(item);

  if(!poly_item && !soup_item && !starset_item)
    return false;

  std::ofstream out(fileinfo.filePath().toUtf8());

  if(poly_item)
    poly_item->save(out);
  else if(soup_item)
    soup_item->save(out);
  else
    starset_item->save(out);

  return true;
}

#include <QtPlugin>
Q_EXPORT_PLUGIN2(Io_off_plugin, Io_off_plugin)
#include "Io_off_plugin.moc"
