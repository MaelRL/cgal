//
// Author(s)     : Jane Tournois
//

#include <CGAL_demo/Plugin_interface.h>
#include <CGAL_demo/Plugin_helper.h>
#include "Implicit_surface_type.h"
#include "constrained_surfaces_implicit/Constrained_surface_implicit_interface.h"
#include "Scene_constrained_surface_implicit_item.h"
#include "ui_Function_dialog.h"

#include <iostream>

#include <QAction>
#include <QMainWindow>
#include <QPluginLoader>
#include <QDir>
#include <QApplication>
#include <QMenu>
#include <QList>
#include <QLibrary>

class Io_constrained_surface_implicit_plugin :
  public QObject,
  protected Plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(Plugin_interface)
  
public:
  Io_constrained_surface_implicit_plugin();
  virtual ~Io_constrained_surface_implicit_plugin() {}
  
  using Plugin_helper::init;
  virtual void init(QMainWindow* mainWindow, Scene_interface* scene_interface);
  
  QList<QAction*> actions() const
  {
    return QList<QAction*>();
  }

  bool applicable() const { return true; }
  
public slots:
  void load_function() const;
  
private:
  void load_function_plugins();
 
private:
  QList<Constrained_surface_implicit_interface*> functions_;
};



Io_constrained_surface_implicit_plugin::
Io_constrained_surface_implicit_plugin()
{
  load_function_plugins();
}


void
Io_constrained_surface_implicit_plugin::
init(QMainWindow* mainWindow, Scene_interface* scene_interface)
{
  this->scene = scene_interface;
  this->mw = mainWindow;
  
  QAction* actionLoad = new QAction("Load Constrained Implicit function (aniso)", mw);
  if( NULL != actionLoad )
    connect(actionLoad, SIGNAL(triggered()), this, SLOT(load_function()));
  
  QMenu* menuFile = mw->findChild<QMenu*>("menuFile");
  if ( NULL != menuFile )
  {
    QList<QAction*> menuFileActions = menuFile->actions();
    
    // Look for action just after "Load..." action
    QAction* actionAfterLoad = NULL;
    for ( QList<QAction*>::iterator it_action = menuFileActions.begin(), 
         end = menuFileActions.end() ; it_action != end ; ++ it_action ) //Q_FOREACH( QAction* action, menuFileActions)
    {
      if ( NULL != *it_action && (*it_action)->text().contains("Load") )
      {
        ++it_action;
        if ( it_action != end && NULL != *it_action )
        {
          actionAfterLoad = *it_action;
        }
      }
    }
    
    // Insert "Load implicit function" action
    if ( NULL != actionAfterLoad )
    {
      menuFile->insertAction(actionAfterLoad,actionLoad);      
    }   
  }
}

void
Io_constrained_surface_implicit_plugin::
load_function() const
{
  QDialog dialog(mw);
  Ui::FunctionDialog ui;
  ui.setupUi(&dialog);
  
  connect(ui.buttonBox, SIGNAL(accepted()), &dialog, SLOT(accept()));
  connect(ui.buttonBox, SIGNAL(rejected()), &dialog, SLOT(reject()));
  
  // Add loaded functions to the dialog
  int i=0;
  Q_FOREACH( Constrained_surface_implicit_interface* f, functions_ )
  {
    ui.functionList->insertItem(i++, f->name());
  }

  // Open window
  int return_code = dialog.exec();
  if(return_code == QDialog::Rejected) { return; }
  
  // Get selected function
  i = ui.functionList->currentIndex();
  Constrained_surface_implicit_interface* function = functions_[i];
  QString param = ui.parametersList->text();
  std::cout << param.toStdString() << std::endl;

  function->set_parameters(param.toStdString());//if param is valid

  // Create Scene_implicit_function object and add it to the framework
  Scene_constrained_surface_implicit_item* item =
    new Scene_constrained_surface_implicit_item(function->surface());
  
  item->setName(tr("%1").arg(function->name()));
  item->setRenderingMode(FlatPlusEdges);

  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  scene->itemChanged(index);
    
  Scene_interface::Item_id new_item_id = scene->addItem(item);
  scene->setSelectedItem(new_item_id);
}

void
Io_constrained_surface_implicit_plugin::
load_function_plugins()
{
  QDir pluginsDir(qApp->applicationDirPath());
  QString dirname = pluginsDir.dirName();
  if ( !pluginsDir.cd("constrained_surfaces_implicit") ) { 
    // In that case, dirname may be "Debug" or "Release" and one has to
    // search in ../implicit_functions/Debug or
    // ../implicit_functions/Release
    QString newDir = QString("../constrained_surfaces_implicit/") + dirname;
    if( !pluginsDir.cd(newDir) ) return; 
  }
  
  Q_FOREACH (QString fileName, pluginsDir.entryList(QDir::Files))
  {
    if( fileName.contains("plugin") && QLibrary::isLibrary(fileName) )
    {
      qDebug("    + Loading Function \"%s\"...", fileName.toUtf8().data());
      QPluginLoader loader;
      loader.setFileName(pluginsDir.absoluteFilePath(fileName));
      QObject *function_plugin = loader.instance();
      if ( NULL != function_plugin )
      {
        Constrained_surface_implicit_interface* function =
          qobject_cast<Constrained_surface_implicit_interface*>(function_plugin);
        if ( NULL != function )
        {
          functions_ << function;
        }
      }
    }
  }
}


#include <QtPlugin>
Q_EXPORT_PLUGIN2(Io_constrained_surface_implicit_plugin, Io_constrained_surface_implicit_plugin)
#include "Io_constrained_surface_implicit_plugin.moc"
