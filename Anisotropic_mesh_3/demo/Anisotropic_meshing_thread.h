//******************************************************************************
// File Description : 
//******************************************************************************

#ifndef CGAL_DEMO_ANISO_MESH_3_MESHING_THREAD_H
#define CGAL_DEMO_ANISO_MESH_3_MESHING_THREAD_H

#include <QThread>
#include <QObject>
#include <QStringList>
#include <QString>
#include <QTimer>


class Scene_starset3_item;

//enum Meshing_Dimension {SURFACE = 2, VOLUME = 3, SURFACE_AND_VOLUME = 5};

class Anisotropic_mesh_function_interface
{
public:
  virtual ~Anisotropic_mesh_function_interface() {}
  
  // Launch
  virtual void launch() = 0;
  
  // Stop
  virtual void stop() = 0;
  
  // Logs
  virtual QStringList parameters_log() const = 0;
  virtual QString status(double time_period) const = 0;
};


class Anisotropic_meshing_thread : public QThread
{
  Q_OBJECT
public:
  // Constructor / Destructor
  Anisotropic_meshing_thread(Anisotropic_mesh_function_interface* f, 
                             Scene_starset3_item* item);
  virtual ~Anisotropic_meshing_thread();


public:
  // Scene item
  Scene_starset3_item* item() const { return item_; }
  
  // Infos about meshing
  double time() const { return time_; }
  
  // Logs
  QStringList parameters_log() const { return f_->parameters_log(); }
  
public slots:
  // Stop
  void stop();
  
private slots:
  // emit signal status report
  void emit_status();
  
signals:
  // Emitted at the end of the process
  void done(Anisotropic_meshing_thread*);
  // Informs about status of meshing
  void status_report(QString);
  
protected:
  // Overload of QThread function
  virtual void run();
  
private:
  Anisotropic_mesh_function_interface* f_;
  Scene_starset3_item* item_;
  double time_; // in seconds
  QTimer* timer_;
  double timer_period_;
};

#endif // CGAL_DEMO_ANISO_MESH_3_MESHING_THREAD_H
