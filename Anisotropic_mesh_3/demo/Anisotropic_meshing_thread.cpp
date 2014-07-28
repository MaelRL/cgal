//
// Author(s)     : Jane Tournois
//
//******************************************************************************
// File Description : 
//******************************************************************************


#include <QTime>
#include <QApplication>

#include "Anisotropic_meshing_thread.h"
#include "Scene_starset3_item.h"

Anisotropic_meshing_thread::
Anisotropic_meshing_thread(Anisotropic_mesh_function_interface* f, 
                           Scene_starset3_item* item)
  : f_(f)
  , item_(item)
  , time_(0)
  , timer_(new QTimer(this))
  , timer_period_(1)
{
  connect(timer_, SIGNAL(timeout()),
          this,   SLOT(emit_status()));
  
  timer_->start(static_cast<int>(timer_period_*1000));
  item_->moveToThread(this);
}


Anisotropic_meshing_thread::
~Anisotropic_meshing_thread()
{
  delete f_;
  delete timer_;
  QApplication::restoreOverrideCursor();
}


void
Anisotropic_meshing_thread::
run()
{
  QTime timer;
  timer.start();
  
  f_->launch();
  time_ = double(timer.elapsed()) / 1000;
  item_->moveToThread(QApplication::instance()->thread());
  emit done(this);
}


void
Anisotropic_meshing_thread::
stop()
{
  f_->stop();
  
  // Block application until thread is deleted
  QApplication::setOverrideCursor(Qt::WaitCursor);
}


void
Anisotropic_meshing_thread::
emit_status()
{
  emit (status_report(f_->status(timer_period_)));
}


#include "Anisotropic_meshing_thread.moc"
