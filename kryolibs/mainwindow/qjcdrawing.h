/*****************************************************************************************
                            qjcdrawing.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef QJCDRAWING_H
#define QJCDRAWING_H

#include <vector>

#include <qwidget.h>
#include <qimage.h>
#include <qpicture.h>
//Added by qt3to4:
#include <QResizeEvent>
#include <QMouseEvent>
#include <QMenu>
#include <QPaintEvent>

#include "fidarray.h"
#include "kryomolcore_export.h"

struct   zoom
{
  zoom() { firstpoint=0; lastpoint=0; ratio=1.0f; offset=fleft=0.0f; fright=1.0f; }
  int npoints; int firstpoint; int lastpoint; float ratio; float offset;float fleft;float fright;

  void SetPoints(int npoints)
  {
    firstpoint=(int)(npoints*fleft);
    lastpoint=(int)(npoints*fright);
    if(lastpoint>(npoints-1))
      lastpoint=npoints-1;
  }

  void SetPoints(int npoints, int left,int right)
  {
    fleft=(float)left/(float)npoints-1;
    fright=(float)right/(float)npoints-1;
    SetPoints(npoints);
  }

  void reset(int size)
  {
    fleft=0.0f;
    fright=1.0f;
    firstpoint=0;
    lastpoint=size-1;
  }

} ;

class QRubberBand;
class QMenu;
class QwtScaleDraw;
class QwtLinearScaleEngine;
class KRYOMOLCORE_EXPORT QJCDrawing : public QWidget
{
  Q_OBJECT
public:
  QJCDrawing(QWidget *parent = 0);

  ~QJCDrawing();
  enum eVisualMode { eNORMAL, eZOOM };
  enum Orientation { LeftRight, RightLeft };
  enum BaseLinePosition { Top, Middle, Bottom };
  void Export(const QString& s, const QString& filter);
  void SetVisualMode(const eVisualMode& e) { m_visualmode=e; SetCursor(); }
  void SetOrientation(const QJCDrawing::Orientation& o) { m_orientation=o; }
  void SetBaseLinePosition (BaseLinePosition b );
  void SetZoom(float l, float r);
  void UpdateDrawing(bool autoscale=false);

protected:
  void UpdateDrawing(bool autoscale,const QRect& rect);
  void paintEvent(QPaintEvent* e);
  QPicture GetPicture();
  QImage GetImage();
  void mousePressEvent(QMouseEvent* e);
  void mouseMoveEvent(QMouseEvent* e);
  void mouseReleaseEvent(QMouseEvent* e);
  void resizeEvent(QResizeEvent* e);
  void SetCursor();
  virtual QSize sizeHint() const { return QSize(400,300); }
private:
  void Draw(QPainter* pDC, QRect rect, bool bprint=false);
  void DrawData(QPainter* p);
  void DrawPrintData(QPainter* p);
  fidarray* m_data;
  QRect m_rect;
  zoom m_zoom;
  QColor m_color;
  QColor m_bkcolor;
  double m_factor; //scaling factor;
  QwtScaleDraw* m_scale;
  QwtLinearScaleEngine* m_scaleengine;
  QPalette* m_scalepalette;
  float m_min;
  float m_max;
  float m_shift;
  int m_npoints;
  eVisualMode m_visualmode;
  QRubberBand* m_tracker;
  QMenu* m_optionsmenu;
  QJCDrawing::Orientation m_orientation;
  BaseLinePosition m_baselineposition;
  std::vector<QPoint> m_points;
public slots:
  void SetData(fidarray* data,float max,float min, float shift, int npoints);
  void OnIncrease();
  void OnDecrease();
  void OnZoom(bool);
  void ResetZoom();
  void OnCopy();
  void OnPrint();
  void OnSnapshot();
protected slots:
  void OnSetLimits();
signals:
  void limits(float,float);
};

#endif
