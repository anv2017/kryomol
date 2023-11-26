#ifndef KRYOVISOROPTICAL_H
#define KRYOVISOROPTICAL_H

#include "kryovisor.h"

class QColor;

namespace kryomol
{
class KRYOMOLCORE_EXPORT KryoVisorOptical : public KryoVisor
{
    Q_OBJECT
public:
    KryoVisorOptical(World* world,QWidget* parent=0, const QGLWidget* shareWidget=0, Qt::WindowFlags f=0 );

public slots:
    void OnShowElectricDipole(bool );
    void OnShowMagneticDipole(bool );
    void OnShowVelocityDipole(bool );
    void RenderMolecule(size_t index,GLenum mode=GL_RENDER);
    void RenderDipole(const kryomol::Coordinate& c, const QColor& color);
    void OnSetActiveTransition(int transition) {m_activetransition = transition;}

protected:
    void RenderScreenText();

private:
    bool m_bshowelectricdipole;
    bool m_bshowmagneticdipole;
    bool m_bshowvelocitydipole;
    bool m_bshowtransitionchanges;
    int m_activetransition;

};
}

#endif // KRYOVISOROPTICAL_H
