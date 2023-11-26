#ifndef ORCADIALOG_H
#define ORCADIALOG_H

#include "orcaengine.h"
#include <QDialog>

namespace kryomol {
class Molecule;
}
class QTextEdit;
class QGridLayout;
class QPushButton;
class OrcaDialog : public QDialog
{
    Q_OBJECT
public:
    OrcaDialog(QWidget* parent=nullptr);
    void SetMolecule(const kryomol::Molecule* mol);
public slots:
    void OnJobFinished(int exitcode);
signals:
    void outputfile(QString );
private:
    void InitTextEdit();
private slots:
    void OnRunButton();
    void OnSaveButton();
    void OnLoadButton();
private:
    OrcaEngine m_engine;
    const kryomol::Molecule* m_molecule;
    QTextEdit* m_textedit;
    QGridLayout* m_layout;
    QPushButton* m_runbutton;
    QPushButton* m_savebutton;
    QPushButton* m_loadbutton;
};

#endif // ORCADIALOG_H
