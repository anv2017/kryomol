#include "orcadialog.h"

#include <QPushButton>
#include <QGridLayout>
#include <QTextEdit>
#include <QTextStream>
#include <QFileDialog>
#include <QTemporaryFile>
#include <QMessageBox>

#include "molecule.h"
#include "atom.h"

OrcaDialog::OrcaDialog(QWidget* parent) : QDialog(parent)
{
    m_layout = new QGridLayout(this);
    m_textedit = new QTextEdit();
    m_layout->addWidget( m_textedit );
    m_runbutton = new QPushButton(this);
    m_savebutton = new QPushButton(this);
    m_loadbutton = new QPushButton(this);
    m_runbutton->setText("Run");
    m_savebutton->setText("Save to disk");
    m_loadbutton->setText("Load into Kryomol");
    m_savebutton->setEnabled(false);
    m_loadbutton->setEnabled(false);
    m_layout->addWidget(m_runbutton);
    m_layout->addWidget(m_savebutton);
    m_layout->addWidget(m_loadbutton);
    connect(m_runbutton,SIGNAL(clicked()),this,SLOT(OnRunButton()));
    connect(m_savebutton,SIGNAL(clicked()),this,SLOT(OnSaveButton()));
    connect(m_loadbutton,SIGNAL(clicked()),this,SLOT(OnLoadButton()));
    connect(&m_engine,SIGNAL(finished(int )),this,SLOT(OnJobFinished(int )));

}

void OrcaDialog::SetMolecule(const kryomol::Molecule *mol)
{
    m_molecule=mol;
    InitTextEdit();
}

void OrcaDialog::InitTextEdit()
{
    QString str;
    QTextStream stream(&str);
    //stream << "!PAL2" << endl;
    stream << "!HF-3c" << endl;
    stream << "*xyz 0 1 " << endl;
    const std::vector<kryomol::Atom>& atoms=m_molecule->Atoms();
    const std::vector<kryomol::Coordinate>& xyz=m_molecule->CurrentFrame().XYZ();
    size_t natoms=atoms.size();
    for(size_t i=0;i<natoms;++i)
    {
        stream << QString::fromStdString(atoms[i].Symbol()) << "  " << xyz.at(i).x() << "  " << xyz.at(i).y() << "  " << xyz.at(i).z() <<endl;
    }
    stream << "*" << endl;
    m_textedit->setText(str);

    this->update();

}

void OrcaDialog::OnRunButton()
{
    m_engine.SetInputFile(m_textedit->toPlainText());
    m_engine.Run();
    this->setCursor(QCursor(Qt::WaitCursor));

}

void OrcaDialog::OnSaveButton()
{
    QString filename=QFileDialog::getSaveFileName(this);
    QFile ofile(filename);
    ofile.open(QFile::WriteOnly);
    QTextStream ostream(&ofile);
    ostream << m_textedit->toPlainText();
    ofile.close();
}

void OrcaDialog::OnLoadButton()
{
    QTemporaryFile* tfile = new QTemporaryFile();
    tfile->open();
    QTextStream tstream(tfile);
    tstream << m_textedit->toPlainText();
    qDebug() << "tempfile=" << tfile->fileName();
    emit outputfile(tfile->fileName());

}

void OrcaDialog::OnJobFinished(int exitcode)
{
    /*QDialog*  d= new QDialog(this);
    //QGridLayout* lay = new QGridLayout(d);
    QTextEdit* textedit= new QTextEdit(d);
    m_layout->addWidget(textedit);*/
    //QFile ofile(m_engine.OutputFileName());
    //ofile.open(QFile::ReadOnly);
    QMessageBox::information(this,"","Orca job finished");
    this->unsetCursor();
    m_textedit->setText(m_engine.OutputString());
    m_savebutton->setEnabled(true);
    m_loadbutton->setEnabled(true);

    //emit outputfile(m_engine.OutputFileName());
    //this->close();
}


