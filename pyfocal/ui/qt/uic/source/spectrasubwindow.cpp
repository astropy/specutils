#include "spectrasubwindow.h"
#include "ui_spectrasubwindow.h"

SpectraSubWindow::SpectraSubWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::SpectraSubWindow)
{
    ui->setupUi(this);
}

SpectraSubWindow::~SpectraSubWindow()
{
    delete ui;
}
