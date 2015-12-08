#ifndef SPECTRASUBWINDOW_H
#define SPECTRASUBWINDOW_H

#include <QMainWindow>

namespace Ui {
class SpectraSubWindow;
}

class SpectraSubWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit SpectraSubWindow(QWidget *parent = 0);
    ~SpectraSubWindow();

private:
    Ui::SpectraSubWindow *ui;
};

#endif // SPECTRASUBWINDOW_H
