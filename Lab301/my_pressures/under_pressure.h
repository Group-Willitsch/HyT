#ifndef UNDER_PRESSURE_H
#define UNDER_PRESSURE_H

#include <QMainWindow>

namespace Ui {
class under_pressure;
}

class under_pressure : public QMainWindow
{
    Q_OBJECT

public:
    explicit under_pressure(QWidget *parent = 0);
    ~under_pressure();

private:
    Ui::under_pressure *ui;
};

#endif // UNDER_PRESSURE_H
