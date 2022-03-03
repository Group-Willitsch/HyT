#include "under_pressure.h"
#include "ui_under_pressure.h"

under_pressure::under_pressure(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::under_pressure)
{
    ui->setupUi(this);
}

under_pressure::~under_pressure()
{
    delete ui;
}
