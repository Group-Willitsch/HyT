#include "under_pressure.h"
#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    under_pressure w;
    w.show();

    return a.exec();
}
