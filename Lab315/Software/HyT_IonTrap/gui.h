#ifndef GUI_H
#define GUI_H

#include <QMainWindow>

namespace Ui {
class gui;
}

class gui : public QMainWindow
{
    Q_OBJECT

public:
    explicit gui(QWidget *parent = 0);
    ~gui();

private:
    Ui::gui *ui;
};

#endif // GUI_H
