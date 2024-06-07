#include "MainWindow.h"
#include <QApplication>

int main(int argc, char *argv[])
{

    QApplication a(argc, argv);

    QFile qss(qApp->applicationDirPath() + "/Resources/QSS/soap.qss");

    if (qss.open(QFile::ReadOnly | QFile::Text)) {

        a.setStyleSheet(QLatin1String(qss.readAll()));

        qss.close();
    }

    MainWindow w;

    w.show();

    return a.exec();
}
