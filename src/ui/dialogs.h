#pragma once

#include <QDialog>
#include <QLineEdit>
#include <QTextEdit>

class brush_dialog : public QDialog {

    Q_OBJECT

public:
    brush_dialog(QWidget* parent = 0);

private:
    QLineEdit* name_box_;
    QTextEdit* code_box_;
};