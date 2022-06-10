#include "dialogs.h"
#include <QtWidgets>

namespace {
    QWidget* spacer() {
        auto spacer = new QWidget();
        spacer->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Preferred);
        return spacer;
    }
}

brush_dialog::brush_dialog(QWidget* parent) :
        QDialog(parent) {
    this->setWindowTitle("New brush...");
    auto layout = new QVBoxLayout(this);
    layout->addWidget( new QLabel("name") );
    layout->addWidget( name_box_ = new QLineEdit() );
    layout->addWidget( new QLabel("code") );
    layout->addWidget( code_box_ = new QTextEdit() );
    QDialogButtonBox* btns;
    layout->addWidget( btns = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel) );

    QPushButton* parse_btn = new QPushButton("View");
    btns->addButton(parse_btn, QDialogButtonBox::ActionRole);

    connect(btns, &QDialogButtonBox::accepted, this, &QDialog::accept);
    connect(btns, &QDialogButtonBox::rejected, this, &QDialog::reject);
}