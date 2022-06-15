#include "dialogs.h"
#include "brush_viewer.h"
#include <QtWidgets>

namespace {
    QWidget* spacer(int n) {
        auto spacer = new QWidget();
        spacer->setFixedSize(QSize(n, n));
        return spacer;
    }

    QMainWindow* get_mainwindow() {
        foreach(QWidget * w, qApp->topLevelWidgets())
            if (QMainWindow* mainWin = qobject_cast<QMainWindow*>(w))
                return mainWin;
        return nullptr;
    }

    void initialize_dlg(QDialog* dlg, const std::string& title) {
        dlg->setWindowTitle(title.c_str());
        auto main_wnd = get_mainwindow();
        QSize main_wnd_sz = main_wnd->size();
        dlg->resize(0.7 * main_wnd_sz);
    }
}

ui::brush_dialog::brush_dialog(QWidget* parent) :
        QDialog(parent) {

    initialize_dlg( this, "New brush...");

    auto layout = new QVBoxLayout(this);
    layout->addWidget( new QLabel("name") );
    layout->addWidget( name_box_ = new QLineEdit() );
    layout->addWidget(spacer(5));
    layout->addWidget( new QLabel("code") );
    layout->addWidget( code_box_ = new QTextEdit() );
    QPushButton* parse_btn;
    layout->addWidget(parse_btn = new QPushButton("Parse"));
    layout->addWidget(spacer(10));
    layout->addWidget( btns_ = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel) );

    view_btn_ = new QPushButton("View");
    btns_->addButton(view_btn_, QDialogButtonBox::ActionRole);

    connect(btns_, &QDialogButtonBox::accepted, this, &QDialog::accept);
    connect(btns_, &QDialogButtonBox::rejected, this, &QDialog::reject);
    connect(parse_btn, &QPushButton::clicked, this, &brush_dialog::parse_brush_code);
    connect(view_btn_, &QPushButton::clicked, this, &brush_dialog::launch_brush_viewer);

    btns_->button(QDialogButtonBox::Ok)->setEnabled(false);
    view_btn_->setEnabled(false);
}

void ui::brush_dialog::parse_brush_code() {
    auto result = ch::brush_language_to_expr(code_box_->toPlainText().toStdString());
    if (std::holds_alternative<std::string>(result)) {
        QMessageBox mb;
        mb.setText(std::get<std::string>(result).c_str());
        mb.exec();
        btns_->button(QDialogButtonBox::Ok)->setEnabled(false);
        view_btn_->setEnabled(false);
    } else {
        brush_ = std::get<ch::brush_expr_ptr>(result);
        btns_->button(QDialogButtonBox::Ok)->setEnabled(true);
        view_btn_->setEnabled(true);
    }
}

void ui::brush_dialog::launch_brush_viewer() {
    auto func = std::get<ch::brush_fn>(brush_->eval());
    auto viewer = new brush_viewer(name_box_->text().toStdString(), func, nullptr);
    viewer->exec();
}

std::string ui::brush_dialog::brush_name() const {
    return name_box_->text().toStdString();
}

ch::brush_expr_ptr ui::brush_dialog::brush_expr() const {
    return brush_;
}

std::optional<std::tuple<std::string, ch::brush_expr_ptr>> ui::brush_dialog::create_brush()
{
    std::unique_ptr<ui::brush_dialog> dlg = std::make_unique<ui::brush_dialog>();
    if (dlg->exec() == QDialog::Accepted) {
        return { { dlg->brush_name(), dlg->brush_expr() } };
    } else {
        return {};
    }
}


ch::brush_expr_ptr ui::brush_dialog::edit_brush(const std::string& name, const std::string& code) {
    return {};
}

ui::add_layer_dialog::add_layer_dialog(const std::vector<std::string>& brushes) : 
        QDialog(nullptr) {
    setWindowTitle("Add layer...");
    auto layout = new QVBoxLayout(this);
    layout->addWidget(new QLabel("Brush"));
    layout->addWidget(brush_box_ = new QComboBox());
    layout->addWidget(new QLabel("End of layer value"));
    layout->addWidget(value_edit_ = new QLineEdit());
    value_edit_->setValidator(new QDoubleValidator(0, 100, 2, this));

    for (const auto& str : brushes) {
        brush_box_->addItem(str.c_str());
    }
    QDialogButtonBox* btns;
    layout->addWidget(btns = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel));
    connect(btns, &QDialogButtonBox::accepted, this, &QDialog::accept);
    connect(btns, &QDialogButtonBox::rejected, this, &QDialog::reject);
}

std::string ui::add_layer_dialog::brush_name() const {
    return brush_box_->currentText().toStdString();
}

double ui::add_layer_dialog::value() const {
    return value_edit_->text().toDouble();
}

std::optional<std::tuple<std::string, double>> ui::add_layer_dialog::create_layer_item(const std::vector<std::string>& brushes, bool is_initial_layer)
{
    std::unique_ptr<ui::add_layer_dialog> dlg = std::make_unique<ui::add_layer_dialog>(brushes);
    if (dlg->exec() == QDialog::Accepted) {
        return { { dlg->brush_name(), dlg->value() } };
    } else {
        return {};
    }
}
