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
    connect(name_box_, &QLineEdit::textChanged, this, &brush_dialog::update_btn_enabled_state);

    update_btn_enabled_state();
}

void ui::brush_dialog::parse_brush_code() {
    auto result = ch::brush_language_to_expr(code_box_->toPlainText().toStdString());
    if (std::holds_alternative<std::string>(result)) {
        QMessageBox mb;
        mb.setText(std::get<std::string>(result).c_str());
        mb.exec();
    } else {
        brush_ = std::get<ch::brush_expr_ptr>(result);
    }
    update_btn_enabled_state();
}

void ui::brush_dialog::launch_brush_viewer() {
    auto func = std::get<ch::brush_fn>(brush_->eval());
    auto viewer = new brush_viewer(name_box_->text().toStdString(), func, nullptr);
    viewer->exec();
}

void ui::brush_dialog::update_btn_enabled_state() {
    bool has_valid_brush = brush_ != nullptr;
    bool has_valid_name = !name_box_->text().isEmpty(); // TODO: need to test for name collisions.
    view_btn_->setEnabled(has_valid_brush);
    btns_->button(QDialogButtonBox::Ok)->setEnabled(has_valid_brush && has_valid_name);
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

ui::layer_dialog::layer_dialog(const std::vector<std::string>& brushes, const ui::layer_param& current_layer_params) :
        QDialog(nullptr) {
    
    auto layout = new QVBoxLayout(this);
    layout->addWidget(new QLabel("Brush"));
    layout->addWidget(brush_box_ = new QComboBox());
    layout->addWidget(new QLabel("End of layer value"));
    layout->addWidget(value_edit_ = new QLineEdit());

    auto validator = new QDoubleValidator(0, 1.0, 2, this);
    validator->setNotation(QDoubleValidator::StandardNotation);
    value_edit_->setValidator(validator);

    for (const auto& str : brushes) {
        brush_box_->addItem(str.c_str());
    }
    layout->addWidget(btns_ = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel));
    connect(btns_, &QDialogButtonBox::accepted, this, &QDialog::accept);
    connect(btns_, &QDialogButtonBox::rejected, this, &QDialog::reject);
    connect(value_edit_, &QLineEdit::textChanged, this, &layer_dialog::update_btn_enabled_state);

    if (std::holds_alternative<bool>(current_layer_params)) {
        init_create_layer_dialog(std::get<bool>(current_layer_params));
    } else {
        const auto& [brush, val] = std::get<std::tuple<std::string, double>>(current_layer_params);
        init_edit_layer_dialog(brush, val);
    }
}

std::string ui::layer_dialog::brush_name() const {
    return brush_box_->currentText().toStdString();
}

double ui::layer_dialog::value() const {
    return value_edit_->text().toDouble();
}

std::optional<std::tuple<std::string, double>> ui::layer_dialog::create_layer_item(const std::vector<std::string>& brushes, bool is_initial)
{
    std::unique_ptr<ui::layer_dialog> dlg = std::make_unique<ui::layer_dialog>(brushes, layer_param{ is_initial });
    if (dlg->exec() == QDialog::Accepted) {
        return { { dlg->brush_name(), dlg->value() } };
    } else {
        return {};
    }
}

std::optional<std::tuple<std::string, double>> ui::layer_dialog::edit_layer_item(
        const std::vector<std::string>& brushes, const std::string& curr_brush, double curr_end_val) {

    std::unique_ptr<ui::layer_dialog> dlg = std::make_unique<ui::layer_dialog>(brushes, std::tuple(curr_brush, curr_end_val));
    if (dlg->exec() == QDialog::Accepted) {
        return { { dlg->brush_name(), dlg->value() } };
    } else {
        return {};
    }
}

void ui::layer_dialog::init_create_layer_dialog(bool initial) {
    setWindowTitle("Create crosshatching layer...");
    if (initial) {
        value_edit_->setText("1.0");
        value_edit_->setEnabled(false);
    }
}

void ui::layer_dialog::init_edit_layer_dialog(const std::string& brush, double val) {
    setWindowTitle("Edit crosshatching layer...");
    int curr_item = brush_box_->findText(brush.c_str());
    if (curr_item == -1) {
        curr_item = 0; // TODO: this is actually a bug if it happens
    }
    brush_box_->setCurrentIndex(curr_item);
    value_edit_->setText( std::to_string(val).c_str() );
}

void ui::layer_dialog::update_btn_enabled_state() {
    auto val = value();
    bool valid_value = val >= 0.0 && val <= 1.0;
    btns_->button(QDialogButtonBox::Ok)->setEnabled(valid_value);
}

ui::rect_in_image_selector::rect_in_image_selector(QWidget* parent) :
        QLabel(parent) {
    selectionStarted = false;
}

void ui::rect_in_image_selector::paintEvent(QPaintEvent* e)
{
    QLabel::paintEvent(e);

    QPainter painter(this);
    QBrush brush(QColor(0, 0, 0, 180));
    QPen pen(brush, 1, Qt::DashLine);
    painter.setPen(pen);
    painter.setBrush(QBrush(QColor(255, 255, 255, 120)));

    painter.drawRect(selectionRect);
}

void  ui::rect_in_image_selector::mousePressEvent(QMouseEvent* e) {
    selectionStarted = true;
    selectionRect.setTopLeft(e->pos());
    selectionRect.setBottomRight(e->pos());
}

void ui::rect_in_image_selector::mouseMoveEvent(QMouseEvent* e)
{
    if (selectionStarted)
    {
        selectionRect.setBottomRight(e->pos());
        repaint();
    }
}

void ui::rect_in_image_selector::mouseReleaseEvent(QMouseEvent* e)
{
    selectionStarted = false;
}

ui::test_swatch_picker::test_swatch_picker(cv::Mat img) :
        src_img_(img) {
    rect_in_image_selector* ss = new rect_in_image_selector(this);
    ss->setPixmap(QPixmap::fromImage(QImage(img.data, img.cols, img.rows, img.step, QImage::Format_BGR888)));
    ss->show();
}

cv::Mat  ui::test_swatch_picker::test_swatch() const {
    return test_swatch_;
}

cv::Mat ui::test_swatch_picker::get_test_swatch(cv::Mat src_img) {
    auto dlg = std::make_unique<ui::test_swatch_picker>(src_img);
    if (dlg->exec() == QDialog::Accepted) {
        return dlg->test_swatch();
    } else {
        return {};
    }
}
