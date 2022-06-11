#pragma once

#include "../crosshatching/brush_language.hpp"
#include <QDialog>
#include <QLineEdit>
#include <QTextEdit>
#include <QPushButton>
#include <QDialogButtonBox>
#include <tuple>
#include <string>

namespace ui {

    class brush_dialog : public QDialog {

        Q_OBJECT

    public:
        brush_dialog(QWidget* parent = 0);

        std::string brush_name() const;
        ch::brush_expr_ptr brush_expr() const;

        static std::tuple<std::string, ch::brush_expr_ptr> create_brush();
        static ch::brush_expr_ptr edit_brush(const std::string& name, const std::string& code);

    private:

        void parse_brush_code();

        QLineEdit* name_box_;
        QTextEdit* code_box_;
        QDialogButtonBox* btns_;
        QPushButton* view_btn_;
        ch::brush_expr_ptr brush_;
    };

}