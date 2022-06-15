#pragma once

#include "../crosshatching/brush_language.hpp"
#include <QDialog>
#include <QLineEdit>
#include <QTextEdit>
#include <QPushButton>
#include <QDialogButtonBox>
#include <QComboBox>
#include <tuple>
#include <string>

namespace ui {

    class brush_dialog : public QDialog {

        Q_OBJECT

    public:
        brush_dialog(QWidget* parent = 0);

        std::string brush_name() const;
        ch::brush_expr_ptr brush_expr() const;

        static std::optional<std::tuple<std::string, ch::brush_expr_ptr>> create_brush();
        static ch::brush_expr_ptr edit_brush(const std::string& name, const std::string& code);

    private:

        void parse_brush_code();
        void launch_brush_viewer();

        QLineEdit* name_box_;
        QTextEdit* code_box_;
        QDialogButtonBox* btns_;
        QPushButton* view_btn_;
        ch::brush_expr_ptr brush_;
    };

    class add_layer_dialog : public QDialog {

        Q_OBJECT

    public:
        add_layer_dialog(const std::vector<std::string>& brushes);
        std::string brush_name() const;
        double value() const;

        static std::optional<std::tuple<std::string, double>> create_layer_item(const std::vector<std::string>& brushes, bool is_initial_layer);
    private:
        QComboBox* brush_box_;
        QLineEdit* value_edit_;
    };
}