#pragma once

#include "../crosshatching/brush_language.hpp"
#include "../crosshatching/drawing.hpp"
#include "cv_image_box.h"
#include "drawing_worker.h"
#include <opencv2/core.hpp>
#include <QDialog>
#include <QLineEdit>
#include <QTextEdit>
#include <QPushButton>
#include <QProgressBar>
#include <QPlainTextEdit>
#include <QDialogButtonBox>
#include <QComboBox>
#include <QLabel>
#include <QMenu>
#include <tuple>
#include <variant>
#include <any>
#include <string>
#include <QStackedWidget>
#include <memory>
#include <QThread>

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
        void update_btn_enabled_state();

        QLineEdit* name_box_;
        QTextEdit* code_box_;
        QDialogButtonBox* btns_;
        QPushButton* view_btn_;
        ch::brush_expr_ptr brush_;
    };

    using layer_param = std::variant<bool, std::tuple<std::string, double>>;
    class layer_dialog : public QDialog {

        Q_OBJECT

    public:
        layer_dialog(const std::vector<std::string>& brushes, const layer_param& current_layer_params);
        std::string brush_name() const;
        double value() const;

        static std::optional<std::tuple<std::string, double>> create_layer_item(const std::vector<std::string>& brushes, bool is_initial);
        static std::optional<std::tuple<std::string, double>> edit_layer_item(const std::vector<std::string>& brushes, const std::string& curr_brush, double curr_end_val);
    private:

        void init_create_layer_dialog(bool initial);
        void init_edit_layer_dialog(const std::string& brush, double val);
        void update_btn_enabled_state();

        QComboBox* brush_box_;
        QLineEdit* value_edit_;
        QDialogButtonBox* btns_;
    };

    class test_swatch_picker : public QDialog {

        Q_OBJECT

    public:
        test_swatch_picker(cv::Mat img);
        cv::Mat test_swatch() const;

        static cv::Mat get_test_swatch(cv::Mat src);

    private:
        void update_btn_enabled_state();

        cv_image_box* selector_;
        QDialogButtonBox* btns_;
        cv::Mat src_img_;
    };

    class drawing_progress : public QDialog {

        Q_OBJECT

    public:
        drawing_progress(const ch::crosshatching_job& job);
        void generate_drawing();
        void update_progress(double pcnt);
        void set_status_line(const std::string& str);
        void log_message(const std::string& str);
        void on_finished();
    private:
        ch::crosshatching_job job_;
        QLabel* status_;
        QProgressBar* progress_;
        QPlainTextEdit* log_;
        QPushButton* button_;
        QStackedWidget* btn_stack_;
        QPushButton* view_drawing_btn_;
        QPushButton* save_drawing_btn_;
        QPushButton* exit_btn_;
    };

    class progress : public QDialog {

        Q_OBJECT

    public:
        progress();
        void update_progress(double pcnt);
        std::function<void(double)> progress_func() const;
        std::any run(const std::string& title, std::function<std::any()> job);

    private:
        QProgressBar* progress_;
        std::unique_ptr<worker> worker_;
        std::unique_ptr<QThread> thread_;
    };
}