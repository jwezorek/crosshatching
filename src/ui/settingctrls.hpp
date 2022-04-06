#pragma once

#include <opencv2/core.hpp>
#include <QtWidgets>
#include <tuple>
#include <functional>
namespace ui {

    using val_to_pos_fn = std::function<int(double val, double range, double min, double max)>;
    using pos_to_val_fn = std::function<double(int val, double range, double min, double max)>;

    class labeled_slider : public QWidget {

        Q_OBJECT

    public:
        static int linear_value_to_position(double val, double range, double min, double max);
        static double linear_position_to_value(int pos, double range, double min, double max);

        labeled_slider(const QString& txt, double min, double max, double init_val, int range = 1000,
            val_to_pos_fn val_to_pos = linear_value_to_position, pos_to_val_fn pos_to_val = linear_position_to_value);
        void set(double value);
        double value() const;

    signals:
        void value_changed(double new_val);
        void slider_released();

    private:
        val_to_pos_fn val_to_pos_;
        pos_to_val_fn pos_to_val_;
        double min_;
        double max_;
        QSlider* slider_;
        QLabel* lbl_val_;

        double range() const;
        int value_to_position(double val) const;
        double position_to_value(int pos) const;

        void handle_position_change(int pos);
        void handle_released();

    };

    class preprocess_settings : public QWidget {
        Q_OBJECT
    private:
        labeled_slider* scale_slider_;
        labeled_slider* contrast_slider_;
        labeled_slider* thresh_slider_;

    public:
        
        cv::Mat scaled_image;
        cv::Mat contrast_applied_image;

        preprocess_settings();
        void initialize();
        double scale() const;
        double beta() const;
        double sigma() const;
        std::tuple<double, double> contrast_params() const;

    signals:
        void contrast_changed(std::tuple<double, double>);
        void scale_changed(double);
    };

    class shock_filter_settings : public QWidget {
        Q_OBJECT
    private:
        labeled_slider* sigma_slider_;
        labeled_slider* str_sigma_slider_;
        labeled_slider* blend_slider_;
        labeled_slider* iter_slider_;
    public:

        cv::Mat filtered_image;

        shock_filter_settings();
        void initialize();

    signals:
        void changed(std::tuple<double, double>);
    };

}