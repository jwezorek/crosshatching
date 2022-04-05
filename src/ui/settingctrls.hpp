#pragma once

#include <QtWidgets>
#include <tuple>

namespace ui {

    class labeled_slider : public QWidget {

        Q_OBJECT

    public:
        labeled_slider(const QString& txt, double min, double max, double init_val, int range = 100);
        void set(double value);
        double value() const;

    signals:
        void value_changed(double new_val);
        void slider_released();

    private:
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
}