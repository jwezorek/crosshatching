#pragma once

#include <QtWidgets>

namespace ui {
    using val_to_pos_fn = std::function<int(double val, double range, double min, double max)>;
    using pos_to_val_fn = std::function<double(int val, double range, double min, double max)>;

    int exp_val_to_pos(double val, double range, double min, double max);
    double exp_pos_to_val(int pos, double range, double min, double max);

    class float_value_slider : public QWidget {

        Q_OBJECT

    public:
        static int linear_value_to_position(double val, double range, double min, double max);
        static double linear_position_to_value(int pos, double range, double min, double max);

        float_value_slider(const QString& txt, double min, double max, double init_val, int range = 1000,
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

    class int_value_slider : public QWidget {

        Q_OBJECT

    private:
        QSlider* slider_;
        QLabel* lbl_val_;

        void handle_position_change(int pos);
        void handle_released();
    public:
        int_value_slider(const QString& txt, int min, int max, int init_val);
        void set(int val);
        int value() const;
    signals:
        void value_changed(int new_val);
        void slider_released();
    };

}
