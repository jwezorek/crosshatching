#pragma once

#include <opencv2/core.hpp>
#include <QtWidgets>
#include <tuple>
#include <functional>
namespace ui {

    using val_to_pos_fn = std::function<int(double val, double range, double min, double max)>;
    using pos_to_val_fn = std::function<double(int val, double range, double min, double max)>;

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

    class image_processing_pipeline_item : public QWidget {
        Q_OBJECT
    protected:
        int index_;
        QWidget* content_area_;
        cv::Mat output_;
    public:
        image_processing_pipeline_item(QString title, int index);
        int index() const;
        cv::Mat output() const;
        void clear_output();
        virtual void populate();
        virtual void initialize();
        virtual cv::Mat process_image(cv::Mat input);
        virtual bool is_on() const;
    signals:
        void changed(int index);
    };

    class scale_and_contrast : public image_processing_pipeline_item {
        Q_OBJECT
    private:
        float_value_slider* scale_slider_;
        float_value_slider* contrast_slider_;
        float_value_slider* thresh_slider_;

    public:
        scale_and_contrast();
        void populate() override;
        void initialize() override;
        cv::Mat process_image(cv::Mat input) override;
        bool is_on() const override;

    };

    class anisotropic_diffusion_filter : public image_processing_pipeline_item {
        Q_OBJECT
    private:
        float_value_slider* alpha_slider_;
        float_value_slider* k_slider_;
        int_value_slider* iters_slider_;

    public:
        anisotropic_diffusion_filter();
        void populate() override;
        void initialize() override;
        cv::Mat process_image(cv::Mat input) override;
        bool is_on() const override;
    };

    class shock_filter : public image_processing_pipeline_item {
        Q_OBJECT
    private:
        int_value_slider* sigma_slider_;
        int_value_slider* str_sigma_slider_;
        float_value_slider* blend_slider_;
        int_value_slider* iter_slider_;
    public:
        shock_filter();
        void populate() override;
        void initialize() override;
        cv::Mat process_image(cv::Mat input) override;
        bool is_on() const override;
    };

    class mean_shift_segmentation : public image_processing_pipeline_item {
        Q_OBJECT
    private:
        int_value_slider* sigmaS_slider_;
        float_value_slider* sigmaR_slider_;
        int_value_slider* min_size_slider_;
    public: 
        mean_shift_segmentation();
        void populate() override;
        void initialize() override;
        cv::Mat process_image(cv::Mat input) override;
        bool is_on() const override;
    };
}
