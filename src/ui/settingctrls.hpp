#pragma once

#include <QtWidgets>
#include "float_value_slider.h"
#include <opencv2/core.hpp>
#include <tuple>
#include <functional>

namespace ui {

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
