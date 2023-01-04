#pragma once

#include <QtWidgets>

#include "float_value_slider.h"
#include "../crosshatching/geometry.hpp"
#include "../crosshatching/raster_to_vector.hpp"
#include "../crosshatching/util.hpp"
#include <variant>
#include <memory>
#include <opencv2/core.hpp>
#include <tuple>
#include <functional>

namespace ui {

    struct vector_graphics {
        std::vector<ch::colored_polygon> polygons;
        ch::dimensions<int> sz;
    };

    using vector_graphics_ptr = std::shared_ptr<vector_graphics>;
    using pipeline_output = std::variant<std::monostate, cv::Mat, vector_graphics_ptr>;

    class image_processing_pipeline_item : public QWidget {
        Q_OBJECT
    protected:
        int index_;
        QWidget* content_area_;
        pipeline_output output_;
    public:
        image_processing_pipeline_item(QString title, int index);
        int index() const;
        pipeline_output output() const;
        void clear_output();
        virtual void populate();
        virtual void initialize();
        virtual pipeline_output process_image(pipeline_output input);
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
        pipeline_output process_image(pipeline_output input) override;
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
        pipeline_output process_image(pipeline_output input) override;
        bool is_on() const override;
    };

    class edge_preserving_filter : public image_processing_pipeline_item {
        Q_OBJECT
    private:
        float_value_slider* sigma_s_;
        float_value_slider* sigma_r_;
        QComboBox* flag_;
    public:
        edge_preserving_filter();
        void populate() override;
        void initialize() override;
        pipeline_output process_image(pipeline_output input) override;
        bool is_on() const override;
    };

    class stylization_filter : public image_processing_pipeline_item {
        Q_OBJECT
    private:
        float_value_slider* sigma_s_;
        float_value_slider* sigma_r_;

    public:
        stylization_filter();
        void populate() override;
        void initialize() override;
        pipeline_output process_image(pipeline_output input) override;
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
        pipeline_output process_image(pipeline_output input) override;
        bool is_on() const override;
    };

    class mean_shift_segmentation : public image_processing_pipeline_item {
        Q_OBJECT
    private:
        int_value_slider* sigmaS_slider_;
        float_value_slider* sigmaR_slider_;
        int_value_slider* min_size_slider_;
        cv::Mat label_image_;
    public: 
        mean_shift_segmentation();
        void populate() override;
        void initialize() override;
        pipeline_output process_image(pipeline_output input) override;
        cv::Mat labels() const;
        bool is_on() const override;
    };

    class raster_to_vector : public image_processing_pipeline_item {
        Q_OBJECT
    private:
        float_value_slider* param_slider_;
    public:
        raster_to_vector();
        void populate() override;
        void initialize() override;
        pipeline_output process_image(pipeline_output input) override;
        bool is_on() const override;
    };
}
