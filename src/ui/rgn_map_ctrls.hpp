#pragma once

#include "../crosshatching/geometry.hpp"
#include "image_box.h"
#include "settingctrls.hpp"
#include <QWidget>
#include <unordered_set>

namespace ui {

    class flow_direction_panel : public image_box {

        Q_OBJECT

    public:
        flow_direction_panel();
    };

    class rgn_map_ctrl : public QWidget {

        Q_OBJECT

    private:

        using poly_tree_val = std::pair<int, const ch::polygon*>;
        struct poly_getter {
            const ch::polygon& operator()(const poly_tree_val& v) const {
                return *v.second;
            }
        };
        using poly_tree = ch::polygon_tree<poly_tree_val, poly_getter>;

        poly_tree tree_;
        std::optional<ch::point> curr_loc_;
        ch::dimensions<int> base_sz_;
        double scale_;
        std::unordered_set<int> selected_;
        std::vector<uchar> colors_;
        std::vector<ch::polygon> scaled_regions_;
        bool selecting_;
        int cursor_radius_;

        template<typename T>
        std::vector<int> find_intersecting_polys(const T& some_shape) const {
            namespace r = ranges;
            namespace rv = ranges::views;
            std::vector<poly_tree_val> results = tree_.query(some_shape);
            return results |
                rv::transform(
                    [](auto&& index_poly) {
                        return index_poly.first;
                    }
            ) | r::to_vector;
        }
 
        std::vector<int> polys_at_point(const ch::point& pt) const;
        ch::polygon cursor_poly(double scale, const ch::point& pt, int sz_index) const;

    protected:

        void paintEvent(QPaintEvent* event) override;
        void keyPressEvent(QKeyEvent* event) override;
        void mousePressEvent(QMouseEvent* event) override;
        void mouseMoveEvent(QMouseEvent* event) override;
        void mouseReleaseEvent(QMouseEvent* event) override;

    public:
        rgn_map_ctrl();
        void set_regions(vector_graphics_ptr gfx);
        void set_scale(double sc);
    };

    class main_window;

    class rgn_map_panel : public QWidget {

        Q_OBJECT

        void populate_ctrls();

    private:
        main_window* parent_;
        QComboBox* layer_cbo_;
        QCheckBox* brush_cbx_;
        QCheckBox* flow_cbx_;
        QComboBox* curr_brush_cbo_;
        flow_direction_panel* flow_ctrl_;
    
    public:
        rgn_map_panel(main_window* parent);
    };

}