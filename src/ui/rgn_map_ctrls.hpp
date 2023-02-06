#pragma once

#include "../crosshatching/geometry.hpp"
#include "../crosshatching/brush_lang.hpp"
#include "../crosshatching/ink_layers.hpp"
#include "image_box.h"
#include "settingctrls.hpp"
#include <QWidget>
#include <unordered_set>
#include <unordered_map>

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
        std::vector<ch::color> colors_;
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
        void set_regions(const ch::dimensions<int>& sz, ch::ink_layer* layer);
        void set_scale(double sc);
        const std::unordered_set<int>& selected() const;

    signals:
        void selection_changed();
    };

    class main_window;

    class rgn_map_panel : public QWidget {

        Q_OBJECT

    private:
        main_window* parent_;
        QComboBox* layer_cbo_;
        QCheckBox* brush_cbx_;
        QCheckBox* flow_cbx_;
        QComboBox* curr_brush_cbo_;
        flow_direction_panel* flow_ctrl_;
        QStackedWidget* stack_;
 
        std::unordered_map<std::string, ch::brush_expr_ptr> brush_to_name_;
        std::unordered_map<ch::brush_expr*, std::string> name_to_brush_;
        std::vector<ch::brush_expr_ptr> default_brushes_;

        std::vector<ch::brush_expr_ptr> get_brush_defaults() const;
        void handle_selection_change();

    public:
        rgn_map_panel(main_window* parent, QStackedWidget* stack);
        void repopulate_ctrls();
        void set_layers(ch::ink_layers* layers);
    };

    struct rgn_map_tools {
        QStackedWidget* rgn_map_stack;
        rgn_map_panel* rgn_props;

        void set_rgn_map_scale(double scale);
    };
}