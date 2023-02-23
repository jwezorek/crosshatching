#pragma once

#include <QWidget>
#include <QtWidgets>
#include <unordered_set>
#include <span>
#include "../crosshatching/ink_layers.hpp"
#include <range/v3/all.hpp>

/*------------------------------------------------------------------------------------------------*/

namespace ui {

    class rgn_selection {
        std::unordered_set<int> old_selection_;
        std::unordered_set<int> selection_;
    public:
        rgn_selection();

        const std::unordered_set<int>& selected_ids() const;
        bool has_selection() const;
        void clear();
        void reselect();
        bool is_selected_id(int id) const;
        void insert(std::span<int> ids);
        void remove(std::span<int> ids);

        auto selected_layer_items(ch::ink_layer* layer) const {
            namespace r = ranges;
            namespace rv = ranges::views;
            return selection_ |
                rv::transform(
                    [layer](auto index)->ch::ink_layer_item* {
                        return &(layer->at(index));
                    }
            );
        }
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

        ch::ink_layer* layer_;
        poly_tree tree_;
        std::optional<ch::point> curr_loc_;
        ch::dimensions<int> base_sz_;
        double scale_;
        ch::brush_expr_ptr def_brush_;
        rgn_selection selection_;
        std::vector<std::unordered_set<int>> set_brushes_;
        std::vector<ch::color> colors_;
        std::vector<ch::polygon> scaled_regions_;
        bool selecting_;
        int cursor_radius_;
        const std::vector<ch::brush_expr_ptr>* brushes_;
        struct brush_info {
            ch::color color;
            std::unordered_set<int> items;
        };
        std::unordered_map<const ch::brush_expr*, brush_info> brush_info_;
        bool is_showing_brushes_;
        bool is_showing_flow_;
        std::optional<ch::point> flow_pt1_;
        std::optional<ch::point> flow_pt2_;
        std::unordered_map<int, QPixmap> patterns_;

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

        void paintEvent(QPaintEvent* event) override;
        void keyPressEvent(QKeyEvent* event) override;
        void mousePressEvent(QMouseEvent* event) override;
        void mouseMoveEvent(QMouseEvent* event) override;
        void mouseReleaseEvent(QMouseEvent* event) override;

    public:
        rgn_map_ctrl(const std::vector<ch::brush_expr_ptr>* brushes);
        void set(ch::brush_expr_ptr default_brush, double scale, const ch::dimensions<int>& sz, ch::ink_layer* layer);
        void set_scale(double scale);
        const rgn_selection& current_selection() const;
        void set_brush_of_selection(ch::brush_expr_ptr br);
        void set_flow_of_selection(double flow);
        void show_brushes(bool v);
        void show_flow(bool v);
        bool is_showing_brushes() const;
        bool is_showing_flow() const;
        ch::ink_layer* layer() const;

        void handle_flow_drag(bool dragging, const ch::point& pt);

    signals:
        void selection_changed();
        void flow_assigned(double theta);
    };
}