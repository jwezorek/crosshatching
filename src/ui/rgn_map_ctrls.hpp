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

    class select_button : public QPushButton
    {
        Q_OBJECT

    public:
        select_button(const std::string& txt);
        void set_items(std::span<const std::string> items);
        void show_popup();
        void item_selected(QListWidgetItem* item);

    signals:
        void item_clicked(QString str);

    private:
        QListWidget* list_;
        QWidget* popup_;
    };

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

        ch::ink_layer* layer_;
        poly_tree tree_;
        std::optional<ch::point> curr_loc_;
        ch::dimensions<int> base_sz_;
        double scale_;
        ch::brush_expr_ptr def_brush_;
        std::unordered_set<int> selected_;
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
        rgn_map_ctrl(const std::vector<ch::brush_expr_ptr>* brushes);
        void set(ch::brush_expr_ptr default_brush, double scale, const ch::dimensions<int>& sz, ch::ink_layer* layer);
        const std::unordered_set<int>& selected() const;
        bool has_selection() const;
        void set_brush_of_selection(ch::brush_expr_ptr br);
        
        auto selected_layer_items() const {
            namespace r = ranges;
            namespace rv = ranges::views;
            return selected_ |
                rv::transform(
                    [this](auto index)->ch::ink_layer_item* {
                        return &(layer_->at(index));
                    }
                );
        }

        void show_brushes(bool v);
        void show_flow(bool v);
        bool is_showing_brushes() const;
        bool is_showing_flow() const;

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
        QLabel* brush_lbl_;
        select_button* select_brush_btn_;
        flow_direction_panel* flow_ctrl_;
        QStackedWidget* stack_;

        std::unordered_map<std::string, ch::brush_expr_ptr> name_to_brush_;
        std::unordered_map<ch::brush_expr*, std::string> brush_to_name_;
        std::vector<ch::brush_expr_ptr> default_brushes_;
        std::vector<ch::brush_expr_ptr> brushes_;

        std::vector<ch::brush_expr_ptr> get_brush_defaults() const;
        void handle_selection_change();
        rgn_map_ctrl* current_rgn_map() const;
        void handle_brush_change(QString str);

    public:
        rgn_map_panel(main_window* parent, QStackedWidget* stack);
        void repopulate_ctrls();
        void set_layers(double scale, ch::ink_layers* layers);
    };

    struct rgn_map_tools {
        QStackedWidget* rgn_map_stack;
        rgn_map_panel* rgn_props;
    };
}