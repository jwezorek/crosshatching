#pragma once

#include "../crosshatching/geometry.hpp"
#include "../crosshatching/brush_lang.hpp"
#include "../crosshatching/ink_layers.hpp"
#include "image_box.h"
#include "settingctrls.hpp"
#include "rgn_map_ctrl.h"
#include <QWidget>
#include <unordered_set>
#include <unordered_map>
#include <memory>

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

    class flow_direction_panel : public QWidget {

        Q_OBJECT

        std::optional<double> direction_;

        float radius() const;
        QPointF center() const;

    public:
        flow_direction_panel();
        void clear();
        void set_direction(double theta);
        std::optional<double> direction() const;
        void paintEvent(QPaintEvent* event) override;
    };

    class main_window;

    class rgn_tool_panel : public QWidget {

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

        std::vector<ch::brush_expr_ptr> get_brush_defaults() const;
        void handle_selection_change();
        rgn_map_ctrl* current_rgn_map() const;
        void handle_brush_change(QString str);
        void handle_flow_assigned(double theta);
        void handle_selection_change_brush();
        void handle_selection_change_flow();
        void handle_brush_name_change(std::string old_name, std::string new_name);
        void repopulate_brush_tables();

    public:
        rgn_tool_panel(main_window* parent, QStackedWidget* stack);
        void repopulate_ctrls();
        void set_layers(double scale, ch::ink_layers* layers);
        void set_scale(double scale);
    };

    class rgn_map_tools {
        QStackedWidget* rgn_map_stack_;
        rgn_tool_panel* rgn_props_;
    public:
        bool has_rgn_maps() const;
        rgn_map_tools();
        void clear();
        void populate(main_window* parent);
        void repopulate_ctrls();
        void set_layers(double scale, ch::ink_layers* layers);
        void set_scale(double scale);
        QStackedWidget* rgn_map_stack() const;
        rgn_tool_panel* rgn_props() const;
    };
}