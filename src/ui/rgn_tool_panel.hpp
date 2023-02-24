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

/*------------------------------------------------------------------------------------------------*/

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
    class rgn_map_tools;

    class rgn_tool_panel : public QWidget {

        Q_OBJECT

    private:
        QComboBox* layer_cbo_;
        QCheckBox* brush_cbx_;
        QCheckBox* flow_cbx_;
        QLabel* brush_lbl_;
        select_button* select_brush_btn_;
        flow_direction_panel* flow_ctrl_;

    public:
        rgn_tool_panel();
        void set_brush_names( const std::vector<std::string>& brushes);
        void set_layers(int n);

        void set_brush_name(const std::string& brush_name = {});
        void set_hetero_brush();
        void set_flow(std::optional<float> flow);
        void connect_to_tools(rgn_map_tools* tools);
    };

    class rgn_map_tools {
        main_window* parent_;
        QStackedWidget* rgn_map_stack_;
        rgn_tool_panel* rgn_props_;

        std::vector<ch::brush_expr_ptr> get_brush_defaults() const;
        std::unordered_map<ch::brush_expr*, std::string> brush_to_name_dict() const;

    public:
        rgn_map_tools(main_window* parent);
        bool has_rgn_maps() const;
        void clear();
        void populate();
        void set_layers(double scale, ch::ink_layers* layers);
        QStackedWidget* rgn_map_stack() const;
        rgn_tool_panel* rgn_props() const;
        void set_scale(double scale);
        ui::rgn_map_ctrl* current_rgn_map() const;
        void handle_brush_assigned(QString str);
        void handle_flow_assigned(double theta);
        void handle_selection_change_flow();
        void handle_selection_change_brush(); 
        void handle_selection_change();
        void handle_brush_name_change();
        void handle_brush_change(ch::brush_expr_ptr old_brush, ch::brush_expr_ptr new_brush);
    };
}