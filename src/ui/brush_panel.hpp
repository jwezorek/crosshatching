#pragma once

#include "layer_panel.hpp"
#include "../crosshatching/brush_lang.hpp"

namespace ui {

    class brush_panel : public ui::tree_panel {

        Q_OBJECT

    public:

        brush_panel(layer_panel& layers);
        std::vector<std::string> brush_names() const;
        std::vector<ch::brush_expr_ptr> brushes() const;
        std::unordered_map<std::string, ch::brush_expr_ptr> brush_dictionary() const;
        std::unordered_map<ch::brush_expr*, std::string> brush_name_dictionary() const;
        ch::json to_json() const;
        void from_json(const ch::json& json);
        void sync_layer_panel();

    private:

        layer_panel& layer_panel_;

        class brush_item : public QTreeWidgetItem {
        public:
            brush_item(const std::string& name, ch::brush_expr_ptr expr);
            brush_item(ch::brush_expr_ptr expr);

            ch::brush_expr_ptr brush_expression;
            bool is_toplevel;

            bool is_leaf() const;
        };

        static void insert_brush_item(brush_item* parent, brush_item* item);
        static void insert_toplevel_item(QTreeWidget* tree, const std::string& name, ch::brush_expr_ptr expr);
        void add_brush_node();
        void delete_brush_node();
        void handle_double_click(QTreeWidgetItem* item, int column);
    };

}