#include "brush_panel.hpp"
#include "dialogs.h"

using json = nlohmann::json;

/*------------------------------------------------------------------------------------------*/

ui::brush_panel::brush_panel(layer_panel& layers) :
    tree_panel("brushes", [&]() {this->add_brush_node(); },
        [&]() {this->delete_brush_node(); }),
    layer_panel_(layers) {
    connect(tree(), &QTreeWidget::itemDoubleClicked, this, &brush_panel::handle_double_click);

    auto interlocking_diagonal = std::get<ch::brush_expr_ptr>(
        ch::parse(
            R"(
                    ( 
                    brush
                        (define lines 
                            (quote
                                (brush 
                                    (horz_strokes
                                        (norm_rnd (lerp 100 1600) (lerp 10 100))
                                        (norm_rnd (lerp 200 0) (lerp 20 0.05))
                                        (norm_rnd (lerp 7 0.5) (lerp 0.5 0.05))
                                    )
                                    (dis (ramp 0.20 false true))
                                    (jiggle (norm_rnd 0.0 0.5))
                                )
                            )
                        )
                        (rot 45)
                        lines
                        (rot -45)
                        lines
                    )
            )"
        )
        );
    auto horz = std::get<ch::brush_expr_ptr>(
        ch::parse(
            R"(
                (brush
                    (horz_strokes
                        (norm_rnd (lerp 100 1600) (lerp 10 100))
                        (norm_rnd (lerp 200 0) (lerp 20 0.05))
                        (norm_rnd (lerp 7 0.5) (lerp 0.5 0.05))
                    )
                    (dis (ramp 0.2 0 1))
                    (jiggle (norm_rnd 0 0.5))
                )
            )"
        )
        );
    auto vert = std::get<ch::brush_expr_ptr>(
        ch::parse(
            R"(
                (brush
                    (rot 90)
                    (horz_strokes
                        (norm_rnd (lerp 100 1600) (lerp 10 100))
                        (norm_rnd (lerp 200 0) (lerp 20 0.05))
                        (norm_rnd (lerp 7 0.5) (lerp 0.5 0.05))
                    )
                    (dis (ramp 0.2 0 1))
                    (jiggle (norm_rnd 0 0.5))
                )
            )"
        )
        );
    auto solid = std::get<ch::brush_expr_ptr>(ch::parse("(solid value)"));

    insert_toplevel_item(tree(), "interlocking-diagonal", interlocking_diagonal);
    insert_toplevel_item(tree(), "horizontal", horz);
    insert_toplevel_item(tree(), "vertical", vert);
    insert_toplevel_item(tree(), "solid", solid);
}

std::vector<std::string> ui::brush_panel::brush_names() const {
    std::vector<std::string> brushes(tree()->topLevelItemCount());
    for (int i = 0; i < tree()->topLevelItemCount(); ++i) {
        QTreeWidgetItem* item = tree()->topLevelItem(i);
        brushes[i] = item->text(0).toStdString();
    }
    return brushes;
}


std::vector<ch::brush_expr_ptr> ui::brush_panel::brushes() const {
    std::vector<ch::brush_expr_ptr> brushes(tree()->topLevelItemCount());
    for (int i = 0; i < tree()->topLevelItemCount(); ++i) {
        brush_item* bi = static_cast<brush_item*>(tree()->topLevelItem(i));
        brushes[i] = bi->brush_expression;
    }
    return brushes;
}

std::unordered_map<std::string, ch::brush_expr_ptr> ui::brush_panel::brush_dictionary() const {
    std::unordered_map<std::string, ch::brush_expr_ptr> dictionary;
    for (int i = 0; i < tree()->topLevelItemCount(); ++i) {
        QTreeWidgetItem* item = tree()->topLevelItem(i);
        std::string name = item->text(0).toStdString();
        brush_item* bi = static_cast<brush_item*>(item);
        dictionary[name] = bi->brush_expression;
    }
    return dictionary;
}

std::unordered_map<ch::brush_expr*, std::string> ui::brush_panel::brush_name_dictionary() const {
    std::unordered_map<ch::brush_expr*, std::string> dict;
    for (int i = 0; i < tree()->topLevelItemCount(); ++i) {
        QTreeWidgetItem* item = tree()->topLevelItem(i);
        std::string name = item->text(0).toStdString();
        brush_item* bi = static_cast<brush_item*>(item);
        dict[bi->brush_expression.get()] = name;
    }
    return dict;
}

json ui::brush_panel::to_json() const {
    auto brushes = brush_dictionary();
    json brushes_json = json::array();
    for (const auto& [name, expr] : brushes) {
        json pair = json::array();
        pair.push_back(name);
        pair.push_back(ch::pretty_print(*expr));
        brushes_json.push_back(pair);
    }
    return brushes_json;
}

void ui::brush_panel::from_json(const ch::json& json) {
    tree()->clear();
    for (const auto& pair : json) {
        auto brush_name = pair[0].get<std::string>();
        auto brush_expr = std::get<ch::brush_expr_ptr>(ch::parse(pair[1].get<std::string>()));
        insert_toplevel_item(tree(), brush_name, brush_expr);
    }
}

ui::brush_panel::brush_item::brush_item(const std::string& name, ch::brush_expr_ptr expr) :
    QTreeWidgetItem(static_cast<QTreeWidget*>(nullptr), QStringList(QString(name.c_str()))),
    brush_expression(expr),
    is_toplevel(true)
{}

ui::brush_panel::brush_item::brush_item(ch::brush_expr_ptr expr) :
    QTreeWidgetItem(
        static_cast<QTreeWidget*>(nullptr),
        QStringList(QString(expr->short_string().c_str()))),
    brush_expression(expr),
    is_toplevel(false)
{}

bool ui::brush_panel::brush_item::is_leaf() const {
    return brush_expression->children().empty();
}

void ui::brush_panel::insert_brush_item(brush_item* parent, brush_item* item) {
    parent->addChild(item);
    auto expr = item->brush_expression;
    for (auto child : expr->children()) {
        insert_brush_item(item, new brush_item(child));
    }
}

void ui::brush_panel::insert_toplevel_item(QTreeWidget* tree,
    const std::string& name, ch::brush_expr_ptr expr) {
    auto toplevel_item = new brush_item(name, expr);
    tree->addTopLevelItem(toplevel_item);
    if (expr) {
        for (auto child : expr->children()) {
            insert_brush_item(toplevel_item, new brush_item(child));
        }
    }
}

void ui::brush_panel::add_brush_node() {
    auto result = ui::brush_dialog::create_brush();
    if (result) {
        const auto& [name, brush] = *result;
        insert_toplevel_item(tree(), name, brush);
    }
}

void ui::brush_panel::delete_brush_node() {
    if (!tree()->selectedItems().empty()) {
        auto item = tree()->selectedItems().first();
        if (!item->parent()) {
            tree()->removeItemWidget(item, 0);
            delete item;
        }
    }
}

QTreeWidgetItem* toplevel_parent(QTreeWidgetItem* twi) {
    while (twi->parent()) {
        twi = twi->parent();
    }
    return twi;
}

void ui::brush_panel::handle_double_click(QTreeWidgetItem* item, int column) {
    brush_item* bi = static_cast<brush_item*>(item);
    ch::brush_expr_ptr old_brush;
    ch::brush_expr_ptr new_brush;
    if (bi->is_toplevel) {
        old_brush = bi->brush_expression;
        auto old_name = bi->text(0).toStdString();
        auto result = brush_dialog::edit_brush(
            bi->text(0).toStdString(),
            ch::pretty_print(*bi->brush_expression)
        );
        if (result) {
            auto [new_name, expr] = result.value();
            tree()->removeItemWidget(item, 0);
            delete item;
            insert_toplevel_item(tree(), new_name, expr);
            if (old_name != new_name) {
                emit brush_name_changed(old_name, new_name);
            }
            new_brush = expr;
        }
    } else {
        auto result = brush_dialog::edit_brush_expr(ch::pretty_print(*bi->brush_expression));
        if (result) {
            auto toplevel_item = toplevel_parent(item);
            auto toplevel_brush_item = static_cast<brush_item*>(toplevel_item);
            old_brush = toplevel_brush_item->brush_expression;
            auto parent = static_cast<brush_item*>(item->parent());
            parent->brush_expression->replace_child(bi->brush_expression, result);

            auto toplevel_parent_body = toplevel_brush_item->brush_expression;
            auto toplevel_parent_name = toplevel_brush_item->text(0).toStdString();
            tree()->removeItemWidget(toplevel_item, 0);
            delete toplevel_item;
            insert_toplevel_item(tree(), toplevel_parent_name, toplevel_parent_body);
            new_brush = toplevel_parent_body;
        }
    }
    if (old_brush && (new_brush != old_brush)) {
        emit brush_changed(old_brush, new_brush);
    }
}
