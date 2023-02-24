#include "layer_panel.hpp"
#include "dialogs.h"

namespace r = ranges;
namespace rv = ranges::views;
using json = nlohmann::json;

/*------------------------------------------------------------------------------------------------------*/

ui::layer_panel::layer_panel() :
    ui::list_panel("layers", 3, [&]() { this->add_layer(); }, [&]() { this->delete_layer(); }) {

    //layers_->horizontalHeader()->setStretchLastSection(true);
    list()->horizontalHeader()->setSectionResizeMode(0, QHeaderView::Stretch);
    list()->horizontalHeader()->setSectionResizeMode(1, QHeaderView::Fixed);
    list()->horizontalHeader()->setSectionResizeMode(2, QHeaderView::Fixed);
    list()->setColumnWidth(1, 75);
    list()->setColumnWidth(2, 75);
    //list()->setSelectionMode(QAbstractItemView::NoSelection);
    list()->setSelectionBehavior(QAbstractItemView::SelectRows);
    list()->connect(list(), &QTableWidget::cellDoubleClicked, this, &layer_panel::cell_double_clicked);
    list()->setHorizontalHeaderLabels(QStringList{ "brush", "from value", "to value" });

    setEnabled(false);
}

void ui::layer_panel::set_brush_names(const std::vector<std::string>& brush_names) {
    brush_names_ = brush_names;
    setEnabled(!brush_names_.empty());
}

std::vector<std::tuple<std::string, double>> ui::layer_panel::layers() const {
    return layers_ |
        rv::transform(
            [](const auto& p)->std::tuple<std::string, double> {
                const auto& [val, name] = p;
                return { name, val };
            }
    ) | r::to_vector;
}

bool ui::layer_panel::has_content() const {
    return !layers_.empty();
}

std::tuple<std::string, double> ui::layer_panel::row(int n) const {
    auto brush = list()->item(n, 0)->text().toStdString();
    auto val = list()->item(n, 2)->text().toDouble();
    return { brush,val };
}

void ui::layer_panel::add_layer() {
    auto result = ui::layer_dialog::create_layer_item(brush_names_, list()->rowCount() == 0);
    if (!result) {
        return;
    }
    auto [brush, end_of_range] = *result;
    insert_layer(brush, end_of_range);
    emit layers_changed();
}

void ui::layer_panel::cell_double_clicked(int r, int col) {
    auto [brush, val] = row(r);
    auto result = ui::layer_dialog::edit_layer_item(brush_names_, brush, val);
    if (result) {
        layers_.erase(val);
        auto [brush, end_of_range] = *result;
        insert_layer(brush, end_of_range);
        emit layers_changed();
    }
}

void ui::layer_panel::delete_layer() {
    auto current_row = list()->currentRow();
    if (current_row >= 0) {
        auto [brush, val] = row(current_row);
        layers_.erase(val);
        sync_layers_to_ui();
        emit layers_changed();
    }
}

void ui::layer_panel::insert_layer(const std::string& brush, double end_of_range) {
    layers_[end_of_range] = brush;
    sync_layers_to_ui();
    emit layers_changed();
}

void ui::layer_panel::setRowText(int row, const std::string& brush,
    const std::string& from, const std::string& to) {
    std::array<QTableWidgetItem*, 3> items = {
        new QTableWidgetItem(brush.c_str()),
        new QTableWidgetItem(from.c_str()),
        new QTableWidgetItem(to.c_str())
    };
    for (auto [col, item] : rv::enumerate(items)) {
        item->setFlags(item->flags() & ~Qt::ItemIsEditable);
        list()->setItem(row, col, item);
    }
}

void ui::layer_panel::sync_layers_to_ui() {
    list()->setRowCount(layers_.size());
    int row = 0;
    std::string prev = "0.0";
    for (const auto& [val, name] : layers_) {
        std::string curr = std::to_string(val);
        setRowText(row, name, prev, curr);
        prev = curr;
        ++row;
    }
}

json ui::layer_panel::to_json() const {
    json js_ary = json::array();
    for (const auto& [val, brush_name] : layers_) {
        json pair = json::array();
        pair.push_back(val);
        pair.push_back(brush_name);
        js_ary.push_back(pair);
    }
    return js_ary;
}

void ui::layer_panel::from_json(const json& json) {
    layers_.clear();
    for (const auto& pair : json) {
        auto thresh = pair[0].get<double>();
        auto brush_name = pair[1].get<std::string>();
        layers_[thresh] = brush_name;
        sync_layers_to_ui();
    }
}