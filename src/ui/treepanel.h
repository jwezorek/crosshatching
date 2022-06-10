#pragma once

#include <QTreeWidget>
#include <QTreeWidgetItem>
#include <QHeaderView>
#include <QPushButton>
#include <functional>

namespace ui {

    using tree_panel_callback = std::function<void(QTreeWidget*, QTreeWidgetItem*)>;

    class tree_panel : public QWidget
    {
        Q_OBJECT

    public:
        tree_panel( tree_panel_callback add_cb, tree_panel_callback delete_cb);

        const QTreeWidget* tree() const;
        QTreeWidget* tree();

        const QPushButton* add_btn() const;
        QPushButton* add_btn();

        const QPushButton* delete_btn() const;
        QPushButton* delete_btn();

    private:
        QTreeWidget* tree_;
        QPushButton* add_btn_;
        QPushButton* delete_btn_;
        tree_panel_callback add_fn_;
        tree_panel_callback delete_fn_;
    };

}