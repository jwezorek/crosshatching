#pragma once

#include <QTreeWidget>
#include <QTreeWidgetItem>
#include <QHeaderView>
#include <QPushButton>
#include <QTableWidget>
#include <functional>

namespace ui {

    using tree_panel_callback = std::function<void(QTreeWidget*, QTreeWidgetItem*)>;

    class tree_panel : public QWidget
    {
        Q_OBJECT

    public:
        tree_panel( const std::string& title, tree_panel_callback add_cb, tree_panel_callback delete_cb );

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

    using list_panel_callback = std::function<void(QTableWidget*)>;
    class list_panel : public QWidget
    {
        Q_OBJECT

    public:
        list_panel(const std::string& title, int columns, list_panel_callback add_cb, list_panel_callback delete_cb);

        const QTableWidget* list() const;
        QTableWidget* list();

        const QPushButton* add_btn() const;
        QPushButton* add_btn();

        const QPushButton* delete_btn() const;
        QPushButton* delete_btn();

    private:
        QTableWidget* list_;
        QPushButton* add_btn_;
        QPushButton* delete_btn_;
        list_panel_callback add_fn_;
        list_panel_callback delete_fn_;
    };

}