#pragma once

#include <QTreeWidget>
#include <QTreeWidgetItem>
#include <QHeaderView>
#include <QPushButton>
#include <QTableWidget>
#include <functional>

namespace ui {

    using callback_fn = std::function<void()>;

    class tree_panel : public QWidget
    {
        Q_OBJECT

    public:
        tree_panel( const std::string& title, callback_fn add_cb, callback_fn delete_cb );

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
        callback_fn add_fn_;
        callback_fn delete_fn_;
    };

    class list_panel : public QWidget
    {
        Q_OBJECT

    public:
        list_panel(const std::string& title, int columns, callback_fn add_cb, callback_fn delete_cb);

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
        callback_fn add_fn_;
        callback_fn delete_fn_;
    };

}