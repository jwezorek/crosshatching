#include "treepanel.h"
#include <QVBoxLayout>
#include <QHBoxLayout>

namespace {

    constexpr int k_btn_sz = 30;
    constexpr int k_btn_bar_hgt = k_btn_sz + 8;

    QTreeWidgetItem* selected_item(QTreeWidget* tree) {
        auto list = tree->selectedItems();
        if (list.isEmpty()) {
            return nullptr;
        }
        return list.first();
    }
}

ui::tree_panel::tree_panel(tree_panel_callback add_fn, tree_panel_callback delete_fn) :
        QWidget(0),
        add_fn_(add_fn),
        delete_fn_(delete_fn) {

    QVBoxLayout* vlayout = new QVBoxLayout(this);
    QWidget* button_bar = new QWidget();
    button_bar->setFixedHeight(k_btn_bar_hgt);
    QHBoxLayout* hlayout = new QHBoxLayout(button_bar);
    hlayout->setSpacing(0);
    hlayout->addWidget(add_btn_ = new QPushButton("+"));
    add_btn_->setFixedSize(QSize(k_btn_sz, k_btn_sz));
    hlayout->addWidget(delete_btn_ = new QPushButton("-"));
    delete_btn_->setFixedSize(QSize(k_btn_sz, k_btn_sz));
    hlayout->addStretch();

    vlayout->addWidget(tree_ = new QTreeWidget());
    vlayout->addWidget(button_bar);
    vlayout->setSpacing(0);

    connect(add_btn_, &QPushButton::clicked, 
        [=]() {
            add_fn_(tree_,selected_item(tree_));
        }
    );

    connect(delete_btn_, &QPushButton::clicked,
        [=]() {
            delete_fn_(tree_,selected_item(tree_));
        }
    );

}

const QTreeWidget* ui::tree_panel::tree() const {
    return tree_;
}

QTreeWidget* ui::tree_panel::tree() {
    return const_cast<QTreeWidget*>(const_cast<const tree_panel*>(this)->tree());
}

const QPushButton* ui::tree_panel::add_btn() const {
    return add_btn_;
}
QPushButton* ui::tree_panel::add_btn() {
    return const_cast<QPushButton*>(const_cast<const tree_panel*>(this)->add_btn());
}

const QPushButton* ui::tree_panel::delete_btn() const {
    return delete_btn_;
}

QPushButton* ui::tree_panel::delete_btn() {
    return const_cast<QPushButton*>(const_cast<const tree_panel*>(this)->delete_btn());
} 
