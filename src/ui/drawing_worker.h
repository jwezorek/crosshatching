#pragma once

#include <functional>
#include <any>
#include <QObject>
#include "../crosshatching/drawing.hpp"

namespace ui {

    class drawing_worker : public QObject
    {
        Q_OBJECT

    public:
        drawing_worker(const ch::crosshatching_job& job);
        void process();
        ch::drawing output() const;

        void update_progress(double val);
        void set_status_line(const std::string& str);
        void log_prog(const std::string& str);

        ~drawing_worker();
    signals:
        void progress(double);
        void status(std::string str);
        void log(std::string str);
        void finished();

    private:
        ch::crosshatching_job job_;
        ch::drawing output_;
    };

    class worker : public QObject
    {
        Q_OBJECT

    private:
        std::function<std::any()> job_;
        std::any output_;
    public:
        worker();
        void set(std::function<std::any()> job);
        std::function<void(double)> progress_function();
        void update_progress(double pcnt);
        void process();
        std::any output() const;

    signals:
        void progress(double);
        void finished();
    };
}