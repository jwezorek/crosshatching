#pragma once

#include <QObject>
#include "../crosshatching/drawing.hpp"

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
