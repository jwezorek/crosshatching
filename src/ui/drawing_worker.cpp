#include "drawing_worker.h"
#include <QThread>
#include <qdebug.h>

drawing_worker::drawing_worker(const ch::crosshatching_job& job) : 
    QObject(nullptr),
    job_(job) 
{}

void drawing_worker::process()
{
    ch::prog_callbacks cbs{
        [this](double pcnt) { this->update_progress(pcnt); },
        [this](const std::string& msg) { this->set_status_line(msg); },
        [this](const std::string& msg) { this->log(msg); }
    };
    
    output_ = ch::generate_crosshatched_drawing(job_, cbs );
    emit finished();
}

ch::drawing drawing_worker::output() const {
    return output_;
}

void drawing_worker::update_progress(double val) {
    emit progress(val);
}
void drawing_worker::set_status_line(const std::string& str) {
    emit status(str);
}

void drawing_worker::log_prog(const std::string& str) {
    emit log(str);
}

drawing_worker::~drawing_worker()
{
}
