#include "drawing_worker.h"
#include <QThread>
#include <qdebug.h>

ui::drawing_worker::drawing_worker(const ch::crosshatching_job& job) : 
    QObject(nullptr),
    job_(job) 
{}

void ui::drawing_worker::process()
{
    ch::callbacks cbs{
        [this](double pcnt) { this->update_progress(pcnt); },
        [this](const std::string& msg) { this->set_status_line(msg); },
        [this](const std::string& msg) { this->log(msg); }
    };
    
    output_ = ch::generate_crosshatched_drawing(job_, cbs );
    emit finished();
}

ch::drawing ui::drawing_worker::output() const {
    return output_;
}

void ui::drawing_worker::update_progress(double val) {
    emit progress(val);
}
void ui::drawing_worker::set_status_line(const std::string& str) {
    emit status(str);
}

void ui::drawing_worker::log_prog(const std::string& str) {
    emit log(str);
}

ui::drawing_worker::~drawing_worker()
{
}

/*--------------------------------------------------------------------------------------------------*/

ui::worker::worker() :
    QObject(nullptr)
{}

void ui::worker::set(std::function<std::any()> job) {
    job_ = job;
}

std::function<void(double)> ui::worker::progress_function() {
    return [this](double val)->void {
        this->update_progress(val);
    };
}

void ui::worker::update_progress(double pcnt) {
    emit progress(pcnt);
}

void ui::worker::process()
{
    output_ = job_();
    emit finished();
}

std::any ui::worker::output() const {
    return output_;
}
