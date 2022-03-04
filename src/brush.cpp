#include "brush.hpp"

namespace {

    class pipeline_visitor {
    private:
        double t_;
        ch::hatching_range hatching_;
    public:
        pipeline_visitor(double t) : t_(t) {
        }

        void operator()(const ch::param_adapter_fn& param_adapter) {
            t_ = param_adapter(t_);
        }

        void operator()(const ch::brush_adapter_fn& brush_adapter) {
            hatching_ = brush_adapter(hatching_, t_);
        }

        void operator()(const ch::brush_fn& brush) {
            hatching_ = brush(t_);
        }

        ch::hatching_range output() const {
            return hatching_;
        }
    };

}

ch::hatching_range ch::run_brush_pipeline(const brush_pipeline& pipeline, double t)
{
    pipeline_visitor visitor(t);
    for (const auto& item : pipeline) {
        std::visit(visitor, item);
    }
    return visitor.output();
}
