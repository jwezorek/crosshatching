#include "drawing.hpp"
#include "geometry.hpp"
#include "point_set.hpp"
#include "brush_lang.hpp"
#include "segmentation.hpp"
#include "qdebug.h"
#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/imgproc.hpp>
#include <range/v3/all.hpp>
#include <stdexcept>
#include <unordered_set>
#include <array>
#include <chrono>
#include <memory>
#include <fstream>

#include "brush_lang.hpp"
#include "raster_to_vector.hpp"

/*------------------------------------------------------------------------------------------------*/

namespace r = ranges;
namespace rv = ranges::views;

namespace {

    constexpr int k_typical_number_of_strokes = 10000;

    class progress {
    private:
        std::string job_name_;
        std::function<void(double)> prog_fn_;
        std::function<void(const std::string&)> stat_fn_;
        std::function<void(const std::string&)> log_fn_;
        size_t total_blobs;
        size_t curr_blob;
        size_t layer_blob;
        size_t total_layer_blobs;

        static bool completed_new_quantile(size_t running_count, size_t total_count, int quantile) {
            if (running_count == 0) {
                return true;
            }
            auto curr_quantile = (running_count * quantile) / total_count;
            auto prev_quantile = ((running_count - 1) * quantile) / total_count;

            return curr_quantile > prev_quantile;
        }

        static int pcnt_complete(size_t running_count, size_t total_count) {
            return static_cast<int>((100.0 * running_count) / total_count);
        }

    public:
        progress(const std::string& job_name, const ch::callbacks& cbs) :
            job_name_(job_name),
            prog_fn_(cbs.update_progress_cb),
            stat_fn_(cbs.update_status_cb),
            log_fn_(cbs.log_message_cb),
            total_blobs(0),
            curr_blob(0)
        {};

        void set_polygon_count(size_t n) {
            total_blobs = n;
        }

        void start_new_layer(size_t n) {
            layer_blob = 0;
            total_layer_blobs = n;
        }

        void tick() {
            if (total_blobs == 0) {
                return;
            }
            ++curr_blob;
            ++layer_blob;
            if (prog_fn_) {
                prog_fn_(static_cast<double>(curr_blob) / static_cast<double>(total_blobs));
            }
            if (log_fn_) {
                if (completed_new_quantile(layer_blob, total_layer_blobs, 10)) {
                    log(std::string("        ") +
                        std::to_string(pcnt_complete(layer_blob, total_layer_blobs)) +
                        "% complete...");
                }
            }
        }

        void log(const std::string& msg) {
            if (log_fn_) {
                log_fn_(msg);
            }
        }

        void status(const std::string& msg) {
            if (stat_fn_) {
                stat_fn_(msg);
                log(std::string("----") + msg + std::string("----"));
            }
        }
    };

    std::tuple<size_t, size_t> drawing_job_stats(const ch::ink_layers& layers) {
        size_t poly_count = r::accumulate(
            layers.content | rv::transform([](auto&& l) {return l.content.size(); }),
            size_t{ 0 }
        );
        size_t vert_count = r::accumulate(
            layers.content | rv::transform(
                [](auto&& l) {
                    return r::accumulate(
                        l.content | rv::transform(
                            [](auto&& ili) {return ch::vert_count(ili.poly); }
                        ),
                        size_t{ 0 }
                    );
                }
            ),
            size_t{ 0 }
        );
        return { poly_count, vert_count };
    }

    using swatch_table = std::unordered_map<ch::brush_token, ch::bkgd_swatches>;
    swatch_table create_swatch_table(int n, int dim) {
        swatch_table tbl;
        tbl[0] = std::vector<cv::Mat>(n, ch::blank_monochrome_bitmap(dim));
        return tbl;
    }

    std::tuple<std::vector<ch::stroke_group>, swatch_table> draw_ink_layer(
            const ch::ink_layer& layer, swatch_table& tbl,
            const ch::parameters& params, progress& prog) {

        size_t n = layer.content.size();
        prog.start_new_layer(n);
        std::vector<ch::stroke_group> output;
        std::unordered_map<ch::brush_token, ch::brush_ptr> brush_table;
        swatch_table output_table = create_swatch_table(
            ch::brush::num_samples(),
            ch::brush::swatch_dim()
        );

        for (const auto& blob: layer.content) {
            auto parent_tok = blob.parent_token();
            auto iter = brush_table.find(parent_tok);
            ch::brush_ptr current_brush;
            if (iter != brush_table.end()) {
                current_brush = iter->second;
            } else {
                auto bkgds = tbl.at(parent_tok);
                current_brush = std::make_shared<ch::brush>(
                    layer.brush,
                    params.epsilon,
                    ch::dimensions(static_cast<double>(params.swatch_sz)),
                    bkgds
                );
                current_brush->build_n(20);
                brush_table[parent_tok] = current_brush;
            }

            auto value = static_cast<double>(blob.value) / 255.0;
            auto strokes = current_brush->draw_strokes(blob.poly, value);
            std::copy(strokes->begin(), strokes->end(), std::back_inserter(output));

            auto tok = blob.token();
            if (output_table.find(tok) == output_table.end()) {
                output_table[tok] = current_brush->render_swatches(value);
            }
            prog.tick();
        }
        output.shrink_to_fit();

        return { std::move(output), std::move(output_table) };
    }

    ch::drawing draw(const ch::ink_layers& inp_layers, const ch::parameters& params, progress& prog) {

        auto [poly_count, vert_count] = drawing_job_stats(inp_layers);
        prog.log("  (job contains " + std::to_string(poly_count) + " polygons with " +
            std::to_string(vert_count) + " total vertices");
        prog.set_polygon_count(poly_count);

        auto layers = scale(inp_layers, params.scale);
        
        std::vector<std::vector<ch::stroke_group>> layer_strokes(layers.content.size());
        swatch_table tok_to_bkgd = create_swatch_table(
            ch::brush::num_samples(),
            ch::brush::swatch_dim()
        );
        for (auto [index, layer] : rv::enumerate(layers.content | rv::reverse)) {
            prog.log(std::string("  - layer ") + std::to_string(index));
            std::tie(layer_strokes[index], tok_to_bkgd) = draw_ink_layer(
                layer, tok_to_bkgd, params, prog);
        }

        return ch::drawing{
            layer_strokes | rv::join | r::to_vector,
            layers.sz
        };
    }
}

ch::drawing ch::generate_crosshatched_drawing(const ch::crosshatching_job& job, const callbacks& cbs) {
    progress prog(job.title, cbs);

    auto start_time = std::chrono::high_resolution_clock().now();

    prog.status(std::string("drawing ") + job.title);
    auto result = draw(job.layers, job.params, prog);
    prog.status(std::string("complete.")); // TODO: or error

    std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock().now() - start_time;
    prog.log(std::string("( ") + std::to_string(elapsed.count()) + " seconds)");

    return result;
}


size_t ch::drawing::stroke_count() const {
    return r::accumulate(
        content |
            rv::transform(
                [](const auto& group)->size_t {
                    return group.strokes.size();
                }
            )
        ,
       size_t{ 0 }
    );
}

void ch::write_to_svg(const std::string& filename, const drawing& d,
        std::function<void(double)> update_progress) {
    std::ofstream outfile(filename);

    outfile << svg_header(static_cast<int>(d.sz.wd), static_cast<int>(d.sz.hgt));

    auto n = static_cast<double>(d.stroke_count());
    for (const auto& [index, poly_info] :rv::enumerate(d.strokes())) {
        const auto& [poly, thickness] = poly_info;
        outfile << polyline_to_svg(poly, thickness) << std::endl;
        if (update_progress) {
            update_progress(index / n);
        }
    }
    outfile << "</svg>" << std::endl;
    outfile.close();
}

cv::Mat ch::paint_drawing(const ch::drawing& d, std::function<void(double)> update_progress_cb) {
    cv::Mat mat(static_cast<int>(d.sz.hgt), static_cast<int>(d.sz.wd), CV_8U, 255);
    return ch::strokes_to_mat(d.content, mat);
}

void ch::debug_drawing() {

}
