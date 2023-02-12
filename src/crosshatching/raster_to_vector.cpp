#include "raster_to_vector.hpp"
#include "point_set.hpp"
#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/imgproc.hpp>
#include <range/v3/all.hpp>

#include "correct.hpp"

namespace r = ranges;
namespace rv = ranges::views;

/*------------------------------------------------------------------------------------------------------*/

namespace {

	using int_polyline = std::vector<cv::Point>;

	struct find_contour_output {
		std::vector<int_polyline> contours;
		std::vector<cv::Vec4i> hierarchy;
	};

	template<typename T>
	std::vector<std::tuple<T, find_contour_output>> perform_find_contours(
		const cv::Mat& input) {
		auto planes = color_masks<T>(input);
		return planes |
			rv::transform(
				[](const auto& pair)->std::tuple<T, find_contour_output> {
					const auto& [value, mat] = pair;
					find_contour_output output;
					cv::findContours(
						mat, output.contours, output.hierarchy, 
						cv::RETR_CCOMP, cv::CHAIN_APPROX_NONE
					);

					return { value, std::move(output) };
				}
		) | r::to_vector;
	}

	enum class direction {
		N, NE, E, SE, S, SW, W, NW
	};

	struct pixel_crawler {
		cv::Point loc;
		direction head;
	};

	pixel_crawler roll(const pixel_crawler& pc) {
		static std::unordered_map<direction, direction> dir_to_next_dir = {
			{direction::NW, direction::SW},
			{direction::SW, direction::SE},
			{direction::SE, direction::NE},
			{direction::NE, direction::NW}
		};
		return { pc.loc, dir_to_next_dir[pc.head] };
	}

	cv::Point normalize_offset(const cv::Point& pt) {
		auto offset = pt;
		offset /= std::max(std::abs(pt.x), std::abs(pt.y));
		return offset;
	}

	direction direction_to(const cv::Point& from_pt, const cv::Point& to_pt) {
		static ch::int_point_map<direction> offset_to_direction = {
			{{0,-1}, direction::N },
			{{1,-1}, direction::NE},
			{{1,0},  direction::E },
			{{1,1},  direction::SE},
			{{0,1},  direction::S },
			{{-1,1}, direction::SW},
			{{-1,0}, direction::W },
			{{-1,-1},direction::NW}
		};

		return offset_to_direction[normalize_offset(to_pt - from_pt)];
	}

	cv::Point2f get_vertex(const pixel_crawler& pc) {
		static std::unordered_map<direction, cv::Point2f> dir_to_vert_offset = {
			{direction::NW, {0,0} },
			{direction::NE, {1,0}},
			{direction::SE, {1,1}},
			{direction::SW, {0,1}}
		};
		return cv::Point2f(pc.loc) + dir_to_vert_offset[pc.head];
	}

	bool is_cardinal_dir(direction dir) {
		static std::unordered_set<direction> nesw = {
			direction::N,
			direction::E,
			direction::S,
			direction::W
		};
		return nesw.find(dir) != nesw.end();
	}

	bool is_ordinal_dir(direction dir) {
		return !is_cardinal_dir(dir);
	}

	pixel_crawler flip(const pixel_crawler& pc, const cv::Point& to_pt) {
		auto dir = direction_to(pc.loc, to_pt);
		static std::unordered_map<direction, direction> dir_to_flip_dir = {
			{direction::N , direction::S },
			{direction::NE, direction::SW},
			{direction::E , direction::W },
			{direction::SE, direction::NW},
			{direction::S , direction::N },
			{direction::SW, direction::NE},
			{direction::W , direction::E },
			{direction::NW, direction::SE}
		};
		auto flipped_crawler = pixel_crawler{ to_pt, dir_to_flip_dir[pc.head] };
		if (is_ordinal_dir(dir)) {
			return flipped_crawler;
		} else {
			return  roll(flipped_crawler);
		}
	}

	using dir_pair = std::tuple<direction, direction>;
	std::tuple<direction, direction> get_shared_vert_directions(
			const cv::Point& from_pt, const cv::Point& to_pt) {
		static std::unordered_map<direction, dir_pair> dir_to_shared_verts = {
			{direction::N , {direction::NW, direction::NE}},
			{direction::NW, {direction::NW, direction::NW}},
			{direction::W , {direction::NW, direction::SW} },
			{direction::SW, {direction::SW, direction::SW}},
			{direction::S , {direction::SW, direction::SE} },
			{direction::SE, {direction::SE, direction::SE}},
			{direction::E , {direction::NE, direction::SE} },
			{direction::NE, {direction::NE, direction::NE}}
		};
		return dir_to_shared_verts[direction_to(from_pt, to_pt)];
	}

	bool can_flip(const pixel_crawler& pc, const cv::Point& to_pt) {
		if (pc.loc == to_pt) {
			return false;
		}
		auto [shared_vert_1, shared_vert_2] = get_shared_vert_directions(pc.loc, to_pt);
		return pc.head == shared_vert_1 || pc.head == shared_vert_2;
	}

	int find_northwest_index(const int_polyline& ip) {
		if (ip.size() == 1) {
			return 0;
		}
		int north_west_index = 0;
		for (int i = 1; i < static_cast<int>(ip.size()); ++i) {
			const auto& current_min = ip[north_west_index];
			const auto& pt = ip[i];
			if (pt.y < current_min.y) {
				north_west_index = i;
			} else if (pt.y == current_min.y && pt.x < current_min.x) {
				north_west_index = i;
			}
		}
		return north_west_index;
	}

	std::tuple<pixel_crawler, int> initialize_contour_crawl(
			const int_polyline& ip, bool counter_clockwise) {
		int northwest_index = find_northwest_index(ip);
		return { 
			pixel_crawler{
				ip[northwest_index], 
				counter_clockwise ? direction::NW : direction::SW
			},
			northwest_index 
		};
	}

	void push_if_new(ch::ring& poly, const ch::point& pt) {
		if (!poly.empty() && poly.back() == pt) {
			return;
		}
		poly.push_back(pt);
	}

	template<typename R>
	auto rotate_view(R rng, int pivot) {
		return ranges::views::concat(
			rng | ranges::views::drop(pivot),
			rng | ranges::views::take(pivot)
		);
	}

	auto canonicalized_cyclic_contour_view(const int_polyline& ip) {
		int start_index = find_northwest_index(ip);
		auto start_point = ip.at(start_index);
		return rotate_view(rv::all(ip), start_index) | rv::cycle;
	}

	template<typename T>
	auto neighbor_view(T rng) {
		return rng | rv::sliding(2) |
			rv::transform(
				[](auto pair)->std::tuple<cv::Point, cv::Point> {
					return { pair[0],pair[1] };
				}
		);
	}

	ch::ring contour_to_ring(const int_polyline& ip, bool counter_clockwise = true) {
		ch::ring poly;
		int n = static_cast<int>(ip.size());
		poly.reserve(n);

		auto contour = canonicalized_cyclic_contour_view(ip);
		auto crawler = pixel_crawler{
			contour[0], counter_clockwise ? direction::NW : direction::SW
		};
		auto neighbors = neighbor_view(contour);
		auto iter = neighbors.begin();

		do {
			const auto& [current, next] = *iter;
			auto current_vert = get_vertex(crawler);
			push_if_new(poly, current_vert);
			if (can_flip(crawler, next)) {
				crawler = flip(crawler, next);
				++iter;
			}
			crawler = roll(crawler);
		} while (get_vertex(crawler) != poly.front());

		poly.shrink_to_fit();
		return poly;
	}

	std::vector<ch::polygon> contour_info_to_polygons(
		const find_contour_output& contours) {

		std::unordered_map<int, int> contour_index_to_poly_index;
		std::vector<ch::polygon> blobs;
		for (int i = 0; i < contours.hierarchy.size(); ++i) {
			const cv::Vec4i& h = contours.hierarchy[i];
			int parent = h[3];
			const auto& contour = contours.contours[i];
			if (parent == -1) { // contour is an outer polygon
				int poly_index = static_cast<int>(blobs.size());
				blobs.emplace_back(
					ch::make_polygon(contour_to_ring(contour, true), {})
				);
				contour_index_to_poly_index[i] = poly_index;
			} else { // contour is a hole, so find its parent...
				auto iter = contour_index_to_poly_index.find(parent);
				if (iter == contour_index_to_poly_index.end()) {
					throw std::runtime_error("contour hierearchy had a child contour before its parent");
				}
				blobs[iter->second].inners().push_back(contour_to_ring(contour, false));
			}
		}
		return blobs;
	}

	template<typename T>
	std::vector<std::tuple<T, cv::Mat>> color_masks(const cv::Mat& input) {
		auto levels = ch::unique_values<T>(input);
		return
			levels |
			rv::transform(
				[&input](auto val)->std::tuple<T, cv::Mat> {
					cv::Mat output;
					cv::inRange(input, val, val, output);
					return { val, output };
				}
		) | r::to_vector;
	}

	using int_polyline = std::vector<cv::Point>;


	template<typename T>
	using blobs_per_gray_t = std::tuple<T, std::vector<ch::polygon>>;

	template<typename T>
	std::vector<std::tuple<T, ch::polygon>> blobs_per_gray_to_layer(
		const std::vector< blobs_per_gray_t<T>>& blobs_per_gray) {
		using blob = std::tuple< T, ch::polygon>;
		return rv::join(
			blobs_per_gray |
			rv::transform(
				[](const auto& bpg) {
					const auto& [value, polygons] = bpg;
					return polygons |
						rv::transform(
							[value](const auto& poly)->blob {
								return {
									value,
									poly
								};
							}
					);
				}
			)
		) | r::to_vector;
	}

	template<typename T>
	std::vector<std::tuple<T, ch::polygon>> image_to_polygons(const cv::Mat& img) {

		auto contours = perform_find_contours<T>(img);
		std::vector<blobs_per_gray_t<T>> blobs_per_gray = contours |
			rv::transform(
				[](const auto& tup)->blobs_per_gray_t<T> {
					const auto& [value, contour_info] = tup;
					auto polygons = contour_info_to_polygons(contour_info);
					return { value, std::move(polygons) };
				}
			) |
			r::to_vector;

		return blobs_per_gray_to_layer(blobs_per_gray);
	}

	int update_vertex_usage_count(ch::point_map<int>& vertex_count, const ch::point& pt,
			const ch::dimensions<int>& dims) {
		bool not_in_table = vertex_count.find(pt) == vertex_count.end();
		if (not_in_table && (pt.x == 0 || pt.x == dims.wd)) {
			++vertex_count[pt];
		}
		if (not_in_table && (pt.y == 0 || pt.y == dims.hgt)) {
			++vertex_count[pt];
		}
		return ++vertex_count[pt];
	}

	r::any_view<ch::point> update_vertex_usage_counts(const ch::ring& ring,
		ch::point_map<int>& vertex_count, const ch::dimensions<int>& dims) {
		return ring | rv::remove_if(
			[&](const ch::point& pt)->bool {
				int count = update_vertex_usage_count(vertex_count, pt, dims);
				return count != 3;
			}
		);
	}

	r::any_view<ch::point> update_vertex_usage_counts(const ch::polygon& poly,
		ch::point_map<int>& vertex_count, const ch::dimensions<int>& dims) {
		return rv::concat(
			rv::single(update_vertex_usage_counts(poly.outer(), vertex_count, dims)),
			poly.inners() |
			rv::transform(
				[&](const ch::ring& r) {
					return update_vertex_usage_counts(r, vertex_count, dims);
				}
			)
		) | rv::join;
	}

	bool is_critical_point_free(const ch::point_set& set, std::span<const ch::point> pts) {
		for (ch::point pt : pts) {
			if (set.find(pt) != set.end()) {
				return false;
			}
		}
		return true;
	}

	void insert_island_critical_points(ch::point_set& set, const ch::ring& r) {
		if (is_critical_point_free(set, r)) {
			set.insert(r.front());
			set.insert(ch::southeast_most_point(r));
		}
	}

	void insert_island_critical_points(
		ch::point_set& set, std::span<const ch::polygon> polygons) {
		for (const auto& poly : polygons) {
			insert_island_critical_points(set, poly.outer());
			for (const auto& hole : poly.inners()) {
				insert_island_critical_points(set, hole);
			}
		}
	}

	ch::dimensions<int> get_dimensions(std::span<const ch::polygon> polygons) {
		auto [x1, y1, x2, y2] = ch::bounding_rectangle(polygons);
		return { x2 - x1 - 1, y2 - y1 - 1 };
	}

	ch::point_set generate_critical_points_set(
		std::span<const ch::polygon> polygons)
	{
		ch::dimensions<int> dims = get_dimensions(polygons);
		ch::point_map<int> vertex_usage_counts;
		vertex_usage_counts.reserve(polygons.size() * 10);

		auto set = polygons |
			rv::transform(
				[&vertex_usage_counts, &dims](const auto& poly) {
					return update_vertex_usage_counts(poly, vertex_usage_counts, dims);
				}
			) |
			rv::join |
			r::to_<ch::point_set>();

		insert_island_critical_points(set, polygons);
		return set;
	}

	using edge = std::tuple<ch::ring::const_iterator, ch::ring::const_iterator>;

	std::vector<edge> get_shared_edges(
		const ch::ring& verts, const ch::point_set& critical_points) {

		auto iter_list = rv::iota(0, static_cast<int>(verts.size())) |
			rv::transform(
				[&verts](int index)->ch::ring::const_iterator {
					return verts.cbegin() + index;
				}
		) | r::to_vector;

		auto critical_pts_ring = iter_list | rv::remove_if(
			[&critical_points](const auto& i) {
				return critical_points.find(*i) == critical_points.end();
			}
		) | r::to_vector;

		return critical_pts_ring |
			rv::cycle |
			rv::sliding(2) |
			rv::take(critical_pts_ring.size()) |
			rv::transform(
				[](auto pair)->edge {
					return { pair[0] , pair[1] };
				}
		) | r::to_vector;
	}

	std::vector<ch::point> expand_edge(const edge& e, const ch::ring& ring) {
		auto [u, v] = e;
		std::vector<ch::point> path;
		if (u < v) {
			std::copy(u, v, std::back_inserter(path));
		} else {
			std::copy(u, ring.end(), std::back_inserter(path));
			std::copy(ring.begin(), v, std::back_inserter(path));
		}
		path.push_back(*v);
		return path;
	}

	std::span<const ch::point> lookup_path_or_expand_and_simplify(
			const edge& e, const ch::ring& ring, const ch::point_set& crit_pts,
			ch::path_table& paths, double param) {

		auto edge = expand_edge(e, ring);
		auto pid = ch::make_path_id(edge);
		if (paths.contains(pid)) {
			return paths.find(pid);
		}

		auto simplified_edge = ch::perform_douglas_peucker_simplification(edge, param);
		auto path_in_table = paths.insert(edge, simplified_edge);
		paths.insert(
			edge | rv::reverse | r::to_vector,
			simplified_edge | rv::reverse | r::to_vector
		);
		return path_in_table;
	}

	ch::ring simplify_ring(const ch::ring& ring, const ch::point_set& crit_pts,
			ch::path_table& paths, double param) {
		auto edges = get_shared_edges(ring, crit_pts);
		return edges |
			rv::transform(
				[&](const auto& e) {
					auto simplified_edge = lookup_path_or_expand_and_simplify(
						e, ring, crit_pts, paths, param
					);
					auto n = simplified_edge.size();
					return simplified_edge | rv::take(n - 1);
				}
			) |
			rv::join |
			r::to_<ch::ring>;
	}

	template<typename T>
	bool is_degenerate_poly_tuple(const std::tuple<T, ch::polygon>& tup) {
		return is_degenerate_ring(std::get<1>(tup).outer());
	}

	ch::polygon simplify_polygon(const ch::polygon& poly, const ch::point_set& crit_pts,
			ch::path_table& paths, double param) {
		return ch::make_polygon(
			simplify_ring(poly.outer(), crit_pts, paths, param),
			poly.inners() |
				rv::transform(
					[&](const auto& hole) {
						return simplify_ring(hole, crit_pts, paths, param);
					}
				) | 
                rv::remove_if( ch::is_degenerate_ring ) | 
                r::to_vector
		);
	}

	std::vector<ch::polygon> simplify_polygons(
		std::span<const ch::polygon> dissection, double param) {

		const auto critical_points = generate_critical_points_set(dissection);
		ch::path_table paths;
		return dissection |
			rv::transform(
				[&](const ch::polygon& poly) {
					return simplify_polygon(poly, critical_points, paths, param);
				}
		) | r::to_vector;
	}

	template<typename... Args>
	auto first(const std::tuple<Args...>& tup)->decltype(auto) {
		return std::get<0>(tup);
	}

	template<typename... Args>
	auto second(const std::tuple<Args...>& tup)->decltype(auto) {
		return std::get<1>(tup);
	}

    template<typename T>
    std::vector<std::tuple<T, ch::polygon>> correct_polygons(
            const std::vector<std::tuple<T, ch::polygon>>& inp) {
        
        std::vector<std::tuple<T, ch::polygon>> output;
        output.reserve(inp.size() + 100); 
        for (auto&& [col, poly] : inp) {
            if (ch::is_degenerate_ring(poly.outer())) {
                continue;
            }

            if (boost::geometry::is_valid(poly)) {
                output.emplace_back(col, std::move(poly));
            }

            // fix self-intersections, spikes, etc.
            double remove_spike_threshold = 1E-7;
            ch::polygons result;
            geometry::correct<ch::point, ch::polygon, ch::ring, ch::polygons>(
                ch::polygons{ poly }, result, remove_spike_threshold
                );

            for (auto&& p : result) {
                if (ch::is_degenerate_ring(p.outer())) {
                    continue;
                }
                output.emplace_back(col, std::move(p));
            }
        }
        output.shrink_to_fit();
        return output;
    }

	template<typename T>
	std::vector<std::tuple<T, ch::polygon>> simplify_colored_polygons(
		std::span< std::tuple<T, ch::polygon>> blobs, double param) {
		namespace r = ranges;
		namespace rv = ranges::views;

		auto polys = blobs |
			rv::transform(second<T, ch::polygon>) |
			r::to_vector;
		
		polys = simplify_polygons(polys, param);
        std::vector<std::tuple<T, ch::polygon>> corrected_polys =
            rv::zip(blobs | rv::transform(first<T, ch::polygon>), polys) |
            rv::transform(
                [](auto&& pair)->std::tuple<T, ch::polygon> {
                    return { pair.first, pair.second };
                }
            ) | r::to_vector;

        auto old_n = 0;
        while (old_n == corrected_polys.size()) {
            old_n = corrected_polys.size();
            corrected_polys = correct_polygons(corrected_polys);
        }

		return corrected_polys;
	}

	struct color_hasher {
		size_t operator()(const ch::color& c) const
		{
			std::size_t seed = 0;
			boost::hash_combine(seed, c[0]);
			boost::hash_combine(seed, c[1]);
			boost::hash_combine(seed, c[2]);

			return seed;
		}
	};

	using color_set = std::unordered_set<ch::color, color_hasher>;

	template<typename T>
	std::vector<std::tuple<T,ch::polygon>> raster_to_vector_aux(cv::Mat mat, double param) {
		auto colored_polygons = image_to_polygons<T>(mat);
		return simplify_colored_polygons<T>(colored_polygons, param);
	}

}

std::vector<ch::gray_polygon> ch::raster_to_vector_grayscale(cv::Mat mat, double param) {
	return raster_to_vector_aux<uchar>(mat, param);
}

std::vector<ch::colored_polygon> ch::raster_to_vector(cv::Mat mat, double param) {
    return raster_to_vector_aux<color>(mat, param);
}

std::vector<ch::point> ch::perform_douglas_peucker_simplification(
	const std::vector<ch::point>& input, double param) {

	std::vector<ch::point> output;
	bool closed = input.front() == input.back();
	cv::approxPolyDP(input, output, param, false);
	if (closed && output.front() != output.back()) {
		output.push_back(output.front());
	}
	return output;
}

std::vector<uchar> ch::detail::unique_1channel_values(const cv::Mat& input) {
	if (input.channels() != 1) {
		throw std::runtime_error("called unique_1channel_values on color image");
	}
	std::array<bool, 256> grays = {};
	input.forEach<uchar>([&grays](uchar gray, const int* pos) { grays[gray] = true;  });
	return rv::iota(0) |
		rv::take(256) |
		rv::filter([&grays](int g) {return grays[g]; }) |
		r::to<std::vector<uchar>>() |
		r::action::reverse;
}

std::vector<ch::color> ch::detail::unique_3channel_values(const cv::Mat& input) {
	if (input.channels() != 3) {
		throw std::runtime_error("called unique_3channel_values on color image");
	}
	color_set colors;
	for (int y = 0; y < input.rows; ++y) {
		for (int x = 0; x < input.cols; ++x) {
			auto pixel = input.at<ch::color>(y, x);
			colors.insert(pixel);
		}
	}
	return colors | r::to_vector;
}

std::vector<ch::gray_polygon> ch::to_monochrome(
		const std::vector<colored_polygon>& polys, bool invert) {
	auto tup_to_mono = [invert](const colored_polygon& tup)->gray_polygon {
		const auto& [col, poly] = tup;
		auto gray_val = color_to_monochrome(col);
		if (invert) {
			gray_val = 255 - gray_val;
		}
		return { gray_val , poly };
	};
	return polys | rv::transform(tup_to_mono) | r::to_vector;
}