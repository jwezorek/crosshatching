#include "point_set.hpp"

ch::path_table::path_table(){
}

ch::path_table::path_ref ch::path_table::insert(
        path_ref path, path_ref path_value) {
    path_key key_1(make_path_id(path));
    auto& path_in_tbl = impl_[key_1];
    path_in_tbl.assign(path_value.begin(), path_value.end());
    return path_in_tbl;
}

ch::path_table::path_ref ch::path_table::find(const path_id& pid) const {
    path_key key_1(pid);
    auto iter = impl_.find(key_1);
    if (iter == impl_.end()) {
        return {};
    } else {
        return iter->second;
    }
}

bool ch::path_table::contains(const path_id& pid) const {
    return impl_.find(path_key(pid)) != impl_.end();
}

ch::path_table::path_id ch::make_path_id(std::span<const ch::point> path) 
{
    if (path.size() >= 3) {
        return { path.front(), path[1], path.back() };
    } else {
        if (path.size() < 2) {
            throw std::runtime_error("bad path");
        }
        return { path.front(), { 0,0 }, path.back() };
    }
}

std::tuple<ch::point, ch::point> ch::detail::sort_points(const std::tuple<ch::point, ch::point>& tup) {
    auto u = std::get<0>(tup);
    auto v = std::get<1>(tup);
    const auto& [u_x, u_y] = u;
    const auto& [v_x, v_y] = v;

    if (v_x < u_x) {
        return { v,u };
    } else if (v_x > u_x) {
        return { u,v };
    }
    return (v_y < u_y) ? std::tuple{v, u} : std::tuple{u, v};
}