#include "point_set.hpp"

ch::path_table::path_table(){
}

ch::path_table::path_ref ch::path_table::insert(path_ref path) {
    path_key key(path);
    auto& path_in_tbl = impl_[key];
    path_in_tbl.assign(path.begin(), path.end());
    return path_in_tbl;
}

ch::path_table::path_ref ch::path_table::find(const path_id& pid) const {
    path_key key(pid);
    auto iter = impl_.find(key);
    if (iter == impl_.end()) {
        return {};
    } else {
        return iter->second;
    }
}

bool ch::path_table::contains(const path_id& pid) const {
    return impl_.find(path_key(pid)) != impl_.end();
}