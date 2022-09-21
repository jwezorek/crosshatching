#pragma once

#include "geometry.hpp"
#include <unordered_map>
#include <unordered_set>
#include <boost/functional/hash.hpp>
#include <limits>

namespace ch {

    namespace detail {

        template<int P>
        class fp_value
        {
            static constexpr double calc_scaling_factor(int digits_of_precision)
            {
                return (digits_of_precision == 1) ? 10.0 : 10.0 * calc_scaling_factor(digits_of_precision - 1);
            }

            static constexpr double scaling_factor = calc_scaling_factor(P);

        public:
            fp_value(double val) :
                impl_(static_cast<int64_t>(std::llround(scaling_factor* val)))
            {}

            bool operator==(fp_value<P> fpv) const
            {
                return impl_ == fpv.impl_;
            }

            double to_float() const
            {
                return impl_ / scaling_factor;
            }

            int64_t impl() const {
                return impl_;
            }

        private:
            int64_t impl_;
        };

        template<int P>
        class path_key {

            template<int P>
            friend struct path_key_hasher;

        public:

            path_key(const ch::point& u, const ch::point& m, const ch::point& v) :
                u_x_(u.x), u_y_(u.y), m_x_(m.x), m_y_(m.y), v_x_(v.x), v_y_(v.y)
            {}

            path_key(const std::tuple<point, point, point>& pid) :
                path_key(std::get<0>(pid), std::get<1>(pid), std::get<2>(pid))
            {}

            bool operator==(const path_key<P>& other) const
            {
                return u_x_ == other.u_x_ && u_y_ == other.u_y_ &&
                    m_x_ == other.m_x_ && m_y_ == other.m_y_ &&
                    v_x_ == other.v_x_ && v_y_ == other.v_y_;
            }

        private:
            fp_value<P> u_x_;
            fp_value<P> u_y_;
            fp_value<P> m_x_;
            fp_value<P> m_y_;
            fp_value<P> v_x_;
            fp_value<P> v_y_;
        };

        template<int P>
        struct path_key_hasher
        {
            size_t operator()(const path_key<P>& key) const {
                size_t seed = 0;

                boost::hash_combine(seed, key.u_x_.impl());
                boost::hash_combine(seed, key.u_y_.impl());
                boost::hash_combine(seed, key.m_x_.impl());
                boost::hash_combine(seed, key.m_y_.impl());
                boost::hash_combine(seed, key.v_x_.impl());
                boost::hash_combine(seed, key.v_y_.impl());

                return seed;
            }
        };

        template<int P>
        class fp_point
        {
            template<int P>
            friend struct fp_point_hasher;

        public:
            fp_point(const ch::point& pt) :
                x_(pt.x), y_(pt.y)
            {}

            bool operator==(const fp_point<P>& other) const
            {
                return x_ == other.x_ && y_ == other.y_;
            }

            ch::point to_point() const {
                return { x_.to_float(), y_.to_float() };
            }

        private:
            fp_value<P> x_;
            fp_value<P> y_;
        };

        template<int P>
        struct fp_point_hasher
        {
            size_t operator()(const fp_point<P>& key) const {
                size_t seed = 0;

                boost::hash_combine(seed, key.x_.impl());
                boost::hash_combine(seed, key.y_.impl());

                return seed;
            }
        };


        template<typename V, int P>
        using path_map = std::unordered_map< path_key<P>, V, path_key_hasher<P>>;

        struct int_point_hasher
        {
            size_t operator()(const int_point& p) const
            {
                std::size_t seed = 0;
                boost::hash_combine(seed, p.x);
                boost::hash_combine(seed, p.y);

                return seed;
            }
        };

        std::tuple<ch::point, ch::point> sort_points(const std::tuple<ch::point, ch::point>& tup);

        template<int P>
        class edge_key {

            template<int P>
            friend struct edge_key_hasher;

        public:

            edge_key() :
                u_x_(-std::numeric_limits<double>::max()),
                u_y_(-std::numeric_limits<double>::max()),
                v_x_(-std::numeric_limits<double>::max()),
                v_y_(-std::numeric_limits<double>::max())
            {
            }

            edge_key(const std::tuple<point, point>& pid) :
                edge_key() {
                auto [pt1, pt2] = detail::sort_points(pid);
                u_x_ = pt1.x;
                u_y_ = pt1.y;
                v_x_ = pt2.x;
                v_y_ = pt2.y;
            }

            edge_key(const ch::point& u, const ch::point& v) :
                edge_key(std::tuple<ch::point, ch::point>{u, v}) {
            }

            bool operator==(const edge_key<P>& other) const {
                return u_x_ == other.u_x_ && u_y_ == other.u_y_ &&
                    v_x_ == other.v_x_ && v_y_ == other.v_y_;
            }

        private:
            detail::fp_value<P> u_x_;
            detail::fp_value<P> u_y_;
            detail::fp_value<P> v_x_;
            detail::fp_value<P> v_y_;
        };

        template<int P>
        struct edge_key_hasher
        {
            size_t operator()(const edge_key<P>& key) const {
                size_t seed = 0;

                boost::hash_combine(seed, key.u_x_.impl());
                boost::hash_combine(seed, key.u_y_.impl());
                boost::hash_combine(seed, key.v_x_.impl());
                boost::hash_combine(seed, key.v_y_.impl());

                return seed;
            }
        };

        struct int_edge_hasher {
            size_t operator()(const std::tuple<int, int>& ie) const {
                size_t seed = 0;
                auto [u, v] = ie;
                if (u > v) {
                    std::swap(u, v);
                }
                boost::hash_combine(seed, u);
                boost::hash_combine(seed, v);
                return seed;
            }
        };
    }

    constexpr int k_point_set_prec = 2;

    using point_set_member = detail::fp_point<k_point_set_prec>;

    using point_set = 
        std::unordered_set<point_set_member, detail::fp_point_hasher<k_point_set_prec>>;

    template <typename V>
    using point_map =
        std::unordered_map<point_set_member, V, detail::fp_point_hasher<k_point_set_prec>>;

    template<typename V>
    using int_point_map = std::unordered_map<int_point, V, detail::int_point_hasher>;

    using int_point_set = std::unordered_set<int_point, detail::int_point_hasher>;

    class path_table {
    private:
        detail::path_map<std::vector<ch::point>, k_point_set_prec> impl_;
    public:
        using path_key = detail::path_key< k_point_set_prec>;
        using path_ref = std::span<const point>;
        using path_id = std::tuple<point, point, point>;

        path_table();
        path_ref insert(path_ref path, path_ref path_value);
        path_ref find(const path_id& pid) const;
        bool contains(const path_id& pid) const;
    };

    path_table::path_id make_path_id(std::span<const point> path);


    constexpr int k_edge_table_prec = 4;
    using edge_key = detail::edge_key<k_edge_table_prec>;

    template<typename V>
    using edge_table = std::unordered_map<edge_key, V, detail::edge_key_hasher<k_edge_table_prec>>;

    using int_edge_set = std::unordered_set<std::tuple<int, int>, detail::int_edge_hasher>;

}