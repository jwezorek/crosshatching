#pragma once

#include "geometry.hpp"
#include <unordered_map>
#include <unordered_set>
#include <boost/functional/hash.hpp>

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

            path_key(const ch::point& u, const ch::point& v) :
                x1_(u.x), y1_(u.y), x2_(v.x), y2_(v.y)
            {}

            path_key(const std::tuple<point, point>& pid) :
                path_key(std::get<0>(pid), std::get<1>(pid))
            {}

            path_key(const std::span<const point>& path) :
                path_key(path.front(), path.back())
            {}

            bool operator==(const path_key<P>& other) const
            {
                return x1_ == other.x1_ && y1_ == other.y1_ &&
                    x2_ == other.x2_ && y2_ == other.y2_;
            }

        private:
            fp_value<P> x1_;
            fp_value<P> y1_;
            fp_value<P> x2_;
            fp_value<P> y2_;
        };

        template<int P>
        struct path_key_hasher
        {
            size_t operator()(const path_key<P>& key) const {
                size_t seed = 0;
                boost::hash_combine(seed, key.x1_.impl());
                boost::hash_combine(seed, key.y1_.impl());
                boost::hash_combine(seed, key.x2_.impl());
                boost::hash_combine(seed, key.y2_.impl());
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
        using path_id = std::tuple<point, point>;

        path_table();
        path_ref insert(path_ref path);
        path_ref find(const path_id& pid) const;
        bool contains(const path_id& pid) const;
    };

}