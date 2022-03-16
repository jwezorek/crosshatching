#include "util.hpp"
#include <sstream>
#include <random>

namespace {
	std::random_device rd;
	std::mt19937 random(rd());
}

std::string ch::svg_header(int wd, int hgt)
{
    std::stringstream ss;

    ss << "<?xml version=\"1.0\" standalone=\"no\"?>\n";
    ss << "<svg width=\"" + std::to_string(wd) + "px\" height=\"" + std::to_string(hgt) + "px\"  xmlns = \"http://www.w3.org/2000/svg\" version = \"1.1\">\n";
    ss << "<rect width=\"" + std::to_string(wd) + "\" height=\"" + std::to_string(hgt) + "\"  fill=\"white\"/>\n";

    return ss.str();
}

double ch::normal_rnd(double mean, double stddev)
{
	std::normal_distribution<double> nd(mean, stddev);
	return nd(random);
}

double ch::uniform_rnd(double lower_bound, double upper_bound)
{
	std::uniform_real_distribution<> ud(lower_bound, upper_bound);
	return ud(random);
}

ch::rnd_fn ch::normal_rnd_fn(double mean, double stddev)
{
	return [mean, stddev]() {return normal_rnd(mean, stddev); };
}

ch::rnd_fn ch::const_rnd_fn(double val)
{
	return [val]() {return val; };
}