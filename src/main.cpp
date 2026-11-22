#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
constexpr double pi = 3.14159265358979323846;
constexpr double radius_max = 14;
constexpr double separation = 0.4;
constexpr double shift_max = 0.1;
constexpr double impact_max = 20;
constexpr std::size_t nucleons = 208;
constexpr std::size_t protons = 82;
constexpr unsigned radial_intervals = 4096;

struct Point { double x, y, z; };
using Nucleus = std::array<Point, nucleons>;
struct Event { std::uint64_t trial; double b; unsigned npart, ncoll; };

// Fix the conversion to a uniform double rather than relying on a
// standard-library distribution's engine consumption.
double uniform(std::mt19937_64& random) {
    return static_cast<double>(random() >> 11) * 0x1.0p-53;
}

class RadialProfile {
    std::array<double, radial_intervals + 1> cumulative_{};

public:
    RadialProfile(double half_density, double diffuseness) {
        const auto density = [=](double r) {
            const double correction = 1.00866 - 0.000461484*r - 0.000203571*r*r;
            return r*r / ((1 + std::exp((r-half_density)/diffuseness)) * correction);
        };
        constexpr double step = radius_max / radial_intervals;
        for (unsigned i = 1; i <= radial_intervals; ++i)
            cumulative_[i] = cumulative_[i-1] + step*(density(i*step)+density((i-1)*step))/2;
    }

    double sample(std::mt19937_64& random) const {
        const double value = uniform(random) * cumulative_.back();
        const auto upper = std::upper_bound(cumulative_.begin(), cumulative_.end(), value);
        const auto i = upper - cumulative_.begin();
        const double frac = (value-cumulative_[i-1])/(cumulative_[i]-cumulative_[i-1]);
        return (i-1+frac)*radius_max/radial_intervals;
    }
};

class Model {
    RadialProfile proton_{6.68, 0.447};
    RadialProfile neutron_{6.69, 0.56};
    std::mt19937_64 random_;
    double collision_distance2_;

    Nucleus sample_nucleus() {
        // Use a fresh random species assignment per generated nucleus.
        std::array<bool, nucleons> is_proton{};
        std::fill_n(is_proton.begin(), protons, true);
        std::shuffle(is_proton.begin(), is_proton.end(), random_);
        Nucleus points{};
        for (unsigned retry = 0; retry < 100000; ++retry) {
            Point sum{};
            for (std::size_t i = 0; i < nucleons; ++i) {
                const auto& profile = is_proton[i] ? proton_ : neutron_;
                for (;;) {
                    const double r = profile.sample(random_);
                    const double cosine = 2*uniform(random_)-1;
                    const double phi = 2*pi*uniform(random_);
                    const double transverse = r*std::sqrt(1-cosine*cosine);
                    Point p{transverse*std::cos(phi), transverse*std::sin(phi), r*cosine};
                    const auto overlaps = [&](Point previous) {
                        const double dx = p.x-previous.x;
                        if (dx*dx >= separation*separation) return false;
                        const double dy = p.y-previous.y, dz = p.z-previous.z;
                        return dx*dx+dy*dy+dz*dz < separation*separation;
                    };
                    if (std::any_of(points.begin(), points.begin()+i, overlaps)) continue;
                    points[i] = p;
                    sum.x += p.x;
                    sum.y += p.y;
                    sum.z += p.z;
                    break;
                }
            }
            const Point centre{sum.x/nucleons, sum.y/nucleons, sum.z/nucleons};
            if (centre.x*centre.x+centre.y*centre.y+centre.z*centre.z > shift_max*shift_max)
                continue;
            for (auto& p : points) {
                p.x -= centre.x;
                p.y -= centre.y;
                p.z -= centre.z;
            }
            return points;
        }
        throw std::runtime_error("could not generate a recentered Pb nucleus");
    }

public:
    Model(std::uint64_t seed, double sigma_mb)
        : random_(seed), collision_distance2_(sigma_mb/(10*pi)) {}

    Event generate(std::uint64_t trial) {
        const double b = impact_max*std::sqrt(uniform(random_));
        const auto left = sample_nucleus(), right = sample_nucleus();
        std::array<bool, nucleons> left_part{}, right_part{};
        unsigned pairs = 0;
        for (std::size_t i = 0; i < nucleons; ++i) {
            const double shifted_x = left[i].x-b;
            for (std::size_t j = 0; j < nucleons; ++j) {
                const double dx = shifted_x-right[j].x, dy = left[i].y-right[j].y;
                if (dx*dx+dy*dy >= collision_distance2_) continue;
                ++pairs;
                left_part[i] = right_part[j] = true;
            }
        }
        const unsigned participants = std::accumulate(left_part.begin(), left_part.end(), 0u)
                                    + std::accumulate(right_part.begin(), right_part.end(), 0u);
        return {trial, b, participants, pairs};
    }
};

std::uint64_t integer(const std::string& text) {
    std::uint64_t value = 0;
    const auto end = text.data()+text.size();
    const auto parsed = std::from_chars(text.data(), end, value);
    if (parsed.ptr != end || parsed.ec != std::errc{}) throw std::invalid_argument("invalid unsigned integer: "+text);
    return value;
}

std::ofstream output(const std::filesystem::path& path) {
    std::ofstream stream(path);
    stream.exceptions(std::ios::badbit | std::ios::failbit);
    stream << std::setprecision(12);
    return stream;
}

void write_run(const std::filesystem::path& directory, const std::vector<Event>& events,
               std::uint64_t trials, std::uint64_t seed, double sigma, unsigned width) {
    const auto count = events.size();
    const unsigned bins = 100/width;
    if (count < 2*bins) throw std::runtime_error("too few inelastic events for centrality summaries");
    if (!std::filesystem::create_directory(directory))
        throw std::runtime_error("output directory already exists");
    auto event_csv = output(directory/"events.csv");
    auto summary_csv = output(directory/"centrality.csv");
    event_csv << "trial,b_fm,npart,ncoll,taa_mb_inverse,centrality_low_pct,centrality_high_pct\n";
    summary_csv << "centrality_low_pct,centrality_high_pct,events,b_low_fm,b_high_fm,b_mean_fm,"
                   "npart_mean,ncoll_mean,taa_mean_mb_inverse,taa_fixed_bin_se_mb_inverse\n";
    for (unsigned bin = 0; bin < bins; ++bin) {
        const auto first = count*bin/bins, last = count*(bin+1)/bins;
        const auto lo = bin*width, hi = lo+width;
        const double lower = first ? (events[first-1].b+events[first].b)/2 : 0;
        const double upper = last < count ? (events[last-1].b+events[last].b)/2 : impact_max;
        double sum_b = 0, sum_part = 0, sum_pairs = 0, square_pairs = 0;
        for (std::size_t i = first; i < last; ++i) {
            const auto& e = events[i];
            event_csv << e.trial << ',' << e.b << ',' << e.npart << ',' << e.ncoll << ','
                      << e.ncoll/sigma << ',' << lo << ',' << hi << '\n';
            sum_b += e.b;
            sum_part += e.npart;
            sum_pairs += e.ncoll;
            square_pairs += static_cast<double>(e.ncoll)*e.ncoll;
        }
        const double n = static_cast<double>(last-first);
        const double se = std::sqrt(std::max(0., square_pairs-sum_pairs*sum_pairs/n)/(n*(n-1)))/sigma;
        summary_csv << lo << ',' << hi << ',' << last-first << ',' << lower << ',' << upper << ','
                    << sum_b/n << ',' << sum_part/n << ',' << sum_pairs/n << ',' << sum_pairs/(n*sigma)
                    << ',' << se << '\n';
    }
    event_csv.close();
    summary_csv.close();
    const double fraction = static_cast<double>(count)/trials;
    const double area_barn = pi*impact_max*impact_max/100;
    auto metadata = output(directory/"metadata.json");
    metadata << "{\n"
             << "  \"schema\": \"glauber-pbpb-v1\",\n"
             << "  \"system\": \"Pb208-Pb208\",\n"
             << "  \"profile\": \"Pbpnrw, continuous placement\",\n"
             << "  \"rng\": \"mt19937_64\",\n"
             << "  \"seed\": " << seed << ",\n"
             << "  \"trials\": " << trials << ",\n"
             << "  \"inelastic_events\": " << count << ",\n"
             << "  \"sigma_nn_inelastic_mb\": " << sigma << ",\n"
             << "  \"sigma_pbpb_inelastic_barn\": " << area_barn*fraction << ",\n"
             << "  \"sigma_pbpb_inelastic_se_barn\": " << area_barn*std::sqrt(fraction*(1-fraction)/trials) << ",\n"
             << "  \"proton_radius_fm\": 6.68,\n  \"proton_diffuseness_fm\": 0.447,\n"
             << "  \"neutron_radius_fm\": 6.69,\n  \"neutron_diffuseness_fm\": 0.56,\n"
             << "  \"radial_reweight_divisor_coefficients\": [1.00866, -0.000461484, -0.000203571],\n"
             << "  \"radial_cutoff_fm\": " << radius_max << ",\n"
             << "  \"radial_grid_intervals\": " << radial_intervals << ",\n"
             << "  \"min_nucleon_separation_fm\": " << separation << ",\n"
             << "  \"max_pre_recenter_shift_fm\": " << shift_max << ",\n"
             << "  \"max_impact_parameter_fm\": " << impact_max << ",\n"
             << "  \"centrality_bin_width_pct\": " << width << ",\n"
             << "  \"centrality_definition\": \"impact-parameter quantiles of accepted inelastic events\",\n"
             << "  \"taa_definition\": \"Ncoll divided by sigma_nn_inelastic_mb\"\n}\n";
    metadata.close();
}
} // namespace

int main(int argc, char** argv) {
    if (argc != 4 && argc != 6) {
        std::cerr << "Usage: glauber-pbpb OUTPUT_DIR TRIALS SEED [SIGMA_NN_MB BIN_WIDTH_PCT]\n";
        return 2;
    }
    try {
        const auto trials = integer(argv[2]), seed = integer(argv[3]);
        std::size_t parsed = 0;
        const double sigma = argc == 6 ? std::stod(argv[4], &parsed) : 67.6;
        const auto width = argc == 6 ? integer(argv[5]) : 5;
        if (trials < 2 || trials > 10000000 || !std::isfinite(sigma) || sigma <= 0
            || (argc == 6 && parsed != std::string(argv[4]).size())
            || !width || width > 100 || 100%width)
            throw std::invalid_argument("use 2–10000000 trials, positive sigma and a bin width dividing 100");
        if (std::filesystem::exists(argv[1])) throw std::invalid_argument("output directory already exists");
        Model model(seed, sigma);
        std::vector<Event> events;
        for (std::uint64_t i = 1; i <= trials; ++i) {
            auto event = model.generate(i);
            if (event.ncoll) events.push_back(event);
        }
        std::sort(events.begin(), events.end(), [](const Event& a, const Event& b) { return a.b < b.b; });
        write_run(argv[1], events, trials, seed, sigma, static_cast<unsigned>(width));
        std::cout << "Wrote " << events.size() << " inelastic events from " << trials << " proposals\n";
    } catch (const std::exception& error) {
        std::cerr << "glauber-pbpb: " << error.what() << '\n';
        return 1;
    }
}
