#include <cmath>
#include "../include/HyperLogLog.hpp"
#include "../nthash/nthash.hpp"

HyperLogLog::HyperLogLog(uint8_t kmer_length, uint8_t msb_length)
    : k(kmer_length), b(msb_length)
{
    table.resize(static_cast<std::size_t>(1) << b);
}

void
HyperLogLog::add(char const * const seq, std::size_t length) noexcept
{
    auto clz = [](uint64_t x) {return x ? __builtin_clzll(x) : 64;};
    nthash::NtHash hasher(seq, length, 1, k, 0);
    while(hasher.roll()) {
        const auto idx = hasher.hashes()[0] >> b;
        const auto lsb = hasher.hashes()[0] & ~((static_cast<std::size_t>(1) << b) - 1);
        const auto v = clz(lsb) + 1;
        if (v > table.at(idx)) table[idx] = v;
    }
}

std::size_t
HyperLogLog::count() const noexcept
{
    const double alpha_m = 0.7213 / (1 + 1.079 / table.size());
    const std::size_t raw_estimate = alpha_m * table.size() * table.size() * harmonic_mean();
    return static_cast<std::size_t>(bias_correction(raw_estimate));
}

HyperLogLog
HyperLogLog::operator+(const HyperLogLog& other) const
{
    if (not compatible(other)) throw std::runtime_error("Merging incompatible sketch");
    HyperLogLog toRet(k, b);
    for (std::size_t i = 0; i < table.size(); ++i) {
        toRet.table[i] = std::max(table.at(i), other.table.at(i));
    }
    return toRet;
}

bool 
HyperLogLog::compatible(const HyperLogLog& other) const noexcept
{
    bool same_k = k == other.k;
    bool same_b = b == other.b;
    bool same_size = table.size() == other.table.size();
    return same_k and same_b and same_size;
}

double 
HyperLogLog::harmonic_mean() const noexcept
{
    double sum_of_inverses = 0;
    for (auto c : table) sum_of_inverses += std::pow(2, -c);
    return static_cast<double>(table.size()) / sum_of_inverses;
}

double
HyperLogLog::bias_correction(double raw_estimate) const noexcept
{
    if (raw_estimate <= 2.5 * table.size()) { // linear counting
        std::size_t count = 0;
        for (auto c : table) if (c == 0) ++count;
        if (count != 0) return table.size() * std::log(static_cast<double>(table.size()) / count);
    } else if (raw_estimate > (static_cast<uint64_t>(1) << 32) / 30) { // large range correction
        return -std::pow(2, 32) * std::log(1 - raw_estimate / std::pow(2, 32));
    } 
    return raw_estimate;
}
