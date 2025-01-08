#include <cmath>
#include "../include/HyperLogLog.hpp"
#include "../nthash/nthash.hpp"

HyperLogLog::HyperLogLog(uint8_t kmer_length, uint8_t msb_length)
    : k(kmer_length), b(msb_length), shift(8 * sizeof(hash_t) - b), mask((static_cast<std::size_t>(1) << shift) - 1)
{
    registers.resize(static_cast<std::size_t>(1) << b);
    for (auto& r : registers) r = 0;
    alpha_m = 0.7213 / (1 + 1.079 / registers.size());
}

void
HyperLogLog::add(char const * const seq, std::size_t length) noexcept
{
    auto clz = [](hash_t x) {return x ? __builtin_clzll(x) : 64;};
    nthash::NtHash hasher(seq, length, 1, k, 0);
    while(hasher.roll()) {
        const auto idx = hasher.hashes()[0] >> shift;
        const auto lsb = hasher.hashes()[0] & mask;
        const std::size_t v = clz(lsb) + 1 - b;
        if (v > registers.at(idx)) registers[idx] = v;
    }
}

std::size_t
HyperLogLog::count() const noexcept
{
    const std::size_t raw_estimate = alpha_m * harmonic_mean() * (registers.size() * registers.size());
    return static_cast<std::size_t>(bias_correction(raw_estimate));
}

double
HyperLogLog::standard_error() const noexcept
{
    return static_cast<double>(1.04) / sqrt(registers.size());
}

HyperLogLog
HyperLogLog::operator+(const HyperLogLog& other) const
{
    if (not compatible(other)) throw std::runtime_error("Merging incompatible sketch");
    HyperLogLog toRet(k, b);
    for (std::size_t i = 0; i < registers.size(); ++i) {
        toRet.registers[i] = std::max(registers.at(i), other.registers.at(i));
    }
    return toRet;
}

bool 
HyperLogLog::compatible(const HyperLogLog& other) const noexcept
{
    bool same_k = k == other.k;
    bool same_b = b == other.b;
    bool same_size = registers.size() == other.registers.size();
    return same_k and same_b and same_size;
}

double 
HyperLogLog::harmonic_mean() const noexcept
{
    double sum_of_inverses = 0;
    for (auto r : registers) {
        auto toadd = 1.0 / (static_cast<std::size_t>(1) << r);
        // std::cerr << "adding " << toadd << "\n";
        sum_of_inverses += toadd;
    }
    return 1.0 / sum_of_inverses;
}

double
HyperLogLog::bias_correction(double raw_estimate) const noexcept
{
    if (raw_estimate <= 2.5 * registers.size()) { // linear counting
        std::size_t count = 0;
        for (auto r : registers) if (r == 0) ++count;
        if (count != 0) return registers.size() * std::log(static_cast<double>(registers.size()) / count);
    }
    if constexpr (sizeof(hash_t) == sizeof(uint32_t)) {
        if (raw_estimate > (static_cast<uint64_t>(1) << 32) / 30) { // large range correction
            return -std::pow(2, 32) * std::log(1 - raw_estimate / std::pow(2, 32));
        }
    }
    return raw_estimate;
}
