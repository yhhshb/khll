#include <stdexcept>
#include <cmath>
#include <cassert>
#include "../include/HyperLogLog.hpp"
#include "../nthash/nthash.hpp"

#include <iostream>

#define BITS_IN_BYTE 8

namespace sketching {

HyperLogLog::HyperLogLog() 
    : k(0), b(0), shift(0), mask(0), total_seen_kmers(0)
{
    sanitize_endianness();
}

HyperLogLog::HyperLogLog(const uint8_t kmer_length, const double error_rate)
    : k(kmer_length), total_seen_kmers(0)
{
    if (error_rate < 0 or error_rate > 1) throw std::invalid_argument("error rate should be in (0, 1)");
    auto x = double(1.04) / error_rate;
    auto l = static_cast<std::size_t>(std::ceil(std::log2l(x * x)));
    if (l > std::numeric_limits<uint8_t>::max()) throw std::invalid_argument("error rate too low (too many buckets)");
    b = static_cast<uint8_t>(l);
    sanitize_endianness();
    sanitize_kmer_length(k);
    sanitize_b(b);
    init();
    clear();
}

HyperLogLog::HyperLogLog(const uint8_t kmer_length, const uint8_t msb_length)
    : k(kmer_length), b(msb_length), total_seen_kmers(0)
{
    sanitize_endianness();
    sanitize_kmer_length(k);
    sanitize_b(b);
    init();
    clear();
}

HyperLogLog::HyperLogLog(std::istream& istrm)
{
    // load k, b;
    istrm.read(reinterpret_cast<char*>(&k), sizeof(k));
    sanitize_kmer_length(k);
    istrm.read(reinterpret_cast<char*>(&b), sizeof(b));
    sanitize_b(b);
    istrm.read(reinterpret_cast<char*>(&total_seen_kmers), sizeof(total_seen_kmers)); // TODO fix endianess
    init();
    // load registers
    istrm.read(reinterpret_cast<char*>(&registers[0]), registers.size()); // uint8_t so no need to endianess nor sizeof
}

void
HyperLogLog::naive_add(char const * const seq, const std::size_t length) noexcept
{
    nthash::NtHash hasher(
        seq, 
        length, 
        std::max(static_cast<std::size_t>(sizeof(hash_t) / sizeof(uint64_t)), static_cast<std::size_t>(1)), 
        k, 
        0
    );
    while(hasher.roll()) {
        hash_t hval = *reinterpret_cast<hash_t const*>(hasher.hashes());
        const auto idx = hval >> shift;
        const auto lsb = hval & mask;
        const std::size_t v = clz(lsb) + 1 - b;
        if (v > registers.at(idx)) registers[idx] = v;
        ++total_seen_kmers;
    }
}

void
HyperLogLog::buffered_add(char const * const seq, const std::size_t length, std::vector<uint64_t>& buffer) noexcept
{
    const std::size_t pack_shift = BITS_IN_BYTE * (sizeof(uint64_t) - sizeof(register_t));
    const uint64_t pack_mask = (static_cast<std::size_t>(1) << pack_shift) - 1;
    nthash::NtHash hasher(
        seq, 
        length, 
        std::max(static_cast<std::size_t>(sizeof(hash_t) / sizeof(uint64_t)), static_cast<std::size_t>(1)), 
        k, 
        0
    );
    buffer.clear();
    while(hasher.roll()) {
        const hash_t hval = *reinterpret_cast<hash_t const*>(hasher.hashes());
        const auto idx = hval >> shift;
        const auto lsb = hval & mask;
        const std::size_t v = clz(lsb) + 1 - b;
        assert(v < BITS_IN_BYTE * sizeof(hash_t));
        assert(idx < pack_mask); // idx must fit into 56 bits
        assert(v < std::numeric_limits<register_t>::max()); // v must fit into registers
        // std::cerr << uint64_t(hval >> 64) << uint64_t(hval) <<  "\n";
        buffer.push_back(v << pack_shift | idx);
    }
    for (auto pack : buffer) {
        register_t v = pack >> pack_shift;
        std::size_t idx = pack & pack_mask;
        if (v > registers.at(idx)) registers[idx] = v;
        ++total_seen_kmers;
    }
}

void 
HyperLogLog::clear() noexcept
{
    for (auto& r : registers) r = 0;
}

std::size_t
HyperLogLog::size() const noexcept
{
    return total_seen_kmers;
}

std::size_t
HyperLogLog::count() const noexcept
{
    // !!! DO NOT compute registers.size()**2 first because if b >= 32 we have an integer overflow
    auto hmean = harmonic_mean();
    const std::size_t raw_estimate = (alpha_m * hmean * registers.size()) * registers.size(); // !!!
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
    if (not compatible(other)) throw std::runtime_error("[operator+] Adding two incompatible sketches");
    HyperLogLog toRet(k, b);
    for (std::size_t i = 0; i < registers.size(); ++i) {
        toRet.registers[i] = std::max(registers.at(i), other.registers.at(i));
    }
    toRet.total_seen_kmers = total_seen_kmers + other.total_seen_kmers;
    return toRet;
}

HyperLogLog& 
HyperLogLog::operator+=(const HyperLogLog& other)
{
    if (not compatible(other)) throw std::runtime_error("[operator+=] Merging incompatible sketches");
    for (std::size_t i = 0; i < registers.size(); ++i) {
        registers[i] = std::max(registers.at(i), other.registers.at(i));
    }
    total_seen_kmers += other.total_seen_kmers;
    return *this;
}

void 
HyperLogLog::store(std::ostream& ostrm) const
{
    // save k, b, total_seen_kmers;
    ostrm.write(reinterpret_cast<const char*>(&k), sizeof(k));
    ostrm.write(reinterpret_cast<const char*>(&b), sizeof(b));
    ostrm.write(reinterpret_cast<const char*>(&total_seen_kmers), sizeof(total_seen_kmers));
    ostrm.write(reinterpret_cast<const char*>(registers.data()), registers.size());
}

void 
HyperLogLog::store(std::string const& sketch_file) const
{
    std::ofstream ostrm(sketch_file, std::ios::binary);
    return store(ostrm);
}

HyperLogLog 
HyperLogLog::load(std::string const& sketch_filename) {
    std::ifstream istrm(sketch_filename, std::ios::binary);
    return HyperLogLog(istrm);
}

void
HyperLogLog::init()
{
    registers.resize(static_cast<hash_t>(1) << b);
    alpha_m = 0.7213 / (1 + 1.079 / registers.size());
    shift = (BITS_IN_BYTE * sizeof(hash_t) - b);
    mask = (static_cast<hash_t>(1) << shift) - 1;
}

void
HyperLogLog::sanitize_endianness() const
{
    int n = 1;
    bool little_endian = (*(char *)&n == 1);
    if (not little_endian) {
        throw std::runtime_error("This implementation of HLL sketches only works on Little Endan architectures");
    }
}

void 
HyperLogLog::sanitize_kmer_length(const std::size_t kmer_length) const
{
    if (kmer_length > (BITS_IN_BYTE * sizeof(hash_t) / 2)) throw std::invalid_argument(std::string("k-mer length should be 0 < k <= ") + std::to_string(BITS_IN_BYTE * sizeof(hash_t) / 2));
}

void 
HyperLogLog::sanitize_b(const std::size_t bval) const 
{
    if (bval >= BITS_IN_BYTE * sizeof(hash_t)) throw std::invalid_argument(std::string("Number of indexing bits should be < ") + std::to_string(BITS_IN_BYTE * sizeof(hash_t)));
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
        auto toadd = 1.0 / (static_cast<hash_t>(1) << r);
        sum_of_inverses += toadd;
    }
    return 1.0 / sum_of_inverses;
}

double
HyperLogLog::bias_correction(const double raw_estimate) const noexcept
{
    if (raw_estimate <= 2.5 * registers.size()) { // linear counting
        std::size_t count = 0;
        for (auto r : registers) if (r == 0) ++count;
        if (count != 0) return registers.size() * std::log(static_cast<double>(registers.size()) / count);
    }
    if constexpr (sizeof(hash_t) == sizeof(uint32_t)) {
        if (raw_estimate > (static_cast<hash_t>(1) << 32) / 30) { // large range correction
            return -std::pow(2, 32) * std::log(1 - raw_estimate / std::pow(2, 32));
        }
    }
    return raw_estimate;
}

int 
HyperLogLog::clz(const uint32_t x) const noexcept
{
    return x ? __builtin_clz(x) : (BITS_IN_BYTE * sizeof(sketching::HyperLogLog::hash_t));
}

int 
HyperLogLog::clz(const uint64_t x) const noexcept
{
    return x ? __builtin_clzll(x) : (BITS_IN_BYTE * sizeof(sketching::HyperLogLog::hash_t));
}

int 
HyperLogLog::clz(const __uint128_t x) const noexcept
{
    uint64_t const* view = reinterpret_cast<uint64_t const*>(&x);
    if (view[1] == 0) return 64 + clz(view[0]);
    return clz(view[1]);
}

} // namespace sketching
