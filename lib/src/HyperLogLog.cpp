#include <stdexcept>
#include <cmath>
#include "../include/HyperLogLog.hpp"
#include "../nthash/nthash.hpp"

namespace sketching {

HyperLogLog::HyperLogLog() 
    : k(0), b(0), shift(0), mask(0), total_seen_kmers(0)
{}

HyperLogLog::HyperLogLog(uint8_t kmer_length, double error_rate)
    : k(kmer_length), total_seen_kmers(0)
{
    if (error_rate < 0 or error_rate > 1) throw std::invalid_argument("error rate should be in (0, 1)");
    auto x = double(1.04) / error_rate;
    auto l = static_cast<std::size_t>(std::ceil(std::log2l(x * x)));
    if (l > std::numeric_limits<uint8_t>::max()) throw std::invalid_argument("error rate too low (too many buckets)");
    b = static_cast<uint8_t>(l);
    sanitize_kmer_length(k);
    sanitize_b(b);
    init();
    clear();
}

HyperLogLog::HyperLogLog(uint8_t kmer_length, uint8_t msb_length)
    : k(kmer_length), b(msb_length), total_seen_kmers(0)
{
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
HyperLogLog::add(char const * const seq, std::size_t length) noexcept
{
    auto clz = [](hash_t x) {return x ? __builtin_clzll(x) : 64;};
    nthash::NtHash hasher(seq, length, 1, k, 0);
    while(hasher.roll()) {
        const auto idx = hasher.hashes()[0] >> shift;
        const auto lsb = hasher.hashes()[0] & mask;
        const std::size_t v = clz(lsb) + 1 - b;
        if (v > registers.at(idx)) registers[idx] = v;
        ++total_seen_kmers;
    }
}

void
HyperLogLog::add_fast(char const * const seq, std::size_t length, std::vector<hash_t>& buffer) noexcept
{
    auto clz = [](hash_t x) {return x ? __builtin_clzll(x) : 64;};
    nthash::NtHash hasher(seq, length, 1, k, 0);
    buffer.clear();
    while(hasher.roll()) buffer.push_back(hasher.hashes()[0]);
    for (auto hash : buffer) {
        const auto idx = hash >> shift;
        const auto lsb = hash & mask;
        const std::size_t v = clz(lsb) + 1 - b;
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
    ostrm.write(reinterpret_cast<const char*>(&total_seen_kmers), sizeof(total_seen_kmers)); // TODO fix endianess
    ostrm.write(reinterpret_cast<const char*>(registers.data()), registers.size()); // uint8_t so no need to endianess nor sizeof
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
    registers.resize(static_cast<std::size_t>(1) << b);
    alpha_m = 0.7213 / (1 + 1.079 / registers.size());
    shift = (8 * sizeof(hash_t) - b);
    mask = (static_cast<std::size_t>(1) << shift) - 1;
}

void 
HyperLogLog::sanitize_kmer_length(std::size_t kmer_length) const
{
    if (kmer_length > 32) throw std::invalid_argument("k-mer length should be 0 < k <= 32");
}

void 
HyperLogLog::sanitize_b(std::size_t bval) const 
{
    if (bval >= 64) throw std::invalid_argument("Number of indexing bits should be < 64");
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

} // namespace sketching
