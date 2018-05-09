#ifndef FILE_DEBUG
#define FILE_DEBUG

void pauli_check(const int neib[][4], const std::vector<std::vector<int>>& k_link);

void flux_check(const int neib[][4], const std::vector<std::vector<int>>& k_link);

void plaquette_check(const int neib[][4], const std::vector<std::vector<int>>& k_link, const std::vector<int>& plaq_occ);

#endif
