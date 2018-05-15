#ifndef FILE_DEBUG
#define FILE_DEBUG

void site_occupation_check(const std::vector<std::vector<int>>& neib, const std::vector<int>& s_site, const std::vector<std::vector<int>>& k_link);

void pauli_check(const std::vector<std::vector<int>>& neib, const std::vector<std::vector<int>>& k_link);

void flux_check(const std::vector<std::vector<int>>& neib, const std::vector<std::vector<int>>& k_link);

void plaquette_check(const std::vector<std::vector<int>>& neib, const std::vector<std::vector<int>>& k_link, const std::vector<int>& plaq_occ);

#endif
