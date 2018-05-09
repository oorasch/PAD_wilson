#ifndef FILE_UPDATES
#define FILE_UPDATES

void do_updates(const int n, const int neib[][4], const std::vector<double>& I_bessel, double& detM, std::vector<int>& s_site, std::vector<int>& plaq_occ, std::vector<std::vector<int>>& k_link, std::vector<std::vector<int>>& sig_link);

void do_loop_updates(const int neib[][4], const std::vector<double>& I_bessel, double& detM, std::vector<int>& s_site, std::vector<int>& plaq_occ, std::vector<std::vector<int>>& k_link, std::vector<std::vector<int>>& sig_link);

void plaquette_ins_update(const int neib[][4], const std::vector<int>& coords, const std::vector<std::vector<int>>& sig_link, const std::vector<double>& I_bessel, double& detM, std::vector<int>& s_site, std::vector<int>& plaq_occ, std::vector<std::vector<int>>& k_link);

void plaquette_del_update(const int neib[][4], const std::vector<int>& coords, const std::vector<std::vector<int>>& sig_link, const std::vector<double>& I_bessel, double& detM, std::vector<int>& s_site, std::vector<int>& plaq_occ, std::vector<std::vector<int>>& k_link);

void loop_exp_update(const int neib[][4], const std::vector<int>& coords, const std::vector<std::vector<int>>& sig_link, const std::vector<double>& I_bessel_p, double& detM, std::vector<int>& s_site, std::vector<int>& plaq_occ, std::vector<std::vector<int>>& k_link);

void loop_col_update(const int neib[][4], const std::vector<int>& coords, const std::vector<std::vector<int>>& sig_link, const std::vector<double>& I_bessel, double& detM, std::vector<int>& s_site, std::vector<int>& plaq_occ, std::vector<std::vector<int>>& k_link);

void loop_rot_update(const std::vector<int>& coords, const std::vector<double>& I_bessel, std::vector<int>& plaq_occ, std::vector<std::vector<int>>& k_link);

void do_blanket_update(const std::vector<double>& I_bessel, std::vector<int>& plaq_occ);

void do_sigma_update(const int neib [][4], const std::vector<int> s_site, double detM, std::vector<std::vector<int>>& sig_link);


#endif
