// Export the initial 1000-grain orientation set used by the texture validation.
#include "cp_model/RVE.hpp"
#include <fstream>

int main()
{
    RVE rve("ver_PPP/rve_texture.txt");
    std::ofstream out("ver_PPP/texture_initial.txt");
    out << "# id qw qx qy qz\n";
    for (size_t i = 0; i < rve.num_grains(); i++) {
        Eigen::Quaterniond q = rve.get_grain(i).get_orientation();
        q.normalize();
        out << i << " " << q.w() << " " << q.x() << " " << q.y() << " " << q.z() << "\n";
    }
    return 0;
}
