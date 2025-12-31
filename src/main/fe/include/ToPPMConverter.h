#include <vector>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <data_type.h>
#include <model_struct.h>

class ToPPMConverter {
private:
    struct RGB {
        unsigned char r, g, b;
    };

    static RGB coolToWarmColormap(float value) {
        // value should be in range [0, 1]
        value = std::clamp(value, 0.0f, 1.0f);
        
        RGB color;
        if (value < 0.5f) {
            // Cool (blue) to white
            float t = value * 2.0f;
            color.r = static_cast<unsigned char>(255 * t);
            color.g = static_cast<unsigned char>(255 * t);
            color.b = 255;
        } else {
            // White to warm (red)
            float t = (value - 0.5f) * 2.0f;
            color.r = 255;
            color.g = static_cast<unsigned char>(255 * (1.0f - t));
            color.b = static_cast<unsigned char>(255 * (1.0f - t));
        }
        return color;
    }

public:
    static void convert(const std::string& filename, std::shared_ptr<model::ModelApi<float, int>> m_mesh, ARRAY_REAL_VIEW data, int i1, int plane,
                std::array<float, 3UL> position, int *nb_nodes) {
        int  width, height;
        size_t size = m_mesh->getNumberOfNodes();
        switch (plane)
        {
        case 0: // xy
            width = nb_nodes[0];
            height = nb_nodes[1];
            break;
        case 1: // yz
            width = nb_nodes[1];
            height = nb_nodes[2];
            break;
        case 2: // xz
            width = nb_nodes[0];
            height = nb_nodes[2];
            break;
        default:
            break;
        }

        // Find min and max for normalization
        float minVal = data(0, i1), maxVal = data(0, i1);
        for (int nodeIndex = 0; nodeIndex < size; nodeIndex++) {
            minVal = std::min(minVal, data(nodeIndex, i1));
            maxVal = std::max(maxVal, data(nodeIndex, i1));
        }

        // Write PPM file
        std::ofstream file(filename, std::ios::binary);
        file << "P6\n" << width << " " << height << "\n255\n";

        for (int nodeIndex = 0; nodeIndex < size; nodeIndex++) {
            float x = m_mesh->nodeCoord(nodeIndex, 0);
            float y = m_mesh->nodeCoord(nodeIndex, 1);
            float z = m_mesh->nodeCoord(nodeIndex, 2);
            float val;
            if (plane == 0 && z == position[2] ||
                plane == 1 && x == position[0] ||
                plane == 2 && y == position[1]) {
                    val = data(nodeIndex, i1);
            } else {
                continue;
            }
            float normalized = (maxVal == minVal) ? 0.5f : (val - minVal) / (maxVal - minVal);
            RGB color = coolToWarmColormap(normalized);
            file.put(color.r);
            file.put(color.g);
            file.put(color.b);
        }

        file.close();
    }
};