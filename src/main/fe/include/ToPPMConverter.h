#include <vector>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <data_type.h>

class ToPPMConverter {
private:
    struct RGB {
        unsigned char r, g, b;
    };

    RGB coolToWarmColormap(float value) {
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
    void convert(std::shared_ptr<model::ModelApi<float, int>> m_mesh, ARRAY_REAL_VIEW data, const std::string& filename) {
        int height = data.size();
        int width = data[0].size();

        // Find min and max for normalization
        float minVal = data[0][0], maxVal = data[0][0];
        for (const auto& row : data) {
            for (float val : row) {
                minVal = std::min(minVal, val);
                maxVal = std::max(maxVal, val);
            }
        }

        // Write PPM file
        std::ofstream file(filename, std::ios::binary);
        file << "P6\n" << width << " " << height << "\n255\n";

        for (const auto& row : data) {
            for (float val : row) {
                float normalized = (maxVal == minVal) ? 0.5f : (val - minVal) / (maxVal - minVal);
                RGB color = coolToWarmColormap(normalized);
                file.put(color.r);
                file.put(color.g);
                file.put(color.b);
            }
        }

        file.close();
    }
};