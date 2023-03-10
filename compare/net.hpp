//
// Created by ljz on 22-11-18.
//

#ifndef COMPARE_NET_HPP
#define COMPARE_NET_HPP
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include "wyhash.h"
#define activate(x)     (x / (1 + (((int)(x > 0) << 1) - 1) * x))
#define act_back(x)     ((1 - (((int)(x > 0) << 1) - 1) * x) * (1 - (((int)(x > 0) << 1) - 1) * x))
using namespace std;

static inline float loss_function(const float *first, float *last, uint32_t &feature) {
    float result = 0.0, wh = 1 / sqrtf(feature);
    for (size_t i = 0; i < feature; ++i) {
        result += 0.5f * (last[i] - first[i]) * (last[i] - first[i]);
        last[i] -= first[i];
//        last[i] *= wh;          // big_data single_cell
    }
    return result;
}

struct NeuralNetwork {
    uint32_t input_neuron_num, output_neuron_num;
    float learning_rate;
    bool isLinear;
    float wh, *w, *output, *gra;
    NeuralNetwork() = default;
    NeuralNetwork(uint32_t input_num, uint32_t output_num, float rate, bool linear) :
            input_neuron_num(input_num), output_neuron_num(output_num), learning_rate(rate), isLinear(linear) {
        wh = 1 / sqrtf(input_neuron_num);
        w = new float[output_neuron_num * input_neuron_num];
        output = new float[output_neuron_num];
        gra = new float[input_neuron_num];
        uint64_t rng = wyhash64(input_neuron_num, output_neuron_num);
        for (size_t i = 0; i < output_neuron_num * input_neuron_num; ++i) w[i] = wy2gau(wyrand(&rng));     // 初始化参数
    }
//    ~NeuralNetwork() {delete [] w; delete [] output; delete [] gra;}
    void forw(const float *input) const {                                           // 前向传播
        for (size_t i = 0; i < output_neuron_num; ++i) {
            float temp = 0.0,   *weight = w + i * input_neuron_num;
            for (size_t j = 0; j < input_neuron_num; ++j) temp += input[j] * weight[j];
            output[i] = temp * wh;
        }
        if (!isLinear) {
            for (size_t i = 0; i < output_neuron_num; ++i) output[i] = activate(output[i]);
        }
    }
    void back(const float *input, const float *gin) const {                           // 后向传播
        memset(gra, 0, input_neuron_num * sizeof(float));
        if (!isLinear) {
            for (size_t i = 0; i < output_neuron_num; ++i) {
                float *weight = w + i * input_neuron_num, s = gin[i] * act_back(output[i]) * wh;
                for (size_t j = 0; j < input_neuron_num; ++j) {
                    gra[j] += s * weight[j];                            // 输入层 的偏导
                    weight[j] -= s * input[j] * learning_rate;          // 更新 w
                }
            }
        } else {
            for(size_t i = 0; i < output_neuron_num; ++i) {
                float *weight = w + i * input_neuron_num, s = gin[i] * wh;
                for (size_t j = 0; j < input_neuron_num; ++j) {
                    gra[j] += s * weight[j];                            // 输入层 的偏导
                    weight[j] -= s * input[j] * learning_rate;          // 更新 w
                }
            }
        }
    }
};

struct IOFile {
    vector<vector<float>> data;
    uint32_t sample, feature;
    const string file;
    IOFile(const string f) : file(f) { load(); }
    ~IOFile() = default;
    bool load() {
        ifstream INFILE(file);
        float temp; string line;
        vector<float> temv;
        if (!INFILE) return false;
        getline(INFILE, line);  istringstream LIN(line);
        while (LIN) {
            LIN >> temp;
            temv.push_back(temp);
        }
        temv.pop_back();
        data.push_back(temv);
        feature = temv.size();
        while (INFILE) {
            for (size_t i = 0; i != feature; ++i) {
                INFILE >> temp;
                temv[i] = temp;
            }
            data.push_back(temv);
        }
        data.pop_back();
        sample = data.size();
        return true;
    }
    void normalization_z() {                                  // z-score 归一化
        for (uint32_t i = 0; i != feature; ++i) {
            float sum_x = 0.0, sum_xx = 0.0;
            for (uint32_t j = 0; j != sample; ++j) {
                sum_x += data[j][i];
                sum_xx += data[j][i] * data[j][i];
            }
            sum_x /= sample;    sum_xx = sum_xx / sample - sum_x * sum_x;   sum_xx = 1 / sqrt(sum_xx);
            for (uint32_t j = 0; j != sample; ++j) data[j][i] = (data[j][i] - sum_x) * sum_xx;
        }
    }
};

struct pca {
    NeuralNetwork n0;
    NeuralNetwork n1;
    NeuralNetwork n2;
    pca() = default;
    pca(uint32_t feature, float rate) {
        n0 = NeuralNetwork(2, 32, rate, false);
        n1 = NeuralNetwork(32, 32, rate, false);
        n2 = NeuralNetwork(32, feature, rate,true);
    }
    void forw(float *input) {
        n0.forw(input);
        n1.forw(n0.output);
        n2.forw(n1.output);
    }
    void back(float *input, float *gin) {
        n2.back(n1.output, gin);
        n1.back(n0.output, n2.gra);
        n0.back(input, n1.gra);
    }
    float pca_mlp(float *input, float *data) {
        forw(input);
        float loss = loss_function(data, n2.output, n2.output_neuron_num);
        back(input, n2.output);
        return loss;
    }
};
#endif //COMPARE_NET_HPP
