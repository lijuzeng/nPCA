#include <iostream>
#include <vector>
#include <algorithm>
#include "npca.hpp"
#include <cfloat>
#include <ctime>
#include <unistd.h>
#include <Eigen/Eigen>
using namespace Eigen;

void npca_output(float *w, vector<vector<float>> &data, ofstream &out) {
    uint32_t sample = data.size(), feature = data[0].size();
    float temp = 0.0;
    for (uint32_t i = 0; i < sample; ++i) {
        for (uint32_t j = 0; j < 2; ++j) {
            float *weight = w + j * feature;
            for (uint32_t z = 0; z < feature; ++z) {
                temp += data[i][z] * weight[z];
            }
            temp *= 1 / sqrtf(feature);
            out << temp << (j ? '\n' : '\t');
            temp = 0.0;
        }
    }
}

void pca_output(vector<vector<float>> &data, float *w, ofstream &out) {
    uint32_t sample = data.size(),  feature = data[0].size();
    /* Decompose the data by SVD and take the V matrix as the eigenvector of the covariance matrix */
    Matrix<float, -1, -1> m(feature, sample);
    for (size_t i = 0; i < sample; ++i) for (size_t j = 0; j < feature; ++j) m(j, i) = data[i][j];
    BDCSVD<Matrix<float, -1, -1>> svd(m, ComputeThinU | ComputeThinV);
    RowVectorXf pc1 = svd.matrixU().col(0);
    RowVectorXf pc2 = svd.matrixU().col(1);
    for (size_t i = 0; i < sample; ++i) out << pc1 * m.col(i) << '\t' << pc2 * m.col(i) << '\n';
    for (size_t i = 0; i < feature; ++i) {w[i] = pc1(0, i); w[i + feature] = pc2(0, i);}
}

void document(void) {
    cerr << "usage :    npca input [option]\n";
    cerr << "\t-e   learning rate = 0.001\n";
    cerr << "\t-g   max of gr = 4\n";
    exit(0);
}

int main(int argc, char *argv[]) {
    size_t t0 = time(nullptr);
    string ofile("./");
    int opt;
    float learning_rate = 0.001;
    uint64_t gr_max = 4;
    while ((opt = getopt(argc, argv, "e:g:")) >= 0) {
        switch (opt) {
            case 'e' : learning_rate = atof(optarg);    break;
            case 'g' : gr_max = atoi(optarg);   break;
            default : document();
        }
    }
    if (argc < optind + 1) document();
    ofstream out;
    IOFile io(argv[optind]);     io.normalization_z();   npca npc(io.feature, learning_rate);
    string pca_file = ofile + "pc.pca";     out.open(pca_file);      pca_output(io.data, npc.n0.w, out);   out.close();
    uint64_t gr = 0, count = 0, rng = wyhash64(0, 0);
    float loss0_npca = FLT_MAX, loss_npca = 0.0;
    size_t *shuff = new size_t[io.sample];
    for (size_t i = 0; i < io.sample; ++i) shuff[i] = i;
    cerr << "data:\t" << io.sample << " * " << io.feature << "\nlearning_rate = " << learning_rate << endl;
    cout << unitbuf;    cout << "loss : \t";    cout.precision(5);  cout.setf(ios::fixed);
    while (gr < gr_max) {
        while (count < 0x100000) {
            for (size_t i = io.sample - 1; i; --i) swap(shuff[i], shuff[wyrand(&rng) % (i + 1)]);      // fisher shuff
            for (size_t i = 0; i < io.sample; ++i) loss_npca += npc.npca_mlp(io.data[shuff[i]].data());
            count += io.sample;
        }
        if (loss_npca > loss0_npca) ++gr;
        cout << loss_npca/count/io.feature << '\t';
        loss0_npca = loss_npca;   loss_npca = 0.0;  count = 0;
    }
    delete [] shuff;
    string npca_file = ofile + "pc.npca";       out.open(npca_file);     npca_output(npc.n0.w, io.data, out);   out.close();
    cout << "\ntime:\t" << time(nullptr) - t0 << "sec"<< endl;
    return 0;
}
