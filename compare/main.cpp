#include <iostream>
#include <vector>
#include "net.hpp"
#include <cfloat>
#include <ctime>
#include <unistd.h>
#include <Eigen/Eigen>
using namespace Eigen;

vector<vector<float>> pca_output(vector<vector<float>> &data) {
    uint32_t sample = data.size(),  feature = data[0].size();
    vector<vector<float>> result(sample);
    /* Decompose the data by SVD and take the V matrix as the eigenvector of the covariance matrix */
    Matrix<float, -1, -1> m(feature, sample);
    for (long i = 0; i != sample; ++i) for (long j = 0; j != feature; ++j) m(j, i) = data[i][j];
    BDCSVD<Matrix<float, -1, -1>> svd(m, ComputeThinU | ComputeThinV);
    RowVectorXf pc1 = svd.matrixU().col(0);
    RowVectorXf pc2 = svd.matrixU().col(1);
    for (long i = 0; i != sample; ++i) {
        result[i].push_back(pc1 * m.col(i));
        result[i].push_back(pc2 * m.col(i));
    }
    return result;
}

void document(void) {
    cerr << "usage :    compare input\n";
    cerr << "\te:   learning rate = 0.001\n";
    cerr << "\tg:   max of gr = 4\n";
    exit(0);
}

int main(int argc, char *argv[]) {
    size_t t0 = time(nullptr);
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
    IOFile io(argv[optind]);     io.normalization_z();
    vector<vector<float>> pca_data;
    pca_data = pca_output(io.data);
    pca pc(io.feature, learning_rate);
    uint64_t gr = 0, count = 0, rng = wyhash64(0, 0);
    float loss0_pca = FLT_MAX, loss_pca = 0.0;
    size_t *shuff = new size_t[io.sample];  for (size_t i = 0; i < io.sample; ++i) {shuff[i] = i;}
    cerr << "data:\t" << io.sample << " * " << io.feature << endl;
    cout << unitbuf;    cout << "loss : \t";    cout.precision(5);  cout.setf(ios::fixed);
    while (gr < gr_max) {
        while (count < 0x100000) {
            for (size_t i = io.sample - 1; i; --i) swap(shuff[i], shuff[wyrand(&rng) % (i + 1)]);      // fisher shuff
            for (size_t i = 0; i < io.sample; ++i) loss_pca += pc.pca_mlp(pca_data[shuff[i]].data(), io.data[shuff[i]].data());
            count += io.sample;
        }
        if (loss_pca > loss0_pca) ++gr;
        cout << loss_pca/count/io.feature << '\t';
        loss0_pca = loss_pca;   loss_pca = 0.0, count = 0;
    }
    cout << "\ntime:\t" << time(nullptr) - t0 << "sec"<< endl;
    return 0;
}
