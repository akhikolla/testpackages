
#include "DBN.h"

using namespace std;

DBN::DBN(int size, int n_i, int *hls, int n_o, int n_l) {
    int input_size;

    N = size;
    n_ins = n_i;
    hidden_layer_sizes = hls;
    n_outs = n_o;
    n_layers = n_l;

    sigmoid_layers = new HiddenLayer*[n_layers];
    rbm_layers = new RBM*[n_layers];

    for(int i=0; i<n_layers; i++) {
        if(i == 0) {
            input_size = n_ins;
        } else {
            input_size = hidden_layer_sizes[i-1];
        }

        sigmoid_layers[i] = new HiddenLayer(N, input_size, hidden_layer_sizes[i], NULL, NULL);

        rbm_layers[i] = new RBM(N, input_size, hidden_layer_sizes[i], sigmoid_layers[i]->W, sigmoid_layers[i]->b, NULL);
    }

    log_layer = new LogisticRegression(N, hidden_layer_sizes[n_layers-1], n_outs);
}

DBN::~DBN() {
    delete log_layer;

    for(int i=0; i<n_layers; i++) {

        delete rbm_layers[i];
    }
    delete[] sigmoid_layers;
    delete[] rbm_layers;
}


void DBN::pretrain(int **input, double lr, int k, int epochs) {
    int *layer_input;
    int prev_layer_input_size;
    int *prev_layer_input;

    int *train_X = new int[n_ins];

    for(int i=0; i<n_layers; i++) {
        for(int epoch=0; epoch<epochs; epoch++) {
            for(int n=0; n<N; n++) {
                for(int m=0; m<n_ins; m++)
                    train_X[m] = input[n][m];

                for(int l=0; l<=i; l++) {
                    if(l == 0) {
                        layer_input = new int[n_ins];
                        for(int j=0; j<n_ins; j++) layer_input[j] = train_X[j];
                    } else {
                        if(l == 1) prev_layer_input_size = n_ins;
                        else prev_layer_input_size = hidden_layer_sizes[l-2];

                        prev_layer_input = new int[prev_layer_input_size];

                        for(int j=0; j<prev_layer_input_size; j++)
                            prev_layer_input[j] = layer_input[j];

                        delete[] layer_input;

                        layer_input = new int[hidden_layer_sizes[l-1]];
                        memset(layer_input, 0, sizeof(int) * hidden_layer_sizes[l-1]);

                        sigmoid_layers[l-1]->sample_h_given_v(prev_layer_input, layer_input);

                        delete[] prev_layer_input;
                    }
                }

                rbm_layers[i]->contrastive_divergence(layer_input, lr, k);
                delete[] layer_input;
            }

        }
    }

    delete[] train_X;
    //delete[] layer_input;
}

void DBN::finetune(int **input, int **label, double lr, int epochs) {
    int *layer_input;

    int *prev_layer_input;

    int *train_X = new int[n_ins];
    int *train_Y = new int[n_outs];

    for(int epoch=0; epoch<epochs; epoch++) {
        for(int n=0; n<N; n++) {
            for(int m=0; m<n_ins; m++)  train_X[m] = input[n][m];
            for(int m=0; m<n_outs; m++) train_Y[m] = label[n][m];

            for(int i=0; i<n_layers; i++) {
                if(i == 0) {
                    prev_layer_input = new int[n_ins];
                    for(int j=0; j<n_ins; j++) prev_layer_input[j] = train_X[j];
                } else {
                    prev_layer_input = new int[hidden_layer_sizes[i-1]];
                    for(int j=0; j<hidden_layer_sizes[i-1]; j++) prev_layer_input[j] = layer_input[j];
                    delete[] layer_input;
                }


                layer_input = new int[hidden_layer_sizes[i]];
                sigmoid_layers[i]->sample_h_given_v(prev_layer_input, layer_input);
                delete[] prev_layer_input;
            }

            log_layer->train(layer_input, train_Y, lr);
            delete[] layer_input;
        }
    }

    //delete[] layer_input;
    delete[] train_X;
    delete[] train_Y;
}

void DBN::predict(int *x, double *y) {
    double *layer_input;
    double *prev_layer_input;
    double linear_output;

    prev_layer_input = new double[n_ins];
    for(int j=0; j<n_ins; j++) prev_layer_input[j] = x[j];

    for(int i=0; i<n_layers; i++) {
        layer_input = new double[sigmoid_layers[i]->n_out];

        for(int k=0; k<sigmoid_layers[i]->n_out; k++) {
            linear_output = 0.0;

            for(int j=0; j<sigmoid_layers[i]->n_in; j++) {
                linear_output += sigmoid_layers[i]->W[k][j] * prev_layer_input[j];
            }
            linear_output += sigmoid_layers[i]->b[k];
            layer_input[k] = sigmoid(linear_output);
        }
        delete[] prev_layer_input;

        if(i < n_layers-1) {
            prev_layer_input = new double[sigmoid_layers[i]->n_out];
            for(int j=0; j<sigmoid_layers[i]->n_out; j++) prev_layer_input[j] = layer_input[j];
            delete[] layer_input;
        }
    }

    for(int i=0; i<log_layer->n_out; i++) {
        y[i] = 0;
        for(int j=0; j<log_layer->n_in; j++) {
            y[i] += log_layer->W[i][j] * layer_input[j];
        }
        y[i] += log_layer->b[i];
    }

    log_layer->softmax(y);

    delete[] layer_input;
}
