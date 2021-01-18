#include<Rcpp.h>
#include <vector>
#include <math.h>
#include "../inst/include/dense.h"
#include "../inst/include/operateurs.h"

using namespace std;

/***********************************/
Dense::Dense(unsigned long _n_neurons, string _activation /*=sigmoid*/, double learning_rate_init_ /*= 0.01*/, bool bias_ /*= 1*/,  const string & alg /*= "sgd"*/)
    :n_neurons (_n_neurons),
    activation (_activation),
    learning_rate_init (learning_rate_init_),
    output_layer (false),
    algo (alg),
    beta_1 (0.9),
    beta_2 (0.999)
{

    /* convert bias (bool) to unsigned */
    bias = bias_?1:0;

    /* Allocate memory for the vectors of the layer */
    this->net. reserve (n_neurons);
    this->E. reserve (n_neurons);
    this->O. reserve (n_neurons);
    this->alpha. resize (n_neurons);
    this->M. resize(n_neurons);
    this->V. resize (n_neurons);
    this->W. resize (n_neurons);
    this->DeltaW. resize (n_neurons);
}

/***********************************/
void Dense::set_input_dim (vector<unsigned long> in_dim)
{
    if (in_dim. size () > 1)
    {
        Rcpp::Rcout << "Error in input dimension for Dense layer, the input should be a 1D vector.\n";
        Rcpp::stop ("\n.");
    }
    this->input_dim = in_dim[0];
    this->input.reserve (in_dim[0] + bias);
    // set uniform weights for each neuron

    // init  parameters (including adam parameters)
    for (unsigned long i = 0; i < this->n_neurons; ++i)
    {
        //this->alpha[i]. resize (input_dim + unsigned (this->bias), this->learning_rate_init);
        for (unsigned j = 0; j < this->input_dim + this->bias; ++j){
          this->alpha[i]. push_back (this->learning_rate_init);}

        this->DeltaW[i]. resize (this->input_dim + this->bias, 0.0);
        this->W[i] =  random_vector (this->input_dim + this->bias);
        this->M[i]. resize (this->input_dim + this->bias, 0.0);
        this->V[i]. resize (this->input_dim + this->bias, 0.0);
    }
}

bool Dense::contains_bias() {return bool( bias);}

bool Dense::is_output(){return output_layer;}

void Dense::set_output_layer(bool last){output_layer = last;}

vector<unsigned long> Dense::get_output_dim(){return {n_neurons};}

vector<unsigned long> Dense::get_input_dim(){return {input_dim};}
/***********************************/
MatD Dense::simulate (const MatD & input_, bool store)
{
    if (input_. size () > 1)
    {
        Rcpp::Rcout << "Input of the dense layer is not correct. Matrix of 1 row is required. \n";
        Rcpp::Rcout << "The input matrix contains " << input_. size () << ".\n";
        Rcpp::stop ("\n.");
    }

    if (input_[0]. size () != this->input_dim)
    {
        Rcpp::Rcout << "      The input of the dense layer is not correct.. \n";
        Rcpp::Rcout << "      The input dim must be: " << this->input_dim << ".\n";
        Rcpp::Rcout << "      The input line is of size: " <<input_. size () << ".\n";
        Rcpp::stop ("\n.");
    }

    VectD output;
    VectD input__;
    input__ = input_[0];

    if (this->bias)
        input__.insert (input__. begin (), 1);

    // output without activation function
    output = matrix_dot (this->W, input__);

    if (store == true)
    {
        net.clear ();
        this->input. clear ();
        net = output;
        this->input = input__;
    }
    // output with activation function
    output = vect_activation (output, activation);

    if (store == true)
    {
        this->O.clear ();
        this->O = output;
    }

    return {output};
}

/***********************************/
void Dense::computeErrors(const MatD & nextErrors)
{
    if (nextErrors. size () > 1)
    {
        Rcpp::Rcout << "Error to backpropagate to the dense layer is not correct. Matrix of 1 row is required. \n";
        Rcpp::Rcout << "The output errors matrix contains " << nextErrors. size () << ".\n";
        Rcpp::stop ("\n.");
    }

    /* nextErrors: errors of the next layer, and if the layer is an output layer, it is : (expected out - output) */
    if (nextErrors[0]. size () != this->n_neurons)
    {
        Rcpp::Rcout << "Error in computing the error, output dimensions are not correct.\n";
        Rcpp::Rcout << "Expecting " <<  this->n_neurons << " as output dimensions \n";
        Rcpp::Rcout << "While, the given errors are of size: " << nextErrors. size ();
    }

    this->E. clear ();

    VectD diff_net =  diff_activation (net, activation);

    E = matrix_dot (nextErrors[0], diff_net);
}

/***********************************/
void Dense::updateWeights (unsigned numb_iter)
{
    double momentum (0);
    double m_hat (0), v_hat (0), delta_alpha (0);
    if (algo == "sgd")
         momentum = 0.9;
    else
         momentum = 0.0;

    unsigned n_input = this->input_dim + this->bias;

    for (unsigned j = 0; j < this->n_neurons; ++j)
    {
        for (unsigned i = 0; i < n_input; ++i)
        {
            //break;
            this->DeltaW[j][i] =  momentum * this->DeltaW[j][i] + (1 - momentum) * this->input[i] * this->E[j];

            // when using the momentum, some time DeltaW[j][i] decrease a lot, so we fix 0.000001 as a min
            /*if (this->DeltaW[j][i] < 0.000001)
                this->DeltaW[j][i] = 0.0;*/

            this->W[j][i] -= alpha[j][i]  * this->DeltaW[j][i];
        }
      }

    // update learning rate
    if (algo == "adam")
    {
        for (unsigned j = 0; j < this->n_neurons; ++j)
        {
             for (unsigned i = 0; i < n_input; ++i)
                {
                    this->M[j][i] = (this->beta_1 * this->M[j][i]) + ((1 - this->beta_1) * this->DeltaW[j][i]) ;
                    this->V[j][i] = (this->beta_2 * this->V[j][i]) + ((1 - this->beta_2) * this->DeltaW[j][i] * this->DeltaW[j][i]) ;

                    m_hat = this->M[j][i]   / (1.0 - double (pow (this->beta_1, numb_iter + 1)));
                    v_hat = this->V[j][i]  / (1.0 - double (pow (this->beta_2, numb_iter + 1)));

                    // the adapted learning rate shouldn't be less that 0.00001 and sould not be greater that initial rate
                    if ((this->alpha[j][i] - delta_alpha) > 0.00001 and (this->alpha[j][i] - delta_alpha) <= learning_rate_init)
                    delta_alpha = (0.001 * m_hat) / (sqrt (v_hat) + 0.00000001);

                    // the adapted learning rate shouldn't be less that 0.001
                    if (  delta_alpha < (this->alpha[j][i] - 0.001) )
                        this->alpha[j][i] = this->alpha[j][i] - delta_alpha;
                }
          }
    }
}

/*********************************/
// return the errors to backpropagate to the previous layer
MatD Dense::get_errors ()
{
    VectD A = matrix_dot (Transpose (this->W), this->E);

    if (this->bias == 1)
        A. erase (A. begin ());
    return {A};
}

MatD Dense::get_weights(){return this->W;}

string Dense::getType(){return "dense";}
