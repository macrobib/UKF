#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
    
    VectorXd rmse(4);
    auto est_size = estimations.size();
    if(est_size != ground_truth.size())
        std::cout<< "Estimation and ground truth size don't match."<<std::endl;
    else if(!est_size){
        std::cout<<"Estimation size is empty."<< std::endl;
    }
    else{
        for(auto i = 0; i < est_size; i++){
            VectorXd residual = estimations[i] - ground_truth[i];
            residual = residual.array() * residual.array();
            rmse += residual;
        }
        rmse = rmse/estimations.size();
        rmse = rmse.array().sqrt();
    }
    std::cout<<"RMSE Calculated."<<std::endl;
    return rmse;
}
