//
// Created by balint.galgoczi on 2021.08.10..
//

#ifndef ABC_SMC_ALGORITHM_UTIL_H
#define ABC_SMC_ALGORITHM_UTIL_H

#include <iostream>
#include <vector>
#include <random>
#include <map>
#include "../include/third-party-library/Eigen/Core"
Eigen::VectorXd getDistribution(double x_mu, double x_sigma, size_t numberOfValues);

std::vector<std::vector<double>> transpose(const std::vector<std::vector<double>> &m);

double arithmetic_mean(const std::vector<double> &vector);

int searchVector(const std::vector<std::pair<int, double>> &vec, double &item);

std::vector<std::pair<int, double>> sortByAscending(std::map<int, double>& M);

void histogram(std::vector<double> data, int nbins = 20 );

Eigen::VectorXi randomWeightedIndices(int draws, Eigen::VectorXd weights);

std::vector<double> resampling(std::vector<double> vector, std::vector<int> indices);

Eigen::MatrixXd resampling(Eigen::MatrixXd vector, Eigen::VectorXi indices);

double getUniformRandomNumber();

Eigen::VectorXd staticPeakModel(Eigen::VectorXd x, Eigen::VectorXd params);

void populationStatistics(Eigen::MatrixXd population);

void populationStatistics(std::vector<std::vector<double>> population);

//template<class T> void reorder(std::vector<T> &v, std::vector<size_t> const &order ){
//    for ( int s = 1, d; s < order.capacity(); ++ s ) {
//        for ( d = order[s]; d < s; d = order[d] ) ;
//        if ( d == s ) while ( d = order[d], d != s ) swap( v[s], v[d] );
//    }
//}

//template< typename order_iterator, typename value_iterator >
//void reorder( order_iterator order_begin, order_iterator order_end, value_iterator v )  {
//    typedef typename std::iterator_traits<value_iterator>::value_type value_t;
//    typedef typename std::iterator_traits<order_iterator>::value_type index_t;
//    typedef typename std::iterator_traits<order_iterator>::difference_type diff_t;
//
//    diff_t remaining = order_end - 1 - order_begin;
//    for ( index_t s = index_t(), d; remaining > 0; ++ s ) {
//        for ( d = order_begin[s]; d > s; d = order_begin[d] ) ;
//        if ( d == s ) {
//            -- remaining;
//            value_t temp = v[s];
//            while ( d = order_begin[d], d != s ) {
//                swap( temp, v[d] );
//                -- remaining;
//            }
//            v[s] = temp;
//        }
//    }
//}
#endif //ABC_SMC_ALGORITHM_UTIL_H
