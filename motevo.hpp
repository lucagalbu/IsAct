#include <iostream>
#ifndef __CMOTEVO
#define __CMOTEVO

#include <vector>
#include <limits.h>
#include <numeric>
#include <cmath>
#include <bitset>
#include <cstring>
#include <boost/dynamic_bitset.hpp>

using namespace boost;
using namespace std;


class Cmotevo{
public:
  //Main constructor
  Cmotevo(const vector<dynamic_bitset<>> sites_sigma_, const vector<dynamic_bitset<>> sites_lexA_,
	  const vector<double> weight_sites_sigma_, const vector<double> weight_sites_lexA_):
    sites_sigma(sites_sigma_),
    sites_lexA(sites_lexA_),
    weight_sites_sigma(weight_sites_sigma_),
    weight_sites_lexA(weight_sites_lexA_)
  {return;}
  
  //This is the main function
  void generate_W(void);
  
  //Print all the possible configurations of sigma70 and LexA, excluding the ones which are not physically possible
  void print_configurations(void);

  //Print the W  matrices
  void print_W(void);

  //Print the W matrices in a R format
  void print_W_R(bool active=true); //print activ eor inactive?

  //Print the W matrices in a Mathematica format
  void print_W_Mathematica(bool active=true);
  
  //Compute the probability to be active
  double compute_P(double C_sigma, double C_lexA);
  
private:
  vector<dynamic_bitset<>> sites_sigma, sites_lexA;
  vector<double> weight_sites_sigma, weight_sites_lexA;
	
  int num_config_sigma, num_config_lexA;
  //I need a vector since dynamic_bitset has a constructor with arguments
  vector<dynamic_bitset<>> config_sigma, config_lexA;  
  vector<int> num_sites_config_sigma, num_sites_config_lexA;
  
  double **W_active, **W_inactive;
  vector<double> weight_config_sigma, weight_config_lexA; 

  //Count the number of trailing (right) zeros
  int count_trailing_zeros(dynamic_bitset<> bits);

  //Check if two binding sites overlap and so they cannot be put in the same configuration
  bool is_overlapping(const dynamic_bitset<> config1, const dynamic_bitset<> config2);

  //Check if a given pair of configuration of sigma70 and LexA is active
  int is_active(const dynamic_bitset<> config_sigma, const dynamic_bitset<> config_lexA);
  
  //Return the number of valid configurations found
  void get_configurations(const vector<dynamic_bitset<>> sites, vector<double>& weight_config,
			 const int num_sites, vector<dynamic_bitset<>>& configurations,
			 const vector<double> weight_sites, vector<int> &num_sites_in_config);
}; 


#endif
