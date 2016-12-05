#include "motevo.hpp"

void Cmotevo::generate_W(void){
  const int num_sites_sigma = sites_sigma.size(); //number of sites
  num_config_sigma = 1<<num_sites_sigma; //2^num_sites
  config_sigma.reserve(num_sites_sigma);
  //weight_config_sigma = new double[num_config_sigma];
  get_configurations(sites_sigma, weight_config_sigma, num_sites_sigma,
		     config_sigma, weight_sites_sigma, num_sites_config_sigma);

  const int num_sites_lexA = sites_lexA.size(); //number of sites
  num_config_lexA = 1<<num_sites_lexA; //2^num_sites
  config_lexA.reserve(num_sites_lexA);
  //weight_config_lexA = new double[num_config_lexA];
  get_configurations(sites_lexA, weight_config_lexA, num_sites_lexA,
		     config_lexA, weight_sites_lexA, num_sites_config_lexA);


  //ALLOCATE THE W MATRICES
  W_active = new double*[num_sites_sigma+1];
  for(int i=0; i<num_sites_sigma+1; i++) W_active[i] = new double[num_sites_lexA+1];

  W_inactive = new double*[num_sites_sigma+1];
  for(int i=0; i<num_sites_sigma+1; i++) W_inactive[i] = new double[num_sites_lexA+1];

  //Set elements of both W to zero
  for(int i=0; i<num_sites_sigma+1; i++){
    for(int j=0; j<num_sites_lexA+1; j++){
      W_active[i][j]=1;
      W_inactive[i][j]=1;
    }
  }


  //LOOP OVER ALL POSSIBLE PAIR OF CONFIGURATIONS OF SIGMA AND LEXA AND FILL W
  int current_num_sigma, current_num_lexA; //number of sigmas or lexAs bound in the current config
  for(int i=0; i<config_sigma.size(); i++){
    current_num_sigma = num_sites_config_sigma[i];
    for(int j=0; j<config_lexA.size(); j++){
      current_num_lexA = num_sites_config_lexA[j];
      //cout << i << "   " << j << endl;
      if(is_active(config_sigma[i], config_lexA[j])==1){
	//cout << "Active with  " << current_num_sigma << " sigma and " << current_num_lexA << " lexA" << endl;
	//cout << config_sigma[i] <<  endl;
	//cout << config_lexA[j]  << endl;
	//cout << "Weight: " << weight_config_sigma[i]*weight_config_lexA[j] << endl << endl;
	W_active[current_num_sigma][current_num_lexA]*=weight_config_sigma[i]*weight_config_lexA[j];
      }
      else if(is_active(config_sigma[i], config_lexA[j])==0){
	//cout << "Inactive with  " << current_num_sigma << " sigma and " << current_num_lexA << " lexA" << endl;
	//cout << config_sigma[i] <<  endl;
	//cout << config_lexA[j]  << endl;
	//cout << "Weight: " << weight_config_sigma[i]*weight_config_lexA[j] << endl << endl;
	W_inactive[current_num_sigma][current_num_lexA]*=weight_config_sigma[i]*weight_config_lexA[j];
      }
      //else cout << "Sigma70 and LexA overlapping" << endl;
    }
  }

  //If some elements of W are 1 it means that no sites have been multiplied, so set them to 0
   //Set elements of both W to zero
  for(int i=0; i<num_sites_sigma+1; i++){
    for(int j=0; j<num_sites_lexA+1; j++){
      if(abs(W_active[i][j]-1)<0.0000001) W_active[i][j]=0;
      if(abs(W_inactive[i][j]-1)<0.0000001) W_inactive[i][j]=0;
    }
  }


}



inline int Cmotevo::count_trailing_zeros(dynamic_bitset<> bits){
  if(bits.any()==0) return(bits.size());

  int count = 0;
  while( bits.test(0)==0 ){
    bits >>= 1;
    count++;
  }

  return(count);
}


inline bool Cmotevo::is_overlapping(const dynamic_bitset<> config1, const dynamic_bitset<> config2){
  if( (config1 & config2).any() ) return(true);
  return(false);
}


inline int Cmotevo::is_active(const dynamic_bitset<> config_sigma, const dynamic_bitset<> config_lexA){
  if(is_overlapping(config_sigma, config_lexA)){
    return(-1); //invalid configuration
  }

  //I create the mask to find the active sigmas (the ones downstream rightmost lexA)
  int tr_zeros_lexA = count_trailing_zeros(config_lexA);
  dynamic_bitset<> mask_active =  (dynamic_bitset<>(config_lexA.size()).set() << tr_zeros_lexA).flip();

  //The configuration is active if the masked config_sigma is different from 0
  //(i.e. if some site have survived the mask, i.e. if some sites are downstream lexA)
  if( (config_sigma & mask_active).any() )   return(1);
  return(0);
}


void Cmotevo::get_configurations(const vector<dynamic_bitset<>> sites, vector<double>& weight_config,
				const int num_sites, vector<dynamic_bitset<>>& configurations,
				const vector<double> weight_sites, vector<int> &num_sites_in_config){
  const int num_config = 1<<num_sites; //2^num_sites = total number of combinations for the sites (allowed and not)
  const size_t length_config = sites[0].size(); //length of one configuration
  vector<dynamic_bitset<>> seq;

  //Create the combinations pattern
  for(int i=0; i<num_config; i++) seq.emplace_back(num_sites, i);

  dynamic_bitset<> bit(length_config);
  dynamic_bitset<> config(length_config);
  bool is_valid; //used to know if a config is valid and so must be stored
  double weight_config_tmp; //stores weight of the configuration before checking if the config is valid or not

  for(int j=0; j<num_config; j++){
    config.reset();
    is_valid=true;
    weight_config_tmp = 1;
    for(int i=0; i<num_sites; i++){
      bit = (seq[j][i]) ? dynamic_bitset<>(length_config).flip() : dynamic_bitset<>(length_config);
      if( is_overlapping(config, bit & sites[i]) ){ // if the binding sites overlap, don't create the config
	is_valid=false;
	break;
      }
      config |= (bit & sites[i]);
      weight_config_tmp *= pow(weight_sites[i], seq[j][i]); //Weight is added only if bit!=0
    }
    if(is_valid){
      configurations.push_back(config);
      num_sites_in_config.push_back(seq[j].count());
      weight_config.push_back(weight_config_tmp);
    }
  }

  return;
}


void Cmotevo::print_configurations(void){
  cerr << endl << "---------------------------------------------------------" << endl;
  cerr << "SIGMA70 SITES and weights: " << endl;
  for(int i=0; i<sites_sigma.size(); i++)
    cerr << sites_sigma[i] << " - " << weight_sites_sigma[i] << endl;
  cerr << "LexA SITES and weights: " << endl;
  for(int i=0; i<sites_lexA.size(); i++)
    cerr << sites_lexA[i] << " - " << weight_sites_lexA[i] << endl;

  cerr << endl << "SIGMA70 CONFIGURATIONS, weights and number of sites in the config:" << endl;
  for(int i=0; i<config_sigma.size(); i++)
    cerr << config_sigma[i] << " - " << weight_config_sigma[i] << " - " << num_sites_config_sigma[i] << endl;
  cerr << endl << "LexA CONFIGURATIONS and weights and number of sites in the config:" << endl;
  for(int i=0; i<config_lexA.size(); i++)
    cerr << config_lexA[i] << " - " << weight_config_lexA[i] << " - " << num_sites_config_lexA[i] << endl;
  cerr << "---------------------------------------------------------" << endl;

}


void Cmotevo::print_W(void){
  cerr << endl << "---------------------------------------------------------" << endl;
  cerr << "W MATRICES" << endl;
  cerr << "Active:" << endl;
  for(int i=0; i<sites_sigma.size()+1; i++){
    for(int j=0; j<sites_lexA.size()+1; j++){
      cout << W_active[i][j] << "   ";
    }
    cout << endl;
  }

  cout << endl;
  cerr << "Inactive:" << endl;
  for(int i=0; i<sites_sigma.size()+1; i++){
    for(int j=0; j<sites_lexA.size()+1; j++){
      cout << W_inactive[i][j] << "   ";
    }
    cout << endl;
  }
  cerr << "---------------------------------------------------------" << endl;
}


void Cmotevo::print_W_R(bool active){
  if(active){
    cout << "W_active <- matrix(c(" << W_active[0][0];
    for(int j=0; j<sites_lexA.size()+1; j++){
      for(int i=0; i<sites_sigma.size()+1; i++){
	if(j==0 & i==0) continue;
	cout << ", " << W_active[i][j];
      }
    }
  }
  else{
    cout << "W_inactive <- matrix(c(" << W_inactive[0][0];
    for(int j=0; j<sites_lexA.size()+1; j++){
      for(int i=0; i<sites_sigma.size()+1; i++){
	if(j==0 & i==0) continue;
	cout << ", " << W_inactive[i][j];
      }
    }
  }

  cout << "), nrow=" << sites_sigma.size()+1 << ")" << endl;

}

void Cmotevo::print_W_Mathematica(bool active){
  if(active){
    cout << "W_active = {{" << W_active[0][0];
    for(int j=0; j<sites_lexA.size()+1; j++){
      for(int i=0; i<sites_sigma.size()+1; i++){
	if(j==0 & i==0) continue;
	cout << "," << W_active[i][j];
      }
      cout << "},{";
    }
  }
  else{
    cout << "W_inactive = {{" << W_inactive[0][0];
    for(int j=0; j<sites_lexA.size()+1; j++){
      for(int i=0; i<sites_sigma.size()+1; i++){
	if(j==0 & i==0) continue;
	cout << "," << W_inactive[i][j];
      }
      cout << "},{";
    }
  }

  cout << "}" << endl;

}

double Cmotevo::compute_P(double C_sigma, double C_lexA){
  //Create vectors with powers of concentrations
  vector<double> sigma, lexA;
  for(int i=0; i<sites_sigma.size()+1; i++) sigma.push_back(pow(C_sigma, i));
  for(int i=0; i<sites_lexA.size()+1; i++) lexA.push_back(pow(C_lexA, i));

  vector<double> WLexA_active, WLexA_inactive;
  for(int i=0; i<sites_sigma.size()+1; i++){
      double tmp_active=0, tmp_inactive=0;
    for(int j=0; j<sites_lexA.size(); j++){
      tmp_active+= W_active[i][j]*lexA[j];
      tmp_inactive+= W_inactive[i][j]*lexA[j];
    }
    WLexA_active.push_back(tmp_active);
    WLexA_inactive.push_back(tmp_inactive);
  }

  double active=0, inactive=0;
  for(int i=0; i<sites_sigma.size()+1; i++){
    active+=WLexA_active[i]*sigma[i];
    inactive+=WLexA_inactive[i]*sigma[i];
  }

  return( exp( log(active)-log(active+inactive)) );


}

