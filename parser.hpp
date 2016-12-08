/*ToDo:
-one easy improvement is to remove all the zeros at the beginning of the sequences
*/

#ifndef __PARSER
#define __PARSER

#include <iostream>
#include <fstream>
#include <cstring>
#include <vector>
#include <boost/program_options.hpp>
#include <boost/dynamic_bitset.hpp>

using namespace boost::program_options;
using namespace boost;
using namespace std;


class CParser{
public:
  struct Sfile_params{
    vector<int> sites_sigma_start, sites_sigma_end;
    vector<int> sites_lexA_start, sites_lexA_end;
    vector<double> energy_sigma;
    vector<double> energy_lexA;
    int seq_length;
    string output_name;
  }file_params; //Will hold the parameters extracted from motevo output file

  string output_name; //name to append to eventually outputted files

  bool verbose;
  
  vector<dynamic_bitset<>> sites_sigma;
  vector<double> weight_sites_sigma;

  vector<dynamic_bitset<>> sites_lexA;
  vector<double> weight_sites_lexA;

  
  CParser(int argc_, char **argv_): argc(argc_), argv(argv_){  //constructor: takes cmd line arguments
    return;
  }
  void start_parsing(void); //this is the main function to call to start the parsing process
    
private:
  //Comd line arguments
  int argc;
  char **argv;
  
  //These vars are read from the command line
  string filename, sigma_name, lexA_name, sigma_strand;
  int sigma_right, sigma_left, lexA_right, lexA_left;
  double thres_sigma, thres_lexA;

  //Functions to print info
  void print_cmdline(void); //print cmd line options
  void print_params(void);  //print params extracted from motevo file

  //Functions to parse the file
  void ParseOptions(void);
  void ParseFile(void); //The filename is extracted from the command line
  void CreateSites(vector<dynamic_bitset<>>& sites, vector<double>& weights, const int seq_length,
		   const vector<double>energy, const vector<int> sites_start, const vector<int> sites_end);
}; 


#endif
