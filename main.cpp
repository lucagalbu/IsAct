/*Example:
./a.out -f sites_reduced --thres_sigma -1000 --thres_lexA -10000
*/

#include <iomanip>
#include "parser.hpp"
#include "motevo.hpp"

int main(int argc, char** argv){
  CParser parser(argc, argv);
  parser.start_parsing();

  Cmotevo motevo(parser.sites_sigma, parser.sites_lexA, parser.weight_sites_sigma, parser.weight_sites_lexA);
  motevo.generate_W();

  if(parser.verbose){
    motevo.print_configurations();
    motevo.print_W();
  }

  cerr << endl << "---------------------------------------------------------" << endl;
  //cerr << "W matrix in R format: " << endl;
  motevo.print_W_R(true);
  motevo.print_W_R(false);

  //motevo.print_W();
  //motevo.save_W(parser.output_name);
  
  typedef std::numeric_limits< double > dbl;
  cerr  << endl << "Probability with [sigma]=1.5, [lexA]=0.5: " << 
    setprecision(dbl::max_digits10) << motevo.compute_P(1.5,0.5) << endl;
  
  return(0);
}
