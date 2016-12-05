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
  cerr << "W matrix in R format: " << endl;
  motevo.print_W_R(true);
  motevo.print_W_R(false);

  typedef std::numeric_limits< double > dbl;
  cerr  << endl << "Probability with [sigma]=1, [lexA]=1e-18: " << 
    setprecision(dbl::max_digits10) << motevo.compute_P(1,1e-18) << endl;

  return(0);
}
