#include "parser.hpp"

void CParser::ParseOptions(void){
  try{
    options_description desc{"Options"};
    desc.add_options()
      ("help,h", "Help screen")
      ("file,f", value<string>(&filename)->required(), "Filename with the results from Motevo")
      ("sigma", value<string>(&sigma_name)->default_value("Sigma70_spacer"), "Identifier of sigma70")
      ("lexA", value<string>(&lexA_name)->default_value("LexA"), "Identifier of LexA")
      ("verbose", value<bool>(&verbose)->default_value(true), "Print information on the standard error output")
      ("strand", value<string>(&sigma_strand)->default_value("+"), "Direction of sigma transcription")
      // ("Csigma", value<double>(&C_sigma)->default_value(1), "Concentration of sigma70")
      // ("ClexA",  value<double>(&C_lexA)->default_value(1), "Concentration of lexA")
      ("TSS",    value<int>(&tss)->default_value(0), "Position of the Tss")
      ("sigma_right", value<int>(&sigma_right)->default_value(0), "How much sigma 70 sites are extended to the right")
      ("sigma_left",  value<int>(&sigma_left)->default_value(0), "How much sigma 70 sites are extended to the left")
      ("lexA_right", value<int>(&lexA_right)->default_value(0), "How much lexA sites are extended to the right")
      ("lexA_left",  value<int>(&lexA_left)->default_value(0), "How much lexA sites are extended to the left")
      ("thres_sigma", value<double>(&thres_sigma)->default_value(4), "Minimal energy for keeping a sigma 70 site")
      ("thres_lexA",  value<double>(&thres_lexA)->default_value(6), "Minimal energy for keeping a lexA site");
      
    variables_map vm;
    store(parse_command_line(argc, argv, desc), vm);
    if (vm.count("help")){
      cerr << desc << "\n";
      exit(0);
    }
    notify(vm);
  }
  catch (const error &ex){
    cerr << ex.what() << '\n';
    exit(0);
  }

  if(verbose) print_cmdline();
}

void CParser::print_cmdline(void){
  cerr << endl;
  cerr << "---------------------------------------------------------" << endl;
  cerr << "PARAMETERS: " << endl;
  cerr << "Filename: " << filename << endl;
  cerr << "Identifier of sigma70: " << sigma_name << endl;
  cerr << "Identifier of lexA: " << lexA_name << endl;
  //  cerr << "Sigma concentration: " << C_sigma << endl;
  // cerr << "LexA concentration: " << C_lexA << endl;
  cerr << "Sigma transcription direction: " << sigma_strand << endl;
  cerr << "Transcription starting site: " << tss << endl;
  cerr << "Sigma binding site extensions [left, right]: [" << sigma_left << ", " << sigma_right << "]" << endl;
  cerr << "LexA binding site extensions [left, right]: [" << lexA_left << ", " << lexA_right << "]" << endl;
  cerr << "Discard all sigma binding sites whose energy is lower than: " << thres_sigma << endl;
  cerr << "Discard all LexA binding sites whose energy is lower than: " << thres_lexA << endl;
  cerr << "---------------------------------------------------------" << endl << endl;
}


void CParser::ParseFile(void){
  ifstream  file(filename);
  if(file.fail()){
    cerr << "Error opening the file" << endl;
    exit(0);
  }

  string tmp_str;
  int site_start, site_end;
  double posterior, energy;
  string alignment, promoter, WM, sequence, strand;
  
  while(file>>tmp_str){
    //Site start and end
    size_t next;
    site_start = stoi(tmp_str, &next);
    site_end = stoi(tmp_str.substr(++next));
    
    //Strand, Posterior, Promoter, WM, Sequence, Energy
    file >> strand >> posterior >> alignment >> WM >> sequence >> energy >> promoter;
    
    //Fill the arrays
    if(alignment.find(sigma_name)!=string::npos){ //If this is a sigma70
      //I don't want sigma on reverse strain or with too low E
      if(strand.compare(sigma_strand) || energy<thres_sigma) continue; 
      file_params.sites_sigma_start.push_back(site_start);
      file_params.sites_sigma_end.push_back(site_end);
      file_params.energy_sigma.push_back(energy);
    }
    else if(alignment.find(lexA_name)!=string::npos){ // If this is a LexA
      //I don't want lexA with too low energy
      if(energy<thres_lexA) continue; 
      file_params.sites_lexA_start.push_back(site_start);
      file_params.sites_lexA_end.push_back(site_end);
      file_params.energy_lexA.push_back(energy);
    }
    else{
      cerr << "ERROR: FILE READING STOPPED AT:" << endl;
      cerr << "Binding site: " << site_start << " - " << site_end << endl;
      cerr << "Strand: " << strand << endl;
      cerr << "Posterior: " << posterior << endl;
      cerr << "Alignment: " << alignment << endl;
      cerr << "WM: " << WM << endl;
      cerr << "Sequence: " << sequence << endl;
      cerr << "Energy: " << energy << endl;
      cerr << "Promoter: " << promoter << endl << endl;
      exit(1);
    }
  }

  //Find length of the sequence (motevo returns ordered sites, so last position of the seq is just the last number in
  //site_end
  file_params.seq_length = site_end; 
  
  file.close();
  return;
}


void CParser::CreateSites(vector<dynamic_bitset<>>& sites, vector<double>& weights, const int seq_length,
		 const vector<double>energy, const vector<int> sites_start, const vector<int> sites_end){
  //Reserve enough memory to store all sites.
  //This is not required, but since I know in advance the dimension of the vector
  //I can avoid a reallocation of the memory in case the number of items exceed the default minimum allocation memory
  sites.reserve(energy.size());
  weights.reserve(energy.size());

  //Populate the vector with bitsets.
  //I use vectors so I can use a non-default constructor (I need it to specify the size of each bitset,
  //i.e. the length of the seq)
  for (int i = 0; i < energy.size(); ++i)  sites.emplace_back(seq_length+1); //+1 to contain also element 0

  //In every bitset store one site
  for(int i=0; i<sites.size(); i++){
    weights.push_back(exp(energy[i])); //PUT EXPONENTIAL IN FINAL VERSION!!!
    for(int j=sites_start[i]; j<=sites_end[i]; j++){
      sites[i].set(seq_length-j); //seq_length-j because sets from right
    }
  }
}

void CParser::print_params(void){
  cerr << endl << "---------------------------------------------------------" << endl;
  cerr << "INFORMATION EXTRACTED FROM MOTEVO OUTPUT: " << endl;
  cerr << "Length of the sequence: " << file_params.seq_length << endl;
  cerr << "Binding sites sigma: " << endl;
  for(auto i = file_params.sites_sigma_start.begin(); i!=file_params.sites_sigma_start.end(); i++)
    cerr << *i << ' ';
  cerr << endl;
  for(auto i = file_params.sites_sigma_end.begin(); i!=file_params.sites_sigma_end.end(); i++)
    cerr << *i << ' ';
  cerr << endl << "Energy of sigma70 binding sites: " << endl;
  for(auto i = file_params.energy_sigma.begin(); i!=file_params.energy_sigma.end(); i++)
    cerr << *i << ' ';
  
  cerr << endl << "Binding sites lexA: " << endl;
  for(auto i = file_params.sites_lexA_start.begin(); i!=file_params.sites_lexA_start.end(); i++)
    cerr << *i << ' ';
  cerr << endl;
  for(auto i = file_params.sites_lexA_end.begin(); i!=file_params.sites_lexA_end.end(); i++)
    cerr << *i << ' ';
  cerr << endl << "Energy of lexA binding sites: " << endl;
  for(auto i = file_params.energy_lexA.begin(); i!=file_params.energy_lexA.end(); i++)
    cerr << *i << ' ';
  cerr << endl << "---------------------------------------------------------" << endl;
}

void CParser::start_parsing(void)
{ 
  ParseOptions();
  ParseFile();

  if(verbose) print_params();

  //Create the binary vector with the sites for sigma70 and LexA
  CreateSites(sites_sigma, weight_sites_sigma, file_params.seq_length,
	      file_params.energy_sigma, file_params.sites_sigma_start, file_params.sites_sigma_end);

  CreateSites(sites_lexA, weight_sites_lexA, file_params.seq_length,
	      file_params.energy_lexA, file_params.sites_lexA_start, file_params.sites_lexA_end);
  /*
  if(verbose){
    cerr << "Sites sigma70:" << endl;
    for(auto i = sites_sigma.begin(); i!=sites_sigma.end(); i++)
      cerr << *i << endl;
    
    cerr << endl << "Sites lexA:" << endl;
    for(auto i = sites_lexA.begin(); i!=sites_lexA.end(); i++)
      cerr << *i << endl;
  }
  */
  return;
}
 
