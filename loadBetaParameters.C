#include <fstream>
#include <iostream>
#include <sstream>
#include <map>
#include <vector>
#include <string>

std::vector<string> parseLine(string line){

  vector<string> tokenizedLine;
  string buffer;

  istringstream tokens(line);
  while(!tokens.eof())
    {
      string buffer;
      tokens >> buffer;
      tokenizedLine.push_back(buffer);
    }

  return tokenizedLine;
}

std::map<int, std::vector<double> > loadBetaParameters(std::string inputFilename ){

  std::map<int, std::vector<double> > betaParametersOut;

  ifstream file(inputFilename.c_str());
  
  string currentLine;
  int counter = 0;
  std::vector< double > fit_parameters;

  while(getline(file, currentLine)){
    // Trying to skip whitespace lines
    if (currentLine.size() > 10 && currentLine[0] != '#'){
      

      vector<string> line = parseLine(currentLine);
      
      int numberOfValues = line.size()-3;
      std::vector< double > fit_para_sect;
      //std::cout << " >> " << numberOfValues<< std::endl;
      for (int i=0; i<numberOfValues; i++){ 
	fit_parameters.push_back( atof(line[i+3].c_str()) ); 
	//std::cout << " >> " <<  atof(line[i+3].c_str()) << std::endl;
      } // Setting values	  
      
      betaParametersOut[counter] = fit_parameters;
      counter++;
      //std::cout << " COUNTER " << counter << " FIT PARAMETERS VECTOR SIZE " << fit_parameters.size() << std::endl;
      fit_parameters.clear();

    }
  }
  file.close();


  return betaParametersOut;
}


