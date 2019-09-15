// Sammie Walker
// 2315055
// swalker@chapman.edu
// CPSC350 - 03
// Assignment 1 - DNA
/*
Sources:
-http://www.cplusplus.com/doc/tutorial/files/
-https://stackoverflow.com/questions/33268513/calculating-standard-deviation-variance-in-c
-http://www.cplusplus.com/reference/string/string/npos/
-https://stackoverflow.com/questions/12225985/how-to-calculate-probability-of-a-letter-in-a-string
*/

#include <iostream>
#include <fstream> // read files in and out
#include <cmath>
#include <algorithm>
using namespace std;

class assign1{

// header functions
public:
  assign1();
  ~assign1();
  int sum(string DNAfile);
  double mean(string DNAfile);
  double variance(string DNAfile, double mean);
  double standardDeviation(double variance);
  double probability(string DNAfile, string nucleotide, int sum);
  string random(double standardDeviation, double mean, double probabilityA, double probabilityC, double probabilityG, double probabilityT);
};

// constructor and destructor
assign1::assign1(){}
assign1::~assign1(){}

// reads in a file, takes the length of each line, returning an int
// that tells the length of the DNA sequence.
int assign1::sum(string DNAfile){
  ifstream myFile; // read input from myFile
  myFile.open("assign1DNA.txt");
  if (myFile.fail()){
    cout << "File could not be opened. " << endl;
    return 0;
  }

  int sumOfDNA = 0;
  int countA = 0;
  int countC = 0;
  int countT = 0;
  int countG = 0;
  string line;

    while(true){
      getline(myFile,line);
      if (myFile.eof()) break; // if my file ends, exit

      sumOfDNA += line.size();

      for(int i = 0; i < line.size(); ++i){
        char temp  =  toupper(line[i]);
        if(temp == 'A'){
          countA++;
        }
        else if(temp == 'C'){
          countC++;
        }
        else if(temp == 'G'){
          countG++;
        }
        else{
          countT++;
        }

      int checkTotalCount = countA + countC + countG + countT; // this should be the same as sumOfDNA
      }
  }

  // have these print out to double check that it's counting correctly
  cout << "Total length of DNA sequence: " << sumOfDNA << endl;
  cout << " A COUNT: " << countA << endl;
  cout << " C COUNT: " << countC << endl;
  cout << " G COUNT: " << countG << endl;
  cout << " T COUNT: " << countT << endl;

  return sumOfDNA;

  myFile.close();
}

// reads in a file, returns double for mean
// calculates amount of lines in file and returns average length of DNA string.
double assign1::mean(string DNAfile){
  ifstream myFile; // read input from myFile
  myFile.open("assign1DNA.txt");
  if (myFile.fail()){
    cout << "File could not be opened. " << endl;
    return 0;
  }

  string line;
  double totalLines = 0;
  double sumOfDNA = 0;

  while(true){
    myFile >> line;
    if (myFile.eof()) break; // iterate through while loop until there is nothing left in file
    totalLines++;
    sumOfDNA += line.size();
    break;
  }
  return sumOfDNA/totalLines;

  cout << "MEAN: " << sumOfDNA/totalLines << endl;
  myFile.close();
}

// calculates the variance using the mean
// mean subtracted from each value, use pow to square
double assign1::variance(string DNAfile, double mean){
  ifstream myFile; // read input from myFile
  myFile.open("assign1DNA.txt");
  if (myFile.fail()){
    cout << "File could not be opened. " << endl;
    return 0;
  }
  string line;
  int totalLines = 0;
  int sumOfDNA = 0;
  double sumofValues = 0;

  while(true){
    myFile >> line;
    if(myFile.eof()) break;
      getline(myFile,line);
      totalLines++;
      sumofValues+= pow (line.size() - mean, 2.0); // pow = get power of given number
  }
  return sumofValues/totalLines;
  cout << "VARIANCE: " << sumofValues/totalLines;
  myFile.close();
}

// takes the variance calculates the standard deviation, which is variance squared
double assign1::standardDeviation(double variance){
  ifstream myFile; // read input from myFile
  myFile.open("assign1DNA.txt");
  if (myFile.fail()){
    cout << "File could not be opened. " << endl;
    return 0;
  }
  return sqrt(variance);
}

// reads in a file, finds the relative porbability of a each nucelotide and nucleotide bigram
double assign1::probability(string DNAfile, string nucleotide, int sum){
  ifstream myFile; // read input from myFile
  myFile.open("assign1DNA.txt");
  if (myFile.fail()){
    cout << "File could not be opened. " << endl;
    return 0;
  }

  double sumOfDNA = 0;
  string line;
  string DNAList;
  locale var;

  while(true){
    myFile >> line;
    if(myFile.eof()) break;

    for(int i = 0; i < line.size(); ++i){
      DNAList += toupper(line[i], var);
    }
  }
  myFile.close();
  int index = DNAList.find(nucleotide);
  while(index != string::npos){ // npos = "until the end of the string"; while we havevnt reached the end
    sumOfDNA++;
    index = DNAList.find(nucleotide, index+1);
  }
  return sumOfDNA/sum;

  char character_A = 'A';
  char character_C = 'C';
  char character_G = 'G';
  char character_T = 'T';

  // calculate the occurence of each character in the sequence, then the probability
  // use this to check work
  int occurrences_A = std::count(line.begin(), line.end(), character_A);
  double probabilityA = (double) occurrences_A/line.size();

  int occurrences_C = std::count(line.begin(), line.end(), character_C);
  double probabilityC = (double) occurrences_C/line.size();

  int occurrences_G = std::count(line.begin(), line.end(), character_G);
  double probabilityG = (double) occurrences_G/line.size();

  int occurrences_T = std::count(line.begin(), line.end(), character_T);
  double probabilityT = (double) occurrences_T/line.size();
  cout << "The probability of A in the sequence is: " << probabilityA << endl;
  cout << "The probability of C in the sequence is: " << probabilityC << endl;
  cout << "The probability of G in the sequence is: " << probabilityG << endl;
  cout << "The probability of T in the sequence is: " << probabilityT << endl;
  myFile.close();
}


// generates 1000 DNA sequences
// using stats calculated above, create Gaussian distribution
// use probability to determine how many times each nucleotide should show up
string assign1::random(double standardDeviation, double mean, double probabilityA, double probabilityC, double probabilityG, double probabilityT){
  ifstream myFile; // read input from myFile
  myFile.open("assign1DNA.txt");
  if (myFile.fail()){
    cout << "File could not be opened. " << endl;
    return 0;
  }

  // open outputFile so we can print results there
  fstream outputFile;
  outputFile.open("SammieWalker.out.txt");

  if (outputFile.is_open()){
    outputFile << "1000 generated DNA sequences: " << endl;
    outputFile << "\n";

    int countA = 0;
    int countC = 0;
    int countG = 0;
    int countT = 0;
    char aChar = 'A';
    char cChar = 'C';
    char tChar = 'T';
    char charG = 'G';
    string aa = "AA";
    string ac = "AC";
    string at = "AT";
    string ag = "AG";
    string ta = "TA";
    string tc = "TC";
    string tt = "TT";
    string tg = "TG";
    string ca = "CA";
    string cc = "CC";
    string ct = "CT";
    string cg = "CG";
    string ga = "GA";
    string gc = "GC";
    string gt = "GT";
    string gg = "GG";

    // steps for Gaussian calculation
    double calcA = (probabilityA * (1000 * mean)) + .5; // add .5 so that it rounds to the right number when cast to an int, cast will truncate it
    double calcC = (probabilityC * (1000 * mean)) + .5;
    double calcG = (probabilityG * (1000 * mean)) + .5;
    double calcT = (probabilityT * (1000 * mean)) + .5;

    for (int i = 0; i < 1000; ++i){

      // generate random numbers between 0-1 to use in equation
      double aVariable = ((double) rand() / (RAND_MAX + 1.0));
      double bVariable = ((double) rand() / (RAND_MAX + 1.0));

      // calculate equations
      double constant = sqrt((-2) * (log(aVariable))) * cos(2 * M_PI * bVariable);
      cout << "A: " << aVariable << "B: " << bVariable << endl;
      cout << "CONSTANT: " << constant << endl;
      double dVariable = (standardDeviation * constant) + mean;
      cout << "D VARIABLE: " << dVariable << endl;

      int equationCast;

      if (dVariable <= .5){
        equationCast = 1;
      }

      else{
        int equationCast = static_cast<int>(dVariable + .5);
      }

      string currLine = " ";
      int numOfNucleotides = 0;

      // create random nums between 1 and 4 to correlate to nucleotides to generate
      while(numOfNucleotides < equationCast){
        int nucleotideNumber = rand()%4; // random numb between 1 and 4
        if(nucleotideNumber == 0){
          if (countA < calcA){ // keep producing A's, add them to print out and add to count
            currLine += "A";
            numOfNucleotides ++;
            countA ++;
          }
        }
        else if(nucleotideNumber == 1){
          if(countC < calcC){
            currLine += "C";
            numOfNucleotides ++;
            countC ++;
          }
        }
        else if(nucleotideNumber == 2){
          if (countG < calcG){
            currLine += "G";
            numOfNucleotides ++;
            countG ++;
          }
        }
        else{
          if(countT < calcT){
            currLine += "T";
            numOfNucleotides ++;
            countT ++;
          }
        }
      break;
      }
      outputFile << "*" << currLine << endl;
  }

  }
  myFile.close();
}

// main method to run program
int main(int argc, char** argv){ // argc = number of paramaters passed, argv = variables passed

  if(argc > 1){ // if there is a paramater
    string inputString = argv[1]; // string being tested is the first variable passed

    // give user option to enter a file to test(even though we clafiry the file in the params for functions)
    cout << "Enter the file you would like to test: " << endl;
    cin >> inputString;

    bool noError = true;
    assign1 a ;
    while(noError){ // while there is no input error
      cout << "The results can be found in 'SammieWalker.out'. " << endl;
      ofstream outputFile;
      outputFile.open("SammieWalker.out.txt");

      // using an instance of the class, pass the params through the functions we created to get results
      int sumResults = a.sum(inputString);
      int meanResults = a.mean(inputString);
      double varianceResults = a.variance(inputString, meanResults);
      double standardDeviationResults = a.standardDeviation(varianceResults);

      string aString = "A";
      string cString = "C";
      string tString = "T";
      string stringG = "G";
      string aaString = "AA";
      string acString = "AC";
      string atString = "AT";
      string agString = "AG";
      string taString = "TA";
      string tcString = "TC";
      string ttString = "TT";
      string tgString = "TG";
      string caString = "CA";
      string ccString = "CC";
      string ctString = "CT";
      string cgString = "CG";
      string gaString = "GA";
      string gcString = "GC";
      string gtString = "GT";
      string ggString = "GG";

      double probabilityAResults = a.probability(inputString, aString, sumResults) * 100;
      double probabilityCRestults = a.probability(inputString,cString, sumResults) * 100;
      double probabilityGResults = a.probability(inputString, stringG, sumResults) * 100;
      double probabilityTResults = a.probability(inputString, tString, sumResults) * 100;
      double probabilityAAResults = a.probability(inputString, aaString, sumResults) * 100;
      double probabilityATResults = a.probability(inputString, atString, sumResults) * 100;
      double probabilityACResults = a.probability(inputString, acString, sumResults) * 100;
      double probabilityAGResults = a.probability(inputString, agString, sumResults) * 100;
      double probabilityTAResults = a.probability(inputString, taString, sumResults) * 100;
      double probabilityTTResults = a.probability(inputString, ttString, sumResults) * 100;
      double probabilityTCResults = a.probability(inputString, tcString, sumResults) * 100;
      double probabilityTGResults = a.probability(inputString, tgString, sumResults) * 100;
      double probabilityCAResults = a.probability(inputString, caString, sumResults) * 100;
      double probabilityCTResults = a.probability(inputString, ctString, sumResults) * 100;
      double probabilityCCResults = a.probability(inputString, ccString, sumResults) * 100;
      double probabilityCGResults = a.probability(inputString, cgString, sumResults) * 100;
      double probabilityGAResults = a.probability(inputString, gaString, sumResults) * 100;
      double probabilityGTResults = a.probability(inputString, gtString, sumResults) * 100;
      double probabilityGCResults = a.probability(inputString, gcString, sumResults) * 100;
      double probabilityGGResults = a.probability(inputString, ggString, sumResults) * 100;

      string randomResults = a.random(standardDeviationResults, meanResults, probabilityAResults, probabilityCRestults, probabilityGResults, probabilityTResults);


      outputFile << "Sammie Walker" << endl;
      outputFile << "2315055" << endl;
      outputFile << "swalker@chapman.edu" << endl;
      outputFile << "CPSC 350 - 03" << endl;
      outputFile << "Assignment #1" << endl;

      outputFile << "\n";

      outputFile << "SUM: " << sumResults << endl;
      outputFile << "MEAN: " << meanResults << endl;
      outputFile << "VARIANCE: " << varianceResults << endl;
      outputFile << "STANDARD DEVIATION: " << standardDeviationResults << endl;

      outputFile << "\n";

      outputFile << "Probabilities for nucleotides: " << endl;
      outputFile << aString << ": " << probabilityAResults << "%" << endl;
      outputFile << cString << ": " << probabilityCRestults << "%" << endl;
      outputFile << stringG << ": " << probabilityGResults << "%" << endl;
      outputFile << tString << ": " << probabilityTResults << "%" << endl;

      outputFile << "\n";

      outputFile << "Probabilities of nucleotide bigrams: " << endl;
      outputFile << aaString << ": " << probabilityAAResults << "%" << endl;
      outputFile << acString << ": " << probabilityACResults << "%" << endl;
      outputFile << atString << ": " << probabilityATResults << "%" << endl;
      outputFile << agString << ": " << probabilityAGResults << "%" << endl;
      outputFile << caString << ": " << probabilityCAResults << "%" << endl;
      outputFile << ccString << ": " << probabilityCCResults << "%" << endl;
      outputFile << ctString << ": " << probabilityCTResults << "%" << endl;
      outputFile << cgString << ": " << probabilityCGResults << "%" << endl;
      outputFile << gaString << ": " << probabilityGAResults << "%" << endl;
      outputFile << gcString << ": " << probabilityGCResults << "%" << endl;
      outputFile << gtString << ": " << probabilityGTResults << "%" << endl;
      outputFile << ggString << ": " << probabilityGGResults << "%" << endl;
      outputFile << taString << ": " << probabilityTAResults << "%" << endl;
      outputFile << tcString << ": " << probabilityTCResults << "%" << endl;
      outputFile << ttString << ": " << probabilityTTResults << "%" << endl;
      outputFile << tgString << ": " << probabilityTGResults << "%" << endl;

      outputFile << "\n";

      outputFile << "Here are the 1000 random DNA sequences: " << endl;

      outputFile << endl;

      outputFile << "\n";

      // ask the user if they would like to test another file
      cout << "Would you like to continue? (y/n): " << endl;
      string yesOrNo;
      cin >> yesOrNo;
      if(yesOrNo == "n"){
        cout << "End of program." << endl;
      }
      else if(yesOrNo == "y"){
        string testFile ;
        cout << "Enter the file you would liek to test: " << endl;
        cin >> testFile;

        // run the program again with the new file
        int newSumResults = a.sum(testFile);
        int newMeanResults = a.mean(testFile);
        double newVarianceResults = a.variance(testFile, meanResults);
        double newStandardDeviationResults = a.standardDeviation(varianceResults) * 100;
        double newProbabilityAResults = a.probability(testFile, aString, sumResults) * 100;
        double newProbabilityCRestults = a.probability(testFile,cString, sumResults) * 100;
        double newProbabilityGResults = a.probability(testFile, stringG, sumResults) * 100;
        double newProbabilityTResults = a.probability(testFile, tString, sumResults) * 100 ;
        double newProbabilityAAResults = a.probability(testFile, aaString, sumResults) * 100;
        double newProbabilityATResults = a.probability(testFile, atString, sumResults) * 100;
        double newProbabilityACResults = a.probability(testFile, acString, sumResults) * 100;
        double newProbabilityAGResults = a.probability(testFile, agString, sumResults) * 100;
        double newProbabilityTAResults = a.probability(testFile, taString, sumResults) * 100;
        double newProbabilityTTResults = a.probability(testFile, ttString, sumResults) * 100;
        double newProbabilityTCResults = a.probability(testFile, tcString, sumResults) * 100;
        double newProbabilityTGResults = a.probability(testFile, tgString, sumResults) * 100;
        double newProbabilityCAResults = a.probability(testFile, caString, sumResults) * 100;
        double newProbabilityCTResults = a.probability(testFile, ctString, sumResults) * 100;
        double newProbabilityCCResults = a.probability(testFile, ccString, sumResults) * 100;
        double newProbabilityCGResults = a.probability(testFile, cgString, sumResults) * 100;
        double newProbabilityGAResults = a.probability(testFile, gaString, sumResults) * 100;
        double newProbabilityGTResults = a.probability(testFile, gtString, sumResults) * 100;
        double newProbabilityGCResults = a.probability(testFile, gcString, sumResults) * 100;
        double newProbabilityGGResults = a.probability(testFile, ggString, sumResults) * 100;

        outputFile << "\n";

        outputFile << "SUM: " << sumResults << endl;
        outputFile << "MEAN: " << meanResults << endl;
        outputFile << "VARIANCE: " << varianceResults << endl;
        outputFile << "STANDARD DEVIATION: " << standardDeviationResults << endl;

        outputFile << "\n";

        outputFile << "Probabilities for nucleotides: " << endl;
        outputFile << aString << ": " << probabilityAResults << "%" << endl;
        outputFile << cString << ": " << probabilityCRestults << "%" << endl;
        outputFile << stringG << ": " << probabilityGResults << "%" << endl;
        outputFile << tString << ": " << probabilityTResults << "%" << endl;

        outputFile << "\n";

        outputFile << "Probabilities of nucleotide bigrams: " << endl;
        outputFile << aaString << ": " << probabilityAAResults << "%" << endl;
        outputFile << acString << ": " << probabilityACResults << "%" << endl;
        outputFile << atString << ": " << probabilityATResults << "%" << endl;
        outputFile << agString << ": " << probabilityAGResults << "%" << endl;
        outputFile << caString << ": " << probabilityCAResults << "%" << endl;
        outputFile << ccString << ": " << probabilityCCResults << "%" << endl;
        outputFile << ctString << ": " << probabilityCTResults << "%" << endl;
        outputFile << cgString << ": " << probabilityCGResults << "%" << endl;
        outputFile << gaString << ": " << probabilityGAResults << "%" << endl;
        outputFile << gcString << ": " << probabilityGCResults << "%" << endl;
        outputFile << gtString << ": " << probabilityGTResults << "%" << endl;
        outputFile << ggString << ": " << probabilityGGResults << "%" << endl;
        outputFile << taString << ": " << probabilityTAResults << "%" << endl;
        outputFile << tcString << ": " << probabilityTCResults << "%" << endl;
        outputFile << ttString << ": " << probabilityTTResults << "%" << endl;
        outputFile << tgString << ": " << probabilityTGResults << "%" << endl;
      }
      break;

    }
  }
}
