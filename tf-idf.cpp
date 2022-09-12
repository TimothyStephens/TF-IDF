#include <cstdlib>
#include <iostream>
#include <vector>
#include <cstdio>
#include <sstream>
#include <algorithm>
#include <string.h>
#include <ctime>
#include <fstream>
#include <algorithm>
#include <map>
#include <chrono>
using namespace std;

// int const bufsize    = 2000000;

bool   usePvalue = true;
double slevel = 0.5; // significent level

// Get current date/time, format is YYYY-MM-DD HH:mm:ss
// Code from https://stackoverflow.com/questions/997946/how-to-get-current-time-and-date-in-c
const std::string currentDateTime() {
    time_t now = time(0);
    struct tm tstruct;
    char buf[80];
    tstruct = *localtime(&now);
    // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
    // for more information about date/time format
    strftime(buf, sizeof(buf), "%Y-%m-%d %X", &tstruct);
    return buf;
}

template <class T>
void printVector(vector<T>& a)
{
    for(int i = 0; i < a.size(); i ++) cout << a[i] << "  ";
    cout << endl;
}

void operator += (vector<int>& a, vector<int>& b)
{
    if (a.size() != b.size()) cout << "error!" << endl;
    for(int i = 0; i < a.size(); i ++) a[i] += b[i];
}

// find the idx for max value
vector<int> maxValue2DMatrix (vector< vector<double> >& rMatrix)
{
    double maxValue = -1;
    vector<int> maxIdx (2, 0);
    for(int i = 0; i < rMatrix.size(); i ++)
        for(int j = 0; j < rMatrix[i].size(); j ++)
            if(maxValue < rMatrix[i][j])
            {
                maxValue = rMatrix[i][j];
                maxIdx[0] = i, maxIdx[1] = j;
            }

    return maxIdx;
}

// big  number of sum
double average2DMatrix (vector< vector<double> >& rMatrix)
{
    double average = 0;
    for(int i = 0; i < rMatrix.size(); i ++)
        for(int j = 0; j < rMatrix[i].size(); j ++)
            average += rMatrix[i][j];
    return average / rMatrix.size() / (rMatrix[0].size() -1);
}

// find which species a sequence belongs to
int belongingSpecies(int seqIdx, vector< vector<int> >& speciesClass)
{
    for(int i = 0; i < speciesClass.size(); i ++)
        if (std::find(speciesClass[i].begin(), speciesClass[i].end(), seqIdx + 1) != speciesClass[i].end()) return i;
    return -1;
}

/*************************************************
// a charactoristic value defined on a class of sequences
// number_of_pattern stands for the overall number of each
//   subseqence in a species.
// pattern_type_sum stands for the number of different type
//   of subsequences.
// tau is a sorted vector of frequencies for each k-mer.
***************************************************/
template <typename T>
vector<int>* getTau (int spIdx, int numOfSub, T& spFreqMatrix, int& number_of_pattern )
{
    int pattern_type_sum = 0;
    vector<int>* tau = new vector<int>();
    for(int i = 0; i < numOfSub; i ++)
    {
        if(spFreqMatrix[spIdx][i] > 0)
        {
            tau->push_back(spFreqMatrix[spIdx][i]);
        }
        number_of_pattern  += spFreqMatrix[spIdx][i];
    }
    sort (tau->begin(), tau->begin());
    //pattern_type_sum  += spFreqMatrix[spIdx][i] > 0, number_of_pattern  += spFreqMatrix[spIdx][i];

    return tau;
}

double getProb(double freq, vector<int>& dist, int sum)
{
    int l_sum =0;
    for(unsigned int i =0; i < dist.size(); i ++)
        if(dist[i] < freq) l_sum += dist[i];
        else break;
    return 1.0 * l_sum / sum;

}

// recode the sequence into the a freqence vector for the accurrance of each of the subseq in a give species
template <typename T>
vector<int> subSeq2Freq(vector<int>& aSeq, T& aSpFreq)
{
    vector<int> freq4subseq(aSeq.size(), 0);
    for(int i =0; i < aSeq.size(); i ++)
        freq4subseq[i] = aSpFreq[aSeq[i] -1];
    return freq4subseq;
}

/******************************************************************************
// delta is defined on the substring of the sequences.
// Substrings are splitted from the frequence sequence from the correlated speices,
//   if the continuing zeros exceeding the number of the size of a predefined windowsize * 2
//    a string is splitted into 2.
// Delta is the sum of the frequences in the belonging spieces for a substring.
//
// Input:
//   aSeq is the sequence with a significant correlation with a species.
//   correlatedSpFreq  is the frequence vector defined on the correlated specie.
//   belongingSpFreq  is the frequence vector defined on the belonging specie.
// Output :
//   An array of three elements <startingNo, length, delta>
//
*******************************************************************************/
template <typename T>
vector<vector<int> > getDelta(vector<int>& aSeq, T& correlatedSpFreq, T& belongingSpFreq, int windowSize)
{
    vector<int> correlatedFreq4subseq = subSeq2Freq(aSeq, correlatedSpFreq);
    vector<int> belongingFreq4subseq  = subSeq2Freq(aSeq, belongingSpFreq);
    vector<vector<int> > deltas ;
    int numOfContinueZeros = 0;
    int aDelta   = 0;
    int aLengths = 0;
    bool startingZero = true;

    for(int i = 0; i < correlatedFreq4subseq.size(); i ++ )
    {
        if(correlatedFreq4subseq[i] == 0)
        {
            // ignoring start zeros
            if(startingZero) continue;
            numOfContinueZeros ++;
        }
        else
        {
            startingZero = false;

            if (numOfContinueZeros >= 2 * windowSize)
            {
                // save a delta, & the length
                vector<int> tmp(3, 0);
                tmp[0] = i - aLengths;
                tmp[1] = aLengths - numOfContinueZeros;
                tmp[2] = aDelta;
                deltas.push_back(tmp);
                // clear the value
                aDelta =0, aLengths = 0;
            }

            numOfContinueZeros = 0;
            aDelta += belongingFreq4subseq[i] ;
            //aDelta += correlatedFreq4subseq[i] ;
        }
        aLengths ++;
    }

    // deal the last or the only sequence
    vector<int> tmp(3, 0);
    tmp[0] = correlatedFreq4subseq.size() - aLengths ;
    tmp[1] = aLengths - numOfContinueZeros;
    tmp[2] = aDelta;
    deltas.push_back(tmp);
    //cout << "deltas : " << endl;
    //for(int i = 0; i < deltas.size(); i ++) printVector(deltas[i]);
    return deltas;
}

// update the rMatrix
void updating(vector< vector<int> >& rMatrix, int belongingSp, vector<int>& otherSeq, int delta)
{
    for(int i = 0; i < otherSeq.size(); i ++)
    {
        if(rMatrix[otherSeq[i] -1][belongingSp] == 0) continue;
        rMatrix[otherSeq[i] -1][belongingSp] -= delta;
    }

}

template <typename T>
void findSignificant(vector< vector<double> >& rMatrix, vector< vector<int> >& allSeq, vector< vector<int> >& speciesClass, T& spFreqMatrix, int numOfSub, int windowSize, bool usePvalue, string ofilename)
{
    cout << "[" << currentDateTime() << "] ## Step 5 of 5 - Finding significant LGT events" << endl;

    bool continueFlag = false;
    int theSeq = -1;
    int otherSpecies = -1;
    double aveFreq = average2DMatrix(rMatrix);
    cout << "[" << currentDateTime() << "] Average value of relationship matrix: " << aveFreq << endl;

    // Open results ouput file.
    fstream fout(ofilename , fstream::out  );

    int iter = 0;
    do
    {
        //cout << iter ++ << endl;

        vector<int> maxFreqIdx = maxValue2DMatrix(rMatrix);
        int theSeq = maxFreqIdx[0];
        int otherSpecies = maxFreqIdx[1];
        cout << "[" << currentDateTime() << "] "
	  << printf("Selected - species relationship: %f for sequence %d and sequence %d", rMatrix[theSeq][otherSpecies], theSeq, otherSpecies)
          << endl;

        // Check if we need to continue printing results
        continueFlag = rMatrix[theSeq][otherSpecies] > aveFreq;

        // Find the species the sequence belongs to
        int spIdx  = belongingSpecies(theSeq, speciesClass);

        if(spIdx  == -1)
        {
            cerr << "[" << currentDateTime() << "] ERROR: Could not find the belonging species for " << theSeq << ", " << spIdx << endl;
            exit(1);
        }

        // pattern_type_sum divide number_of_pattern in the species class.
        int sum =0;
        vector<int>* tau = getTau(spIdx, numOfSub, spFreqMatrix, sum);

        vector<vector<int> > deltas = getDelta(allSeq[theSeq], spFreqMatrix[otherSpecies], spFreqMatrix[spIdx], windowSize);
        for(int i = 0; i < deltas.size(); i ++)
        {
            double freq = deltas[i][2] / double(deltas[i][1]);
            vector<int>* dist = (vector<int> *)(tau);
            if(getProb(freq, *dist, sum) < slevel)  // pValue < significant level
            {
                vector<int> tmp(4, 0);
                tmp[0] = theSeq +1;     // since the sequence ID starts from 0 in the internal logics
                tmp[1] = otherSpecies +1;
                tmp[2] = deltas[i][0];
                tmp[3] = deltas[i][1];
                for(int i = 0; i < tmp.size(); i ++ )
                {   
                    fout << tmp[i] << " ";
                }
                fout << endl;
            }
        }
        rMatrix[theSeq][otherSpecies] = 0;
        delete(tau);
    }
    while (continueFlag);
    
    // Clode results file
    fout.close();
    
    cout << "[" << currentDateTime() << "] ## Step 5 of 5 - Finished finding significant LGT events!" << endl;
}

/******************************************************************************
*     Input :
*        speciesClass : containing the grouping of the sequences into speices
*            each row stands for a speices.
*        codedSeqs    : containing the coded sequences, each row of which is a coded sequence.
*        uniqueStr    : a database of <unique sequence, seqence ID>
*     logics :
*       freqMatrix  -->  spFreqMatrix --> rMatrix
*       freqMatrix is defined as the frequence vector of each sequence.
*            each sequence vector is calculated as the occurance of all possible substrings.
*       spFreqMatrix is calculated by suming up the frequence vectors of the belonging sequences.
*       rMatrix describes the sequences, species relations.
*     Output :
*       the LGT sequences information which is in the  format of <sequenceNo, correlated species, starting char position, sequence length>
*
******************************************************************************/
template <typename T>
T** createMatrix(int row, int col)
{
    cout << "[" << currentDateTime() << "] Allocating a matrix with an area of row:" << row << " col:" << col << endl;
    T ** c=NULL;
    if (( c = ( int** )malloc( row *sizeof( T* ))) == NULL )
    {
        cout << "[" << currentDateTime() << "] ERROR: Fail to allocate memory" << endl; exit(-1);
    }

    for (int  i = 0; i < row; i++ )
    {
        if (( c[i] = ( int* )malloc( col *sizeof( T ) )) == NULL )
        {
            cout << "[" << currentDateTime() << "] ERROR: Fail to allocate memory" << endl; exit(-1);
        }
    }
    cout << "[" << currentDateTime() << "]   - Done allocating matrix" << endl;
    return c;

}
int** createMatrix(int row, int col)
{
    cout << "[" << currentDateTime() << "] Allocating a matrix with an area of row:" << row << " col:" << col << endl;
    int ** c=NULL;
    if (( c = ( int** )malloc( row *sizeof( int* ))) == NULL )
    {
        cout << "[" << currentDateTime() << "] ERROR: Fail to allocate memory" << endl; exit(-1);
    }

    for (int  i = 0; i < row; i++ )
    {
        if (( c[i] = ( int* )malloc( col *sizeof( int ) )) == NULL )
        {
            cout << "[" << currentDateTime() << "] ERROR: Fail to allocate memory" << endl; exit(-1);
        }
    }
    cout << "[" << currentDateTime() << "]   - Done allocating matrix" << endl;
    return c;

}

void findLGT(vector< vector<int> > speciesClass, vector< vector<int> >& codedSeqs, map<string, int>& uniqueStr, int windowSize, string ofilename)
{
    cout << "[" << currentDateTime() << "] ## Step 4 of 5 - Finding LGT" << endl;

    // number of the substring types, number of the sequences, number of the species
    int numOfSub = 0;
    int numOfSeq = 0;
    int numOfSp = 0;
    numOfSub = uniqueStr.size();
    numOfSeq = codedSeqs.size();
    numOfSp  = speciesClass.size();

    /// build the frequence matrix.
    /// It is a 2D matrix storing the frequence of a particular subseq pattern in a sequence.
    cout << "[" << currentDateTime() << "] Building frequence matrix" << endl;
    int ** freqMatrix = createMatrix(numOfSeq, numOfSub);  // Error arises if fail to allocate more memory.

    cout << "[" << currentDateTime() << "] Populating frequence matrix" << endl;
    for(int i =0; i < numOfSeq; i ++)
    {
        cout << "[" << currentDateTime() << "]   - Loading sequence " << i << " of " << numOfSeq << endl;
        for (int j =0; j < codedSeqs[i].size(); j ++)
        {
            freqMatrix [i][codedSeqs[i][j] -1] ++;
        }
    }

    /// Build the species frequence matrix.
    /// It is a 2D matrix storing the frequence of a particular subseq pattern in a particular species.
    cout << "[" << currentDateTime() << "] Building species frequence matrix" << endl;
    int** spFreqMatrix = createMatrix(numOfSp, numOfSub);

    cout << "[" << currentDateTime() << "] Populating species frequence matrix" << endl;
    for (int i=0; i < numOfSp; i ++)
    {
        cout << "[" << currentDateTime() << "]   - Loading group " << i << " of " << numOfSp-1 << endl;
        for (int j =0; j < speciesClass[i].size(); j ++)
        {
            int seq_idx = speciesClass[i][j] -1;  //  to get the interal sequence ID
            for(int k=0; k < numOfSub; k ++)
                spFreqMatrix [i][k]  += freqMatrix[seq_idx][k];
        }
    }

    // Generate relationship matrix
    cout << "[" << currentDateTime() << "] Generate relationship matrix" << endl;
    vector< vector<double> > rMatrix(codedSeqs.size(), vector<double> (speciesClass.size(), 0));

    // Each sequence
    for(int i = 0; i < numOfSeq; i ++)
    {
        cout << "[" << currentDateTime() << "]   - Loading sequence " << i << " of " << numOfSeq-1 << endl;
        // for each subclass
        for(int j = 0; j < numOfSp; j ++)
        {
            cout << "[" << currentDateTime() << "]     - Loading group " << j << " of " << numOfSp-1 << endl;
            // avoid the calucation of seq with its own belonging species
            if( belongingSpecies(i, speciesClass) == j ) continue;
            // each subsequence pattern
            for(int m = 0; m < numOfSub; m ++)
            {
                // if a subsequence pattern does not exist in a sequence
                if (freqMatrix[i][m] == 0) continue;
                rMatrix[i][j] += spFreqMatrix[j][m] * freqMatrix[i][m];
            }
            // normalizaion origin/(speciesNO*seqLength)
            rMatrix[i][j] = rMatrix[i][j] / (speciesClass[j]).size()  / (codedSeqs[i]).size();
        }
    }

    // Print the relation Matrix
    //cout << "Relation Matrix:" << endl;
    //for(int i=0; i < rMatrix.size(); i ++)
    //{
    //    printVector(rMatrix[i]);
    //}

    // Release memory
    for(int i = 0; i < numOfSeq; i++)
    {
        free(freqMatrix[i]);
    }
    free(freqMatrix);

    cout << "[" << currentDateTime() << "] ## Step 4 of 5 - Finished finding LGT!" << endl;

    /*********************************************************
    *      detecting the significant correlations in rMatrix
    **********************************************************/
    findSignificant(rMatrix, codedSeqs, speciesClass, spFreqMatrix, numOfSub, windowSize, usePvalue, ofilename);

    // Release memory
    for(int i = 0; i < numOfSp; i++)
    {
        free(spFreqMatrix[i]);
    }
    free(spFreqMatrix);
}


int main(int argc, char* argv[])
{
    vector<string> seqs;
    time_t aTimer = std::time(NULL);

    if (argc != 4)
    {
        cout << "\n-----------------                  INFO                   -----------------" << endl;
        cout << "CMD should be in the format of:\n"  << endl;
        cout << "./tf-idf seqfile.txt k-mer_size sigLevel < speciesInfo.txt\n" << endl;
	cout << "  - k-mer_size is the k-mer size to use for analysis.\n" << endl;
        cout << "  - sigLevel is the significant level in the significant test." << endl;
        cout << "      E.g. 0.1 or 0.05.\n" << endl;
        cout << "  - speciesInfo should be in the format of :" << endl;
        cout << "      number of species in the first line. Followed by this are the " << endl;
        cout << "      the sequence order in each of the species in a single line.\n" << endl;
        cout << "      e.g.   " << endl;
        cout << "      2" << endl;
        cout << "      1,2" << endl;
        cout << "      3,4,5\n" << endl;
        cout << "      This says there are 5 sequences in the seqfile.txt. These sequences " << endl;
        cout << "      group in two species one contains 1,2, the other one includes 3,4,5" << endl;
        cout << "*******************************************************************************\n" << endl;
        exit(0);
    }
    cout << "[" << currentDateTime() << "] ---- TF-IDF starting ----" << endl;

    // Start reading arguments from command line
    ifstream readSequence(argv[1]);

    string ofilename = string(argv[1]);
    unsigned found = ofilename.find_last_of(".");
    if (found)
    {
        ofilename = ofilename.substr(0,found);
    }
    ofilename = (ofilename + ".res").c_str();

    int k = atoi(argv[2]); // k-mer size
    cout << "[" << currentDateTime() << "] K-mer: " << k << endl;

    int z = 1; // z size
    cout << "[" << currentDateTime() << "] Z size: " << z << endl;

    int windowSize = k*z;
    cout << "[" << currentDateTime() << "] Window size: " << windowSize << endl;

    double sigLevel = atof(argv[3]);
    cout << "[" << currentDateTime() << "] Significance level: " << sigLevel << endl;
    slevel = sigLevel;

    // Read sequences
    cout << "[" << currentDateTime() << "] ## Step 1 of 5 - Reading sequences from file" << endl;

    string seq;
    while (getline(readSequence, seq))
    {
        seqs.push_back(seq);
    }
    cout << "[" << currentDateTime() << "] Loaded " << seqs.size() << " sequences from file" << endl;

    cout << "[" << currentDateTime() << "] ## Step 1 of 5 - Finished reading sequences from file!" << endl;

    /*******************************************************************
    *        1. divide a sequence into substrings
    *        2. code substrings into single numbers
    *           rules: if a substring has ID value in uniqueStr
    *                use the ID value, otherwise assign an unique ID.
    *        3. uniqueStr.size() has the total number of unique substrings
    *           codedSeqs saves the coded sequences.
    *******************************************************************/
    cout << "[" << currentDateTime() << "] ## Step 2 of 5 - Generating k-mers" << endl;
    
    map<string, int> uniqueStr;
    map<string, int>::iterator iter;
    vector< vector<int> > codedSeqs;

    for (int i = 0; i < seqs.size(); i++)
    {
        cout << "[" << currentDateTime() << "]   - Processing sequence " << i << " of " << seqs.size()-1 << endl;
        vector<string> tmp;
        vector<int> tmpSeq;
        int kn = seqs[i].length() -k + 1;
	for (int j = 0; j < kn; j++)
	{
            string aSubStr = seqs[i].substr(j,k);
            int code = -1;
            iter = uniqueStr.find(aSubStr);
            if(iter == uniqueStr.end())
            {
                uniqueStr[aSubStr] = uniqueStr.size() +1;
                code = uniqueStr.size();
            }
            else
            {
                code = iter->second;
            }
            tmpSeq.push_back(code);
	}
    codedSeqs.push_back(tmpSeq);
    }
    cout << "[" << currentDateTime() << "] Number of unique k-mers: " << uniqueStr.size() << endl;

    cout << "[" << currentDateTime() << "] ## Step 2 of 5 - Finished generating k-mers!" << endl;

    /***********************************************
    *         read species information
    *     Sample Format :
    *        2  # the number of species
    *        1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16          # the ID number of sequence belonging to the 1st species
    *        17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 # the ID number of sequence belonging to the 2nd species
    *************************************************/
    cout << "[" << currentDateTime() << "] ## Step 3 of 5 - Reading species information from file" << endl;
    
    vector<vector<int> > speciesClass;
    int numOfSp = 0;
    cin >> numOfSp;
    int const MAX_CHAR_INLINE = 2000000;
    char seqIDsInStr [MAX_CHAR_INLINE];
    cin.getline(seqIDsInStr, MAX_CHAR_INLINE); // For get the newline in the first line.
    
    for (int i = 0; i < numOfSp; i ++)
    {
        char seqIDsInStr [MAX_CHAR_INLINE];
        vector<int> eachClass;
        cin.getline(seqIDsInStr, MAX_CHAR_INLINE);
        //cout << seqIDsInStr <<endl;
        stringstream ss(seqIDsInStr); // Insert the string into a stream
        string buf;
        while (ss >> buf)
            eachClass.push_back(atoi(buf.c_str()));
        cout << "[" << currentDateTime() << "] Loaded species " << i << " which contains sequences:";
        printVector(eachClass) ;
        speciesClass.push_back(eachClass);
    }

    cout << "[" << currentDateTime() << "] Loaded " << speciesClass.size() << " groups from file" << endl;

    cout << "[" << currentDateTime() << "] ## Step 3 of 5 - Finished reading species information from file!" << endl;

    /****************************************************************
    *
    *               main logic
    ****************************************************************/
    findLGT(speciesClass, codedSeqs, uniqueStr, windowSize, ofilename);

    time_t endOfTime = std::time(NULL);
    cout << "[" << currentDateTime() << "] ---- TF-IDF finished after " << endOfTime - aTimer   << " seconds ----" << endl;
}




