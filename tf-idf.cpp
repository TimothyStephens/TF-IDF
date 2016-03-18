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
using namespace std;

// int const bufsize    = 2000000;

bool   usePvalue = true;
double slevel = 0.5; // significent level

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
void* getTau (int spIdx, int numOfSub, T& spFreqMatrix, int& number_of_pattern )
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

    return (void *)(tau);
}
template <typename T>
void* getTau (int spIdx, int numOfSub, T& spFreqMatrix)
{
    int pattern_type_sum = 0, number_of_pattern = 0;
    double* res = new double();
    for(int i = 0; i < numOfSub; i ++)
    {
        pattern_type_sum  += spFreqMatrix[spIdx][i] > 0, number_of_pattern  += spFreqMatrix[spIdx][i];
    }
    *res = (   double(number_of_pattern) / pattern_type_sum);
    return (void *)(res);
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
vector< vector<int> > findSignificant(vector< vector<double> >& rMatrix, vector< vector<int> >& allSeq, vector< vector<int> >& speciesClass, T& spFreqMatrix, int numOfSub, int windowSize, bool usePvalue )
{

    bool continueFlag = false;
    int theSeq = -1,  otherSpecies = -1;
    vector< vector<int> > resultingSubstr;
    double      aveFreq = average2DMatrix(rMatrix);
		cout << "Average value of relationship matrix: " << aveFreq << endl;

    int times = 0;
    do
    {
        //cout << times ++ << endl;

        vector<int> maxFreqIdx = maxValue2DMatrix(rMatrix);

        int theSeq = maxFreqIdx[0],  otherSpecies = maxFreqIdx[1];
        continueFlag = rMatrix[theSeq][otherSpecies] > aveFreq;
	printf("select :  %f for seq: %d, group: %d. \n ", rMatrix[theSeq][otherSpecies], theSeq, otherSpecies);
        // find the species the sequence belongs to
        int spIdx  = belongingSpecies(theSeq, speciesClass);

        if(spIdx  == -1) {cerr <<  " [msg] Could not find the belonging species for " << theSeq << ", " << spIdx << endl;exit(0);}

        // pattern_type_sum divide  number_of_pattern in the species class.
        int sum =0;
        //vector<int> tau = getTau(spIdx, numOfSub, spFreqMatrix, sum);
        void* tau = NULL;
        if(usePvalue == true)
        {
            tau = getTau(spIdx, numOfSub, spFreqMatrix, sum);
        }
        else
        {
            tau = getTau(spIdx, numOfSub, spFreqMatrix);
        }

        vector<vector<int> > deltas = getDelta(allSeq[theSeq], spFreqMatrix[otherSpecies], spFreqMatrix[spIdx], windowSize);
        for(int i = 0; i < deltas.size(); i ++)
        {
            if (!usePvalue)  /// tau > delta / Length
            {
                //cout << "Using tau > delta / Length " << endl;
                double* tau1 = ((double*)(tau));
                if(*tau1 >= deltas[i][2] / double(deltas[i][1]) )
                {
                    vector<int> tmp(4, 0);
                    tmp[0] = theSeq +1;     // since the sequence ID starts from 0 in the internal logics
                    tmp[1] = otherSpecies +1;
                    tmp[2] = deltas[i][0];
                    tmp[3] = deltas[i][1];
                    resultingSubstr.push_back(tmp);
                    // update the rMatrix
                    //updating(rMatrix, spIdx, speciesClass[otherSpecies], deltas[i][2]);
                }

            }
            else  /// P-value of having a k-mer at the frequency of (delta / Length)
            {
                cout << "Using Pvalue at significant level of " << slevel << "."<< endl;
                double freq = deltas[i][2] / double(deltas[i][1]);
                vector<int>* dist = (vector<int> *)(tau);
                if(getProb(freq, *dist, sum) < slevel)  // pValue < significant level
                {
                    vector<int> tmp(4, 0);
                    tmp[0] = theSeq +1;     // since the sequence ID starts from 0 in the internal logics
                    tmp[1] = otherSpecies +1;
                    tmp[2] = deltas[i][0];
                    tmp[3] = deltas[i][1];
                    resultingSubstr.push_back(tmp);
                    // update the rMatrix
                    //updating(rMatrix, spIdx, speciesClass[otherSpecies], deltas[i][2]);
                }
            }
        }
        rMatrix[theSeq][otherSpecies] = 0;
        free(tau);


    }
    while (continueFlag);

    return resultingSubstr;
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

    T ** c=NULL;
    if (( c = ( int** )malloc( row *sizeof( T* ))) == NULL )
    {
        cout << "[msg] Fail to allocate memory" << endl; exit(-1);
    }

    for (int  i = 0; i < row; i++ )
    {
        if (( c[i] = ( int* )malloc( col *sizeof( T ) )) == NULL )
        {
            cout << "[msg] Fail to allocate memory." << endl; exit(-1);
        }
    }
    return c;

}
int** createMatrix(int row, int col)
{

    int ** c=NULL;
    if (( c = ( int** )malloc( row *sizeof( int* ))) == NULL )
    {
        cout << "[msg] Fail to allocate memory" << endl; exit(-1);
    }

    for (int  i = 0; i < row; i++ )
    {
        if (( c[i] = ( int* )malloc( col *sizeof( int ) )) == NULL )
        {
            cout << "[msg] Fail to allocate memory." << endl; exit(-1);
        }
    }
    return c;

}
vector< vector<int> > findLGT(vector< vector<int> > speciesClass, vector< vector<int> >& codedSeqs, map<string, int>& uniqueStr, int windowSize)
{
    cout << "Find LGT!" << endl;
    // number of the substring types, number of the sequences, number of the species
    int numOfSub =0, numOfSeq =0, numOfSp =0;
    numOfSub = uniqueStr.size();
    numOfSeq = codedSeqs.size();
    numOfSp  = speciesClass.size();
    /// build the frequence matrix.
    /// It is a 2D matrix storing the frequence of a particular subseq pattern in a sequence.


    int ** freqMatrix = createMatrix(numOfSeq, numOfSub);  // Error arises if fail to allocate more memory.
    cout << "sucessfully alloc" << endl;


    cout << "number of k-mer types : "<< numOfSub << endl;
    for(int i =0; i < codedSeqs.size(); i ++)
        for (int j =0; j < codedSeqs[i].size(); j ++)
            freqMatrix [i][codedSeqs[i][j] -1] ++;

    cout << "Overall k-mer frequencies in different sequences in a single species." << endl;

    /// Build the species frequence matrix.
    /// It is a 2D matrix storing the frequence of a particular subseq pattern in a particular species.
    int** spFreqMatrix = createMatrix(speciesClass.size(), numOfSub);

    for (int i=0; i < speciesClass.size(); i ++)
    {
        cout << "species " << i  << ":";
        for (int j =0; j < speciesClass[i].size(); j ++)
        {
            int seq_idx = speciesClass[i][j] -1;  //  to get the interal sequence ID


            for(int k=0; k < numOfSub; k ++)
                spFreqMatrix [i][k]  += freqMatrix[seq_idx][k];


            cout << seq_idx << ", ";
        }
        cout << endl;
        //if ( spFreqMatrix.size() <= i) cout << "error" << endl;
    }



cout << "Build the species frequence matrix" << endl;
    /// generate relationship matrix
    vector< vector<double> > rMatrix(codedSeqs.size(), vector<double> (speciesClass.size(), 0));
    //each sequence
    for(int i = 0; i < codedSeqs.size(); i ++)
    {
        // for each subclass
        for(int j = 0; j < speciesClass.size(); j ++)
        {
            // avoid the calucation of seq with its own belonging species
            if( belongingSpecies(i,speciesClass) == j) continue;
            // each subsequence pattern
            for(int m = 0; m < numOfSub; m ++)
            {
                // if a subsequence pattern does not exist in a sequence
                if (freqMatrix[i][m] == 0) continue;
                rMatrix[i][j] += spFreqMatrix[j][m] * freqMatrix[i][m];

            }
            rMatrix[i][j] = rMatrix[i][j] / (speciesClass[j]).size()  / (codedSeqs[i]).size();
        }
    }
    /// print the relation Matrix
    for(int i=0; i < rMatrix.size(); i ++)
    {
        printVector(rMatrix[i]);
    }


    free(freqMatrix);  // release memory  (important)


/*********************************************************
*      detecting the significant correlations in rMatrix
**********************************************************/

    vector< vector<int> > res=findSignificant(rMatrix, codedSeqs, speciesClass, spFreqMatrix, numOfSub, windowSize, usePvalue);

    free(spFreqMatrix); // release memory (important)
    return res;
}



vector< vector<int> > backfindLGT(vector< vector<int> > speciesClass, vector< vector<int> >& codedSeqs, map<string, int>& uniqueStr, int windowSize)
{
    // number of the substring types, number of the sequences, number of the species
    int numOfSub =0, numOfSeq =0, numOfSp =0;
    numOfSub = uniqueStr.size();
    numOfSeq = codedSeqs.size();
    numOfSp  = speciesClass.size();
    /// build the frequence matrix.
    /// It is a 2D matrix storing the frequence of a particular subseq pattern in a sequence.
    vector< vector<int> > freqMatrix(codedSeqs.size(), vector<int> (numOfSub, 0) );
    cout << "number of k-mer types : "<< numOfSub << endl;
    for(int i =0; i < codedSeqs.size(); i ++)
    {
        vector<int> aFreq(numOfSub, 0);
        for (int j =0; j < codedSeqs[i].size(); j ++)
            aFreq[codedSeqs[i][j] -1] ++;
        freqMatrix [i] = aFreq;
    }
    cout << "Overall k-mer frequencies in different sequences in a single species." << endl;
    /// Build the species frequence matrix.
    /// It is a 2D matrix storing the frequence of a particular subseq pattern in a particular species.
    vector< vector<int> > spFreqMatrix(speciesClass.size(), vector<int> (numOfSub, 0) );
    for (int i=0; i < speciesClass.size(); i ++)
    {
        cout << "species " << i  << ":";
        vector<int> aFreq(numOfSub, 0);
        for (int j =0; j < speciesClass[i].size(); j ++)
        {
            int seq_idx = speciesClass[i][j] -1;  //  to get the interal sequence ID
            aFreq  += freqMatrix[seq_idx];
            cout << seq_idx << ", ";
        }
        cout << endl;
        if ( spFreqMatrix.size() <= i) cout << "error" << endl;
        spFreqMatrix [i] = aFreq;
    }
cout << "Build the species frequence matrix" << endl;
    /// generate relationship matrix
    vector< vector<double> > rMatrix(codedSeqs.size(), vector<double> (speciesClass.size(), 0));
    //each sequence
    for(int i = 0; i < codedSeqs.size(); i ++)
    {
        // for each subclass
        for(int j = 0; j < speciesClass.size(); j ++)
        {
            // avoid the calucation of seq with its own belonging species
            if( belongingSpecies(i,speciesClass) == j) continue;
            // each subsequence pattern
            for(int m = 0; m < numOfSub; m ++)
            {
                // if a subsequence pattern does not exist in a sequence
                if (freqMatrix[i][m] == 0) continue;
                rMatrix[i][j] += spFreqMatrix[j][m] * freqMatrix[i][m];

            }
        }
    }
    /// print the relation Matrix
    for(int i=0; i < rMatrix.size(); i ++)
    {
        printVector(rMatrix[i]);
    }
/*********************************************************
*      detecting the significant correlations in rMatrix
**********************************************************/
    return findSignificant(rMatrix, codedSeqs, speciesClass, spFreqMatrix, numOfSub, windowSize, usePvalue);
}

int main(int argc, char* argv[])
{
    vector<string> seqs;
time_t aTimer = std::time(NULL);

    if (argc < 3)
    {
        cout << "\n------------                  INFO                   ------------" << endl;
        cout << "CMD should be in the format of : main.exe seqfile.txt [sigLevel] < speciesInfo.txt" << endl;
        cout << "  speciesInfo should be in the format of :" << endl;
        cout << "    number of species in the first line. Followed by this are the " << endl;
        cout << "    the sequence order in each of the species in a single line.\n" << endl;
        cout << "    [sigLevel] is the significant level in the significant test." << endl;
        cout << "    E.g. 0.1 or 0.05. " << endl;
        cout << "    when [sigLevel] =0, the default average method is performed." << endl;
        cout << "e.g.   " << endl;
        cout << "2\n1,2\n3,4,5\n" << endl;
        cout << " This says there are 5 sequences in the seqfile.txt. These sequences " << endl;
        cout << "group in two species one contains 1,2, the other one includes 3,4,5" << endl;
        cout << "*********************************************************************" << endl;

        exit(0);
    }



	ifstream readSequence(argv[1]);
	double sigLevel = atof(argv[2]);

	if (sigLevel > 1e-20)
	{
	    cout << "Will be using significant test." << endl;

        usePvalue = true;
        slevel    = sigLevel;
	}
	else
	{
	    cout << "Will be using average." << endl;
        usePvalue = false;
	}

	string seq;
	int k = 40; // k-mer size
	cout << "k is " << k << endl;
	int z = 1; // z size
	int windowSize = k*z;
	cout << "z is " << z << endl;
    // reading sequences
	while (getline(readSequence, seq))
    {
        seqs.push_back(seq);
	}



/*******************************************************************
*        1. divide a sequence into substrings
*        2. code substrings into single numbers
*           rules: if a substring has ID value in uniqueStr
*                use the ID value, otherwise assign an unique ID.
*        3. uniqueStr.size() has the total number of unique substrings
*           codedSeqs saves the coded sequences.
*******************************************************************/
    map<string, int> uniqueStr;
	map<string, int>::iterator iter;
	vector< vector<int> > codedSeqs;

    for (int i = 0; i < seqs.size(); i++)
	{
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
    cout << uniqueStr.size() << endl;

/***********************************************
*         read species information
*     Sample Format :
*        2  # the number of species
*        1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16          # the ID number of sequence belonging to the 1st species
*        17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 # the ID number of sequence belonging to the 2nd species
*************************************************/
    vector<vector<int> > speciesClass;
    int numOfSp = 0;
	cin >> numOfSp;
	int const MAX_CHAR_INLINE = 100000;
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
        //printVector(eachClass) ;
        speciesClass.push_back(eachClass);

    }

    cout << "# finish reading the sequences. Main logic!" <<endl;
    time_t endOfTime = std::time(NULL);




    //exit(0);
/****************************************************************
*
*               main logic
****************************************************************/
    vector< vector<int> > resultingSubstr = findLGT(speciesClass, codedSeqs, uniqueStr, windowSize);

/**************************************************************
*
*              Output results to file
*
***************************************************************/

    string ofilename = string(argv[1]);
    unsigned found = ofilename.find_last_of(".");
    if (found)
    {
         ofilename = ofilename.substr(0,found);

    }

    fstream fout((ofilename + ".res").c_str() , fstream::out  );
    for(int i = 0; i < resultingSubstr.size(); i ++)
    {
        for(int j = 0; j < resultingSubstr[i].size(); j ++ )
        {

                fout << resultingSubstr[i][j] << " ";
        }
        fout << endl;
    }
    fout.close();

cout << endOfTime - aTimer   << "s elapsed." << endl;

}



// int multiSeq(int times)
// {
//     vector<string> seqs;
// 	ifstream readSequence("F84_LoSeqs1000_SpecSt1_SubRtA0.seq");
// 
// 	string seq;
// 	int k = 40; // k-mer size
// 	int windowSize = k;
// 
// 	while (getline(readSequence, seq))
//     {
//         seqs.push_back(seq);
// 	}
// 	for(int i = 0; i < seqs.size(); i ++)
//     {
//         for(int j = 0; j < times; j ++)
//         {
//             cout << seqs[i];
//         }
//         cout << endl;
//     }
//     return 0;
// }
// 
// 
// 
// 
// int forDebug(int argc, char *argv[])
// {
//     /// a sample input with full anotation
//     /*
// # numberOfSubstringTypes numberOfSeq NumOfSpecies
// 4 3 3
// > 1
// 1
// > 2
// 2
// > 3
// 3
// # string matrix
// 1 1 2 3 4 4
// 2 2 2 2 1 4
// 3 3 1 1
//     */
//     int windowSize = 50;
// 
//     int numOfSub =0, numOfSeq =0, numOfSp =0;
//     cin >> numOfSub  >> numOfSeq >> numOfSp;
//     vector< vector<int> > speciesClass;
// 
// 
//     char seqIDsInStr [256];
//     cin.getline(seqIDsInStr, 256);
//     seqIDsInStr [0] = 0;
//     /******************************
//     /// read species information
//     ********************************/
//     for (int i = 0; i < numOfSp; i ++)
//     {
//         char seqIDsInStr [256];
//         vector<int> eachClass;
//         cin.getline(seqIDsInStr, 256);
// 
//         stringstream ss(seqIDsInStr); // Insert the string into a stream
//         string buf;
//         while (ss >> buf)
//             eachClass.push_back(atoi(buf.c_str()));
// 
// 
//         speciesClass.push_back(eachClass);
//     }
// 
//     /********************************
//     /// sequence information
//     ********************************/
//     vector< vector<int> > allSeq;
//     for (int i = 0; i < numOfSeq; i ++)
//     {
//         int bufsize = 10000;
//         char seqIDsInStr [bufsize];
//         vector<int> eachSeq;
//         cin.getline(seqIDsInStr, bufsize);
//         stringstream ss(seqIDsInStr); // Insert the string into a stream
//         string buf;
//         while (ss >> buf)
//             eachSeq.push_back(atoi(buf.c_str()));
//         allSeq.push_back(eachSeq);
//     }
//     // read
//     //[allSeq.size()][numOfSub][1];
// 
// 
//     /// build the frequence matrix.
//     /// It is a 2D matrix storing the frequence of a particular subseq pattern in a sequence.
//     vector< vector<int> > freqMatrix(allSeq.size(), vector<int> (numOfSub, 0) );
// 
// 
//     for(int i =0; i < allSeq.size(); i ++)
//     {
//         vector<int> aFreq(numOfSub, 0);
//         for (int j =0; j < allSeq[i].size(); j ++)
//             aFreq[allSeq[i][j] -1] ++;
//         freqMatrix [i] = aFreq;
//     }
//     /* print
//     for(int i =0; i < freqMatrix.size(); i ++)
//     {
//         printVector(freqMatrix[i]);
//         cout << endl;
//     }
//     */
// 
//     /// Build the species frequence matrix.
//     /// It is a 2D matrix storing the frequence of a particular subseq pattern in a particular species.
//     vector< vector<int> > spFreqMatrix(speciesClass.size(), vector<int> (numOfSub, 0) );
//     for (int i=0; i < speciesClass.size(); i ++)
//     {
//         vector<int> aFreq(numOfSub, 0);
//         for (int j =0; j < speciesClass[i].size(); j ++)
//         {
//             int seq_idx = speciesClass[i][j] -1;
//             aFreq += freqMatrix[seq_idx];
//         }
// 
//         spFreqMatrix [i] = aFreq;
//     }
//     // print
// 
//     //for(int i =0; i < spFreqMatrix.size(); i ++)
//     //{
//     //    printVector(spFreqMatrix[i]);
//     //}
// 
// 
//     /// generate relationship matrix
// 
//     vector< vector<int> > rMatrix(allSeq.size(), vector<int> (speciesClass.size(), 0) );
//     //each sequence
//     for(int i = 0; i < allSeq.size(); i ++)
//     {
//         // for each subclass
//         for(int j = 0; j < speciesClass.size(); j ++)
//         {
//             // avoid the calucation of seq with its own belonging species
//             if( belongingSpecies(i,speciesClass) == j) continue;
//             // each subsequence pattern
//             for(int m = 0; m < numOfSub; m ++)
//             {
//                 // if a subsequence pattern does not exist in a sequence
//                 if (freqMatrix[i][m] == 0) continue;
//                 rMatrix[i][j] += spFreqMatrix[j][m] * freqMatrix[i][m];
//                 //printVector(spFreqMatrix[j]);
//                 //cout << m << ", " << j << ", ";
//                 //cout << spFreqMatrix[j][m -1] << endl;
//             }
// 
//         }
//     }
// 
//     for(int i=0; i < rMatrix.size(); i ++)
//     {
//         printVector(rMatrix[i]);
//     }
//     findSignificant(rMatrix, allSeq, speciesClass, spFreqMatrix, numOfSub, windowSize, usePvalue );
//     system("PAUSE");
//     return EXIT_SUCCESS;
// }
// 
// 
// 
