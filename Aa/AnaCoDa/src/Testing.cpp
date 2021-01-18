#include "include/Testing.h"
#include <cstring>

#ifndef STANDALONE
#include <Rcpp.h>
using namespace Rcpp;
#endif


/* testUtility (RCPP EXPOSED)
 * Arguments: None
 * Performs Unit Testing on functions within Utility.h
 * Returns 0 if successful, 1 if error found.
 * Note: This should be the only function besides those in Utility.h that uses cout and cerr.
*/
int testUtility()
{
    int error = 0;
    int globalError = 0;
    int outError = 0;
    int errError = 0;

    error = my_print("Testing my_print, no argument.\n");

    if (error)
    {
#ifndef STANDALONE
        Rcpp::Rcerr << "Error in my_print, no argument.\n";
#else
        std::cerr << "Error in my_print, no argument.\n";
#endif
        outError = 1;
        globalError = 1;
    }

    error = my_print("Testing my_print, one argument: %.\n", 0);

    if (error)
    {
#ifndef STANDALONE
        Rcpp::Rcerr << "Error in my_print, single argument.\n";
#else
        std::cerr << "Error in my_print, single argument.\n";
#endif
        outError = 1;
        globalError = 1;
    }

    error = my_print("Testing my_print, multiple arguments: %, %, %.\n", "String", 0, 0.5);

    if (error)
    {
#ifndef STANDALONE
        Rcpp::Rcerr << "Error in my_print, multiple arguments.\n";
#else
        std::cerr << "Error in my_print, multiple arguments.\n";
#endif
        outError = 1;
        globalError = 1;
    }

    if (!outError)
        my_print("\nUtility my_print --- Pass\n");

    error = my_printError("Testing my_printError, no argument.\n");

    if (error)
    {
#ifndef STANDALONE
        Rcpp::Rcerr << "Error in my_printError, no argument.\n";
#else
        std::cerr << "Error in my_printError, no argument.\n";
#endif
        errError = 1;
        globalError = 1;
    }

    error = my_printError("Testing my_printError, one argument: %.\n", 0);

    if (error)
    {
#ifndef STANDALONE
        Rcpp::Rcerr << "Error in my_printError, single argument.\n";
#else
        std::cerr << "Error in my_printError, single argument.\n";
#endif
        errError = 1;
        globalError = 1;
    }

    error = my_printError("Testing my_printError, multiple arguments: %, %, %.\n", "String", 0, 0.5);

    if (error)
    {
#ifndef STANDALONE
        Rcpp::Rcerr << "Error in my_printError, multiple arguments.\n";
#else
        std::cerr << "Error in my_printError, multiple arguments.\n";
#endif
        errError = 1;
        globalError = 1;
    }

    // Possible TODO: Create intentional errors involving number of arguments (%), error check.

    if (!errError)
        my_print("Utility my_printError --- Pass\n");

    return globalError;
}


/* testSequenceSummary (RCPP EXPOSED)
 * Arguments: None
 * Performs Unit Testing on functions within SequenceSummary.cpp
 * that are not exposed to RCPP already.
 * Returns 0 if successful, 1 if error found.
*/
int testSequenceSummary()
{
    SequenceSummary SS("ATGCTCATTCTCACTGCTGCCTCGTAG");
    std::vector <unsigned> *uVectStar;
    std::vector <unsigned> iVect;
    std::vector <unsigned> uVect;
    int error = 0;
    int globalError = 0;

    //----------------------------//
    //------ Clear Function ------//
    //----------------------------//
    SS.clear();
    for (unsigned i = 0; i < 64; i++)
    {
        if (0 != SS.getCodonCountForCodon(i))
        {
            my_printError("Problem with Sequence Summary \"clear\" function.\n Problem at codon index %.\n", i);
            error = 1;
            globalError = 1;
        }
    }
    for (unsigned i = 0; i < 22; i++)
    {
        if (0 != SS.getAACountForAA(i))
        {
            my_printError("Problem with Sequence Summary \"clear\" function.\n Problem at amino acid index %.\n", i);
            error = 1;
            globalError = 1;
        }
    }

    if (!error)
        my_print("Sequence Summary clear --- Pass\n");
    else
        error = 0; //Reset for next function.

    //--------------------------------------//
    //------ ProcessSequence Function ------//
    //--------------------------------------//
    SS.processSequence("ATGCTCATTCTCACTGCTGCCTCGTAG");

    if (1 != SS.getAACountForAA("I"))
    {
        my_printError("Problem with Sequence Summary \"processSequence\" function.\n");
        my_printError("Problem with amino acid \"I\". I is in the sequence once, but is returning %.\n",
                      SS.getAACountForAA("I"));
        error = 1;
        globalError = 1;
    }

    if (1 != SS.getAACountForAA("T"))
    {
        my_printError("Problem with Sequence Summary \"processSequence\" function.\n");
        my_printError("Problem with amino acid \"T\". T is in the sequence once, but is returning %.\n",
                      SS.getAACountForAA("T"));
        error = 1;
        globalError = 1;
    }

    std::string codon = "ATT";
    if (1 != SS.getCodonCountForCodon(codon))
    {
        my_printError("Problem with Sequence Summary \"processSequence\" function.\n");
        my_printError("Problem with codon \"ATT\". ATT is in the sequence once, but is returning %.\n",
                      SS.getCodonCountForCodon(codon));
        error = 1;
        globalError = 1;
    }

    codon = "ACT";
    if (1 != SS.getCodonCountForCodon(codon))
    {
        my_printError("Problem with Sequence Summary \"processSequence\" function.\n");
        my_printError("Problem with codon \"ACT\". ACT is in the sequence once, but is returning %.\n",
                      SS.getCodonCountForCodon(codon));
        error = 1;
        globalError = 1;
    }

    codon = "GCT";
    if (1 != SS.getCodonCountForCodon(codon))
    {
        my_printError("Problem with Sequence Summary \"processSequence\" function.\n");
        my_printError("Problem with codon \"GCT\". GCT is in the sequence once, but is returning %.\n",
                      SS.getCodonCountForCodon(codon));
        error = 1;
        globalError = 1;
    }

    codon = "GCC";
    if (1 != SS.getCodonCountForCodon(codon))
    {
        my_printError("Problem with Sequence Summary \"processSequence\" function.\n");
        my_printError("Problem with codon \"GCC\". GCC is in the sequence once, but is returning %.\n",
                      SS.getCodonCountForCodon(codon));
        error = 1;
        globalError = 1;
    }

    uVectStar = SS.getCodonPositions("CTC");
    if ((1 != uVectStar -> at(0)) && (3 != uVectStar -> at(1)))
    {
        my_printError("Codon CTC should be found at position 1 and 3(zero indexed), but is found at these locations:\n");
        for (unsigned i = 0; i < uVectStar -> size(); i++)
        {
            my_printError("%\n", uVectStar -> at(i));
        }
        error = 1;
        globalError = 1;
    }

    uVectStar= SS.getCodonPositions("ATT");
    if (2 != uVectStar -> at(0))
    {
        my_printError("Codon ATT should be found at position 2(zero indexed), but is found at these locations:\n");
        for (unsigned i = 0; i < uVectStar -> size(); i++)
        {
            my_printError("%\n", uVectStar -> at(i));
        }
        error = 1;
        globalError = 1;
    }

    if (!error)
        my_print("Sequence Summary processSequence --- Pass\n");
    else
        error = 0; //Reset for next function.

    //------------------------------------------//
    //------ complimentNucleotide Function------//
    //------------------------------------------//
    if ('T' != SequenceSummary::complimentNucleotide('A'))
    {
        my_printError("The compliment of A should be T\n");
        error = 1;
        globalError = 1;
    }

    if ('A' != SequenceSummary::complimentNucleotide('T'))
    {
        my_printError("The compliment of T should be A\n");
        error = 1;
        globalError = 1;
    }

    if ('G' != SequenceSummary::complimentNucleotide('C'))
    {
        my_printError("The compliment of C should be G\n");
        error = 1;
        globalError = 1;
    }

    if ('C' != SequenceSummary::complimentNucleotide('G'))
    {
        my_printError("The compliment of G should be C\n");
        error = 1;
        globalError = 1;
    }

    if ('C' != SequenceSummary::complimentNucleotide('Q'))
    {
        my_printError("The compliment of Q should be C\n");
        error = 1;
        globalError = 1;
    }

    if (!error)
        my_print("Sequence Summary complimentNucleotide --- Pass\n");
    else
        error = 0; //Reset for next function.

    //----------------------------------------------//
    //------ getAACountForAA(string) Function ------//
    //----------------------------------------------//
    if (1 != SS.getAACountForAA("M"))
    {
        my_printError("Error with getAACountForAA(string) for amino acid M. Should return 1, returns %.\n",
                      SS.getAACountForAA("M"));
        error = 1;
        globalError = 1;
    }

    if (2 != SS.getAACountForAA("L"))
    {
        my_printError("Error with getAACountForAA(string) for amino acid L. Should return 2, returns %.\n",
                      SS.getAACountForAA("L"));
        error = 1;
        globalError = 1;
    }

    if (1 != SS.getAACountForAA("I"))
    {
        my_printError("Error with getAACountForAA(string) for amino acid I. Should return 1, returns %.\n",
                      SS.getAACountForAA("I"));
        error = 1;
        globalError = 1;
    }

    if (1 != SS.getAACountForAA("T"))
    {
        my_printError("Error with getAACountForAA(string) for amino acid T. Should return 1, returns %.\n",
                      SS.getAACountForAA("T"));
        error = 1;
        globalError = 1;
    }

    if (2 != SS.getAACountForAA("A"))
    {
        my_printError("Error with getAACountForAA(string) for amino acid A. Should return 2, returns %\n",
                      SS.getAACountForAA("A"));
        error = 1;
        globalError = 1;
    }

    if (1 != SS.getAACountForAA("S"))
    {
        my_printError("Error with getAACountForAA(string) for amino acid S. Should return 1, returns %\n",
                      SS.getAACountForAA("S"));
        error = 1;
        globalError = 1;
    }

    if (1 != SS.getAACountForAA("X"))
    {
        my_printError("Error with getAACountForAA(string) for amino acid X. Should return 1, returns %\n",
                      SS.getAACountForAA("X"));
        error = 1;
        globalError = 1;
    }

    if (0 != SS.getAACountForAA("G"))
    {
        my_printError("Error with getAACountForAA(string) for amino acid G. Should return 0, returns %\n",
                      SS.getAACountForAA("G"));
        error = 1;
        globalError = 1;
    }

    if (!error)
        my_print("Sequence Summary getAACountForAA(string) --- Pass\n");
    else
        error = 0; //Reset for next function.

    //---------------------------------------------//
    //------ getAACountForAA(index) Function ------//
    //---------------------------------------------//
    if (1 != SS.getAACountForAA(10))
    {
        my_printError("Error with getAACountForAA(index) for amino acid M (index 10).\n Should return 1, returns %\n",
                      SS.getAACountForAA(10));
        error = 1;
        globalError = 1;
    }

    if (2 != SS.getAACountForAA(9))
    {
        my_printError("Error with getAACountForAA(index) for amino acid L (index 9).\n Should return 2, returns %\n",
                      SS.getAACountForAA(9));
        error = 1;
        globalError = 1;
    }

    if (1 != SS.getAACountForAA(7))
    {
        my_printError("Error with getAACountForAA(index) for amino acid I (index 7).\n Should return 1, returns %\n",
                      SS.getAACountForAA(7));
        error = 1;
        globalError = 1;
    }

    if (1 != SS.getAACountForAA(16))
    {
        my_printError("Error with getAACountForAA(index) for amino acid T (index 16).\n Should return 1, returns %\n",
                      SS.getAACountForAA(16));
        error = 1;
        globalError = 1;
    }

    if (2 != SS.getAACountForAA(0))
    {
        my_printError("Error with getAACountForAA(index) for amino acid A (index 0).\n Should return 2, returns %\n",
                      SS.getAACountForAA(0));
        error = 1;
        globalError = 1;
    }

    if (1 != SS.getAACountForAA(15))
    {
        my_printError("Error with getAACountForAA(index) for amino acid S (index 15).\n Should return 1, returns %\n",
                      SS.getAACountForAA(15));
        error = 1;
        globalError = 1;
    }

    if (1 != SS.getAACountForAA(21))
    {
        my_printError("Error with getAACountForAA(index) for amino acid X (index 21).\n Should return 1, returns %\n",
                      SS.getAACountForAA(21));
        error = 1;
        globalError = 1;
    }

    if (0 != SS.getAACountForAA(2))
    {
        my_printError("Error with getAACountForAA(index) for amino acid D (index 2).\n Should return 0, returns %\n",
                      SS.getAACountForAA(2));
        error = 1;
        globalError = 1;
    }

    if (!error)
        my_print("Sequence Summary getAACountForAA(index) --- Pass\n");
    else
        error = 0; //Reset for next function.

    //--------------------------------------------//
    //------ getCodonCountsForCodon(string) ------//
    //--------------------------------------------//
    codon = "ATG";
    if (1 != SS.getCodonCountForCodon(codon))
    {
        my_printError("Error with getCodonCountForCodon(string) for %.\n Should return 1, but returns %.\n",
                      codon, SS.getCodonCountForCodon(codon));
        error = 1;
        globalError = 1;
    }

    codon = "CTC";
    if (2 != SS.getCodonCountForCodon(codon))
    {
        my_printError("Error with getCodonCountForCodon(string) for %.\n Should return 2, but returns %.\n",
                      codon, SS.getCodonCountForCodon(codon));
        error = 1;
        globalError = 1;
    }

    codon = "ATT";
    if (1 != SS.getCodonCountForCodon(codon))
    {
        my_printError("Error with getCodonCountForCodon(string) for %.\n Should return 1, but returns %.\n",
                      codon, SS.getCodonCountForCodon(codon));
        error = 1;
        globalError = 1;
    }

    codon = "ACT";
    if (1 != SS.getCodonCountForCodon(codon))
    {
        my_printError("Error with getCodonCountForCodon(string) for %.\n Should return 1, but returns %.\n",
                      codon, SS.getCodonCountForCodon(codon));
        error = 1;
        globalError = 1;
    }

    codon = "GCT";
    if (1 != SS.getCodonCountForCodon(codon))
    {
        my_printError("Error with getCodonCountForCodon(string) for %.\n Should return 1, but returns %.\n",
                      codon, SS.getCodonCountForCodon(codon));
        error = 1;
        globalError = 1;
    }

    codon = "GCC";
    if (1 != SS.getCodonCountForCodon(codon))
    {
        my_printError("Error with getCodonCountForCodon(string) for %.\n Should return 1, but returns %.\n",
                      codon, SS.getCodonCountForCodon(codon));
        error = 1;
        globalError = 1;
    }

    codon = "TCG";
    if (1 != SS.getCodonCountForCodon(codon))
    {
        my_printError("Error with getCodonCountForCodon(string) for %.\n Should return 1, but returns %.\n",
                      codon, SS.getCodonCountForCodon(codon));
        error = 1;
        globalError = 1;
    }

    codon = "TAG";
    if (1 != SS.getCodonCountForCodon(codon))
    {
        my_printError("Error with getCodonCountForCodon(string) for %.\n Should return 1, but returns %.\n",
                      codon, SS.getCodonCountForCodon(codon));
        error = 1;
        globalError = 1;
    }

    codon = "AAA";
    if (0 != SS.getCodonCountForCodon(codon))
    {
        my_printError("Error with getCodonCountForCodon(string) for %.\n Should return 0, but returns %.\n",
                      codon, SS.getCodonCountForCodon(codon));
        error = 1;
        globalError = 1;
    }

    if (!error)
        my_print("Sequence Summary getCodonCountsForCodon(string) --- Pass\n");
    else
        error = 0; //Reset for next function.

    //----------------------------------------------------//
    //------ getCodonCountsForCodon(index) Function ------//
    //----------------------------------------------------//
    if (1 != SS.getCodonCountForCodon(29))
    {
        my_printError("Error with getCodonCountForCodon(index) for codon \"ATG\" (index 29).\n");
        my_printError("Should return 1, but returns %\n", SS.getCodonCountForCodon(29));
        error = 1;
        globalError = 1;
    }

    if (2 != SS.getCodonCountForCodon(24))
    {
        my_printError("Error with getCodonCountForCodon(index) for codon \"CTC\" (index 24).\n");
        my_printError("Should return 2, but returns %\n", SS.getCodonCountForCodon(24));
        error = 1;
        globalError = 1;
    }

    if (1 != SS.getCodonCountForCodon(20))
    {
        my_printError("Error with getCodonCountForCodon(index) for codon \"ATT\" (index 20).\n");
        my_printError("Should return 1, but returns %\n", SS.getCodonCountForCodon(20));
        error = 1;
        globalError = 1;
    }

    if (1 != SS.getCodonCountForCodon(51))
    {
        my_printError("Error with getCodonCountForCodon(index) for codon \"ACT\" (index 51).\n");
        my_printError("Should return 1, but returns %\n", SS.getCodonCountForCodon(51));
        error = 1;
        globalError = 1;
    }

    if (1 != SS.getCodonCountForCodon(3))
    {
        my_printError("Error with getCodonCountForCodon(index) for codon \"GCT\" (index 3).\n");
        my_printError("Should return 1, but returns %\n", SS.getCodonCountForCodon(3));
        error = 1;
        globalError = 1;
    }

    if (1 != SS.getCodonCountForCodon(1))
    {
        my_printError("Error with getCodonCountForCodon(index) for codon \"GCC\" (index 1).\n");
        my_printError("Should return 1, but returns %\n", SS.getCodonCountForCodon(1));
        error = 1;
        globalError = 1;
    }

    if (1 != SS.getCodonCountForCodon(46))
    {
        my_printError("Error with getCodonCountForCodon(index) for codon \"TCG\" (index 46).\n");
        my_printError("Should return 1, but returns %\n", SS.getCodonCountForCodon(46));
        error = 1;
        globalError = 1;
    }

    if (1 != SS.getCodonCountForCodon(62))
    {
        my_printError("Error with getCodonCountForCodon(index) for codon \"TAG\" (index 62).\n");
        my_printError("Should return 1, but returns %\n", SS.getCodonCountForCodon(62));
        error = 1;
        globalError = 1;
    }

    if (0 != SS.getCodonCountForCodon(2))
    {
        my_printError("Error with getCodonCountForCodon(index) for codon \"AAA\" (index 2).\n");
        my_printError("Should return 1, but returns %\n", SS.getCodonCountForCodon(2));
        error = 1;
        globalError = 1;
    }

    if (!error)
        my_print("Sequence Summary getCodonCountsForCodon(index) --- Pass\n");
    else
        error = 0; //Reset for next function.

    //------------------------------------------------//
    //------ getCodonPositions(string) Function ------//
    //------------------------------------------------//
    uVectStar = SS.getCodonPositions("ATG");
    if (uVectStar -> at(0) != 0 || uVectStar -> size() != 1)
    {
        my_printError("Error with getCodonPositions(string) for codon \"ATG\".\n Should return 0, but returns:\n");
        for (unsigned i = 0; i < uVectStar -> size(); i++)
        {
            my_printError("%\n", uVectStar -> at(i));
        }
        error = 1;
        globalError = 1;
    }

    uVectStar = SS.getCodonPositions("CTC");
    if (uVectStar -> at(0) != 1 || uVectStar -> at(1) != 3|| uVectStar -> size() != 2)
    {
        my_printError("Error with getCodonPositions(string) for codon \"CTC\".\n Should return 1 and 3, but returns:\n");
        for (unsigned i = 0; i < uVectStar -> size(); i++)
        {
            my_printError("%\n", uVectStar -> at(i));
        }
        error = 1;
        globalError = 1;
    }

    uVectStar = SS.getCodonPositions("ATT");
    if (uVectStar -> at(0) != 2 || uVectStar -> size() != 1)
    {
        my_printError("Error with getCodonPositions(string) for codon \"ATT\".\n Should return 2, but returns:\n");
        for (unsigned i = 0; i < uVectStar -> size(); i++)
        {
            my_printError("%\n", uVectStar -> at(i));
        }
        error = 1;
        globalError = 1;
    }

    uVectStar = SS.getCodonPositions("ACT");
    if (uVectStar -> at(0) != 4 || uVectStar -> size() != 1)
    {
        my_printError("Error with getCodonPositions(string) for codon \"ACT\".\n Should return 4, but returns:\n");
        for (unsigned i = 0; i < uVectStar -> size(); i++)
        {
            my_printError("%\n", uVectStar -> at(i));
        }
        error = 1;
        globalError = 1;
    }


    uVectStar = SS.getCodonPositions("GCT");
    if (uVectStar -> at(0) != 5 || uVectStar -> size() != 1)
    {
        my_printError("Error with getCodonPositions(string) for codon \"GCT\".\n Should return 5, but returns:\n");
        for (unsigned i = 0; i < uVectStar -> size(); i++)
        {
            my_printError("%\n", uVectStar -> at(i));
        }
        error = 1;
        globalError = 1;
    }

    uVectStar = SS.getCodonPositions("GCC");
    if (uVectStar -> at(0) != 6 || uVectStar -> size() != 1)
    {
        my_printError("Error with getCodonPositions(string) for codon \"GCC\".\n Should return 6, but returns:\n");
        for (unsigned i = 0; i < uVectStar -> size(); i++)
        {
            my_printError("%\n", uVectStar -> at(i));
        }
        error = 1;
        globalError = 1;
    }

    uVectStar = SS.getCodonPositions("TCG");
    if (uVectStar -> at(0) != 7 || uVectStar -> size() != 1)
    {
        my_printError("Error with getCodonPositions(string) for codon \"TCG\".\n Should return 7, but returns:\n");
        for (unsigned i = 0; i < uVectStar -> size(); i++)
        {
            my_printError("%\n", uVectStar -> at(i));
        }
        error = 1;
        globalError = 1;
    }

    uVectStar = SS.getCodonPositions("TAG");
    if (uVectStar -> at(0) != 8 || uVectStar -> size() != 1)
    {
        my_printError("Error with getCodonPositions(string) for codon \"TAG\".\n Should return 8, but returns:\n");
        for (unsigned i = 0; i < uVectStar -> size(); i++)
        {
            my_printError("%\n", uVectStar -> at(i));
        }
        error = 1;
        globalError = 1;
    }

    uVectStar = SS.getCodonPositions("GTG");
    if (uVectStar -> size() != 0)
    {
        my_printError("Error with getCodonPositions(string) for codon \"GTG\".\n");
        my_printError("Should return an empty vector, but returns:\n");
        for (unsigned i = 0; i < uVectStar -> size(); i++)
        {
            my_printError("%\n", uVectStar -> at(i));
        }
        error = 1;
        globalError = 1;
    }

    if (!error)
        my_print("Sequence Summary getCodonPositions(string) --- Pass\n");
    else
        error = 0; //Reset for next function.

    //-----------------------------------------------//
    //------ getCodonPositions(index) Function ------//
    //-----------------------------------------------//
    uVectStar = SS.getCodonPositions(29);
    if (uVectStar -> at(0) != 0 || uVectStar -> size() != 1)
    {
        my_printError("Error with getCodonPositions(index) for codon index 29.\n Should return 0, but returns:\n");
        for (unsigned i = 0; i < uVectStar -> size(); i++)
        {
            my_printError("%\n", uVectStar -> at(i));
        }
        error = 1;
        globalError = 1;
    }

    uVectStar = SS.getCodonPositions(24);
    if (uVectStar -> at(0) != 1 || uVectStar -> at(1) != 3 || uVectStar -> size() != 2)
    {
        my_printError("Error with getCodonPositions(index) for codon index 24.\n Should return 1 and 3, but returns:\n");
        for (unsigned i = 0; i < uVectStar -> size(); i++)
        {
            my_printError("%\n", uVectStar -> at(i));
        }
        error = 1;
        globalError = 1;
    }

    uVectStar = SS.getCodonPositions(20);
    if (uVectStar -> at(0) != 2 || uVectStar -> size() != 1)
    {
        my_printError("Error with getCodonPositions(index) for codon index 20.\n Should return 2, but returns:\n");
        for (unsigned i = 0; i < uVectStar -> size(); i++)
        {
            my_printError("%\n", uVectStar -> at(i));
        }
        error = 1;
        globalError = 1;
    }

    uVectStar = SS.getCodonPositions(51);
    if (uVectStar -> at(0) != 4 || uVectStar -> size() != 1)
    {
        my_printError("Error with getCodonPositions(index) for codon index 51.\n Should return 4, but returns:\n");
        for (unsigned i = 0; i < uVectStar -> size(); i++)
        {
            my_printError("%\n", uVectStar -> at(i));
        }
        error = 1;
        globalError = 1;
    }

    uVectStar = SS.getCodonPositions(3);
    if (uVectStar -> at(0) != 5 || uVectStar -> size() != 1)
    {
        my_printError("Error with getCodonPositions(index) for codon index 3.\n Should return 4, but returns:\n");
        for (unsigned i = 0; i < uVectStar -> size(); i++)
        {
            my_printError("%\n", uVectStar -> at(i));
        }
        error = 1;
        globalError = 1;
    }

    uVectStar = SS.getCodonPositions(1);
    if (uVectStar -> at(0) != 6 || uVectStar -> size() != 1)
    {
        my_printError("Error with getCodonPositions(index) for codon index 1.\n Should return 4, but returns:\n");
        for (unsigned i = 0; i < uVectStar -> size(); i++)
        {
            my_printError("%\n", uVectStar -> at(i));
        }
        error = 1;
        globalError = 1;
    }

    uVectStar = SS.getCodonPositions(46);
    if (uVectStar -> at(0) != 7 || uVectStar -> size() != 1)
    {
        my_printError("Error with getCodonPositions(index) for codon index 46.\n Should return 7, but returns:\n");
        for (unsigned i = 0; i < uVectStar -> size(); i++)
        {
            my_printError("%\n", uVectStar -> at(i));
        }
        error = 1;
        globalError = 1;
    }

    uVectStar = SS.getCodonPositions(62);
    if (uVectStar -> at(0) != 8 || uVectStar -> size() != 1)
    {
        my_printError("Error with getCodonPositions(index) for codon index 62.\n Should return 8, but returns:\n");
        for (unsigned i = 0; i < uVectStar -> size(); i++)
        {
            my_printError("%\n", uVectStar -> at(i));
        }
        error = 1;
        globalError = 1;
    }

    uVectStar = SS.getCodonPositions(54);
    if (uVectStar -> size() != 0)
    {
        my_printError("Error with getCodonPositions(index) for codon index 54.\n");
        my_printError("Should return an empty vector, but returns:\n");
        for (unsigned i = 0; i < uVectStar -> size(); i++)
        {
            my_printError("%\n", uVectStar -> at(i));
        }
        error = 1;
        globalError = 1;
    }

    if (!error)
        my_print("Sequence Summary getCodonPositions(index) --- Pass\n");
    else
        error = 0; //Reset for next function.

    //----------------------------------------------------------------//
    //------ init/get/setRFPCount / getSingleRFPCount Functions ------//
    //----------------------------------------------------------------//
    SS.initRFPCount(1);
    iVect = SS.getRFPCount(0);

    if (0 != iVect.size())
    {
        my_printError("Error with initRFPCount or getRFPCount. Function should return an empty vector but returns:\n");
        for (unsigned i = 0; i < iVect.size(); i++)
        {
            my_printError("%\n", iVect[i]);
        }
        error = 1;
        globalError = 1;
    }

    iVect = {1, 2, 3, 4, 5};
    SS.setRFPCount(iVect, 0);

    if (SS.getRFPCount(0) != iVect)
    {
        my_printError("Error in initRFPCount, getRFPCount or setRFPCount.\n");
        my_printError("Function should return 1, 2, 3, 4, 5, but returns:\n");
        for (unsigned i = 0; i < iVect.size(); i++)
        {
            my_printError("%\n", iVect[i]);
        }
        error = 1;
        globalError = 1;
    }

    for (unsigned i = 0; i < 5; i++)
    {
        int tmp = SS.getSingleRFPCount(i, 0);
        if (iVect[i] != tmp)
        {
            my_printError("Error in initRFPCount, getSingleRFPCount or setRFPCount.\n");
            my_printError("Function should return %, but returns %\n", iVect[i], tmp);
            error = 1;
            globalError = 1;
        }
    }

    if (!error)
        my_print("Sequence Summary init/get/setRFPCount / getSingleRFPCount --- Pass\n");
    else
        error = 0; //Reset for next function.

    //-----------------------------------------------//
    //------ init/get/setSumRFPCount Functions ------//
    //-----------------------------------------------//
    SS.initSumRFPCount(1);
    std::array <unsigned, 64> u64Array = SS.getSumRFPCount(0);
    std::array <unsigned, 64> empty;
    empty.fill(0);

    if (u64Array != empty)
    {
        my_printError("Error with initSumRFPCount or getSumRFPCount.\n");
        my_printError("Function should return an array filled with zeroes but returns:\n");
        for (unsigned i = 0; i < u64Array.size(); i++)
        {
            if (u64Array[i] != 0) my_printError("For index %, value %\n", i, u64Array[i]);
        }
        error = 1;
        globalError = 1;
    }

    u64Array[4] = 5;
    u64Array[16] = 9;
    u64Array[24] = 6;
    u64Array[18] = 2;
    u64Array[47] = 0;

    SS.setSumRFPCount(u64Array, 0);

    if (SS.getSumRFPCount(0) != u64Array)
    {
        my_printError("Error in initSumRFPCount, getSumRFPCount or setSumRFPCount\n.");
        my_printError("Function should return the index, value pairs: 4, 5; 16, 9; 24, 6; 18, 2; 47, 0 but returns:\n");
        my_printError("4, %\n", u64Array[4]);
        my_printError("16, %\n", u64Array[16]);
        my_printError("24, %\n", u64Array[24]);
        my_printError("18, %\n", u64Array[18]);
        my_printError("27, %\n", u64Array[47]);
        error = 1;
        globalError = 1;
    }

    if (!error)
        my_print("Sequence Summary init/get/setSumRFPCount --- Pass\n");
    else
        error = 0; //Reset for next function.

    //----------------------------------------------//
    //------ get/setPositionCodonID Functions ------//
    //----------------------------------------------//
    uVect = SS.getPositionCodonID();

    if (0 != uVect.size())
    {
        my_printError("Error with getPositionCodonID. Function should return an empty vector but returns:\n");
        for (unsigned i = 0; i < uVect.size(); i++)
        {
            my_printError("%\n", uVect[i]);
        }
        error = 1;
        globalError = 1;
    }

    uVect = {4, 7, 16, 32};
    SS.setPositionCodonID(uVect);

    if (SS.getPositionCodonID() != uVect)
    {
        my_printError("Error in getPositionCodonID or setPositionCodonID.\n");
        my_printError("Function should return 4, 7, 16, 32, but returns:\n");
        for (unsigned i = 0; i < uVect.size(); i++)
        {
            my_printError("%\n", uVect[i]);
        }
        error = 1;
        globalError = 1;
    }

    if (!error)
        my_print("Sequence Summary get/setPositionCodonID --- Pass\n");
    else
        error = 0; //Reset for next function.

    //---------------------------------------------------------//
    //------ getCodonSpecificSumRFPCount(string)---------------//
    //-------setCodonSpecificSumRFPCount Functions ------------//
    //---------------------------------------------------------//
    SS.setCodonSpecificSumRFPCount(4, 35);
    SS.setCodonSpecificSumRFPCount(16,45);
    SS.setCodonSpecificSumRFPCount(54,2);
    SS.setCodonSpecificSumRFPCount(45,0);

    unsigned tmp = SS.getCodonSpecificSumRFPCount("TGC");
    if (35 != tmp)
    {
        my_printError("Error in getCodonSpecificSumRFPCount(string) or setCodonSpecificSumRFPCount for codon \"TGC\".\n");
        my_printError("Should return 35, but returns %\n", tmp);
        error = 1;
        globalError = 1;
    }

    tmp = SS.getCodonSpecificSumRFPCount("CAC");
    if (45 != tmp)
    {
        my_printError("Error in getCodonSpecificSumRFPCount(string) or setCodonSpecificSumRFPCount for codon \"CAC\".\n");
        my_printError("Should return 45, but returns %\n", tmp);
        error = 1;
        globalError = 1;
    }

    tmp = SS.getCodonSpecificSumRFPCount("GTG");
    if (2 != tmp)
    {
        my_printError("Error in getCodonSpecificSumRFPCount(string) or set RFPValue for codon \"GTG\".\n");
        my_printError("Should return 2, but returns %\n", tmp);
        error = 1;
        globalError = 1;
    }

    tmp = SS.getCodonSpecificSumRFPCount("TCC");
    if (0 != tmp)
    {
        my_printError("Error in getCodonSpecificSumRFPCount(string) or setCodonSpecificSumRFPCount for codon \"TCC\".\n");
        my_printError("Should return 0, but returns %\n", tmp);
        error = 1;
        globalError = 1;
    }

    if (!error)
        my_print("Sequence Summary getCodonSpecificSumRFPCount(string) / setCodonSpecificSumRFPCount --- Pass\n");
    else
        error = 0; //Reset for next function.

    //-----------------------------------------//
    //------ getCodonSpecificSumRFPCount(index) Function ------//
    //-----------------------------------------//
    SS.setCodonSpecificSumRFPCount(0,45);
    SS.setCodonSpecificSumRFPCount(1,52);
    SS.setCodonSpecificSumRFPCount(2,63);
    SS.setCodonSpecificSumRFPCount(60,23);

    tmp = SS.getCodonSpecificSumRFPCount(0);
    if (45 != tmp)
    {
        my_printError("Error with getCodonSpecificSumRFPCount(index) for codon index 0.\n Should return 45, but returns %\n", tmp);
        error = 1;
        globalError = 1;
    }

    tmp = SS.getCodonSpecificSumRFPCount(1);
    if (52 != tmp)
    {
        my_printError("Error with getCodonSpecificSumRFPCount(index) for codon index 1.\n Should return 52, but returns %\n", tmp);
        error = 1;
        globalError = 1;
    }

    tmp = SS.getCodonSpecificSumRFPCount(2);
    if (63 != tmp)
    {
        my_printError("Error with getCodonSpecificSumRFPCount(index) for codon index 2.\n Should return 63, but returns %\n", tmp);
        error = 1;
        globalError = 1;
    }

    tmp = SS.getCodonSpecificSumRFPCount(60);
    if (23 != tmp)
    {
        my_printError("Error with getCodonSpecificSumRFPCount(index) for codon index 60.\n should return 23, but returns %\n", tmp);
        error = 1;
        globalError = 1;
    }

    if (!error)
        my_print("Sequence Summary getCodonSpecificSumRFPCount(index) --- Pass\n");
    // No need to reset error

    return globalError;
}


/* testGene (RCPP EXPOSED)
 * Arguments: None
 * Performs Unit Testing on functions within Gene.cpp
 * that are not exposed to RCPP already.
 * Returns 0 if successful, 1 if error found.
*/
int testGene()
{
    Gene testGene;
    int error = 0;
    int globalError = 0;

    //---------------------------------//
    //------ get/setId Functions ------//
    //---------------------------------//
    testGene.setId("testGene");

    if (testGene.getId() != "testGene")
    {
        my_printError("Error in testGene: setId or getId.\n");
        globalError = 1;
    }
    else
        my_print("Gene get/setId --- Pass\n");

    //------------------------------------------//
    //------ get/setDescription Functions ------//
    //------------------------------------------//
    testGene.setDescription("Test Gene for Unit Testing");

    if (testGene.getDescription() != "Test Gene for Unit Testing")
    {
        my_printError("Error in testGene: setDescription or getDescription.\n");
        globalError = 1;
    }
    else
        my_print("Gene get/setDescription --- Pass\n");

    //---------------------------------------//
    //------ get/setSequence Functions ------//
    //---------------------------------------//
    testGene.setSequence("ATGCTCATTCTCACTGCTGCCTCGTAG");

    if (testGene.getSequence() != "ATGCTCATTCTCACTGCTGCCTCGTAG")
    {
        my_printError("Error in testGene: setSequence or getSequence.\n");
        globalError = 1;
    }
    else
        my_print("Gene get/setSequence --- Pass\n");

    //---------------------------------------------//
    //------ init/get/setRFPCount Functions ------//
    //---------------------------------------------//
    testGene.initRFPCount(1);
    std::vector <unsigned> RFPCounts = testGene.getRFPCount(0);

    if (0 != RFPCounts.size())
    {
        my_printError("Error in testGene: initRFPCount or getRFPCount.\n");
        my_printError("Function should return an empty vector but returns:\n");
        for (unsigned i = 0; i < RFPCounts.size(); i++)
        {
            my_printError("%\n", RFPCounts[i]);
        }
        error = 1;
        globalError = 1;
    }

    RFPCounts = {0, 1, 1};
    testGene.setRFPCount(RFPCounts, 0);

    if (testGene.getRFPCount(0) != RFPCounts)
    {
        my_printError("Error in testGene: initRFPCount, getRFPCount or setRFPCount.\n");
        my_printError("Function should return 0, 1, 1, but returns:\n");
        for (unsigned i = 0; i < RFPCounts.size(); i++)
        {
            my_printError("%\n", RFPCounts[i]);
        }
        error = 1;
        globalError = 1;
    }

    if (!error)
        my_print("Gene init/get/setRFPCount --- Pass\n");
    else
        error = 0; //Reset for next function.

    //------------------------------------------------//
    //------ init/get/setSumRFPCount Functions ------//
    //------------------------------------------------//
    testGene.initSumRFPCount(1);
    std::array <unsigned, 64> sumRFPCounts = testGene.getSumRFPCount(0);
    std::array <unsigned, 64> empty;
    empty.fill(0);

    if (sumRFPCounts != empty)
    {
        my_printError("Error in testGene: initSumRFPCount or getSumRFPCount.\n");
        my_printError("Function should return an array filled with zeroes but returns:\n");
        for (unsigned i = 0; i < sumRFPCounts.size(); i++)
        {
            if (sumRFPCounts[i] != 0) my_printError("For index %, value %\n", i, sumRFPCounts[i]);
        }
        error = 1;
        globalError = 1;
    }

    sumRFPCounts[4] = 5;
    sumRFPCounts[16] = 9;
    sumRFPCounts[24] = 6;
    sumRFPCounts[18] = 2;
    sumRFPCounts[47] = 0;

    testGene.setSumRFPCount(sumRFPCounts, 0);

    if (testGene.getSumRFPCount(0) != sumRFPCounts)
    {
        my_printError("Error in testGene: initSumRFPCount, getSumRFPCount or setSumRFPCount\n.");
        my_printError("Function should return the index, value pairs: 4, 5; 16, 9; 24, 6; 18, 2; 47, 0 but returns:\n");
        my_printError("4, %\n", sumRFPCounts[4]);
        my_printError("16, %\n", sumRFPCounts[16]);
        my_printError("24, %\n", sumRFPCounts[24]);
        my_printError("18, %\n", sumRFPCounts[18]);
        my_printError("27, %\n", sumRFPCounts[47]);
        error = 1;
        globalError = 1;
    }

    if (!error)
        my_print("Gene init/get/setSumRFPCount --- Pass\n");

    //-----------------------------------------//
    //------ getSequenceSummary Function ------//
    //-----------------------------------------//
    SequenceSummary SS("ATGCTCATTCTCACTGCTGCCTCGTAG");
    SequenceSummary *GeneSS = testGene.getSequenceSummary();
    for (unsigned i = 0; i < 64; i++)
    {
        if (SS.getCodonCountForCodon(i) != GeneSS->getCodonCountForCodon(i))
        {
            my_printError("Error in testGene: getSequenceSummary. Codon counts are incorrect for codon %, %.\n",
                          i, SequenceSummary::codonArray[i]);
            my_printError("Should return %, but returns %\n", SS.getCodonCountForCodon(i), GeneSS->getCodonCountForCodon(i));
            error = 1;
            globalError = 1;
        }
    }

    /*
    for (unsigned i = 0; i < 64; i++)
    {
        if (SS.getCodonSpecificSumRFPCount(i) != GeneSS->getCodonSpecificSumRFPCount(i))
        {
            my_printError("Error in testGene: getSequenceSummary. RFP observed is incorrect for codon %.\n", i);
            my_printError("Should return %, but returns %\n", SS.getCodonSpecificSumRFPCount(i), GeneSS->getCodonSpecificSumRFPCount(i));
            error = 1;
            globalError = 1;
        }
    }
     */

    //This fails because this returns pointers to vectors and they need to be compared differently.
    std::vector <unsigned> *SSvec;
    std::vector <unsigned> *Gvec;
    for (unsigned i = 0; i < 64; i++)
    {
        SSvec = SS.getCodonPositions(i);
        Gvec = GeneSS->getCodonPositions(i);
        if (SSvec->size() != Gvec->size())
        {
            my_printError("Error in testGene: getSequenceSummary. Codon positions are incorrect.\n");
            my_printError("Information in compared vectors are not of equal size.\n");
            error = 1;
            globalError = 1;
        }
        else
        {
            for (unsigned j = 0; j < SSvec->size(); j++)
            {
                if (SSvec->at(j) != Gvec->at(j))
                {
                    my_printError("Error in testGene: getSequenceSummary. Codon positions are incorrect for codon %.\n", i);
                    my_printError("Should return %, but returns %\n", SSvec->at(j), Gvec->at(j));
                    error = 1;
                    globalError = 1;
                }
            }
        }
    }

    unsigned AAListSize = (unsigned)SequenceSummary::aminoAcids().size();
    for (unsigned i = 0; i < AAListSize; i++)
    {
        if (SS.getAACountForAA(i) != GeneSS->getAACountForAA(i))
        {
            my_printError("Error in testGene: getSequenceSummary. AA counts are incorrect for amino acid %.\n", i);
            my_printError("Should return %, but returns %\n", SS.getAACountForAA(i), GeneSS->getAACountForAA(i));
            error = 1;
            globalError = 1;
        }
    }
    if (!error)
        my_print("Gene getSequenceSummary --- Pass\n");
    else
        error = 0; //Reset for next function.

    //----------------------------------------------------------//
    //------ get/setObservedSynthesisRateValues Functions ------//
    //----------------------------------------------------------//
    std::vector <double> tmp;
    tmp = testGene.getObservedSynthesisRateValues();

    if (0 != tmp.size())
    {
        my_printError("Error in testGene: getObservedSynthesisRateValues.\n");
        my_printError("Function should return an empty vector but returns:\n");
        for (unsigned i = 0; i < tmp.size(); i++)
        {
            my_printError("%\n", tmp[i]);
        }
        error = 1;
        globalError = 1;
    }

    tmp = {2.34, 3.234, 0.123};
    testGene.setObservedSynthesisRateValues(tmp);

    if (testGene.getObservedSynthesisRateValues() != tmp)
    {
        my_printError("Error in testGene: getObservedSynthesisRateValues or setObservedSynthesisRateValues.\n");
        my_printError("Function should return 2.34, 3.234, 0.123, but returns:\n");
        for (unsigned i = 0; i < tmp.size(); i++)
        {
            my_printError("%\n", tmp[i]);
        }
        error = 1;
        globalError = 1;
    }

    if (!error)
        my_print("Gene get/setObservedSynthesisRateValues --- Pass\n");
    else
        error = 0; //Reset for next function.

    //--------------------------------------------------//
    //------ getNumObservedSynthesisSets Function ------//
    //--------------------------------------------------//
    if (3 != testGene.getNumObservedSynthesisSets())
    {
        my_printError("Error in testGene: getNumObservedSynthesisSets. Function should return 3, but returns %.\n",
                      testGene.getNumObservedSynthesisSets());
        globalError = 1;
    }
    else
        my_print("Gene getNumObservedSynthesisSets --- Pass\n");

    //-----------------------------------------------//
    //------ getObservedSynthesisRate Function ------//
    //-----------------------------------------------//

    // Declared above: tmp = {2.34, 3.234, 0.123}

    for (unsigned i = 0; i < 3; i++)
    {
        if (testGene.getObservedSynthesisRate(i) != tmp[i])
        {
            my_printError("Error in testGene: getObservedSynthesisRate. Function should return % at index %, but returns %.\n",
                          tmp[i], i, testGene.getObservedSynthesisRate(i));
            error = 1;
            globalError = 1;
        }
    }

    if (!error)
        my_print("Gene getObservedSynthesisRate --- Pass\n");
    else
        error = 0; //Reset for next function.

    //--------------------------------------//
    //------ getNucleotideAt Function ------//
    //--------------------------------------//
    if ('A' != testGene.getNucleotideAt(0))
    {
        my_printError("Error in testGene: getNucleotideAt. At index 0, the return value should be 'A', but is %.\n",
                      testGene.getNucleotideAt(0));
        error = 1;
        globalError = 1;
    }
    if ('T' != testGene.getNucleotideAt(1))
    {
        my_printError("Error in testGene: getNucleotideAt. At index 1, the return value should be 'T', but is %.\n",
                      testGene.getNucleotideAt(1));
        error = 1;
        globalError = 1;
    }
    if ('G' != testGene.getNucleotideAt(2))
    {
        my_printError("Error in testGene: getNucleotideAt. At index 2, the return value should be 'G', but is %.\n",
                      testGene.getNucleotideAt(2));
        error = 1;
        globalError = 1;
    }
    if ('C' != testGene.getNucleotideAt(3))
    {
        my_printError("Error in testGene: getNucleotideAt. At index 3, the return value should be 'C', but is %.\n",
                      testGene.getNucleotideAt(3));
        error = 1;
        globalError = 1;
    }
    if ('T' != testGene.getNucleotideAt(10))
    {
        my_printError("Error in testGene: getNucleotideAt. At index 10, the return value should be 'T', but is %.\n",
                      testGene.getNucleotideAt(10));
        error = 1;
        globalError = 1;
    }
    if ('G' != testGene.getNucleotideAt(23))
    {
        my_printError("Error in testGene: getNucleotideAt. At index 23, the return value should be 'G', but is %.\n",
                      testGene.getNucleotideAt(23));
        error = 1;
        globalError = 1;
    }
    if ('G' != testGene.getNucleotideAt(26))
    {
        my_printError("Error in testGene: getNucleotideAt. At index 26, the return value should be 'G', but is %.\n",
                      testGene.getNucleotideAt(26));
        error = 1;
        globalError = 1;
    }
    if (!error)
        my_print("Gene getNucleotideAt --- Pass\n");
    else
        error = 0; //Reset for next function.

    //Todo: consider out of range test case?

    //-----------------------------//
    //------ length Function ------//
    //-----------------------------//
    if (testGene.length() == strlen("ATGCTCATTCTCACTGCTGCCTCGTAG"))
        my_print("Gene length --- Pass\n");
    else
    {
        my_printError("Error in testGene: length. Should return % but returns: %.\n",
                      strlen("ATGCTCATTCTCACTGCTGCCTCGTAG"), testGene.length());
        globalError = 1;
    }

    //----------------------------------------//
    //------ reverseComplement Function ------//
    //----------------------------------------//
    Gene tmpGene;
    tmpGene = testGene.reverseComplement();
    if ("CTACGAGGCAGCAGTGAGAATGAGCAT" == tmpGene.getSequence())
        my_print("Gene reverseComplement --- Pass\n");
    else
    {
        my_printError("Error in testGene: reverseComplement. Should return \"CTACGAGGCAGCAGTGAGAATGAGCAT\" but returns: %.\n",
                      tmpGene.getSequence());
        globalError = 1;
    }

    //-----------------------------------//
    //------ toAASequence Function ------//
    //-----------------------------------//
    if ("MLILTAASX" == testGene.toAASequence())
        my_print("Gene toAASequence --- Pass\n");
    else
    {
        my_printError("Error in testGene: toAASequence. Should return \"MLILTAASX\", but returns: %\n",
                      testGene.toAASequence());
        globalError = 1;
    }

    //----------------------------//
    //------ clear Function ------//
    //----------------------------//
    testGene.clear();
    if ("" != testGene.getId())
    {
        my_printError("Error in testGene: clear. Gene Id should be blank, but is %.\n", testGene.getId());
        error = 1;
        globalError = 1;
    }
    if ("" != testGene.getDescription())
    {
        my_printError("Error in testGene: clear. Gene description should be blank, but is %.\n",
                      testGene.getDescription());
        error = 1;
        globalError = 1;
    }
    if ("" != testGene.getSequence())
    {
        my_printError("Error in testGene: clear. Gene sequence should be blank, but is %.\n", testGene.getSequence());
        error = 1;
        globalError = 1;
    }

    if (!error)
        my_print("Gene clear --- Pass\n");
    // No need to reset error

    return globalError;
}


/* testGenomePAHelper (NOT EXPOSED)
 * Arguments: The genome to add hard-coded PA-formatted genes into, a boolean signifying if genes are simulated or not
 * Creates and adds the hard-coded PA-formatted genes (simulated or not) into the genome specified by the argument
 * Used in testGenome for convenience.
 * If simulated is true, we will only count up categories for long.
 */
void testGenomePAHelper(Genome* genome, bool simulated)
{
    // All values here are derived from readRFPData.csv's hardcoded values
    // or from readSimulatedGenome.csv

    std::string genomeString1 = simulated ? "TTTTTTATTCTTGCTGGG" : "CTTGCTATTTTTTTTGGG";
    std::string genomeString2 = simulated ? "TGGTGGATTCCTGTA" : "CCTGTAATTTGGTGG";

    Gene panse1(genomeString1, "TEST001", "No description for PA(NSE) Model");
    Gene panse2(genomeString2, "TEST002", "No description for PA(NSE) Model");

    // RFPCount for TEST001: value[position] = RFPCount
    std::vector <unsigned> test1Cat1 = {0, 0, 2, 0, 1, 1};
    std::vector <unsigned> test1Cat2 = {0, 17, 0, 1, 0, 0};

    // RFPCount for TEST002: value[position] = RFPCount
    std::vector <unsigned> test2Cat1 = {1, 1, 0, 0, 1};
    std::vector <unsigned> test2Cat2 = {2, 0, 2, 3, 6};

    std::string codon;
    unsigned index1, index2, index3, index4, index5, index6, index7, index8;

    // sumRFPCount for TEST001
    std::array <unsigned, 64> sumTest1Cat1;
    std::array <unsigned, 64> sumTest1Cat2;
    sumTest1Cat1.fill(0);
    sumTest1Cat2.fill(0);

    // sumRFPCount for TEST002
    std::array <unsigned, 64> sumTest2Cat1;
    std::array <unsigned, 64> sumTest2Cat2;
    sumTest2Cat1.fill(0);
    sumTest2Cat2.fill(0);

    codon = "CTT";
    index1 = SequenceSummary::codonToIndex(codon);

    codon = "GCT";
    index2 = SequenceSummary::codonToIndex(codon);
    sumTest1Cat2[index2] = 17;

    codon = "ATT";
    index3 = SequenceSummary::codonToIndex(codon);
    sumTest1Cat1[index3] = 2;
    sumTest2Cat2[index3] = 2;

    codon = "TTT";
    index4 = SequenceSummary::codonToIndex(codon);
    sumTest1Cat1[index4] = 1;
    sumTest1Cat2[index4] = 1;

    codon = "CCT";
    index5 = SequenceSummary::codonToIndex(codon);
    sumTest2Cat1[index5] = 1;
    sumTest2Cat2[index5] = 2;

    codon = "GTA";
    index6 = SequenceSummary::codonToIndex(codon);
    sumTest2Cat1[index6] = 1;

    codon = "TGG";
    index7 = SequenceSummary::codonToIndex(codon);
    sumTest2Cat1[index7] = 1;
    sumTest2Cat2[index7] = 9;

    codon = "GGG";
    index8 = SequenceSummary::codonToIndex(codon);
    sumTest1Cat1[index8] = 1;

    if (!simulated)
    {
        panse1.geneData.setPositionCodonID({index1, index2, index3, index4, index4, index8});
        panse2.geneData.setPositionCodonID({index5, index6, index3, index7, index7});

        panse1.initRFPCount(2);
        panse2.initRFPCount(2);
        panse1.setRFPCount(test1Cat1, 0);
        panse1.setRFPCount(test1Cat2, 1);
        panse2.setRFPCount(test2Cat1, 0);
        panse2.setRFPCount(test2Cat2, 1);
    }

    if (!simulated)
    {
        panse1.initSumRFPCount(2);
        panse2.initSumRFPCount(2);

        panse1.setSumRFPCount(sumTest1Cat2, 1);
        panse2.setSumRFPCount(sumTest2Cat2, 1);
    }
    else
    {
        panse1.initSumRFPCount(1);
        panse2.initSumRFPCount(1);
    }

    panse1.setSumRFPCount(sumTest1Cat1, 0);
    panse2.setSumRFPCount(sumTest2Cat1, 0);

    genome->addGene(panse1, simulated);
    genome->addGene(panse2, simulated);

    if (!simulated)
    {
        genome->addRFPCountColumnName("long");
        genome->addRFPCountColumnName("short");
    }
}


/* testGenomeSimulatedEqualityHelper (NOT EXPOSED)
 * Arguments: Two genomes to check if they are equal to one another, with genome1 created by the simulateGenome function.
 * Compares if two genomes are equivalent under the assumption that genome1 was created by the simulateGenome function,
 * and therefore this genome lacks position-based data and therefore cannot be checked by a simple == operator.
 *
 * Thus, compared to an == statement, the following do not need to get checked: simulatedGenes (nothing is stored),
 * numGenesWithPhi (nothing for both), RFPCountColumnNames (only one RFPCountColumn is stored, and it does not need a name [yet]),
 * seq (indeterminable because position is not stored), id and description (not handled in testing) observedSynthesisRates
 * (nothing for both), codonPositions (used in FONSE), naa (not handled in testing), RFPCount (this is a position-based value,
 * as opposed to sumRFPCount), and positionCodonID (position-based).
 *
 * Used in testGenome for convenience.
 * Returns true if successful, false if error is found.
 */
bool testGenomeSimulatedPAEqualityHelper(Genome genome1, Genome genome2)
{
    std::vector <Gene> genes1, genes2;

    genes1 = genome1.getGenes(false);
    genes2 = genome2.getGenes(false);

    // Check if number of genes in genomes are equal
    if (genes1.size() != genes2.size()) return false;

    for (unsigned i = 0u; i < genes1.size(); i++)
    {
        SequenceSummary seq1 = genes1[i].geneData;
        SequenceSummary seq2 = genes2[i].geneData;

        // Check if sumRFPCount is equal
        if (seq1.getSumRFPCount() != seq2.getSumRFPCount()) return false;

        // Check if nCodons is equal for all 64 codons.
        for (unsigned j = 0u; j < 64; j++)
        {
            if (seq1.getCodonCountForCodon(j) != seq2.getCodonCountForCodon(j))
                return false;
        }
    }

    return true;
}


/* testGenome (RCPP EXPOSED)
 * Arguments: string of the name of the directory in which testing files are found for reading and writing files
 * Performs Unit Testing on functions within Genome.cpp
 * that are not exposed to RCPP already.
 * Returns 0 if successful, 1 if error found.
 */
int testGenome(std::string testFileDir)
{
    Genome genome1;
    Genome genome2;
    Gene g1("ATGGCCACTATTGGGTCTTAG", "TEST001", "TEST001 Test Gene");
    Gene g2("TGGGATTACCAA", "TEST002", "TEST002 Test Gene");
    Gene g3("TTGGAAACCACA", "TEST003", "TEST003 Test Gene");
    Gene g4("TGGGATTACCCC", "TEST004", "TEST004 Test Gene");
    Gene s1("TGGGATTACCAA", "TEST011", "TEST011 Test Gene"); //simulated gene
    int error = 0;
    int globalError = 0;

    /* Section 1:
     * Testing / Gene / Other Functions:
     * getGene, addGene, getGenes,
     * getNumGenesWithPhi, setNumGenesWithPhi, getNumGenesWithPhiForIndex,
     * getGenomeSize, getCodonCountsPerGene, get/addRFPCountColumnNames, clear
    */

    //TODO: should improper input be given (bad id/index)?

    //-----------------------------------//
    //------ get/addGene Functions ------//
    //-----------------------------------//
    genome1.addGene(g1, false);
    genome1.addGene(s1, true); //add the simulated gene s1

    Gene test = genome1.getGene("TEST001", false);
    Gene test2 = genome1.getGene(0, false);
    Gene test3 = genome1.getGene("TEST011", true);
    Gene test4 = genome1.getGene(0, true);

    if (!(test == g1 && test2 == g1)) //checking both by string and index
    {
        my_printError("Error in testGenome: addGene or getGene with genes.\n");
        error = 1;
        globalError = 1;
    }

    if (!(test3 == s1 && test4 == s1)) //checking both by string and index
    {
        my_printError("Error in testGenome: addGene or getGene with simulated genes.\n");
        error = 1;
        globalError = 1;
    }

    if (!error)
        my_print("Genome get/addGene --- Pass\n");
    else
        error = 0; //Reset for next function.

    //-------------------------------//
    //------ getGenes Function ------//
    //-------------------------------//
    std::vector<Gene> testVec;
    testVec.push_back(g1);

    if (!(testVec == genome1.getGenes(false)))
    {
        my_printError("Error in testGenome: getGenes(false).\n");
        error = 1;
        globalError = 1;
    }

    testVec.clear();
    testVec.push_back(s1);

    if (!(testVec == genome1.getGenes(true)))
    {
        my_printError("Error in testGenome: getGenes(true).\n");
        error = 1;
        globalError = 1;
    }

    if (!error)
        my_print("Genome getGenes --- Pass\n");
    else
        error = 0; //Reset for next function.

    //----------------------------------------------//
    //------ get/setNumGenesWithPhi Functions ------//
    //----------------------------------------------//
    genome1.setNumGenesWithPhi({0, 1, 2, 3});

    std::vector<unsigned> uVector = {0, 1, 2, 3};

    if (genome1.getNumGenesWithPhi() == uVector)
        my_print("Genome get/setNumGenesWithPhi --- Pass\n");
    else
    {
        my_printError("Error in testGenome: setNumGenesWithPhi or getNumGenesWithPhi.\n");
        globalError = 1;
    }

    //-------------------------------------------------//
    //------ getNumGenesWithPhiForIndex Function ------//
    //-------------------------------------------------//
    for (unsigned i = 1; i < 4; i++)
    {
        if (genome1.getNumGenesWithPhiForIndex(i) != i)
        {
            my_printError("Error in testGenome: getNumGenesWithPhiForIndex with index %. Should return %, but instead returns %.\n",
                          i, i, genome1.getNumGenesWithPhiForIndex(i));
            error = 1;
            globalError = 1;
        }
    }

    if (!error)
        my_print("Genome getNumGenesWithPhiForIndex --- Pass\n");
    else
        error = 0; //Reset for next function.

    //------------------------------------//
    //------ getGenomeSize Function ------//
    //------------------------------------//
    if (1 != genome1.getGenomeSize(false))
    {
        my_printError("Error in testGenome: getGenomeSize(false). Should return 1, but returns %.\n",
                      genome1.getGenomeSize(false));
        error = 1;
        globalError = 1;
    }

    if (1 != genome1.getGenomeSize(true))
    {
        my_printError("Error in testGenome: getGenomeSize(true). Should return 1, but returns %.\n",
                      genome1.getGenomeSize(true));
        error = 1;
        globalError = 1;
    }

    if (!error)
        my_print("Genome getGenomeSize --- Pass\n");
    else
        error = 0; //Reset for next function.

    //--------------------------------------------//
    //------ getCodonCountsPerGene Function ------//
    //--------------------------------------------//

    //reuse generic vector of unsigned integers
    uVector = {1};
    if (uVector != genome1.getCodonCountsPerGene("ATG"))
    {
        my_printError("Error in testGenome: getCodonCountsPerGene with a single gene.\n");
        error = 1;
        globalError = 1;
    }

    genome1.addGene(g2);
    genome1.addGene(g4);

    uVector = {0, 1, 1};

    if (uVector != genome1.getCodonCountsPerGene("GAT"))
    {
        my_printError("Error in testGenome: getCodonCountsPerGene with three genes.\n");
        error = 1;
        globalError = 1;
    }

    if (!error)
        my_print("Genome getCodonCountsPerGene --- Pass\n");
    else
        error = 0; //Reset for next function.

    //-------------------------------------------------//
    //------ get/addRFPCountColumnNames Function ------//
    //-------------------------------------------------//

    std::vector<std::string> sVector = genome1.getRFPCountColumnNames();

    if (sVector.size() != 0)
    {
        my_printError("Error in testGenome: getRFPCountColumnNames. Function should return an empty vector but returns:\n");
        for (unsigned i = 0; i < sVector.size(); i++)
        {
            my_printError("%\n", sVector[i]);
        }
        error = 1;
        globalError = 1;
    }

    genome1.addRFPCountColumnName("short");
    genome1.addRFPCountColumnName("long");

    sVector = genome1.getRFPCountColumnNames();
    std::vector<std::string> sVector2 = {"short", "long"};

    if (sVector != sVector2)
    {
        my_printError("Error in testGenome: getRFPCountColumnNames or addRFPCountColumnNames. Function should return 'short', 'long', but returns:\n");
        for (unsigned i = 0; i < sVector.size(); i++)
        {
            my_printError("%\n", sVector[i]);
        }
        error = 1;
        globalError = 1;
    }

    if (!error)
        my_print("Genome get/addRFPCountColumnNames --- Pass\n");
    else
        error = 0; //Reset for next function.

    //----------------------------//
    //------ Clear Function ------//
    //----------------------------//

    // Empty Genome as a control variable
    Genome empty;

    // Test adding ObservedSynthesisRateValues
    Gene clear1("TTGATCGGGCAT", "TEST005", "TEST005 Test Gene");
    clear1.setObservedSynthesisRateValues({1, 2, 3, 4});
    genome1.addGene(clear1, false);

    genome1.clear();

    if (genome1 == empty)
        my_print("Genome clear --- Pass\n");
    else
    {
        my_printError("Error in testGenome: clear. Genome is not empty.\n");
        globalError = 1;
    }

    /* Section 2:
     * Other and File I/O Functions:
     * getGenomeForGeneIndices
     * readFasta
     * readRFPData
     * readObservedPhiValues
     */

    //-----------------------------------------------//
    //------ getGenomeForGeneIndices Function ------//
    //-----------------------------------------------//

    // add more simulated and non-simulated genes
    genome1.addGene(g1, false);
    genome1.addGene(g2, false);
    genome1.addGene(g3, false);
    genome1.addGene(g4, false);

    //reuse generic vector of unsigned integers
    uVector = {0, 1, 2, 3};
    //uVector = {1, 2, 3, 4};

    if (!(genome1 == genome1.getGenomeForGeneIndices(uVector, false)))
    {
        my_printError("Error in testGenome: getGenomeForGeneIndices with genes.\n");
        error = 1;
        globalError = 1;
    }

    genome1.clear();

    Gene s2("TAGCATGATCCA", "TEST012", "TEST002 Test Gene"); //simulated gene
    Gene s3("TCATCAGGATTC", "TEST013", "TEST002 Test Gene"); //simulated gene
    Gene s4("AAACATGTCACG", "TEST014", "TEST002 Test Gene"); //simulated gene

    genome1.addGene(s1, true);
    genome1.addGene(s2, true);
    genome1.addGene(s3, true);
    genome1.addGene(s4, true);

    if (!(genome1 == genome1.getGenomeForGeneIndices(uVector, true)))
    {
        my_printError("Error in testGenome: getGenomeForGeneIndices with simulated genes.\n");
        error = 1;
        globalError = 1;
    }

    if (!error)
        my_print("Genome getGenomeForGeneIndices --- Pass\n");
    else
        error = 0; //Reset for next function.

    //--------------------------------//
    //------ readFasta Function ------//
    //--------------------------------//
    std::string file = testFileDir + "/" + "readFasta.fasta";
    genome1.readFasta(file, false);

    Gene fasta1("ATGACCGTAATTTTTTACTAG", "TEST002", "TEST002 Test Gene");
    Gene fasta2("ATGGTCTACTTTCTGACATAG", "TEST003", "TEST003 Test Gene");

    genome2.addGene(g1, false);
    genome2.addGene(fasta1, false);
    genome2.addGene(fasta2, false);

    if (genome1 == genome2)
        my_print("Genome readFasta --- Pass\n");
    else
    {
        my_printError("Error in testGenome: readFasta. Genomes are not equivalent.\n");
        globalError = 1;
    }

    //---------------------------------//
    //------ writeFasta Function ------//
    //---------------------------------//

    // Now write a genome described above in readFasta to a file, read it in again, and then compare its validity again.
    file = testFileDir + "/" + "writeFasta.fasta";
    genome1.writeFasta(file, false);
    genome2.readFasta(file, false);

    if (!(genome1 == genome2))
    {
        my_printError("Error in testGenome: writeFasta with genes. Genomes are not equivalent.\n");
        error = 1;
        globalError = 1;
    }

    // Now, re-do writing check but with simulated genes.
    genome1.clear();

    genome1.addGene(g1, true);
    genome1.addGene(fasta1, true);
    genome1.addGene(fasta2, true);

    genome1.writeFasta(file, true);

    /* Note that while these genes were originally simulated, they are printed
     * as non-simulated genes.
     * It is up to the user to know that they were simulated, but they will
     * now be read in as non-simulated genes (and Unit Testing will compare their validity as such)
     */
    genome2.readFasta(file, false);

    genome1.clear();
    genome1.addGene(g1, false);
    genome1.addGene(fasta1, false);
    genome1.addGene(fasta2, false);

    if (!(genome1 == genome2))
    {
        my_printError("Error in testGenome: writeFasta with simulated genes. Genomes are not equivalent.\n");
        error = 1;
        globalError = 1;
    }

    if (!error)
        my_print("Genome writeFasta --- Pass\n");
    else
        error = 0; //Reset for next function.

    //----------------------------------//
    //------ readRFPData Function ------//
    //----------------------------------//
    genome2.clear();

    file = testFileDir + "/" + "readRFPData.csv";
    genome1.readRFPData(file, false);

    testGenomePAHelper(&genome2, false);

    if (genome1 == genome2)
        my_print("Genome readRFPData --- Pass\n");
    else
    {
        my_printError("Error in testGenome: readRFPData. Genomes are not equivalent.\n");
        globalError = 1;
    }

    //-----------------------------------//
    //------ writeRFPData Function ------//
    //-----------------------------------//

    // Now write a genome described above in readRFPData to file2, read it in again, and then compare its validity again.

    std::string file2 = testFileDir + "/" + "writeRFPData.csv";

    genome1.writeRFPData(file2, false);
    genome2.readRFPData(file2, false);

    if (genome1 == genome2)
        my_print("Genome writeRFPData --- Pass\n");
    else
    {
        my_printError("Error in testGenome: writeRFPData. Genomes are not equivalent.\n");
        globalError = 1;
    }

    //-----------------------------------------------------//
    //------ readSimulatedGenomeFromPAModel Function ------//
    //-----------------------------------------------------//
    file = testFileDir + "/" + "readSimulatedGenome.csv";

    // First, check if the function works compared to a PA-formatted read-in genome (genome 1).
    testGenomePAHelper(&genome1, true);
    genome2.readSimulatedGenomeFromPAModel(file);

    if (!(testGenomeSimulatedPAEqualityHelper(genome1, genome2)))
    {
        my_printError("Error in testGenome: readSimulatedGenomeFromPAModel with non-simulated genes. Genomes are not equivalent.\n");
        error = 1;
        globalError = 1;
    }

    file2 = testFileDir + "/" + "writeSimulatedGenome.csv";

    // Then, check if the function works compared to an RFPData-formatted read-in genome (genome 1 again).
    genome2.writeRFPData(file2, true);
    genome1.readSimulatedGenomeFromPAModel(file2);

    if (!(testGenomeSimulatedPAEqualityHelper(genome1, genome2)))
    {
        my_printError("Error in testGenome: readSimulatedGenomeFromPAModel with simulated genes. Genomes are not equivalent.\n");
        error = 1;
        globalError = 1;
    }

    if (!error)
        my_print("Genome readSimulatedGenomeFromPAModel --- Pass\n");
    else
        error = 0; //Reset for next function.

    /* readObservedPhiValues Testing Function
     *
     * Compares a genome with the readObservedPhiValues function's created genome.
     * Reads in "readObservedPhiValues.csv" and "readObservedPhiValuesError.csv" twice each,
     * once for byID and once for byIndex.
     *
     * Significant standard error output is produced by design: both files exhibit some errors.
     */

    /* TODO NOTE: This testing function has been disabled for tester convenience.
     * It purposefully spits out a half-dozen glaring red error messages, which may
     * be confusing when checking for actual red error messages.
     */
    /*
    //--------------------------------------------//
    //------ readObservedPhiValues Function ------//
    //--------------------------------------------//
    genome1.clear();
    genome2.clear();
    std::vector <double> emptyVec; // Empty vector used to clear ObservedSynthesisRateValues

    // Test 1: Test non-error file vs by ID readObservedPhiValues function
    genome1.addGene(g1, false);
    g1.setObservedSynthesisRateValues({1, 2, 3, 4});
    genome2.addGene(g1, false);

    genome1.addGene(g2, false);
    g2.setObservedSynthesisRateValues({4, 3, 2, 1});
    genome2.addGene(g2, false);

    genome1.addGene(g3, false);
    g3.setObservedSynthesisRateValues({-1, -1, 4, 2});
    genome2.addGene(g3, false);

    genome1.addGene(g4, false);
    g4.setObservedSynthesisRateValues({2, 1, 4, -1});
    genome2.addGene(g4, false);

    genome2.setNumGenesWithPhi({3, 3, 4, 3});

    file = testFileDir + "/" + "readObservedPhiValues.csv";
    genome1.readObservedPhiValues(file, true);

    if (!(genome1 == genome2))
    {
        my_printError("Error in testGenome comparing genomes: readObservedPhiValues.csv ");
        my_printError("by ID produces a different genome than expected.\n");
        error = 1;
        globalError = 1;
    }
    genome1.clear();

    // Test 2: Test non-error file vs by index readObservedPhiValues function
    // Re-input genome as it was in the previous test, then run it by index instead
    g1.setObservedSynthesisRateValues(emptyVec);
    genome1.addGene(g1, false);
    g2.setObservedSynthesisRateValues(emptyVec);
    genome1.addGene(g2, false);
    g3.setObservedSynthesisRateValues(emptyVec);
    genome1.addGene(g3, false);
    g4.setObservedSynthesisRateValues(emptyVec);
    genome1.addGene(g4, false);

    genome1.readObservedPhiValues(file, false);

    if (!(genome1 == genome2))
    {
        my_printError("Error in testGenome comparing genomes: readObservedPhiValues.csv ");
        my_printError("by index produces a different genome than expected.\n");
        error = 1;
        globalError = 1;
    }
    genome1.clear();
    genome2.clear();

    // Test 3: Test error file vs by ID readObservedPhiValues function
    // Since this file has an error in number of phi values, the ObservedSynthesisRateValues are cleared
    genome1.addGene(g1, false);
    genome2.addGene(g1, false);

    genome1.addGene(g2, false);
    genome2.addGene(g2, false);

    genome1.addGene(g3, false);
    genome2.addGene(g3, false);

    genome1.addGene(g4, false);
    genome2.addGene(g4, false);

    // As discussed in the documentation, however, NumGenesWithPhi is still initialized with 0's despite the error
    genome2.setNumGenesWithPhi({0, 0, 0, 0});

    file = testFileDir + "/" + "readObservedPhiValuesError.csv";
    genome1.readObservedPhiValues(file, true);

    if (!(genome1 == genome2))
    {
        my_printError("Error in testGenome comparing genomes: readObservedPhiValuesError.csv ");
        my_printError("by ID produces a different genome than expected.\n");
        error = 1;
        globalError = 1;
    }

    genome1.clear();

    // Test 4: Test error file vs by index readObservedPhiValues function
    // Re-input genome as it was in the previous test, then run it by index instead
    genome1.addGene(g1, false);
    genome1.addGene(g2, false);
    genome1.addGene(g3, false);
    genome1.addGene(g4, false);

    genome1.readObservedPhiValues(file, false);

    if (!(genome1 == genome2))
    {
        my_printError("Error in testGenome comparing genomes: readObservedPhiValuesError.csv ");
        my_printError("by index produces a different genome than expected.\n");
        error = 1;
        globalError = 1;
    }

    // If any errors are produced, reset variable for next function
    if (!error)
        my_print("Genome readObservedPhiValues --- Pass\n");
    // No need to reset error
    */

    return globalError;
}


/* testParameter (RCPP EXPOSED)
 * Arguments: None
 * Performs Unit Testing on functions within Parameter.cpp
 * that are not exposed to RCPP already.
 * Returns 0 if successful, 1 if error found.
*/
int testParameter(std::string testFileDir)
{
    Parameter parameter;
    int error = 0;
    int globalError = 0;

    /* Section 1: 26 functions tested in total.
     * initParameterSet Function
     * and related get/set functions as a consequence of the function setup:
     * get/setMixtureAssignment,
     * getMutationSelectionState, getNumParam, getNumMixtureElements
     * get/setStdDevSynthesisRate, getCurrentStdDevSynthesisRateProposalWidth
     * getNumAcceptForStdDevSynthesisRate, getStdCspForIndex, getNumAcceptForCspForIndex
     * getNumMutationCategories, getNumSelectionCategories, getNumSynthesisRateCategories
     * getMutationCategory, getSelectionCategory
     * getMixtureElementsOfMutationCategory, getMixtureElementsOfSelectionCategory
     * get/setCategoryProbability
     * get/setSynthesisRate, getSynthesisRateCategory
     * getNumAcceptForSynthesisRate, getSynthesisRateProposalWidth, getCurrentSynthesisRateProposalWidth
    */

    //---------------------------------------//
    //------ initParameterSet Function ------//
    //---------------------------------------//

    /* Initialize parameter:
     * Arguments: vector <double> stdDevSynthesisRate, unsigned numMixtures, vector <unsigned> geneAssignment,
     *           vector <vector <unsigned>> mixtureDefinitionMatrix, bool splitSer, string mutationSelectionState
     *
     * Thus, let:
    */
    Genome genome;
    std::string file = testFileDir + "/" + "simulatedAllUniqueR.fasta";
    genome.readFasta(file);

    unsigned numMixtures = 3;
    std::vector<double> stdDev(numMixtures, 1);
    unsigned numGenes = genome.getGenomeSize();
    std::vector<unsigned> geneAssignment(numGenes);
    if (numMixtures == 1)
    {
        for (unsigned i = 0u; i < numGenes; i++)
        {
            geneAssignment[i] = 0u;
        }
    }
    else if (numMixtures == 3)
    {
        for (unsigned i = 0u; i < numGenes; i++)
        {
            if (i < 961) geneAssignment[i] = 0u;
            else if (i < 1418) geneAssignment[i] = 1u;
            else geneAssignment[i] = 2u;
        }
    }
    std::vector<std::vector<unsigned>> mixtureDefinitionMatrix;
    bool splitSer = true;
    std::string mutationSelectionState = Parameter::allUnique;

    parameter.initParameterSet(stdDev, numMixtures, geneAssignment, mixtureDefinitionMatrix, splitSer, mutationSelectionState);
    unsigned initParameterSetError = 0;

    /* This call changes many variables in parameter that must now be checked.
     * Thus, unit testing is done in order of variable changed.
     * See initParameterSet in Parameter.cpp or the table of contents above.
     * This also introduces a level of uncertainty in what may be wrong, and thus an error in the following
     * unit testing checks may be a result of the checking function or initParameterSet.
    */

    //------------------------------------------//
    //------ getMixtureAssignment Function------//
    //------------------------------------------//

    for (unsigned i = 0u; i < numGenes; i++)
    {
        // Note: This section of code is because vectors in R are 1-indexed (i.e. for getMixtureAssignment)
#ifndef STANDALONE
        if (parameter.getMixtureAssignment(i) != geneAssignment[i] - 1)
        {
            my_printError("Error in initParameterSet or getMixtureAssignment for index %.", i);
            my_printError(" Value should be %, but is instead %.\n", geneAssignment[i], parameter.getMixtureAssignment(i + 1));
            error = 1;
            globalError = 1;
            initParameterSetError = 1;
        }
#else
        if (parameter.getMixtureAssignment(i) != geneAssignment[i])
        {
            my_printError("Error in initParameterSet or getMixtureAssignment for index %.", i);
            my_printError(" Value should be %, but is instead %.\n", i, geneAssignment[i], parameter.getMixtureAssignment(i));
            error = 1;
            globalError = 1;
            initParameterSetError = 1;
        }
#endif
    }

    if (!error)
        my_print("Parameter getMixtureAssignment --- Pass\n");
    else
        error = 0; //Reset for next function.

    //-------------------------------------------//
    //------ setMixtureAssignment Functions------//
    //-------------------------------------------//
    for (unsigned i = 0u; i < numGenes; i++)
    {
        parameter.setMixtureAssignment(i, i);
        if (parameter.getMixtureAssignment(i) != i)
        {
            my_printError("Error in setMixtureAssignment for index %. Value should be %, but is instead %.\n",
                          i, i, parameter.getMixtureAssignment(i));
            error = 1;
            globalError = 1;
        }
    }

    if (!error)
        my_print("Parameter setMixtureAssignment --- Pass\n");
    else
        error = 0; //Reset for next function.

    //------------------------------------------------//
    //------ getMutationSelectionState Function ------//
    //------------------------------------------------//
    if (parameter.getMutationSelectionState() != mutationSelectionState)
    {
        my_printError("Error in initParameterSet or getMutationSelectionState. Value should be %, but is instead %.\n",
                      mutationSelectionState, parameter.getMutationSelectionState());
        globalError = 1;
        initParameterSetError = 1;
    }
    else
        my_print("Parameter getMutationSelectionState --- Pass\n");

    //----------------------------------//
    //------ getNumParam Function ------//
    //----------------------------------//

    unsigned int numParam = ((splitSer) ? 40 : 41);

    if (parameter.getNumParam() != numParam)
    {
        my_printError("Error in initParameterSet or getNumParam. Value should be %, but is instead %.\n",
                      numParam, parameter.getNumParam());
        globalError = 1;
        initParameterSetError = 1;
    }
    else
        my_print("Parameter getNumParam --- Pass\n");

    //--------------------------------------------//
    //------ getNumMixtureElements Function ------//
    //--------------------------------------------//
    if (parameter.getNumMixtureElements() != numMixtures)
    {
        my_printError("Error in initParameterSet or getNumMixtureElements. Value should be %, but is instead %.\n",
                      numMixtures, parameter.getNumMixtureElements());
        globalError = 1;
        initParameterSetError = 1;
    }
    else
        my_print("Parameter getNumMixtureElements --- Pass\n");

    //--------------------------------------------//
    //------ getStdDevSynthesisRate Function------//
    //--------------------------------------------//

    // Check proposed StdDevSynthesisRate
    for (unsigned i = 0u; i < numMixtures; i++)
    {
        if (parameter.getStdDevSynthesisRate(i, true) != stdDev[i])
        {
            my_printError("Error in initParameterSet or getStdDevSynthesisRate(proposed) for index %.", i);
            my_printError(" Value should be %, but is instead %.\n", stdDev[i], parameter.getStdDevSynthesisRate(i, true));
            error = 1;
            globalError = 1;
            initParameterSetError = 1;
        }
    }

    // Check non-proposed StdDevSynthesisRate
    for (unsigned i = 0u; i < numMixtures; i++)
    {
        if (parameter.getStdDevSynthesisRate(i, false) != stdDev[i])
        {
            my_printError("Error in initParameterSet or getStdDevSynthesisRate(non-proposed) for index %.", i);
            my_printError(" Value should be %, but is instead %.\n", stdDev[i], parameter.getStdDevSynthesisRate(i, false));
            error = 1;
            globalError = 1;
            initParameterSetError = 1;
        }
    }

    if (!error)
        my_print("Parameter getStdDevSynthesisRate --- Pass\n");
    else
        error = 0; //Reset for next function.

    // Not part of initParameterSet, but added here for convenience.
    //--------------------------------------------//
    //------ setStdDevSynthesisRate Function------//
    //--------------------------------------------//
    for (unsigned i = 0u; i < numMixtures; i++)
    {
        parameter.setStdDevSynthesisRate(i, i);
        if (parameter.getStdDevSynthesisRate(i) != i)
        {
            my_printError("Error in setStdDevSynthesisRate for index %. Value should be %, but is instead %.\n",
                          i, i, parameter.getStdDevSynthesisRate(i));
            error = 1;
            globalError = 1;
        }
    }

    if (!error)
        my_print("Parameter setStdDevSynthesisRate --- Pass\n");
    else
        error = 0; //Reset for next function.

    //-----------------------------------------------------------------//
    //------ getCurrentStdDevSynthesisRateProposalWidth Function ------//
    //-----------------------------------------------------------------//

    // This value is initialized to 0.1 in initParameterSet.
    if (parameter.getCurrentStdDevSynthesisRateProposalWidth() != 0.1)
    {
        my_printError("Error in initParameterSet or getCurrentStdDevSynthesisRateProposalWidth.");
        my_printError(" Value should be 0.1, but is instead %.\n", parameter.getCurrentStdDevSynthesisRateProposalWidth());
        globalError = 1;
        initParameterSetError = 1;
    }
    else
        my_print("Parameter getCurrentStdDevSynthesisRateProposalWidth --- Pass\n");

    // For unit testing only.
    //---------------------------------------------------------//
    //------ getNumAcceptForStdDevSynthesisRate Function ------//
    //---------------------------------------------------------//

    // This value is initialized to 0 in initParameterSet.
    if (parameter.getNumAcceptForStdDevSynthesisRate() != 0)
    {
        my_printError("Error in initParameterSet or getNumAcceptForStdDevSynthesisRate.");
        my_printError(" Value should be 0, but is instead %.\n", parameter.getNumAcceptForStdDevSynthesisRate());
        globalError = 1;
        initParameterSetError = 1;
    }
    else
        my_print("Parameter getNumAcceptForStdDevSynthesisRate --- Pass\n");

    // For unit testing only.
    //----------------------------------------//
    //------ getStdCspForIndex Function ------//
    //----------------------------------------//
    std::vector <double> tmpStd_Csp(numParam, 0.1);

    for (unsigned i = 0u; i < numParam; i++)
    {
        if (parameter.getStdCspForIndex(i) != tmpStd_Csp[i])
        {
            my_printError("Error in initParameterSet or getStdCspForIndex for index %.", i);
            my_printError(" Value should be %, but is instead %.\n", tmpStd_Csp[i], parameter.getStdCspForIndex(i));
            error = 1;
            globalError = 1;
            initParameterSetError = 1;
        }
    }

    if (!error)
        my_print("Parameter getStdCspForIndex --- Pass\n");
    else
        error = 0; //Reset for next function.

    // For unit testing only.
    //-------------------------------------------------//
    //------ getNumAcceptForCspForIndex Function ------//
    //-------------------------------------------------//

    // Because default constructor was used, maxGrouping is = 22.
    // TODO: Make this more dynamic / tested with other settings?
    unsigned maxGrouping = 22;
    std::vector <unsigned> tmpNumAcceptForCsp(maxGrouping, 0u);

    for (unsigned i = 0u; i < maxGrouping; i++)
    {
        if (parameter.getNumAcceptForCspForIndex(i) != tmpNumAcceptForCsp[i])
        {
            my_printError("Error in initParameterSet or getNumAcceptForCspForIndex for index %.", i);
            my_printError(" Value should be %, but is instead %.\n", tmpNumAcceptForCsp[i], parameter.getNumAcceptForCspForIndex(i));
            error = 1;
            globalError = 1;
            initParameterSetError = 1;
        }
    }

    if (!error)
        my_print("Parameter getNumAcceptForCspForIndex --- Pass\n");
    else
        error = 0; //Reset for next function.

    /* TODO NOTE: setNumMutationSelectionValues is accessed in initParameterSet
     * and through this function numMutationCategories and numSelectionCategories is changed.
     * This function may itself need to be tested as a middleman function.
    */

    //-----------------------------------------------//
    //------ getNumMutationCategories Function ------//
    //-----------------------------------------------//

    // Because mutationSelectionState is allUnique, by initParameterSet the numMutationCategories should be = numMixtures.
    // TODO: Make this more dynamic / tested with other settings?
    if (parameter.getNumMutationCategories() != numMixtures)
    {
        my_printError("Error in initParameterSet or getNumMutationCategories. Value should be %, but is instead %.\n",
                      numMixtures, parameter.getNumMutationCategories());
        globalError = 1;
        initParameterSetError = 1;
    }
    else
        my_print("Parameter getNumMutationCategories --- Pass\n");

    //------------------------------------------------//
    //------ getNumSelectionCategories Function ------//
    //------------------------------------------------//

    // For use in future unit testing as well
    unsigned numSelectionCategories = parameter.getNumSelectionCategories();

    // Because mutationSelectionState is allUnique, by initParameterSet the numSelectionCategories should be = numMixtures.
    // TODO: Make this more dynamic / tested with other settings?
    if (numSelectionCategories != numMixtures)
    {
        my_printError("Error in initParameterSet or getNumSelectionCategories. Value should be %, but is instead %.\n",
                      numMixtures, numSelectionCategories);
        globalError = 1;
        initParameterSetError = 1;
    }
    else
        my_print("Parameter getNumSelectionCategories --- Pass\n");

    //----------------------------------------------------//
    //------ getNumSynthesisRateCategories Function ------//
    //----------------------------------------------------//

    // Because mutationSelectionState is allUnique, by initParameterSet the numSynthesisRateCategories should be = numMixtures.
    // TODO: Make this more dynamic / tested with other settings?
    if (parameter.getNumSynthesisRateCategories() != numMixtures)
    {
        my_printError("Error in initParameterSet or getNumSynthesisRateCategories.");
        my_printError(" Value should be %, but is instead %.\n", numMixtures, parameter.getNumSynthesisRateCategories());
        globalError = 1;
        initParameterSetError = 1;
    }
    else
        my_print("Parameter getNumSynthesisRateCategories--- Pass\n");

    /* TODO NOTE: initCategoryDefinitions is accessed in initParameterSet
     * and through this function categories.delM, categories.delEta,
     * mutationIsInMixture, and selectionIsInMixture are changed.
     * This function may itself need to be tested as a middleman function.
    */

    /* Note: While the categories' .delM and .delEta values are checked,
     * The categories variable itself is not checked. Possible TODO.
    */

    //-----------------------------------------//
    //------ getMutationCategory Function------//
    //-----------------------------------------//

    // Because mutationSelectionState is allUnique,
    // by initParameterSet each category's delM = numMixtures
    // TODO: Make this more dynamic / tested with other settings?

    for (unsigned i = 0u; i < numMixtures; i++)
    {
        if (parameter.getMutationCategory(i) != i)
        {
            my_printError("Error in initParameterSet or getMutationCategory for index %.", i);
            my_printError(" Value should be %, but is instead %.\n", i, parameter.getMutationCategory(i));
            error = 1;
            globalError = 1;
            initParameterSetError = 1;
        }
    }

    if (!error)
        my_print("Parameter getMutationCategory --- Pass\n");
    else
        error = 0; //Reset for next function.

    //------------------------------------------//
    //------ getSelectionCategory Function------//
    //------------------------------------------//

    // Because mutationSelectionState is allUnique,
    // by initParameterSet each category's delEta = numMixtures
    // TODO: Make this more dynamic / tested with other settings?

    for (unsigned i = 0u; i < numMixtures; i++)
    {
        if (parameter.getSelectionCategory(i) != i)
        {
            my_printError("Error in initParameterSet or getSelectionCategory for index %.", i);
            my_printError(" Value should be %, but is instead %.\n", i, parameter.getSelectionCategory(i));
            error = 1;
            globalError = 1;
            initParameterSetError = 1;
        }
    }

    if (!error)
        my_print("Parameter getSelectionCategory --- Pass\n");
    else
        error = 0; //Reset for next function.

    //----------------------------------------------------------//
    //------ getMixtureElementsOfMutationCategory Function------//
    //----------------------------------------------------------//

    /* Because mutationSelectionState is allUnique, by initParameterSet
     * the mutationIsInMixture should be equal to a corresponding vector of unsigned integers
     * defined in the loop below -- but not for other mutation selection states.
    */
    // TODO: Make this more dynamic / tested with other settings?
    std::vector <std::vector <unsigned> > tmp(numMixtures);

    for (unsigned i = 0u; i < numMixtures; i++)
    {
        tmp[i].push_back(i);
    }

    for (unsigned i = 0u; i < numMixtures; i++)
    {
        if (parameter.getMixtureElementsOfMutationCategory(i) != tmp[i])
        {
            my_printError("Error in initParameterSet or getMixtureElementsOfMutationCategory for index %.\n", i);
            error = 1;
            globalError = 1;
            initParameterSetError = 1;
        }
    }

    if (!error)
        my_print("Parameter getMixtureElementsOfMutationCategory --- Pass\n");
    else
        error = 0; //Reset for next function.

    //-----------------------------------------------------------//
    //------ getMixtureElementsOfSelectionCategory Function------//
    //-----------------------------------------------------------//

    /* Because mutationSelectionState is allUnique, by initParameterSet
     * the mutationIsInMixture should be equal to the vector of unsigned integers
     * defined above -- but not for other mutation selection states.
    */
    // TODO: Make this more dynamic / tested with other settings?

    for (unsigned i = 0u; i < numMixtures; i++)
    {
        if (parameter.getMixtureElementsOfSelectionCategory(i) != tmp[i])
        {
            my_printError("Error in initParameterSet or getMixtureElementsOfSelectionCategory for index %.\n", i);
            error = 1;
            globalError = 1;
            initParameterSetError = 1;
        }
    }

    if (!error)
        my_print("Parameter getMixtureElementsOfSelectionCategory --- Pass\n");
    else
        error = 0; //Reset for next function.

    //--------------------------------------------//
    //------ getCategoryProbability Function------//
    //--------------------------------------------//

    // Each value is initialized to 1.0 / numMixtures in initParameterSet.
    for (unsigned i = 0u; i < numMixtures; i++)
    {
        if (parameter.getCategoryProbability(i) != 1.0/(double)numMixtures)
        {
            my_printError("Error in initParameterSet or getCategoryProbability for index %.", i);
            my_printError(" Value should be %, but is instead %.\n", 1.0/(double)numMixtures, parameter.getCategoryProbability(i));
            error = 1;
            globalError = 1;
            initParameterSetError = 1;
        }
    }

    if (!error)
        my_print("Parameter getCategoryProbability --- Pass\n");
    else
        error = 0; //Reset for next function.

    // Not part of initParameterSet, but added here for convenience.
    //--------------------------------------------//
    //------ setCategoryProbability Function------//
    //--------------------------------------------//
    for (unsigned i = 0u; i < numMixtures; i++)
    {
        parameter.setCategoryProbability(i, (double)i/(double)numMixtures);
        if (parameter.getCategoryProbability(i) != (double)i/(double)numMixtures)
        {
            my_printError("Error in setCategoryProbability for index %. Value should be %, but is instead %.\n",
                          i, (double)i/(double)numMixtures, parameter.getCategoryProbability(i));
            error = 1;
            globalError = 1;
        }
    }

    if (!error)
        my_print("Parameter setCategoryProbability --- Pass\n");
    else
        error = 0; //Reset for next function.

    //--------------------------------------//
    //------ getSynthesisRate Function------//
    //--------------------------------------//

    /* Each value is initialized to 0.0 in initParameterSet.
     * Note that numSelectionCategories = numMixtures since allUnique; therefore,
     * the outer loop only iterates once.
    */

    // Check proposed SynthesisRate
    for (unsigned i = 0u; i < numSelectionCategories; i++)
    {
        for (unsigned j = 0u; j < numGenes; j++)
        {
            if (parameter.getSynthesisRate(j, i, true) != 0.0)
            {
                my_printError("Error in initParameterSet or getSynthesisRate(proposed) for index % of mixture %.", j, i);
                my_printError(" Value should be 0.0, but is instead %.\n", parameter.getSynthesisRate(j, i, true));
                error = 1;
                globalError = 1;
                initParameterSetError = 1;
            }
        }
    }

    // Check non-proposed SynthesisRate
    for (unsigned i = 0u; i < numSelectionCategories; i++)
    {
        for (unsigned j = 0u; j < numGenes; j++)
        {
            if (parameter.getSynthesisRate(j, i, false) != 0.0)
            {
                my_printError("Error in initParameterSet or getSynthesisRate(non-proposed) for index % of mixture %.", j, i);
                my_printError(" Value should be 0.0, but is instead %.\n", parameter.getSynthesisRate(j, i, false));
                error = 1;
                globalError = 1;
                initParameterSetError = 1;
            }
        }
    }

    if (!error)
        my_print("Parameter getSynthesisRate --- Pass\n");
    else
        error = 0; //Reset for next function.

    // Not part of initParameterSet, but added here for convenience.
    //--------------------------------------//
    //------ setSynthesisRate Function------//
    //--------------------------------------//

    for (unsigned i = 0u; i < numSelectionCategories; i++)
    {
        for (unsigned j = 0u; j < numGenes; j++)
        {
            parameter.setSynthesisRate(j, j, i);
            if (parameter.getSynthesisRate(j, i, false) != j)
            {
                my_printError("Error in setSynthesisRate for index % of mixture %. Value should be %, but is instead %.\n",
                              j, i, j, parameter.getSynthesisRate(j, i, false));
                error = 1;
                globalError = 1;
            }
        }
    }

    if (!error)
        my_print("Parameter setSynthesisRate --- Pass\n");
    else
        error = 0; //Reset for next function.

    //----------------------------------------------//
    //------ getSynthesisRateCategory Function------//
    //----------------------------------------------//

    // Because mutationSelectionState is allUnique,
    // by initParameterSet each category's delEta = numMixtures
    // TODO: Make this more dynamic / tested with other settings?

    for (unsigned i = 0u; i < numMixtures; i++)
    {
        if (parameter.getSynthesisRateCategory(i) != i)
        {
            my_printError("Error in initParameterSet or getSynthesisRateCategory for index %.", i);
            my_printError(" Value should be %, but is instead %.\n", i, parameter.getSynthesisRateCategory(i));
            error = 1;
            globalError = 1;
            initParameterSetError = 1;
        }
    }

    if (!error)
        my_print("Parameter getSynthesisRateCategory --- Pass\n");
    else
        error = 0; //Reset for next function.

    //--------------------------------------------------//
    //------ getNumAcceptForSynthesisRate Function------//
    //--------------------------------------------------//

    // Each value is initialized to 0 in initParameterSet.

    for (unsigned i = 0u; i < numMixtures; i++)
    {
        unsigned expressionCategory = parameter.getSynthesisRateCategory(i);
        for (unsigned j = 0u; j < numGenes; j++)
        {
            if (parameter.getNumAcceptForSynthesisRate(expressionCategory, j) != 0)
            {
                my_printError("Error in initParameterSet or getNumAcceptForSynthesisRate");
                my_printError(" for index % of expression category %. Value should be 0, but is instead %.\n",
                              j, expressionCategory, parameter.getNumAcceptForSynthesisRate(expressionCategory, j));
                error = 1;
                globalError = 1;
                initParameterSetError = 1;
            }
        }
    }

    if (!error)
        my_print("Parameter getNumAcceptForSynthesisRate --- Pass\n");
    else
        error = 0; //Reset for next function.

    //---------------------------------------------------//
    //------ getSynthesisRateProposalWidth Function------//
    //---------------------------------------------------//

    /* Each value is initialized to 5 in initParameterSet.
     * Note that numSelectionCategories = numMixtures since allUnique; therefore,
     * the outer loop only iterates once.
    */
    for (unsigned i = 0u; i < numSelectionCategories; i++)
    {
        for (unsigned j = 0u; j < numGenes; j++)
        {
            if (parameter.getSynthesisRateProposalWidth(j, i) != 5)
            {
                my_printError("Error in initParameterSet or getSynthesisRateProposalWidth for index % of mixture %.", j, i);
                my_printError(" Value should be 5, but is instead %.\n", parameter.getSynthesisRateProposalWidth(j, i));
                error = 1;
                globalError = 1;
                initParameterSetError = 1;
            }
        }
    }

    if (!error)
        my_print("Parameter getSynthesisRateProposalWidth --- Pass\n");
    else
        error = 0; //Reset for next function.

    //----------------------------------------------------------//
    //------ getCurrentSynthesisRateProposalWidth Function------//
    //----------------------------------------------------------//

    // Each value is initialized to 5 in initParameterSet.

    for (unsigned i = 0u; i < numMixtures; i++)
    {
        unsigned expressionCategory = parameter.getSynthesisRateCategory(i);
        for (unsigned j = 0u; j < numGenes; j++)
        {
            if (parameter.getCurrentSynthesisRateProposalWidth(expressionCategory, j) != 5)
            {
                my_printError("Error in initParameterSet or getCurrentSynthesisRateProposalWidth");
                my_printError(" for index % of expression category %. Value should be 5, but is instead %.\n",
                              j, expressionCategory, parameter.getCurrentSynthesisRateProposalWidth(expressionCategory, j));
                error = 1;
                globalError = 1;
                initParameterSetError = 1;
            }
        }
    }

    if (!error)
        my_print("Parameter getCurrentSynthesisRateProposalWidth --- Pass\n");
    else
        error = 0; //Reset for next function.

    //------------------------------------------------------------------------//
    //------ End of Unit Testing for initParameterSet-related Functions ------//
    //------------------------------------------------------------------------//

    if (!initParameterSetError)
        my_print("Parameter initParameterSet --- Pass\n");

    /* Section 2:
     * Group List functions: 4 functions tested in total.
     * get/setGroupList, getGroupListSize, getGrouping
    */

    std::vector <std::string> fullGroupList = {"GCA", "GCC", "GCG", "GCT", "TGC", "TGT", "GAC", "GAT", "GAA", "GAG",
                                               "TTC", "TTT", "GGA", "GGC", "GGG", "GGT", "CAC", "CAT", "ATA", "ATC",
                                               "ATT", "AAA", "AAG", "CTA", "CTC", "CTG", "CTT", "TTA", "TTG", "ATG",
                                               "AAC", "AAT", "CCA", "CCC", "CCG", "CCT", "CAA", "CAG", "AGA", "AGG",
                                               "CGA", "CGC", "CGG", "CGT", "TCA", "TCC", "TCG", "TCT", "ACA", "ACC",
                                               "ACG", "ACT", "GTA", "GTC", "GTG", "GTT", "TGG", "TAC", "TAT", "AGC",
                                               "AGT"};

    std::vector <std::string> twoCodonGroupList = {"TTT", "TTC", "TGT", "TGC", "TAT", "TAC", "CAA", "CAG", "AAT",
                                                   "AAC", "CAT", "CAC", "GAA", "GAG", "GAT", "GAC", "AAA", "AAG"};

    //----------------------------------------//
    //------ get/setGroupList Functions ------//
    //----------------------------------------//
    parameter.setGroupList(fullGroupList);

    if (parameter.getGroupList() != fullGroupList)
    {
        my_printError("Error in setGroupList or getGroupList with a full group list.\n");
        error = 1;
        globalError = 1;
    }

    parameter.setGroupList(twoCodonGroupList);
    if (parameter.getGroupList() != twoCodonGroupList)
    {
        my_printError("Error in setGroupList or getGroupList with only two-codon AAs in the group list.\n");
        error = 1;
        globalError = 1;
    }

    if (!error)
        my_print("Parameter get/setGroupList --- Pass\n");
    else
        error = 0; //Reset for next function.

    //----------------------------------------//
    //------ getGroupListSize Functions ------//
    //----------------------------------------//
    if (parameter.getGroupListSize() != 18)
    {
        my_printError("Error in getGroupListSize. Size should be 18, but is instead %.\n", parameter.getGroupListSize());
        globalError = 1;
    }
    else
        my_print("Parameter getGroupListSize --- Pass\n");

    //-----------------------------------//
    //------ getGrouping Functions ------//
    //-----------------------------------//
    // We will only test 6 indices

    if (parameter.getGrouping(0) != "TTT")
    {
        my_printError("Error in getGrouping at index 0. Value should be TTT, but is instead %.\n", parameter.getGrouping(0));
        error = 1;
        globalError = 1;
    }
    if (parameter.getGrouping(9) != "AAC")
    {
        my_printError("Error in getGrouping at index 9. Value should be AAC, but is instead %.\n", parameter.getGrouping(9));
        error = 1;
        globalError = 1;
    }
    if (parameter.getGrouping(14) != "GAT")
    {
        my_printError("Error in getGrouping at index 14. Value should be GAT, but is instead %.\n", parameter.getGrouping(14));
        error = 1;
        globalError = 1;
    }
    if (parameter.getGrouping(4) != "TAT")
    {
        my_printError("Error in getGrouping at index 4. Value should be TAT, but is instead %.\n", parameter.getGrouping(28));
        error = 1;
        globalError = 1;
    }
    if (parameter.getGrouping(17) != "AAG")
    {
        my_printError("Error in getGrouping at index 17. Value should be AAG, but is instead %.\n", parameter.getGrouping(31));
        error = 1;
        globalError = 1;
    }
    if (parameter.getGrouping(1) != "TTC")
    {
        my_printError("Error in getGrouping at index 1. Value should be TTC, but is instead %.\n", parameter.getGrouping(60));
        error = 1;
        globalError = 1;
    }

    if (!error)
        my_print("Parameter getGrouping --- Pass\n");
    else
        error = 0; //Reset for next function.

    /* Section 3:
     * InitializeSynthesisRate Function
     */

    //----------------------------------------------//
    //------ InitializeSynthesisRate Function ------//
    //----------------------------------------------//
    parameter.InitializeSynthesisRate(genome, stdDev[0]);

    // This call changes currentSynthesisRateLevel, std_phi, and numAcceptForSynthesisRate.
    // These functions must now be checked.
    for (unsigned category = 0u; category < numSelectionCategories; category++)
    {
        for (unsigned j = 0u; j < numGenes; j++)
        {
            // TODO: Check if currentSynthesisRateLevel is set correctly
            // TODO: Check the follow functions: calculateSCUO, Parameter::randLogNorm, quickSortPair, pivotPair
            // currentSynthesisRateLevel[category][index[j]] = expression[j];

            // Check if std_phi = 0.1 as set in InitializeSynthesisRate
            if (parameter.getCurrentSynthesisRateProposalWidth(category, j) != 0.1)
            {
                my_printError("Error in InitializeSynthesisRate -- std_phi is not set correctly.\n");
                my_printError(" Value at index % of expression category % should be 0.1, but is instead %.\n",
                              j, category, parameter.getCurrentSynthesisRateProposalWidth(category, j));
                error = 1;
                globalError = 1;
            }

            //Check if numAcceptForSynthesisRate = 0 as set in InitializeSynthesisRate
            if (parameter.getNumAcceptForSynthesisRate(category, j) != 0)
            {
                my_printError("Error in InitializeSynthesisRate -- numAcceptForSynthesisRate is not set correctly.\n");
                my_printError(" Value at index % of expression category % should be 0, but is instead %.\n",
                              j, category, parameter.getNumAcceptForSynthesisRate(category, j));
                error = 1;
                globalError = 1;
            }
        }
    }

    if (!error)
        my_print("Parameter InitializeSynthesisRate --- Pass\n");
    else
        error = 0; //Reset for next function.

    /* Section 4:
     * Other functions: 4 functions tested in total.
     * get/setLastIteration, get/setNumObservedPhiSets
    */

    //--------------------------------------------//
    //------ get/setLastIteration Functions ------//
    //--------------------------------------------//

    //lastIteration should be initialized to 0 in the constructor.
    if (parameter.getLastIteration() != 0)
    {
        my_printError("Error in getLastIteration. Value should be 0, but is instead %.\n", parameter.getLastIteration());
        error = 1;
        globalError = 1;
    }

    parameter.setLastIteration(9);
    if (parameter.getLastIteration() != 9)
    {
        my_printError("Error in setLastIteration or getLastIteration. Value should be 9, but is instead %.\n",
                      parameter.getLastIteration());
        error = 1;
        globalError = 1;
    }

    if (!error)
        my_print("Parameter get/setLastIteration --- Pass\n");
    else
        error = 0; //Reset for next function.

    //-------------------------------------------------//
    //------ get/setNumObservedPhiSets Functions ------//
    //-------------------------------------------------//

    //obsPhiSets should be initialized to 0 in the constructor.
    if (parameter.getNumObservedPhiSets() != 0)
    {
        my_printError("Error in getNumObservedPhiSets. Value should be 0, but is instead %.\n",
                      parameter.getNumObservedPhiSets());
        error = 1;
        globalError = 1;
    }

    parameter.setNumObservedPhiSets(22);
    if (parameter.getNumObservedPhiSets() != 22)
    {
        my_printError("Error in setNumObservedPhiSets or getNumObservedPhiSets. Value should be 22, but is instead %.\n",
                      parameter.getNumObservedPhiSets());
        error = 1;
        globalError = 1;
    }

    if (!error)
        my_print("Parameter get/setNumObservedPhiSets --- Pass\n");
    else
        error = 0; //Reset for next function.

    return globalError;
}


/* TODO: Rework or remove!
int testParameterWithFile(std::string filename)
{
    Parameter parameter;

    parameter.initBaseValuesFromFile(filename);

    return 0;
}
*/


/* testCovarianceMatrix (RCPP EXPOSED)
 * Arguments: None
 * Performs Unit Testing on functions within CovarianceMatrix.cpp
 * that are not exposed to RCPP already.
 * Returns 0 if successful, 1 if error found.
*/
int testCovarianceMatrix()
{
    CovarianceMatrix covM; //Default constructor sets numVariates to 2.
    int error = 0;
    int globalError = 0;

    //-----------------------------------------------------------//
    //------ getCovMatrix & initCovarianceMatrix Functions ------//
    //-----------------------------------------------------------//

    std::vector <double> covM2 = {0.0025, 0, 0, 0, \
                                  0, 0.0025, 0, 0, \
                                  0, 0, 0.0025, 0, \
                                  0, 0, 0, 0.0025};

    covM.initCovarianceMatrix(4);

    for (unsigned i = 0u; i < 16; i++)
    {
        // Compare, for each position of the vector of doubles, if the statically created covM2
        // equals the initialized and extracted covM
        if (covM2[i] != covM.getCovMatrix()->at(i))
        {
            my_printError("Error in getCovMatrix or initCovarianceMatrix:");
            my_printError("at index % matrix extracted should return % but instead returns %.\n",
                          i, covM2[i], covM.getCovMatrix()->at(i));
            error = 1;
            globalError = 1;
        }
    }

    if (!error)
        my_print("CovarianceMatrix getCovMatrix & initCovarianceMatrix --- Pass\n");
    else
        error = 0; //Reset for next function.

    //------------------------------//
    //------ setDiag Function ------//
    //------------------------------//
    covM.setDiag(3.14);

    covM2 = {3.14, 0, 0, 0, \
             0, 3.14, 0, 0, \
             0, 0, 3.14, 0, \
             0, 0, 0, 3.14};

    for (unsigned i = 0u; i < 16; i++)
    {
        // Compare, for each position of the vector of doubles, if the statically created covM2
        // equals the set covM
        if (covM2[i] != covM.getCovMatrix()->at(i))
        {
            my_printError("Error in setDiag: at index % matrix extracted should return % but instead returns %.\n",
                          i, covM2[i], covM.getCovMatrix()->at(i));
            error = 1;
            globalError = 1;
        }
    }

    if (!error)
        my_print("CovarianceMatrix setDiag --- Pass\n");
    else
        error = 0; //Reset for next function.

    //-----------------------------------------------------------------//
    //------ getCholeskyMatrix & choleskyDecomposition Functions ------//
    //-----------------------------------------------------------------//
    covM.choleskyDecomposition();

    // Perform Cholesky decomposition on covM2
    // (This should be the same code as used in Cholesky decomposition in CovarianceMatrix.cpp)
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < (i + 1); j++)
        {
            double LsubstractSum = 0.0;
            for (int k = 0; k < j; k++)
            {
                LsubstractSum += covM2[i * 4 + k] * covM2[j * 4 + k];
            }
            covM2[i * 4 + j] = (i == j) ? std::sqrt(covM2[i * 4 + i] - LsubstractSum) :
                               (1.0 / covM2[j * 4 + j]) * (covM2[i * 4 + j] - LsubstractSum);
        }
    }

    for (unsigned i = 0u; i < 16; i++)
    {
        // Compare, for each position of the vector of doubles, if the statically created covM2
        // equals the set covM for its Cholesky matrix
        if (covM2[i] != covM.getCholeskyMatrix()->at(i))
        {
            my_printError("Error in getCholeskyMatrix or choleskyDecomposition:");
            my_printError(" at index % matrix extracted should return % but instead returns %.\n",
                          i, covM2[i], covM.getCholeskyMatrix()->at(i));
            error = 1;
            globalError = 1;
        }
    }

    if (!error)
        my_print("CovarianceMatrix getCholeskyMatrix & choleskyDecomposition --- Pass\n");
    else
        error = 0; //Reset for next function.

    //-------------------------------------//
    //------ getNumVariates Function ------//
    //-------------------------------------//
    covM.getNumVariates();

    if (covM.getNumVariates() != 4)
    {
        my_printError("Error in getNumVariates. Function should return 4, but returns %.\n", covM.getNumVariates());
        globalError = 1;
    }
    else
        my_print("CovarianceMatrix getNumVariates --- Pass\n");

    //----------------------------------//
    //------ == Operator Function ------//
    //----------------------------------//
    CovarianceMatrix covMcp; //Default constructor sets numVariates to 2.

    covMcp = covM;

    if (!(covMcp == covM))
    {
        my_printError("Error in CovarianceMatrix == operator. Fails to be equivalent.\n");
        globalError = 1;
    }
    else
        my_print("CovarianceMatrix == operator --- Pass\n");

    //TODO: Test these final two functions.
    //-------------------------------------------------------------//
    //------ transformIidNumbersIntoCovaryingNumbers Function ------//
    //-------------------------------------------------------------//
    //TODO: Should implement Parameter unit testing for RandNorm before doing this!

    //Size is numCodons * (numMutationCategories + numSelectionCategories for ROCParameter
    //Each value is based on randNorm(0.0, 1.0)
    //covM.transformIidNumbersIntoCovaryingNumbers(iidTest);

    //------------------------------------------------//
    //------ calculateSampleCovariance Function ------//
    //------------------------------------------------//

    return globalError;
}


/* TODO: Rework or remove!
int testPATrace()
{
    Trace RFP; //initialize with 0 categories, 2 codon-specific parameter types
    Trace ROC;
    Trace FONSE;
    Trace PANSE;
    //int samples = 10;
    //int thinning = 10;
    //int useSamples = 100;
    int globalError = 0;

    //-----------------------------------------//
    //------ initializePATrace Function ------//
    //-----------------------------------------//
    //RFP.initializePATrace();

    //-----------------------------------------//
    //------ initializeROCTrace Function ------//
    //-----------------------------------------//
    //ROC.initializeROCTrace();

    //-------------------------------------------//
    //------ initializeFONSETrace Function ------//
    //-------------------------------------------//
    //FONSE.initializeFONSETrace();

    //-------------------------------------------//
    //------ initializePANSETrace Function ------//
    //-------------------------------------------//
    //PANSE.initializePANSETrace();

    return globalError;
}
*/


/* TODO: Rework or remove!
int testPAParameter()
{
    int error = 0;
    int globalError = 0;

     * Section 1: 1 function tested in total.
     * initPAParameterSet Function
     * and related get/set functions as a consequence of the function setup:


    //------------------------------------------//
    //------ initPAParameterSet Function ------//
    //------------------------------------------//

     * Initialize parameter:
     * Arguments: vector <double> stdDevSynthesisRate, unsigned numMixtures, vector <unsigned> geneAssignment,
     *           vector <vector <unsigned>> mixtureDefinitionMatrix, bool splitSer, string mutationSelectionState
     *
     * Thus, let:


    Genome genome;
    //genome.readRFPData("/Users/hollisbui/RibModelDev/data/rfp/rfp.counts.by.codon.and.gene.GSE63789.wt.csv");
    unsigned numMixtures = 3;
    std::vector <double> stdDev(numMixtures, 1);
    unsigned numGenes = genome.getGenomeSize();
    std::vector <unsigned> geneAssignment(numGenes);
    if (numMixtures == 1)
    {
        for (unsigned i = 0u; i < numGenes; i++)
        {
            geneAssignment[i] = 0u;
        }
    }
    else if (numMixtures == 3)
    {
        for (unsigned i = 0u; i < numGenes; i++)
        {
            if (i < 961) geneAssignment[i] = 0u;
            else if (i < 1418) geneAssignment[i] = 1u;
            else geneAssignment[i] = 0u;
        }
    }
    std::vector<std::vector<unsigned>> mixtureDefinitionMatrix;
    bool splitSer = true;
    std::string mutationSelectionState = Parameter::allUnique;

    PAParameter parameter(stdDev, numMixtures, geneAssignment, mixtureDefinitionMatrix, splitSer, mutationSelectionState);

     * This constructor in turn calls two functions: initParameterSet() and initPAParameterSet().
     * initParameterSet should have been tested in testParameter(), above, but we must now
     * test initPAParameterSet
     *
     * Thus, unit testing is done in order of variable changed:
     * numParam, currentCodonSpecificParameter, proposedCodonSpecificParameter, std_csp, and groupList.
     * This also introduces a level of uncertainty in what may be wrong, and thus an error in the following
     * unit testing checks may be a result of the checking function or initParameterSet.


    // numParam is set to 61 in initPAParameterSet.
    unsigned numParam = parameter.getNumParam();
    if (numParam != 61)
    {
        my_printError("Error in initPAParameterSet -- numParam is not set correctly.");
        my_printError(" Value should be 61 but is instead %.\n", numParam);
        error = 1;
        globalError = 1;
    }

    // TODO: check the changed:
    // currentCodonSpecificParameter
    // proposedCodonSpecificParameter

    // std_csp is set to 0.1 for each index in initPAParameterSet.
    for (unsigned i = 0u; i < numParam; i++)
    {
        if (parameter.getStdCspForIndex(i) != 0.1)
        {
            my_printError("Error in InitializeSynthesisRate -- std_csp is not set correctly.");
            my_printError(" Value at index % should be 0.1, but is instead %.\n", i, parameter.getStdCspForIndex(i));
            error = 1;
            globalError = 1;
        }
    }

    // groupList is set to the same as this temporary group list in initPAParameterSet.
    std::vector <std::string> tmpGroupList = {"GCA", "GCC", "GCG", "GCT", "TGC", "TGT", "GAC", "GAT", "GAA", "GAG",
                                              "TTC", "TTT", "GGA", "GGC", "GGG", "GGT", "CAC", "CAT", "ATA", "ATC",
                                              "ATT", "AAA", "AAG", "CTA", "CTC", "CTG", "CTT", "TTA", "TTG", "ATG",
                                              "AAC", "AAT", "CCA", "CCC", "CCG", "CCT", "CAA", "CAG", "AGA", "AGG",
                                              "CGA", "CGC", "CGG", "CGT", "TCA", "TCC", "TCG", "TCT", "ACA", "ACC",
                                              "ACG", "ACT", "GTA", "GTC", "GTG", "GTT", "TGG", "TAC", "TAT", "AGC",
                                              "AGT"};

    if (parameter.getGroupList() != tmpGroupList)
    {
        my_printError("Error in initPAParameterSet -- groupList is not set correctly.\n");
        error = 1;
        globalError = 1;
    }

    if (!error)
        my_print("PAParameter initPAParameterSet --- Pass\n");
    else
        error = 0; //Reset for next function.

    //parameter.InitializeSynthesisRate(genome, stdDev[0]);

    return globalError;
}
*/


/* testMCMCAlgorithm (RCPP EXPOSED)
 * Arguments: None
 * Performs Unit Testing on functions within MCMCAlgorithm.cpp
 * that are not exposed to RCPP already.
 * Returns 0 if successful, 1 if error found.
*/
int testMCMCAlgorithm()
{
    unsigned samples = 10;
    unsigned thinning = 10;
    int error = 0;
    int globalError = 0;

    MCMCAlgorithm mcmc = MCMCAlgorithm(samples, thinning, 10, true, true, true);

    if (!mcmc.isEstimateSynthesisRate())
    {
        my_printError("Error in isEstimateSynthesisRate. Function should return true, but returns false.\n");
        error = 1;
        globalError = 1;
    }
    my_print("checked mcmc.isEstimateSynthesisRate(default)\n");

    mcmc.setEstimateSynthesisRate(false);
    if (mcmc.isEstimateSynthesisRate())
    {
        my_printError("Error in isEstimateSynthesisRate or setEstimateSynthesisRate.");
        my_printError(" Function should return false, but returns true.\n");
        error = 1;
        globalError = 1;
    }
    my_print("checked mcmc.isEstimateSynthesisRate(false)\n");

    mcmc.setEstimateSynthesisRate(true);
    if (!mcmc.isEstimateSynthesisRate())
    {
        my_printError("Error in isEstimateSynthesisRate or setEstimateSynthesisRate.");
        my_printError(" Function should return true, but returns false.\n");
        error = 1;
        globalError = 1;
    }
    my_print("checked mcmc.isEstimateSynthesisRate(true)\n");

    if (!error)
        my_print("MCMCAlgorithm is/setEstimateSynthesisRate --- Pass\n");
    else
        error = 0; //Reset for next function.

    //------------------------------------------------------------//
    //------ is/setEstimateCodonSpecificParameter Functions ------//
    //------------------------------------------------------------//
    if (!mcmc.isEstimateCodonSpecificParameter())
    {
        my_printError("Error in isEstimateCodonSpecificParameter. Function should return true, but returns false.\n");
        error = 1;
        globalError = 1;
    }

    mcmc.setEstimateCodonSpecificParameter(false);
    if (mcmc.isEstimateCodonSpecificParameter())
    {
        my_printError("Error in isEstimateCodonSpecificParameter or setEstimateCodonSpecificParameter.");
        my_printError(" Function should return false, but returns true.\n");
        error = 1;
        globalError = 1;
    }

    mcmc.setEstimateCodonSpecificParameter(true);
    if (!mcmc.isEstimateCodonSpecificParameter())
    {
        my_printError("Error in isEstimateCodonSpecificParameter or setEstimateCodonSpecificParameter.");
        my_printError(" Function should return true, but returns false.\n");
        error = 1;
        globalError = 1;
    }

    if (!error)
        my_print("MCMCAlgorithm is/setEstimateCodonSpecificParameter --- Pass\n");
    else
        error = 0; //Reset for next function.

    //----------------------------------------------------//
    //------ is/setEstimateHyperParameter Functions ------//
    //----------------------------------------------------//
    if (!mcmc.isEstimateHyperParameter())
    {
        my_printError("Error in isEstimateHyperParameter. Function should return true, but returns false.\n");
        error = 1;
        globalError = 1;
    }

    mcmc.setEstimateHyperParameter(false);
    if (mcmc.isEstimateHyperParameter())
    {
        my_printError("Error in isEstimateHyperParameter or setEstimateHyperParameter.");
        my_printError(" Function should return false, but returns true.\n");
        error = 1;
        globalError = 1;
    }

    mcmc.setEstimateHyperParameter(true);
    if (!mcmc.isEstimateHyperParameter())
    {
        my_printError("Error in isEstimateHyperParameter or setEstimateHyperParameter.");
        my_printError(" Function should return true, but returns false.\n");
        error = 1;
        globalError = 1;
    }

    if (!error)
        my_print("MCMCAlgorithm is/setEstimateHyperParameter --- Pass\n");
    else
        error = 0; //Reset for next function.

    //-------------------------------------------------------//
    //------ is/setEstimateMixtureAssignment Functions ------//
    //-------------------------------------------------------//

    /* NOTE: By default, both constructors initialize estimateMixtureAssignment to true,
    // although it is not one of the arguments. */

    if (!mcmc.isEstimateMixtureAssignment())
    {
        my_printError("Error in isEstimateMixtureAssignment. Function should return true, but returns false.\n");
        error = 1;
        globalError = 1;
    }

    mcmc.setEstimateMixtureAssignment(false);
    if (mcmc.isEstimateMixtureAssignment())
    {
        my_printError("Error in isEstimateMixtureAssignment or setEstimateMixtureAssignment.");
        my_printError(" Function should return false, but returns true.\n");
        error = 1;
        globalError = 1;
    }

    mcmc.setEstimateMixtureAssignment(true);
    if (!mcmc.isEstimateMixtureAssignment())
    {
        my_printError("Error in isEstimateMixtureAssignment or setEstimateMixtureAssignment.");
        my_printError(" Function should return true, but returns false.\n");
        error = 1;
        globalError = 1;
    }

    if (!error)
        my_print("MCMCAlgorithm is/setEstimateMixtureAssignment --- Pass\n");
    else
        error = 0; //Reset for next function.

    //-------------------------------------------//
    //------ get/setStepsToAdapt Functions ------//
    //-------------------------------------------//

    // NOTE: By default, both constructors initialize stepsToAdapt to -1
    if (mcmc.getStepsToAdapt() != -1)
    {
        my_printError("Error in getStepsToAdapt. Function should return -1, but returns %.\n", mcmc.getStepsToAdapt());
        error = 1;
        globalError = 1;
    }

    mcmc.setStepsToAdapt(52);
    if (mcmc.getStepsToAdapt() != 52)
    {
        my_printError("Error in getStepsToAdapt or setStepsToAdapt. Function should return 52, but returns %.\n",
                      mcmc.getStepsToAdapt());
        error = 1;
        globalError = 1;
    }

    // Intentional error checking: Should print an error message with no change to stepsToAdapt
    mcmc.setStepsToAdapt(101);
    if (mcmc.getStepsToAdapt() != 52)
    {
        my_printError("Error in getStepsToAdapt or setStepsToAdapt.");
        my_printError(" Function should return 52, with no change, but returns %.\n", mcmc.getStepsToAdapt());
        error = 1;
        globalError = 1;
    }

    if (!error)
        my_print("MCMCAlgorithm get/setStepsToAdapt --- Pass\n");
    else
        error = 0; //Reset for next function.

    //-------------------------------------------//
    //------ getLogPosteriorTrace Function ------//
    //-------------------------------------------//

    // NOTE: By default, both constructors initialize likelihoodTrace to a vector with
    // a size of samples + 1 zeroes.
    std::vector <double> posteriorTrace = mcmc.getLogPosteriorTrace();
    std::vector <double> tmp;
    tmp.resize(samples+1);

    if (tmp != posteriorTrace)
    {
        my_printError("Error in getLogPosteriorTrace. Function should return a vector of % + 1 zeroes.\n", samples);
        globalError = 1;
    }
    else
        my_print("MCMCAlgorithm getLogPosteriorTrace --- Pass\n");

    //----------------------------------------------------//
    //------ getLogLikelihoodPosteriorMean Function ------//
    //----------------------------------------------------//
    //TODO. LikelihoodTrace is *not* set in this file, and therefore testing will
    //require implementation of models beforehand.

    return globalError;
}





//-----------------------------------------------------------------------------------------------------//
//---------------------------------------- R SECTION --------------------------------------------------//
//-----------------------------------------------------------------------------------------------------//

#ifndef STANDALONE
//---------------------------------//
//---------- RCPP Module ----------//
//---------------------------------//

RCPP_MODULE(Test_mod)
{
	function("testUtility", &testUtility);
	function("testSequenceSummary", &testSequenceSummary);
	function("testGene", &testGene);
	function("testGenome", &testGenome);
	function("testParameter", &testParameter);
	function("testCovarianceMatrix", &testCovarianceMatrix);
	//function("testPAParameter", &testPAParameter);
	function("testMCMCAlgorithm", &testMCMCAlgorithm);
}
#endif
