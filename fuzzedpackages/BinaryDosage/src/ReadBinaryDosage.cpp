#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include <vector>

extern const std::ios_base::openmode READWRITEBINARY;
extern const std::ios_base::openmode READBINARY;
extern const std::ios_base::openmode WRITEBINARY;

extern const int NUMBEROFBASES;
extern const unsigned short USBASE[];
extern const double DBASE[];

//***************************************************************************//
//                        Reading the header                                 //
//***************************************************************************//
// These functions read the headers to the binary dosage files

//  ************************ Constants **************************************//
// Magic word for binary dosage files
extern const int MAGICWORD;
// Format ID stored in file
extern const std::vector<std::vector<int> > FORMAT;

// Reads the base header for a binary dosage file
// Parameter filename - Name of binary dosage file
// Return - list with format and subformat
// [[Rcpp::export]]
Rcpp::List ReadBinaryDosageBaseHeader(std::string &filename) {
  std::ifstream infile;
  int magicWord;
  int formatInt;
  unsigned int ui, uj;
  int format, subformat;
  Rcpp::List retVal;

  // Open the file - if file already exists, truncates to size 0.
  // Only opens for output
  infile.open(filename.c_str(), READBINARY);
  if (!infile.good())
    return Rcpp::List::create(Rcpp::Named("error") = "Unable to open binary dosage file");

  infile.read((char *)&magicWord, sizeof(int));
  infile.read((char *)&formatInt, sizeof(int));
  if (magicWord != MAGICWORD) {
    infile.close();
    return Rcpp::List::create(Rcpp::Named("error") = "File does not appear to be a binary dosage file");
  }

  format = 0;
  subformat = 0;
  for (ui = 0; ui < FORMAT.size(); ++ui) {
    for (uj = 0; uj < FORMAT[ui].size(); ++uj) {
      if(FORMAT[ui][uj] == formatInt) {
        format = ui + 1;
        subformat = uj + 1;
        ui = FORMAT.size();
        break;
      }
    }
  }

  infile.close();
  if (format == 0)
    return Rcpp::List::create(Rcpp::Named("error") = "Unknown binary dosage file fromat");

  return Rcpp::List::create(Rcpp::Named("format") = format,
                            Rcpp::Named("subformat") = subformat);
}

// Writes the additional header info for formats 3.1 and 3.2
// Parameter filename - Name of binary dosage file
// Return - number of subjects in data
// [[Rcpp::export]]
Rcpp::List ReadBinaryDosageHeader3A(std::string &filename) {
  std::ifstream infile;
  int numSubjects = 0;
  Rcpp::List retVal;

  // Open the file for reading
  // Only opens for output
  infile.open(filename.c_str(), READBINARY);

  infile.seekg(8);
  infile.read((char *)&numSubjects, sizeof(int));

  infile.close();
  return Rcpp::List::create(Rcpp::Named("numsub") = numSubjects);
}

// Writes the additional header info for formats 3.3 and 3.4
// Parameter filename - Name of binary dosage file
// Return - StringVector with MD5 hashes of family and map files
// [[Rcpp::export]]
Rcpp::List ReadBinaryDosageHeader3B(std::string &filename) {
  std::ifstream infile;
  Rcpp::StringVector md5(2);
  Rcpp::List retVal;
  char md5hash[33];

  // Open the file for appending
  // Only opens for output
  infile.open(filename.c_str(), READBINARY);

  infile.seekg(8);
  infile.read(md5hash, 32);
  md5hash[32] = 0;
  md5[0] = md5hash;
  infile.read(md5hash, 32);
  md5[1] = md5hash;

  infile.close();
  return Rcpp::List::create(Rcpp::Named("md5") = md5);
}

std::string ReadBDString(std::ifstream &infile, int length) {
  char *stringRead = NULL;
  std::string outString = "";

  if (length != 0) {
    stringRead = new char[length];
    infile.read(stringRead, length);
    outString = stringRead;
  }
  return outString;
}

std::vector<int> ReadBDInteger(std::ifstream &infile, int length) {
  std::vector<int> retVal;

  if (length > 0) {
    retVal.resize(length);
    infile.read((char *)retVal.data(), length * sizeof(int));
  }
  return retVal;
}

std::vector<double> ReadBDNumeric(std::ifstream &infile, int length) {
  std::vector<double> retVal;

  if (length > 0) {
    retVal.resize(length);
    infile.read((char *)retVal.data(), length * sizeof(double));
  }
  return retVal;
}

Rcpp::List ReadBDSubjects(std::ifstream &infile) {
  int sidsize;
  int fidsize;
  std::string sidString;
  std::string fidString;

  infile.read((char *)&sidsize, sizeof(int));
  infile.read((char *)&fidsize, sizeof(int));

  sidString = ReadBDString(infile, sidsize);
  fidString = ReadBDString(infile, fidsize);

  return Rcpp::List::create(Rcpp::Named("sidsize") = sidsize,
                            Rcpp::Named("fidsize") = fidsize,
                            Rcpp::Named("sidstring") = sidString,
                            Rcpp::Named("fidstring") = fidString);
}

Rcpp::List ReadBDSNPs(std::ifstream &infile, int numSNPs, int numGroups, int snpoptions) {
  int snpSize, chrSize, refSize, altSize;
  std::string snpString, chrString, refString, altString;
  std::vector<int> location;
  int aafoffset, mafoffset, avgcalloffset, rsqoffset;
  std::vector<double> aaf, maf, avgCall, rsq;

  infile.read((char *)&snpSize, sizeof(int));
  infile.read((char *)&chrSize, sizeof(int));
  infile.read((char *)&refSize, sizeof(int));
  infile.read((char *)&altSize, sizeof(int));

  snpString = ReadBDString(infile, snpSize);
  chrString = ReadBDString(infile, chrSize);
  if ((snpoptions & 0x0010) != 0)
    location = ReadBDInteger(infile, numSNPs);
  refString = ReadBDString(infile, refSize);
  altString = ReadBDString(infile, altSize);

  aafoffset = 0;
  if ((snpoptions & 0x0080) != 0) {
    aafoffset = infile.tellg();
    aaf = ReadBDNumeric(infile, numSNPs * numGroups);
  }
  mafoffset = 0;
  if ((snpoptions & 0x0100) != 0) {
    mafoffset = infile.tellg();
    maf = ReadBDNumeric(infile, numSNPs * numGroups);
  }
  avgcalloffset = 0;
  if ((snpoptions & 0x0200) != 0) {
    avgcalloffset = infile.tellg();
    avgCall = ReadBDNumeric(infile, numSNPs * numGroups);
  }
  rsqoffset = 0;
  if ((snpoptions & 0x0400) != 0) {
    rsqoffset = infile.tellg();
    rsq = ReadBDNumeric(infile, numSNPs * numGroups);
  }

  return Rcpp::List::create(Rcpp::Named("snpsize") = snpSize,
                            Rcpp::Named("chrsize") = chrSize,
                            Rcpp::Named("refsize") = refSize,
                            Rcpp::Named("altsize") = altSize,
                            Rcpp::Named("snpstring") = snpString,
                            Rcpp::Named("chrstring") = chrString,
                            Rcpp::Named("location") = location,
                            Rcpp::Named("refstring") = refString,
                            Rcpp::Named("altstring") = altString,
                            Rcpp::Named("aafoffset") = aafoffset,
                            Rcpp::Named("mafoffset") = mafoffset,
                            Rcpp::Named("avgcallfoffset") = avgcalloffset,
                            Rcpp::Named("rsqoffset") = rsqoffset,
                            Rcpp::Named("aaf") = aaf,
                            Rcpp::Named("maf") = maf,
                            Rcpp::Named("avgcall") = avgCall,
                            Rcpp::Named("rsq") = rsq);
}

// Writes the additional header info for formats 4.1 and 4.2
// Parameter filename - Name of binary dosage file
// Return - list of values in header
// [[Rcpp::export]]
Rcpp::List ReadBinaryDosageHeader4A(std::string &filename) {
  std::ifstream infile;
  int numSubjects;
  int numSNPs;
  int numGroups;
  int subjectOptions;
  int snpOptions;
  int subjectOffset;
  int snpOffset;
  int dosageOffset;
  std::vector<int> groupSizes;
  Rcpp::List samples;
  Rcpp::List snps;
  Rcpp::List retVal;

  // Open the file for appending
  // Only opens for output
  infile.open(filename.c_str(), READBINARY);

  infile.seekg(8);
  infile.read((char *)&numSubjects, sizeof(int));
  infile.read((char *)&numSNPs, sizeof(int));
  infile.read((char *)&numGroups, sizeof(int));
  infile.read((char *)&subjectOptions, sizeof(int));
  infile.read((char *)&snpOptions, sizeof(int));
  infile.read((char *)&subjectOffset, sizeof(int));
  infile.read((char *)&snpOffset, sizeof(int));
  infile.read((char *)&dosageOffset, sizeof(int));

  groupSizes.resize(numGroups);
  infile.read((char *)groupSizes.data(), numGroups * sizeof(int));

  samples = ReadBDSubjects(infile);
  snps = ReadBDSNPs(infile, numSNPs, numGroups, snpOptions);

  infile.close();

  return Rcpp::List::create(Rcpp::Named("numsub") = numSubjects,
                            Rcpp::Named("numSNPs") = numSNPs,
                            Rcpp::Named("numgroups") = numGroups,
                            Rcpp::Named("suboptions") = subjectOptions,
                            Rcpp::Named("snpoptions") = snpOptions,
                            Rcpp::Named("subjectoffset") = subjectOffset,
                            Rcpp::Named("snpoffset") = snpOffset,
                            Rcpp::Named("dosageoffset") = dosageOffset,
                            Rcpp::Named("groups") = groupSizes,
                            Rcpp::Named("samples") = samples,
                            Rcpp::Named("snps") = snps);
}


// Writes the additional header info for formats 4.2 and 4.3
// Parameter filename - Name of binary dosage file
// Return - 0 successful, 1 error
// [[Rcpp::export]]
Rcpp::List ReadBinaryDosageHeader4B (std::string &filename) {
  std::ifstream infile;
  int subjectOffset;
  int snpOffset;
  int indexOffset;
  int dosageOffset;
  int numGroups;
  std::vector<int> groupSizes;
  int snpOptions;
  int numSub;
  Rcpp::List samples;
  int numSNPs;
  Rcpp::List snps;
  Rcpp::List retVal;

  // Open the file for appending
  // Only opens for output
  infile.open(filename.c_str(), READBINARY);

  infile.seekg(8);
  infile.read((char *)&subjectOffset, sizeof(int));
  infile.read((char *)&snpOffset, sizeof(int));
  infile.read((char *)&indexOffset, sizeof(int));
  infile.read((char *)&dosageOffset, sizeof(int));

  infile.read((char *)&numGroups, sizeof(int));
  groupSizes.resize(numGroups);
  infile.read((char *)groupSizes.data(), numGroups * sizeof(int));

  infile.read((char *)&numSub, sizeof(int));
  samples = ReadBDSubjects(infile);

  infile.read((char *)&numSNPs, sizeof(int));
  infile.read((char *)&snpOptions, sizeof(int));
  snps = ReadBDSNPs(infile, numSNPs, numGroups, snpOptions);

  infile.close();

  return Rcpp::List::create(Rcpp::Named("suboffset") = subjectOffset,
                            Rcpp::Named("snpoffset") = snpOffset,
                            Rcpp::Named("indexoffset") = indexOffset,
                            Rcpp::Named("dosageoffset") = dosageOffset,
                            Rcpp::Named("numgroups") = numGroups,
                            Rcpp::Named("groups") = groupSizes,
                            Rcpp::Named("numsub") = numSub,
                            Rcpp::Named("samples") = samples,
                            Rcpp::Named("numSNPs") = numSNPs,
                            Rcpp::Named("snpoptions") = snpOptions,
                            Rcpp::Named("snps") = snps);
}

// [[Rcpp::export]]
Rcpp::List ReadBDIndices3C(std::string filename,
                                    int numSNPs,
                                    int indexStart) {
  std::ifstream infile;
  Rcpp::IntegerVector datasize(numSNPs);
  Rcpp::NumericVector indices(numSNPs);
  int ds;
  Rcpp::List retval;

  infile.open(filename.c_str(), READBINARY);

  infile.seekg(indexStart);
  for (int i = 0; i < numSNPs; ++i) {
    infile.read((char *)&ds, sizeof(int));
    datasize[i] = ds;
    indices[i] = infile.tellg();
    infile.seekg(ds, std::ios_base::cur);
  }

  infile.close();

  return Rcpp::List::create(Rcpp::Named("datasize") = datasize,
                            Rcpp::Named("indices") = indices);
}

// [[Rcpp::export]]
Rcpp::List ReadBDIndices4C(std::string filename,
                                    int numSNPs,
                                    int headersize) {
  int indexstart;
  std::ifstream infile;
  Rcpp::IntegerVector datasize(numSNPs);
  Rcpp::NumericVector indices(numSNPs);
  Rcpp::List retval;

  infile.open(filename.c_str(), READBINARY);

  indexstart = headersize - numSNPs * sizeof(int);
  infile.seekg(indexstart);
  infile.read((char *)&datasize[0], numSNPs * sizeof(int));

  infile.close();

  indices[0] = headersize;
  for (int i = 1; i < numSNPs; ++i)
    indices[i] = indices[i - 1] + datasize[i - 1];

  return Rcpp::List::create(Rcpp::Named("datasize") = datasize,
                            Rcpp::Named("indices") = indices);
}

// Routine to convert short values to double
// Parameter us - is the vector of shorts to convert
// Parameter x - vector of unsigned doubles to store converted values
// Parameter base - Index of USBASE to use as base
// NOTE: missing values are coded as 0xffff for shorts
void UShortToDouble(Rcpp::IntegerVector &us,
                    Rcpp::NumericVector &x,
                    const int numsub,
                    const int base) {
  unsigned short *ps1;
  double *d;
  int i;

  ps1 = (unsigned short *)&us[0];
  d = (double *)&x[0];
  for (i = 0; i < numsub; ++i, ++ps1, ++d) {
    if (*ps1 == 0xffff) {
      // Missing
      *d = NA_REAL;
    } else {
      *d = DBASE[base] * *ps1;
    }
  }
}

// [[Rcpp::export]]
int ReadBinaryDosageDataC(std::string &filename,
                          int headersize,
                          int numsub,
                          int snp,
                          Rcpp::NumericVector &dosage,
                          Rcpp::IntegerVector &us,
                          int base) {
  std::ifstream infile;
  std::streampos loc;

  infile.open(filename.c_str(), READBINARY);

  loc = headersize + 2 * (snp - 1) * numsub;
  infile.seekg(loc);
  infile.read((char *)&us[0], numsub * sizeof(short));
  UShortToDouble(us, dosage, numsub, base - 1);
  infile.close();
  return 0;
}

// [[Rcpp::export]]
int ReadBinaryDosageDataP1P2(std::string &filename,
                             int headersize,
                             int numsub,
                             int snp,
                             Rcpp::NumericVector &dosage,
                             Rcpp::NumericVector &p0,
                             Rcpp::NumericVector &p1,
                             Rcpp::NumericVector &p2,
                             Rcpp::IntegerVector &us,
                             int base) {
  std::ifstream infile;
  std::streampos loc;
  int readsize;

  infile.open(filename.c_str(), READBINARY);

  loc = headersize + 4 * (snp - 1) * numsub;
  readsize = numsub * sizeof(short);

  infile.seekg(loc);
  infile.read((char *)&us[0], readsize);
  UShortToDouble(us, p1, numsub, base - 1);
  infile.read((char *)&us[0], readsize);
  UShortToDouble(us, p2, numsub, base - 1);
  dosage = p1 + p2 + p2;
  p0 = 1. - p1 - p2;
  for (int i = 0; i < numsub; ++i) {
    if (dosage[i] > 2.)
      dosage[i] = 2.;
    if (p0[i] < 0.)
      p0[i] = 0.;
  }
  infile.close();
  return 0;
}

// [[Rcpp::export]]
int ReadBinaryDosageDataCompressed(std::string &filename,
                                   double index,
                                   double datasize,
                                   int numsub,
                                   Rcpp::NumericVector &dosage,
                                   Rcpp::NumericVector &p0,
                                   Rcpp::NumericVector &p1,
                                   Rcpp::NumericVector &p2,
                                   Rcpp::IntegerVector &us) {
  unsigned short *usbase, *usadd;
  std::ifstream infile;

  infile.open(filename.c_str(), READBINARY);

  usbase = (unsigned short *)&us[0];
  usadd = usbase + numsub;

  infile.seekg(index);
  infile.read((char *)usbase, datasize);
  for (int i = 0; i < numsub; ++i, ++usbase) {
    if (*usbase == 0xffff) {
      dosage[i] = NA_REAL;
      p0[i] = NA_REAL;
      p1[i] = NA_REAL;
      p2[i] = NA_REAL;
      continue;
    }
    if (*usbase & 0x8000) {
      dosage[i] = (*usbase & 0x7fff) / 10000.;
      if (*usadd == 0xffff) {
        p0[i] = NA_REAL;
        p1[i] = NA_REAL;
        p2[i] = NA_REAL;
        ++usadd;
      } else if (*usadd & 0x8000) {
        p1[i] = (*usadd & 0x7fff) / 10000.;
        ++usadd;
        p0[i] = *usadd / 10000.;
        ++usadd;
        p2[i] = *usadd / 10000.;
        ++usadd;
      } else {
        p1[i] = *usadd / 10000.;
        ++usadd;
        p2[i] = (dosage[i] - p1[i]) / 2.;
        p0[i] = 1. - p2[i] - p1[i];
      }
    } else {
      dosage[i] = *usbase / 10000.;
      if (dosage[i] > 1.) {
        p2[i] = dosage[i] - 1.;
        p1[i] = dosage[i] - p2[i] - p2[i];
        p0[i] = 0.;
      } else {
        p0[i] = 1. - dosage[i];
        p1[i] = dosage[i];
        p2[i] = 0.;
      }
    }
  }

  infile.close();
  return 0;
}
