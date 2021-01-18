#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include <vector>

//***************************************************************************//
//                                                                           //
//                        Writing the header                                 //
//                                                                           //
//***************************************************************************//
// These functions write the headers to the binary dosage files
// NOTE: There is no error checking done here. It is assumed this was done
// prior to calling these routines

//***************************************************************************//
//                           Constants                                       //
//***************************************************************************//

// Magic word for binary dosage files
extern const int MAGICWORD = 0x65736f62;
// Format ID stored in file
extern const std::vector<std::vector<int> > FORMAT = {
  { 0x01000100, 0x02000100},
  { 0x01000200, 0x02000200},
  { 0x01000300, 0x02000300, 0x03000300, 0x04000300},
  { 0x01000400, 0x02000400, 0x03000400, 0x04000400}
};
// Various modes of opening the binary dosage file
extern const std::ios_base::openmode NEWBINARY = std::ios_base::out | std::ios_base::binary;
extern const std::ios_base::openmode READBINARY = std::ios_base::in | std::ios_base::binary;
extern const std::ios_base::openmode APPENDBINARY = std::ios_base::out | std::ios_base::binary | std::ios_base::app;
extern const std::ios_base::openmode READWRITEBINARY = std::ios_base::in | std::ios_base::out | std::ios_base::binary | std::ios_base::ate;
// zero values for filling initial values of binary dosage file to 0
const char CHARZERO = 0x0;
const int INTZERO = 0;
//const double DOUBLEZERO = 0.;

// Tolerance for writing data
const double WRITINGTOLERANCE = 1e-6;

struct OFFSETS {
  enum offsets { numsubjects, numsnps, numgroups, suboptions, snpoptions,
                 subjectoffset, snpoffset, indexoffset, dosageoffset};
};

//***************************************************************************//
//                      Support function                                     //
//***************************************************************************//

// ******************** Writing vectors *************************************//

// Write a NULL terminated string to the header
int WriteBDString(std::fstream &outfile, std::string &outstring) {
  if (outstring.length() > 0) {
    outfile.write(outstring.c_str(), outstring.length());
    outfile.write(&CHARZERO, 1);
  }
  return 0;
}

// Write a integer vector to the header
int WriteBDInteger(std::fstream &outfile, Rcpp::IntegerVector &outvector) {
  if (outvector.length() > 0)
    outfile.write((char *)&outvector[0], outvector.length() * sizeof(int));
  return 0;
}

// Write a numeric vector
int WriteBDNumeric(std::fstream &outfile, Rcpp::NumericVector &outvector) {
  if (outvector.length() > 0)
    outfile.write((char *)&outvector[0], outvector.length() * sizeof(double));
  return 0;
}

//  ********** Writing addition information for formats 4.x *****************//

// Write the group size information
int WriteBDGroups(std::fstream &outfile,
                  Rcpp::IntegerVector &groups,
                  int numGroupLoc,
                  int subjectOffsetLoc) {
  int numGroups;
  int subjectOffset;

  numGroups = groups.size();
  if (numGroupLoc >= 0)
    outfile.seekp(numGroupLoc);
  outfile.write((char *)&numGroups, sizeof(int));

  outfile.seekp(0, std::ios_base::end);
  outfile.write((char *)&groups[0], numGroups * sizeof(int));

  subjectOffset = (int)outfile.tellp();
  outfile.seekp(subjectOffsetLoc);
  outfile.write((char *)&subjectOffset, sizeof(int));

  return 0;
}

// Write a the family information
int WriteBDFamilyInfo(std::fstream &outfile,
                      int numSub,
                      std::string &sid,
                      std::string &fid,
                      int numSubLoc,
                      int suboffsetLoc,
                      int snpoffsetLoc) {
  int sidsize, fidsize;
  int suboffset, snpoffset;

  sidsize = sid.length() + 1;
  fidsize = fid.length();
  if (fidsize > 0)
    ++fidsize;

  outfile.seekg(suboffsetLoc);
  outfile.read((char *)&suboffset, sizeof(int));
  if (numSubLoc < 0) {
    outfile.seekp(suboffset);
    outfile.write((char *)&numSub, sizeof(int));
  } else {
    outfile.seekp(numSubLoc);
    outfile.write((char *)&numSub, sizeof(int));
    outfile.seekp(suboffset);
  }

  outfile.write((char *)&sidsize, sizeof(int));
  outfile.write((char *)&fidsize, sizeof(int));
  WriteBDString(outfile, sid);
  WriteBDString(outfile, fid);

  snpoffset = outfile.tellp();
  outfile.seekp(snpoffsetLoc);
  outfile.write((char *)&snpoffset, sizeof(int));

  return 0;
}

// Write a the SNP information for format 4
int WriteBDSNPInfo(std::fstream &outfile,
                  int numSNPs,
                  std::string &snpid,
                  std::string &chromosome,
                  Rcpp::IntegerVector &location,
                  std::string &reference,
                  std::string &alternate,
                  Rcpp::NumericVector &aaf,
                  Rcpp::NumericVector &maf,
                  Rcpp::NumericVector &avgCall,
                  Rcpp::NumericVector &rsq,
                  int numSNPloc,
                  int snpOptionsLoc,
                  int snpOffsetLoc,
                  int nextOffsetLoc) {
  int snpoffset, nextoffset;
  int stringSize[4];
  int snpOptions = 0;

  if (snpid.length() != 0)
    snpOptions |= 0x02;
  if (chromosome.size() !=  0) {
    snpOptions |= 0x04;
    if (chromosome.find('\t') == std::string::npos)
      snpOptions |= 0x08;
  }
  if (location.length() != 0)
    snpOptions |= 0x10;
  if (reference.length() != 0)
    snpOptions |= 0x20;
  if (alternate.length() != 0)
    snpOptions |= 0x40;
  if (aaf.length() != 0)
    snpOptions |= 0x80;
  if (maf.length() != 0)
    snpOptions |= 0x100;
  if (avgCall.length() != 0)
    snpOptions |= 0x200;
  if (rsq.length() != 0)
    snpOptions |= 0x400;

  outfile.seekg(snpOffsetLoc);
  outfile.read((char *)&snpoffset, sizeof(int));

  if (snpOptionsLoc < 0) {
    outfile.seekp(snpoffset);
    outfile.write((char *)&numSNPs, sizeof(int));
    outfile.write((char *)&snpOptions, sizeof(int));
  } else {
    outfile.seekp(numSNPloc);
    outfile.write((char *)&numSNPs, sizeof(int));
    outfile.seekp(snpOptionsLoc);
    outfile.write((char *)&snpOptions, sizeof(int));
    outfile.seekp(snpoffset);
  }

  stringSize[0] = snpid.length() == 0 ? 0 : snpid.length() + 1;
  stringSize[1] = chromosome.length() == 0 ? 0 : chromosome.length() + 1;
  stringSize[2] = reference.length() == 0 ? 0 : reference.length() + 1;
  stringSize[3] = alternate.length()  == 0 ? 0 : alternate.length() + 1;
  outfile.write((char *)&stringSize[0], sizeof(stringSize));

  WriteBDString(outfile, snpid);
  WriteBDString(outfile, chromosome);
  WriteBDInteger(outfile, location);
  WriteBDString(outfile, reference);
  WriteBDString(outfile, alternate);
  WriteBDNumeric(outfile, aaf);
  WriteBDNumeric(outfile, maf);
  WriteBDNumeric(outfile, avgCall);
  WriteBDNumeric(outfile, rsq);

  nextoffset = outfile.tellp();
  outfile.seekp(nextOffsetLoc);
  outfile.write((char *)&nextoffset, sizeof(int));

  return 0;
}

// Write the initial indices values of 0
int WriteBDIndices(std::fstream &outfile, int numIndices,
                   int indexOffsetLoc, int dosageOffsetLoc) {
  int indexoffset, dosageoffset;

  outfile.seekg(indexOffsetLoc);
  outfile.read((char *)&indexoffset, sizeof(int));
  outfile.seekp(indexoffset);
  for (int i = 0; i < numIndices; ++i)
    outfile.write((char *)&INTZERO, sizeof(int));
  dosageoffset = outfile.tellp();
  outfile.seekp(dosageOffsetLoc);
  outfile.write((char *)&dosageoffset, sizeof(int));
  return 0;
}

//***************************************************************************//
//                 Open the binary dosage file                               //
//***************************************************************************//

// Open a new file for writing
// If file already exists, data is lost
int OpenBDFileNewWrite(std::ofstream &outfile, std::string &filename) {
  outfile.open(filename.c_str(), NEWBINARY);
  if (!outfile.good())
    return 1;
  return 0;
}

// Open a file for appending
// The file is opened for writing with the initial writing pointer
// set to the end of the file.
int OpenBDFileAppend(std::ofstream &outfile, std::string &filename) {
  outfile.open(filename.c_str(), APPENDBINARY);
  if (!outfile.good())
    return 1;
  return 0;
}

// Open a file for reading and writing
int OpenBDFileReadWrite(std::fstream &outfile, std::string &filename) {
  outfile.open(filename.c_str(), READWRITEBINARY);
  if (!outfile.good())
    return 1;
  return 0;
}

//***************************************************************************//
//                 Writing the header                                        //
//***************************************************************************//

//  ************** Write the base header for all formats*********************//

// Writes the base header for a binary dosage file
// This is the complete header for formats 1.x and 2.x
// Parameter filename - Name of binary dosage file
// Parameter format - Foramt of the binary dosage file
// Parameter subformat - Subformat of the binary dosage file
// Return - 0 successful, 1 error
// [[Rcpp::export]]
int WriteBinaryDosageBaseHeader(std::string &filename, int format, int subformat) {
  std::ofstream outfile;

  if (OpenBDFileNewWrite(outfile, filename) != 0)
    Rcpp::stop("Unable to create output file");

  outfile.write((char *)&MAGICWORD, sizeof(int));
  outfile.write((char *)&FORMAT[format][subformat], sizeof(int));

  outfile.close();
  return 0;
}

// Write the header for formats 3.1 and 3.2
// Parameter filename - Name of binary dosage file
// Parameter numSubjects - number of subjects in data
// Return - 0 successful, 1 error
// [[Rcpp::export]]
int WriteBinaryDosageHeader3A(std::string &filename,
                              int numSubjects) {
  std::ofstream outfile;

  // Open the file for appending
  if (OpenBDFileAppend(outfile, filename) != 0)
    Rcpp::stop("Unable to open file for appending");

  outfile.write((char *)&numSubjects, sizeof(int));

  outfile.close();
  return 0;
}

// Write header for formats 3.3 and 3.4
// Parameter filename - Name of binary dosage file
// Parameter md5samples - MD5 hash for the samples data frame
// Parameter md5SNPs - MD5 hash for the SNPs data frame
// Paramter numIndices - Number of indices to write
// Return - 0 successful, 1 error
// [[Rcpp::export]]
int WriteBinaryDosageHeader3B(std::string &filename,
                              std::string &md5samples,
                              std::string &md5SNPs,
                              int numIndices) {
  std::ofstream outfile;

  // Open the file for appending
  if (OpenBDFileAppend(outfile, filename) != 0)
    Rcpp::stop("Unable to open file for appending");

  outfile.write(md5samples.c_str(), 32);
  outfile.write(md5SNPs.c_str(), 32);
  for (int i = 0; i < numIndices; ++i)
    outfile.write((char *)&INTZERO, sizeof(int));

  outfile.close();
  return 0;
}

// Write the header for formats 4.1 and 4.2
// Parameter filename - Name of binary dosage file
// Parameter numSubjects - number of subjects in data
// Parameter numSubjects - number of SNPs in data
// Return - 0 successful, 1 error
// [[Rcpp::export]]
int WriteBinaryDosageHeader4A(std::string &filename,
                              int headerEntries,
                              int numSubjects,
                              int numSNPs,
                              Rcpp::IntegerVector &groups,
                              std::string &sid,
                              std::string &fid,
                              std::string &snpid,
                              std::string &chromosome,
                              Rcpp::IntegerVector &location,
                              std::string &reference,
                              std::string &alternate,
                              Rcpp::NumericVector &aaf,
                              Rcpp::NumericVector &maf,
                              Rcpp::NumericVector &avgCall,
                              Rcpp::NumericVector &rsq,
                              Rcpp::IntegerVector &offsets,
                              int numIndices) {
  std::fstream outfile;

  if (OpenBDFileReadWrite(outfile, filename) != 0)
    Rcpp::stop("Unable to open file for read/write");
  outfile.seekp(8);

  // Zero out the rest of the data. It will be filled in later
  for (int i = 0; i < headerEntries; ++i)
    outfile.write((char *)&INTZERO, sizeof(int));

  WriteBDGroups(outfile,
                groups,
                offsets[OFFSETS::offsets::numgroups],
                offsets[OFFSETS::offsets::subjectoffset]);
  WriteBDFamilyInfo(outfile,
                    numSubjects,
                    sid,
                    fid,
                    offsets[OFFSETS::offsets::numsubjects],
                    offsets[OFFSETS::offsets::subjectoffset],
                    offsets[OFFSETS::offsets::snpoffset]);
  WriteBDSNPInfo(outfile,
                numSNPs,
                snpid,
                chromosome,
                location,
                reference,
                alternate,
                aaf,
                maf,
                avgCall,
                rsq,
                offsets[OFFSETS::offsets::numsnps],
                offsets[OFFSETS::offsets::snpoptions],
                offsets[OFFSETS::offsets::snpoffset],
                offsets[OFFSETS::offsets::indexoffset]);
  WriteBDIndices(outfile,
                 numIndices,
                 offsets[OFFSETS::offsets::indexoffset],
                 offsets[OFFSETS::offsets::dosageoffset]);
  outfile.close();
  return 0;
}

//***************************************************************************//
//                                                                           //
//                        Writing the data                                   //
//                                                                           //
//***************************************************************************//
// These functions appends data to the end of the binary dosage file
// NOTE: There is no error checking done here. It is assumed this was done
// prior to calling these routines

//  ************************ Constants **************************************//

// Values that dosages and genetic probabilities are multiplied by to change
// them to short integers

extern const int NUMBEROFBASES = 3;
// 0x7ffe is 32,767 or 2^16 - 1
// 0xfff3 is 65,534 or 2^32 - 1
// 0x2710 is 10,000
extern const unsigned short USBASE[NUMBEROFBASES] = {
  0x7ffe, // Used for format 1.1
  0xfffe, // Used for format 1.2
  0x2710  // Used for all other formats
};

// Values the short integers are multiplied by to get dosage and genetic
// probabilities
extern const double DBASE[NUMBEROFBASES] = {
  1. / USBASE[0],
  1. / USBASE[1],
  1. / USBASE[2]
};

// Routine to convert a double value to an unsigned short
// Parameter x - is the vector of doubles to convert
// Parameter us - vector of unsigned shorts to store converted values
// Parameter base - Index of USBASE to use as base
// NOTE: missing values are coded as 0xffff
unsigned short DoubleToUShort(const double x, const int base) {
  unsigned short r1, r2;
  double x1, x2;

  r1 = r2 = 0;
  x1 = x2 = 0.;

// Can never happen. Missingness was already checked.
// Left in case something changes
//  if (x != x) {
    // Missing
//    return 0xffff;
//  }
  r1 = (unsigned short)(x * USBASE[base]);
  // The following section checks if r1 or (r1 -1) or (r1 + 1)
  // gives the closest value to the double passed
  x1 = r1 * DBASE[base];
  if (r1 * DBASE[base] < x)
    r2 = r1 + 1;
  else
    r2 = r1 - 1;
  x2 = r2 * DBASE[base];

  return (fabs(x - x1) < fabs(x - x2)) ? r1 : r2;
}

// Routine to convert a vector of double values to a vector of unsigned shorts
// Parameter x - is the vector of doubles to convert
// Parameter us - vector of unsigned shorts to store converted values
// Parameter base - Index of USBASE to use as base
// NOTE: missing values are coded as 0xffff
void DoubleToUShort(Rcpp::NumericVector &x,
                    Rcpp::IntegerVector &us,
                    const int base) {
  unsigned short r1, r2;
  unsigned short *ps1;
  double x1, x2;
  int i;

  r1 = r2 = 0;
  x1 = x2 = 0.;
  ps1 = (unsigned short *)&us[0];
  for (i = 0; i < x.size(); ++i, ++ps1) {
    if (x[i] != x[i]) {
      // Missing
      *ps1 = 0xffff;
    } else {
      r1 = (unsigned short)(x[i] * USBASE[base]);
      // The following section checks if r1 or (r1 -1) or (r1 + 1)
      // gives the closest value to the double passed
      x1 = r1 * DBASE[base];
      if (r1 * DBASE[base] < x[i])
        r2 = r1 + 1;
      else
        r2 = r1 - 1;
      x2 = r2 * DBASE[base];

      *ps1 = (fabs(x[i] - x1) < fabs(x[i] - x2)) ? r1 : r2;
    }
  }
}

// Write only the dosage data
// Paramater filename - name of bindary dosage file
// Parameter dosage - vector of dosage values to write
// Parameter usdosage - vector used to store the converted dosage values. This
//                      is passed to avoid allocating and dealllocating memory
// Parameter base - Index of USBASE to use as base
// [[Rcpp::export]]
int WriteBinaryDosageDataC(std::string &filename,
                          Rcpp::NumericVector &dosage,
                          Rcpp::IntegerVector &us,
                          int base) {
  std::ofstream outfile;

  // Opens file for appending
  if (OpenBDFileAppend(outfile, filename) != 0)
    Rcpp::stop("Unable to open file for appending");

  DoubleToUShort(dosage, us, base - 1);
  outfile.write((char *)&us[0], dosage.size() * sizeof(short));

  outfile.close();
  return 0;
}

// Write only the p1 and p2 data
// Paramater filename - name of bindary dosage file
// Parameter p1 - vector of P(g=1) to write
// Parameter p2 - vector of P(g=2) to write
// Parameter us - vector used to store the converted values. This
//                is passed to avoid allocating and dealllocating memory
// Parameter base - Index of USBASE to use as base
// [[Rcpp::export]]
int WriteBinaryP1P2Data(std::string &filename,
                        Rcpp::NumericVector &p1,
                        Rcpp::NumericVector &p2,
                        Rcpp::IntegerVector &us,
                        int base) {
  std::ofstream outfile;

  // Opens file for appending
  if (OpenBDFileAppend(outfile, filename) != 0)
    Rcpp::stop("Unable to open file for appending");

  DoubleToUShort(p1, us, base - 1);
  outfile.write((char *)&us[0], p1.size() * sizeof(short));
  DoubleToUShort(p2, us, base - 1);
  outfile.write((char *)&us[0], p2.size() * sizeof(short));

  outfile.close();
  return 0;
}

// Write compressed data for formats 3 and 4
// Paramater filename - name of bindary dosage file
// Parameter dosage - vector of dosages to write
// Parameter p0 - vector of P(g=0) to write
// Parameter p1 - vector of P(g=1) to write
// Parameter p2 - vector of P(g=2) to write
// Parameter snpnumber - The 0-based number of snp being written
// Parameter datasize - vector of number of bytes used to store SNP data
// Parameter us - vector used to store the converted values. This
//                is passed to avoid allocating and dealllocating memory
// [[Rcpp::export]]
int WriteBinaryCompressed(std::string &filename,
                          Rcpp::NumericVector &dosage,
                          Rcpp::NumericVector &p0,
                          Rcpp::NumericVector &p1,
                          Rcpp::NumericVector &p2,
                          Rcpp::IntegerVector &snpnumber,
                          Rcpp::IntegerVector &datasize,
                          Rcpp::IntegerVector &us) {
  int i;
  int additionallength;
  unsigned short *usdose, *usadd;
  unsigned int writesize;
  std::ofstream outfile;

  // Opens file for appending
  if (OpenBDFileAppend(outfile, filename) != 0)
    Rcpp::stop("Unable to open file for appending");

  DoubleToUShort(dosage, us, 2);

  additionallength = 0;

  usdose = (unsigned short *)&us[0];
  usadd = usdose + dosage.size();
  for (i = 0; i < dosage.length(); ++i, ++usdose) {
    if (dosage[i] != dosage[i])
      continue;
    if (p0[i] != p0[i] || p1[i] != p1[i] || p2[i] != p2[i]) {
      *usdose |= 0x8000;
      *usadd = 0xffff;
      ++usadd;
      ++additionallength;
    } else if (fabs(p0[i] + p1[i] + p2[i] - 1.) < WRITINGTOLERANCE
                 && fabs(p1[i] + p2[i] + p2[i] - dosage[i]) < WRITINGTOLERANCE) {
      if (p0[i] != 0 && p2[i] != 0) {
        *usdose |= 0x8000;
        *usadd = DoubleToUShort(p1[i], 2);
        ++usadd;
        ++additionallength;
      }
    } else {
      *usdose |= 0x8000;
      *usadd = 0x8000 | DoubleToUShort(p1[i], 2);
      ++usadd;
      *usadd = DoubleToUShort(p0[i], 2);
      ++usadd;
      *usadd = DoubleToUShort(p2[i], 2);
      ++usadd;
      additionallength += 3;
    }
  }

  writesize = (dosage.size() + additionallength) * sizeof(short);
  if (snpnumber[0] < 0) {
    outfile.write((char *)&writesize, sizeof(int));
  } else {
    datasize[snpnumber[0]] = writesize;
    ++snpnumber[0];
  }
  outfile.write((char *)&us[0], writesize);

  outfile.close();
  return 0;
}

// [[Rcpp::export]]
int WriteBinaryDosageIndicesC(std::string &filename, int headersize, Rcpp::IntegerVector &datasize) {
  int indexloc;
  std::fstream outfile;

  if (OpenBDFileReadWrite(outfile, filename) != 0)
    Rcpp::stop("Unable to open file for read/write");

  indexloc = headersize - datasize.size() * sizeof(int);
  outfile.seekp(indexloc);
  outfile.write((char *)&datasize[0], datasize.size() * sizeof(int));
  outfile.close();
  return 0;
}

// [[Rcpp::export]]
int updatesnpinfo(std::string &filename, int offset, Rcpp::NumericVector &value) {
  std::fstream outfile;

  if (OpenBDFileReadWrite(outfile, filename) != 0)
    Rcpp::stop("Unable to open file for read/write");

  outfile.seekp(offset);
  outfile.write((char *)&value[0], value.size() * sizeof(double));
  outfile.close();
  return 0;
}
