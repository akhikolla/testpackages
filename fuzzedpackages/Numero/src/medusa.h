/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#ifndef medusa_INCLUDED
#define medusa_INCLUDED

#include <string>
#include <vector>

namespace medusa {

  /* Set precision. */
  typedef unsigned char mdbyte; /* small unsigned integer */
  typedef unsigned int mdsize;  /* unsigned integer */
  typedef double mdreal;  /* floating point */

  /* Defined constants. */
  extern mdbyte bmin();
  extern mdbyte bmax();
  extern mdbyte bnan();
  extern mdreal rlim();
  extern mdreal rnan();
  extern mdsize snan();
  extern long lnan();

  /*
   * File object for text spreadsheets. 
   */
  class File {
  private:
    void* buffer;
  public:
    File();
    File(const File&); /* disabled */
    void operator=(const File&); /* disabled */
    ~File();

    /* Return file name if an active stream. */
    std::string active() const;
    
    /* Close the file stream. */
    void close();
    
    /* Non-empty if the last call to a member function
       resulted in an error. */
    std::string error() const;

    /* Return a message on the current file name and how many bytes
       have been read and written since the file was opened. */
    std::string info() const;

    /* Move the access position. A negative value moves the cursor
       towards the beginning, and a positive value towards the end. */
    bool jump(const long);
    
    /* Open a new file stream. The second argument sets the mode
       of activity, please see fopen() in <cstdio> for details. */
    bool open(const std::string&, const std::string&);

    /* Return the current access position in the file, or a negative
       value if position cannot be established. */
    long position() const;

    /* Read bytes from the file stream until a newline character. */
    std::string read();

    /* Read bytes from the file stream until a newline character.
       The line is split at bytes matching the first input. The
       second input sets the minimum number of elements in the
       output (filled with empty strings if necessary). */
    std::vector<std::string> read(const char, const mdsize);

    /* Return the number of read/written bytes. */
    unsigned long size() const;

    /* Write a single string in the file stream. Returns the
       number of bytes written. */
    unsigned long write(const std::string&); 

    /* Write a sequence of delimited strings in the file stream.
       Returns the number of bytes written. */
    unsigned long write(const std::vector<std::string>&, const char); 

    /* Check if a file pointed by the path can be opened.
       If failed, returns an error message. */
    static std::string check(const std::string&, const std::string&);
  };

  /* Buffered file functions. */
  extern FILE* openfile(const std::string&, const std::string&);
  extern bool closefile(FILE*);

  /* Utility function to split lines into strings. */
  extern std::vector<std::string> splitstr(const std::string&, const char);
  
  /*
   * Data structure to hold sparse spreadsheets of strings.
   */
  class Table {
  public:
    void* buffer;
  public:
    Table();
    Table(const Table&);
    void operator=(const Table&);
    ~Table();

    /* Insert or replace a value at the cell indicated by the
       row (1st input) and column rank (2nd input). Returns true
       if the input was stored in the cell. */
    bool insert(const mdsize, const mdsize, const std::string&);

    /* Return and erase table cell contents. */
    std::string remove(const mdsize, const mdsize);

    /* Return the values on a row. Only elements up to the last
       non-empty value are returned. */
    std::vector<std::string> row(const mdsize) const;

    /* Return the value in a specific cell. */
    std::string value(const mdsize, const mdsize) const;
  };

  /* Output format for binary search. The bounds contain the positions
     of the smaller and larger adjacent element (negative if none), and
     the weights contain the interpolation coefficients for the precise
     search key if within two existing positions. */
  struct Site {
    unsigned short usable; /* num of adjacent positions */
    std::pair<mdsize, mdsize> bounds;
    std::pair<mdreal, mdreal> weights;
  };

  /* Sorting and searching functions. */
  extern Site binsearch(const std::vector<mdreal>&, const mdreal);
  extern std::vector<mdsize> sortreal(std::vector<mdreal>&, const int);
  extern std::vector<mdsize> sortsize(std::vector<mdsize>&, const int);
  extern std::vector<mdsize> sortstr(std::vector<std::string>&, const int);

  /* Set functions. */
  extern std::vector<mdsize> uniqreal(std::vector<mdreal>&);
  extern std::vector<mdsize> uniqsize(std::vector<mdsize>&);
  extern std::vector<mdsize> uniqstr(std::vector<std::string>&);
  extern std::vector<mdsize> combine(const std::vector<mdsize>&,
				     const std::vector<mdsize>&,
				     const int);
  extern medusa::mdsize match(std::vector<medusa::mdsize>&,
			      std::vector<medusa::mdsize>&,
			      const std::vector<std::string>&,
			      const std::vector<std::string>&);

  /* Conversion functions. */
  extern long string2long(const std::string&);
  extern mdsize string2size(const std::string&);
  extern mdreal string2real(const std::string&);
  extern std::string string2safe(const std::string&, const mdsize);
  extern std::string long2string(const long);
  extern std::string long2text(const long);
  extern std::string real2string(const mdreal);
  extern std::string real2text(const mdreal);
  extern std::string time2text(const mdreal);
 
  /* Notifications. */
  extern void panic(const std::string&, const char*, const int);
  extern void worry(const std::string&, const char*);

  /* Time and date functions. */
  extern std::string currtime();

  /* Version information. */
  extern std::string version();
}

#endif /* medusa_INCLUDED */
