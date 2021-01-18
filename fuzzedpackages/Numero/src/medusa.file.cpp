/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "medusa.local.h"

/*
 *
 */
File::File() {
  this->buffer = new FileBuffer();
}

/*
 *
 */
File::File(const File& t) {
  panic("Copy constructor not available.", __FILE__, __LINE__);
}

/*
 *
 */
void
File::operator=(const File& t) {
  panic("Copy operator not available.", __FILE__, __LINE__);
}

/*
 *
 */
File::~File() {
  FileBuffer* p = (FileBuffer*)buffer;
  delete p;
}

/*
 *
 */
string
File::active() const {
  FileBuffer* p = (FileBuffer*)buffer;
  if(p->fid != NULL) return p->name;
  return "";
}

/*
 *
 */
void
File::close() {
  FileBuffer* p = (FileBuffer*)buffer; delete p;
  this->buffer = new FileBuffer();
}

/*
 *
 */
string
File::error() const {
  FileBuffer* p = (FileBuffer*)buffer;
  return p->errtext;
}

/*
 *
 */
string
File::info() const {
  FileBuffer* p = (FileBuffer*)buffer;
  string fn = p->name;
  mdsize len = fn.size();
  if(len > 32) fn = ("[..]" + fn.substr(len - 28));
  if(p->nread > 0) return (fn + " -> " + long2text(p->nread) + " bytes");
  if(p->nwritten > 0) return (long2text(p->nwritten) + " bytes -> " + fn);
  return fn;
}

/*
 *
 */
bool
File::jump(const long offset) {
  FileBuffer* p = (FileBuffer*)buffer;
  if(p->fid == NULL) {
    p->abort("No file stream.");
    return false;
  }
  return (fseek(p->fid, offset, SEEK_CUR) == 0);
}

/*
 *
 */
bool
File::open(const string& fname, const string& prm) {
  FileBuffer* p = (FileBuffer*)buffer;

  /* Discard previous contents. */
  p->clear();

  /* Open a new file stream. */
  p->fid = fopen(fname.c_str(), prm.c_str());
  p->name = fname;

  /* Set I/O buffer. */
  if(p->fid == NULL) p->abort("Could not open '" + fname + "'.");
  else setvbuf(p->fid, p->iobuf, _IOFBF, IOBUFCAP_medusa);
  return (p->fid != NULL);
}

/*
 *
 */
long
File::position() const {
  FileBuffer* p = (FileBuffer*)buffer;
  if(p->fid == NULL) {
    p->abort("No file stream.");
    return -2;
  }
  return ftell(p->fid);
}

/*
 *
 */
string
File::read() {
  FileBuffer* p = (FileBuffer*)buffer;

  /* Check for previous errors. */
  if((p->errtext).size() > 0) return "";
  if(p->fid == NULL) return "";

  /* Read a line from the file stream. */
  char* data = fgets(p->bytes, IOBUFCAP_medusa, p->fid);
  if(data == NULL) {
    p->abort("No data.");
    return "";
  }

  /* Check that line was fully read. */
  mdsize ndone = strlen(data);
  p->nread += ndone;
  if(ndone >= IOBUFCAP_medusa) {
    p->abort("Line data exceeded buffer capacity.");
    return "";
  }
  
  /* Remove trailing line breaks. */
  while(ndone > 0) {
    if((data[ndone] != '\r') &&
       (data[ndone] != '\n')) break;
    data[ndone] = '\0';
    ndone--;
  }
  return string(data);
}

/*
 *
 */
vector<string>
File::read(const char delim, const mdsize nmin) {
  FileBuffer* p = (FileBuffer*)buffer;
  vector<string> fields(nmin);

  /* Check for previous errors. */
  if((p->errtext).size() > 0) return fields;
  if(p->fid == NULL) return fields;

  /* Read a line from the file stream. */
  char* data = fgets(p->bytes, IOBUFCAP_medusa, p->fid);
  if(data == NULL) {
    p->abort("No data.");
    return fields;
  }

  /* Replace carriage returns and delimiters. */
  mdsize ndone = 0;
  mdsize nbytes = 0;
  bool spaceprev = false;
  bool spaceflag = (delim == '\0');
  for(mdsize i = 0; data[i] != '\0'; i++, ndone++) {
    char c = data[i];
    if(c == '\r') continue;
    if(spaceflag && isspace(c)) {
      if(spaceprev) continue;
      c = '\0';
    }
    spaceprev = (c == delim);
    if(spaceprev) c = '\0';
    data[nbytes] = c;
    nbytes++;
  }

  /* Check that line was fully read. */
  p->nread += ndone;
  if(ndone >= IOBUFCAP_medusa) {
    p->abort("Line data exceeded buffer capacity.");
    return fields;
  }

  /* Ensure correct termination. */
  if(data[nbytes-1] == '\n') nbytes--;
  if(nbytes < 1) return fields;
  data[nbytes] = '\0';

  /* Collect string segments. */
  char* ptr = data;
  mdsize nfields = 0;
  for(mdsize i = 0; i <= nbytes; i++) {
    if(data[i] != '\0') continue;
    if(nfields < nmin) fields[nfields] = string(ptr);
    else fields.push_back(string(ptr));
    ptr = (char*)(data + i + 1);
    nfields++;
  }
  return fields;
}

/*
 *
 */
unsigned long
File::size() const {
  FileBuffer* p = (FileBuffer*)buffer;
  return (p->nread + p->nwritten);
}

/*
 *
 */
unsigned long
File::write(const string& s) {
  FileBuffer* p = (FileBuffer*)buffer;
  if(p->fid == NULL) {
    p->abort("No file stream.");
    return 0;
  }
  unsigned long n = fprintf(p->fid, "%s", s.c_str());
  if(n < s.size()) p->abort("Write failed.");
  p->nwritten += n;
  return n;
}

/*
 *
 */
unsigned long
File::write(const vector<string>& fields, const char delim) {
  FileBuffer* p = (FileBuffer*)buffer;
  unsigned long n = 0;

  /* Check for errors. */
  if(fields.size() < 1) return 0;
  if(p->fid == NULL) {
    p->abort("No file stream.");
    return 0;
  }

  /* Print fields. */
  unsigned long nbytes = fields[0].size();
  n += fprintf(p->fid, "%s", fields[0].c_str());
  for(mdsize j = 1; j < fields.size(); j++) {
    n += fprintf(p->fid, "%c%s", delim, fields[j].c_str());
    nbytes += (fields[j].size() + 1);
  }

  /* Check byte count. */
  n += fprintf(p->fid, "\n"); nbytes++;
  if(n < nbytes) p->abort("Write failed.");
  p->nwritten += n;
  return n;
}

/*
 *
 */
string
File::check(const string& fname, const string& prm) {
  if(fname == "") return "Empty file name.";
  FILE* fid = fopen(fname.c_str(), prm.c_str());
  if(fid == NULL) return ("File '" + fname + "' is inaccessible.");
  fclose(fid);
  return "";
}
