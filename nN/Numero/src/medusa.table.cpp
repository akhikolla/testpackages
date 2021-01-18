/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "medusa.local.h"

/*
 *
 */
Table::Table() {
  this->buffer = new TableBuffer();
}

/*
 *
 */
Table::Table(const Table& t) {
  this->buffer = new TableBuffer(t.buffer);
}

/*
 *
 */
void
Table::operator=(const Table& t) {
  if(this == &t) return;
  TableBuffer* p = (TableBuffer*)buffer; delete p;
  this->buffer = new TableBuffer(t.buffer);
}

/*
 *
 */
Table::~Table() {
  TableBuffer* p = (TableBuffer*)buffer;
  delete p;
}

/*
 *
 */
bool
Table::insert(const mdsize r, const mdsize c, const string& s) {
  TableBuffer* p = (TableBuffer*)buffer;
  mdsize sznan = medusa::snan();

  /* Check inputs. */
  if(r == sznan) return false;
  if(c == sznan) return false;
  if(s.size() < 1) return false;

  /* Create a new row. */
  unordered_map<mdsize, unordered_map<mdsize, mdsize> >::iterator rpos;
  unordered_map<mdsize, unordered_map<mdsize, mdsize> >& data = p->data;
  if((rpos = data.find(r)) == data.end()) {
    data[r].clear();
    rpos = data.find(r);
  }
  
  /* Find or create the table cell. */
  unordered_map<mdsize, mdsize>::iterator cpos;
  unordered_map<mdsize, mdsize>& rowdata = rpos->second;
  if((cpos = rowdata.find(c)) == rowdata.end()) {
    rowdata[c] = sznan;
    cpos = rowdata.find(c);
  }

  /* Update an empty cell. */
  if(cpos->second == sznan) {
    cpos->second = p->grow(s);
    return true;
  }

  /* Update a cell with different contents. */
  string& prev = p->words[cpos->second];
  if(s != prev) {
    p->shrink(cpos->second);
    cpos->second = p->grow(s);
  }
  return true;
}

/*
 *
 */
string
Table::remove(const mdsize r, const mdsize c) {
  TableBuffer* p = (TableBuffer*)buffer;

  /* Find the row. */
  unordered_map<mdsize, unordered_map<mdsize, mdsize> >& data = p->data;
  unordered_map<mdsize, unordered_map<mdsize, mdsize> >::iterator rpos;
  if((rpos = data.find(r)) == data.end()) return "";

  /* Find the table cell. */
  unordered_map<mdsize, mdsize>::iterator cpos;
  unordered_map<mdsize, mdsize>& rowdata = rpos->second;
  if((cpos = rowdata.find(c)) == rowdata.end()) return "";

  /* Copy cell contents. */
  string s = p->words[cpos->second];

  /* Erase the cell. */
  rowdata.erase(c);
  if(rowdata.size() < 1) data.erase(r);
  p->shrink(cpos->second);
  return s;
}

/*
 *
 */
vector<string>
Table::row(const mdsize r) const {
  TableBuffer* p = (TableBuffer*)buffer;
  vector<string> array;

  /* Find the row. */
  unordered_map<mdsize, unordered_map<mdsize, mdsize> >& data = p->data;
  unordered_map<mdsize, unordered_map<mdsize, mdsize> >::iterator rpos;
  if((rpos = data.find(r)) == data.end()) return array;

  /* Return the row contents. */
  mdsize ndata = 0;
  unordered_map<mdsize, string>& words = p->words;
  unordered_map<mdsize, mdsize>& rowdata = rpos->second;
  for(mdsize c = 0; ndata < rowdata.size(); c++) {
    unordered_map<mdsize, mdsize>::iterator cpos;
    if((cpos = rowdata.find(c)) == rowdata.end()) continue;
    else array.resize(c);
    array.push_back(words[cpos->second]);
    ndata++;
  }
  return array;
}

/*
 *
 */
string
Table::value(const mdsize r, const mdsize c) const {
  TableBuffer* p = (TableBuffer*)buffer;

  /* Find the row. */
  unordered_map<mdsize, unordered_map<mdsize, mdsize> >& data = p->data;
  unordered_map<mdsize, unordered_map<mdsize, mdsize> >::iterator rpos;
  if((rpos = data.find(r)) == data.end()) return "";

  /* Find the table cell. */
  unordered_map<mdsize, mdsize>::iterator cpos;
  unordered_map<mdsize, mdsize>& rowdata = rpos->second;
  if((cpos = rowdata.find(c)) == rowdata.end()) return "";

  /* Return cell contents. */
  return p->words[cpos->second];
}

/*
 *
 */
mdsize
TableBuffer::grow(const string& w) {
  if(w.size() < 1) panic("Empty string.", __FILE__, __LINE__);

  /* Create a new word. */
  mdsize rank = medusa::snan();
  unordered_map<string, pair<mdsize, mdsize> >::iterator pos;
  if((pos = word2rank.find(w)) == word2rank.end()) {
    
    /* Find the next available rank. */
    mdsize down = word2rank.size();
    mdsize up = (word2rank.size() + 1);
    while(true) {
      if(words.count(down) < 1) {rank = down; break;}
      if(words.count(up) < 1) {rank = up; break;}
      if(down > 0) down--;
      up++;
    }

    /* Update data structures. */
    word2rank[w] = pair<mdsize, mdsize>(rank, 1);
    words[rank] = w;
    return rank;
  }

  /* Update existing word. */
  (pos->second).second += 1;
  return (pos->second).first;
}

/*
 *
 */
mdsize
TableBuffer::shrink(const mdsize wrank) {

  /* Check that word exists. */
  unordered_map<mdsize, string>::iterator pos;
  if((pos = words.find(wrank)) == words.end())
    panic("Unusable input.", __FILE__, __LINE__);

  /* Update counter. */
  pair<mdsize, mdsize>& entry = word2rank[pos->second];
  entry.second -= 1;

  /* Delete unused word. */
  if(entry.second < 1) {
    word2rank.erase(pos->second);
    words.erase(wrank);
    return 0;
  }
  return entry.second;
}
