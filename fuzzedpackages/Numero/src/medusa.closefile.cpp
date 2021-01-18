/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "medusa.local.h"

/*
 *
 */
bool
medusa::closefile(FILE* fid) {
  if((FileBuffer::handles).count(fid) < 1)
    panic("Invalid stream.\n", __FILE__, __LINE__);
  fclose(fid);
  free(FileBuffer::handles[fid]);
  (FileBuffer::handles).erase(fid);
  return true;
}
