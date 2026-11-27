#ifndef PROJECT_PATHS_HPP
#define PROJECT_PATHS_HPP

#include <TString.h>
#include <TSystem.h>

namespace Paths {
// The 73Ge project root, from the path of a file one directory below it
// (czt/, cdte/): pass the macro's __FILE__.
inline TString ProjectRootOf(const char *file) {
  TString path = file;
  path.ReplaceAll("/./", "/");
  TString macros_dir = gSystem->DirName(path);
  return gSystem->DirName(macros_dir);
}
} // namespace Paths

#endif
