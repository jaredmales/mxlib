/** \file gitRepo.hpp
  * \author Jared R. Males
  * \brief Interrogate the current state of a git repository (declarations)
  * \ingroup utils_files
  */

//***********************************************************************//
// Copyright 2021 Jared R. Males (jaredmales@gmail.com)
//
// This file is part of mxlib.
//
// mxlib is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// mxlib is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with mxlib.  If not, see <http://www.gnu.org/licenses/>.
//***********************************************************************//

#ifndef gitRepo_hpp
#define gitRepo_hpp

#include <iostream>
#include <set>

namespace mx
{
namespace sys
{
   
/// Interrogate the current state of a git repository
/** Once the target directory is set, either on construction
  * or using the dir() member function, the repo name, branch,
  * commit hash, and modification status are interrogated with 
  * calls to git.  This includes a list of uncommitted changes,
  * including untracked files.
  *
  * Examples:
  * \code
    mx::sys::gitRepo gr( "/path/to/repo/name");
    std::cout << gr.modified() << '\n'; //will be 1 if modified, 0 if not modified
    std::cout << gr.isNotCommitted("filename") << '\n'; //will be 1 if this file is not committed.  0 otherwise.
    \endcode
  *
  */ 
class gitRepo
{
protected:

   //Set by user:
   std::string m_dir; ///< The directory of the git repository

   //Found using git:
   std::string m_name;      ///< The repo name
   std::string m_branch;    ///< The current branch
   std::string m_hash;      ///< The complete commit hash
   bool m_modified {false}; ///< The modification status, true or false.
   
   std::set<std::string> m_modifiedFiles;  ///< Files which git lists as modified
   std::set<std::string> m_deletedFiles;   ///< Files which git lists as deleted
   std::set<std::string> m_renamedFiles;   ///< Files which git lists as renamed-from
   std::set<std::string> m_renamedFiles2;  ///< Files which git lists as renamed-to
   std::set<std::string> m_untrackedFiles; ///< Files which git lists as untracked

   /// Get the name of the git repo
   /** Called whenever m_dir is set.
     * 
     * \returns 0 on success
     * \returns -1 on error
     */ 
   int getGitName();
   
   /// Get the name of the current commit hash
   /** Called whenever m_dir is set.
     * 
     * \returns 0 on success
     * \returns -1 on error
     */
   int getGitHash();
   
   /// Get the modification status of the repo
   /** Called whenever m_dir is set.
     * 
     * \returns 0 on success
     * \returns -1 on error
     */
   int getGitModified();
   
   /// Get the list of modified files, and the branch name.
   /** Called whenever m_dir is set.
     * 
     * \returns 0 on success
     * \returns -1 on error
     */
   int getGitFileState();
      
public:
   
   /// Default c'tor
   gitRepo();
   
   /// Constructor which sets the directory.
   /** This results in the git repo status being interrogated.
     */
   gitRepo(const std::string & d);
   
   /// Set the directory
   /** This results in the git repo status being interrogated.
     */
   void dir(const std::string & d);
   
   /// Get the current directory
   /** \returns the git repo directory
     */
   std::string dir();
   
   /// Get the current repo's .git directory
   /** \returns the directory plus "/.git"
     */
   std::string gitDir();
   
   /// Get the repo's name
   /** \returns the repo's name
     */
   std::string name();
   
   /// Get the current branch
   /** \returns the current branch name
     */
   std::string branch();
   
   /// Get the current commit hash
   /** \returns the current value of the hash
     */
   std::string hash();
   
   /// Get whether the repo is modified
   /** \returns true is modified
     * \returns false if not modified 
     */
   bool modified();
   
   /// Check whether a file is listed as not committed
   /** Not committed means modified, deleted, renamed (from or to), or untracked.
     *
     * \returns true if not committed
     * \returns false otherwise
     */ 
   bool isNotCommitted(const std::string & file);
   
};

}
}

#endif //gitRepo_hpp
