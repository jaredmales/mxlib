#ifndef MXLIB_UNCOMP_VERSION_H
#define MXLIB_UNCOMP_VERSION_H

#define MXLIB_UNCOMP_CURRENT_SHA1 "626b6cf3b3564a0a9b631cedd77dea6fd8be53c9"
#define MXLIB_UNCOMP_REPO_MODIFIED  1


#if MXLIB_UNCOMP_REPO_MODIFIED == 1
  #ifndef GITHEAD_NOWARNING
    #pragma message ("********************************")
    #pragma message ("*                              *")
    #pragma message ("* WARNING: repository modified *")
    #pragma message ("*  changes not committed for   *")
    #pragma message ("*    mxlib    *")
    #pragma message ("*                              *")
    #pragma message ("********************************")
  #endif
#endif


#endif
