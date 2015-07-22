#ifndef MXLIB_UNCOMP_VERSION_H
#define MXLIB_UNCOMP_VERSION_H

#define MXLIB_UNCOMP_CURRENT_SHA1 "cf3a627212f31a72aa54ce1e3e48bbdd5219116b"
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
