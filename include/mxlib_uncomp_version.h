#ifndef MXLIB_UNCOMP_VERSION_H
#define MXLIB_UNCOMP_VERSION_H

#define MXLIB_UNCOMP_CURRENT_SHA1 "9597a2407a375a08b80a568e5ebf8a3c6b1651b3"
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
