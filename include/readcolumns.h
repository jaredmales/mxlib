/************************************************************
 *    readcolumns.h
 *
 * Author: Jared R. Males (jrmales@email.arizona.edu)
 *
 * Declarations for the readcolumns functions
 *
 * Developed as part of the Magellan Adaptive Optics system.
 ************************************************************/

/** \file readcolumns.h
 * \author Jared R. Males
 * \brief Declarations for the readcolumns functions
 *
 */


#include <iostream>
#include <fstream>

#include <stdarg.h>
#include <cstring>

#include <cstdlib>

#ifndef __readcolumns_h__
#define __readcolumns_h__

/// Read columns of data from a text file.
/** You provide pointers to pointers of the appropriate type.  These are allocated to the correct size
  * based on the number of lines which match the format.
  *
  * The format of the columns is of the form "cslfd", for char, short, long, float, and double.
  * \todo change this to match normal c formatting.
  * We should admit the above form or "%lf%lf%i", that is if you use % you use it the whole time
  * And then specify using a MMatrix<type, 1> with capitals, "DID", "%LF%LF%I".
  *
  * \param fname is the name of the file to read
  * \param delims specifies the column delimiters.  If empty then " " and "," are defaults.
  * \param comment specifies the comment character,  If 0 then "#" is the default.
  * \param format specifies the number of columns and the data type of each column
  * \param ... the remaining arguments are pointers to pointers,which are allocated and filled with the data.
  * \retval -1 on error
  * \retval n the number of columns read. 
  */  
int readcolumns(std::string fname,  std::string delims, char comment, std::string format, ...);

#endif //__readcolumns_h__
