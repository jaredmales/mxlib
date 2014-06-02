/** \file msgq.h
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Declarations for the mxlib c message queue facility
  * \ingroup IPC_msgq
  * \ingroup IPC
  * 
*/


#ifndef __mx_msgq_h__
#define __mx_msgq_h__

#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/msg.h>

#include "IPC.h"


#ifdef __cplusplus
extern "C"
{
#endif

   
#define MX_MSGQ_CREAT (0666|IPC_CREAT)

#define MX_MSGQ_TEXTMSG_MAXLEN  1024

/** \addtogroup IPC_msgq
  *  @{
  */

///A simple message c structure.
/** See \ref msgq for how this is used for IPC with a message queue.
  * 
  */
typedef struct 
{
   ///The message type
   long mtype;
   
   ///The message, which has maximum length MX_MSGQ_TEXTMSG_MAXLEN.
   char mtext[MX_MSGQ_TEXTMSG_MAXLEN];
   
} msgq_textmsg;


///Initialize a text message 
/**
  */ 
void msgq_textmsg_initialize(msgq_textmsg *m);

///Set the values of a text message
/**
  */ 
void msgq_textmsg_set(msgq_textmsg *m, long type, const char * str);

///Calculate the length of the text message after the string is populated.
/**
  */ 
size_t msgq_textmsg_length(msgq_textmsg *m);



///A c structure to manage a System V message queue
/** Steps to using this structure:
  * \code
  * msgq q;
  * msgq_initialize(&q);
  * msgq_set_key(&q, "path/to/key/file", 1);
  * msgq_connect(&q);
  * 
  * //Now send a message
  * msgq_textmsg msg;
  * msgq_textmsg_initialize(&msg);
  * msgq_textmsg_set(&msg, 1, "hi")"
  * msgq_send(&msgq, &msg, msgq_textmsg_length(&msg), 0);
  * 
  * //Now receive a message
  * msgq_set_listen_type(&q, 1); //Must always set the message listen_type before calling receive
  * msgq_receive(&q, &msg, sizeof(msgq_textmsg), 0); //this blocks until receipt unless you change msgflag
  * \endcode
  * 
  * \sa Functions for working with msgq are \ref msgq_initialize, \ref msgq_set_key, \ref msgq_set_listen_type, 
  * \ref msgq_connect, \ref msgq_send, and \ref msgq_receive.
  * 
  */
typedef struct
{

   char key_path[MX_IPC_KEYLEN];
   int key_id;

   key_t key;

   int msgq_id;
   
   long listen_type;
   
} msgq;


///Initialize a \ref msgq structure
/** 
  * \param q is a pointer to the \ref msgq to initialize
  * 
  */
void msgq_initialize(msgq * q);

///Set the key for a \ref msgq structure
/** 
  * \param q is a pointer to a \ref msgq
  * \param path is the full path to use in a call to ftok
  * \param id is the id number to use in a call to ftok
  * 
  * \returns the key value, which is also set in the msgq
  * 
  */
key_t msgq_set_key(msgq * q, const char * path, const int id);

///Set the listen message type for a \ref msgq structure
/** 
  * \param q is a pointer to a \ref msgq
  * \param t is the message type
  * 
  */
void msgq_set_listen_type(msgq *q, const long t);

///Connect a \ref msgq structure to a message queue.
/** The key must have already been set by a call to \ref msgq_set_key
  * 
  * \param q is a pointer to the \ref msgq to connect
  * \param msgflag is the message flag argument to pass to msgget
  * \returns the message queue id (>=0) on success
  * \retval -1 on an error
  *
  */
int msgq_connect(msgq * q, int msgflag);

///Send a message using the \ref msgq structure
/** The \ref msgq must be connected.
  * 
  * \param q is a pointer to a \ref msgq
  * \param msg is the message to msg to send 
  * \param msgsz is the size of msg
  * \param msgflag is the msessage flag to pass to msgsnd 
  * 
  * \returns the return value of msgsnd, 0 on success, -1 on failure.
  * 
  */
int msgq_send(msgq * q, const void *msg, size_t msgsz, int msgflag);


///Receive a message using the \ref msgq structure
/** The \ref msgq must be connected, and the listen_type set.
  * 
  * \param q is a pointer to a \ref msgq
  * \param msg is the message structure to populate 
  * \param msgsz is the max size of msg
  * \param msgflag is the message flag to pass to msgrcv
  * 
  * \returns on success, returns the number of bytes read.
  * \returns on failure, (size_t)-1
  * 
  */
size_t msgq_receive(msgq * q, void *msg, size_t msgsz, int msgflag);

/** @}
 */

#ifdef __cplusplus
} //extern "C"
#endif



#endif //__mx_msgq_h__

