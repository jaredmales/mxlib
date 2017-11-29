/** \file msgq.c
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Definitions for the mxlib c message queue facility
  * 
*/

#include "msgq.h"

void msgq_textmsg_initialize(msgq_textmsg *m)
{
   m->mtype = 0;
   m->mtext[0] = '\0';
}

void msgq_textmsg_set(msgq_textmsg *m, long type, const char * str)
{
   m->mtype = type;
   strncpy(m->mtext, str, MX_MSGQ_TEXTMSG_MAXLEN);
}

size_t msgq_textmsg_length(msgq_textmsg *m)
{
   return sizeof(long) + (strlen(m->mtext)+1)*sizeof(char);
}


void msgq_initialize(msgq * q)
{
   q->key_id = 0;
   q->key = -1;
   q->msgq_id = -1;
}


key_t msgq_set_key(msgq * q, const char * path, const int id)
{
   strncpy(q->key_path, path, MX_IPC_KEYLEN);
   q->key_id = id;
   
   q->key = ftok(q->key_path, q->key_id);

   return q->key_id;
}

void msgq_set_listen_type(msgq *q, const long t)
{
   q->listen_type = t;
}

int msgq_connect(msgq * q, int msgflag)
{
   if(q->key < 0)
   {
      errno = EINVAL;
      return -1;
   }

   q->msgq_id = msgget(q->key, msgflag);
   
   return q->msgq_id;
}

int msgq_send(msgq * q, const void *msg, size_t msgsz, int msgflag)
{
   return msgsnd(q->msgq_id, msg, msgsz, msgflag);
}

size_t msgq_receive(msgq * q, void *msg, size_t msgsz, int msgflag)
{
   return msgrcv(q->msgq_id, msg, msgsz, q->listen_type, msgflag);
}

