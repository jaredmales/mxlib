/** \file ds9_interface.c
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Definitions for the mxlib c ds9 interface
  * 
*/



#include "ds9_interface.h"


int ds9_interface_init(ds9_interface * ds9i)
{
   ds9_interface_set_title(ds9i, "ds9");
   
   ds9i->port = 0;

   ds9i->nsegs = 0;
   ds9i->segs = 0;

   return 0;
}

int ds9_interface_set_title(ds9_interface * ds9i, const char * new_title)
{
   strncpy(ds9i->title, new_title, DS9_TITLE_LENGTH);

   return 0;
}

int ds9_interface_spawn(ds9_interface * ds9i)
{
   int i;
   char resp[32];
   char cmd[DS9_CMD_MAX_LENGTH];
      
   snprintf(cmd, DS9_CMD_MAX_LENGTH, "xpaaccess %s", ds9i->title);
   
   resp[0] = 0;
   if( command_response(cmd, resp, 128) ) return -1;
   
   //Don't respawn if it already exists.
   if(strcmp(resp, "yes\n") == 0) return 0;
   
   system("xpans &");
   
   snprintf(cmd, DS9_CMD_MAX_LENGTH, "ds9 -title %s -port %i &", ds9i->title, ds9i->port);
   
   system(cmd);
   
   //Now wait for ds9 to respond or timeout.
   snprintf(cmd, DS9_CMD_MAX_LENGTH, "xpaaccess %s", ds9i->title);
   
   
   for(i=0;i<10000;i++)
   {
      resp[0] = 0;
      if( command_response(cmd, resp, 128) ) return -1;
      if(strcmp(resp, "yes\n") == 0) return 0;
      usleep(100);
   }
   
   
   return -1; //a timeout
}
   
int ds9_interface_addsegment(ds9_interface *ds9i, int frame)
{
   int i;
   size_t curr_n;
   
   if(frame-1 < ds9i->nsegs) return 0;
   
   curr_n = ds9i->nsegs;
   
   ds9i->nsegs = frame;
   
   ds9i->segs = (sharedmem_segment *) realloc( ds9i->segs, sizeof(sharedmem_segment)*ds9i->nsegs);
      
   for(i = curr_n; i<ds9i->nsegs; i++)
   {
      sharedmem_segment_initialize(&ds9i->segs[i]);
      sharedmem_segment_set_key(&ds9i->segs[i], 0, IPC_PRIVATE);
   }
   
}

int ds9_interface_addframe(ds9_interface *ds9i, int frame)
{
   char cmd[DS9_CMD_MAX_LENGTH];
   
   ds9_interface_addsegment(ds9i, frame);

   snprintf(cmd, DS9_CMD_MAX_LENGTH, "xpaset -p %s frame %i", ds9i->title, frame);
   system(cmd);
   
   
}

int ds9_interface_display_raw(ds9_interface *ds9i, int frame, void *im, size_t dim1, size_t dim2, size_t dim3, int bitpix)
{
   size_t i, tot_size;
   char cmd[DS9_CMD_MAX_LENGTH];
   size_t pixsz;
   
   if(frame < 1)
   {
      fprintf(stderr, "frame must >= 1\n");
      return -1;
   }
   
   ds9_interface_spawn(ds9i);
   
   switch(bitpix)
   {
      case 8:
         pixsz = sizeof(char);
         break;
      case 16:
         pixsz = sizeof(short);
         break;
      case 32:
         pixsz = sizeof(long);
         break;
      case 64:
         pixsz = sizeof(long long);
         break;
      case -32:
         pixsz = sizeof(float);
         break;
      case -64:
         pixsz = sizeof(double);
         break;
      default:
         pixsz = sizeof(double);
   }
   
   //If needed add a shared memory segment for this frame
   if(frame-1 >= ds9i->nsegs)
   {
      ds9_interface_addsegment(ds9i, frame);
   }

   //Calculate total size
   tot_size= pixsz;
   tot_size*=dim1;
   tot_size*=dim2;
   tot_size*=dim3;

   //Re-allocate shared memory if necessary
   if(tot_size > ds9i->segs[frame-1].size)
   {
      if(ds9i->segs[frame-1].size > 0)
      {
         sharedmem_segment_detach(&ds9i->segs[frame-1]);
      }
      sharedmem_segment_create(&ds9i->segs[frame-1], tot_size);
   }
   
   memcpy(ds9i->segs[frame-1].addr, im, tot_size);
   
   
   snprintf(cmd, DS9_CMD_MAX_LENGTH, "xpaset -p %s frame %i", ds9i->title, frame);
   system(cmd);
   
   //Handle single image so that the cube dialog doesn't open up
   if(dim3 == 1)
   {
      snprintf(cmd, DS9_CMD_MAX_LENGTH, "xpaset -p %s shm array shmid %i [xdim=%zu,ydim=%zu,bitpix=%i] &", 
                                      ds9i->title, ds9i->segs[frame-1].shmemid,
                                     dim1, dim2, bitpix);
   }
   else
   {
      snprintf(cmd, DS9_CMD_MAX_LENGTH, "xpaset -p %s shm array shmid %i [xdim=%zu,ydim=%zu,zdim=%zu,bitpix=%i] &", 
                                      ds9i->title, ds9i->segs[frame].shmemid,
                                     dim1, dim2, dim3, bitpix);
   }
   
   system(cmd);
                
   return 0;
   
}

int ds9_interface_shutdown(ds9_interface *ds9i)
{
   size_t i;
   
   for(i=0; i<ds9i->nsegs; i++) sharedmem_segment_detach(&ds9i->segs[i]);
   
   free(ds9i->segs);
   
   ds9i->nsegs = 0;
   
   return 0;
}

ds9_interface * static_ds9(int shutdown)
{
   static int inited = 0;
   static ds9_interface ds9;
   
   if(shutdown)
   {
      if(!inited) return 0;
      
      ds9_interface_shutdown(&ds9);
      
      inited = 0;
      return 0;
   }
   
   if(inited) return &ds9;
   
   ds9_interface_init(&ds9);
   inited = 1;
   
   return &ds9;
}
   
int ds9_display(int frame, void *im, size_t dim1, size_t dim2, size_t dim3, int bitpix)
{
   return ds9_interface_display_raw(static_ds9(0), frame, im, dim1, dim2, dim3, bitpix);
}

int ds9_display_shutdown()
{
   static_ds9(1);
   
   return 0;
}
