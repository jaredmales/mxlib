#include "gnuplot_interface.h"


int gnuplot_interface_init(gnuplot_interface * gpi)
{
   gpi->gnuplot = 0;
   gpi->gperr = 0;
   gpi->gpout = 0;
   
//    if (mkfifo(MX_GNUPLOT_ERR, 0600)) 
//    {
//       if (errno != EEXIST) 
//       {
//          fprintf(stderr, "could not create gnuplot error fifo\n");
//          return -1;
//       }
//    }
   
   if (mkfifo(MX_GNUPLOT_OUT, 0600)) 
   {
      if (errno != EEXIST) 
      {
         fprintf(stderr, "could not create gnuplot output fifo\n");
         return -1;
      }
   }

   return 0;
}

int gnuplot_interface_connect(gnuplot_interface *gpi)
{
   char cmd[512];
   snprintf(cmd, 512, "gnuplot");
   
   printf("%s\n", cmd);
   
   if ((gpi->gnuplot = popen(cmd,"w")) == NULL) 
   {
      fprintf(stderr, "could not connect to gnuplot\n");
      return -1;
   }
  
   
//    if ((gpi->gpout = fopen(MX_GNUPLOT_OUT,"r")) == NULL) 
//    {
//       fprintf(stderr, "could not open gnuplot output fifo\n");
//       pclose(gpi->gnuplot);
//       return -1;
//    }
   
   //fprintf(gpi->gnuplot, "set print \"%s\"\n\n", MX_GNUPLOT_OUT);
   //fflush(gpi->gnuplot);
   
   return 0;
}

int gnuplot_interface_geterror(gnuplot_interface *gpi, char *str, size_t len)
{
   size_t nb;
   
   if ((gpi->gperr = fopen(MX_GNUPLOT_ERR,"r")) == NULL) 
   {
      return -1;
   }
   
   nb = fread(str, sizeof(char), len, gpi->gperr);
   
   fclose(gpi->gperr);
   
   return nb;
}
int gnuplot_interface_command(gnuplot_interface *gpi, const char * cmd)
{
   size_t nb;
   char buf[1024];
   
   fprintf(gpi->gnuplot, "%s\n", cmd);
   fflush(gpi->gnuplot);

   return 0;
}

int gnuplot_interface_shutdown(gnuplot_interface *gpi)
{
   pclose(gpi->gnuplot);
   //fclose(gpi->gpout);
   
   remove(MX_GNUPLOT_ERR);
   remove(MX_GNUPLOT_OUT);
   
   return 0;
}

