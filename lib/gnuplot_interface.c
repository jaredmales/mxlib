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

   gpi->plotno = 0;
   
   return 0;
}

int gnuplot_interface_connect(gnuplot_interface *gpi)
{
   char cmd[512];
   snprintf(cmd, 512, "gnuplot");
   
   
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

int gnuplot_interface_plot(gnuplot_interface *gpi, double *x, double *y, int n)
{
   FILE * f;
   char fname[256];
   char cmd[256]; 
   
   snprintf(fname, 256, "/tmp/gnuplot_interface_out%i.dat", gpi->plotno);
   gpi->plotno ++;
   
   snprintf(cmd, 256, "plot \"%s\" using 1:2 notitle\n", fname);
            
   f = fopen(fname, "w");
  
   for(int i=0;i<n;i++)
   {
      fprintf(f, "%f %f\n", x[i], y[i]);
   }
   
   fclose(f);
   
   
   return gnuplot_interface_command(gpi, cmd);
}
 
int gnuplot_interface_replot(gnuplot_interface *gpi, double *x, double *y, int n)
{
   FILE * f;
    char fname[256];
   char cmd[256]; 
   
   snprintf(fname, 256, "/tmp/gnuplot_interface_out%i.dat", gpi->plotno);
   gpi->plotno++;
   
   snprintf(cmd, 256, "replot \"%s\" using 1:2 notitle\n", fname);
   
   f = fopen(fname, "w");
  
   for(int i=0;i<n;i++)
   {
      fprintf(f, "%f %f\n", x[i], y[i]);
   }
   
   fclose(f);
   
   
   return gnuplot_interface_command(gpi, cmd);
} 
   


int gnuplot_interface_shutdown(gnuplot_interface *gpi)
{
   pclose(gpi->gnuplot);
   //fclose(gpi->gpout);
   
   remove(MX_GNUPLOT_ERR);
   remove(MX_GNUPLOT_OUT);
   
   return 0;
}


gnuplot_interface * static_gnuplot(int shutdown)
{
   static int inited = 0;
   static gnuplot_interface gnuplot;
   
   if(shutdown)
   {
      if(!inited) return 0;
      
      gnuplot_interface_shutdown(&gnuplot);
      
      inited = 0;
      return 0;
   }
   
   if(inited) return &gnuplot;
   
   gnuplot_interface_init(&gnuplot);
   gnuplot_interface_connect(&gnuplot);
   inited = 1;
   
   return &gnuplot;
}
 
int gnuplot_command(const char *cmd)
{
   return gnuplot_interface_command(static_gnuplot(0), cmd);
}

int gnuplot_plot(double *x, double *y, int n)
{
   return gnuplot_interface_plot(static_gnuplot(0), x, y, n);
}

int gnuplot_replot(double *x, double *y, int n)
{
   return gnuplot_interface_replot(static_gnuplot(0), x, y, n);
}

int gnuplot_display_shutdown()
{
   static_gnuplot(1);
   
   return 0;
}
