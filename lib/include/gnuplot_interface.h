
#ifndef __gnuplot_interface_h__
#define __gnuplot_interface_h__

#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

#ifdef __cplusplus
extern "C"
{
#endif
   
#define MX_GNUPLOT_OUT    "./gpout"
#define MX_GNUPLOT_ERR    "./gperr"

typedef struct
{
   FILE * gnuplot;
   FILE * gperr;
   FILE * gpout;
   
   int plotno;
   
} gnuplot_interface;

int gnuplot_interface_init(gnuplot_interface * gpi);

int gnuplot_interface_connect(gnuplot_interface *gpi);

int gnuplot_interface_geterror(gnuplot_interface *gpi, char *str, size_t len);

int gnuplot_interface_command(gnuplot_interface *gpi, const char * cmd);

int gnuplot_interface_plot(gnuplot_interface *gpi, double *x, double *y, int n);
int gnuplot_interface_replot(gnuplot_interface *gpi, double *x, double *y, int n);

int gnuplot_interface_shutdown(gnuplot_interface *gpi);


int gnuplot_command(const char *cmd);
int gnuplot_plot(double *x, double *y, int n);
int gnuplot_replot(double *x, double *y, int n);

int gnuplot_display_shutdown();


#ifdef __cplusplus
} //extern "C"
#endif
   
#endif // __gnuplot_interface_h__
