

#include "readcolumns.h"
#include "astroconstants.h"

int readcolumns(std::string fname,  std::string delims, char comment, std::string format, ...)
{
   va_list ap;
   
   char ** charptrptr;
   short ** shortptrptr;
   int ** intptrptr;
   long ** longptrptr;
   float ** floatptrptr;
   double ** dblptrptr;
   
   char tmpline[1024];
   int pos, opos;
   std::string argstr, tmpstr;
   
   std::ifstream fin;

   int nargs, nlines;

   nargs = format.length();

   fin.open(fname.c_str(), std::ifstream::in);

   if(!fin.good())
   {
      std::cerr << "Error in file. ";
      if(fin.bad()) std::cerr << "(bad)\n";
      if(fin.fail()) std::cerr << "(fail)\n";
      if(fin.eof()) std::cerr << "(eof)\n";
      return -1;
   }
   
   //First figure out how many lines we have
   nlines = 0;
   while(!fin.eof())
   {
      fin.getline(tmpline, 1024);
      //here test for delimiters
      tmpstr = tmpline;
      pos = tmpstr.find(comment,0);
      opos = tmpstr.find_first_not_of(delims+comment, 0);
      if(pos >=0)
      {
         if(opos == -1 || opos > pos) continue;
      }
      if(opos < 0) continue;

      nlines++;
   }

   
   va_start(ap, format);

   for(int i = 0; i < nargs; i++)
   {
      switch(format[i])
      {
         case 'c':
            charptrptr = va_arg(ap, char**);
            (*charptrptr) = new char[nlines];
            break;
         case 'i':
            intptrptr = va_arg(ap, int**);
            (*intptrptr) = new int[nlines];
            break;            
         case 'd':
            dblptrptr = va_arg(ap, double **);
            (*dblptrptr) = new double[nlines];
            break;
         default:
            std::cerr << "Unrecognized format specifier in readcol\n";
            return -1;
      }
   }

   va_end(ap);

   //fin.seekg(0);
   fin.close();
   fin.open(fname.c_str());
   
   

   delims += "\n";
   
   for(int i=0; i < nlines; i++)
   {
      fin.getline(tmpline, 1024);
      tmpstr = tmpline;
      
      //First check for comments;
      pos = tmpstr.find(comment,0);
      opos = tmpstr.find_first_not_of(delims+comment, 0);
      
      if(pos >=0)
      {
         if(opos == -1 || opos > pos) 
         {
            i--;
            continue;
         }
      }
      if(opos < 0) 
      {
         i--;
         continue;
      }
      //Now find first argument

      //Skip leading delims.
      opos = tmpstr.find_first_not_of(delims+comment, 0);
     
      va_start(ap, format);
     
      for(int j = 0; j < nargs; j++)
      {
         pos = tmpstr.find_first_of(delims+'\n', opos);
         
         if(pos == -1) pos = tmpstr.length();
         
         argstr = tmpstr.substr(opos, pos-opos);
         
         
         opos = tmpstr.find_first_not_of(delims+'\n', pos+1);
         if(opos == -1) opos = tmpstr.length();
         
   
         switch(format[j])
         {
            case 'c':
               charptrptr = va_arg(ap, char**);
               (*charptrptr)[i] = atoi(argstr.c_str());
               break;
            case 'i':
               intptrptr = va_arg(ap, int**);
               (*intptrptr)[i] = atoi(argstr.c_str());
               break;            
            case 'd':
               dblptrptr = va_arg(ap, double **);
               (*dblptrptr)[i] = strtod(argstr.c_str(), 0);
               break;
            default:
               std::cerr << "Unrecognized format specifier in readcol\n";
               return -1;
         }

      }

      va_end(ap);
   }
   fin.close();
   
   return nlines;
}

/*
int main()
{
   int * i1, *i2;
   double *d1, *d2, *d3;
   int nlines;
   
   nlines = readcolumns("testcols.txt",  " ,\t", '#', "ididd", &i1, &d1, &i2, &d2, &d3);
            
   std::cout << "nlines: " << nlines << "\n";
   
   for(int j=0; j<nlines;j++)
   {
      std::cout << i1[j] << " " << d1[j] << " " << i2[j] << " " << d2[j] << " " << d3[j] << "\n";
   }
   
   delete[] i1;
   delete[] i2;
   delete[] d1;
   delete[] d2;
   delete[] d3;
   
   return 0;
}
   
/**/

      
   