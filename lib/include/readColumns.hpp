
template<char delim=' ', char eol='\n'> 
void readcol(fstream & fin, int i)
{
   return;
}

template<char delim=' ', char eol='\n', typename arrT, typename... arrTs> 
void readcol(fstream & fin, int i, arrT * array, arrTs*... arrays)
{
   static const unsigned short int nargs = sizeof...(arrTs);
   
   if(nargs > 0) 
   {
      std::string str;
      fin << str;

      array[i] = convertFromString<arrT>(str);

      ++i;
   }
      
   readcol<delim,eol>(fin, i, arrays...);
 
}

template<char delim=' ', char eol='\n', typename... arrTs> 
void readcol(const std::string & fname, int n_lines, arrTs*... arrays)
{
   //open file
   std::ifstream fin;
   fin.open(fname);
   
   for(int j = 0; j< n_lines; j++)
   {
      readcol<delim,eol>(fin, 0, arrays...);
   }
   
   fin.close();
}

