

#include "../include/readColumns.hpp"
#include <iostream>

#include <iostream>
#include <vector>
#include <string>


int main()
{
   std::vector<std::string> errnos;

   mx::readColumns("c++11_errnos.txt", errnos);


   for(int i=0; i< errnos.size(); ++i)
   {
      if(errnos[i] == "ENOTSUP")
      {
         std::cout << "      #ifdef " << errnos[i] << "\n";
         std::cout << "      #if (ENOTSUP != EOPNOTSUPP)\n";
         std::cout << "      case " << errnos[i] << ":\n";
         std::cout << "          return \"" << errnos[i] << "\";\n";
         std::cout << "      #endif\n";
         std::cout << "      #endif \n";      
         
         continue;
      }
      
      if(errnos[i] == "EOPNOTSUP")
      {
         std::cout << "      #ifdef " << errnos[i] << "\n";
         std::cout << "      case " << errnos[i] << ":\n";
         std::cout << "          #if (ENOTSUP == EOPNOTSUPP)\n";
         std::cout << "             return \"ENOTSUP / EOPNOTSUPP\";\n";
         std::cout << "          #else\n";
         std::cout << "             return \"" << errnos[i] << "\";\n";
         std::cout << "          #endif\n";
         std::cout << "      #endif \n";      
         
         continue;
      }
      
      if(errnos[i] == "EWOULDBLOCK")
      {
         std::cout << "      #ifdef " << errnos[i] << "\n";
         std::cout << "      #if (EWOULDBLOCK != EAGAIN)\n";
         std::cout << "      case " << errnos[i] << ":\n";
         std::cout << "          return \"" << errnos[i] << "\";\n";
         std::cout << "      #endif\n";
         std::cout << "      #endif \n";      
         
         continue;
      }
      
      if(errnos[i] == "EAGAIN")
      {
         std::cout << "      #ifdef " << errnos[i] << "\n";
         std::cout << "      case " << errnos[i] << ":\n";
         std::cout << "          #if (EWOULDBLOCK == EAGAIN)\n";
         std::cout << "             return \"EAGIAN / EWOULDBLOCK\";\n";
         std::cout << "          #else\n";
         std::cout << "             return \"" << errnos[i] << "\";\n";
         std::cout << "          #endif\n";
         std::cout << "      #endif \n";      
         
         continue;
      }
      
      std::cout << "      #ifdef " << errnos[i] << "\n";
      std::cout << "      case " << errnos[i] << ":\n";
      std::cout << "          return \"" << errnos[i] << "\";\n";
      std::cout << "      #endif \n";
   }

}
