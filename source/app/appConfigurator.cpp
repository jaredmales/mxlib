/** \file appConfigurator.cpp
 * \author Jared R. Males
 * \brief Implementation of An application configuration manager
 *
 * \ingroup mxApp_files
 *
 */

//***********************************************************************//
// Copyright 2021 Jared R. Males (jaredmales@gmail.com)
//
// This file is part of mxlib.
//
// mxlib is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// mxlib is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with mxlib.  If not, see <http://www.gnu.org/licenses/>.
//***********************************************************************//

#include "app/appConfigurator.hpp"

namespace mx
{
namespace app 
{

void appConfigurator::clear()
{
   m_targets.clear();
   clOnlyTargets.clear();
   nonOptions.clear();
   m_unusedConfigs.clear();
}

void appConfigurator::add( const configTarget & tgt )
{

   //First check for duplicate name and command line only
   if(m_targets.count(tgt.name) > 0 && tgt.section == "" && tgt.keyword == "")
   {
      clOnlyTargets.push_back(tgt);
      clOnlyTargets.back().orderAdded = nAdded;
   }
   else
   {
      std::pair<targetIterator, bool> res = m_targets.insert({tgt.name, tgt});
      res.first->second.orderAdded = nAdded;
   }

   ++nAdded;

}

void appConfigurator::add( const std::string &n,
                           const std::string &so,
                           const std::string &lo,
                           int clt,
                           const std::string & s,
                           const std::string & kw,
                           bool isReq,
                           const std::string & ht,
                           const std::string & he )
{
   add( configTarget(n,so,lo,clt, s, kw, isReq, ht, he) );
}

void appConfigurator::parseCommandLine( int argc,
                                        char ** argv,
                                        const std::string & oneTarget
                                      )
{
   if(argc == 0) return;

   clOptions clOpts;

   //First load the options into the clOptions parser
   targetIterator it;
   for(it = m_targets.begin(); it != m_targets.end(); ++it)
   {
      if(it->second.shortOpt == "" && it->second.longOpt == "") continue; //No command line opts specified.

      clOpts.add(it->second.name,it->second.shortOpt.c_str(),it->second.longOpt.c_str(), it->second.clType);
   }

   //Then load the command-line-only options.
   clOnlyTargetIterator cloit;
   for(cloit = clOnlyTargets.begin(); cloit != clOnlyTargets.end(); ++cloit)
   {
      if(cloit->shortOpt == "" && cloit->longOpt == "") continue; //Nothing to add?

      clOpts.add(cloit->name,cloit->shortOpt.c_str(),cloit->longOpt.c_str(), cloit->clType);
   }

   //Now we parse
   clOpts.parse(argc, argv, &nonOptions);

   if(clOpts.numUnknown() > 0)
   {
      std::vector<std::string> unk;
      clOpts.unknown(unk);
      
      for (size_t n=0;n<unk.size();++n)
      {
         //Insert or update existing
         m_unusedConfigs[unk[n]].name = unk[n];
         if(m_sources) m_unusedConfigs[unk[n]].sources.push_back("command line");
         m_unusedConfigs[unk[n]].set = true;
      }
   }
      
   //If nothing more to do, get out
   if(clOpts.nOpts == 0)
   {
      return;
   }
   //And then load the results in the config target map.
   for(it = m_targets.begin(); it != m_targets.end(); ++it)
   {
      if(oneTarget != "" && it->second.name != oneTarget) continue;

      if(clOpts.optSet(it->second.name))
      {
         std::vector<std::string> args;

         clOpts.getAll(args, it->second.name);
         
         it->second.values.insert( it->second.values.end(), args.begin(), args.end());
         if(m_sources)
         {
            for(size_t n=0; n < args.size(); ++n) it->second.sources.push_back("command line");
         }
         
         it->second.verbosity = clOpts.count(it->second.name);
         it->second.set = true;
      }
   }


}

int appConfigurator::readConfig( const std::string & fname,
                                 bool reportFileNotFound
                               )
{
   //Handle empty string quietly
   if(fname == "") return 0;

   iniFile iF;

   ///\todo update error handling to include >0 (line numer of parse error) and -2 memory allocation error.
   int prv = iF.parse(fname); 
   
   if( prv == -1)
   {
      if(!reportFileNotFound) return -1;

      mxError("appConfigurator: ", MXE_FILENOTFOUND, "The file " + fname + " was not found");
      return -1;
   }
   
   if( prv == -2)
   {
      mxError("appConfigurator: ", MXE_ALLOCERR, "Memory allocation error in config file parser");
      return -1;
   }
   
   if( prv > 0)
   {
      mxError("appConfigurator: ", MXE_PARSEERR, "Parsing error in " + fname + " at line " + std::to_string(prv));
      return -1;
   }

   targetIterator it;

   for(it = m_targets.begin(); it != m_targets.end(); ++it)
   {
      if(iF.count(it->second.section, it->second.keyword) > 0)
      {
         it->second.values.push_back(iF(it->second.section, it->second.keyword));
         if(m_sources)
         {
            it->second.sources.push_back(fname);
         }
         it->second.set = true;
         
         iF.erase(it->second.section, it->second.keyword); //Erase it from the iniFile map.
      }
   }
   
   //here set aside non-deleted iF entries
   for( auto iFit = iF.names.begin(); iFit != iF.names.end(); ++iFit)
   {
      //Insert or update existing
      m_unusedConfigs[iFit->first].name = iFit->first;
      
      std::string sect, nam;
      iniFile::parseKey(sect, nam, iFit->first);
      m_unusedConfigs[iFit->first].section = sect;
      m_unusedConfigs[iFit->first].keyword = nam;
      
      m_unusedConfigs[iFit->first].values.push_back(iFit->second);
      
      if(m_sources) m_unusedConfigs[iFit->first].sources.push_back(fname);
      
      m_unusedConfigs[iFit->first].set = true;
   }
   
   return 0;
}

bool appConfigurator::isSet( const std::string & name,
                             std::unordered_map<std::string, configTarget> & targets
                           )
{
   if(targets.count(name) == 0) return false;
   return targets[name].set;
}

bool appConfigurator::isSet( const std::string & name )
{
   return isSet(name, m_targets);
}

int appConfigurator::count( const std::string & name,
                            std::unordered_map<std::string, configTarget> & targets
                          )
{
   return targets[name].values.size();
}

int appConfigurator::count( const std::string & name )
{
   return count(name, m_targets);
}

int appConfigurator::verbosity( const std::string & name,
                                std::unordered_map<std::string, configTarget> & targets
                              )
{
   return targets[name].verbosity;
}

int appConfigurator::verbosity( const std::string & name )
{
   return verbosity(name, m_targets);
}

//+++++++++++++++++++++++++++++++++++++
// Explicit instants:
//+++++++++++++++++++++++++++++++++++++

template
int appConfigurator::get<char>( char & v, const std::string & name, size_t i, std::unordered_map<std::string, configTarget> & targets);

template
int appConfigurator::get<char16_t>( char16_t & v, const std::string & name, size_t i, std::unordered_map<std::string, configTarget> & targets);

template
int appConfigurator::get<char32_t>( char32_t & v, const std::string & name, size_t i, std::unordered_map<std::string, configTarget> & targets);

template
int appConfigurator::get<wchar_t>( wchar_t & v, const std::string & name, size_t i, std::unordered_map<std::string, configTarget> & targets);

template
int appConfigurator::get<signed char>( signed char & v, const std::string & name, size_t i, std::unordered_map<std::string, configTarget> & targets);

template
int appConfigurator::get<short>( short & v, const std::string & name, size_t i, std::unordered_map<std::string, configTarget> & targets);

template
int appConfigurator::get<int>( int & v, const std::string & name, size_t i, std::unordered_map<std::string, configTarget> & targets);

template
int appConfigurator::get<long>( long & v, const std::string & name, size_t i, std::unordered_map<std::string, configTarget> & targets);

template
int appConfigurator::get<long long>( long long & v, const std::string & name, size_t i, std::unordered_map<std::string, configTarget> & targets);

template
int appConfigurator::get<unsigned char>( unsigned char & v, const std::string & name, size_t i, std::unordered_map<std::string, configTarget> & targets);

template
int appConfigurator::get<unsigned short>( unsigned short & v, const std::string & name, size_t i, std::unordered_map<std::string, configTarget> & targets);

template
int appConfigurator::get<unsigned int>( unsigned int & v, const std::string & name, size_t i, std::unordered_map<std::string, configTarget> & targets);

template
int appConfigurator::get<unsigned long>( unsigned long & v, const std::string & name, size_t i, std::unordered_map<std::string, configTarget> & targets);

template
int appConfigurator::get<unsigned long long>( unsigned long long & v, const std::string & name, size_t i, std::unordered_map<std::string, configTarget> & targets);

template
int appConfigurator::get<float>( float & v, const std::string & name, size_t i, std::unordered_map<std::string, configTarget> & targets);

template
int appConfigurator::get<double>( double & v, const std::string & name, size_t i, std::unordered_map<std::string, configTarget> & targets);

template
int appConfigurator::get<long double>( long double & v, const std::string & name, size_t i, std::unordered_map<std::string, configTarget> & targets);

#ifdef HASQUAD
template
int appConfigurator::get<__float128>( __float128 & v, const std::string & name, size_t i, std::unordered_map<std::string, configTarget> & targets);
#endif

template
int appConfigurator::get<bool>( bool & v, const std::string & name, size_t i, std::unordered_map<std::string, configTarget> & targets);

template
int appConfigurator::get<std::string>( std::string & v, const std::string & name, size_t i, std::unordered_map<std::string, configTarget> & targets);

//+++++++++++++++++++++++++++++++++++++

template
int appConfigurator::get<char>( char & v, const std::string & name, size_t i);

template
int appConfigurator::get<char16_t>( char16_t & v, const std::string & name, size_t i);

template
int appConfigurator::get<char32_t>( char32_t & v, const std::string & name, size_t i);

template
int appConfigurator::get<wchar_t>( wchar_t & v, const std::string & name, size_t i);

template
int appConfigurator::get<signed char>( signed char & v, const std::string & name, size_t i);

template
int appConfigurator::get<short>( short & v, const std::string & name, size_t i);

template
int appConfigurator::get<int>( int & v, const std::string & name, size_t i);

template
int appConfigurator::get<long>( long & v, const std::string & name, size_t i);

template
int appConfigurator::get<long long>( long long & v, const std::string & name, size_t i);

template
int appConfigurator::get<unsigned char>( unsigned char & v, const std::string & name, size_t i);

template
int appConfigurator::get<unsigned short>( unsigned short & v, const std::string & name, size_t i);

template
int appConfigurator::get<unsigned int>( unsigned int & v, const std::string & name, size_t i);

template
int appConfigurator::get<unsigned long>( unsigned long & v, const std::string & name, size_t i);

template
int appConfigurator::get<unsigned long long>( unsigned long long & v, const std::string & name, size_t i);

template
int appConfigurator::get<float>( float & v, const std::string & name, size_t i);

template
int appConfigurator::get<double>( double & v, const std::string & name, size_t i);

template
int appConfigurator::get<long double>( long double & v, const std::string & name, size_t i);

#ifdef HASQUAD
template
int appConfigurator::get<__float128>( __float128 & v, const std::string & name, size_t i);
#endif

template
int appConfigurator::get<bool>( bool & v, const std::string & name, size_t i);

template
int appConfigurator::get<std::string>( std::string & v, const std::string & name, size_t i);

//+++++++++++++++++++++++++++++++++++++

template
int appConfigurator::get<char>( char & v, const std::string & name, std::unordered_map<std::string, configTarget> & targets);

template
int appConfigurator::get<char16_t>( char16_t & v, const std::string & name, std::unordered_map<std::string, configTarget> & targets);

template
int appConfigurator::get<char32_t>( char32_t & v, const std::string & name, std::unordered_map<std::string, configTarget> & targets);

template
int appConfigurator::get<wchar_t>( wchar_t & v, const std::string & name, std::unordered_map<std::string, configTarget> & targets);

template
int appConfigurator::get<signed char>( signed char & v, const std::string & name, std::unordered_map<std::string, configTarget> & targets);

template
int appConfigurator::get<short>( short & v, const std::string & name, std::unordered_map<std::string, configTarget> & targets);

template
int appConfigurator::get<int>( int & v, const std::string & name, std::unordered_map<std::string, configTarget> & targets);

template
int appConfigurator::get<long>( long & v, const std::string & name, std::unordered_map<std::string, configTarget> & targets);

template
int appConfigurator::get<long long>( long long & v, const std::string & name, std::unordered_map<std::string, configTarget> & targets);

template
int appConfigurator::get<unsigned char>( unsigned char & v, const std::string & name, std::unordered_map<std::string, configTarget> & targets);

template
int appConfigurator::get<unsigned short>( unsigned short & v, const std::string & name, std::unordered_map<std::string, configTarget> & targets);

template
int appConfigurator::get<unsigned int>( unsigned int & v, const std::string & name, std::unordered_map<std::string, configTarget> & targets);

template
int appConfigurator::get<unsigned long>( unsigned long & v, const std::string & name, std::unordered_map<std::string, configTarget> & targets);

template
int appConfigurator::get<unsigned long long>( unsigned long long & v, const std::string & name, std::unordered_map<std::string, configTarget> & targets);

template
int appConfigurator::get<float>( float & v, const std::string & name, std::unordered_map<std::string, configTarget> & targets);

template
int appConfigurator::get<double>( double & v, const std::string & name, std::unordered_map<std::string, configTarget> & targets);

template
int appConfigurator::get<long double>( long double & v, const std::string & name, std::unordered_map<std::string, configTarget> & targets);

#ifdef HASQUAD
template
int appConfigurator::get<__float128>( __float128 & v, const std::string & name, std::unordered_map<std::string, configTarget> & targets);
#endif

template
int appConfigurator::get<bool>( bool & v, const std::string & name, std::unordered_map<std::string, configTarget> & targets);

template
int appConfigurator::get<std::string>( std::string & v, const std::string & name, std::unordered_map<std::string, configTarget> & targets);

//+++++++++++++++++++++++++++++++++++++

template
int appConfigurator::get<char>( char & v, const std::string & name);

template
int appConfigurator::get<char16_t>( char16_t & v, const std::string & name);

template
int appConfigurator::get<char32_t>( char32_t & v, const std::string & name);

template
int appConfigurator::get<wchar_t>( wchar_t & v, const std::string & name);

template
int appConfigurator::get<signed char>( signed char & v, const std::string & name);

template
int appConfigurator::get<short>( short & v, const std::string & name);

template
int appConfigurator::get<int>( int & v, const std::string & name);

template
int appConfigurator::get<long>( long & v, const std::string & name);

template
int appConfigurator::get<long long>( long long & v, const std::string & name);

template
int appConfigurator::get<unsigned char>( unsigned char & v, const std::string & name);

template
int appConfigurator::get<unsigned short>( unsigned short & v, const std::string & name);

template
int appConfigurator::get<unsigned int>( unsigned int & v, const std::string & name);

template
int appConfigurator::get<unsigned long>( unsigned long & v, const std::string & name);

template
int appConfigurator::get<unsigned long long>( unsigned long long & v, const std::string & name);

template
int appConfigurator::get<float>( float & v, const std::string & name);

template
int appConfigurator::get<double>( double & v, const std::string & name);

template
int appConfigurator::get<long double>( long double & v, const std::string & name);

#ifdef HASQUAD
template
int appConfigurator::get<__float128>( __float128 & v, const std::string & name);
#endif

template
int appConfigurator::get<bool>( bool & v, const std::string & name);

template
int appConfigurator::get<std::string>( std::string & v, const std::string & name);

//++++++++++++++++++++++++++++++++++++++++++++++

template
int appConfigurator::get<char>( std::vector<char> & v, const std::string & name, size_t i, std::unordered_map<std::string, configTarget> & targets);

template
int appConfigurator::get<char16_t>( std::vector<char16_t> & v, const std::string & name, size_t i, std::unordered_map<std::string, configTarget> & targets);

template
int appConfigurator::get<char32_t>( std::vector<char32_t> & v, const std::string & name, size_t i, std::unordered_map<std::string, configTarget> & targets);

template
int appConfigurator::get<wchar_t>( std::vector<wchar_t> & v, const std::string & name, size_t i, std::unordered_map<std::string, configTarget> & targets);

template
int appConfigurator::get<signed char>( std::vector<signed char> & v, const std::string & name, size_t i, std::unordered_map<std::string, configTarget> & targets);

template
int appConfigurator::get<short>( std::vector<short> & v, const std::string & name, size_t i, std::unordered_map<std::string, configTarget> & targets);

template
int appConfigurator::get<int>( std::vector<int> & v, const std::string & name, size_t i, std::unordered_map<std::string, configTarget> & targets);

template
int appConfigurator::get<long>( std::vector<long> & v, const std::string & name, size_t i, std::unordered_map<std::string, configTarget> & targets);

template
int appConfigurator::get<long long>( std::vector<long long> & v, const std::string & name, size_t i, std::unordered_map<std::string, configTarget> & targets);

template
int appConfigurator::get<unsigned char>( std::vector<unsigned char> & v, const std::string & name, size_t i, std::unordered_map<std::string, configTarget> & targets);

template
int appConfigurator::get<unsigned short>( std::vector<unsigned short> & v, const std::string & name, size_t i, std::unordered_map<std::string, configTarget> & targets);

template
int appConfigurator::get<unsigned int>( std::vector<unsigned int> & v, const std::string & name, size_t i, std::unordered_map<std::string, configTarget> & targets);

template
int appConfigurator::get<unsigned long>( std::vector<unsigned long> & v, const std::string & name, size_t i, std::unordered_map<std::string, configTarget> & targets);

template
int appConfigurator::get<unsigned long long>( std::vector<unsigned long long> & v, const std::string & name, size_t i, std::unordered_map<std::string, configTarget> & targets);

template
int appConfigurator::get<float>( std::vector<float> & v, const std::string & name, size_t i, std::unordered_map<std::string, configTarget> & targets);

template
int appConfigurator::get<double>( std::vector<double> & v, const std::string & name, size_t i, std::unordered_map<std::string, configTarget> & targets);

template
int appConfigurator::get<long double>( std::vector<long double> & v, const std::string & name, size_t i, std::unordered_map<std::string, configTarget> & targets);

#ifdef HASQUAD
template
int appConfigurator::get<__float128>( std::vector<__float128> & v, const std::string & name, size_t i, std::unordered_map<std::string, configTarget> & targets);
#endif

template
int appConfigurator::get<bool>( std::vector<bool> & v, const std::string & name, size_t i, std::unordered_map<std::string, configTarget> & targets);

template
int appConfigurator::get<std::string>( std::vector<std::string> & v, const std::string & name, size_t i, std::unordered_map<std::string, configTarget> & targets);

//+++++++++++++++++++++++++++++++++++++

template
int appConfigurator::get<char>( std::vector<char> & v, const std::string & name, size_t i);

template
int appConfigurator::get<char16_t>( std::vector<char16_t> & v, const std::string & name, size_t i);

template
int appConfigurator::get<char32_t>( std::vector<char32_t> & v, const std::string & name, size_t i);

template
int appConfigurator::get<wchar_t>( std::vector<wchar_t> & v, const std::string & name, size_t i);

template
int appConfigurator::get<signed char>( std::vector<signed char> & v, const std::string & name, size_t i);

template
int appConfigurator::get<short>( std::vector<short> & v, const std::string & name, size_t i);

template
int appConfigurator::get<int>( std::vector<int> & v, const std::string & name, size_t i);

template
int appConfigurator::get<long>( std::vector<long> & v, const std::string & name, size_t i);

template
int appConfigurator::get<long long>( std::vector<long long> & v, const std::string & name, size_t i);

template
int appConfigurator::get<unsigned char>( std::vector<unsigned char> & v, const std::string & name, size_t i);

template
int appConfigurator::get<unsigned short>( std::vector<unsigned short> & v, const std::string & name, size_t i);

template
int appConfigurator::get<unsigned int>( std::vector<unsigned int> & v, const std::string & name, size_t i);

template
int appConfigurator::get<unsigned long>( std::vector<unsigned long> & v, const std::string & name, size_t i);

template
int appConfigurator::get<unsigned long long>( std::vector<unsigned long long> & v, const std::string & name, size_t i);

template
int appConfigurator::get<float>( std::vector<float> & v, const std::string & name, size_t i);

template
int appConfigurator::get<double>( std::vector<double> & v, const std::string & name, size_t i);

template
int appConfigurator::get<long double>( std::vector<long double> & v, const std::string & name, size_t i);

#ifdef HASQUAD
template
int appConfigurator::get<__float128>( std::vector<__float128> & v, const std::string & name, size_t i);
#endif

template
int appConfigurator::get<bool>( std::vector<bool> & v, const std::string & name, size_t i);

template
int appConfigurator::get<std::string>( std::vector<std::string> & v, const std::string & name, size_t i);

//+++++++++++++++++++++++++++++++++++++++++++++++
template
int appConfigurator::get<char>( std::vector<char> & v, const std::string & name, std::unordered_map<std::string, configTarget> & targets);

template
int appConfigurator::get<char16_t>( std::vector<char16_t> & v, const std::string & name, std::unordered_map<std::string, configTarget> & targets);

template
int appConfigurator::get<char32_t>( std::vector<char32_t> & v, const std::string & name, std::unordered_map<std::string, configTarget> & targets);

template
int appConfigurator::get<wchar_t>( std::vector<wchar_t> & v, const std::string & name, std::unordered_map<std::string, configTarget> & targets);

template
int appConfigurator::get<signed char>( std::vector<signed char> & v, const std::string & name, std::unordered_map<std::string, configTarget> & targets);

template
int appConfigurator::get<short>( std::vector<short> & v, const std::string & name, std::unordered_map<std::string, configTarget> & targets);

template
int appConfigurator::get<int>( std::vector<int> & v, const std::string & name, std::unordered_map<std::string, configTarget> & targets);

template
int appConfigurator::get<long>( std::vector<long> & v, const std::string & name, std::unordered_map<std::string, configTarget> & targets);

template
int appConfigurator::get<long long>( std::vector<long long> & v, const std::string & name, std::unordered_map<std::string, configTarget> & targets);

template
int appConfigurator::get<unsigned char>( std::vector<unsigned char> & v, const std::string & name, std::unordered_map<std::string, configTarget> & targets);

template
int appConfigurator::get<unsigned short>( std::vector<unsigned short> & v, const std::string & name, std::unordered_map<std::string, configTarget> & targets);

template
int appConfigurator::get<unsigned int>( std::vector<unsigned int> & v, const std::string & name, std::unordered_map<std::string, configTarget> & targets);

template
int appConfigurator::get<unsigned long>( std::vector<unsigned long> & v, const std::string & name, std::unordered_map<std::string, configTarget> & targets);

template
int appConfigurator::get<unsigned long long>( std::vector<unsigned long long> & v, const std::string & name, std::unordered_map<std::string, configTarget> & targets);

template
int appConfigurator::get<float>( std::vector<float> & v, const std::string & name, std::unordered_map<std::string, configTarget> & targets);

template
int appConfigurator::get<double>( std::vector<double> & v, const std::string & name, std::unordered_map<std::string, configTarget> & targets);

template
int appConfigurator::get<long double>( std::vector<long double> & v, const std::string & name, std::unordered_map<std::string, configTarget> & targets);

#ifdef HASQUAD
template
int appConfigurator::get<__float128>( std::vector<__float128> & v, const std::string & name, std::unordered_map<std::string, configTarget> & targets);
#endif

template
int appConfigurator::get<bool>( std::vector<bool> & v, const std::string & name, std::unordered_map<std::string, configTarget> & targets);

template
int appConfigurator::get<std::string>( std::vector<std::string> & v, const std::string & name, std::unordered_map<std::string, configTarget> & targets);
//+++++++++++++++++++++++++++++++++++++

template
int appConfigurator::get<char>( std::vector<char> & v, const std::string & name);

template
int appConfigurator::get<char16_t>( std::vector<char16_t> & v, const std::string & name);

template
int appConfigurator::get<char32_t>( std::vector<char32_t> & v, const std::string & name);

template
int appConfigurator::get<wchar_t>( std::vector<wchar_t> & v, const std::string & name);

template
int appConfigurator::get<signed char>( std::vector<signed char> & v, const std::string & name);

template
int appConfigurator::get<short>( std::vector<short> & v, const std::string & name);

template
int appConfigurator::get<int>( std::vector<int> & v, const std::string & name);

template
int appConfigurator::get<long>( std::vector<long> & v, const std::string & name);

template
int appConfigurator::get<long long>( std::vector<long long> & v, const std::string & name);

template
int appConfigurator::get<unsigned char>( std::vector<unsigned char> & v, const std::string & name);

template
int appConfigurator::get<unsigned short>( std::vector<unsigned short> & v, const std::string & name);

template
int appConfigurator::get<unsigned int>( std::vector<unsigned int> & v, const std::string & name);

template
int appConfigurator::get<unsigned long>( std::vector<unsigned long> & v, const std::string & name);

template
int appConfigurator::get<unsigned long long>( std::vector<unsigned long long> & v, const std::string & name);

template
int appConfigurator::get<float>( std::vector<float> & v, const std::string & name);

template
int appConfigurator::get<double>( std::vector<double> & v, const std::string & name);

template
int appConfigurator::get<long double>( std::vector<long double> & v, const std::string & name);

#ifdef HASQUAD
template
int appConfigurator::get<__float128>( std::vector<__float128> & v, const std::string & name);
#endif

template
int appConfigurator::get<bool>( std::vector<bool> & v, const std::string & name);

template
int appConfigurator::get<std::string>( std::vector<std::string> & v, const std::string & name);

//+++++++++++++++++++++++++++++++++++++++

template
int appConfigurator::operator()<char>( char & v, const std::string & name);

template
int appConfigurator::operator()<char16_t>( char16_t & v, const std::string & name);

template
int appConfigurator::operator()<char32_t>( char32_t & v, const std::string & name);

template
int appConfigurator::operator()<wchar_t>( wchar_t & v, const std::string & name);

template
int appConfigurator::operator()<signed char>( signed char & v, const std::string & name);

template
int appConfigurator::operator()<short>( short & v, const std::string & name);

template
int appConfigurator::operator()<int>( int & v, const std::string & name);

template
int appConfigurator::operator()<long>( long & v, const std::string & name);

template
int appConfigurator::operator()<long long>( long long & v, const std::string & name);

template
int appConfigurator::operator()<unsigned char>( unsigned char & v, const std::string & name);

template
int appConfigurator::operator()<unsigned short>( unsigned short & v, const std::string & name);

template
int appConfigurator::operator()<unsigned int>( unsigned int & v, const std::string & name);

template
int appConfigurator::operator()<unsigned long>( unsigned long & v, const std::string & name);

template
int appConfigurator::operator()<unsigned long long>( unsigned long long & v, const std::string & name);

template
int appConfigurator::operator()<float>( float & v, const std::string & name);

template
int appConfigurator::operator()<double>( double & v, const std::string & name);

template
int appConfigurator::operator()<long double>( long double & v, const std::string & name);

#ifdef HASQUAD
template
int appConfigurator::operator()<__float128>( __float128 & v, const std::string & name);
#endif

template
int appConfigurator::operator()<bool>( bool & v, const std::string & name);

template
int appConfigurator::operator()<std::string>( std::string & v, const std::string & name);

//++++++++++++++++++++++++++++++++++++++++++++++

template
int appConfigurator::configUnused<char>( char & v, const std::string & key);

template
int appConfigurator::configUnused<char16_t>( char16_t & v, const std::string & key);

template
int appConfigurator::configUnused<char32_t>( char32_t & v, const std::string & key);

template
int appConfigurator::configUnused<wchar_t>( wchar_t & v, const std::string & key);

template
int appConfigurator::configUnused<signed char>( signed char & v, const std::string & key);

template
int appConfigurator::configUnused<short>( short & v, const std::string & key);

template
int appConfigurator::configUnused<int>( int & v, const std::string & key);

template
int appConfigurator::configUnused<long>( long & v, const std::string & key);

template
int appConfigurator::configUnused<long long>( long long & v, const std::string & key);

template
int appConfigurator::configUnused<unsigned char>( unsigned char & v, const std::string & key);

template
int appConfigurator::configUnused<unsigned short>( unsigned short & v, const std::string & key);

template
int appConfigurator::configUnused<unsigned int>( unsigned int & v, const std::string & key);

template
int appConfigurator::configUnused<unsigned long>( unsigned long & v, const std::string & key);

template
int appConfigurator::configUnused<unsigned long long>( unsigned long long & v, const std::string & key);

template
int appConfigurator::configUnused<float>( float & v, const std::string & key);

template
int appConfigurator::configUnused<double>( double & v, const std::string & key);

template
int appConfigurator::configUnused<long double>( long double & v, const std::string & key);

#ifdef HASQUAD
template
int appConfigurator::configUnused<__float128>( __float128 & v, const std::string & key);
#endif

template
int appConfigurator::configUnused<bool>( bool & v, const std::string & key);

template
int appConfigurator::configUnused<std::string>( std::string & v, const std::string & key);

//+++++++++++++++++++++++++++++++++++++++++++++++++++++

template
int appConfigurator::configUnused<char>( char & v, const std::string & section, const std::string & keyword);

template
int appConfigurator::configUnused<char16_t>( char16_t & v, const std::string & section, const std::string & keyword);

template
int appConfigurator::configUnused<char32_t>( char32_t & v, const std::string & section, const std::string & keyword);

template
int appConfigurator::configUnused<wchar_t>( wchar_t & v, const std::string & section, const std::string & keyword);

template
int appConfigurator::configUnused<signed char>( signed char & v, const std::string & section, const std::string & keyword);

template
int appConfigurator::configUnused<short>( short & v, const std::string & section, const std::string & keyword);

template
int appConfigurator::configUnused<int>( int & v, const std::string & section, const std::string & keyword);

template
int appConfigurator::configUnused<long>( long & v, const std::string & section, const std::string & keyword);

template
int appConfigurator::configUnused<long long>( long long & v, const std::string & section, const std::string & keyword);

template
int appConfigurator::configUnused<unsigned char>( unsigned char & v, const std::string & section, const std::string & keyword);

template
int appConfigurator::configUnused<unsigned short>( unsigned short & v, const std::string & section, const std::string & keyword);

template
int appConfigurator::configUnused<unsigned int>( unsigned int & v, const std::string & section, const std::string & keyword);

template
int appConfigurator::configUnused<unsigned long>( unsigned long & v, const std::string & section, const std::string & keyword);

template
int appConfigurator::configUnused<unsigned long long>( unsigned long long & v, const std::string & section, const std::string & keyword);

template
int appConfigurator::configUnused<float>( float & v, const std::string & section, const std::string & keyword);

template
int appConfigurator::configUnused<double>( double & v, const std::string & section, const std::string & keyword);

template
int appConfigurator::configUnused<long double>( long double & v, const std::string & section, const std::string & keyword);

#ifdef HASQUAD
template
int appConfigurator::configUnused<__float128>( __float128 & v, const std::string & section, const std::string & keyword);
#endif

template
int appConfigurator::configUnused<bool>( bool & v, const std::string & section, const std::string & keyword);

template
int appConfigurator::configUnused<std::string>( std::string & v, const std::string & section, const std::string & keyword);

//+++++++++++++++++++++++++++++++++++++

int appConfigurator::unusedSections( std::vector<std::string> & sections )
{
   sections.clear();
   
   //Wind through all the targets
   for(auto it = m_unusedConfigs.begin(); it != m_unusedConfigs.end(); ++it)
   {
      std::string sect = it->second.section;
      
      bool add = true;

      //Check if this section is already in the vector -- is there a std::algorithms way to do this?
      for(size_t i=0;i < sections.size(); ++i)
      {
         if(sections[i] == sect)
         {
            add = false;
            break;
         }
      }
      
      //Add it if it wasn't found.
      if(add) sections.push_back(sect);
   }
   
   return 0;
}

int appConfigurator::isSetUnused( const std::string & name )
{
   return isSet( name, m_unusedConfigs);
}

void writeConfigFile( const std::string & fname,                 
                      const std::vector<std::string> & sections, 
                      const std::vector<std::string> & keywords, 
                      const std::vector<std::string> & values    
                    )
{
   std::ofstream fout;
   
   fout.open(fname);
   
   std::string lastSection;
   for(size_t i=0; i< sections.size(); ++i)
   {
      std::string currSection = sections[i];
      
      if( currSection != lastSection && currSection != "")
      {
         fout << "\n[" << currSection << "]\n";
      }
         
      fout << keywords[i] << "=" << values[i] << "\n";
         
      lastSection = currSection;
   }
   
   fout.close();
   
   return;
}

} //namespace app
} //namespace mx

