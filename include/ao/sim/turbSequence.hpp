/** \file turbSequence.hpp
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Manage a turbulence sequence saved as FITS files.
  * \ingroup mxAO_sim_files
  *
  */

#ifndef turbSequence_hpp
#define turbSequence_hpp

#include <vector>
#include <string>

#include <Eigen/Dense>

#include "../../ioutils/fileUtils.hpp"
#include "../../ioutils/fits/fitsFile.hpp"

#include "wavefront.hpp"

namespace mx
{
namespace AO
{
namespace sim
{

template<typename _realT>
struct turbSequence
{
   typedef _realT realT;
   typedef Eigen::Array<realT, Eigen::Dynamic, Eigen::Dynamic> imageT;

   std::vector<std::string> _phaseFnames;
   std::vector<std::string> _ampFnames;

   int _size {0};

   int _nPerCube {0};

   int _files {0};

   int _frames {0};

   int _currFileNo {0};
   int _currFrameNo {0};

   realT _wfPS {1};
   realT _F0Photons {2e11};

   realT _starMag {0};

   realT _pixVal;

   imageT * _pupil;

   bool _loopClosed;

   improc::eigenCube<realT> _currPhase;
   improc::eigenCube<realT> _currAmp;

   bool _phaseOnly {false};

   turbSequence()
   {
   }

   /// Return the size of the turbulence screens
   /**
     * \returns _size
     */
   int size()
   {
      return _size;
   }

   /// Return the number of frames per cube
   /**
     * \returns _nPerCube
     */
   int nPerCube()
   {
      return _nPerCube;
   }

   /// Return the number of files
   /**
     * \returns _files
     */
   int files()
   {
      return _files;
   }

   /// Return the number of frames in each file
   /**
     * \returns _frames
     */
   int frames()
   {
      return _frames;
   }


   realT wfPS()
   {
      return _wfPS;
   }

   void wfPS(realT ps)
   {
      _wfPS = ps;
      calcPixVal();
   }

   realT F0Photons()
   {
      return _F0Photons;
   }

   void F0Photons(realT f0)
   {
      _F0Photons = f0;
      calcPixVal();
   }

   realT starMag()
   {
      return _starMag;
   }

   void starMag(realT sm)
   {
      _starMag = sm;
      calcPixVal();
   }


   realT pixVal()
   {
      return _pixVal;
   }

   void calcPixVal()
   {
      //This is the square root of photons/pixel
      _pixVal = sqrt(_F0Photons)*pow(10., -0.2*_starMag)*_wfPS; //pow(_wfPS,2);
   }


   /// Get the file names of the sequence assuming they are stored in a directory.
   /**
     * \param [in] dir is the directory to search for phase and amplitude files.
     * \param [in] max [optional] specifies the maximum number of files to include the file name list.
     *
     * \retval 0 on success.
     * \retval -1 on error.
     */
   int turbFnames(std::string dir, int max= 0);

   void openPhaseFrame(int fn);

   void nextPhase(wavefront<realT> &wf);

   void nextWF(wavefront<realT> & wf);


};

template<typename realT>
int turbSequence<realT>::turbFnames(std::string dir, int max)
{
   _phaseFnames = ioutils::getFileNames(dir, "", ".pha", ".fits");
   _ampFnames = ioutils::getFileNames(dir, "", ".amp", ".fits");


   if(_phaseFnames.size() == 0)
   {
      mxError("turbSequence", MXE_FILENOTFOUND, "No turbulent phase files found.");
      return -1;
   }

   if(_ampFnames.size() == 0)
   {
      std::cerr << "turbSequence: no turbulent amplitude files found.\n";

      _phaseOnly = true;
   }

   if(max > 0)
   {
      _phaseFnames.erase(_phaseFnames.begin()+max, _phaseFnames.end());
      if(_phaseOnly == false) _ampFnames.erase(_ampFnames.begin()+max, _ampFnames.end());
   }

   _files = _phaseFnames.size();

   openPhaseFrame(0);

    _size = _currPhase.rows();
   _nPerCube = _currPhase.planes();

   _frames = _nPerCube*_phaseFnames.size();

}

template<typename realT>
void turbSequence<realT>::openPhaseFrame(int fn)
{
   fits::fitsFile<float> ff;

   std::cout << _phaseFnames[fn] << "\n";

   ff.read(_phaseFnames[fn], _currPhase);
   ff.close();

   if(_phaseOnly == false)
   {
      ff.read(_ampFnames[fn], _currAmp);
      ff.close();
   }

   _currFileNo = fn;
   _currFrameNo = 0;
}

template<typename realT>
void turbSequence<realT>::nextPhase(wavefront<realT> &wf)
{

   int Npix = _pupil->sum();

   if(_currFrameNo > _nPerCube-1)
   {
      openPhaseFrame(_currFileNo + 1);
   }

   wf.phase.resize(_pupil->rows(),_pupil->cols());
   //mx::cutImage(wf.phase, _currPhase.image(_currFrameNo), _pupil->rows(), _pupil->cols());

   wf.phase = _currPhase.image(_currFrameNo).block( 0.5*( _currPhase.rows()-_pupil->rows()), 0.5*(_currPhase.cols() - _pupil->cols()), _pupil->rows(), _pupil->cols());


   wf.phase = (wf.phase - (wf.phase* (*_pupil)).sum()/Npix)* (*_pupil);

   //wf.amplitude.Constant(_pixVal);

   wf.amplitude = (*_pupil)*_pixVal;

   ++_currFrameNo;
}

template<typename realT>
void turbSequence<realT>::nextWF(wavefront<realT> & wf)
{

   int Npix = _pupil->sum();


   if(_currFrameNo > _nPerCube-1)
   {
      openPhaseFrame(_currFileNo + 1);
   }
   wf.phase.resize(_pupil->rows(),_pupil->cols());
   //mx::cutImage(wf.phase, _currPhase.image(_currFrameNo), _pupil->rows(), _pupil->cols());
   wf.phase = _currPhase.image(_currFrameNo).block( 0.5*( _currPhase.rows()-_pupil->rows()), 0.5*(_currPhase.cols() - _pupil->cols()), _pupil->rows(), _pupil->cols());

   wf.phase = (wf.phase - (wf.phase* (*_pupil)).sum()/Npix)* (*_pupil);

   if(_phaseOnly)
   {
      wf.amplitude = _pixVal*(*_pupil);
   }
   else
   {
      wf.amplitude.resize(_pupil->rows(),_pupil->cols());
      //mx::cutImage(wf.amplitude, _currAmp.image(_currFrameNo), _pupil->rows(), _pupil->cols());
      wf.amplitude = _currAmp.image(_currFrameNo).block( 0.5*( _currAmp.rows()-_pupil->rows()), 0.5*(_currAmp.cols() - _pupil->cols()), _pupil->rows(), _pupil->cols());

      wf.amplitude *= _pixVal*(*_pupil);
   }
   _currFrameNo++;
}

} //namespace sim
} //namespace AO
} //namespace mx

#endif //turbSequence_hpp
