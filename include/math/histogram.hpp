

template<typename realT>
class histogramUniform
{
public:

   realT _min;
   realT _max;
   realT _width;

   std::vector<realT> _freqs;
   
   void setup( realT mn,
               realT mx,
               realT w )
   {
      _min = mn;
      _max = mx;
      _width = w;
      
      reset();
   }
   
   void reset()
   {
      _freqs.resize( (_max - _min)/_width + 1, 0);
   }
   
   void accum( const realT & val )
   {
      int i = (val - _min) / _width;
      if( i < 0) i = 0;
      if( i >= _freqs.size()) i = _freqs.size()-1;
      
      
      ++_freqs[i];
   }
   
   void accum( const std::vector<realT> & vals )
   {
      for(int i=0; i< vals.size(); ++i) accum(vals[i]);
   }
   
   realT freq(int i /**< [in] the bin number */)
   {
      return _freqs[i];
   }
   
   int bins()
   {
      return _freqs.size();
   }
   
   realT binLeft(int i /**< [in] the bin number */)
   {
      return _min + i*_width;
   }
   
   realT binMid(int i /**< [in] the bin number */)
   {
      return _min + i*_width + 0.5*_width;
   }
   
   ///Get the right edge of the i-th bin.
   realT binRight(int i /**< [in] the bin number */)
   {
      return _min + i*_width + 1.0*_width;
   }
   
   void normalize( int excludeTop = 0 )
   {
      realT sum = 0;
      
      for(int i=0; i< _freqs.size()-excludeTop; ++i)
      {
         sum += _freqs[i];
      }
      
      for(int i=0; i< _freqs.size(); ++i)
      {
         _freqs[i] /= (sum*_width);
      }
   }
   
};
