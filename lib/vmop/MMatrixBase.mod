/** \file MMatrixBase
  * \brief Declares and defines the base class for all MMatrices.
  *
  * \author Jared R. Males (jaredmales@gmail.com)
  *
  */

#ifndef __mx_MMatrixBase__
#define __mx_MMatrixBase__

#include <cstdarg>
#include <cstdlib>
#include <complex>

#include <iostream>

#include "MMConfig.h"
//#include "MMatrixIterator"
#include "MMaccessor"
#include "MMfunctional"
#include "MMatrixBound"
#include "MMatrixSlice"

namespace mx
{

///This is the base class for all MMatrix types.
/** \ingroup vmop 
  * This base class provides most of the underlying functionality of the MMatrix system.  It is not
  * intended to be used by itself, rather should be a base class.
  *
  * The dimension of an MMatrix is hard coded via the template parameter \p _dim, which allows the compiler
  * to optimize loops and allows specialization by dimension for further optimization.
  * 
  * \tparam dataT the data type held by the mmatrix.
  * \tparam _dim is the dimension of the mmatrix.
  * \tparam absT is the data type returned by the absolute value operation.  This is to handle complex numbers.
  */
template<class dataT, size_t _dim, class absT>
class MMatrixBase
{

public:
   
   typedef bool is_mmatrix; ///<For use by the is_mmatrix test
   
   typedef dataT data_type; ///<The type stored by the matrix
   typedef absT  abs_type; ///<The type returned by the absolute value operation
   
   static const size_t dimension = _dim; ///<The dimension of the matrix.
      
public:
   /** @name Construction and Destruction
     */
   //@{
   
   ///Default constructor. 
   /** Initializes the non-array members but does not allocate storage.
     */
   MMatrixBase();
   
   ///Construct and allocate the data array.
   /** The size of each dimension is specified in an array.
     * 
     * \param sizes is an array of length _dim specifying the size of each dimension
     */
   explicit MMatrixBase(const size_t *sizes);
            
   ///Copy constructor.
   /**
     * 
     * \param m is an MMatrix to copy
     */
   MMatrixBase( const MMatrixBase & m );

   ///Move constructor.
   /** If \p m owns its storage, it is just re-pointed and ownership is adjusted.
     * 
     * \param m is an MMatrix r-value to copy
     */
   MMatrixBase( MMatrixBase && m );

   ///Destructor.  
   /** If owned, the _data storage is freed.
     */
   virtual ~MMatrixBase();

   //@}


   
protected:
   
   
   /** @name Set-By-User Variables
     * \brief Variables which are set explicitly by a call to a member method.
     */
   //@{
      
   size_t  _size[_dim];///<The raw size of each dimension
   
   size_t _offset[_dim];///<The offsets from the beginning of each dimension
   
   size_t  _stride[_dim]; ///<The stride of each dimension.  Strides are >= 1.
   
   //@}
   
   /** @name Calculated Variables
     * \brief Variables which are calculated by a call to a member method.
     */
   //@{
   
   dataT * _raw_data; ///<The raw data storage array
   size_t _tot_size; ///<The total size of the storage array, independent of stride.
   bool _owner; ///<ownership status of the storage array

   dataT * _data; ///<The current start of the data, affected by offsets
   
   size_t _size_prod[_dim]; ///<The size of 1 increment of each dimension, which is the product of sizes for the dimensions below
   
   bool _isvector; ///<0 if all strides == 1, 1 otherwise.  Allows for direct indexing if no strides > 1 are used.
   
   size_t _length[_dim]; ///<The apparent length of each dimension, taking stride into account.
   
   size_t _tot_length; ///<The total apparent length of the matrix, taking strides into account

   //@}
   
protected:

   /** @name Low-level allocation and initialization
     * \brief These set initial values, allocate, and free memory.
     */
   //@{
      
   ///Initializes _data and _owner
   /** Must be called on construction to ensure free() works properly.
     */
   void initialize();
   
   ///Resizes the storage array.
   /** Sets size based on  _tot_size.
     */
   virtual void allocate();
   
   ///Free the storage array.
   /** Only frees if owned.  
     */
   virtual void free();
   
   //@}
   
public:

   /** @name Resizing the matrix
     */
   //@{
         
   ///Resize the storage array
   /** If size does not change, and the block is already owned, this does nothing.
     * 
     * \param sizes is an array of length _dim
     */  
   virtual void resize(const size_t * sizes);
      
   //@}
   
   /** @name Copying
     */  
   //@{
     
   /// Copy the meta-data from one MMatrixBase to another.
   /** This only copies the meta-data, does not copy the storage array.
     * 
     * \param m is the matrix to copy
     */
   virtual void copy( const MMatrixBase & m );
  
   
   //@}
   
             
public:

   /** @name Member Access
     */
   //@{
   
   ///Get the dimension of the matrix.
   /** 
     * \returns the number of dimensions of the matrix
     */
   size_t dim() const;

   ///Get the total storage size of the matrix.
   /** 
     * \returns the total storage size of the matrix
     */
   size_t total_size() const;

   ///Returns the size array.
   /** 
     * \returns a pointer to an array of size _dim
     */
   const size_t * size() const;
   
   ///Returns the storage size of the specified dimension, does not take stride into account.
   /** 
     * \param d is the dimension desired
     * 
     * \returns the storage size of dimension \p d
     */
   size_t size(const size_t d) const;

   ///Returns the size product of the specified dimension.
   /** 
     * \param d is the dimension desired
     * 
     * \returns the size product of dimension \p d;
     */
   size_t size_prod(size_t d) const;
      
   ///Returns the offset of the specified dimension.
   /** 
     * \param d is the dimension desired
     * 
     * \returns the offset of dimension \p d;
     */
   size_t offset(size_t d) const;
   
   ///Returns the stride of the specified dimension.
   /** 
     * \param d is the dimension desired
     * 
     * \returns the stride of dimension \p d
     */
   size_t stride(const size_t d) const;

   ///Check whether all dimensions have stride==1, in which case the matrix can be treated like a vector
   /** This allows using vectorized loops instead of using tuple indexing.  Whether or not dimensions are
     * have stride is determined once, upon allocation, so this is just accessing a constant.
     *
     * \returns true if all dimensions have stride==1, false otherwise
     */ 
   bool vector() const;
   
   ///Get the total apparent length of the matrix.
   /** 
     * \returns the total apparent length of the matrix
     */
   size_t total_length() const;
  
   ///Returns the length array.
   /** 
     * \returns a pointer to an array of size _dim
     */
   const size_t * length() const;
   
   ///Returns the length of the specified dimension, taking stride into account.
   /** 
     * \param d is the dimension desired
     * 
     * \returns the length of dimension \p d
     */
   size_t length(const size_t d) const;
   
   ///Get the address of the raw storage array
   /**
     * \returns the address of the raw storage array
     */
   dataT * raw_data() const;
   
   ///Get the address of the storage array including the offsets
   /**
     * \returns the address of the storage array with the offsets
     */
   dataT * data() const;
      
   ///Get the ownership status of the data array
   /** 
     *\returns true if this object owns the data array, false otherwise
     */
   bool owner() const;
   
   //@}
      
};

template<class dataT, size_t _dim, class absT>
MMatrixBase<dataT, _dim, absT>::MMatrixBase()
{
   MX_MM_DEBUG_TRACE("MMatrixBase::default ctor");
   
   initialize();
   
}

template<class dataT, size_t _dim, class absT>
MMatrixBase<dataT, _dim, absT>::MMatrixBase( const size_t * sizes )
{
   MX_MM_DEBUG_TRACE("MMatrixBase::ctor size_t*");
   
   initialize();
   
   resize(sizes);
}

template<class dataT, size_t _dim, class absT>
MMatrixBase<dataT, _dim, absT>::MMatrixBase( const MMatrixBase & m )
{
   MX_MM_DEBUG_TRACE("MMatrixBase::copy ctor &");
   
   initialize();

   copy(m);
}

template<class dataT, size_t _dim, class absT>
MMatrixBase<dataT, _dim, absT>::MMatrixBase(MMatrixBase && m)
{
   MX_MM_DEBUG_TRACE("MMatrixBase::copy ctor &&");
   
   initialize();
   
   copy(m);
}

template<class dataT, size_t _dim, class absT>
MMatrixBase<dataT, _dim, absT>::~MMatrixBase()
{
   MX_MM_DEBUG_TRACE("MMatrixBase::dtor");
   
   free();  
}

template<class dataT, size_t _dim, class absT>
void MMatrixBase<dataT, _dim, absT>::initialize()
{  
   MX_MM_DEBUG_TRACE("MMatrixBase::initialize");
   
   _raw_data = 0;
   
   _data = 0;
   
   _owner = false;
   
   _isvector = 0;
   
   _tot_size=0;
   
   _tot_length=0;
   
}//void initialize()

template<class dataT, size_t _dim,  class absT> 
void MMatrixBase<dataT, _dim, absT>::allocate()
{
   MX_MM_DEBUG_TRACE("MMatrixBase::allocate");
   
   //Always check if already allocated.
   free();

   _raw_data = new dataT[_tot_size];
   _data = _raw_data;
   
   _owner = true;
   
}//void allocate()

template<class dataT, size_t _dim,  class absT> 
void MMatrixBase<dataT, _dim, absT>::free()
{
   MX_MM_DEBUG_TRACE("MMatrixBase::free");
   
   if(_data && _owner)
   {
      delete[] _raw_data;
      _owner = false;
      _data = 0;
   }
}


template<class dataT, size_t _dim,  class absT>
void MMatrixBase<dataT, _dim, absT>::resize(const size_t * sizes)
{
   MX_MM_DEBUG_TRACE("MMatrixBase::resize");
   
   size_t old_size = _tot_size;
 
   //Update meta data;
   _tot_size = sizes[0];
   
   _offset[0] = 0;
   _stride[0] = 1;
   _length[0] = sizes[0];
   
   _size_prod[0] = 1;
   
   for(size_t i=1;i<_dim;i++)
   {
      _size[i] = sizes[i];
      
      _offset[i] = 0;
      _stride[i] = 1;
      _length[i] = sizes[i];
      
      _tot_size *= sizes[i];
      
      _size_prod[i] = 1;
   }

   for(size_t i=0; i < _dim; i++)
   {
      for(size_t j=i+1; j<_dim;j++)
      {
         _size_prod[i] *= _size[j];
      }
   }
   
   _tot_length = _tot_size;
   
   //And now allocate the storage array, if necessary
   if(_tot_size != old_size || _owner == false)
   {
      allocate();
   }
   
}//void resize(size_t *)

template<class dataT, size_t _dim, class absT>
MMatrixBase<dataT, _dim, absT> & MMatrixBase<dataT, _dim, absT>::copy( const MMatrixBase<dataT, _dim, absT> & m )
{
   MX_MM_DEBUG_TRACE("MMatrixBase::copy");
      
   _raw_data = m._raw_data;
   _tot_size = m._tot_size;
   
   _owner = false;
   
   _data = m._data;
   
   for(size_t i=0;i<_dim;i++)
   {
      _size[i] = m._size[i];
      _offset[i] = m._offset[i];
      _stride[i] =  m._stride[i];
      _isvector = m._isvector;
      _length[i] = m._length[i];                
   }
      
   _tot_length = m._tot_length;
         
   return *this;
}//MMatrixBase & copy(const MMatrixBase &)

template<class dataT, size_t _dim, class absT>
size_t MMatrixBase<dataT, _dim, absT>::total_size() const
{ 
   return _tot_size;
}

template<class dataT, size_t _dim, class absT>
const size_t * MMatrixBase<dataT, _dim, absT>::size() const
{ 
   return _size; 
}
   
template<class dataT, size_t _dim, class absT>
size_t MMatrixBase<dataT, _dim, absT>::size(const size_t d) const 
{ 
   return _size[d];
}
 
template<class dataT, size_t _dim, class absT>
size_t MMatrixBase<dataT, _dim, absT>::size_prod(const size_t d) const 
{ 
   return _size_prod[d];
}

template<class dataT, size_t _dim, class absT>
size_t MMatrixBase<dataT, _dim, absT>::offset(size_t d) const
{
   return _offset[d];
}
   
template<class dataT, size_t _dim, class absT>
size_t MMatrixBase<dataT, _dim, absT>::stride(const size_t d) const
{ 
   return _stride[d];
}

template<class dataT, size_t _dim, class absT>
bool MMatrixBase<dataT, _dim, absT>::vector() const 
{
   return _isvector;
}

template<class dataT, size_t _dim, class absT>
size_t MMatrixBase<dataT, _dim, absT>::total_length() const
{ 
   return _tot_length;
}
   
template<class dataT, size_t _dim, class absT>
const size_t * MMatrixBase<dataT, _dim, absT>::length() const
{ 
   return _length; 
}
   
template<class dataT, size_t _dim, class absT>
size_t MMatrixBase<dataT, _dim, absT>::length(const size_t d) const 
{ 
   return _length[d];
}

template<class dataT, size_t _dim, class absT>
dataT * MMatrixBase<dataT, _dim, absT>::raw_data() const
{ 
   return _raw_data;
}

template<class dataT, size_t _dim, class absT>
dataT * MMatrixBase<dataT, _dim, absT>::data() const
{ 
   return _data;
}
 
template<class dataT, size_t _dim, class absT>
bool MMatrixBase<dataT, _dim, absT>::owner() const
{ 
   return _owner;
}


} //namespace mx
#endif // __mx_mmatrix_base__



   