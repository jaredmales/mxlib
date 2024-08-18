/** \file changeable.hpp
 * \brief A simple class to track member data changes
 *
 * \author Jared R. Males (jaredmales@gmail.com)
 *
 * \ingroup utils_files
 *
 */

#ifndef mx_base_changeable_hpp
#define mx_base_changeable_hpp

namespace mx
{
namespace base
{

/// A simple class to track member data changes.
/** Maintains a monotonic counter that derived classes increment anytime a member
 * datum is changed.  This allows for, say, re-initializing automatically on a subsequent
 * call to a function.
 *
 * \ingroup base
 */
template <class _derivedT>
class changeable
{

  public:
    typedef _derivedT derivedT;

    /// The integer type of the counter.
    typedef uint64_t changeT;

  private:
    /// The counter itself.
    changeT m_change{ 0 };

    /// A marker for last change, set by \ref setChangePoint
    changeT m_changePoint{ std::numeric_limits<changeT>::max() };

  public:
    /// Increment the counter
    /** Call this function anytime the derived class changes something important.
     *
     */
    void changed()
    {
        ++m_change;
    }

    /// Get the value of the counter
    /** Get this at a checkpoint and compare it later to decide if any
     * action should be taken.
     */
    changeT change()
    {
        return m_change;
    }

    void setChangePoint()
    {
        m_changePoint = m_change;
    }

    bool isChanged()
    {
        return ( m_change != m_changePoint );
    }
};

} // namespace base
} // namespace mx

#endif // mx_base_changeable_hpp
