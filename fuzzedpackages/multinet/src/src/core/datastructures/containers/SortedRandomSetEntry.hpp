/**
 * An entry in a SortedRandomSet
 *
 * History:
 * - 2018.01.01 file adapted from version 1.0 of the multinet library
 */

#ifndef UU_CORE_DATASTRUCTURES_CONTAINERS_SORTEDRANDOMSETENTRY_H_
#define UU_CORE_DATASTRUCTURES_CONTAINERS_SORTEDRANDOMSETENTRY_H_

#include <memory>
#include <vector>

namespace uu {
namespace core {

template <class ELEMENT_TYPE>
class
    SortedRandomSet;

//template <class ELEMENT_TYPE>
//class
//    SortedRandomSetEntry;

//template <class ELEMENT_TYPE> using
//SortedRandomSetEntrySharedPtr =
//    std::shared_ptr<SortedRandomSetEntry<ELEMENT_TYPE> >;

//template <class ELEMENT_TYPE> using
//constSortedRandomSetEntrySharedPtr =
//    std::shared_ptr<const SortedRandomSetEntry<ELEMENT_TYPE> >;

/**
 * An entry in a sorted set, which is implemented as a skip list.
 */
template <class ELEMENT_TYPE>
class
    SortedRandomSetEntry
{

    friend class SortedRandomSet<ELEMENT_TYPE>;

  private:

    /** The object corresponding to this entry. */
    ELEMENT_TYPE value;

    /** An array of pointers to the next entry, for each level in the skip list. */
    std::vector<std::shared_ptr<SortedRandomSetEntry<ELEMENT_TYPE>>> forward;

    /** The number of entries before the next entry on each level, used for positional access. */
    std::vector<int> link_length;

  public:

    /**
     * Constructor, to be used for entries not corresponding to any value (that is, the header).
     * @param level height of the entry in the skip list
     */
    SortedRandomSetEntry(
        int level
    );

    /**
     * Constructor.
     * @param level height of the entry in the skip list
     * @param obj the object corresponding to this entry
     */
    SortedRandomSetEntry(
        int level,
        const ELEMENT_TYPE& obj
    );

    /**
     * Constructor.
     * @param level height of the entry in the skip list
     * @param obj the object corresponding to this entry
     */
    SortedRandomSetEntry(
        int level,
        ELEMENT_TYPE&& obj
    );

    /**
     * This function is used to increase the level of the header.
     * @param skipped_entries number of skipped entries (to be set to the number of entries)
     */
    void
    increment(
        long skipped_entries
    );

};


/* TEMPLATE CODE */

template <class ELEMENT_TYPE>
SortedRandomSetEntry<ELEMENT_TYPE>::
SortedRandomSetEntry(
    int level
)
{
    forward.resize(level+1);
    link_length.resize(level+1);
}

template <class ELEMENT_TYPE>
SortedRandomSetEntry<ELEMENT_TYPE>::
SortedRandomSetEntry(
    int level,
    const ELEMENT_TYPE& value
)
{
    forward.resize(level+1);
    link_length.resize(level+1);
    this->value = value;
}

template <class ELEMENT_TYPE>
SortedRandomSetEntry<ELEMENT_TYPE>::
SortedRandomSetEntry(
    int level,
    ELEMENT_TYPE&& value
)
{
    forward.resize(level+1);
    link_length.resize(level+1);
    this->value = std::move(value);
}

template <class ELEMENT_TYPE>
void
SortedRandomSetEntry<ELEMENT_TYPE>::
increment(
    long skipped_entries
)
{
    int current_size = forward.size();
    forward.resize(current_size+1,nullptr);
    link_length.resize(current_size+1,skipped_entries);
}

}
}

#endif
