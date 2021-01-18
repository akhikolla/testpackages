/**
 * A bag supporting both sorted and random access.
 *
 * A sorted random bag is a class used to store a bag of objects that can be accessed:
 * 1. by key in (average) log time.
 * 2. by index (position) in constant time.
 * 3. by const iterating over its elements.
 *
 * History:
 * - 2019.07.28 file created
 */


#ifndef UU_CORE_DATASTRUCTURES_CONTAINERS_SORTEDRANDOMBAG_H_
#define UU_CORE_DATASTRUCTURES_CONTAINERS_SORTEDRANDOMBAG_H_

#include <functional>
#include <memory>
#include <vector>
#include "core/datastructures/containers/DefaultLT.hpp"
#include "core/datastructures/containers/DefaultEQ.hpp"
#include "core/datastructures/containers/SortedRandomBagEntry.hpp"
#include "core/utils/random.hpp"
#include "core/exceptions/ElementNotFoundException.hpp"

namespace uu {
namespace core {


/**
 * A sorted bag is a class used to store a set of components that:
 * 1. can be accessed by id in log time.
 * 2. can be accessed by index (position) in the set in constant time.
 * 3. can be accessed using a const iterator.
 *
 * A sorted bag is implemented as a skip list.
 */
template <class ELEMENT_TYPE>
class SortedRandomBag
{
  private:
    float P = 0.5;

    std::shared_ptr<SortedRandomBagEntry<ELEMENT_TYPE>> header;

    /* Number of entries for which the sorted set is optimized. */
    size_t capacity = 1;

    /* Current number of entries, that is, distinct objects. */
    size_t num_entries = 0;

    /* Current number of objects. */
    size_t num_objects = 0;

    /* Maximum level */
    size_t MAX_LEVEL = 0;

    /* Current maximum level in use. */
    size_t level;

  public:

    /**
     * Creates a sorted set.
     */
    SortedRandomBag();

    /**
     * Destructs a sorted set, making sure that this does not result in a long
     * sequence of nested calls to entry destructors that would fill in the stack
     */
    ~SortedRandomBag();

    /**
     * Creates a sorted set optimized to store a pre-defined number of entries.
     * @param start_capacity the initial capacity for which the sorted set is optimized
     */
    SortedRandomBag(
        size_t start_capacity
    );

    /** Iterator over the objects in this collection. */
    class
        iterator
    {
        typedef std::forward_iterator_tag iterator_category;

      public:

        iterator();

        /** Returns an iterator pointing at the input object */
        iterator(
            std::shared_ptr<SortedRandomBagEntry<ELEMENT_TYPE>>
        );

        /** Return the object pointed by this iterator */
        const ELEMENT_TYPE&
        operator*(
        );

        /** Moves the iterator to the next object in the collection (prefix) */
        iterator
        operator++(
        );

        /** Moves the iterator to the next object in the collection (postfix) */
        iterator
        operator++(
            int
        );

        /** Checks if this iterator equals the input one */
        bool
        operator==(
            const SortedRandomBag<ELEMENT_TYPE>::iterator& rhs
        );

        /** Checks if this iterator differs from the input one */
        bool
        operator!=(
            const SortedRandomBag<ELEMENT_TYPE>::iterator& rhs
        );

      private:

        /** Entry currently pointed to by this iterator */
        std::shared_ptr<SortedRandomBagEntry<ELEMENT_TYPE>> current;

    };

    /** Returns an iterator to the first object in the collection */
    SortedRandomBag<ELEMENT_TYPE>::
    iterator
    begin(
    ) const;

    /** Returns an iterator after the last object in the collection */
    SortedRandomBag<ELEMENT_TYPE>::
    iterator
    end(
    ) const;

    /** Returns the number of objects in the collection */
    size_t
    size(
    ) const;

    /** Returns the number of distinct objects in the collection */
    size_t
    dsize(
    ) const;

    /** Returns true if an object with the input id is present in the collection */
    template <
        typename KEY = ELEMENT_TYPE,
        typename LT = DefaultLT<ELEMENT_TYPE, KEY>,
        typename EQ = DefaultEQ<ELEMENT_TYPE, KEY>
        >
    bool
    contains(
        const KEY
    ) const;

    /** Returns the position of the input value in the collection, or -1 */
    template <
        typename KEY = ELEMENT_TYPE,
        typename LT = DefaultLT<ELEMENT_TYPE, KEY>,
        typename EQ = DefaultEQ<ELEMENT_TYPE, KEY>
        >
    int
    index_of(
        const KEY
    ) const;

    /** Returns the object at the given position in the collection.
     * @throw ElementNotFoundException if the index is outside the bounds on the set
     */
    const ELEMENT_TYPE&
    at(
        size_t
    ) const;

    /** Returns a random object, uniform probability */
    const ELEMENT_TYPE&
    get_at_random(
    ) const;

    /**
     * Inserts a new object in the collection, or updates its count.
     * @return true if KEY was not already present, false otherwise
     * (in which case the object is updated with the new value)
     */
    bool
    add(
        const ELEMENT_TYPE&
    );

    /**
     * Inserts a new object in the collection, or updates its count.
     * @return true if KEY was not already present, false otherwise
     * (in which case the object is updated with the new value)
     */
    bool
    add(
        ELEMENT_TYPE&&
    );


    /**
     * Reduces the count of the input object, or removes it from the collection
     * if the count is 1.
     * @return true if the object is erased from the collection, false if the object was not present.
     */
    template <
        typename KEY = ELEMENT_TYPE,
        typename LT = DefaultLT<ELEMENT_TYPE, KEY>,
        typename EQ = DefaultEQ<ELEMENT_TYPE, KEY>
        >
    bool
    erase(
        KEY
    );

};

/* TEMPLATE CODE */

template <class ELEMENT_TYPE>
SortedRandomBag<ELEMENT_TYPE>::
SortedRandomBag(
)
{
    header = std::make_shared<SortedRandomBagEntry<ELEMENT_TYPE> >(MAX_LEVEL);
    level = 0;
}

template <class ELEMENT_TYPE>
SortedRandomBag<ELEMENT_TYPE>::
~SortedRandomBag()
{
    std::vector<std::shared_ptr<SortedRandomBagEntry<ELEMENT_TYPE>> > tmp(size());

    std::shared_ptr<SortedRandomBagEntry<ELEMENT_TYPE>> current = header;

    std::shared_ptr<SortedRandomBagEntry<ELEMENT_TYPE>> next = current->forward.at(0);

    tmp.push_back(current);

    while (next)
    {
        tmp.push_back(next);

        for (size_t i=0; i<current->forward.size(); i++)
        {
            current->forward[i] = nullptr;
        }

        current = next;
        next = current->forward.at(0);

    }
}

template <class ELEMENT_TYPE>
SortedRandomBag<ELEMENT_TYPE>::
SortedRandomBag(
    size_t start_capacity
)
{
    capacity = start_capacity;
    MAX_LEVEL = std::ceil(std::log2(capacity));
    header = std::make_shared<SortedRandomBagEntry<ELEMENT_TYPE> >(MAX_LEVEL);
    level = 0;
}

template <class ELEMENT_TYPE>
typename SortedRandomBag<ELEMENT_TYPE>::iterator
SortedRandomBag<ELEMENT_TYPE>::
begin(
) const
{
    return iterator(header->forward[0]);
}

template <class ELEMENT_TYPE>
typename SortedRandomBag<ELEMENT_TYPE>::iterator
SortedRandomBag<ELEMENT_TYPE>::
end(
) const
{
    return iterator(nullptr);
}

template <class ELEMENT_TYPE>
const ELEMENT_TYPE&
SortedRandomBag<ELEMENT_TYPE>::iterator::operator*(
)
{
    return current->value;
}

template <class ELEMENT_TYPE>
SortedRandomBag<ELEMENT_TYPE>::iterator::
iterator(
    std::shared_ptr<SortedRandomBagEntry<ELEMENT_TYPE>> iter
) :
    current(iter)
{
}

template <class ELEMENT_TYPE>
typename SortedRandomBag<ELEMENT_TYPE>::iterator
SortedRandomBag<ELEMENT_TYPE>::iterator::
operator++(
)
{
    // PREFIX
    current=current->forward[0];
    return *this;
}

template <class ELEMENT_TYPE>
typename SortedRandomBag<ELEMENT_TYPE>::iterator
SortedRandomBag<ELEMENT_TYPE>::iterator::
operator++(
    int
)
{
    // POSTFIX
    SortedRandomBag<ELEMENT_TYPE>::iterator tmp(current);
    current=current->forward[0];
    return tmp;
}

template <class ELEMENT_TYPE>
bool
SortedRandomBag<ELEMENT_TYPE>::iterator::
operator==(
    const SortedRandomBag<ELEMENT_TYPE>::iterator& rhs
)
{
    return current == rhs.current;
}

template <class ELEMENT_TYPE>
bool
SortedRandomBag<ELEMENT_TYPE>::iterator::
operator!=(
    const SortedRandomBag<ELEMENT_TYPE>::iterator& rhs
)
{
    return current != rhs.current;
}

template <class ELEMENT_TYPE>
size_t
SortedRandomBag<ELEMENT_TYPE>::size(
) const
{
    return num_objects;
}

template <class ELEMENT_TYPE>
size_t
SortedRandomBag<ELEMENT_TYPE>::dsize(
) const
{
    return num_entries;
}

template <class ELEMENT_TYPE>
template <typename KEY, typename LT, typename EQ>
bool
SortedRandomBag<ELEMENT_TYPE>::
contains(
    const KEY search_value
) const
{
    LT lt = LT();
    EQ eq = EQ();
    std::shared_ptr<const SortedRandomBagEntry<ELEMENT_TYPE>> x = header;

    for (int i = level; i >= 0; i--)
    {
        while (x->forward[i] != nullptr && lt(x->forward[i]->value, search_value))
        {
            x = x->forward[i];
        }
    }

    x = x->forward[0];
    return x != nullptr && eq(x->value, search_value);
}

template <class ELEMENT_TYPE>
template <typename KEY, typename LT, typename EQ>
int
SortedRandomBag<ELEMENT_TYPE>::
index_of(
    const KEY search_value
) const
{
    LT lt = LT();
    EQ eq = EQ();
    std::shared_ptr<const SortedRandomBagEntry<ELEMENT_TYPE>> x = header;
    long so_far=0;

    for (int i = level; i >= 0; i--)
    {
        while (x->forward[i] != nullptr && lt(x->forward[i]->value, search_value))
        {
            so_far+= x->link_length[i];
            x = x->forward[i];
        }
    }

    so_far+= x->link_length[0];
    x = x->forward[0];

    if (x != nullptr && eq(x->value, search_value))
    {
        return so_far-1;
    }

    else
    {
        return -1;
    }
}

template <class ELEMENT_TYPE>
const ELEMENT_TYPE&
SortedRandomBag<ELEMENT_TYPE>::
at(
    size_t pos
) const
{
    if (pos < 0 || pos >= num_entries)
    {
        throw ElementNotFoundException("Index out of bounds");
    }

    std::shared_ptr<const SortedRandomBagEntry<ELEMENT_TYPE>> x = header;
    size_t so_far=0;

    for (int i = level; i >= 0; i--)
    {
        while (x->forward[i] != nullptr && x->link_length[i] + so_far <= pos + 1)
        {
            so_far += x->link_length[i];
            x = x->forward[i];
        }
    }

    return x->value;
}

template <class ELEMENT_TYPE>
const ELEMENT_TYPE&
SortedRandomBag<ELEMENT_TYPE>::
get_at_random(
) const
{
    return at(getRandomInt(dsize()));
}

template <class ELEMENT_TYPE>
bool
SortedRandomBag<ELEMENT_TYPE>::
add(
    const ELEMENT_TYPE& value
)
{
    // the number of objects is always incremented
    num_objects++;

    std::shared_ptr<SortedRandomBagEntry<ELEMENT_TYPE>> x = header;
    std::vector<std::shared_ptr<SortedRandomBagEntry<ELEMENT_TYPE>> > update;
    update.resize(level+1);
    std::vector<size_t> skipped_positions_per_level;
    size_t skipped_positions = 0;
    skipped_positions_per_level.resize(level+1,0);

    for (int i = level; i >= 0; i--)
    {
        skipped_positions_per_level[i] = skipped_positions;

        while (x->forward[i] != nullptr && x->forward[i]->value < value)
        {
            skipped_positions_per_level[i] += x->link_length[i];
            skipped_positions += x->link_length[i];
            x = x->forward[i];
        }

        update[i] = x;
    }

    x = x->forward[0];


    if (x == nullptr || x->value != value)
    {
        num_entries++;

        if (num_entries>capacity)
        {
            // resize the sorted list
            capacity *= 2;
            MAX_LEVEL++;
            header->increment(num_entries);
        }

        size_t lvl = random_level(MAX_LEVEL,P);

        if (lvl > level)
        {
            update.resize(lvl+1);
            skipped_positions_per_level.resize(lvl+1,0);

            for (size_t i = level + 1; i <= lvl; i++)
            {
                update[i] = header;
                update[i]->link_length[i] = num_entries;
            }

            level = lvl;
        }

        x = std::make_shared<SortedRandomBagEntry<ELEMENT_TYPE> >(lvl, value);

        for (size_t i = 0; i <= lvl; i++)
        {
            size_t offset = skipped_positions-skipped_positions_per_level[i];

            x->forward[i] = update[i]->forward[i];

            if (update[i]->forward[i]==nullptr)
            {
                x->link_length[i] = num_entries - skipped_positions;
            }

            else
            {
                x->link_length[i] = update[i]->link_length[i]-offset;
            }

            update[i]->forward[i] = x;
            update[i]->link_length[i] = offset+1;
        }

        for (size_t i = lvl+1; i <= level; i++)
        {
            update[i]->link_length[i]++;
        }
    }

    else
    {
        x->value = value;
        x->counter++;
    }

    return true;
}


template <class ELEMENT_TYPE>
bool
SortedRandomBag<ELEMENT_TYPE>::
add(
    ELEMENT_TYPE&& value
)
{
    // the number of objects is always incremented
    num_objects++;

    std::shared_ptr<SortedRandomBagEntry<ELEMENT_TYPE>> x = header;
    std::vector<std::shared_ptr<SortedRandomBagEntry<ELEMENT_TYPE>> > update;
    update.resize(level+1);
    std::vector<size_t> skipped_positions_per_level;
    size_t skipped_positions = 0;
    skipped_positions_per_level.resize(level+1,0);

    for (int i = level; i >= 0; i--)
    {
        skipped_positions_per_level[i] = skipped_positions;

        while (x->forward[i] != nullptr && x->forward[i]->value < value)
        {
            skipped_positions_per_level[i] += x->link_length[i];
            skipped_positions += x->link_length[i];
            x = x->forward[i];
        }

        update[i] = x;
    }

    x = x->forward[0];


    if (x == nullptr || x->value != value)
    {
        num_entries++;
        //std::cout << num_entries << std::endl;

        if (num_entries>capacity)
        {
            // resize the sorted list
            capacity *= 2;
            MAX_LEVEL++;
            header->increment(num_entries);
        }

        size_t lvl = random_level(MAX_LEVEL,P);

        if (lvl > level)
        {
            update.resize(lvl+1);
            skipped_positions_per_level.resize(lvl+1,0);

            for (size_t i = level + 1; i <= lvl; i++)
            {
                update[i] = header;
                update[i]->link_length[i] = num_entries;
            }

            level = lvl;
        }

        x = std::make_shared<SortedRandomBagEntry<ELEMENT_TYPE> >(lvl, std::move(value));

        for (size_t i = 0; i <= lvl; i++)
        {
            size_t offset = skipped_positions-skipped_positions_per_level[i];

            x->forward[i] = update[i]->forward[i];

            if (update[i]->forward[i]==nullptr)
            {
                x->link_length[i] = num_entries - skipped_positions;
            }

            else
            {
                x->link_length[i] = update[i]->link_length[i]-offset;
            }

            update[i]->forward[i] = x;
            update[i]->link_length[i] = offset+1;
        }

        for (size_t i = lvl+1; i <= level; i++)
        {
            update[i]->link_length[i]++;
        }
    }

    else
    {
        x->value = std::move(value);
        x->counter++;
    }

    return true;
}


template <class ELEMENT_TYPE>
template <typename KEY, typename LT, typename EQ>
bool
SortedRandomBag<ELEMENT_TYPE>::
erase(
    KEY value
)
{

    LT lt = LT();
    EQ eq = EQ();

    std::shared_ptr<SortedRandomBagEntry<ELEMENT_TYPE>> x = header;
    std::vector<std::shared_ptr<SortedRandomBagEntry<ELEMENT_TYPE>> > update;
    update.resize(MAX_LEVEL+1);

    for (int i = level; i >= 0; i--)
    {
        while (x->forward[i] != nullptr && lt(x->forward[i]->value, value))
        {
            x = x->forward[i];
        }

        update[i] = x;
    }

    x = x->forward[0];

    if (x == nullptr)
    {
        return false;
    }

    if (eq(x->value, value))
    {
        // the number of objects will decrease by one in all cases
        num_objects--;

        // if the counter is higher than one, just decrease it
        if (x->counter > 1)
        {
            x->counter--;
            return true;
        }

        // otherwise delete the entry
        for (size_t i = 0; i <= level; i++)
        {
            if (update[i]->forward[i] != x)
            {
                update[i]->link_length[i]--;
            }

            else
            {
                update[i]->forward[i] = x->forward[i];
                update[i]->link_length[i] += x->link_length[i]-1;
            }
        }

        num_entries--;

        while (level > 0 && header->forward[level] == nullptr)
        {
            level--;
        }

        return true;
    }

    else
    {
        return false;
    }
}


/*
// only for debugging
template <class ELEMENT_TYPE>
void sorted_set<ELEMENT_TYPE>::print(int lev) {
    const Entry<ELEMENT_TYPE> *x = header;
    int num = 0;
    double std = 0;
    while (x!=nullptr) {
    	std += x->link_length[lev]*x->link_length[lev];
    	num++;
        x = x->forward[lev];
    }
    std::cout << "l" << lev << " e: " << num << " std: " << std/num << "\n";
}

template <class ELEMENT_TYPE>
void sorted_set<ELEMENT_TYPE>::print() {
std::cout << "cap:" << capacity << " num:" << num_entries << " lev:" << level << " MaxLevel:" << MAX_LEVEL << "\n";
   for (int i=level; i>=0; i--)
	   print(i);
}
*/

}
}

#endif
