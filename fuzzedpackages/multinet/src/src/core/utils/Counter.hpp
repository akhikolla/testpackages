/**
 * Classes to count the number of occurrences of some objects/values.
 *
 * History:
 * - 2018.01.01 file imported from version 1.0 of the multinet library
 */


#ifndef UU_CORE_UTILS_COUNTER_H_
#define UU_CORE_UTILS_COUNTER_H_

#include <unordered_map>

namespace uu {
namespace core {


/**********************************************************************/
/** Counter ***********************************************************/
/**********************************************************************/

/**
 * A Counter for a single value
 */
template <class T>
class Counter
{
  private:
    /** A map where the count is kept */
    std::unordered_map<T,int> values;

  public:
    /**
     * Constructor.
     */
    Counter(
    );

    /**
     * Destructor.
     */
    ~Counter(
    );

    /**
     * This function is used to increase the count of T.
     * @param val value whose count should be incremented
     */
    void
    inc(
        const T& val
    );


    /**
     * This function is used to set a given value for the count of T.
     * @param key value whose number of occurrences should be incremented
     * @param num number of occurrences of key
     */
    void
    set(
        const T& key,
        int num
    );

    /**
     * @param val value whose count should be returned
     * @return the count of T
     */
    int
    count(
        const T& val
    ) const;

    /**
     * @return the key with the highest number of occurrences (in case of tie, any of the maximal keys can be returned)
     */
    T
    max(
    ) const;

    /**
     * Grants access to the internal data structure storing the numbers of occurrences.
     * @return a reference to the map storing the number of occurrences
     */
    const std::unordered_map<T,int>&
    map(
    );
};

/**
 * A Counter for a pair of values
 */
template <class T1, class T2>
class PairCounter
{
  private:
    /** A map where the count is kept */
    std::unordered_map<T1, std::unordered_map<T2,int> > values;

  public:
    /**
     * Constructor.
     */
    PairCounter(
    );

    /**
     * Destructor.
     */
    ~PairCounter(
    );

    /**
     * This function is used to increase the count of the pair (key1,key2).
     * @param key1 first part of the value whose number of occurrences should be incremented
     * @param key2 second part of the value whose number of occurrences should be incremented
     */
    void
    inc(
        const T1& key1,
        const T2& key2
    );

    /**
     * This function is used to set a given value for the count of the pair (key1,key2).
     * @param key1 first part of the value whose number of occurrences should be set
     * @param key2 second part of the value whose number of occurrences should be set
     * @param num number of occurrences of the pair (key1,key2) to be set
     */
    void
    set(
        const T1& key1,
        const T2& key2,
        int num
    );

    /**
     * @param key1 first part of the value whose number of occurrences should be returned
     * @param key2 second part of the value whose number of occurrences should be returned
     * @return count of the pair (key1,key2)
     */
    int
    count(
        const T1& key1,
        const T2& key2
    ) const;

    /**
     * Grants access to the internal data structure storing the numbers of occurrences.
     * @return a reference to the map storing the number of occurrences
     */
    std::unordered_map<T1, std::unordered_map<T2,int> >&
    map(
    );
};

/**
 * A Counter for a triplet of values
 */
template <class T1, class T2, class T3>
class TripletCounter
{
  private:
    /** A map where the count is kept */
    std::unordered_map<T1, std::unordered_map<T2, std::unordered_map<T3,int> > > values;

  public:
    /**
     * Constructor.
     */
    TripletCounter(
    );

    /**
     * Destructor.
     */
    ~TripletCounter(
    );

    /**
     * This function is used to increase the count of the triple (key1,key2,key3).
     * @param key1 first part of the value whose number of occurrences should be incremented
     * @param key2 second part of the value whose number of occurrences should be incremented
     * @param key3 third part of the value whose number of occurrences should be incremented
     */
    void
    inc(
        const T1& key1,
        const T2& key2,
        const T3& key3
    );

    /**
     * This function is used to set a given value for the count of the triple (key1,key2,key3).
     * @param key1 first part of the value whose number of occurrences should be set
     * @param key2 second part of the value whose number of occurrences should be set
     * @param key3 third part of the value whose number of occurrences should be set
     * @param num the number of occurrences of the triple (key1,key2,key3) to be set
     */
    void
    set(
        const T1& key1,
        const T2& key2,
        const T3& key3,
        int num
    );

    /**
     * @param key1 first part of the value whose number of occurrences should be returned
     * @param key2 second part of the value whose number of occurrences should be returned
     * @param key3 third part of the value whose number of occurrences should be returned
     * @return the number of occurrences of the triple (key1,key2,key3)
     */
    int
    count(
        const T1& key1,
        const T2& key2,
        const T3& key3
    ) const;

    /**
     * Grants access to the internal data structure storing the numbers of occurrences.
     * @return a reference to the map storing the number of occurrences
     */
    std::unordered_map<T1, std::unordered_map<T2, std::unordered_map<T3,int> > >&
    map(
    );

};

template <class T>
Counter<T>::Counter() {}

template <class T>
Counter<T>::~Counter()
{
}

template <class T>
void
Counter<T>::inc(
    const T& val
)
{
    if (values.count(val)==0)
    {
        values[val] = 0;
    }

    values[val]++;
}

template <class T>
void
Counter<T>::set(
    const T& val,
    int num
)
{
    values[val] = num;
}

template <class T>
int
Counter<T>::count(
    const T& val
) const
{
    if (values.count(val)==0)
    {
        return 0;
    }

    else
    {
        return values.at(val);
    }
}

template <class T>
T
Counter<T>::max(
) const
{
    int max = -1;
    T max_value = 0;

    for (auto pair: values)
    {
        if (pair.second>max)
        {
            max_value = pair.first;
            max = pair.second;
        }
    }

    return max_value;
}

template <class T>
const std::unordered_map<T, int>&
Counter<T>::map(
)
{
    return values;
}

template <class T1, class T2>
PairCounter<T1, T2>::PairCounter(
)
{

}

template <class T1, class T2>
PairCounter<T1, T2>::~PairCounter(
)
{
}

template <class T1, class T2>
void
PairCounter<T1, T2>::inc(
    const T1& val1,
    const T2& val2
)
{
    if (values.count(val1)==0 || values.at(val1).count(val2)==0)
    {
        values[val1][val2] = 0;
    }

    values[val1][val2]++;
}

template <class T1, class T2>
void
PairCounter<T1, T2>::set(
    const T1& val1,
    const T2& val2,
    int num
)
{
    values[val1][val2] = num;
}

template <class T1, class T2>
int
PairCounter<T1,T2>::count(
    const T1& val1,
    const T2& val2
) const
{
    if (values.count(val1)==0 || values.at(val1).count(val2)==0)
    {
        return 0;
    }

    else
    {
        return values.at(val1).at(val2);
    }
}

template <class T1, class T2>
std::unordered_map<T1, std::unordered_map<T2,int> >&
PairCounter<T1,T2>::map(
)
{
    return values;
}

template <class T1, class T2, class T3>
TripletCounter<T1, T2, T3>::TripletCounter(
)
{

}

template <class T1, class T2, class T3>
TripletCounter<T1, T2, T3>::~TripletCounter(
)
{
}

template <class T1, class T2, class T3>
void
TripletCounter<T1, T2, T3>::inc(
    const T1& val1,
    const T2& val2,
    const T3& val3
)
{
    if (values.count(val1)==0 || values.at(val1).count(val2)==0 || values.at(val1).at(val2).count(val3)==0)
    {
        values[val1][val2][val3] = 0;
    }

    values[val1][val2][val3]++;
}

template <class T1, class T2, class T3>
void
TripletCounter<T1, T2, T3>::set(
    const T1& val1,
    const T2& val2,
    const T3& val3,
    int num
)
{
    values[val1][val2][val3] = num;
}

template <class T1, class T2, class T3>
int
TripletCounter<T1,T2,T3>::count(
    const T1& val1,
    const T2& val2,
    const T3& val3
) const
{
    if (values.count(val1)==0 || values.at(val1).count(val2)==0 || values.at(val1).at(val2).count(val3)==0)
    {
        return 0;
    }

    else
    {
        return values.at(val1).at(val2).at(val3);
    }
}

template <class T1, class T2, class T3>
std::unordered_map<T1, std::unordered_map<T2, std::unordered_map<T3,int> > >&
TripletCounter<T1,T2,T3>::map(
)
{
    return values;
}

} // namespace core
} // namespace uu

#endif /* UU_CORE_UTILS_COUNTER_H_ */
