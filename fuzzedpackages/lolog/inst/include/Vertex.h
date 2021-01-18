#ifndef VERTEXH_
#define VERTEXH_

#include <set>
//#include <ext/hashSet>
#include <vector>
#include <memory>
#include <boost/shared_ptr.hpp>
#include <boost/container/flat_set.hpp>
namespace lolog{


//typedef _GnuCxx::hashSet<int> Set;
//typedef std::set<int> Set;
typedef boost::container::flat_set<int> Set;
typedef boost::shared_ptr< Set > SetPtr;
typedef boost::shared_ptr< const Set > ConstSetPtr;


/*!
 * A class for a vertex (node)
 */
class Vertex {
protected:
    int idNum;
    std::vector<double> contVar;
    std::vector<int> disVar;
    std::vector<bool> contObs;
    std::vector<bool> disObs;

public:
    Vertex() : idNum(-1){}

    virtual ~Vertex(){}

    /*!
     * every vertex in a network has a unique id
     * \returns the id
     */
    inline int id(){
        return idNum;
    }

    /*!
     * set the vertex's id
     * \param newId the id
     */
    inline void setId(int newId){
        idNum = newId;
    }

    /*!
     * gets the value for a continuous nodal variable
     * \param index the index of the continuous variable
     * \returns the value
     */
    inline double continVariable(int index) const{
        return contVar[index];
    }

    /*!
     * sets the value for a continuous nodal variable
     * \param value the new value
     * \param index which variable to set
     *
     */
    inline void setContinVariable(double value,int index){
        contVar[index] = value;
    }
    /*!
     * gets the value for a discrete nodal variable
     * \param index the index of the discrete variable
     * \returns the value
     */
    inline int discreteVariable(int index) const{
        return disVar[index];
    }

    /*!
     * sets the value for a discrete nodal variable
     * \param value the new value
     * \param index which variable to set
     *
     */
    inline void setDiscreteVariable(int value,int index){
        disVar[index] = value;
    }

    /*!
     * removes a continuous variable
     * \param index which to remove
     */
    inline void removeContinVariable(int index){
        contVar.erase(contVar.begin() + index);
        contObs.erase(contObs.begin() + index);
    }

    /*!
     * removes a discrete variable
     * \param index which to remove
     */
    inline void removeDiscreteVariable(int index){
        disVar.erase(disVar.begin() + index);
        disObs.erase(disObs.begin() + index);
    }

    /*!
     * add a continuous variable to the end of the variables
     * \param value the value to add
     */
    inline void addContinVariable(double value){
        contVar.push_back(value);
        contObs.push_back(true);
    }

    /*!
     * add a discrete variable to the end of the variables
     * \param value the value to add
     */
    inline void addDiscreteVariable(int value){
        disVar.push_back(value);
        disObs.push_back(true);
    }

    /*!
     * is observed
     * \param index which variable
     * \returns true if it is observed. false if missing
     */
    inline bool continObserved(int index){
        return contObs[index];
    }

    /*!
     * is observed
     * \param index which variable
     * \returns true if it is observed. false if missing
     */
    inline bool discreteObserved(int index){
        return disObs[index];
    }

    /*!
     * sets missingness
     * \param index which variable
     * \param observed missingness mask
     */
    inline void setDiscreteObserved(int index,bool observed){
        disObs[index] = observed;
    }

    /*!
     * sets missingness
     * \param index which variable
     * \param observed missingness mask
     */
    inline void setContinObserved(int index,bool observed){
        contObs[index] = observed;
    }





};


}

#endif /* VERTEXH_ */
