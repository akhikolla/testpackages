/*--------------------------------------------------------------------*/
/*     Copyright (C) 2011-2015  <MODAL team @INRIA,Lille & U.M.R. C.N.R.S. 6599 Heudiasyc, UTC>

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as
 published by the Free Software Foundation; either version 2 of the
 License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public
 License along with this program; if not, write to the
 Free Software Foundation, Inc.,
 59 Temple Place,
 Suite 330,
 Boston, MA 02111-1307
 USA

 Contact : parmeet.bhatia@inria.fr , bhatia.parmeet@gmail.com
 */

/** @file InputParameters.h
 *  @brief This file declares InputParameters class.
 **/

#ifndef INPUTPARAMETERS_H_
#define INPUTPARAMETERS_H_

#include "../typedefs/typedef.h"
#include <list>
#include "../enumerations/enumerations.h"

/**
 * @brief Structure to store various Strategy parameters.
 */
struct StrategyParameters
{
    //iterations for various stages of XStrategy and between EM algo for Rows and Columns
    int nbinitmax_;
    int nbtry_;
    int nbxem_;
    int nbiter_xem_;
    int nbiter_XEM_;
    //member function pointers
    StopCriteria_poiter Stop_Criteria;
};

/**
 * @brief Structure to store various model parameters.
 */
struct ModelParameters
{
    //iterations and epsilon inside model
    STK::Real epsilon_int_ ;
    int nbiterations_int_;
    STK::Real eps_xem_;
    STK::Real eps_XEM_;

    //iterations and epsilon inside initialization
    int nbinititerations_;
    STK::Real initepsilon_;

    //epsilon set either to eps_xem_ or esp_XEM_ depending on where we are in XEMstrategy
    STK::Real epsilon_;

    //other input options, self explanatory
    int nbrowclust_;
    int nbcolclust_;
    int nbRow_;
    int nbCol_;

    //proportion
    bool fixedproportions_;
};

std::ostream& operator<< (std::ostream& s, ModelParameters const& param);

/**
 * @brief Structure to store various Strategy options.
 */
struct Strategy
{
    Algorithm Algo_;
    StopCriteria stopcriteria_;
    Initialization Init_;
    Model Model_;
    DataType DataType_;
    bool SemiSupervised;
};

/** @brief This class contains all the parameters that are required to run various algorithms.
 *  This class publicly inherits IPage class defined in STKpp
 *  project http://www.stkpp.org/html/classSTK_1_1IPage.html
 *  which allow easy access to reading and writing option values
 *  from text file provided by the user.
 */
#ifdef STK_DMANAGER
class InputParameters : public STK::IPage
{
  public:
    /**Constructor*/
    InputParameters(STK::Integer const& level);
    /**
     * function that must be called to initialize various enumeration mappings
     */
    void initializeParamEnum();
    /**
     * This function read various input parameters from text file
     * @param optionfilename path to text file with input parameters
     */
    void ReadFromOptionFile(std::string optionfilename);
    /**Get Model parameters*/
    const ModelParameters& GetModelparameters() const;
    /**Get Algorithm parameters*/
    const StrategyParameters& GetStrategyparameters() const;
    /**Get Strategy*/
    const Strategy& GetStrategy();
    /**Get filename*/
    std::string GetDatafilename(){return datafilename_;}
    /**Get optionalfilenames*/
    std::list<std::string>& GetOptionalfilenames(){return optionalfilenames_;}
    void setModelparameters(ModelParameters const&);
    /**Set Algorithm parameters*/
    void setAlgoparameters(StrategyParameters const&);
    /**Set strategy*/
    void setStrategy(Strategy const&);
    /**Destructor*/
    ~InputParameters(){};

  protected:
    ModelParameters Mparam_;
    StrategyParameters strategyParam_;
    Strategy strategy_;
    std::string datafilename_;
    std::list<std::string> optionalfilenames_;
    /**String to enumeration mappings*/
    std::map<std::string,Algorithm> S_Algorithm;
    std::map<std::string,StopCriteria> S_StopCriteria;
    std::map<std::string,DataType> S_DataType;
    std::map<std::string,Initialization> S_Init;
    std::map<std::string,Model> S_Model;

};
#else
class InputParameters
{
  public:
    /**Constructor*/
    InputParameters(){initializeParamEnum();}
    /**Get Model parameters*/
    const ModelParameters& GetModelparameters() const;
    /**
     * function that must be called to initialize various enumeration mappings
     */
    void initializeParamEnum();
    /**Get Algorithm parameters*/
    const StrategyParameters& GetStrategyparameters() const;
    /**Get Strategy*/
    const Strategy& GetStrategy();
    /**Get filename*/
    std::string GetDatafilename(){return datafilename_;}
    /**Get optionalfilenames*/
    std::list<std::string>& GetOptionalfilenames(){return optionalfilenames_;}
    /**Set Model parameters*/
    void setModelparameters(ModelParameters const&);
    /**Set Algorithm parameters*/
    void setAlgoparameters(StrategyParameters const&);
    /**Set strategy*/
    void setStrategy(Strategy const&);
    /**Destructor*/
    ~InputParameters(){};

  protected:
    ModelParameters Mparam_;
    StrategyParameters strategyParam_;
    Strategy strategy_;
    std::string datafilename_;
    std::list<std::string> optionalfilenames_;
    std::map<std::string,Algorithm>       S_Algorithm;
    std::map<std::string,StopCriteria>    S_StopCriteria;
    std::map<std::string,DataType>        S_DataType;
    std::map<std::string,Initialization>  S_Init;
    std::map<std::string,Model>           S_Model;

};
#endif


inline const ModelParameters& InputParameters::GetModelparameters() const
{ return Mparam_;}

inline const StrategyParameters& InputParameters::GetStrategyparameters() const
{ return strategyParam_;}

inline const Strategy& InputParameters::GetStrategy()
{ return strategy_;}

inline void InputParameters::setAlgoparameters(StrategyParameters const& Aparam)
{ strategyParam_ = Aparam;}

inline void InputParameters::setModelparameters(ModelParameters const& Mparam)
{ Mparam_ = Mparam;}

inline void InputParameters::setStrategy(Strategy const& strat)
{ strategy_ = strat;}

#endif /* INPUTPARAMETERS_H_ */
