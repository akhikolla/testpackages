

#ifndef STK_PARTITIONHANDLER_H
#define STK_PARTITIONHANDLER_H

#include <Arrays/include/STK_CArrayVector.h>
#include <STatistiK/include/STK_Law_UniformDiscrete.h>
#include <Sdk/include/STK_IRunner.h>

namespace STK
{

/** @ingroup DManager
 *  CvHanler is an utility function for building the submatrix/subvectors
 *  needed when creating learning and test data sets.
 **/
class PartitionHandler: public IRunnerBase
{
  public:
    /** Default constructor.
     *  @param rangeData prop range of the data to partition
     *  @param  prop to set in test data
     **/
    PartitionHandler( Range const& rangeData, Real prop);
    /** destructor */
    inline virtual ~PartitionHandler() {}

    /** @return the number of folds */
    inline Real const& proportion() const { return prop_;}
    /** @return the range of the data */
    inline Range const& rangeData() const { return rangeData_;}
    /** @return the partitions */
    inline CVectorXi const& partitions() const { return partitions_;}

    inline virtual bool run()
    { partition(); hasRun_ = true; return true;}

    inline void setData( Range const& rangeData, Real prop)
    {
      partitions_.clear();
      rangeData_ = rangeData;
      prop_     = prop;
      sizeTest_ = int(rangeData_.size() * prop_);
      hasRun_   = false;
    }
    /** get the data set when setting out fold k and test data set  */
    template<class Data>
    bool getPartitions( Data const& x, Data& xLearn, Data& xTest);
    /** get the data set when setting out fold k and test data set  */
    template<class xData, class yData>
    bool getPartitions( xData const& x, xData& xLearn, xData& xTest
                      , yData const& y, yData& yLearn, yData& yTest);

  protected:
    /** create a random partition */
    void partition();

  private:
    /** Range of the data set (number of rows) */
    Range rangeData_;
    /** proportion */
    Real prop_;
    /** size of the test (should  */
    int sizeTest_;
    /** repartition of the sample into k-folds */
    CVectorXi partitions_;
};

/* Default constructor. nbLearns is set to the number of observation
 *  @param rangeData the range of the data to set
 *  @param nbLearns numbbe of Learns
 **/
inline PartitionHandler::PartitionHandler( Range const& rangeData, Real prop)
                                         : IRunnerBase()
                                         , rangeData_(rangeData), prop_(prop)
                                         , sizeTest_(int(rangeData_.size()*prop_))
                                         , partitions_()
{
  // check proportion
  if (prop_>1)
  { STKRUNTIME_ERROR_1ARG(PartitionHandler::PartitionHandler,prop,prop>1);}
  if (prop_<0)
  { STKRUNTIME_ERROR_1ARG(PartitionHandler::PartitionHandler,prop,prop<0);}
}
/* get the data set when setting out fold k and test data set  */
template<class Data>
bool PartitionHandler::getPartitions( Data const& x,Data& xLearn, Data& xTest)
{
  // check if partitions are determined
  if (!hasRun_)
  { msg_error_ = STKERROR_NO_ARG(PartitionHandler::getPartitions,PartitionHandler has to run);
    return false;
  }
  // check dimensions
  if (x.rows() != rangeData_)
  { msg_error_ = STKERROR_NO_ARG(PartitionHandler::getKLearn,x.rows()!=rangeData_);
    return false;
  }
  // prepare containers
  Range xLearnRows = x.rows();
  xLearnRows.decLast(sizeTest_);
  xLearn.resize(xLearnRows, x.cols());
  xTest.resize(sizeTest_, x.cols());
  // copy data
  int iLearnRow = xLearn.beginRows(), iTestRow = xTest.beginRows();
  for (int i = partitions_.begin(); i < partitions_.end(); ++i)
  {
    if (partitions_[i] == 1)
    {
      xTest.row(iTestRow) = x.row(i);
      ++iTestRow;
    }
    else
    {
      xLearn.row(iLearnRow) = x.row(i);
      ++iLearnRow;
    }
  }
  return true;
}
/* get the data set when setting out fold k and test data set  */
template<class xData, class yData>
bool PartitionHandler::getPartitions( xData const& x, xData& xLearn, xData& xTest
                                    , yData const& y, yData& yLearn, yData& yTest)
{
  // check if partitions are determined
  if (!hasRun_)
  { msg_error_ = STKERROR_NO_ARG(PartitionHandler::getKLearn,PartitionHandler has to run);
    return false;
  }
  // check dimensions
  if (x.rows() != rangeData_)
  { msg_error_ = STKERROR_NO_ARG(PartitionHandler::getKLearn,x.rows()!=rangeData_);
    return false;
  }
  if (y.rows() != rangeData_)
  { msg_error_ = STKERROR_NO_ARG(PartitionHandler::getKLearn,y.rows()!=rangeData_);
    return false;
  }
  // prepare constainers
  Range xLearnRows = x.rows();
  xLearnRows.decLast(sizeTest_);
  xLearn.resize(xLearnRows, x.cols());
  xTest.resize(sizeTest_, x.cols());
  yLearn.resize(xLearnRows, y.cols());
  yTest.resize(sizeTest_, y.cols());
  // copy data
  int iLearnRow = xLearn.beginRows(), iTestRow = xTest.beginRows();
  for (int i = partitions_.begin(); i < partitions_.end(); ++i)
  {
    if (partitions_[i] == 1)
    {
      xTest.row(iTestRow) = x.row(i);
      yTest.row(iTestRow) = y.row(i);
      ++iTestRow;
    }
    else
    {
      xLearn.row(iLearnRow) = x.row(i);
      yLearn.row(iLearnRow) = y.row(i);
      ++iLearnRow;
    }
  }
  return true;
}

/* create a random partition in k folds*/
inline void PartitionHandler::partition()
{
  partitions_.resize(rangeData_);
  int endTest = partitions_.begin()+sizeTest_;
  //fill the container with the index of the partition (1 test, 0 learn)
  for(int i = partitions_.begin() ; i< endTest ;i++) { partitions_[i] = 1;}
  for(int i = endTest ; i< partitions_.end() ;i++)   { partitions_[i] = 0;}
  //make a random rearrangement
  int begin = partitions_.begin();
  for (int i=partitions_.end()-2; i>begin; --i)
  { std::swap(partitions_[i], partitions_[Law::UniformDiscrete::rand(begin, i+1)]);}
}

} // namespace STK

#endif /* STK_PARTITIONHANDLER_H */
