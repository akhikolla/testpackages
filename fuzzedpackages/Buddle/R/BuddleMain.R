
#' Implementing Statistical Classification and Regression.
#'
#' Build a multi-layer feed-forward neural network model for statistical classification and regression analysis with random effects.
#'@param formula.string a formula string or a vector of numeric values. When it is a string, it denotes a classification or regression equation, of the form label ~ predictors or response ~ predictors, where predictors are separated by + operator. If it is a numeric vector, it will be a label or a response variable of a classification or regression equation, respectively.
#'@param data a data frame or a design matrix. When formula.string is a string, data should be a data frame which includes the label (or the response) and the predictors expressed in the formula string. When formula.string is a vector, i.e. a vector of labels or responses, data should be an nxp numeric matrix whose columns are predictors for further classification or regression.
#'@param train.ratio a ratio that is used to split data into training and test sets. When data is an n-by-p matrix, the resulting train data will be a (train.ratio x n)-by-p matrix. The default is 0.7.
#'@param arrange a logical value to arrange data for the classification only (automatically set up to FALSE for regression) when splitting data into training and test sets. If it is true, data will be arranged for the resulting training set to contain the specified ratio (train.ratio) of labels of the whole data. See also Split2TrainTest().
#'@param batch.size a batch size used for training during iterations.
#'@param total.iter a number of iterations used for training.
#'@param hiddenlayer a vector of numbers of nodes in hidden layers.
#'@param batch.norm a logical value to specify whether or not to use the batch normalization option for training. The default is TRUE.
#'@param drop a logical value to specify whether or not to use the dropout option for training. The default is TRUE.
#'@param drop.ratio a ratio for the dropout; used only if drop is TRUE. The default is 0.1.
#'@param lr a learning rate. The default is 0.1.
#'@param init.weight a weight used to initialize the weight matrix of each layer. The default is 0.1.
#'@param activation a vector of activation functions used in all hidden layers. For two hidden layers (e.g., hiddenlayer=c(100, 50)), it is a vector of two activation functions, e.g., c("Sigmoid", "SoftPlus"). The list of available activation functions includes Sigmoid, Relu, LeakyRelu, TanH, ArcTan, ArcSinH, ElliotSig, SoftPlus, BentIdentity, Sinusoid, Gaussian, Sinc, and Identity. For details of the activation functions, please refer to Wikipedia.
#'@param optim an optimization method which is used for training. The following methods are available: "SGD", "Momentum", "AdaGrad", "Adam", "Nesterov", and "RMSprop."
#'@param type a statistical model for the analysis: "Classification" or "Regression."
#'@param rand.eff a logical value to specify whether or not to add a random effect into classification or regression.  
#'@param distr a distribution of a random effect; used only if rand.eff is TRUE. The following distributions are available: "Normal", "Exponential", "Logistic", and "Cauchy."
#'@param disp a logical value which specifies whether or not to display intermediate training results (loss and accuracy) during the iterations.
#'
#'@return A list of the following values: 
#'\describe{
#'\item{lW}{a list of n terms of weight matrices where n is equal to the number of hidden layers plus one.} 
#'
#'\item{lb}{a list of n terms of bias vectors where n is equal to the number of hidden layers plus one.}
#'
#'\item{lParam}{a list of parameters used for the training process.}
#'
#'\item{train.loss}{a vector of loss values of the training set obtained during iterations where its length is eqaul to number of epochs.}
#'
#'\item{train.accuracy}{a vector of accuracy values of the training set obtained during during iterations where its length is eqaul to number of epochs.}
#'
#'\item{test.loss}{a vector of loss values of the test set obtained during the iterations where its length is eqaul to number of epochs.}
#'
#'\item{test.accuracy}{a vector of accuracy values of the test set obtained during the iterations where its length is eqaul to number of epochs.}
#'
#'\item{predicted.softmax}{an r-by-n numeric matrix where r is the number of labels (classification) or 1 (regression), and n is the size of the test set. Its entries are predicted softmax values (classification) or predicted values (regression) of the test sets, obtained by using the weight matrices (lW) and biases (lb).}
#'
#'\item{predicted.encoding}{an r-by-n numeric matrix which is a result of one-hot encoding of the predicted.softmax; valid for classification only.}
#'
#'\item{confusion.matrix}{an r-by-r confusion matrix; valid classification only.}
#'
#'\item{precision}{an (r+1)-by-3 matrix which reports precision, recall, and F1 of each label; valid classification only.}
#'
#'}
#'
#'@examples
#'####################
#'# train.ratio = 0.6                    ## 60% of data is used for training
#'# batch.size = 10     
#'# total.iter = 100
#'# hiddenlayer=c(20,10)                ## Use two hidden layers
#'# arrange=TRUE                         #### Use "arrange" option 
#'# activations = c("Relu","SoftPlus")   ### Use Relu and SoftPlus 
#'# optim = "Nesterov"                   ### Use the "Nesterov" method for the optimization.
#'# type = Classification  
#'# rand.eff = TRUE                      #### Add some random effect
#'# distr="Normal"                       #### The random effect is a normal random variable
#'# disp = TRUE                          #### Display intemeidate results during iterations.
#'
#'
#'data(iris)
#'
#'lst = TrainBuddle("Species~Sepal.Width+Petal.Width", iris, train.ratio=0.6, 
#'          arrange=TRUE, batch.size=10, total.iter = 100, hiddenlayer=c(20, 10), 
#'          batch.norm=TRUE, drop=TRUE, drop.ratio=0.1, lr=0.1, init.weight=0.1, 
#'          activation=c("Relu","SoftPlus"), optim="Nesterov", 
#'          type = "Classification", rand.eff=TRUE, distr = "Normal", disp=TRUE)
#'
#'lW = lst$lW
#'lb = lst$lb
#'lParam = lst$lParam
#'
#'confusion.matrix = lst$confusion.matrix
#'precision = lst$precision
#'
#'confusion.matrix
#'precision
#'
#'
#'### Another classification example 
#'### Using mnist data
#'
#'
#'data(mnist_data)
#'
#'Img_Mat = mnist_data$Images
#'Img_Label = mnist_data$Labels
#'
#'                              ##### Use 100 images
#'
#'X = Img_Mat                   ### X: 100 x 784 matrix
#'Y = Img_Label                 ### Y: 100 x 1 vector
#'
#'lst = TrainBuddle(Y, X, train.ratio=0.6, arrange=TRUE, batch.size=10, total.iter = 100, 
#'                  hiddenlayer=c(20, 10), batch.norm=TRUE, drop=TRUE, 
#'                  drop.ratio=0.1, lr=0.1, init.weight=0.1, 
#'                  activation=c("Relu","SoftPlus"), optim="AdaGrad", 
#'                  type = "Classification", rand.eff=TRUE, distr = "Logistic", disp=TRUE)
#'
#'
#'confusion.matrix = lst$confusion.matrix
#'precision = lst$precision
#'
#'confusion.matrix
#'precision
#'
#'
#'
#'
#'
#'
#'###############   Regression example  
#'
#'
#'n=100
#'p=10

#'X = matrix(rnorm(n*p, 1, 1), n, p)  ## X is a 100-by-10 design matrix 
#'b = matrix( rnorm(p, 1, 1), p,1)
#'e = matrix(rnorm(n, 0, 1), n,1)
#'Y = X %*% b + e                     ### Y=X b + e

#'######### train.ratio=0.7
#'######### batch.size=20
#'######### arrange=FALSE
#'######### total.iter = 100
#'######### hiddenlayer=c(20)
#'######### activation = c("Identity")
#'######### "optim" = "Adam"
#'######### type = "Regression"
#'######### rand.eff=FALSE
#'
#'lst = TrainBuddle(Y, X, train.ratio=0.7, arrange=FALSE, batch.size=20, total.iter = 100, 
#'                  hiddenlayer=c(20), batch.norm=TRUE, drop=TRUE, drop.ratio=0.1, lr=0.1, 
#'                  init.weight=0.1, activation=c("Identity"), optim="AdaGrad", 
#'                  type = "Regression", rand.eff=FALSE, disp=TRUE)
#'
#'
#'
#'
#'@references
#'[1] Geron, A. Hand-On Machine Learning with Scikit-Learn and TensorFlow. Sebastopol: O'Reilly, 2017. Print.
#'@references
#'[2] Han, J., Pei, J, Kamber, M. Data Mining: Concepts and Techniques. New York: Elsevier, 2011. Print.
#'@references
#'[3] Weilman, S. Deep Learning from Scratch. O'Reilly Media, 2019. Print.
#'@export
#'@seealso CheckNonNumeric(), GetPrecision(), FetchBuddle(), MakeConfusionMatrix(), OneHot2Label(), Split2TrainTest()
#'@importFrom Rcpp evalCpp 
#'@useDynLib Buddle


TrainBuddle = function(formula.string, data, train.ratio=0.7, arrange=0, batch.size=10, 
                       total.iter=10000, hiddenlayer=c(100), batch.norm=TRUE, drop=TRUE, drop.ratio=0.1,
                      lr=0.1, init.weight=0.1, activation=c("Sigmoid"), optim="SGD", 
                      type="Classification", rand.eff=FALSE, distr="Normal",  disp=TRUE){
  
  
  ########## Changing R env to C++ env
  Train_ratio = train.ratio
  bArrange = arrange
  nBatch_Size = batch.size
  nTotal_Iterations = total.iter
  HiddenLayer = hiddenlayer
  bBatch = batch.norm
  bDrop = drop  
  drop_ratio = drop.ratio
  Activation = activation
  strOpt = optim
  Type = type
  bRand = rand.eff
  strDist = distr
  bDisp = disp
  
  d_learning_rate = lr
  d_init_weight = init.weight
  
  
  if(Type=="Regression"){
    bArrange=0
  }
  
  
  ######################
  nHiddenLayer =length(HiddenLayer)
  nAct = length(Activation)
  
  if(nAct>nHiddenLayer){
    stop("Length of Activation vector should be equal or smaller than the length of HiddenLayer.")
  }else if(nAct==nHiddenLayer){
    nstrVec=GetStrVec(Activation)
  }else{
    NewActVec=rep("", times=nHiddenLayer)
    NewActVec[1:nAct] = Activation
    NewActVec[(nAct+1):nHiddenLayer] = "Relu"
    nstrVec = GetStrVec(NewActVec)
  }
  
  ######################
  
  
  if(length(formula.string)==1){
    
    lOneHot = OneHotEncoding(formula.string, data)
    
    Y = lOneHot$Y
    X = lOneHot$X                                ### X:n x p
    T_Mat=lOneHot$OneHot                         #### T: rxn
    Label = lOneHot$Label
    
    dimm = dim(X)
    n = dimm[1]
    p = dimm[2]
    
    
  }else{
    Y = formula.string
    X = data                                   #### X : nxp
    dimm = dim(X)
    n = dimm[1]
    p = dimm[2]
    T_Mat = OneHotEncodingSimple(Y, n)         #### T: rxn
  }
  
  
  cn = count(Y)
  Label = cn$x
  
  lCheck = CheckNonNumeric(X)
  
  if(lCheck[[1]]!=0){
    print("There are non-numeric values in the design matrix X.")

    return(lCheck)

  }
  
  
  ###################### Split X and T into train and test

  nTrain = floor(n*Train_ratio)
  
  if(nTrain<=10){
    
    print(paste("The size of the train set is "+ nTrain, ". Increase the train ratio or get more data.", sep="") )
    
  }
  
  
  
  if(bArrange==1){
    
    lYX = Split2TrainTest(Y, X, Train_ratio)
    Y_test = lYX$y.test
    Y_train = lYX$y.train
    nTrain = length(Y_train)
    
    
    Y[1:nTrain] = Y_train
    Y[(nTrain+1):n] = Y_test
    T_Mat = OneHotEncodingSimple(Y, n)
    
    T_train = T_Mat[ , 1:nTrain]
    T_test = T_Mat[ , (nTrain+1):n]
    
    X_test = lYX$x.test
    X_train = lYX$x.train
    
    
  }else{
    
    Y_train = Y[1:nTrain]
    Y_test = Y[(nTrain+1):n]
    
    X_train = X[1:nTrain, ]
    X_test = X[(nTrain+1):n, ]
    
    T_train = T_Mat[ , 1:nTrain]
    T_test = T_Mat[ , (nTrain+1):n]
    
    
  }
  
  dimmT = dim(T_train)
  r = dimmT[1]
  
  nPerEpoch = nTrain/nBatch_Size
  nEpoch = floor(nTotal_Iterations/nPerEpoch)
  
  if(nBatch_Size>=nTrain){
    
    print("The batch size is bigger than the size of the train set.")
    print("The half of the size of the train set will be tried as a new batch size.")
    
    nBatch_Size = floor(nTrain/2)
    
    if(nBatch_Size==0){
      stop("Batch size is 0.")
    }
    
  }else{
    if(nEpoch==0){
      
      print("The number of epoch is zero. Increase total iteration number, reduce the train ratio, or increase the batch size.")
      stop()
    }
    
  }
  
  
  ############################### Start Buddle
  
  lst = Buddle_Main(t(X_train), T_train, t(X_test), T_test, nBatch_Size, 
              nTotal_Iterations, HiddenLayer, bBatch, bDrop, drop_ratio,
              d_learning_rate, d_init_weight,nstrVec, strOpt, Type, bRand, strDist, bDisp)
  
  lW=lst[[1]]; lb=lst[[2]]
  train_loss= lst[[3]][[1]]
  train_accuracy = lst[[3]][[2]]
  test_loss = lst[[3]][[3]]
  test_accuracy = lst[[3]][[4]]
  
  
  nLen = length(test_accuracy)
  plot(1:nLen, test_accuracy, main = "Accuracy: Training vs. Test", 
       ylab="Accuracy", xlab="Epoch", type="l", col="red", ylim=c(0,1))
  lines(1:nLen, train_accuracy, type="l", col="blue")
  legend("topleft", c("Test", "Train"), fill=c("red", "blue"))
  
  
  lParam = list(label = r, hiddenLayer=HiddenLayer, batch=bBatch, drop=bDrop, drop.ratio=drop_ratio,
                lr = d_learning_rate, init.weight = d_init_weight,  activation = nstrVec, 
                optim=strOpt, type = Type, rand.eff = bRand, distr = strDist, disp = bDisp)

  
  lst2 = Buddle_Predict(t(X_test), lW, lb, lParam)
  Predicted_SoftMax = lst2[[1]]
  Predicted_OneHotEconding = lst2[[2]]
  
  if(type == "Classification"){
    Predicted_Label = OneHot2Label(Predicted_OneHotEconding, Label)
    CM = MakeConfusionMatrix(Predicted_Label, Y_test, Label)
    Precision = GetPrecision(CM)
    
  }else{
    Predicted_Label = Predicted_SoftMax
    CM = NA
    Precision = NA
  }

  lResult = list(lW=lW, lb=lb, 
                 lParam = lParam,
                 train.loss=train_loss,
                 train.accuracy = train_accuracy,
                 test.loss=test_loss,
                 test.accuracy = test_accuracy,
                 predicted.softmax = t(Predicted_SoftMax),
                 predicted.encoding = t(Predicted_OneHotEconding),
                 predicted.label = Predicted_Label,
                 confusion.matrix = CM, precision=Precision)
  
  
  
  
  return(lResult)
}




#' Predicting Classification and Regression.
#'
#' Yield prediction (softmax value or value) for regression and classification for given data based on the results of training.
#'@param X a matrix of real values which will be used for predicting classification or regression.
#'@param lW a list of weight matrices obtained after training.
#'@param lb a list of bias vectors obtained after training.
#'@param lParam a list of parameters used for training. It includes: label, hiddenlayer, batch, drop, drop.ratio, lr, init.weight, activation, optim, type, rand.eff, distr, and disp.
#'
#'@return A list of the following values: 
#'\describe{
#'\item{predicted}{predicted real values (regression) or softmax values (classification).}
#'
#'\item{One.Hot.Encoding}{one-hot encoding values of the predicted softmax values for classification. For regression, a zero matrix will be returned. To convert the one-hot encoding values to labels, use OneHot2Label().}
#'}
#'
#'@examples
#'
#'### Using mnist data again
#'
#'data(mnist_data)
#'
#'X1 = mnist_data$Images       ### X1: 100 x 784 matrix
#'Y1 = mnist_data$Labels       ### Y1: 100 x 1 vector
#'
#'
#'
#'############################# Train Buddle
#'
#'lst = TrainBuddle(Y1, X1, train.ratio=0.6, arrange=TRUE, batch.size=10, total.iter = 100, 
#'                  hiddenlayer=c(20, 10), batch.norm=TRUE, drop=TRUE, 
#'                  drop.ratio=0.1, lr=0.1, init.weight=0.1, 
#'                  activation=c("Relu","SoftPlus"), optim="AdaGrad", 
#'                  type = "Classification", rand.eff=TRUE, distr = "Logistic", disp=TRUE)

#'
#'lW = lst[[1]]
#'lb = lst[[2]]
#'lParam = lst[[3]]
#'
#'
#'X2 = matrix(rnorm(20*784,0,1), 20,784)  ## Genderate a 20-by-784 matrix
#'
#'lst = FetchBuddle(X2, lW, lb, lParam)   ## Pass X2 to FetchBuddle for prediction
#'
#'
#'
#'
#'
#'@references
#'[1] Geron, A. Hand-On Machine Learning with Scikit-Learn and TensorFlow. Sebastopol: O'Reilly, 2017. Print.
#'@references
#'[2] Han, J., Pei, J, Kamber, M. Data Mining: Concepts and Techniques. New York: Elsevier, 2011. Print.
#'@references
#'[3] Weilman, S. Deep Learning from Scratch. O'Reilly Media, 2019. Print.
#'@export
#'@seealso CheckNonNumeric(), GetPrecision(), MakeConfusionMatrix(), OneHot2Label(), Split2TrainTest(), TrainBuddle()
#'@importFrom Rcpp evalCpp 
#'@useDynLib Buddle



FetchBuddle = function(X, lW, lb, lParam){
  
  if(is.matrix(X)==FALSE){
    p = length(X)
    tmpX = matrix(0, 1, p)
    for(i in 1:p){
      tmpX[1, i] = X[i]
    }
    rm(X)
    X = tmpX
  }else{
    dimm = dim(X)
    n = dimm[1]
    p = dimm[2]
  }
  
  lst = Buddle_Predict(t(X), lW, lb, lParam)
  lResult = list(predicted = lst[[1]], One.Hot.Encoding = lst[[2]])
  
  return(lResult)
  
}




#' Detecting Non-numeric Values.
#'
#' Check whether or not an input matrix includes any non-numeric values (NA, NULL, "", character, etc) before being used for training. If any non-numeric values exist, then TrainBuddle() or FetchBuddle() will return non-numeric results.
#'@param X an n-by-p matrix.
#'
#'@return A list of (n+1) values where n is the number of non-numeric values. The first element of the list is n, and all other elements are entries of X where non-numeric values occur. For example, when the (1,1)th and the (2,3)th entries of a 5-by-5 matrix X are non-numeric, then the list returned by CheckNonNumeric() will contain 2, (1,1), and (2,3).
#'
#'@examples
#'
#'n = 5;
#'p = 5;
#'X = matrix(0, n, p)       #### Generate a 5-by-5 matrix which includes two NA's. 
#'X[1,1] = NA
#'X[2,3] = NA
#'
#'lst = CheckNonNumeric(X)
#'
#'lst
#'
#'@export
#'@seealso GetPrecision(), FetchBuddle(), MakeConfusionMatrix(), OneHot2Label(), Split2TrainTest(), TrainBuddle()
#'@importFrom Rcpp evalCpp 
#'@useDynLib Buddle


CheckNonNumeric = function(X){
  
  dimm = dim(X)
  n = dimm[1]
  p = dimm[2]
  
  
  nInc = 0
  lst = list()
  nIndex=2
  
  for(i in 1:n){
    
    for(j in 1:p){
      
      val = X[i, j]
      if((is.na(val)==TRUE) || is.null(val)==TRUE || is.numeric(val)==FALSE){
        nInc = nInc+1
        lst[[nIndex]] = c(i,j)
        nIndex=nIndex+1
      }
    }
  }
  
  lst[[1]] = nInc
  
  return(lst)
  
}



#' Splitting Data into Training and Test Sets.
#'
#' Convert data into training and test sets so that the training set contains approximately the specified ratio of all labels. 
#'@param Y an n-by-1 vector of responses or labels.
#'@param X an n-by-p design matrix of predictors.
#'@param train.ratio a ratio of the size of the resulting training set to the size of data.
#'
#'
#'@return A list of the following values:
#'\describe{
#'
#'\item{y.train}{the training set of Y.}
#'\item{y.test}{the test set of Y.}
#'\item{x.train}{the training set of X.}
#'\item{x.test}{the test set of X.}
#'
#'}
#'@examples
#'
#'data(iris)
#'
#'Label = c("setosa", "versicolor", "virginica")
#'
#'
#'train.ratio=0.8
#'Y = iris$Species 
#'X = cbind( iris$Sepal.Length, iris$Sepal.Width, iris$Petal.Length, iris$Petal.Width)
#'
#'lst = Split2TrainTest(Y, X, train.ratio)
#'
#'Ytrain = lst$y.train
#'Ytest = lst$y.test
#'
#'length(Ytrain)
#'length(Ytest)
#'
#'length(which(Ytrain==Label[1]))
#'length(which(Ytrain==Label[2]))
#'length(which(Ytrain==Label[3]))
#'
#'length(which(Ytest==Label[1]))
#'length(which(Ytest==Label[2]))
#'length(which(Ytest==Label[3]))
#'
#'
#'
#'
#'@export
#'@seealso CheckNonNumeric(), GetPrecision(), FetchBuddle(), MakeConfusionMatrix(), OneHot2Label(), TrainBuddle()
#'@importFrom Rcpp evalCpp
#'@useDynLib Buddle


Split2TrainTest=function(Y, X, train.ratio){
  
  Train_ratio = train.ratio
  
  dimm = dim(X)
  n=dimm[1];p=dimm[2]
  
  cn = count(Y)
  cnx = cn$x
  newcnf = floor(cn$freq * Train_ratio )
  nLev = length(newcnf)
  for(i in 1:nLev){
    if(newcnf[i]==0){newcnf[i]=1}
    
  }
  
  nTrain = sum(newcnf)
  nTest = n-nTrain
  
  
  YTrain = rep(Y[1], times=nTrain)
  XTrain = X
  
  YTest = rep(Y[1], times=nTest)
  XTest = X
  
  nIncYTr=1;nIncYTst=1;nIncXTr=1;nIncXTst=1;
  
  for(i in 1:nLev){
    
    val = cnx[i]
    nMany = newcnf[i]       #### How many for train
    Wh = which(Y==val)     
    nLenWh = length(Wh) 
    Ord = Wh[1:nMany]         ############ Train index
    nLenOrd = length(Ord)
    
    for(i in 1:nLenOrd){
      nIndex = Ord[i]
      YTrain[nIncYTr] = Y[nIndex]
      XTrain[nIncYTr, ]= X[nIndex, ]
      nIncYTr = nIncYTr+1
      
    }
    
    if(nLenWh != nMany){
      
      NOrd = Wh[(nMany+1):nLenWh]      ############ Test index
      nLenNOrd = length(NOrd)
      
      for(i in 1:nLenNOrd){
        nIndex = NOrd[i]
        YTest[nIncYTst] = Y[nIndex]
        XTest[nIncYTst, ]= X[nIndex, ]
        nIncYTst = nIncYTst+1
        
      }
      
      
    }
    
    
  }
  
  nIncYTr = nIncYTr-1
  nIncYTst = nIncYTst-1
  
  XTrain = XTrain[1:nIncYTr,]
  XTest = XTest[1:nIncYTst,]
  
  lst = list(y.train=YTrain, y.test=YTest, x.train=XTrain, x.test=XTest)
  
  return(lst)
  
  
}



#' Obtaining Labels
#'
#' Convert a one-hot encoding matrix to a vector of labels.
#'@param OHE an r-by-n one-hot encoding matrix.
#'@param Label an r-by-1 vector of values or levels which a label can take.
#'
#'
#'@return An n-by-1 vector of labels.  
#'
#'@examples
#'
#'Label = c("setosa", "versicolor", "virginica")
#'r = length(Label)
#'
#'n=10
#'OHE = matrix(0, r, n)           ### Generate a random one-hot encoding matrix     
#'for(i in 1:n){
  
#'  if(i%%r==0){
#'    OHE[i, 3] = 1
#'  }else if(i\%\%r==1){
#'    OHE[i, 1] = 1
#'  }else{
#'    OHE[i, 2] = 1
#'  }
#'  
#'}
#'
#'pred.label = OneHot2Label(OHE, Label)
#'
#'pred.label
#'
#'
#'
#'
#'@export
#'@seealso CheckNonNumeric(), GetPrecision(), FetchBuddle(), MakeConfusionMatrix(), Split2TrainTest(), TrainBuddle()
#'@importFrom Rcpp evalCpp 
#'@useDynLib Buddle



OneHot2Label = function(OHE, Label){
  
  T_Mat = OHE
  dimm=dim(T_Mat)
  p = dimm[1]
  n = dimm[2]
  
  AnswerKey = as.character(Label)
  ans = rep("", times=n)
  
  for(i in 1:n){
    nIndex = which(T_Mat[,i]==1)
    ans[i] = AnswerKey[nIndex]
  }
  
  return(ans)
}





#' Making a Confusion Matrix.
#'
#' Create a confusion matrix from two vectors of labels: predicted label obtained from FetchBuddle() as a result of prediction and true label of a test set.
#'@param predicted.label a vector of predicted labels.
#'@param true.label a vector of true labels.
#'@param Label a vector of all possible values or levels which a label can take.
#'
#'@return An r-by-r confusion matrix where r is the length of Label.
#'
#'@examples
#'
#'
#'data(iris)
#'
#'Label = c("setosa", "versicolor", "virginica")
#'
#'predicted.label = c("setosa", "setosa",    "virginica", "setosa", "versicolor", "versicolor")
#'true.label      = c("setosa", "virginica", "versicolor","setosa", "versicolor", "virginica")
#'
#'confusion.matrix = MakeConfusionMatrix(predicted.label, true.label, Label)
#'precision = GetPrecision(confusion.matrix)
#'
#'confusion.matrix
#'precision
#'
#'
#'
#'
#'
#'@export
#'@seealso CheckNonNumeric(), GetPrecision(), FetchBuddle(), OneHot2Label(), Split2TrainTest(), TrainBuddle()
#'@importFrom Rcpp evalCpp 
#'@useDynLib Buddle


MakeConfusionMatrix = function(predicted.label, true.label, Label){
  
  predicted = as.character(predicted.label)
  answerkey = as.character(true.label)
  Label = as.character(Label)
  
  nLen = length(predicted)
  
  n = length(Label)
  
  CM = matrix(0,n,n)
  
  colnames(CM, do.NULL = TRUE)
  colnames(CM) = as.character(Label)
  
  rownames(CM, do.NULL = TRUE)
  rownames(CM) = as.character(Label)
  
  for(i in 1:nLen){
    
    val = predicted[i]
    ans = answerkey[i]
    
    nval = which(Label==val)
    nans = which(Label==ans)
    
    CM[nval, nans] = CM[nval, nans]+1
    
    
  }
  
  return(CM)
  
  
  
}




#' Obtaining Accuracy.
#'
#' Compute measures of accuracy such as precision, recall, and F1 from a given confusion matrix. 
#'@param confusion.matrix a confusion matrix.
#'
#'@return An (r+1)-by-3 matrix when the input is an r-by-r confusion matrix. 
#'
#'@examples
#'
#'data(iris)
#'
#'Label = c("setosa", "versicolor", "virginica")
#'
#'predicted.label = c("setosa", "setosa",    "virginica", "setosa", "versicolor", "versicolor")
#'true.label      = c("setosa", "virginica", "versicolor","setosa", "versicolor", "virginica")
#'
#'confusion.matrix = MakeConfusionMatrix(predicted.label, true.label, Label)
#'precision = GetPrecision(confusion.matrix)
#'
#'confusion.matrix
#'precision
#'
#'
#'@export
#'@seealso CheckNonNumeric(), FetchBuddle(), MakeConfusionMatrix(), OneHot2Label(), Split2TrainTest(), TrainBuddle()
#'@importFrom Rcpp evalCpp 
#'@useDynLib Buddle


GetPrecision = function(confusion.matrix){
  
  
  CM = confusion.matrix
  
  RT = rownames(CM)
  
  dimm = dim(CM)
  nClass = dimm[1]
  
  MeasureMatrix=matrix(0,(nClass+1), 3)
  
  
  
  PrecisionVec= rep(0,times=nClass)
  RecallVec= rep(0,times=nClass)
  F1Vec= rep(0,times=nClass)
  
  ConsVec = rep(0, times=nClass)
  
  nAccuracy = 0
  for(i in 1:nClass){
    TP = CM[i,i]
    FP = sum(CM[,i])- TP
    FN = sum(CM[i,])- TP
    
    precisionVal = TP/(TP+FP)
    recallVal = TP/(TP+FN)
    f1Val = (2*precisionVal*recallVal)/(precisionVal+recallVal)
    
    PrecisionVec[i] = precisionVal
    RecallVec[i] = recallVal
    F1Vec[i] = f1Val
    
    ConsVec[i] = sum(CM[i,])
    nAccuracy = nAccuracy+TP
  }
  
  cTitle=c("Precision", "Recall", "F1")
  
  colnames(MeasureMatrix, do.NULL = TRUE)
  colnames(MeasureMatrix)=cTitle
  rTitle = rep("", times=(nClass+1))
  rTitle[1:nClass+1] = RT
  rTitle[(nClass+1)] = "Total"
  rownames(MeasureMatrix, do.NULL = TRUE)
  rownames(MeasureMatrix) = rTitle
  
  MeasureMatrix[(1:nClass),1]= PrecisionVec
  MeasureMatrix[(1:nClass),2]= RecallVec
  MeasureMatrix[(1:nClass),3]= F1Vec
  
  nInstance = sum(ConsVec)
  
  MeasureMatrix[(nClass+1),1] = t(ConsVec)%*% PrecisionVec /  nInstance
  MeasureMatrix[(nClass+1),2] = t(ConsVec)%*% RecallVec /  nInstance
  MeasureMatrix[(nClass+1),3] = t(ConsVec)%*% F1Vec /  nInstance
  
  out = list()
  out[[1]] = MeasureMatrix
  out[[2]] = nAccuracy/nInstance
  return(out)
}





ListVar = function(Str){
  
  fstr = formula(Str)
  
  Response = fstr[[2]]
  lPredictor = fstr[[3]]
  
  lst = list()
  lst[[1]] = Response
  
  nIter = 100000
  
  nInc=2
  for(i in 1:nIter){
    
    len = length(lPredictor)
    
    if(len==1){
      
      lst[[nInc]] = lPredictor
      break
    }else{
      lst[[nInc]] = lPredictor[[3]]
      nInc=nInc+1
      lPredictor = lPredictor[[2]]
    }
    
    
    
  }
  
  
  return(lst)
  
}


SplitVariable = function(Str){
  
  if(class(Str)=="character"){
    fstr = formula(Str)
  }else{
    fstr = Str
  }
  
  Response = fstr[[2]]
  lPredictor = fstr[[3]]
  
  lst = list()
  lst[[1]] = Response
  
  nIter = 100000
  
  nInc=2
  for(i in 1:nIter){
    
    len = length(lPredictor)
    
    if(len==1){
      
      lst[[nInc]] = lPredictor
      break
    }else{
      lst[[nInc]] = lPredictor[[3]]
      nInc=nInc+1
      lPredictor = lPredictor[[2]]
    }
    
    
    
  }
  
  
  return(lst)
  
}





OneHotEncodingSimple = function(y, n){
  
  ############################ Make T mat
  
  cn = count(y)
  
  #lev = as.numeric( as.character( cn$x ) ) 
  lev = cn$x
  
  # if(is.na(lev[1])==TRUE){
  #   lev = as.character( cn$x )
  # }
  nLen = length(lev)
  
  T_Mat = matrix(0, nLen, n)
  
  for(i in 1:n){
    yVal = y[i]
    nIndex = which(lev==yVal)
    T_Mat[nIndex, i]=1
    
  }
  
  
  return(T_Mat)
  
}


OneHotEncoding = function(Str, DM){
  
  dimm = dim(DM)
  n = dimm[1]
  
  lst = ListVar(Str)
  
  nLen = length(lst)
  nPred = nLen-1
  
  ColVec = rep("", nPred)
  yname = as.character( lst[[1]])
  
  for(i in 1:nPred){
    ColVec[i] = as.character( lst[[nLen+1-i]] )
  }
  
  ############################ Make T mat
  
  y = DM[, yname]
  
  cn = count(y)
  lev = cn$x
  # lev = as.numeric( as.character( cn$x ) ) 
  # 
  # if(is.na(lev[1])==TRUE){
  #   lev = as.character( cn$x )
  # }
  nLen = length(lev)
  
  T_Mat = matrix(0, nLen, n)
  
  for(i in 1:n){
    yVal = y[i]
    nIndex = which(lev==yVal)
    T_Mat[nIndex, i]=1
    
  }
  
  FreqVec = as.numeric(cn$freq)
  ############################ Make X
  X = matrix(0, n, nPred)
  
  for(i in 1:nPred){
    Var = ColVec[i]
    X[, i] = DM[, Var]
  }
  
  lst = list(Y=y, X=X, OneHotMatrix=T_Mat, Label = lev)
  return(lst)
  
}











RActiveStr2Num = function(Str){
  
  if(Str=="Sigmoid"){
    return(1)
  }else if(Str=="Relu"){
    return(2)
  }else if(Str=="LeakyRelu"){
    return(3)
  }else if(Str=="TanH"){
    return(4)
  }else if(Str=="ArcTan"){
    return(5)
  }else if(Str=="ArcSinH"){
    return(6)
  }else if(Str=="ElliotSig"){
    return(7)
  }else if(Str=="SoftPlus"){
    return(8)
  }else if(Str=="BentIdentity"){
    return(9)
  }else if(Str=="Sinusoid"){
    return(10)
  }else if(Str=="Gaussian"){
    return(11)
  }else if(Str=="Sinc"){
    return(12)
  }else{
    return(2)
  }
  
}

GetStrVec = function(ActVec){
  
  n = length(ActVec)
  ans = rep(0, times=n)
  for(i in 1:n){
    ans[i] = RActiveStr2Num(ActVec[i])
  }
  
  return(ans)
}





















