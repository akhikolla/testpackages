#' @title hiddenBlock
#'
#' @description This function is to construct the hidden block .
#'
#' @param inlayer input layer, keras layer  
#'
#' @param  nodes list of integers, list of the number of nodes for all the hidden layers 
#'
#' @param  acts list of strings, list of activations for hidden layers 
#'
#' @param  idepth integer, index of the hidden layer 
#'
#' @param  orginlayer keras layer, original layer to be added to decoding layer (default: NULL)
#'
#' @param  reg string, regularization (default: NULL)
#'
#' @param  dropout double, dropout rate for the target hidden layer (default 0) 
#'
#' @param  batchnorm bool, flag to conduct batch normalization (default: TRUE) 
#'
#' @return keras layer, block of a hidden layer (with addtion of actionvation or/and batch normalization) 
#'
#' @export hiddenBlock


hiddenBlock=function(inlayer, nodes,acts,idepth, orginlayer=NULL,reg=NULL,dropout=0,batchnorm=TRUE){
  layer=inlayer %>% keras::layer_dense(units =nodes[idepth], kernel_initializer ="glorot_uniform", bias_initializer = "zeros",
                                       kernel_regularizer=reg)
  layer=layer %>% keras::layer_activation(acts[idepth])
  if(batchnorm){
    layer=layer %>% keras::layer_batch_normalization()
  }
  if(!is.null(dropout) && dropout!=0){
    layer=layer %>% keras::layer_dropout(dropout)
  }
  if(!is.null(orginlayer)){
    layer=keras::layer_add(inputs=c(orginlayer,layer))
    layer=layer %>% keras::layer_activation(acts[idepth])
    if(batchnorm){
      layer=layer %>% keras::layer_batch_normalization()
    }
  }
  return(layer)
}

#' @title AutoEncoderModel
#'
#' @description This function is to construct a residual autoencoder-based deep network.
#'
#' @param nfea integer, Number of features 
#'
#' @param  nout integer, Number of output units 
#'
#' @param  nodes list of integers, list of the number of nodes for the hidden layers in encoding component  
#'
#' @param  acts list of strings, list of activation function names 
#'
#' @param  mdropout double, dropout rate of the coding (middle) layer (default:0)
#'
#' @param  reg string, regularization string (default: NULL)
#'
#' @param  batchnorm bool, flag to conduct batch normalization (default:TRUE)
#'
#' @param  isres bool, flag to conduct the residual connections (default: TRUE)
#'
#' @param  outtype integer, output type, 0 indicating nout outputs and 1 nout+nfea outputs (default: 0)
#'
#' @param  fact string, activation for output layer (default:"linear")
#'
#' @return keras model, model of (residual) autoencoder-based deep network  
#'
#' @export AutoEncoderModel
#' 
AutoEncoderModel=function(nfea,nout,nodes,acts,mdropout=0,reg=NULL,batchnorm=TRUE,isres=TRUE,outtype=0,fact="linear"){
  S = rstack::stack$new()
  inlayer=keras::layer_input(shape=nfea,name='feat')
  layer=inlayer
  for(i in c(1:length(nodes))){
    layer=layer %>% keras::layer_dense(units =nodes[i], kernel_initializer ="glorot_uniform", bias_initializer = "zeros",
                                       kernel_regularizer=reg)
    layer=layer %>% keras::layer_activation(acts[i])
    if(batchnorm){
      layer=layer %>% keras::layer_batch_normalization()
    }
    if(i<length(nodes)){
      S$push(layer)
    }else{
      if(!is.null(mdropout) && mdropout!=0){
        layer=layer %>% keras::layer_dropout(mdropout)
      }
    }
  }
  for(i in c((length(nodes)-1):1)){
    originlayer=S$pop()
    layer = hiddenBlock(layer,nodes,acts,i,originlayer,reg,NULL,batchnorm)
  }
  layer=layer %>% keras::layer_dense(units =nfea, kernel_initializer ="glorot_uniform", bias_initializer = "zeros",
                                     kernel_regularizer=reg)
  if(isres){
    layer=keras::layer_add(inputs=c(inlayer,layer))
    layer=layer %>% keras::layer_activation(acts[1])
    if(batchnorm){
      layer=layer %>% keras::layer_batch_normalization() 
    }
  }
  pnout=nout
  if(outtype==1){
    pnout=pnout+nfea
  }
  outlayer =layer %>% keras::layer_dense(units =pnout, kernel_initializer ="glorot_uniform", bias_initializer = "zeros",
                                         activation=fact, kernel_regularizer=reg)
  model = keras::keras_model(inputs = inlayer, outputs = outlayer)
  return(model)
}
