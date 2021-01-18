s2netR <- function(data, params, loss = "default", frame = "ExtJT", proj = "auto", fista = NULL, S3 = TRUE){
  
  switch (loss,
    logit = {type_loss = 1},
    linear = {type_loss = 0},
    {
      switch (data$type,
        regression = {type_loss = 0},
        classification = {type_loss = 1}
      )
    }
  )
  
  switch (frame,
    JT = {type_frame = 0},
    {type_frame = 1}
  )
  
  switch (proj,
    no = {type_proj = 0},
    yes = {type_proj = 1},
    {type_proj = 2}
  )
  
  if(!class(data)=="s2Data")stop("[data] must be a s2Data object")
  if(!class(params)=="s2Params")stop("[params] must be a s2Params object")
  
  obj = new(s2net, data, type_loss)
  
  if(class(fista)=="s2Fista")obj$setupFista(fista)
  
  obj$fit(params, type_frame, type_proj)
  
  if(S3){
    ret = list(
      s2Data = data,
      s2Params = params,
      s2Fista = fista,
      loss = loss,
      type_loss = type_loss,
      frame = frame,
      type_frame = type_frame,
      proj = proj,
      type_proj = type_proj,
      beta = obj$beta,
      intercept = obj$intercept
    )
    class(ret) = "s2netR"
  }else{
    ret = obj
  }
  
  return(ret)
}

# predict_response_s2netR = function(object, newX){
#   return(newX %*% object$beta + object$intercept)
# }
# predict_probs_s2netR = function(object, newX){
#   eta = predict_response_s2netR(object, newX)
#   return(1/(1 + exp(-eta)))
# }
# predict_class_s2netR = function(object, newX){
#   return()
# }
predict.s2netR = function(object, newX, type = "default", ...){
  switch (type,
    reponse = {
      type_pred = 1
    },
    probs = {
      type_pred = 2
    },
    class = {
      type_pred = 3
    },
    {type_pred = 0}
  )
  
  obj = new(s2net, object$s2Data, object$type_loss)
  obj$beta = object$beta
  obj$intercept = object$intercept
  obj$predict(newX, type_pred)
}