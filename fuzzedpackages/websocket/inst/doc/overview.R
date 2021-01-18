## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- eval = FALSE------------------------------------------------------------
#  ws <- WebSocket$new("ws://echo.websocket.org/")

## ---- eval = FALSE------------------------------------------------------------
#  ws <- WebSocket$new("ws://echo.websocket.org/", autoConnect = FALSE)
#  # Set up callbacks here...
#  ws$connect()

## ---- eval = FALSE------------------------------------------------------------
#  {
#    ws <- WebSocket$new("ws://echo.websocket.org/")
#    ws$onOpen(function(event) { message("websocket opened") })
#  }

## ---- eval = FALSE------------------------------------------------------------
#  {
#    ws <- WebSocket$new("ws://echo.websocket.org/")
#    ws$onOpen(function(event) {
#      cat("connected\n")
#    })
#  }

## ---- eval = FALSE------------------------------------------------------------
#  {
#    ws <- WebSocket$new("ws://echo.websocket.org/")
#    removeThis <- ws$onMessage(function(event) {
#      cat("this is the last time i'll run\n")
#      removeThis()
#    })
#    ws$onOpen(function(event) {
#      ws$send("one")
#      ws$send("two")
#    })
#  }

