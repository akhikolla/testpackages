## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---- message=FALSE-----------------------------------------------------------
library(Ryacas)

## -----------------------------------------------------------------------------
eq <- "x^2 + 4 + 2*x + 2*x"

## -----------------------------------------------------------------------------
yac_str(eq) # No task was given to yacas, so we simply get the same returned
yac_str(paste0("Simplify(", eq, ")"))
yac_str(paste0("Factor(", eq, ")"))
yac_str(paste0("TeXForm(Factor(", eq, "))"))

## -----------------------------------------------------------------------------
y_fn(eq, "Simplify")
yac_str(y_fn(eq, "Simplify"))
yac_str(y_fn(eq, "Factor"))
yac_str(y_fn(y_fn(eq, "Factor"), "TeXForm"))

## -----------------------------------------------------------------------------
eq %>% y_fn("Simplify")
eq %>% y_fn("Simplify") %>% yac_str()
eq %>% y_fn("Factor") %>% yac_str()
eq %>% y_fn("Factor") %>% y_fn("TeXForm") %>% yac_str()

## -----------------------------------------------------------------------------
eq
eq %>% yac_expr() # Alternative to "yac_expr(eq)"
cmd <- eq %>% y_fn("Factor")
cmd
e <- yac_expr(cmd)
e
eval(e, list(x = 2))

## -----------------------------------------------------------------------------
cmd <- "2*{x, x^2, x^3}"
cmd %>% yac_str()
e <- cmd %>% yac_expr()
e
eval(e, list(x = 1.5))

## -----------------------------------------------------------------------------
cmd <- "{{1, 2}, {3, 4}}"
yac_str(cmd)
y_print(cmd) # Convenience function for prettier print
e <- cmd %>% yac_expr()
e
eval(e)

## -----------------------------------------------------------------------------
cmd %>% y_fn("TeXForm") %>% yac_str()

## -----------------------------------------------------------------------------
cmd1 <- paste0("a * ", cmd, "")
cmd2 <- cmd1 %>% y_fn("Inverse")
cmd2
cmd2 %>% yac_str()
cmd2 %>% y_fn("TeXForm") %>% yac_str()

## -----------------------------------------------------------------------------
paste0(cmd2, "*", cmd1) %>% 
  y_fn("Simplify") %>% 
  yac_str()
e2 <- cmd2 %>% yac_expr()
eval(e2, list(a = 2.2))

## -----------------------------------------------------------------------------
Achr <- matrix(0, nrow = 3, ncol = 3)
diag(Achr) <- 1
Achr[2, 3] <- "a"
Achr[1, 3] <- "a"
Achr

## -----------------------------------------------------------------------------
Ayac <- Achr %>% as_y()
Ayac

## -----------------------------------------------------------------------------
cmd <- Ayac %>% y_fn("Inverse")
cmd %>% yac_str()

## -----------------------------------------------------------------------------
way1 <- cmd %>% yac_str()
way1 %>% y_print()
way2 <- cmd %>% y_fn("TeXForm") %>% yac_str()
way2
way3 <- cmd %>% y_fn("PrettyForm") %>% yac_str()
way3 %>% cat() # Result of PrettyForm() must be printed
way4 <- cmd %>% yac_str()
way4
way4 %>% as_r()
way4 %>% as_r() %>% print(quote = FALSE)

## -----------------------------------------------------------------------------
A_inv_yac <- way4 %>% as_r()
Bchr <- A_inv_yac[2:3, 2:3]
Bchr
Bchr %>% as_y()
Bchr %>% as_y() %>% 
    y_fn("Inverse") %>% 
    yac_str() %>% 
    as_r()

## -----------------------------------------------------------------------------
yac_str("poly := (x-3)*(x+2)")

## -----------------------------------------------------------------------------
yac_silent("poly := (x-3)*(x+2)")

## -----------------------------------------------------------------------------
yac_str("Variables()")

## -----------------------------------------------------------------------------
yac_str("Expand(poly)")
"poly" %>% y_fn("Expand") %>% yac_str()

## -----------------------------------------------------------------------------
yac_str("Sum(k, 0, n, a^k)")

## -----------------------------------------------------------------------------
cmd <- "Limit(n, Infinity) (1+(1/n))^n"
yac_str(cmd)
yac_expr(cmd)

## -----------------------------------------------------------------------------
yac_str("Limit(h, 0) (Sin(x+h)-Sin(x))/h") 

## -----------------------------------------------------------------------------
cmd <- "Solve(poly == 0, x)"
cmd %>% yac_str()

## -----------------------------------------------------------------------------
cmd <- "Solve(poly, x)"
cmd %>% yac_str()

## -----------------------------------------------------------------------------
cmd
cmd %>% y_rmvars() %>% yac_str()
cmd %>% y_rmvars() %>% yac_expr()

## -----------------------------------------------------------------------------
"poly == 0" %>% y_fn("Solve", "x")
"poly" %>% y_fn("Solve", "x") # default is == 0
"poly == 0" %>% y_fn("Solve", "x") %>% y_rmvars() %>% yac_str()

## -----------------------------------------------------------------------------
f <- function(x, y) 2*x^2 + 3*y + y*x^2
f_body <- body(f)
f_body
f_body_chr <- as.character(as.expression(body(f)))
f_body_chr
# or:
# f_body_chr <- "2 * x^2 + 3 * y + y * x^2"
# f <- function(x, y) NULL
# body(f) <- parse(text = f_body_chr, keep.source = FALSE)

## -----------------------------------------------------------------------------
cmd_g <- paste0("{ D(x) ", f_body_chr, ", D(y) ", f_body_chr, " }")
cmd_g
g_body <- yac_expr(cmd_g)
g_body
g <- function(x, y) NULL
body(g) <- g_body
g
g(2, 4)

## -----------------------------------------------------------------------------
cmd_H <- paste0("HessianMatrix(", f_body_chr, ", {x, y})")
cmd_H
H_body <- yac_expr(cmd_H)
H_body
H <- function(x, y) NULL
body(H) <- H_body
H
H(2, 4)

