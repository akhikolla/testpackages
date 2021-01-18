       subroutine rosen(n, x, fval)
       double precision x(n), fval, dx
       integer n, i
       fval = 0.0D0
       do 10 i=1,(n-1)
          dx = x(i + 1) - x(i) * x(i)
          fval = fval + 100.0 * dx * dx
          dx = 1.0 - x(i)
          fval = fval + dx * dx
 10    continue
       return
       end
