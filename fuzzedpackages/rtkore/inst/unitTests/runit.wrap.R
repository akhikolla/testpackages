#-----------------------------------------------------------------------
#     Copyright (C) 2012-2017  Serge Iovleff, University Lille 1, Inria
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as
#    published by the Free Software Foundation; either version 2 of the
#    License, or (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public
#    License along with this program; if not, write to the
#    Free Software Foundation, Inc.,
#    59 Temple Place,
#    Suite 330,
#    Boston, MA 02111-1307
#    USA
#
#    Contact : S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
#
#-----------------------------------------------------------------------
# test wrap method (mostly of myself, could be enhanced...)
#
# This file cannot be set in the tests directory (require the inline package)
#
if (require("inline"))
{
  body <- '
    RMatrix<double> m_rtab(tab1);
    RVector<double> v_rtab(tab2);

    // create stk arrays and test assignement
    ArrayXX                m_stk; m_stk = m_rtab;
    VectorX                v_stk; v_stk = v_rtab;
    PointX                 p_stk; p_stk = v_rtab;
    ArrayDiagonalX         d_stk; d_stk = v_rtab;
    CSquareX               m_stk2; m_stk2 = m_rtab;
    ArrayLowerTriangularXX l_stk; l_stk = m_rtab.lowerTriangularize();
    ArrayUpperTriangularXX u_stk; u_stk = m_rtab.upperTriangularize();

    // return wrapped arrays
    List ret;
    ret("RMatrix")    = STK::wrap(m_rtab);
    ret("RVector")    = STK::wrap(v_rtab);
    ret("STKMatrix")  = STK::wrap(m_stk);
    ret("STKSquare")  = STK::wrap(m_stk2);
    ret("STKVector")  = STK::wrap(v_stk);
    ret("STKPoint")   = STK::wrap(p_stk);
    ret("STKDiagonal")= STK::wrap(d_stk);
    ret("STKUpper")   = STK::wrap(u_stk);
    ret("STKLower")   = STK::wrap(l_stk);
    return ret;
  '

  wrap <- cxxfunction( signature(tab1 = "matrix", tab2 = "vector")
                   , body
                   , plugin = "rtkore"
                   , verbose = TRUE
                   )
  tab1 <- matrix(rnorm(64), ncol=8)
  tab2 <- rnorm(8)
  res_wrap <- wrap( tab1, tab2)
  res_wrap
}else
{
  print("Install the package inline in order to launch this test")
}
