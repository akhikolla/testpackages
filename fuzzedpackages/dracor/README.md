
<!-- README.md is generated from README.Rmd. Please edit that file -->

# dracor

<!-- badges: start -->

[![Travis build
status](https://travis-ci.com/jefferis/dracor.svg?branch=master)](https://travis-ci.com/jefferis/dracor)
[![CRAN
status](https://www.r-pkg.org/badges/version/dracor)](https://CRAN.R-project.org/package=dracor)
[![GitHub](https://img.shields.io/github/v/release/jefferis/dracor)](https://github.com/jefferis/dracor/releases/)
[![Downloads](http://cranlogs.r-pkg.org/badges/dracor?color=brightgreen)](https://www.r-pkg.org/pkg/dracor)
<!-- badges: end -->

The goal of **dracor** is to allow decoding of the Draco compressed
meshes in R. This is done by wrapping the
[Draco](https://github.com/google/draco) C++ decoding library with the
assistance of the [Rcpp
package](https://cran.r-project.org/package=Rcpp).

The original motivation for **dracor** was decoding
[neuroglancer](https://github.com/google/neuroglancer) meshes of neurons
for example as used by <https://flywire.ai/>.

## Installation

**dracor** is available from
[CRAN](https://cran.r-project.org/package=dracor):

``` r
install.packages('dracor')
```

but you can also install the development version like so:

``` r
remotes::install_github("jefferis/dracor")
```

## Example

This is a basic example using a sample from the draco repository

``` r
library(dracor)
# get sample file from draco repository
carurl='https://github.com/google/draco/blob/master/testdata/car.drc?raw=true'
car.m3d=dracor::draco_decode(carurl)
str(car.m3d)
#> List of 2
#>  $ vb: num [1:4, 1:1856] 1.54 1.65 -1.21 1 1.57 ...
#>  $ it: int [1:3, 1:1744] 1 2 3 3 2 4 4 2 5 5 ...
#>  - attr(*, "class")= chr [1:2] "mesh3d" "shape3d"
```

[rgl](https://cran.r-project.org/package=rgl) is the most widely used R
package for 3D visualisation. By default we return meshes as rgl
`mesh3d` objects, which can then be displayed by `rgl` or manipulated by
a range of R packages including
[Rvcg](https://cran.r-project.org/package=Rvcg).

``` r
# install.packages("rgl")
# convert to rgl mesh3d format
rgl::shade3d(car.m3d, col='red')
# set a nice viewpoint
rgl::view3d(theta = 60, fov=0, zoom=.7)
```

<script>/*
* Copyright (C) 2009 Apple Inc. All Rights Reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions
* are met:
* 1. Redistributions of source code must retain the above copyright
*    notice, this list of conditions and the following disclaimer.
* 2. Redistributions in binary form must reproduce the above copyright
*    notice, this list of conditions and the following disclaimer in the
*    documentation and/or other materials provided with the distribution.
*
* THIS SOFTWARE IS PROVIDED BY APPLE INC. ``AS IS'' AND ANY
* EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
* PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL APPLE INC. OR
* CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
* EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
* PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
* PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
* OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
* OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
* Copyright (2016) Duncan Murdoch - fixed CanvasMatrix4.ortho,
* cleaned up.
*/
/*
CanvasMatrix4 class
This class implements a 4x4 matrix. It has functions which
duplicate the functionality of the OpenGL matrix stack and
glut functions.
IDL:
[
Constructor(in CanvasMatrix4 matrix),           // copy passed matrix into new CanvasMatrix4
Constructor(in sequence<float> array)           // create new CanvasMatrix4 with 16 floats (row major)
Constructor()                                   // create new CanvasMatrix4 with identity matrix
]
interface CanvasMatrix4 {
attribute float m11;
attribute float m12;
attribute float m13;
attribute float m14;
attribute float m21;
attribute float m22;
attribute float m23;
attribute float m24;
attribute float m31;
attribute float m32;
attribute float m33;
attribute float m34;
attribute float m41;
attribute float m42;
attribute float m43;
attribute float m44;
void load(in CanvasMatrix4 matrix);                 // copy the values from the passed matrix
void load(in sequence<float> array);                // copy 16 floats into the matrix
sequence<float> getAsArray();                       // return the matrix as an array of 16 floats
WebGLFloatArray getAsCanvasFloatArray();           // return the matrix as a WebGLFloatArray with 16 values
void makeIdentity();                                // replace the matrix with identity
void transpose();                                   // replace the matrix with its transpose
void invert();                                      // replace the matrix with its inverse
void translate(in float x, in float y, in float z); // multiply the matrix by passed translation values on the right
void scale(in float x, in float y, in float z);     // multiply the matrix by passed scale values on the right
void rotate(in float angle,                         // multiply the matrix by passed rotation values on the right
in float x, in float y, in float z);    // (angle is in degrees)
void multRight(in CanvasMatrix matrix);             // multiply the matrix by the passed matrix on the right
void multLeft(in CanvasMatrix matrix);              // multiply the matrix by the passed matrix on the left
void ortho(in float left, in float right,           // multiply the matrix by the passed ortho values on the right
in float bottom, in float top,
in float near, in float far);
void frustum(in float left, in float right,         // multiply the matrix by the passed frustum values on the right
in float bottom, in float top,
in float near, in float far);
void perspective(in float fovy, in float aspect,    // multiply the matrix by the passed perspective values on the right
in float zNear, in float zFar);
void lookat(in float eyex, in float eyey, in float eyez,    // multiply the matrix by the passed lookat
in float ctrx, in float ctry, in float ctrz,    // values on the right
in float upx, in float upy, in float upz);
}
*/
CanvasMatrix4 = function(m)
{
if (typeof m == 'object') {
if ("length" in m && m.length >= 16) {
this.load(m[0], m[1], m[2], m[3], m[4], m[5], m[6], m[7], m[8], m[9], m[10], m[11], m[12], m[13], m[14], m[15]);
return;
}
else if (m instanceof CanvasMatrix4) {
this.load(m);
return;
}
}
this.makeIdentity();
};
CanvasMatrix4.prototype.load = function()
{
if (arguments.length == 1 && typeof arguments[0] == 'object') {
var matrix = arguments[0];
if ("length" in matrix && matrix.length == 16) {
this.m11 = matrix[0];
this.m12 = matrix[1];
this.m13 = matrix[2];
this.m14 = matrix[3];
this.m21 = matrix[4];
this.m22 = matrix[5];
this.m23 = matrix[6];
this.m24 = matrix[7];
this.m31 = matrix[8];
this.m32 = matrix[9];
this.m33 = matrix[10];
this.m34 = matrix[11];
this.m41 = matrix[12];
this.m42 = matrix[13];
this.m43 = matrix[14];
this.m44 = matrix[15];
return;
}
if (arguments[0] instanceof CanvasMatrix4) {
this.m11 = matrix.m11;
this.m12 = matrix.m12;
this.m13 = matrix.m13;
this.m14 = matrix.m14;
this.m21 = matrix.m21;
this.m22 = matrix.m22;
this.m23 = matrix.m23;
this.m24 = matrix.m24;
this.m31 = matrix.m31;
this.m32 = matrix.m32;
this.m33 = matrix.m33;
this.m34 = matrix.m34;
this.m41 = matrix.m41;
this.m42 = matrix.m42;
this.m43 = matrix.m43;
this.m44 = matrix.m44;
return;
}
}
this.makeIdentity();
};
CanvasMatrix4.prototype.getAsArray = function()
{
return [
this.m11, this.m12, this.m13, this.m14,
this.m21, this.m22, this.m23, this.m24,
this.m31, this.m32, this.m33, this.m34,
this.m41, this.m42, this.m43, this.m44
];
};
CanvasMatrix4.prototype.getAsWebGLFloatArray = function()
{
return new WebGLFloatArray(this.getAsArray());
};
CanvasMatrix4.prototype.makeIdentity = function()
{
this.m11 = 1;
this.m12 = 0;
this.m13 = 0;
this.m14 = 0;
this.m21 = 0;
this.m22 = 1;
this.m23 = 0;
this.m24 = 0;
this.m31 = 0;
this.m32 = 0;
this.m33 = 1;
this.m34 = 0;
this.m41 = 0;
this.m42 = 0;
this.m43 = 0;
this.m44 = 1;
};
CanvasMatrix4.prototype.transpose = function()
{
var tmp = this.m12;
this.m12 = this.m21;
this.m21 = tmp;
tmp = this.m13;
this.m13 = this.m31;
this.m31 = tmp;
tmp = this.m14;
this.m14 = this.m41;
this.m41 = tmp;
tmp = this.m23;
this.m23 = this.m32;
this.m32 = tmp;
tmp = this.m24;
this.m24 = this.m42;
this.m42 = tmp;
tmp = this.m34;
this.m34 = this.m43;
this.m43 = tmp;
};
CanvasMatrix4.prototype.invert = function()
{
// Calculate the 4x4 determinant
// If the determinant is zero,
// then the inverse matrix is not unique.
var det = this._determinant4x4();
if (Math.abs(det) < 1e-8)
return null;
this._makeAdjoint();
// Scale the adjoint matrix to get the inverse
this.m11 /= det;
this.m12 /= det;
this.m13 /= det;
this.m14 /= det;
this.m21 /= det;
this.m22 /= det;
this.m23 /= det;
this.m24 /= det;
this.m31 /= det;
this.m32 /= det;
this.m33 /= det;
this.m34 /= det;
this.m41 /= det;
this.m42 /= det;
this.m43 /= det;
this.m44 /= det;
};
CanvasMatrix4.prototype.translate = function(x,y,z)
{
if (x === undefined)
x = 0;
if (y === undefined)
y = 0;
if (z === undefined)
z = 0;
var matrix = new CanvasMatrix4();
matrix.m41 = x;
matrix.m42 = y;
matrix.m43 = z;
this.multRight(matrix);
};
CanvasMatrix4.prototype.scale = function(x,y,z)
{
if (x === undefined)
x = 1;
if (z === undefined) {
if (y === undefined) {
y = x;
z = x;
}
else
z = 1;
}
else if (y === undefined)
y = x;
var matrix = new CanvasMatrix4();
matrix.m11 = x;
matrix.m22 = y;
matrix.m33 = z;
this.multRight(matrix);
};
CanvasMatrix4.prototype.rotate = function(angle,x,y,z)
{
// angles are in degrees. Switch to radians
angle = angle / 180 * Math.PI;
angle /= 2;
var sinA = Math.sin(angle);
var cosA = Math.cos(angle);
var sinA2 = sinA * sinA;
// normalize
var length = Math.sqrt(x * x + y * y + z * z);
if (length === 0) {
// bad vector, just use something reasonable
x = 0;
y = 0;
z = 1;
} else if (length != 1) {
x /= length;
y /= length;
z /= length;
}
var mat = new CanvasMatrix4();
// optimize case where axis is along major axis
if (x == 1 && y === 0 && z === 0) {
mat.m11 = 1;
mat.m12 = 0;
mat.m13 = 0;
mat.m21 = 0;
mat.m22 = 1 - 2 * sinA2;
mat.m23 = 2 * sinA * cosA;
mat.m31 = 0;
mat.m32 = -2 * sinA * cosA;
mat.m33 = 1 - 2 * sinA2;
mat.m14 = mat.m24 = mat.m34 = 0;
mat.m41 = mat.m42 = mat.m43 = 0;
mat.m44 = 1;
} else if (x === 0 && y == 1 && z === 0) {
mat.m11 = 1 - 2 * sinA2;
mat.m12 = 0;
mat.m13 = -2 * sinA * cosA;
mat.m21 = 0;
mat.m22 = 1;
mat.m23 = 0;
mat.m31 = 2 * sinA * cosA;
mat.m32 = 0;
mat.m33 = 1 - 2 * sinA2;
mat.m14 = mat.m24 = mat.m34 = 0;
mat.m41 = mat.m42 = mat.m43 = 0;
mat.m44 = 1;
} else if (x === 0 && y === 0 && z == 1) {
mat.m11 = 1 - 2 * sinA2;
mat.m12 = 2 * sinA * cosA;
mat.m13 = 0;
mat.m21 = -2 * sinA * cosA;
mat.m22 = 1 - 2 * sinA2;
mat.m23 = 0;
mat.m31 = 0;
mat.m32 = 0;
mat.m33 = 1;
mat.m14 = mat.m24 = mat.m34 = 0;
mat.m41 = mat.m42 = mat.m43 = 0;
mat.m44 = 1;
} else {
var x2 = x*x;
var y2 = y*y;
var z2 = z*z;
mat.m11 = 1 - 2 * (y2 + z2) * sinA2;
mat.m12 = 2 * (x * y * sinA2 + z * sinA * cosA);
mat.m13 = 2 * (x * z * sinA2 - y * sinA * cosA);
mat.m21 = 2 * (y * x * sinA2 - z * sinA * cosA);
mat.m22 = 1 - 2 * (z2 + x2) * sinA2;
mat.m23 = 2 * (y * z * sinA2 + x * sinA * cosA);
mat.m31 = 2 * (z * x * sinA2 + y * sinA * cosA);
mat.m32 = 2 * (z * y * sinA2 - x * sinA * cosA);
mat.m33 = 1 - 2 * (x2 + y2) * sinA2;
mat.m14 = mat.m24 = mat.m34 = 0;
mat.m41 = mat.m42 = mat.m43 = 0;
mat.m44 = 1;
}
this.multRight(mat);
};
CanvasMatrix4.prototype.multRight = function(mat)
{
var m11 = (this.m11 * mat.m11 + this.m12 * mat.m21 +
this.m13 * mat.m31 + this.m14 * mat.m41);
var m12 = (this.m11 * mat.m12 + this.m12 * mat.m22 +
this.m13 * mat.m32 + this.m14 * mat.m42);
var m13 = (this.m11 * mat.m13 + this.m12 * mat.m23 +
this.m13 * mat.m33 + this.m14 * mat.m43);
var m14 = (this.m11 * mat.m14 + this.m12 * mat.m24 +
this.m13 * mat.m34 + this.m14 * mat.m44);
var m21 = (this.m21 * mat.m11 + this.m22 * mat.m21 +
this.m23 * mat.m31 + this.m24 * mat.m41);
var m22 = (this.m21 * mat.m12 + this.m22 * mat.m22 +
this.m23 * mat.m32 + this.m24 * mat.m42);
var m23 = (this.m21 * mat.m13 + this.m22 * mat.m23 +
this.m23 * mat.m33 + this.m24 * mat.m43);
var m24 = (this.m21 * mat.m14 + this.m22 * mat.m24 +
this.m23 * mat.m34 + this.m24 * mat.m44);
var m31 = (this.m31 * mat.m11 + this.m32 * mat.m21 +
this.m33 * mat.m31 + this.m34 * mat.m41);
var m32 = (this.m31 * mat.m12 + this.m32 * mat.m22 +
this.m33 * mat.m32 + this.m34 * mat.m42);
var m33 = (this.m31 * mat.m13 + this.m32 * mat.m23 +
this.m33 * mat.m33 + this.m34 * mat.m43);
var m34 = (this.m31 * mat.m14 + this.m32 * mat.m24 +
this.m33 * mat.m34 + this.m34 * mat.m44);
var m41 = (this.m41 * mat.m11 + this.m42 * mat.m21 +
this.m43 * mat.m31 + this.m44 * mat.m41);
var m42 = (this.m41 * mat.m12 + this.m42 * mat.m22 +
this.m43 * mat.m32 + this.m44 * mat.m42);
var m43 = (this.m41 * mat.m13 + this.m42 * mat.m23 +
this.m43 * mat.m33 + this.m44 * mat.m43);
var m44 = (this.m41 * mat.m14 + this.m42 * mat.m24 +
this.m43 * mat.m34 + this.m44 * mat.m44);
this.m11 = m11;
this.m12 = m12;
this.m13 = m13;
this.m14 = m14;
this.m21 = m21;
this.m22 = m22;
this.m23 = m23;
this.m24 = m24;
this.m31 = m31;
this.m32 = m32;
this.m33 = m33;
this.m34 = m34;
this.m41 = m41;
this.m42 = m42;
this.m43 = m43;
this.m44 = m44;
};
CanvasMatrix4.prototype.multLeft = function(mat)
{
var m11 = (mat.m11 * this.m11 + mat.m12 * this.m21 +
mat.m13 * this.m31 + mat.m14 * this.m41);
var m12 = (mat.m11 * this.m12 + mat.m12 * this.m22 +
mat.m13 * this.m32 + mat.m14 * this.m42);
var m13 = (mat.m11 * this.m13 + mat.m12 * this.m23 +
mat.m13 * this.m33 + mat.m14 * this.m43);
var m14 = (mat.m11 * this.m14 + mat.m12 * this.m24 +
mat.m13 * this.m34 + mat.m14 * this.m44);
var m21 = (mat.m21 * this.m11 + mat.m22 * this.m21 +
mat.m23 * this.m31 + mat.m24 * this.m41);
var m22 = (mat.m21 * this.m12 + mat.m22 * this.m22 +
mat.m23 * this.m32 + mat.m24 * this.m42);
var m23 = (mat.m21 * this.m13 + mat.m22 * this.m23 +
mat.m23 * this.m33 + mat.m24 * this.m43);
var m24 = (mat.m21 * this.m14 + mat.m22 * this.m24 +
mat.m23 * this.m34 + mat.m24 * this.m44);
var m31 = (mat.m31 * this.m11 + mat.m32 * this.m21 +
mat.m33 * this.m31 + mat.m34 * this.m41);
var m32 = (mat.m31 * this.m12 + mat.m32 * this.m22 +
mat.m33 * this.m32 + mat.m34 * this.m42);
var m33 = (mat.m31 * this.m13 + mat.m32 * this.m23 +
mat.m33 * this.m33 + mat.m34 * this.m43);
var m34 = (mat.m31 * this.m14 + mat.m32 * this.m24 +
mat.m33 * this.m34 + mat.m34 * this.m44);
var m41 = (mat.m41 * this.m11 + mat.m42 * this.m21 +
mat.m43 * this.m31 + mat.m44 * this.m41);
var m42 = (mat.m41 * this.m12 + mat.m42 * this.m22 +
mat.m43 * this.m32 + mat.m44 * this.m42);
var m43 = (mat.m41 * this.m13 + mat.m42 * this.m23 +
mat.m43 * this.m33 + mat.m44 * this.m43);
var m44 = (mat.m41 * this.m14 + mat.m42 * this.m24 +
mat.m43 * this.m34 + mat.m44 * this.m44);
this.m11 = m11;
this.m12 = m12;
this.m13 = m13;
this.m14 = m14;
this.m21 = m21;
this.m22 = m22;
this.m23 = m23;
this.m24 = m24;
this.m31 = m31;
this.m32 = m32;
this.m33 = m33;
this.m34 = m34;
this.m41 = m41;
this.m42 = m42;
this.m43 = m43;
this.m44 = m44;
};
CanvasMatrix4.prototype.ortho = function(left, right, bottom, top, near, far)
{
var tx = (left + right) / (left - right);
var ty = (top + bottom) / (bottom - top);
var tz = (far + near) / (near - far);
var matrix = new CanvasMatrix4();
matrix.m11 = 2 / (right - left);
matrix.m12 = 0;
matrix.m13 = 0;
matrix.m14 = 0;
matrix.m21 = 0;
matrix.m22 = 2 / (top - bottom);
matrix.m23 = 0;
matrix.m24 = 0;
matrix.m31 = 0;
matrix.m32 = 0;
matrix.m33 = -2 / (far - near);
matrix.m34 = 0;
matrix.m41 = tx;
matrix.m42 = ty;
matrix.m43 = tz;
matrix.m44 = 1;
this.multRight(matrix);
};
CanvasMatrix4.prototype.frustum = function(left, right, bottom, top, near, far)
{
var matrix = new CanvasMatrix4();
var A = (right + left) / (right - left);
var B = (top + bottom) / (top - bottom);
var C = -(far + near) / (far - near);
var D = -(2 * far * near) / (far - near);
matrix.m11 = (2 * near) / (right - left);
matrix.m12 = 0;
matrix.m13 = 0;
matrix.m14 = 0;
matrix.m21 = 0;
matrix.m22 = 2 * near / (top - bottom);
matrix.m23 = 0;
matrix.m24 = 0;
matrix.m31 = A;
matrix.m32 = B;
matrix.m33 = C;
matrix.m34 = -1;
matrix.m41 = 0;
matrix.m42 = 0;
matrix.m43 = D;
matrix.m44 = 0;
this.multRight(matrix);
};
CanvasMatrix4.prototype.perspective = function(fovy, aspect, zNear, zFar)
{
var top = Math.tan(fovy * Math.PI / 360) * zNear;
var bottom = -top;
var left = aspect * bottom;
var right = aspect * top;
this.frustum(left, right, bottom, top, zNear, zFar);
};
CanvasMatrix4.prototype.lookat = function(eyex, eyey, eyez, centerx, centery, centerz, upx, upy, upz)
{
var matrix = new CanvasMatrix4();
// Make rotation matrix
// Z vector
var zx = eyex - centerx;
var zy = eyey - centery;
var zz = eyez - centerz;
var mag = Math.sqrt(zx * zx + zy * zy + zz * zz);
if (mag) {
zx /= mag;
zy /= mag;
zz /= mag;
}
// Y vector
var yx = upx;
var yy = upy;
var yz = upz;
// X vector = Y cross Z
xx =  yy * zz - yz * zy;
xy = -yx * zz + yz * zx;
xz =  yx * zy - yy * zx;
// Recompute Y = Z cross X
yx = zy * xz - zz * xy;
yy = -zx * xz + zz * xx;
yx = zx * xy - zy * xx;
// cross product gives area of parallelogram, which is < 1.0 for
// non-perpendicular unit-length vectors; so normalize x, y here
mag = Math.sqrt(xx * xx + xy * xy + xz * xz);
if (mag) {
xx /= mag;
xy /= mag;
xz /= mag;
}
mag = Math.sqrt(yx * yx + yy * yy + yz * yz);
if (mag) {
yx /= mag;
yy /= mag;
yz /= mag;
}
matrix.m11 = xx;
matrix.m12 = xy;
matrix.m13 = xz;
matrix.m14 = 0;
matrix.m21 = yx;
matrix.m22 = yy;
matrix.m23 = yz;
matrix.m24 = 0;
matrix.m31 = zx;
matrix.m32 = zy;
matrix.m33 = zz;
matrix.m34 = 0;
matrix.m41 = 0;
matrix.m42 = 0;
matrix.m43 = 0;
matrix.m44 = 1;
matrix.translate(-eyex, -eyey, -eyez);
this.multRight(matrix);
};
// Support functions
CanvasMatrix4.prototype._determinant2x2 = function(a, b, c, d)
{
return a * d - b * c;
};
CanvasMatrix4.prototype._determinant3x3 = function(a1, a2, a3, b1, b2, b3, c1, c2, c3)
{
return a1 * this._determinant2x2(b2, b3, c2, c3) -
b1 * this._determinant2x2(a2, a3, c2, c3) +
c1 * this._determinant2x2(a2, a3, b2, b3);
};
CanvasMatrix4.prototype._determinant4x4 = function()
{
var a1 = this.m11;
var b1 = this.m12;
var c1 = this.m13;
var d1 = this.m14;
var a2 = this.m21;
var b2 = this.m22;
var c2 = this.m23;
var d2 = this.m24;
var a3 = this.m31;
var b3 = this.m32;
var c3 = this.m33;
var d3 = this.m34;
var a4 = this.m41;
var b4 = this.m42;
var c4 = this.m43;
var d4 = this.m44;
return a1 * this._determinant3x3(b2, b3, b4, c2, c3, c4, d2, d3, d4) -
b1 * this._determinant3x3(a2, a3, a4, c2, c3, c4, d2, d3, d4) +
c1 * this._determinant3x3(a2, a3, a4, b2, b3, b4, d2, d3, d4) -
d1 * this._determinant3x3(a2, a3, a4, b2, b3, b4, c2, c3, c4);
};
CanvasMatrix4.prototype._makeAdjoint = function()
{
var a1 = this.m11;
var b1 = this.m12;
var c1 = this.m13;
var d1 = this.m14;
var a2 = this.m21;
var b2 = this.m22;
var c2 = this.m23;
var d2 = this.m24;
var a3 = this.m31;
var b3 = this.m32;
var c3 = this.m33;
var d3 = this.m34;
var a4 = this.m41;
var b4 = this.m42;
var c4 = this.m43;
var d4 = this.m44;
// Row column labeling reversed since we transpose rows & columns
this.m11  =   this._determinant3x3(b2, b3, b4, c2, c3, c4, d2, d3, d4);
this.m21  = - this._determinant3x3(a2, a3, a4, c2, c3, c4, d2, d3, d4);
this.m31  =   this._determinant3x3(a2, a3, a4, b2, b3, b4, d2, d3, d4);
this.m41  = - this._determinant3x3(a2, a3, a4, b2, b3, b4, c2, c3, c4);
this.m12  = - this._determinant3x3(b1, b3, b4, c1, c3, c4, d1, d3, d4);
this.m22  =   this._determinant3x3(a1, a3, a4, c1, c3, c4, d1, d3, d4);
this.m32  = - this._determinant3x3(a1, a3, a4, b1, b3, b4, d1, d3, d4);
this.m42  =   this._determinant3x3(a1, a3, a4, b1, b3, b4, c1, c3, c4);
this.m13  =   this._determinant3x3(b1, b2, b4, c1, c2, c4, d1, d2, d4);
this.m23  = - this._determinant3x3(a1, a2, a4, c1, c2, c4, d1, d2, d4);
this.m33  =   this._determinant3x3(a1, a2, a4, b1, b2, b4, d1, d2, d4);
this.m43  = - this._determinant3x3(a1, a2, a4, b1, b2, b4, c1, c2, c4);
this.m14  = - this._determinant3x3(b1, b2, b3, c1, c2, c3, d1, d2, d3);
this.m24  =   this._determinant3x3(a1, a2, a3, c1, c2, c3, d1, d2, d3);
this.m34  = - this._determinant3x3(a1, a2, a3, b1, b2, b3, d1, d2, d3);
this.m44  =   this._determinant3x3(a1, a2, a3, b1, b2, b3, c1, c2, c3);
};</script>

<script>// To generate the help pages for this library, use
// jsdoc --destination ../../../doc/rglwidgetClass --template ~/node_modules/jsdoc-baseline rglClass.src.js
// To validate, use
// setwd(".../inst/htmlwidgets/lib/rglClass")
// hints <- js::jshint(readLines("rglClass.src.js"))
// hints[, c("line", "reason")]
/**
* The class of an rgl widget
* @class
*/
rglwidgetClass = function() {
this.canvas = null;
this.userMatrix = new CanvasMatrix4();
this.types = [];
this.prMatrix = new CanvasMatrix4();
this.mvMatrix = new CanvasMatrix4();
this.vp = null;
this.prmvMatrix = null;
this.origs = null;
this.gl = null;
this.scene = null;
this.select = {state: "inactive", subscene: null, region: {p1: {x:0, y:0}, p2: {x:0, y:0}}};
this.drawing = false;
};
/**
* Multiply matrix by vector
* @returns {number[]}
* @param M {number[][]} Left operand
* @param v {number[]} Right operand
*/
rglwidgetClass.prototype.multMV = function(M, v) {
return [ M.m11 * v[0] + M.m12 * v[1] + M.m13 * v[2] + M.m14 * v[3],
M.m21 * v[0] + M.m22 * v[1] + M.m23 * v[2] + M.m24 * v[3],
M.m31 * v[0] + M.m32 * v[1] + M.m33 * v[2] + M.m34 * v[3],
M.m41 * v[0] + M.m42 * v[1] + M.m43 * v[2] + M.m44 * v[3]
];
};
/**
* Multiply row vector by Matrix
* @returns {number[]}
* @param v {number[]} left operand
* @param M {number[][]} right operand
*/
rglwidgetClass.prototype.multVM = function(v, M) {
return [ M.m11 * v[0] + M.m21 * v[1] + M.m31 * v[2] + M.m41 * v[3],
M.m12 * v[0] + M.m22 * v[1] + M.m32 * v[2] + M.m42 * v[3],
M.m13 * v[0] + M.m23 * v[1] + M.m33 * v[2] + M.m43 * v[3],
M.m14 * v[0] + M.m24 * v[1] + M.m34 * v[2] + M.m44 * v[3]
];
};
/**
* Euclidean length of a vector
* @returns {number}
* @param v {number[]}
*/
rglwidgetClass.prototype.vlen = function(v) {
return Math.sqrt(this.dotprod(v, v));
};
/**
* Dot product of two vectors
* @instance rglwidgetClass
* @returns {number}
* @param a {number[]}
* @param b {number[]}
*/
rglwidgetClass.prototype.dotprod = function(a, b) {
return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
};
/**
* Cross product of two vectors
* @returns {number[]}
* @param a {number[]}
* @param b {number[]}
*/
rglwidgetClass.prototype.xprod = function(a, b) {
return [a[1]*b[2] - a[2]*b[1],
a[2]*b[0] - a[0]*b[2],
a[0]*b[1] - a[1]*b[0]];
};
/**
* Bind vectors or matrices by columns
* @returns {number[][]}
* @param a {number[]|number[][]}
* @param b {number[]|number[][]}
*/
rglwidgetClass.prototype.cbind = function(a, b) {
if (b.length < a.length)
b = this.repeatToLen(b, a.length);
else if (a.length < b.length)
a = this.repeatToLen(a, b.length);
return a.map(function(currentValue, index, array) {
return currentValue.concat(b[index]);
});
};
/**
* Swap elements
* @returns {any[]}
* @param a {any[]}
* @param i {number} Element to swap
* @param j {number} Other element to swap
*/
rglwidgetClass.prototype.swap = function(a, i, j) {
var temp = a[i];
a[i] = a[j];
a[j] = temp;
};
/**
* Flatten a matrix into a vector
* @returns {any[]}
* @param a {any[][]}
*/
rglwidgetClass.prototype.flatten = function(arr, result) {
var value;
if (typeof result === "undefined") result = [];
for (var i = 0, length = arr.length; i < length; i++) {
value = arr[i];
if (Array.isArray(value)) {
this.flatten(value, result);
} else {
result.push(value);
}
}
return result;
};
/**
* set element of 1d or 2d array as if it was flattened.
* Column major, zero based!
* @returns {any[]|any[][]}
* @param {any[]|any[][]} a - array
* @param {number} i - element
* @param {any} value
*/
rglwidgetClass.prototype.setElement = function(a, i, value) {
if (Array.isArray(a[0])) {
var dim = a.length,
col = Math.floor(i/dim),
row = i % dim;
a[row][col] = value;
} else {
a[i] = value;
}
};
/**
* Transpose an array
* @returns {any[][]}
* @param {any[][]} a
*/
rglwidgetClass.prototype.transpose = function(a) {
var newArray = [],
n = a.length,
m = a[0].length,
i;
for(i = 0; i < m; i++){
newArray.push([]);
}
for(i = 0; i < n; i++){
for(var j = 0; j < m; j++){
newArray[j].push(a[i][j]);
}
}
return newArray;
};
/**
* Calculate sum of squares of a numeric vector
* @returns {number}
* @param {number[]} x
*/
rglwidgetClass.prototype.sumsq = function(x) {
var result = 0, i;
for (i=0; i < x.length; i++)
result += x[i]*x[i];
return result;
};
/**
* Convert a matrix to a CanvasMatrix4
* @returns {CanvasMatrix4}
* @param {number[][]|number[]} mat
*/
rglwidgetClass.prototype.toCanvasMatrix4 = function(mat) {
if (mat instanceof CanvasMatrix4)
return mat;
var result = new CanvasMatrix4();
mat = this.flatten(this.transpose(mat));
result.load(mat);
return result;
};
/**
* Convert an R-style numeric colour string to an rgb vector
* @returns {number[]}
* @param {string} s
*/
rglwidgetClass.prototype.stringToRgb = function(s) {
s = s.replace("#", "");
var bigint = parseInt(s, 16);
return [((bigint >> 16) & 255)/255,
((bigint >> 8) & 255)/255,
(bigint & 255)/255];
};
/**
* Take a component-by-component product of two 3 vectors
* @returns {number[]}
* @param {number[]} x
* @param {number[]} y
*/
rglwidgetClass.prototype.componentProduct = function(x, y) {
if (typeof y === "undefined") {
this.alertOnce("Bad arg to componentProduct");
}
var result = new Float32Array(3), i;
for (i = 0; i<3; i++)
result[i] = x[i]*y[i];
return result;
};
/**
* Get next higher power of two
* @returns { number }
* @param { number } value - input value
*/
rglwidgetClass.prototype.getPowerOfTwo = function(value) {
var pow = 1;
while(pow<value) {
pow *= 2;
}
return pow;
};
/**
* Unique entries
* @returns { any[] }
* @param { any[] } arr - An array
*/
rglwidgetClass.prototype.unique = function(arr) {
arr = [].concat(arr);
return arr.filter(function(value, index, self) {
return self.indexOf(value) === index;
});
};
/**
* Shallow compare of arrays
* @returns { boolean }
* @param { any[] } a - An array
* @param { any[] } b - Another array
*/
rglwidgetClass.prototype.equalArrays = function(a, b) {
return a === b || (a && b &&
a.length === b.length &&
a.every(function(v, i) {return v === b[i];}));
};
/**
* Repeat an array to a desired length
* @returns {any[]}
* @param {any | any[]} arr The input array
* @param {number} len The desired output length
*/
rglwidgetClass.prototype.repeatToLen = function(arr, len) {
arr = [].concat(arr);
while (arr.length < len/2)
arr = arr.concat(arr);
return arr.concat(arr.slice(0, len - arr.length));
};
/**
* Give a single alert message, not to be repeated.
* @param {string} msg  The message to give.
*/
rglwidgetClass.prototype.alertOnce = function(msg) {
if (typeof this.alerted !== "undefined")
return;
this.alerted = true;
alert(msg);
};
rglwidgetClass.prototype.f_is_lit = 1;
rglwidgetClass.prototype.f_is_smooth = 2;
rglwidgetClass.prototype.f_has_texture = 4;
rglwidgetClass.prototype.f_depth_sort = 8;
rglwidgetClass.prototype.f_fixed_quads = 16;
rglwidgetClass.prototype.f_is_transparent = 32;
rglwidgetClass.prototype.f_is_lines = 64;
rglwidgetClass.prototype.f_sprites_3d = 128;
rglwidgetClass.prototype.f_sprite_3d = 256;
rglwidgetClass.prototype.f_is_subscene = 512;
rglwidgetClass.prototype.f_is_clipplanes = 1024;
rglwidgetClass.prototype.f_fixed_size = 2048;
rglwidgetClass.prototype.f_is_points = 4096;
rglwidgetClass.prototype.f_is_twosided = 8192;
rglwidgetClass.prototype.f_fat_lines = 16384;
rglwidgetClass.prototype.f_is_brush = 32768;
/**
* Which list does a particular id come from?
* @returns { string }
* @param {number} id The id to look up.
*/
rglwidgetClass.prototype.whichList = function(id) {
var obj = this.getObj(id),
flags = obj.flags;
if (obj.type === "light")
return "lights";
if (flags & this.f_is_subscene)
return "subscenes";
if (flags & this.f_is_clipplanes)
return "clipplanes";
if (flags & this.f_is_transparent)
return "transparent";
return "opaque";
};
/**
* Get an object by id number.
* @returns { Object }
* @param {number} id
*/
rglwidgetClass.prototype.getObj = function(id) {
if (typeof id !== "number") {
this.alertOnce("getObj id is "+typeof id);
}
return this.scene.objects[id];
};
/**
* Get ids of a particular type from a subscene or the whole scene
* @returns { number[] }
* @param {string} type What type of object?
* @param {number} subscene  Which subscene?  If not given, find in the whole scene
*/
rglwidgetClass.prototype.getIdsByType = function(type, subscene) {
var
result = [], i, self = this;
if (typeof subscene === "undefined") {
Object.keys(this.scene.objects).forEach(
function(key) {
key = parseInt(key, 10);
if (self.getObj(key).type === type)
result.push(key);
});
} else {
ids = this.getObj(subscene).objects;
for (i=0; i < ids.length; i++) {
if (this.getObj(ids[i]).type === type) {
result.push(ids[i]);
}
}
}
return result;
};
/**
* Get a particular material property for an id
* @returns { any }
* @param {number} id  Which object?
* @param {string} property Which material property?
*/
rglwidgetClass.prototype.getMaterial = function(id, property) {
var obj = this.getObj(id), mat;
if (typeof obj.material === "undefined")
console.error("material undefined");
mat = obj.material[property];
if (typeof mat === "undefined")
mat = this.scene.material[property];
return mat;
};
/**
* Is a particular id in a subscene?
* @returns { boolean }
* @param {number} id Which id?
* @param {number} subscene Which subscene id?
*/
rglwidgetClass.prototype.inSubscene = function(id, subscene) {
return this.getObj(subscene).objects.indexOf(id) > -1;
};
/**
* Add an id to a subscene.
* @param {number} id Which id?
* @param {number} subscene Which subscene id?
*/
rglwidgetClass.prototype.addToSubscene = function(id, subscene) {
var thelist,
thesub = this.getObj(subscene),
ids = [id],
obj = this.getObj(id), i;
if (typeof obj != "undefined" && typeof (obj.newIds) !== "undefined") {
ids = ids.concat(obj.newIds);
}
thesub.objects = [].concat(thesub.objects);
for (i = 0; i < ids.length; i++) {
id = ids[i];
if (thesub.objects.indexOf(id) == -1) {
thelist = this.whichList(id);
thesub.objects.push(id);
thesub[thelist].push(id);
}
}
};
/**
* Delete an id from a subscene
* @param { number } id - the id to add
* @param { number } subscene - the id of the subscene
*/
rglwidgetClass.prototype.delFromSubscene = function(id, subscene) {
var thelist,
thesub = this.getObj(subscene),
obj = this.getObj(id),
ids = [id], i;
if (typeof obj !== "undefined" && typeof (obj.newIds) !== "undefined")
ids = ids.concat(obj.newIds);
thesub.objects = [].concat(thesub.objects); // It might be a scalar
for (j=0; j<ids.length;j++) {
id = ids[j];
i = thesub.objects.indexOf(id);
if (i > -1) {
thesub.objects.splice(i, 1);
thelist = this.whichList(id);
i = thesub[thelist].indexOf(id);
thesub[thelist].splice(i, 1);
}
}
};
/**
* Set the ids in a subscene
* @param { number[] } ids - the ids to set
* @param { number } subsceneid - the id of the subscene
*/
rglwidgetClass.prototype.setSubsceneEntries = function(ids, subsceneid) {
var sub = this.getObj(subsceneid);
sub.objects = ids;
this.initSubscene(subsceneid);
};
/**
* Get the ids in a subscene
* @returns {number[]}
* @param { number } subscene - the id of the subscene
*/
rglwidgetClass.prototype.getSubsceneEntries = function(subscene) {
return this.getObj(subscene).objects;
};
/**
* Get the ids of the subscenes within a subscene
* @returns { number[] }
* @param { number } subscene - the id of the subscene
*/
rglwidgetClass.prototype.getChildSubscenes = function(subscene) {
return this.getObj(subscene).subscenes;
};
/**
* Start drawing
* @returns { boolean } Previous state
*/
rglwidgetClass.prototype.startDrawing = function() {
var value = this.drawing;
this.drawing = true;
return value;
};
/**
* Stop drawing and check for context loss
* @param { boolean } saved - Previous state
*/
rglwidgetClass.prototype.stopDrawing = function(saved) {
this.drawing = saved;
if (!saved && this.gl && this.gl.isContextLost())
this.restartCanvas();
};
/**
* Generate the vertex shader for an object
* @returns {string}
* @param { number } id - Id of object
*/
rglwidgetClass.prototype.getVertexShader = function(id) {
var obj = this.getObj(id),
userShader = obj.userVertexShader,
flags = obj.flags,
type = obj.type,
is_lit = flags & this.f_is_lit,
has_texture = flags & this.f_has_texture,
fixed_quads = flags & this.f_fixed_quads,
sprites_3d = flags & this.f_sprites_3d,
sprite_3d = flags & this.f_sprite_3d,
nclipplanes = this.countClipplanes(),
fixed_size = flags & this.f_fixed_size,
is_points = flags & this.f_is_points,
is_twosided = flags & this.f_is_twosided,
fat_lines = flags & this.f_fat_lines,
is_brush = flags & this.f_is_brush,
result;
if (type === "clipplanes" || sprites_3d) return;
if (typeof userShader !== "undefined") return userShader;
result = "  /* ****** "+type+" object "+id+" vertex shader ****** */\n"+
"  attribute vec3 aPos;\n"+
"  attribute vec4 aCol;\n"+
" uniform mat4 mvMatrix;\n"+
" uniform mat4 prMatrix;\n"+
" varying vec4 vCol;\n"+
" varying vec4 vPosition;\n";
if ((is_lit && !fixed_quads && !is_brush) || sprite_3d)
result = result + "  attribute vec3 aNorm;\n"+
" uniform mat4 normMatrix;\n"+
" varying vec3 vNormal;\n";
if (has_texture || type === "text")
result = result + " attribute vec2 aTexcoord;\n"+
" varying vec2 vTexcoord;\n";
if (fixed_size)
result = result + "  uniform vec2 textScale;\n";
if (fixed_quads)
result = result + "  attribute vec2 aOfs;\n";
else if (sprite_3d)
result = result + "  uniform vec3 uOrig;\n"+
"  uniform float uSize;\n"+
"  uniform mat4 usermat;\n";
if (is_twosided)
result = result + "  attribute vec3 aPos1;\n"+
"  attribute vec3 aPos2;\n"+
"  varying float normz;\n";
if (fat_lines) {
result = result +   "  attribute vec3 aNext;\n"+
"  attribute vec2 aPoint;\n"+
"  varying vec2 vPoint;\n"+
"  varying float vLength;\n"+
"  uniform float uAspect;\n"+
"  uniform float uLwd;\n";
}
result = result + "  void main(void) {\n";
if ((nclipplanes || (!fixed_quads && !sprite_3d)) && !is_brush)
result = result + "    vPosition = mvMatrix * vec4(aPos, 1.);\n";
if (!fixed_quads && !sprite_3d && !is_brush)
result = result + "    gl_Position = prMatrix * vPosition;\n";
if (is_points) {
var size = this.getMaterial(id, "size");
result = result + "    gl_PointSize = "+size.toFixed(1)+";\n";
}
result = result + "    vCol = aCol;\n";
if (is_lit && !fixed_quads && !sprite_3d && !is_brush)
result = result + "    vNormal = normalize((normMatrix * vec4(aNorm, 1.)).xyz);\n";
if (has_texture || type == "text")
result = result + "    vTexcoord = aTexcoord;\n";
if (fixed_size)
result = result + "    vec4 pos = prMatrix * mvMatrix * vec4(aPos, 1.);\n"+
"   pos = pos/pos.w;\n"+
"   gl_Position = pos + vec4(aOfs*textScale, 0.,0.);\n";
if (type == "sprites" && !fixed_size)
result = result + "    vec4 pos = mvMatrix * vec4(aPos, 1.);\n"+
"   pos = pos/pos.w + vec4(aOfs, 0., 0.);\n"+
"   gl_Position = prMatrix*pos;\n";
if (sprite_3d)
result = result + "   vNormal = normalize((normMatrix * vec4(aNorm, 1.)).xyz);\n"+
"   vec4 pos = mvMatrix * vec4(uOrig, 1.);\n"+
"   vPosition = pos/pos.w + vec4(uSize*(vec4(aPos, 1.)*usermat).xyz,0.);\n"+
"   gl_Position = prMatrix * vPosition;\n";
if (is_twosided)
result = result + "   vec4 pos1 = prMatrix*(mvMatrix*vec4(aPos1, 1.));\n"+
"   pos1 = pos1/pos1.w - gl_Position/gl_Position.w;\n"+
"   vec4 pos2 = prMatrix*(mvMatrix*vec4(aPos2, 1.));\n"+
"   pos2 = pos2/pos2.w - gl_Position/gl_Position.w;\n"+
"   normz = pos1.x*pos2.y - pos1.y*pos2.x;\n";
if (fat_lines) 
/* This code was inspired by Matt Deslauriers' code in https://mattdesl.svbtle.com/drawing-lines-is-hard */
result = result + "   vec2 aspectVec = vec2(uAspect, 1.0);\n"+
"   mat4 projViewModel = prMatrix * mvMatrix;\n"+
"   vec4 currentProjected = projViewModel * vec4(aPos, 1.0);\n"+
"   currentProjected = currentProjected/currentProjected.w;\n"+
"   vec4 nextProjected = projViewModel * vec4(aNext, 1.0);\n"+
"   vec2 currentScreen = currentProjected.xy * aspectVec;\n"+
"   vec2 nextScreen = (nextProjected.xy / nextProjected.w) * aspectVec;\n"+
"   float len = uLwd;\n"+
"   vec2 dir = vec2(1.0, 0.0);\n"+
"   vPoint = aPoint;\n"+
"   vLength = length(nextScreen - currentScreen)/2.0;\n"+
"   vLength = vLength/(vLength + len);\n"+
"   if (vLength > 0.0) {\n"+
"     dir = normalize(nextScreen - currentScreen);\n"+
"   }\n"+
"   vec2 normal = vec2(-dir.y, dir.x);\n"+
"   dir.x /= uAspect;\n"+
"   normal.x /= uAspect;\n"+
"   vec4 offset = vec4(len*(normal*aPoint.x*aPoint.y - dir), 0.0, 0.0);\n"+
"   gl_Position = currentProjected + offset;\n";
if (is_brush)
result = result + "   gl_Position = vec4(aPos, 1.);\n";
result = result + "  }\n";
// console.log(result);
return result;
};
/**
* Generate the fragment shader for an object
* @returns {string}
* @param { number } id - Id of object
*/
rglwidgetClass.prototype.getFragmentShader = function(id) {
var obj = this.getObj(id),
userShader = obj.userFragmentShader,
flags = obj.flags,
type = obj.type,
is_lit = flags & this.f_is_lit,
has_texture = flags & this.f_has_texture,
fixed_quads = flags & this.f_fixed_quads,
sprites_3d = flags & this.f_sprites_3d,
is_twosided = (flags & this.f_is_twosided) > 0,
fat_lines = flags & this.f_fat_lines,
is_transparent = flags & this.f_is_transparent,
nclipplanes = this.countClipplanes(), i,
texture_format, nlights,
result;
if (type === "clipplanes" || sprites_3d) return;
if (typeof userShader !== "undefined") return userShader;
if (has_texture)
texture_format = this.getMaterial(id, "textype");
result = "/* ****** "+type+" object "+id+" fragment shader ****** */\n"+
"#ifdef GL_ES\n"+
"#ifdef GL_FRAGMENT_PRECISION_HIGH\n"+
"  precision highp float;\n"+
"#else\n"+
"  precision mediump float;\n"+
"#endif\n"+
"#endif\n"+
"  varying vec4 vCol; // carries alpha\n"+
"  varying vec4 vPosition;\n";
if (has_texture || type === "text")
result = result + "  varying vec2 vTexcoord;\n"+
" uniform sampler2D uSampler;\n";
if (is_lit && !fixed_quads)
result = result + "  varying vec3 vNormal;\n";
for (i = 0; i < nclipplanes; i++)
result = result + "  uniform vec4 vClipplane"+i+";\n";
if (is_lit) {
nlights = this.countLights();
if (nlights)
result = result + "  uniform mat4 mvMatrix;\n";
else
is_lit = false;
}
if (is_lit) {
result = result + "   uniform vec3 emission;\n"+
"   uniform float shininess;\n";
for (i=0; i < nlights; i++) {
result = result + "   uniform vec3 ambient" + i + ";\n"+
"   uniform vec3 specular" + i +"; // light*material\n"+
"   uniform vec3 diffuse" + i + ";\n"+
"   uniform vec3 lightDir" + i + ";\n"+
"   uniform bool viewpoint" + i + ";\n"+
"   uniform bool finite" + i + ";\n";
}
}
if (is_twosided)
result = result + "   uniform bool front;\n"+
"   varying float normz;\n";
if (fat_lines)
result = result + "   varying vec2 vPoint;\n"+
"   varying float vLength;\n";
result = result + "  void main(void) {\n";
if (fat_lines) {
result = result + "    vec2 point = vPoint;\n"+
"    bool neg = point.y < 0.0;\n"+
"    point.y = neg ? "+
"      (point.y + vLength)/(1.0 - vLength) :\n"+
"     -(point.y - vLength)/(1.0 - vLength);\n";
if (is_transparent && type == "linestrip")
result = result+"    if (neg && length(point) <= 1.0) discard;\n";
result = result + "    point.y = min(point.y, 0.0);\n"+
"    if (length(point) > 1.0) discard;\n";
}
for (i=0; i < nclipplanes;i++)
result = result + "    if (dot(vPosition, vClipplane"+i+") < 0.0) discard;\n";
if (fixed_quads) {
result = result +   "    vec3 n = vec3(0., 0., 1.);\n";
} else if (is_lit) {
result = result +   "    vec3 n = normalize(vNormal);\n";
}
if (is_twosided) {
result = result +   "    if ((normz <= 0.) != front) discard;\n";
}
if (is_lit) {
result = result + "    vec3 eye = normalize(-vPosition.xyz);\n"+
"   vec3 lightdir;\n"+
"   vec4 colDiff;\n"+
"   vec3 halfVec;\n"+
"   vec4 lighteffect = vec4(emission, 0.);\n"+
"   vec3 col;\n"+
"   float nDotL;\n";
if (!fixed_quads) {
result = result +   "   n = -faceforward(n, n, eye);\n";
}
for (i=0; i < nlights; i++) {
result = result + "   colDiff = vec4(vCol.rgb * diffuse" + i + ", vCol.a);\n"+
"   lightdir = lightDir" + i + ";\n"+
"   if (!viewpoint" + i +")\n"+
"     lightdir = (mvMatrix * vec4(lightdir, 1.)).xyz;\n"+
"   if (!finite" + i + ") {\n"+
"     halfVec = normalize(lightdir + eye);\n"+
"   } else {\n"+
"     lightdir = normalize(lightdir - vPosition.xyz);\n"+
"     halfVec = normalize(lightdir + eye);\n"+
"   }\n"+
"    col = ambient" + i + ";\n"+
"   nDotL = dot(n, lightdir);\n"+
"   col = col + max(nDotL, 0.) * colDiff.rgb;\n"+
"   col = col + pow(max(dot(halfVec, n), 0.), shininess) * specular" + i + ";\n"+
"   lighteffect = lighteffect + vec4(col, colDiff.a);\n";
}
} else {
result = result +   "   vec4 colDiff = vCol;\n"+
"    vec4 lighteffect = colDiff;\n";
}
if (type === "text")
result = result +   "    vec4 textureColor = lighteffect*texture2D(uSampler, vTexcoord);\n";
if (has_texture) {
result = result + {
rgb:            "   vec4 textureColor = lighteffect*vec4(texture2D(uSampler, vTexcoord).rgb, 1.);\n",
rgba:           "   vec4 textureColor = lighteffect*texture2D(uSampler, vTexcoord);\n",
alpha:          "   vec4 textureColor = texture2D(uSampler, vTexcoord);\n"+
"   float luminance = dot(vec3(1.,1.,1.), textureColor.rgb)/3.;\n"+
"   textureColor =  vec4(lighteffect.rgb, lighteffect.a*luminance);\n",
luminance:      "   vec4 textureColor = vec4(lighteffect.rgb*dot(texture2D(uSampler, vTexcoord).rgb, vec3(1.,1.,1.))/3., lighteffect.a);\n",
"luminance.alpha":"    vec4 textureColor = texture2D(uSampler, vTexcoord);\n"+
"   float luminance = dot(vec3(1.,1.,1.),textureColor.rgb)/3.;\n"+
"   textureColor = vec4(lighteffect.rgb*luminance, lighteffect.a*textureColor.a);\n"
}[texture_format]+
"   gl_FragColor = textureColor;\n";
} else if (type === "text") {
result = result +   "    if (textureColor.a < 0.1)\n"+
"     discard;\n"+
"   else\n"+
"     gl_FragColor = textureColor;\n";
} else
result = result +   "   gl_FragColor = lighteffect;\n";
//if (fat_lines)
//  result = result +   "   gl_FragColor = vec4(0.0, abs(point.x), abs(point.y), 1.0);"
result = result + "  }\n";
// console.log(result);
return result;
};
/**
* Call gl functions to create and compile shader
* @returns {Object}
* @param { number } shaderType - gl code for shader type
* @param { string } code - code for the shader
*/
rglwidgetClass.prototype.getShader = function(shaderType, code) {
var gl = this.gl, shader;
shader = gl.createShader(shaderType);
gl.shaderSource(shader, code);
gl.compileShader(shader);
if (!gl.getShaderParameter(shader, gl.COMPILE_STATUS) && !gl.isContextLost())
alert(gl.getShaderInfoLog(shader));
return shader;
};
/**
* Handle a texture after its image has been loaded
* @param { Object } texture - the gl texture object
* @param { Object } textureCanvas - the canvas holding the image
*/
rglwidgetClass.prototype.handleLoadedTexture = function(texture, textureCanvas) {
var gl = this.gl || this.initGL();
gl.pixelStorei(gl.UNPACK_FLIP_Y_WEBGL, true);
gl.bindTexture(gl.TEXTURE_2D, texture);
gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, gl.RGBA, gl.UNSIGNED_BYTE, textureCanvas);
gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR_MIPMAP_NEAREST);
gl.generateMipmap(gl.TEXTURE_2D);
gl.bindTexture(gl.TEXTURE_2D, null);
};
/**
* Get maximum dimension of texture in current browser.
* @returns {number}
*/
rglwidgetClass.prototype.getMaxTexSize = function() {
var gl = this.gl || this.initGL();  
return Math.min(4096, gl.getParameter(gl.MAX_TEXTURE_SIZE));
};
/**
* Load an image to a texture
* @param { string } uri - The image location
* @param { Object } texture - the gl texture object
*/
rglwidgetClass.prototype.loadImageToTexture = function(uri, texture) {
var canvas = this.textureCanvas,
ctx = canvas.getContext("2d"),
image = new Image(),
self = this;
image.onload = function() {
var w = image.width,
h = image.height,
canvasX = self.getPowerOfTwo(w),
canvasY = self.getPowerOfTwo(h),
gl = self.gl || self.initGL(),
maxTexSize = self.getMaxTexSize();
while (canvasX > 1 && canvasY > 1 && (canvasX > maxTexSize || canvasY > maxTexSize)) {
canvasX /= 2;
canvasY /= 2;
}
canvas.width = canvasX;
canvas.height = canvasY;
ctx.imageSmoothingEnabled = true;
ctx.drawImage(image, 0, 0, canvasX, canvasY);
self.handleLoadedTexture(texture, canvas);
self.drawScene();
};
image.src = uri;
};
/**
* Draw text to the texture canvas
* @returns { Object } object with text measurements
* @param { string } text - the text
* @param { number } cex - expansion
* @param { string } family - font family
* @param { number } font - font number
*/
rglwidgetClass.prototype.drawTextToCanvas = function(text, cex, family, font) {
var canvasX, canvasY,
textY,
scaling = 20,
textColour = "white",
backgroundColour = "rgba(0,0,0,0)",
canvas = this.textureCanvas,
ctx = canvas.getContext("2d"),
i, textHeight = 0, textHeights = [], width, widths = [], 
offsetx, offsety = 0, line, lines = [], offsetsx = [],
offsetsy = [], lineoffsetsy = [], fontStrings = [],
maxTexSize = this.getMaxTexSize(),
getFontString = function(i) {
textHeights[i] = scaling*cex[i];
var fontString = textHeights[i] + "px",
family0 = family[i],
font0 = font[i];
if (family0 === "sans")
family0 = "sans-serif";
else if (family0 === "mono")
family0 = "monospace";
fontString = fontString + " " + family0;
if (font0 === 2 || font0 === 4)
fontString = "bold " + fontString;
if (font0 === 3 || font0 === 4)
fontString = "italic " + fontString;
return fontString;
};
cex = this.repeatToLen(cex, text.length);
family = this.repeatToLen(family, text.length);
font = this.repeatToLen(font, text.length);
canvasX = 1;
line = -1;
offsetx = maxTexSize;
for (i = 0; i < text.length; i++)  {
ctx.font = fontStrings[i] = getFontString(i);
width = widths[i] = ctx.measureText(text[i]).width;
if (offsetx + width > maxTexSize) {
line += 1;
offsety = lineoffsetsy[line] = offsety + 2*textHeight;
if (offsety > maxTexSize)
console.error("Too many strings for texture.");
textHeight = 0;
offsetx = 0;
}
textHeight = Math.max(textHeight, textHeights[i]);
offsetsx[i] = offsetx;
offsetx += width;
canvasX = Math.max(canvasX, offsetx);
lines[i] = line;
}
offsety = lineoffsetsy[line] = offsety + 2*textHeight;
for (i = 0; i < text.length; i++) {
offsetsy[i] = lineoffsetsy[lines[i]];
}
canvasX = this.getPowerOfTwo(canvasX);
canvasY = this.getPowerOfTwo(offsety);
canvas.width = canvasX;
canvas.height = canvasY;
ctx.fillStyle = backgroundColour;
ctx.fillRect(0, 0, ctx.canvas.width, ctx.canvas.height);
ctx.textBaseline = "alphabetic";
for(i = 0; i < text.length; i++) {
ctx.font = fontStrings[i];
ctx.fillStyle = textColour;
ctx.textAlign = "left";
ctx.fillText(text[i], offsetsx[i],  offsetsy[i]);
}
return {canvasX:canvasX, canvasY:canvasY,
widths:widths, textHeights:textHeights,
offsetsx:offsetsx, offsetsy:offsetsy};
};
/**
* Set the gl viewport and scissor test
* @param { number } id - id of subscene
*/
rglwidgetClass.prototype.setViewport = function(id) {
var gl = this.gl || this.initGL(),
vp = this.getObj(id).par3d.viewport,
x = vp.x*this.canvas.width,
y = vp.y*this.canvas.height,
width = vp.width*this.canvas.width,
height = vp.height*this.canvas.height;
this.vp = {x:x, y:y, width:width, height:height};
gl.viewport(x, y, width, height);
gl.scissor(x, y, width, height);
gl.enable(gl.SCISSOR_TEST);
};
/**
* Set the projection matrix for a subscene
* @param { number } id - id of subscene
*/
rglwidgetClass.prototype.setprMatrix = function(id) {
var subscene = this.getObj(id),
embedding = subscene.embeddings.projection;
if (embedding === "replace")
this.prMatrix.makeIdentity();
else
this.setprMatrix(subscene.parent);
if (embedding === "inherit")
return;
// This is based on the Frustum::enclose code from geom.cpp
var bbox = subscene.par3d.bbox,
scale = subscene.par3d.scale,
ranges = [(bbox[1]-bbox[0])*scale[0]/2,
(bbox[3]-bbox[2])*scale[1]/2,
(bbox[5]-bbox[4])*scale[2]/2],
radius = Math.sqrt(this.sumsq(ranges))*1.1; // A bit bigger to handle labels
if (radius <= 0) radius = 1;
var observer = subscene.par3d.observer,
distance = observer[2],
FOV = subscene.par3d.FOV, ortho = FOV === 0,
t = ortho ? 1 : Math.tan(FOV*Math.PI/360),
near = distance - radius,
far = distance + radius,
hlen,
aspect = this.vp.width/this.vp.height,
z = subscene.par3d.zoom,
userProjection = subscene.par3d.userProjection;
if (far < 0.0)
far = 1.0;
if (near < far/100.0)
near = far/100.0;
hlen = t*near;
if (ortho) {
if (aspect > 1)
this.prMatrix.ortho(-hlen*aspect*z, hlen*aspect*z,
-hlen*z, hlen*z, near, far);
else
this.prMatrix.ortho(-hlen*z, hlen*z,
-hlen*z/aspect, hlen*z/aspect,
near, far);
} else {
if (aspect > 1)
this.prMatrix.frustum(-hlen*aspect*z, hlen*aspect*z,
-hlen*z, hlen*z, near, far);
else
this.prMatrix.frustum(-hlen*z, hlen*z,
-hlen*z/aspect, hlen*z/aspect,
near, far);
}
this.prMatrix.multRight(userProjection);
};
/**
* Set the model-view matrix for a subscene
* @param { number } id - id of the subscene
*/
rglwidgetClass.prototype.setmvMatrix = function(id) {
var observer = this.getObj(id).par3d.observer;
this.mvMatrix.makeIdentity();
this.setmodelMatrix(id);
this.mvMatrix.translate(-observer[0], -observer[1], -observer[2]);
};
/**
* Set the model matrix for a subscene
* @param { number } id - id of the subscene
*/
rglwidgetClass.prototype.setmodelMatrix = function(id) {
var subscene = this.getObj(id),
embedding = subscene.embeddings.model;
if (embedding !== "inherit") {
var scale = subscene.par3d.scale,
bbox = subscene.par3d.bbox,
center = [(bbox[0]+bbox[1])/2,
(bbox[2]+bbox[3])/2,
(bbox[4]+bbox[5])/2];
this.mvMatrix.translate(-center[0], -center[1], -center[2]);
this.mvMatrix.scale(scale[0], scale[1], scale[2]);
this.mvMatrix.multRight( subscene.par3d.userMatrix );
}
if (embedding !== "replace")
this.setmodelMatrix(subscene.parent);
};
/**
* Set the normals matrix for a subscene
* @param { number } subsceneid - id of the subscene
*/
rglwidgetClass.prototype.setnormMatrix = function(subsceneid) {
var self = this,
recurse = function(id) {
var sub = self.getObj(id),
embedding = sub.embeddings.model;
if (embedding !== "inherit") {
var scale = sub.par3d.scale;
self.normMatrix.scale(1/scale[0], 1/scale[1], 1/scale[2]);
self.normMatrix.multRight(sub.par3d.userMatrix);
}
if (embedding !== "replace")
recurse(sub.parent);
};
self.normMatrix.makeIdentity();
recurse(subsceneid);
};
/**
* Set the combined projection-model-view matrix
*/
rglwidgetClass.prototype.setprmvMatrix = function() {
this.prmvMatrix = new CanvasMatrix4( this.mvMatrix );
this.prmvMatrix.multRight( this.prMatrix );
};
/**
* Count clipping planes in a scene
* @returns {number}
*/
rglwidgetClass.prototype.countClipplanes = function() {
return this.countObjs("clipplanes");
};
/**
* Count lights in a scene
* @returns { number }
*/
rglwidgetClass.prototype.countLights = function() {
return this.countObjs("light");
};
/**
* Count objects of specific type in a scene
* @returns { number }
* @param { string } type - Type of object to count
*/
rglwidgetClass.prototype.countObjs = function(type) {
var self = this,
bound = 0;
Object.keys(this.scene.objects).forEach(
function(key) {
if (self.getObj(parseInt(key, 10)).type === type)
bound = bound + 1;
});
return bound;
};
/**
* Initialize a subscene
* @param { number } id - id of subscene.
*/
rglwidgetClass.prototype.initSubscene = function(id) {
var sub = this.getObj(id),
i, obj;
if (sub.type !== "subscene")
return;
sub.par3d.userMatrix = this.toCanvasMatrix4(sub.par3d.userMatrix);
sub.par3d.userProjection = this.toCanvasMatrix4(sub.par3d.userProjection);
sub.par3d.userProjection.transpose();
sub.par3d.listeners = [].concat(sub.par3d.listeners);
sub.backgroundId = undefined;
sub.subscenes = [];
sub.clipplanes = [];
sub.transparent = [];
sub.opaque = [];
sub.lights = [];
for (i=0; i < sub.objects.length; i++) {
obj = this.getObj(sub.objects[i]);
if (typeof obj === "undefined") {
sub.objects.splice(i, 1);
i--;
} else if (obj.type === "background")
sub.backgroundId = obj.id;
else
sub[this.whichList(obj.id)].push(obj.id);
}
};
/**
* Copy object
* @param { number } id - id of object to copy
* @param { string } reuse - Document id of scene to reuse
*/
rglwidgetClass.prototype.copyObj = function(id, reuse) {
var obj = this.getObj(id),
prev = document.getElementById(reuse);
if (prev !== null) {
prev = prev.rglinstance;
var
prevobj = prev.getObj(id),
fields = ["flags", "type",
"colors", "vertices", "centers",
"normals", "offsets",
"texts", "cex", "family", "font", "adj",
"material",
"radii",
"texcoords",
"userMatrix", "ids",
"dim",
"par3d", "userMatrix",
"viewpoint", "finite",
"pos"],
i;
for (i = 0; i < fields.length; i++) {
if (typeof prevobj[fields[i]] !== "undefined")
obj[fields[i]] = prevobj[fields[i]];
}
} else
console.warn("copyObj failed");
};
/**
* Update the triangles used to display a plane
* @param { number } id - id of the plane
* @param { Object } bbox - bounding box in which to display the plane
*/
rglwidgetClass.prototype.planeUpdateTriangles = function(id, bbox) {
var perms = [[0,0,1], [1,2,2], [2,1,0]],
x, xrow, elem, A, d, nhits, i, j, k, u, v, w, intersect, which, v0, v2, vx, reverse,
face1 = [], face2 = [], normals = [],
obj = this.getObj(id),
nPlanes = obj.normals.length;
obj.bbox = bbox;
obj.vertices = [];
obj.initialized = false;
for (elem = 0; elem < nPlanes; elem++) {
//    Vertex Av = normal.getRecycled(elem);
x = [];
A = obj.normals[elem];
d = obj.offsets[elem][0];
nhits = 0;
for (i=0; i<3; i++)
for (j=0; j<2; j++)
for (k=0; k<2; k++) {
u = perms[0][i];
v = perms[1][i];
w = perms[2][i];
if (A[w] !== 0.0) {
intersect = -(d + A[u]*bbox[j+2*u] + A[v]*bbox[k+2*v])/A[w];
if (bbox[2*w] < intersect && intersect < bbox[1+2*w]) {
xrow = [];
xrow[u] = bbox[j+2*u];
xrow[v] = bbox[k+2*v];
xrow[w] = intersect;
x.push(xrow);
face1[nhits] = j + 2*u;
face2[nhits] = k + 2*v;
nhits++;
}
}
}
if (nhits > 3) {
/* Re-order the intersections so the triangles work */
for (i=0; i<nhits-2; i++) {
which = 0; /* initialize to suppress warning */
for (j=i+1; j<nhits; j++) {
if (face1[i] == face1[j] || face1[i] == face2[j] ||
face2[i] == face1[j] || face2[i] == face2[j] ) {
which = j;
break;
}
}
if (which > i+1) {
this.swap(x, i+1, which);
this.swap(face1, i+1, which);
this.swap(face2, i+1, which);
}
}
}
if (nhits >= 3) {
/* Put in order so that the normal points out the FRONT of the faces */
v0 = [x[0][0] - x[1][0] , x[0][1] - x[1][1], x[0][2] - x[1][2]];
v2 = [x[2][0] - x[1][0] , x[2][1] - x[1][1], x[2][2] - x[1][2]];
/* cross-product */
vx = this.xprod(v0, v2);
reverse = this.dotprod(vx, A) > 0;
for (i=0; i<nhits-2; i++) {
obj.vertices.push(x[0]);
normals.push(A);
for (j=1; j<3; j++) {
obj.vertices.push(x[i + (reverse ? 3-j : j)]);
normals.push(A);
}
}
}
}
obj.pnormals = normals;
};
rglwidgetClass.prototype.getAdj = function (pos, offset, text) {
switch(pos) {
case 1: return [0.5, 1 + offset];
case 2: return [1 + offset/text.length, 0.5];
case 3: return [0.5, -offset];
case 4: return [-offset/text.length, 0.5];
}
}
/**
* Initialize object for display
* @param { number } id - id of object to initialize
*/
rglwidgetClass.prototype.initObj = function(id) {
var obj = this.getObj(id),
flags = obj.flags,
type = obj.type,
is_lit = flags & this.f_is_lit,
is_lines = flags & this.f_is_lines,
fat_lines = flags & this.f_fat_lines,
has_texture = flags & this.f_has_texture,
fixed_quads = flags & this.f_fixed_quads,
is_transparent = obj.is_transparent,
depth_sort = flags & this.f_depth_sort,
sprites_3d = flags & this.f_sprites_3d,
sprite_3d = flags & this.f_sprite_3d,
fixed_size = flags & this.f_fixed_size,
is_twosided = (flags & this.f_is_twosided) > 0,
is_brush = flags & this.f_is_brush,
gl = this.gl || this.initGL(),
polygon_offset,
texinfo, drawtype, nclipplanes, f, nrows, oldrows,
i,j,v,v1,v2, mat, uri, matobj, pass, passes, pmode,
dim, nx, nz, attr;
if (typeof id !== "number") {
this.alertOnce("initObj id is "+typeof id);
}
obj.initialized = true;
if (type === "bboxdeco" || type === "subscene")
return;
if (type === "light") {
obj.ambient = new Float32Array(obj.colors[0].slice(0,3));
obj.diffuse = new Float32Array(obj.colors[1].slice(0,3));
obj.specular = new Float32Array(obj.colors[2].slice(0,3));
obj.lightDir = new Float32Array(obj.vertices[0]);
return;
}
if (type === "clipplanes") {
obj.vClipplane = this.flatten(this.cbind(obj.normals, obj.offsets));
return;
}
if (type === "background" && typeof obj.ids !== "undefined") {
obj.quad = this.flatten([].concat(obj.ids));
return;
}
polygon_offset = this.getMaterial(id, "polygon_offset");
if (polygon_offset[0] != 0 || polygon_offset[1] != 0)
obj.polygon_offset = polygon_offset;
if (is_transparent) {
depth_sort = ["triangles", "quads", "surface",
"spheres", "sprites", "text"].indexOf(type) >= 0;
}
if (is_brush)
this.initSelection(id);
if (typeof obj.vertices === "undefined")
obj.vertices = [];
v = obj.vertices;
obj.vertexCount = v.length;
if (!obj.vertexCount) return;
if (is_twosided) {
if (typeof obj.userAttributes === "undefined")
obj.userAttributes = {};
v1 = Array(v.length);
v2 = Array(v.length);
if (obj.type == "triangles" || obj.type == "quads") {
if (obj.type == "triangles")
nrow = 3;
else
nrow = 4;
for (i=0; i<Math.floor(v.length/nrow); i++)
for (j=0; j<nrow; j++) {
v1[nrow*i + j] = v[nrow*i + ((j+1) % nrow)];
v2[nrow*i + j] = v[nrow*i + ((j+2) % nrow)];
}
} else if (obj.type == "surface") {
dim = obj.dim[0];
nx = dim[0];
nz = dim[1];
for (j=0; j<nx; j++) {
for (i=0; i<nz; i++) {
if (i+1 < nz && j+1 < nx) {
v2[j + nx*i] = v[j + nx*(i+1)];
v1[j + nx*i] = v[j+1 + nx*(i+1)];
} else if (i+1 < nz) {
v2[j + nx*i] = v[j-1 + nx*i];
v1[j + nx*i] = v[j + nx*(i+1)];
} else {
v2[j + nx*i] = v[j + nx*(i-1)];
v1[j + nx*i] = v[j-1 + nx*(i-1)];
}
}
}
}
obj.userAttributes.aPos1 = v1;
obj.userAttributes.aPos2 = v2;
}
if (!sprites_3d) {
if (gl.isContextLost()) return;
obj.prog = gl.createProgram();
gl.attachShader(obj.prog, this.getShader( gl.VERTEX_SHADER,
this.getVertexShader(id) ));
gl.attachShader(obj.prog, this.getShader( gl.FRAGMENT_SHADER,
this.getFragmentShader(id) ));
//  Force aPos to location 0, aCol to location 1
gl.bindAttribLocation(obj.prog, 0, "aPos");
gl.bindAttribLocation(obj.prog, 1, "aCol");
gl.linkProgram(obj.prog);
var linked = gl.getProgramParameter(obj.prog, gl.LINK_STATUS);
if (!linked) {
// An error occurred while linking
var lastError = gl.getProgramInfoLog(obj.prog);
console.warn("Error in program linking:" + lastError);
gl.deleteProgram(obj.prog);
return;
}
}
if (type === "text") {
texinfo = this.drawTextToCanvas(obj.texts,
this.flatten(obj.cex),
this.flatten(obj.family),
this.flatten(obj.family));
}
if (fixed_quads && !sprites_3d) {
obj.ofsLoc = gl.getAttribLocation(obj.prog, "aOfs");
}
if (sprite_3d) {
obj.origLoc = gl.getUniformLocation(obj.prog, "uOrig");
obj.sizeLoc = gl.getUniformLocation(obj.prog, "uSize");
obj.usermatLoc = gl.getUniformLocation(obj.prog, "usermat");
}
if (has_texture || type == "text") {
if (!obj.texture)
obj.texture = gl.createTexture();
obj.texLoc = gl.getAttribLocation(obj.prog, "aTexcoord");
obj.sampler = gl.getUniformLocation(obj.prog, "uSampler");
}
if (has_texture) {
mat = obj.material;
if (typeof mat.uri !== "undefined")
uri = mat.uri;
else if (typeof mat.uriElementId === "undefined") {
matobj = this.getObj(mat.uriId);
if (typeof matobj !== "undefined") {
uri = matobj.material.uri;
} else {
uri = "";
}
} else
uri = document.getElementById(mat.uriElementId).rglinstance.getObj(mat.uriId).material.uri;
this.loadImageToTexture(uri, obj.texture);
}
if (type === "text") {
this.handleLoadedTexture(obj.texture, this.textureCanvas);
}
var stride = 3, nc, cofs, nofs, radofs, oofs, tofs, vnew, fnew,
nextofs = -1, pointofs = -1, alias, colors, key, selection, filter, adj, pos, offset;
obj.alias = undefined;
colors = obj.colors;
j = this.scene.crosstalk.id.indexOf(id);
if (j >= 0) {
key = this.scene.crosstalk.key[j];
options = this.scene.crosstalk.options[j];
colors = colors.slice(0); 
for (i = 0; i < v.length; i++)
colors[i] = obj.colors[i % obj.colors.length].slice(0);
if ( (selection = this.scene.crosstalk.selection) &&
(selection.length || !options.selectedIgnoreNone) )
for (i = 0; i < v.length; i++) {
if (!selection.includes(key[i])) {
if (options.deselectedColor)
colors[i] = options.deselectedColor.slice(0);
colors[i][3] = colors[i][3]*options.deselectedFade;   /* default: mostly transparent if not selected */
} else if (options.selectedColor)
colors[i] = options.selectedColor.slice(0);
}
if ( (filter = this.scene.crosstalk.filter) )
for (i = 0; i < v.length; i++) 
if (!filter.includes(key[i])) {
if (options.filteredColor)
colors[i] = options.filteredColor.slice(0);
colors[i][3] = colors[i][3]*options.filteredFade;   /* default: completely hidden if filtered */
}
}  
nc = obj.colorCount = colors.length;
if (nc > 1) {
cofs = stride;
stride = stride + 4;
v = this.cbind(v, colors);
} else {
cofs = -1;
obj.onecolor = this.flatten(colors);
}
if (typeof obj.normals !== "undefined") {
nofs = stride;
stride = stride + 3;
v = this.cbind(v, typeof obj.pnormals !== "undefined" ? obj.pnormals : obj.normals);
} else
nofs = -1;
if (typeof obj.radii !== "undefined") {
radofs = stride;
stride = stride + 1;
// FIXME:  always concat the radii?
if (obj.radii.length === v.length) {
v = this.cbind(v, obj.radii);
} else if (obj.radii.length === 1) {
v = v.map(function(row, i, arr) { return row.concat(obj.radii[0]);});
}
} else
radofs = -1;
// Add default indices
f = Array(v.length);
for (i = 0; i < v.length; i++)
f[i] = i;
obj.f = [f,f];
if (type == "sprites" && !sprites_3d) {
tofs = stride;
stride += 2;
oofs = stride;
stride += 2;
vnew = new Array(4*v.length);
fnew = new Array(4*v.length);
alias = new Array(v.length);
var rescale = fixed_size ? 72 : 1,
size = obj.radii, s = rescale*size[0]/2;
last = v.length;
f = obj.f[0];
for (i=0; i < v.length; i++) {
if (size.length > 1)
s = rescale*size[i]/2;
vnew[i]  = v[i].concat([0,0,-s,-s]);
fnew[4*i] = f[i];
vnew[last]= v[i].concat([1,0, s,-s]);
fnew[4*i+1] = last++;
vnew[last]= v[i].concat([1,1, s, s]);
fnew[4*i+2] = last++;
vnew[last]= v[i].concat([0,1,-s, s]);
fnew[4*i+3] = last++;
alias[i] = [last-3, last-2, last-1];
}
v = vnew;
obj.vertexCount = v.length;
obj.f = [fnew, fnew];
} else if (type === "text") {
tofs = stride;
stride += 2;
oofs = stride;
stride += 2;
vnew = new Array(4*v.length);
f = obj.f[0];
fnew = new Array(4*f.length);
alias = new Array(v.length);
last = v.length;
adj = this.flatten(obj.adj);
if (typeof obj.pos !== "undefined") {
pos = this.flatten(obj.pos);
offset = adj[0];
}
for (i=0; i < v.length; i++) {
if (typeof pos !== "undefined")
adj = this.getAdj(pos[i % pos.length], offset, obj.texts[i]);
vnew[i]  = v[i].concat([0,-0.5]).concat(adj);
fnew[4*i] = f[i];
vnew[last] = v[i].concat([1,-0.5]).concat(adj);
fnew[4*i+1] = last++;
vnew[last] = v[i].concat([1, 1.5]).concat(adj);
fnew[4*i+2] = last++;
vnew[last] = v[i].concat([0, 1.5]).concat(adj);
fnew[4*i+3] = last++;
alias[i] = [last-3, last-2, last-1];
for (j=0; j < 4; j++) {
v1 = vnew[fnew[4*i+j]];
v1[tofs+2] = 2*(v1[tofs]-v1[tofs+2])*texinfo.widths[i];
v1[tofs+3] = 2*(v1[tofs+1]-v1[tofs+3])*texinfo.textHeights[i];
v1[tofs] = (texinfo.offsetsx[i] + v1[tofs]*texinfo.widths[i])/texinfo.canvasX;
v1[tofs+1] = 1.0-(texinfo.offsetsy[i] -
v1[tofs+1]*texinfo.textHeights[i])/texinfo.canvasY;
vnew[fnew[4*i+j]] = v1;
}
}
v = vnew;
obj.vertexCount = v.length;
obj.f = [fnew, fnew];
} else if (typeof obj.texcoords !== "undefined") {
tofs = stride;
stride += 2;
oofs = -1;
v = this.cbind(v, obj.texcoords);
} else {
tofs = -1;
oofs = -1;
}
obj.alias = alias;
if (typeof obj.userAttributes !== "undefined") {
obj.userAttribOffsets = {};
obj.userAttribLocations = {};
obj.userAttribSizes = {};
for (attr in obj.userAttributes) {
obj.userAttribLocations[attr] = gl.getAttribLocation(obj.prog, attr);
if (obj.userAttribLocations[attr] >= 0) { // Attribute may not have been used
obj.userAttribOffsets[attr] = stride;
v = this.cbind(v, obj.userAttributes[attr]);
stride = v[0].length;
obj.userAttribSizes[attr] = stride - obj.userAttribOffsets[attr];
}
}
}
if (typeof obj.userUniforms !== "undefined") {
obj.userUniformLocations = {};
for (attr in obj.userUniforms)
obj.userUniformLocations[attr] = gl.getUniformLocation(obj.prog, attr);
}
if (sprites_3d) {
obj.userMatrix = new CanvasMatrix4(obj.userMatrix);
obj.objects = this.flatten([].concat(obj.ids));
is_lit = false;
for (i=0; i < obj.objects.length; i++)
this.initObj(obj.objects[i]);
}
if (is_lit && !fixed_quads) {
obj.normLoc = gl.getAttribLocation(obj.prog, "aNorm");
}
nclipplanes = this.countClipplanes();
if (nclipplanes && !sprites_3d) {
obj.clipLoc = [];
for (i=0; i < nclipplanes; i++)
obj.clipLoc[i] = gl.getUniformLocation(obj.prog,"vClipplane" + i);
}
if (is_lit) {
obj.emissionLoc = gl.getUniformLocation(obj.prog, "emission");
obj.emission = new Float32Array(this.stringToRgb(this.getMaterial(id, "emission")));
obj.shininessLoc = gl.getUniformLocation(obj.prog, "shininess");
obj.shininess = this.getMaterial(id, "shininess");
obj.nlights = this.countLights();
obj.ambientLoc = [];
obj.ambient = new Float32Array(this.stringToRgb(this.getMaterial(id, "ambient")));
obj.specularLoc = [];
obj.specular = new Float32Array(this.stringToRgb(this.getMaterial(id, "specular")));
obj.diffuseLoc = [];
obj.lightDirLoc = [];
obj.viewpointLoc = [];
obj.finiteLoc = [];
for (i=0; i < obj.nlights; i++) {
obj.ambientLoc[i] = gl.getUniformLocation(obj.prog, "ambient" + i);
obj.specularLoc[i] = gl.getUniformLocation(obj.prog, "specular" + i);
obj.diffuseLoc[i] = gl.getUniformLocation(obj.prog, "diffuse" + i);
obj.lightDirLoc[i] = gl.getUniformLocation(obj.prog, "lightDir" + i);
obj.viewpointLoc[i] = gl.getUniformLocation(obj.prog, "viewpoint" + i);
obj.finiteLoc[i] = gl.getUniformLocation(obj.prog, "finite" + i);
}
}
obj.passes = is_twosided + 1;
obj.pmode = new Array(obj.passes);
for (pass = 0; pass < obj.passes; pass++) {
if (type === "triangles" || type === "quads" || type === "surface")
pmode = this.getMaterial(id, (pass === 0) ? "front" : "back");
else pmode = "filled";
obj.pmode[pass] = pmode;
}
obj.f.length = obj.passes;
for (pass = 0; pass < obj.passes; pass++) {
f = fnew = obj.f[pass];
pmode = obj.pmode[pass];
if (pmode === "culled")
f = [];
else if (pmode === "points") {
// stay with default
} else if ((type === "quads" || type === "text" ||
type === "sprites") && !sprites_3d) {
nrows = Math.floor(obj.vertexCount/4);
if (pmode === "filled") {
fnew = Array(6*nrows);
for (i=0; i < nrows; i++) {
fnew[6*i] = f[4*i];
fnew[6*i+1] = f[4*i + 1];
fnew[6*i+2] = f[4*i + 2];
fnew[6*i+3] = f[4*i];
fnew[6*i+4] = f[4*i + 2];
fnew[6*i+5] = f[4*i + 3];
}
} else {
fnew = Array(8*nrows);
for (i=0; i < nrows; i++) {
fnew[8*i] = f[4*i];
fnew[8*i+1] = f[4*i + 1];
fnew[8*i+2] = f[4*i + 1];
fnew[8*i+3] = f[4*i + 2];
fnew[8*i+4] = f[4*i + 2];
fnew[8*i+5] = f[4*i + 3];
fnew[8*i+6] = f[4*i + 3];
fnew[8*i+7] = f[4*i];
}
}
} else if (type === "triangles") {
nrows = Math.floor(obj.vertexCount/3);
if (pmode === "filled") {
fnew = Array(3*nrows);
for (i=0; i < fnew.length; i++) {
fnew[i] = f[i];
}
} else if (pmode === "lines") {
fnew = Array(6*nrows);
for (i=0; i < nrows; i++) {
fnew[6*i] = f[3*i];
fnew[6*i + 1] = f[3*i + 1];
fnew[6*i + 2] = f[3*i + 1];
fnew[6*i + 3] = f[3*i + 2];
fnew[6*i + 4] = f[3*i + 2];
fnew[6*i + 5] = f[3*i];
}
}
} else if (type === "spheres") {
// default
} else if (type === "surface") {
dim = obj.dim[0];
nx = dim[0];
nz = dim[1];
if (pmode === "filled") {
fnew = [];
for (j=0; j<nx-1; j++) {
for (i=0; i<nz-1; i++) {
fnew.push(f[j + nx*i],
f[j + nx*(i+1)],
f[j + 1 + nx*(i+1)],
f[j + nx*i],
f[j + 1 + nx*(i+1)],
f[j + 1 + nx*i]);
}
}
} else if (pmode === "lines") {
fnew = [];
for (j=0; j<nx; j++) {
for (i=0; i<nz; i++) {
if (i+1 < nz)
fnew.push(f[j + nx*i],
f[j + nx*(i+1)]);
if (j+1 < nx)
fnew.push(f[j + nx*i],
f[j+1 + nx*i]);
}
}
}
}
obj.f[pass] = fnew;
if (depth_sort) {
drawtype = "DYNAMIC_DRAW";
} else {
drawtype = "STATIC_DRAW";
}
}
if (fat_lines) {
alias = undefined;
obj.nextLoc = gl.getAttribLocation(obj.prog, "aNext");
obj.pointLoc = gl.getAttribLocation(obj.prog, "aPoint");
obj.aspectLoc = gl.getUniformLocation(obj.prog, "uAspect");
obj.lwdLoc = gl.getUniformLocation(obj.prog, "uLwd");
// Expand vertices to turn each segment into a pair of triangles
for (pass = 0; pass < obj.passes; pass++) {
f = obj.f[pass];    
oldrows = f.length;
if (obj.pmode[pass] === "lines") 
break;
}
if (type === "linestrip") 
nrows = 4*(oldrows - 1); 
else
nrows = 2*oldrows;
vnew = new Array(nrows);
fnew = new Array(1.5*nrows);
var fnext = new Array(nrows),
fpt = new Array(nrows), 
pt, start, gap = type === "linestrip" ? 3 : 1;
// We're going to turn each pair of vertices into 4 new ones, with the "next" and "pt" attributes
// added.
// We do this by copying the originals in the first pass, adding the new attributes, then in a 
// second pass add new vertices at the end.
for (i = 0; i < v.length; i++) {
vnew[i] = v[i].concat([0,0,0,0,0]); 
}
nextofs = stride;
pointofs = stride + 3;
stride = stride + 5;
// Now add the extras
last = v.length - 1;
ind = 0;
alias = new Array(f.length);
for (i = 0; i < f.length; i++)
alias[i] = [];
for (i = 0; i < f.length - 1; i++) {
if (type !== "linestrip" && i % 2 == 1)
continue;
k = ++last;
vnew[k] = vnew[f[i]].slice();
for (j=0; j<3; j++)
vnew[k][nextofs + j] = vnew[f[i+1]][j];
vnew[k][pointofs] = -1;
vnew[k][pointofs+1] = -1;
fnew[ind] = k;
last++;
vnew[last] = vnew[k].slice();
vnew[last][pointofs] = 1;
fnew[ind+1] = last;
alias[f[i]].push(last-1, last);
last++;
k = last;
vnew[k] = vnew[f[i+1]].slice();
for (j=0; j<3; j++)
vnew[k][nextofs + j] = vnew[f[i]][j];
vnew[k][pointofs] = -1;
vnew[k][pointofs+1] = 1;
fnew[ind+2] = k;
fnew[ind+3] = fnew[ind+1];
last++;
vnew[last] = vnew[k].slice();
vnew[last][pointofs] = 1;
fnew[ind+4] = last;
fnew[ind+5] = fnew[ind+2];
ind += 6;
alias[f[i+1]].push(last-1, last);
}
vnew.length = last+1;
v = vnew;
obj.vertexCount = v.length;
if (typeof alias !== "undefined" && typeof obj.alias !== "undefined") {  // Already have aliases from previous section?
var oldalias = obj.alias, newalias = Array(obj.alias.length);
for (i = 0; i < newalias.length; i++) {
newalias[i] = oldalias[i].slice();
for (j = 0; j < oldalias[i].length; j++)
Array.prototype.push.apply(newalias[i], alias[oldalias[j]]); // pushes each element 
}
obj.alias = newalias;
} else
obj.alias = alias;
for (pass = 0; pass < obj.passes; pass++)
if (type === "lines" || type === "linestrip" || obj.pmode[pass] == "lines") {
obj.f[pass] = fnew;
}
if (depth_sort) 
drawtype = "DYNAMIC_DRAW";
else
drawtype = "STATIC_DRAW";
}
for (pass = 0; pass < obj.passes; pass++) {
if (obj.vertexCount > 65535) {
if (this.index_uint) {
obj.f[pass] = new Uint32Array(obj.f[pass]);
obj.index_uint = true;
} else
this.alertOnce("Object has "+obj.vertexCount+" vertices, not supported in this browser.");
} else {
obj.f[pass] = new Uint16Array(obj.f[pass]);
obj.index_uint = false;
}
}
if (stride !== v[0].length) {
this.alertOnce("problem in stride calculation");
}
obj.vOffsets = {vofs:0, cofs:cofs, nofs:nofs, radofs:radofs, oofs:oofs, tofs:tofs,
nextofs:nextofs, pointofs:pointofs, stride:stride};
obj.values = new Float32Array(this.flatten(v));
if (type !== "spheres" && !sprites_3d) {
obj.buf = gl.createBuffer();
gl.bindBuffer(gl.ARRAY_BUFFER, obj.buf);
gl.bufferData(gl.ARRAY_BUFFER, obj.values, gl.STATIC_DRAW); //
obj.ibuf = Array(obj.passes);
obj.ibuf[0] = gl.createBuffer();
gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, obj.ibuf[0]);
gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, obj.f[0], gl[drawtype]);
if (is_twosided) {
obj.ibuf[1] = gl.createBuffer();
gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, obj.ibuf[1]);
gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, obj.f[1], gl[drawtype]);
}
}
if (!sprites_3d) {
obj.mvMatLoc = gl.getUniformLocation(obj.prog, "mvMatrix");
obj.prMatLoc = gl.getUniformLocation(obj.prog, "prMatrix");
}
if (fixed_size) {
obj.textScaleLoc = gl.getUniformLocation(obj.prog, "textScale");
}
if (is_lit && !sprites_3d) {
obj.normMatLoc = gl.getUniformLocation(obj.prog, "normMatrix");
}
if (is_twosided) {
obj.frontLoc = gl.getUniformLocation(obj.prog, "front");
}
};
/**
* Set gl depth test based on object's material
* @param { number } id - object to use
*/
rglwidgetClass.prototype.setDepthTest = function(id) {
var gl = this.gl || this.initGL(),
tests = {never: gl.NEVER,
less:  gl.LESS,
equal: gl.EQUAL,
lequal:gl.LEQUAL,
greater: gl.GREATER,
notequal: gl.NOTEQUAL,
gequal: gl.GEQUAL,
always: gl.ALWAYS},
test = tests[this.getMaterial(id, "depth_test")];
gl.depthFunc(test);
};
rglwidgetClass.prototype.mode4type = {points : "POINTS",
linestrip : "LINE_STRIP",
abclines : "LINES",
lines : "LINES",
sprites : "TRIANGLES",
planes : "TRIANGLES",
text : "TRIANGLES",
quads : "TRIANGLES",
surface : "TRIANGLES",
triangles : "TRIANGLES"};
/**
* Sort objects from back to front
* @returns { number[] }
* @param { Object } obj - object to sort
*/
rglwidgetClass.prototype.depthSort = function(obj) {
var n = obj.centers.length,
depths = new Float32Array(n),
result = new Array(n),
compare = function(i,j) { return depths[j] - depths[i]; },
z, w;
for(i=0; i<n; i++) {
z = this.prmvMatrix.m13*obj.centers[i][0] +
this.prmvMatrix.m23*obj.centers[i][1] +
this.prmvMatrix.m33*obj.centers[i][2] +
this.prmvMatrix.m43;
w = this.prmvMatrix.m14*obj.centers[i][0] +
this.prmvMatrix.m24*obj.centers[i][1] +
this.prmvMatrix.m34*obj.centers[i][2] +
this.prmvMatrix.m44;
depths[i] = z/w;
result[i] = i;
}
result.sort(compare);
return result;
};
rglwidgetClass.prototype.disableArrays = function(obj, enabled) {
var gl = this.gl || this.initGL(),
objLocs = ["normLoc", "texLoc", "ofsLoc", "pointLoc", "nextLoc"],
thisLocs = ["posLoc", "colLoc"], i, attr;
for (i = 0; i < objLocs.length; i++) 
if (enabled[objLocs[i]]) gl.disableVertexAttribArray(obj[objLocs[i]]);
for (i = 0; i < thisLocs.length; i++)
if (enabled[thisLocs[i]]) gl.disableVertexAttribArray(this[objLocs[i]]);
if (typeof obj.userAttributes !== "undefined") {
for (attr in obj.userAttribSizes) {  // Not all attributes may have been used
gl.disableVertexAttribArray( obj.userAttribLocations[attr] );
}
}
}
/**
* Draw an object in a subscene
* @param { number } id - object to draw
* @param { number } subsceneid - id of subscene
*/
rglwidgetClass.prototype.drawObj = function(id, subsceneid) {
var obj = this.getObj(id),
subscene = this.getObj(subsceneid),
flags = obj.flags,
type = obj.type,
is_lit = flags & this.f_is_lit,
has_texture = flags & this.f_has_texture,
fixed_quads = flags & this.f_fixed_quads,
is_transparent = flags & this.f_is_transparent,
depth_sort = flags & this.f_depth_sort,
sprites_3d = flags & this.f_sprites_3d,
sprite_3d = flags & this.f_sprite_3d,
is_lines = flags & this.f_is_lines,
fat_lines = flags & this.f_fat_lines,
is_points = flags & this.f_is_points,
fixed_size = flags & this.f_fixed_size,
is_twosided = (flags & this.f_is_twosided) > 0,
gl = this.gl || this.initGL(),
mat,
sphereMV, baseofs, ofs, sscale, i, count, light,
pass, mode, pmode, attr,
enabled = {};
if (typeof id !== "number") {
this.alertOnce("drawObj id is "+typeof id);
}
if (type === "planes") {
if (obj.bbox !== subscene.par3d.bbox || !obj.initialized) {
this.planeUpdateTriangles(id, subscene.par3d.bbox);
}
}
if (!obj.initialized)
this.initObj(id);
if (type === "clipplanes") {
count = obj.offsets.length;
var IMVClip = [];
for (i=0; i < count; i++) {
IMVClip[i] = this.multMV(this.invMatrix, obj.vClipplane.slice(4*i, 4*(i+1)));
}
obj.IMVClip = IMVClip;
return;
}
if (type === "light" || type === "bboxdeco" || !obj.vertexCount)
return;
if (!is_transparent &&
obj.someHidden) {
is_transparent = true;
depth_sort = ["triangles", "quads", "surface",
"spheres", "sprites", "text"].indexOf(type) >= 0;
}        
this.setDepthTest(id);
if (sprites_3d) {
var norigs = obj.vertices.length,
savenorm = new CanvasMatrix4(this.normMatrix);
this.origs = obj.vertices;
this.usermat = new Float32Array(obj.userMatrix.getAsArray());
this.radii = obj.radii;
this.normMatrix = subscene.spriteNormmat;
for (this.iOrig=0; this.iOrig < norigs; this.iOrig++) {
for (i=0; i < obj.objects.length; i++) {
this.drawObj(obj.objects[i], subsceneid);
}
}
this.normMatrix = savenorm;
return;
} else {
gl.useProgram(obj.prog);
}
if (typeof obj.polygon_offset !== "undefined") {
gl.polygonOffset(obj.polygon_offset[0],
obj.polygon_offset[1]);
gl.enable(gl.POLYGON_OFFSET_FILL);
}
if (sprite_3d) {
gl.uniform3fv(obj.origLoc, new Float32Array(this.origs[this.iOrig]));
if (this.radii.length > 1) {
gl.uniform1f(obj.sizeLoc, this.radii[this.iOrig][0]);
} else {
gl.uniform1f(obj.sizeLoc, this.radii[0][0]);
}
gl.uniformMatrix4fv(obj.usermatLoc, false, this.usermat);
}
if (type === "spheres") {
gl.bindBuffer(gl.ARRAY_BUFFER, this.sphere.buf);
} else {
gl.bindBuffer(gl.ARRAY_BUFFER, obj.buf);
}
gl.uniformMatrix4fv( obj.prMatLoc, false, new Float32Array(this.prMatrix.getAsArray()) );
gl.uniformMatrix4fv( obj.mvMatLoc, false, new Float32Array(this.mvMatrix.getAsArray()) );
var clipcheck = 0,
clipplaneids = subscene.clipplanes,
clip, j;
for (i=0; i < clipplaneids.length; i++) {
clip = this.getObj(clipplaneids[i]);
for (j=0; j < clip.offsets.length; j++) {
gl.uniform4fv(obj.clipLoc[clipcheck + j], clip.IMVClip[j]);
}
clipcheck += clip.offsets.length;
}
if (typeof obj.clipLoc !== "undefined")
for (i=clipcheck; i < obj.clipLoc.length; i++)
gl.uniform4f(obj.clipLoc[i], 0,0,0,0);
if (is_lit) {
gl.uniformMatrix4fv( obj.normMatLoc, false, new Float32Array(this.normMatrix.getAsArray()) );
gl.uniform3fv( obj.emissionLoc, obj.emission);
gl.uniform1f( obj.shininessLoc, obj.shininess);
for (i=0; i < subscene.lights.length; i++) {
light = this.getObj(subscene.lights[i]);
if (!light.initialized) this.initObj(subscene.lights[i]);
gl.uniform3fv( obj.ambientLoc[i], this.componentProduct(light.ambient, obj.ambient));
gl.uniform3fv( obj.specularLoc[i], this.componentProduct(light.specular, obj.specular));
gl.uniform3fv( obj.diffuseLoc[i], light.diffuse);
gl.uniform3fv( obj.lightDirLoc[i], light.lightDir);
gl.uniform1i( obj.viewpointLoc[i], light.viewpoint);
gl.uniform1i( obj.finiteLoc[i], light.finite);
}
for (i=subscene.lights.length; i < obj.nlights; i++) {
gl.uniform3f( obj.ambientLoc[i], 0,0,0);
gl.uniform3f( obj.specularLoc[i], 0,0,0);
gl.uniform3f( obj.diffuseLoc[i], 0,0,0);
}
}
if (fixed_size) {
gl.uniform2f( obj.textScaleLoc, 0.75/this.vp.width, 0.75/this.vp.height);
}
gl.enableVertexAttribArray( this.posLoc );
enabled.posLoc = true;
var nc = obj.colorCount;
count = obj.vertexCount;
if (type === "spheres") {
subscene = this.getObj(subsceneid);
var scale = subscene.par3d.scale,
scount = count, indices;
gl.vertexAttribPointer(this.posLoc,  3, gl.FLOAT, false, 4*this.sphere.vOffsets.stride,  0);
gl.enableVertexAttribArray(obj.normLoc );
enabled.normLoc = true;
gl.vertexAttribPointer(obj.normLoc,  3, gl.FLOAT, false, 4*this.sphere.vOffsets.stride,  0);
gl.disableVertexAttribArray( this.colLoc );
var sphereNorm = new CanvasMatrix4();
sphereNorm.scale(scale[0], scale[1], scale[2]);
sphereNorm.multRight(this.normMatrix);
gl.uniformMatrix4fv( obj.normMatLoc, false, new Float32Array(sphereNorm.getAsArray()) );
if (nc == 1) {
gl.vertexAttrib4fv( this.colLoc, new Float32Array(obj.onecolor));
}
if (has_texture) {
gl.enableVertexAttribArray( obj.texLoc );
enabled.texLoc = true;
gl.vertexAttribPointer(obj.texLoc, 2, gl.FLOAT, false, 4*this.sphere.vOffsets.stride,
4*this.sphere.vOffsets.tofs);
gl.activeTexture(gl.TEXTURE0);
gl.bindTexture(gl.TEXTURE_2D, obj.texture);
gl.uniform1i( obj.sampler, 0);
}
if (depth_sort)
indices = this.depthSort(obj);
for (i = 0; i < scount; i++) {
sphereMV = new CanvasMatrix4();
if (depth_sort) {
baseofs = indices[i]*obj.vOffsets.stride;
} else {
baseofs = i*obj.vOffsets.stride;
}
ofs = baseofs + obj.vOffsets.radofs;
sscale = obj.values[ofs];
sphereMV.scale(sscale/scale[0], sscale/scale[1], sscale/scale[2]);
sphereMV.translate(obj.values[baseofs],
obj.values[baseofs+1],
obj.values[baseofs+2]);
sphereMV.multRight(this.mvMatrix);
gl.uniformMatrix4fv( obj.mvMatLoc, false, new Float32Array(sphereMV.getAsArray()) );
if (nc > 1) {
ofs = baseofs + obj.vOffsets.cofs;
gl.vertexAttrib4f( this.colLoc, obj.values[ofs],
obj.values[ofs+1],
obj.values[ofs+2],
obj.values[ofs+3] );
}
gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, this.sphere.ibuf);
gl.drawElements(gl.TRIANGLES, this.sphere.sphereCount, gl.UNSIGNED_SHORT, 0);
}
this.disableArrays(obj, enabled);
if (typeof obj.polygon_offset !== "undefined") 
gl.disable(gl.POLYGON_OFFSET_FILL);
return;
} else {
if (obj.colorCount === 1) {
gl.disableVertexAttribArray( this.colLoc );
gl.vertexAttrib4fv( this.colLoc, new Float32Array(obj.onecolor));
} else {
gl.enableVertexAttribArray( this.colLoc );
enabled.colLoc = true;
gl.vertexAttribPointer(this.colLoc, 4, gl.FLOAT, false, 4*obj.vOffsets.stride, 4*obj.vOffsets.cofs);
}
}
if (is_lit && obj.vOffsets.nofs > 0) {
gl.enableVertexAttribArray( obj.normLoc );
enabled.normLoc = true;
gl.vertexAttribPointer(obj.normLoc, 3, gl.FLOAT, false, 4*obj.vOffsets.stride, 4*obj.vOffsets.nofs);
}
if (has_texture || type === "text") {
gl.enableVertexAttribArray( obj.texLoc );
enabled.texLoc = true;
gl.vertexAttribPointer(obj.texLoc, 2, gl.FLOAT, false, 4*obj.vOffsets.stride, 4*obj.vOffsets.tofs);
gl.activeTexture(gl.TEXTURE0);
gl.bindTexture(gl.TEXTURE_2D, obj.texture);
gl.uniform1i( obj.sampler, 0);
}
if (fixed_quads) {
gl.enableVertexAttribArray( obj.ofsLoc );
enabled.ofsLoc = true;
gl.vertexAttribPointer(obj.ofsLoc, 2, gl.FLOAT, false, 4*obj.vOffsets.stride, 4*obj.vOffsets.oofs);
}
if (typeof obj.userAttributes !== "undefined") {
for (attr in obj.userAttribSizes) {  // Not all attributes may have been used
gl.enableVertexAttribArray( obj.userAttribLocations[attr] );
gl.vertexAttribPointer( obj.userAttribLocations[attr], obj.userAttribSizes[attr],
gl.FLOAT, false, 4*obj.vOffsets.stride, 4*obj.userAttribOffsets[attr]);
}
}
if (typeof obj.userUniforms !== "undefined") {
for (attr in obj.userUniformLocations) {
var loc = obj.userUniformLocations[attr];
if (loc !== null) {
var uniform = obj.userUniforms[attr];
if (typeof uniform.length === "undefined")
gl.uniform1f(loc, uniform);
else if (typeof uniform[0].length === "undefined") {
uniform = new Float32Array(uniform);
switch(uniform.length) {
case 2: gl.uniform2fv(loc, uniform); break;
case 3: gl.uniform3fv(loc, uniform); break;
case 4: gl.uniform4fv(loc, uniform); break;
default: console.warn("bad uniform length");
}
} else if (uniform.length == 4 && uniform[0].length == 4)
gl.uniformMatrix4fv(loc, false, new Float32Array(uniform.getAsArray()));
else
console.warn("unsupported uniform matrix");
}
}
}
for (pass = 0; pass < obj.passes; pass++) {
pmode = obj.pmode[pass];
if (pmode === "culled")
continue;
mode = fat_lines && (is_lines || pmode == "lines") ? "TRIANGLES" : this.mode4type[type];
if (depth_sort && pmode == "filled") {// Don't try depthsorting on wireframe or points
var faces = this.depthSort(obj),
nfaces = faces.length,
frowsize = Math.floor(obj.f[pass].length/nfaces);
if (type !== "spheres") {
var f = obj.index_uint ? new Uint32Array(obj.f[pass].length) : new Uint16Array(obj.f[pass].length);
for (i=0; i<nfaces; i++) {
for (j=0; j<frowsize; j++) {
f[frowsize*i + j] = obj.f[pass][frowsize*faces[i] + j];
}
}
gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, obj.ibuf[pass]);
gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, f, gl.DYNAMIC_DRAW);
}
}
if (is_twosided)
gl.uniform1i(obj.frontLoc, pass !== 0);
if (type !== "spheres") 
gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, obj.ibuf[pass]);
if (type === "sprites" || type === "text" || type === "quads") {
count = count * 6/4;
} else if (type === "surface") {
count = obj.f[pass].length;
}
count = obj.f[pass].length;
if (!is_lines && pmode === "lines" && !fat_lines) {
mode = "LINES";
} else if (pmode === "points") {
mode = "POINTS";
}
if ((is_lines || pmode === "lines") && fat_lines) {
gl.enableVertexAttribArray(obj.pointLoc);
enabled.pointLoc = true;
gl.vertexAttribPointer(obj.pointLoc, 2, gl.FLOAT, false, 4*obj.vOffsets.stride, 4*obj.vOffsets.pointofs);
gl.enableVertexAttribArray(obj.nextLoc );
enabled.nextLoc = true;
gl.vertexAttribPointer(obj.nextLoc, 3, gl.FLOAT, false, 4*obj.vOffsets.stride, 4*obj.vOffsets.nextofs);
gl.uniform1f(obj.aspectLoc, this.vp.width/this.vp.height);
gl.uniform1f(obj.lwdLoc, this.getMaterial(id, "lwd")/this.vp.height);
}
gl.vertexAttribPointer(this.posLoc,  3, gl.FLOAT, false, 4*obj.vOffsets.stride,  4*obj.vOffsets.vofs);
gl.drawElements(gl[mode], count, obj.index_uint ? gl.UNSIGNED_INT : gl.UNSIGNED_SHORT, 0);
this.disableArrays(obj, enabled);
}
if (typeof obj.polygon_offset !== "undefined") 
gl.disable(gl.POLYGON_OFFSET_FILL);
};
/**
* Draw the background for a subscene
* @param { number } id - id of background object
* @param { number } subsceneid - id of subscene
*/
rglwidgetClass.prototype.drawBackground = function(id, subsceneid) {
var gl = this.gl || this.initGL(),
obj = this.getObj(id),
bg, i;
if (!obj.initialized)
this.initObj(id);
if (obj.colors.length) {
bg = obj.colors[0];
gl.clearColor(bg[0], bg[1], bg[2], bg[3]);
gl.depthMask(true);
gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
}
if (typeof obj.quad !== "undefined") {
this.prMatrix.makeIdentity();
this.mvMatrix.makeIdentity();
gl.disable(gl.BLEND);
gl.disable(gl.DEPTH_TEST);
gl.depthMask(false);
for (i=0; i < obj.quad.length; i++)
this.drawObj(obj.quad[i], subsceneid);
}
};
/**
* Draw a subscene
* @param { number } subsceneid - id of subscene
* @param { boolean } opaquePass - is this the opaque drawing pass?
*/
rglwidgetClass.prototype.drawSubscene = function(subsceneid, opaquePass) {
var gl = this.gl || this.initGL(),
sub = this.getObj(subsceneid),
objects = this.scene.objects,
subids = sub.objects,
subscene_has_faces = false,
subscene_needs_sorting = false,
flags, i, obj;
if (sub.par3d.skipRedraw)
return;
for (i=0; i < subids.length; i++) {
obj = objects[subids[i]];
flags = obj.flags;
if (typeof flags !== "undefined") {
subscene_has_faces |= (flags & this.f_is_lit)
& !(flags & this.f_fixed_quads);
obj.is_transparent = (flags & this.f_is_transparent) || obj.someHidden;
subscene_needs_sorting |= (flags & this.f_depth_sort) || obj.is_transparent;
}
}
this.setViewport(subsceneid);
if (typeof sub.backgroundId !== "undefined" && opaquePass)
this.drawBackground(sub.backgroundId, subsceneid);
if (subids.length) {
this.setprMatrix(subsceneid);
this.setmvMatrix(subsceneid);
if (subscene_has_faces) {
this.setnormMatrix(subsceneid);
if ((sub.flags & this.f_sprites_3d) &&
typeof sub.spriteNormmat === "undefined") {
sub.spriteNormmat = new CanvasMatrix4(this.normMatrix);
}
}
if (subscene_needs_sorting)
this.setprmvMatrix();
var clipids = sub.clipplanes;
if (typeof clipids === "undefined") {
console.warn("bad clipids");
}
if (clipids.length > 0) {
this.invMatrix = new CanvasMatrix4(this.mvMatrix);
this.invMatrix.invert();
for (i = 0; i < clipids.length; i++)
this.drawObj(clipids[i], subsceneid);
}
subids = sub.opaque.concat(sub.transparent);
if (opaquePass) {
gl.enable(gl.DEPTH_TEST);
gl.depthMask(true);
gl.disable(gl.BLEND);
for (i = 0; i < subids.length; i++) {
if (!this.getObj(subids[i]).is_transparent) 
this.drawObj(subids[i], subsceneid);
}
} else {
gl.depthMask(false);
gl.blendFuncSeparate(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA,
gl.ONE, gl.ONE);
gl.enable(gl.BLEND);
for (i = 0; i < subids.length; i++) {
if (this.getObj(subids[i]).is_transparent)
this.drawObj(subids[i], subsceneid);
}
}
subids = sub.subscenes;
for (i = 0; i < subids.length; i++) {
this.drawSubscene(subids[i], opaquePass);
}
}
};
/**
* Respond to brush change
*/
rglwidgetClass.prototype.selectionChanged = function() {
var i, j, k, id, subid = this.select.subscene, subscene,
objids, obj,
p1 = this.select.region.p1, p2 = this.select.region.p2,
filter, selection = [], handle, keys, xmin, x, xmax, ymin, y, ymax, z, v,
someHidden;
if (!subid)
return;
subscene = this.getObj(subid);
objids = subscene.objects;
filter = this.scene.crosstalk.filter;
this.setmvMatrix(subid);
this.setprMatrix(subid);
this.setprmvMatrix();
xmin = Math.min(p1.x, p2.x);
xmax = Math.max(p1.x, p2.x);
ymin = Math.min(p1.y, p2.y);
ymax = Math.max(p1.y, p2.y);
for (i = 0; i < objids.length; i++) {
id = objids[i];
j = this.scene.crosstalk.id.indexOf(id);
if (j >= 0) {
keys = this.scene.crosstalk.key[j];
obj = this.getObj(id);
someHidden = false;
for (k = 0; k < keys.length; k++) {
if (filter && filter.indexOf(keys[k]) < 0) {
someHidden = true;
continue;
}
v = [].concat(obj.vertices[k]).concat(1.0);
v = this.multVM(v, this.prmvMatrix);
x = v[0]/v[3];
y = v[1]/v[3];
z = v[2]/v[3];
if (xmin <= x && x <= xmax && ymin <= y && y <= ymax && -1.0 <= z && z <= 1.0) {
selection.push(keys[k]);
} else
someHidden = true;
}
obj.someHidden = someHidden && (filter || selection.length);
obj.initialized = false;
/* Who should we notify?  Only shared data in the current subscene, or everyone? */
if (!this.equalArrays(selection, this.scene.crosstalk.selection)) {
handle = this.scene.crosstalk.sel_handle[j];
handle.set(selection, {rglSubsceneId: this.select.subscene});
}
}
}
};
/**
* Respond to selection or filter change from crosstalk
* @param { Object } event - crosstalk event
* @param { boolean } filter - filter or selection?
*/
rglwidgetClass.prototype.selection = function(event, filter) {
var i, j, ids, obj, keys, crosstalk = this.scene.crosstalk,
selection, someHidden;
// Record the message and find out if this event makes some objects have mixed values:
crosstalk = this.scene.crosstalk;
if (filter) {
filter = crosstalk.filter = event.value;
selection = crosstalk.selection;
} else {  
selection = crosstalk.selection = event.value;
filter = crosstalk.filter;
}
ids = crosstalk.id;
for (i = 0; i < ids.length ; i++) {
obj = this.getObj(ids[i]);
obj.initialized = false;
keys = crosstalk.key[i];
someHidden = false;
for (j = 0; j < keys.length && !someHidden; j++) {
if ((filter && filter.indexOf(keys[j]) < 0) ||
(selection.length && selection.indexOf(keys[j]) < 0))
someHidden = true;
}
obj.someHidden = someHidden;
}
this.drawScene();
};
/**
* Clear the selection brush
* @param { number } except - Subscene that should ignore this request
*/
rglwidgetClass.prototype.clearBrush = function(except) {
if (this.select.subscene != except) {
this.select.state = "inactive";
this.delFromSubscene(this.scene.brushId, this.select.subscene);
}
this.drawScene();
};
/**
* Compute mouse coordinates relative to current canvas
* @returns { Object }
* @param { Object } event - event object from mouse click
*/
rglwidgetClass.prototype.relMouseCoords = function(event) {
var totalOffsetX = 0,
totalOffsetY = 0,
currentElement = this.canvas;
do {
totalOffsetX += currentElement.offsetLeft;
totalOffsetY += currentElement.offsetTop;
currentElement = currentElement.offsetParent;
}
while(currentElement);
var canvasX = event.pageX - totalOffsetX,
canvasY = event.pageY - totalOffsetY;
return {x:canvasX, y:canvasY};
};
/**
* Set mouse handlers for the scene
*/
rglwidgetClass.prototype.setMouseHandlers = function() {
var self = this, activeSubscene, handler,
handlers = {}, drag = 0;
handlers.rotBase = 0;
this.screenToVector = function(x, y) {
var viewport = this.getObj(activeSubscene).par3d.viewport,
width = viewport.width*this.canvas.width,
height = viewport.height*this.canvas.height,
radius = Math.max(width, height)/2.0,
cx = width/2.0,
cy = height/2.0,
px = (x-cx)/radius,
py = (y-cy)/radius,
plen = Math.sqrt(px*px+py*py);
if (plen > 1.e-6) {
px = px/plen;
py = py/plen;
}
var angle = (Math.SQRT2 - plen)/Math.SQRT2*Math.PI/2,
z = Math.sin(angle),
zlen = Math.sqrt(1.0 - z*z);
px = px * zlen;
py = py * zlen;
return [px, py, z];
};
handlers.trackballdown = function(x,y) {
var activeSub = this.getObj(activeSubscene),
activeModel = this.getObj(this.useid(activeSub.id, "model")),
i, l = activeModel.par3d.listeners;
handlers.rotBase = this.screenToVector(x, y);
this.saveMat = [];
for (i = 0; i < l.length; i++) {
activeSub = this.getObj(l[i]);
activeSub.saveMat = new CanvasMatrix4(activeSub.par3d.userMatrix);
}
};
handlers.trackballmove = function(x,y) {
var rotCurrent = this.screenToVector(x,y),
rotBase = handlers.rotBase,
dot = rotBase[0]*rotCurrent[0] +
rotBase[1]*rotCurrent[1] +
rotBase[2]*rotCurrent[2],
angle = Math.acos( dot/this.vlen(rotBase)/this.vlen(rotCurrent) )*180.0/Math.PI,
axis = this.xprod(rotBase, rotCurrent),
objects = this.scene.objects,
activeSub = this.getObj(activeSubscene),
activeModel = this.getObj(this.useid(activeSub.id, "model")),
l = activeModel.par3d.listeners,
i;
for (i = 0; i < l.length; i++) {
activeSub = this.getObj(l[i]);
activeSub.par3d.userMatrix.load(objects[l[i]].saveMat);
activeSub.par3d.userMatrix.rotate(angle, axis[0], axis[1], axis[2]);
}
this.drawScene();
};
handlers.trackballend = 0;
this.clamp = function(x, lo, hi) {
return Math.max(lo, Math.min(x, hi));
};
this.screenToPolar = function(x,y) {
var viewport = this.getObj(activeSubscene).par3d.viewport,
width = viewport.width*this.canvas.width,
height = viewport.height*this.canvas.height,
r = Math.min(width, height)/2,
dx = this.clamp(x - width/2, -r, r),
dy = this.clamp(y - height/2, -r, r);
return [Math.asin(dx/r), Math.asin(-dy/r)];
};
handlers.polardown = function(x,y) {
var activeSub = this.getObj(activeSubscene),
activeModel = this.getObj(this.useid(activeSub.id, "model")),
i, l = activeModel.par3d.listeners;
handlers.dragBase = this.screenToPolar(x, y);
this.saveMat = [];
for (i = 0; i < l.length; i++) {
activeSub = this.getObj(l[i]);
activeSub.saveMat = new CanvasMatrix4(activeSub.par3d.userMatrix);
activeSub.camBase = [-Math.atan2(activeSub.saveMat.m13, activeSub.saveMat.m11),
Math.atan2(activeSub.saveMat.m32, activeSub.saveMat.m22)];
}
};
handlers.polarmove = function(x,y) {
var dragCurrent = this.screenToPolar(x,y),
activeSub = this.getObj(activeSubscene),
activeModel = this.getObj(this.useid(activeSub.id, "model")),
objects = this.scene.objects,
l = activeModel.par3d.listeners,
i, changepos = [];
for (i = 0; i < l.length; i++) {
activeSub = this.getObj(l[i]);
for (j=0; j<2; j++)
changepos[j] = -(dragCurrent[j] - handlers.dragBase[j]);
activeSub.par3d.userMatrix.makeIdentity();
activeSub.par3d.userMatrix.rotate(changepos[0]*180/Math.PI, 0,-1,0);
activeSub.par3d.userMatrix.multRight(objects[l[i]].saveMat);
activeSub.par3d.userMatrix.rotate(changepos[1]*180/Math.PI, -1,0,0);
}
this.drawScene();
};
handlers.polarend = 0;
handlers.axisdown = function(x,y) {
handlers.rotBase = this.screenToVector(x, this.canvas.height/2);
var activeSub = this.getObj(activeSubscene),
activeModel = this.getObj(this.useid(activeSub.id, "model")),
i, l = activeModel.par3d.listeners;
for (i = 0; i < l.length; i++) {
activeSub = this.getObj(l[i]);
activeSub.saveMat = new CanvasMatrix4(activeSub.par3d.userMatrix);
}
};
handlers.axismove = function(x,y) {
var rotCurrent = this.screenToVector(x, this.canvas.height/2),
rotBase = handlers.rotBase,
angle = (rotCurrent[0] - rotBase[0])*180/Math.PI,
rotMat = new CanvasMatrix4();
rotMat.rotate(angle, handlers.axis[0], handlers.axis[1], handlers.axis[2]);
var activeSub = this.getObj(activeSubscene),
activeModel = this.getObj(this.useid(activeSub.id, "model")),
i, l = activeModel.par3d.listeners;
for (i = 0; i < l.length; i++) {
activeSub = this.getObj(l[i]);
activeSub.par3d.userMatrix.load(activeSub.saveMat);
activeSub.par3d.userMatrix.multLeft(rotMat);
}
this.drawScene();
};
handlers.axisend = 0;
handlers.y0zoom = 0;
handlers.zoom0 = 0;
handlers.zoomdown = function(x, y) {
var activeSub = this.getObj(activeSubscene),
activeProjection = this.getObj(this.useid(activeSub.id, "projection")),
i, l = activeProjection.par3d.listeners;
handlers.y0zoom = y;
for (i = 0; i < l.length; i++) {
activeSub = this.getObj(l[i]);
activeSub.zoom0 = Math.log(activeSub.par3d.zoom);
}
};
handlers.zoommove = function(x, y) {
var activeSub = this.getObj(activeSubscene),
activeProjection = this.getObj(this.useid(activeSub.id, "projection")),
i, l = activeProjection.par3d.listeners;
for (i = 0; i < l.length; i++) {
activeSub = this.getObj(l[i]);
activeSub.par3d.zoom = Math.exp(activeSub.zoom0 + (y-handlers.y0zoom)/this.canvas.height);
}
this.drawScene();
};
handlers.zoomend = 0;
handlers.y0fov = 0;
handlers.fovdown = function(x, y) {
handlers.y0fov = y;
var activeSub = this.getObj(activeSubscene),
activeProjection = this.getObj(this.useid(activeSub.id, "projection")),
i, l = activeProjection.par3d.listeners;
for (i = 0; i < l.length; i++) {
activeSub = this.getObj(l[i]);
activeSub.fov0 = activeSub.par3d.FOV;
}
};
handlers.fovmove = function(x, y) {
var activeSub = this.getObj(activeSubscene),
activeProjection = this.getObj(this.useid(activeSub.id, "projection")),
i, l = activeProjection.par3d.listeners;
for (i = 0; i < l.length; i++) {
activeSub = this.getObj(l[i]);
activeSub.par3d.FOV = Math.max(1, Math.min(179, activeSub.fov0 +
180*(y-handlers.y0fov)/this.canvas.height));
}
this.drawScene();
};
handlers.fovend = 0;
handlers.selectingdown = function(x, y) {
var viewport = this.getObj(activeSubscene).par3d.viewport,
width = viewport.width*this.canvas.width,
height = viewport.height*this.canvas.height, 
p = {x: 2.0*x/width - 1.0, y: 2.0*y/height - 1.0};
this.select.region = {p1: p, p2: p};
if (this.select.subscene && this.select.subscene != activeSubscene)
this.delFromSubscene(this.scene.brushId, this.select.subscene);
this.select.subscene = activeSubscene;
this.addToSubscene(this.scene.brushId, activeSubscene);
this.select.state = "changing";
if (typeof this.scene.brushId !== "undefined")
this.getObj(this.scene.brushId).initialized = false;
this.drawScene();
};
handlers.selectingmove = function(x, y) {
var viewport = this.getObj(activeSubscene).par3d.viewport,
width = viewport.width*this.canvas.width,
height = viewport.height*this.canvas.height;
if (this.select.state === "inactive") 
return;
this.select.region.p2 = {x: 2.0*x/width - 1.0, y: 2.0*y/height - 1.0};
if (typeof this.scene.brushId !== "undefined")
this.getObj(this.scene.brushId).initialized = false;
this.drawScene();
};
handlers.selectingend = 0;
this.canvas.onmousedown = function ( ev ){
if (!ev.which) // Use w3c defns in preference to MS
switch (ev.button) {
case 0: ev.which = 1; break;
case 1:
case 4: ev.which = 2; break;
case 2: ev.which = 3;
}
drag = ["left", "middle", "right"][ev.which-1];
var coords = self.relMouseCoords(ev);
coords.y = self.canvas.height-coords.y;
activeSubscene = self.whichSubscene(coords);
var sub = self.getObj(activeSubscene), f;
handler = sub.par3d.mouseMode[drag];
switch (handler) {
case "xAxis":
handler = "axis";
handlers.axis = [1.0, 0.0, 0.0];
break;
case "yAxis":
handler = "axis";
handlers.axis = [0.0, 1.0, 0.0];
break;
case "zAxis":
handler = "axis";
handlers.axis = [0.0, 0.0, 1.0];
break;
}
f = handlers[handler + "down"];
if (f) {
coords = self.translateCoords(activeSubscene, coords);
f.call(self, coords.x, coords.y);
ev.preventDefault();
} else
console.warn("Mouse handler '" + handler + "' is not implemented.");
};
this.canvas.onmouseup = function ( ev ){
if ( drag === 0 ) return;
var f = handlers[handler + "end"];
if (f) {
f.call(self);
ev.preventDefault();
}
drag = 0;
};
this.canvas.onmouseout = this.canvas.onmouseup;
this.canvas.onmousemove = function ( ev ) {
if ( drag === 0 ) return;
var f = handlers[handler + "move"];
if (f) {
var coords = self.relMouseCoords(ev);
coords.y = self.canvas.height - coords.y;
coords = self.translateCoords(activeSubscene, coords);
f.call(self, coords.x, coords.y);
}
};
handlers.wheelHandler = function(ev) {
var del = 1.02, i;
if (ev.shiftKey) del = 1.002;
var ds = ((ev.detail || ev.wheelDelta) > 0) ? del : (1 / del);
if (typeof activeSubscene === "undefined")
activeSubscene = self.scene.rootSubscene;
var activeSub = self.getObj(activeSubscene),
activeProjection = self.getObj(self.useid(activeSub.id, "projection")),
l = activeProjection.par3d.listeners;
for (i = 0; i < l.length; i++) {
activeSub = self.getObj(l[i]);
activeSub.par3d.zoom *= ds;
}
self.drawScene();
ev.preventDefault();
};
this.canvas.addEventListener("DOMMouseScroll", handlers.wheelHandler, false);
this.canvas.addEventListener("mousewheel", handlers.wheelHandler, false);
};
/**
* Find a particular subscene by inheritance
* @returns { number } id of subscene to use
* @param { number } subsceneid - child subscene
* @param { string } type - type of inheritance:  "projection" or "model"
*/
rglwidgetClass.prototype.useid = function(subsceneid, type) {
var sub = this.getObj(subsceneid);
if (sub.embeddings[type] === "inherit")
return(this.useid(sub.parent, type));
else
return subsceneid;
};
/**
* Check whether point is in viewport of subscene
* @returns {boolean}
* @param { Object } coords - screen coordinates of point
* @param { number } subsceneid - subscene to check
*/
rglwidgetClass.prototype.inViewport = function(coords, subsceneid) {
var viewport = this.getObj(subsceneid).par3d.viewport,
x0 = coords.x - viewport.x*this.canvas.width,
y0 = coords.y - viewport.y*this.canvas.height;
return 0 <= x0 && x0 <= viewport.width*this.canvas.width &&
0 <= y0 && y0 <= viewport.height*this.canvas.height;
};
/**
* Find which subscene contains a point
* @returns { number } subscene id
* @param { Object } coords - coordinates of point
*/
rglwidgetClass.prototype.whichSubscene = function(coords) {
var self = this,
recurse = function(subsceneid) {
var subscenes = self.getChildSubscenes(subsceneid), i, id;
for (i=0; i < subscenes.length; i++) {
id = recurse(subscenes[i]);
if (typeof(id) !== "undefined")
return(id);
}
if (self.inViewport(coords, subsceneid))
return(subsceneid);
else
return undefined;
},
rootid = this.scene.rootSubscene,
result = recurse(rootid);
if (typeof(result) === "undefined")
result = rootid;
return result;
};
/**
* Translate from window coordinates to viewport coordinates
* @returns { Object } translated coordinates
* @param { number } subsceneid - which subscene to use?
* @param { Object } coords - point to translate
*/
rglwidgetClass.prototype.translateCoords = function(subsceneid, coords) {
var viewport = this.getObj(subsceneid).par3d.viewport;
return {x: coords.x - viewport.x*this.canvas.width,
y: coords.y - viewport.y*this.canvas.height};
};
/**
* Initialize the sphere object
*/
rglwidgetClass.prototype.initSphere = function() {
var verts = this.scene.sphereVerts,
reuse = verts.reuse, result;
if (typeof reuse !== "undefined") {
var prev = document.getElementById(reuse).rglinstance.sphere;
result = {values: prev.values, vOffsets: prev.vOffsets, it: prev.it};
} else
result = {values: new Float32Array(this.flatten(this.cbind(this.transpose(verts.vb),
this.transpose(verts.texcoords)))),
it: new Uint16Array(this.flatten(this.transpose(verts.it))),
vOffsets: {vofs:0, cofs:-1, nofs:-1, radofs:-1, oofs:-1,
tofs:3, nextofs:-1, pointofs:-1, stride:5}};
result.sphereCount = result.it.length;
this.sphere = result;
};
/**
* Set the vertices in the selection box object
*/
rglwidgetClass.prototype.initSelection = function(id) {
if (typeof this.select.region === "undefined")
return;
var obj = this.getObj(id),
width = this.canvas.width,
height = this.canvas.height, 
p1 = this.select.region.p1,
p2 = this.select.region.p2;
obj.vertices = [[p1.x, p1.y, 0.0],
[p2.x, p1.y, 0.0],
[p2.x, p2.y, 0.0],
[p1.x, p2.y, 0.0],
[p1.x, p1.y, 0.0]];
};
/**
* Do the gl part of initializing the sphere
*/
rglwidgetClass.prototype.initSphereGL = function() {
var gl = this.gl || this.initGL(), sphere = this.sphere;
if (gl.isContextLost()) return;
sphere.buf = gl.createBuffer();
gl.bindBuffer(gl.ARRAY_BUFFER, sphere.buf);
gl.bufferData(gl.ARRAY_BUFFER, sphere.values, gl.STATIC_DRAW);
sphere.ibuf = gl.createBuffer();
gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, sphere.ibuf);
gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, sphere.it, gl.STATIC_DRAW);
return;
};
/**
* Initialize the DOM object
* @param { Object } el - the DOM object
* @param { Object } x - the scene data sent by JSON from R
*/
rglwidgetClass.prototype.initialize = function(el, x) {
this.textureCanvas = document.createElement("canvas");
this.textureCanvas.style.display = "block";
this.scene = x;
this.normMatrix = new CanvasMatrix4();
this.saveMat = {};
this.distance = null;
this.posLoc = 0;
this.colLoc = 1;
if (el) {
el.rglinstance = this;
this.el = el;
this.webGLoptions = el.rglinstance.scene.webGLoptions;
this.initCanvas();
}
if (typeof Shiny !== "undefined") {
var self = this;
Shiny.addCustomMessageHandler("shinyGetPar3d",
function(message) {
var i, param, 
subscene = self.getObj(message.subscene),
parameters = [].concat(message.parameters),
result = {tag: message.tag, subscene: message.subscene};
if (typeof subscene !== "undefined") {
for (i = 0; i < parameters.length; i++) {
param = parameters[i];
result[param] = subscene.par3d[param];
};
} else {
console.log("subscene "+message.subscene+" undefined.")
}
Shiny.setInputValue("par3d:shinyPar3d", result, {priority: "event"});
});
Shiny.addCustomMessageHandler("shinySetPar3d",
function(message) {
var param = message.parameter, 
subscene = self.getObj(message.subscene);
if (typeof subscene !== "undefined") {
subscene.par3d[param] = message.value;
subscene.initialized = false;
self.drawScene();
} else {
console.log("subscene "+message.subscene+" undefined.")
}
})
}
};
/**
* Restart the WebGL canvas
*/
rglwidgetClass.prototype.restartCanvas = function() {
var newcanvas = document.createElement("canvas"),
self = this;
newcanvas.width = this.el.width;
newcanvas.height = this.el.height;
newcanvas.addEventListener("webglcontextrestored",
this.onContextRestored, false);
newcanvas.addEventListener("webglcontextlost",
this.onContextLost, false);
while (this.el.firstChild) {
this.el.removeChild(this.el.firstChild);
}
this.el.appendChild(newcanvas);
this.canvas = newcanvas;
this.setMouseHandlers();
if (this.gl) 
Object.keys(this.scene.objects).forEach(function(key){
self.getObj(parseInt(key, 10)).texture = undefined; 
});
this.gl = null;
};
/**
* Initialize the WebGL canvas
*/
rglwidgetClass.prototype.initCanvas = function() {
this.restartCanvas();
var objs = this.scene.objects,
self = this;
Object.keys(objs).forEach(function(key){
var id = parseInt(key, 10),
obj = self.getObj(id);
if (typeof obj.reuse !== "undefined")
self.copyObj(id, obj.reuse);
});
Object.keys(objs).forEach(function(key){
self.initSubscene(parseInt(key, 10));
});
this.setMouseHandlers();
this.initSphere();
this.onContextRestored = function(event) {
self.initGL();
self.drawScene();
};
this.onContextLost = function(event) {
if (!self.drawing)
this.gl = null;
event.preventDefault();
};
this.initGL0();
this.lazyLoadScene = function() {
if (typeof self.slide === "undefined")
self.slide = self.getSlide();
if (self.isInBrowserViewport()) {
if (!self.gl || self.gl.isContextLost())
self.initGL();
self.drawScene();
}
};
window.addEventListener("DOMContentLoaded", this.lazyLoadScene, false);
window.addEventListener("load", this.lazyLoadScene, false);
window.addEventListener("resize", this.lazyLoadScene, false);
window.addEventListener("scroll", this.lazyLoadScene, false);
this.slide = this.getSlide();
if (this.slide) {
if (typeof this.slide.rgl === "undefined")
this.slide.rgl = [this];
else
this.slide.rgl.push(this);
if (this.scene.context.rmarkdown) 
if (this.scene.context.rmarkdown === "ioslides_presentation") {
this.slide.setAttribute("slideenter", "this.rgl.forEach(function(scene) { scene.lazyLoadScene.call(window);})");
} else if (this.scene.context.rmarkdown === "slidy_presentation") {
// This method would also work in ioslides, but it gets triggered
// something like 5 times per slide for every slide change, so
// you'd need a quicker function than lazyLoadScene.
var MutationObserver = window.MutationObserver || window.WebKitMutationObserver || window.MozMutationObserver,
observer = new MutationObserver(function(mutations) {
mutations.forEach(function(mutation) {
self.slide.rgl.forEach(function(scene) { scene.lazyLoadScene.call(window); });});});
observer.observe(this.slide, { attributes: true, attributeFilter:["class"] });
}
}
};
/**
* Start the writeWebGL scene. This is only used by writeWebGL; rglwidget has
no debug element and does the drawing in rglwidget.js.
*/
rglwidgetClass.prototype.start = function() {
if (typeof this.prefix !== "undefined") {
this.debugelement = document.getElementById(this.prefix + "debug");
this.debug("");
}
this.drag = 0;
this.drawScene();
};
/**
* Display a debug message
* @param { string } msg - The message to display
* @param { Object } [img] - Image to insert before message
*/
rglwidgetClass.prototype.debug = function(msg, img) {
if (typeof this.debugelement !== "undefined" && this.debugelement !== null) {
this.debugelement.innerHTML = msg;
if (typeof img !== "undefined") {
this.debugelement.insertBefore(img, this.debugelement.firstChild);
}
} else if (msg !== "")
alert(msg);
};
/**
* Get the snapshot image of this scene
* @returns { Object } The img DOM element
*/
rglwidgetClass.prototype.getSnapshot = function() {
var img;
if (typeof this.scene.snapshot !== "undefined") {
img = document.createElement("img");
img.src = this.scene.snapshot;
img.alt = "Snapshot";
}
return img;
};
/**
* Initial test for WebGL
*/
rglwidgetClass.prototype.initGL0 = function() {
if (!window.WebGLRenderingContext){
alert("Your browser does not support WebGL. See http://get.webgl.org");
return;
}
};
/**
* If we are in an ioslides or slidy presentation, get the
* DOM element of the current slide
* @returns { Object }
*/
rglwidgetClass.prototype.getSlide = function() {
var result = this.el, done = false;
while (result && !done && this.scene.context.rmarkdown) {
switch(this.scene.context.rmarkdown) {
case "ioslides_presentation":
if (result.tagName === "SLIDE") return result;
break;
case "slidy_presentation":
if (result.tagName === "DIV" && result.classList.contains("slide"))
return result;
break;
default: return null;
}
result = result.parentElement;
}
return null;
};
/**
* Is this scene visible in the browser?
* @returns { boolean }
*/
rglwidgetClass.prototype.isInBrowserViewport = function() {
var rect = this.canvas.getBoundingClientRect(),
windHeight = (window.innerHeight || document.documentElement.clientHeight),
windWidth = (window.innerWidth || document.documentElement.clientWidth);
if (this.scene.context && this.scene.context.rmarkdown !== null) {
if (this.slide)
return (this.scene.context.rmarkdown === "ioslides_presentation" &&
this.slide.classList.contains("current")) ||
(this.scene.context.rmarkdown === "slidy_presentation" &&
!this.slide.classList.contains("hidden"));
}
return (
rect.top >= -windHeight &&
rect.left >= -windWidth &&
rect.bottom <= 2*windHeight &&
rect.right <= 2*windWidth);
};
/**
* Initialize WebGL
* @returns { Object } the WebGL context
*/
rglwidgetClass.prototype.initGL = function() {
var self = this;
if (this.gl) {
if (!this.drawing && this.gl.isContextLost())
this.restartCanvas();
else
return this.gl;
}
// if (!this.isInBrowserViewport()) return; Return what??? At this point we know this.gl is null.
this.canvas.addEventListener("webglcontextrestored",
this.onContextRestored, false);
this.canvas.addEventListener("webglcontextlost",
this.onContextLost, false);
this.gl = this.canvas.getContext("webgl", this.webGLoptions) ||
this.canvas.getContext("experimental-webgl", this.webGLoptions);
this.index_uint = this.gl.getExtension("OES_element_index_uint");
var save = this.startDrawing();
this.initSphereGL();
Object.keys(this.scene.objects).forEach(function(key){
self.initObj(parseInt(key, 10));
});
this.stopDrawing(save);
return this.gl;
};
/**
* Resize the display to match element
* @param { Object } el - DOM element to match
*/
rglwidgetClass.prototype.resize = function(el) {
this.canvas.width = el.width;
this.canvas.height = el.height;
};
/**
* Draw the whole scene
*/
rglwidgetClass.prototype.drawScene = function() {
var gl = this.gl || this.initGL(),
wasDrawing = this.startDrawing();
if (!wasDrawing) {
if (this.select.state !== "inactive")
this.selectionChanged();
gl.enable(gl.DEPTH_TEST);
gl.depthFunc(gl.LEQUAL);
gl.clearDepth(1.0);
gl.clearColor(1,1,1,1);
gl.depthMask(true); // Must be true before clearing depth buffer
gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
this.drawSubscene(this.scene.rootSubscene, true);
this.drawSubscene(this.scene.rootSubscene, false);
}
this.stopDrawing(wasDrawing);
};
/**
* Change the displayed subset
* @param { Object } el - Element of the control; not used.
* @param { Object } control - The subset control data.
*/
rglwidgetClass.prototype.subsetSetter = function(el, control) {
if (typeof control.subscenes === "undefined" ||
control.subscenes === null)
control.subscenes = this.scene.rootSubscene;
var value = Math.round(control.value),
subscenes = [].concat(control.subscenes),
fullset = [].concat(control.fullset),
i, j, entries, subsceneid,
adds = [], deletes = [],
ismissing = function(x) {
return fullset.indexOf(x) < 0;
},
tointeger = function(x) {
return parseInt(x, 10);
};
if (isNaN(value))
value = control.value = 0;
if (control.accumulate)
for (i=0; i <= value; i++)
adds = adds.concat(control.subsets[i]);
else
adds = adds.concat(control.subsets[value]);
deletes = fullset.filter(function(x) { return adds.indexOf(x) < 0; });
for (i = 0; i < subscenes.length; i++) {
subsceneid = subscenes[i];
if (typeof this.getObj(subsceneid) === "undefined")
this.alertOnce("typeof object is undefined");
for (j = 0; j < adds.length; j++)
this.addToSubscene(adds[j], subsceneid);
for (j = 0; j < deletes.length; j++)
this.delFromSubscene(deletes[j], subsceneid);
}
};
/**
* Change the requested property
* @param { Object } el - Element of the control; not used.
* @param { Object } control - The property setter control data.
*/
rglwidgetClass.prototype.propertySetter = function(el, control)  {
var value = control.value,
values = [].concat(control.values),
svals = [].concat(control.param),
direct = values[0] === null,
entries = [].concat(control.entries),
ncol = entries.length,
nrow = values.length/ncol,
properties = this.repeatToLen(control.properties, ncol),
objids = this.repeatToLen(control.objids, ncol),
property, objid = objids[0],
obj = this.getObj(objid),
propvals, i, v1, v2, p, entry, gl, needsBinding,
newprop, newid,
getPropvals = function() {
if (property === "userMatrix")
return obj.par3d.userMatrix.getAsArray();
else if (property === "scale" || property === "FOV" || property === "zoom")
return [].concat(obj.par3d[property]);
else
return [].concat(obj[property]);
};
putPropvals = function(newvals) {
if (newvals.length == 1)
newvals = newvals[0];
if (property === "userMatrix")
obj.par3d.userMatrix.load(newvals);
else if (property === "scale" || property === "FOV" || property === "zoom")
obj.par3d[property] = newvals;
else
obj[property] = newvals;
};
if (direct && typeof value === "undefined")
return;
if (control.interp) {
values = values.slice(0, ncol).concat(values).
concat(values.slice(ncol*(nrow-1), ncol*nrow));
svals = [-Infinity].concat(svals).concat(Infinity);
for (i = 1; i < svals.length; i++) {
if (value <= svals[i]) {
if (svals[i] === Infinity)
p = 1;
else
p = (svals[i] - value)/(svals[i] - svals[i-1]);
break;
}
}
} else if (!direct) {
value = Math.round(value);
}
for (j=0; j<entries.length; j++) {
entry = entries[j];
newprop = properties[j];
newid = objids[j];
if (newprop !== property || newid != objid) {
if (typeof property !== "undefined")
putPropvals(propvals);
property = newprop;
objid = newid;
obj = this.getObj(objid);
propvals = getPropvals();
}
if (control.interp) {
v1 = values[ncol*(i-1) + j];
v2 = values[ncol*i + j];
this.setElement(propvals, entry, p*v1 + (1-p)*v2);
} else if (!direct) {
this.setElement(propvals, entry, values[ncol*value + j]);
} else {
this.setElement(propvals, entry, value[j]);
}
}
putPropvals(propvals);
needsBinding = [];
for (j=0; j < entries.length; j++) {
if (properties[j] === "values" &&
needsBinding.indexOf(objids[j]) === -1) {
needsBinding.push(objids[j]);
}
}
for (j=0; j < needsBinding.length; j++) {
gl = this.gl || this.initGL();
obj = this.getObj(needsBinding[j]);
gl.bindBuffer(gl.ARRAY_BUFFER, obj.buf);
gl.bufferData(gl.ARRAY_BUFFER, obj.values, gl.STATIC_DRAW);
}
};
/**
* Change the requested vertices
* @param { Object } el - Element of the control; not used.
* @param { Object } control - The vertext setter control data.
*/
rglwidgetClass.prototype.vertexSetter = function(el, control)  {
var svals = [].concat(control.param),
j, k, p, a, propvals, stride, ofs, obj, entry,
attrib,
ofss    = {x:"vofs", y:"vofs", z:"vofs",
red:"cofs", green:"cofs", blue:"cofs",
alpha:"cofs", radii:"radofs",
nx:"nofs", ny:"nofs", nz:"nofs",
ox:"oofs", oy:"oofs", oz:"oofs",
ts:"tofs", tt:"tofs"},
pos     = {x:0, y:1, z:2,
red:0, green:1, blue:2,
alpha:3,radii:0,
nx:0, ny:1, nz:2,
ox:0, oy:1, oz:2,
ts:0, tt:1},
values = control.values,
direct = values === null,
ncol,
interp = control.interp,
vertices = [].concat(control.vertices),
attributes = [].concat(control.attributes),
value = control.value, newval, aliases, alias;
ncol = Math.max(vertices.length, attributes.length);
if (!ncol)
return;
vertices = this.repeatToLen(vertices, ncol);
attributes = this.repeatToLen(attributes, ncol);
if (direct)
interp = false;
/* JSON doesn't pass Infinity */
svals[0] = -Infinity;
svals[svals.length - 1] = Infinity;
for (j = 1; j < svals.length; j++) {
if (value <= svals[j]) {
if (interp) {
if (svals[j] === Infinity)
p = 1;
else
p = (svals[j] - value)/(svals[j] - svals[j-1]);
} else {
if (svals[j] - value > value - svals[j-1])
j = j - 1;
}
break;
}
}
obj = this.getObj(control.objid);
// First, make sure color attributes vary in original
if (typeof obj.vOffsets !== "undefined") {
varies = true;
for (k = 0; k < ncol; k++) {
attrib = attributes[k];
if (typeof attrib !== "undefined") {
ofs = obj.vOffsets[ofss[attrib]];
if (ofs < 0) {
switch(attrib) {
case "alpha":
case "red":
case "green":
case "blue":
obj.colors = [obj.colors[0], obj.colors[0]];
break;
}
varies = false;
}
}
}
if (!varies)
this.initObj(control.objid);
}
propvals = obj.values;
aliases = obj.alias;
if (typeof aliases === "undefined")
aliases = [];
for (k=0; k<ncol; k++) {
if (direct) {
newval = value;
} else if (interp) {
newval = p*values[j-1][k] + (1-p)*values[j][k];
} else {
newval = values[j][k];
}       
attrib = attributes[k];
vertex = vertices[k];
alias = aliases[vertex];
if (obj.type === "planes" || obj.type === "clipplanes") {
ofs = ["nx", "ny", "nz", "offset"].indexOf(attrib);
if (ofs >= 0) {
if (ofs < 3) {
if (obj.normals[vertex][ofs] != newval) {  // Assume no aliases here...
obj.normals[vertex][ofs] = newval;
obj.initialized = false;
}
} else {
if (obj.offsets[vertex][0] != newval) {
obj.offsets[vertex][0] = newval;
obj.initialized = false;
}
}
continue;
}
}
// Not a plane setting...
ofs = obj.vOffsets[ofss[attrib]];
if (ofs < 0)
this.alertOnce("Attribute '"+attrib+"' not found in object "+control.objid);
else {
stride = obj.vOffsets.stride;
ofs = ofs + pos[attrib];
entry = vertex*stride + ofs;
propvals[entry] = newval;
if (typeof alias !== "undefined")
for (a = 0; a < alias.length; a++)
propvals[alias[a]*stride + ofs] = newval;
}
}
if (typeof obj.buf !== "undefined") {
var gl = this.gl || this.initGL();
gl.bindBuffer(gl.ARRAY_BUFFER, obj.buf);
gl.bufferData(gl.ARRAY_BUFFER, propvals, gl.STATIC_DRAW);
}
};
/**
* Change the requested vertex properties by age
* @param { Object } el - Element of the control; not used.
* @param { Object } control - The age setter control data.
*/
rglwidgetClass.prototype.ageSetter = function(el, control) {
var objids = [].concat(control.objids),
nobjs = objids.length,
time = control.value,
births = [].concat(control.births),
ages = [].concat(control.ages),
steps = births.length,
j = Array(steps),
p = Array(steps),
i, k, age, j0, propvals, stride, ofs, objid, obj,
attrib, dim, varies, alias, aliases, a, d,
attribs = ["colors", "alpha", "radii", "vertices",
"normals", "origins", "texcoords",
"x", "y", "z",
"red", "green", "blue"],
ofss    = ["cofs", "cofs", "radofs", "vofs",
"nofs", "oofs", "tofs",
"vofs", "vofs", "vofs",
"cofs", "cofs", "cofs"],
dims    = [3,1,1,3,
3,2,2,
1,1,1,
1,1,1],
pos     = [0,3,0,0,
0,0,0,
0,1,2,
0,1,2];
/* Infinity doesn't make it through JSON */
ages[0] = -Infinity;
ages[ages.length-1] = Infinity;
for (i = 0; i < steps; i++) {
if (births[i] !== null) {  // NA in R becomes null
age = time - births[i];
for (j0 = 1; age > ages[j0]; j0++);
if (ages[j0] == Infinity)
p[i] = 1;
else if (ages[j0] > ages[j0-1])
p[i] = (ages[j0] - age)/(ages[j0] - ages[j0-1]);
else
p[i] = 0;
j[i] = j0;
}
}
// First, make sure color attributes vary in original
for (l = 0; l < nobjs; l++) {
objid = objids[l];
obj = this.getObj(objid);
varies = true;
if (typeof obj.vOffsets === "undefined")
continue;
for (k = 0; k < attribs.length; k++) {
attrib = control[attribs[k]];
if (typeof attrib !== "undefined") {
ofs = obj.vOffsets[ofss[k]];
if (ofs < 0) {
switch(attribs[k]) {
case "colors":
case "alpha":
case "red":
case "green":
case "blue":
obj.colors = [obj.colors[0], obj.colors[0]];
break;
}
varies = false;
}
}
}
if (!varies)
this.initObj(objid);
}
for (l = 0; l < nobjs; l++) {
objid = objids[l];
obj = this.getObj(objid);
if (typeof obj.vOffsets === "undefined")
continue;
aliases = obj.alias;
if (typeof aliases === "undefined")
aliases = [];
propvals = obj.values;
stride = obj.vOffsets.stride;
for (k = 0; k < attribs.length; k++) {
attrib = control[attribs[k]];
if (typeof attrib !== "undefined") {
ofs = obj.vOffsets[ofss[k]];
if (ofs >= 0) {
dim = dims[k];
ofs = ofs + pos[k];
for (i = 0; i < steps; i++) {
alias = aliases[i];
if (births[i] !== null) {
for (d=0; d < dim; d++) {
propvals[i*stride + ofs + d] = p[i]*attrib[dim*(j[i]-1) + d] + (1-p[i])*attrib[dim*j[i] + d];
if (typeof alias !== "undefined")
for (a=0; a < alias.length; a++)
propvals[alias[a]*stride + ofs + d] = propvals[i*stride + ofs + d];
}
}
}
} else
this.alertOnce("\'"+attribs[k]+"\' property not found in object "+objid);
}
}
obj.values = propvals;
if (typeof obj.buf !== "undefined") {
gl = this.gl || this.initGL();
gl.bindBuffer(gl.ARRAY_BUFFER, obj.buf);
gl.bufferData(gl.ARRAY_BUFFER, obj.values, gl.STATIC_DRAW);
}
}
};
/**
* Bridge to old style control
* @param { Object } el - Element of the control; not used.
* @param { Object } control - The bridge control data.
*/
rglwidgetClass.prototype.oldBridge = function(el, control) {
var attrname, global = window[control.prefix + "rgl"];
if (global)
for (attrname in global)
this[attrname] = global[attrname];
window[control.prefix + "rgl"] = this;
};
/**
* Set up a player control
* @param { Object } el - The player control element
* @param { Object } control - The player data.
*/
rglwidgetClass.prototype.Player = function(el, control) {
var
self = this,
components = [].concat(control.components),
buttonLabels = [].concat(control.buttonLabels),
Tick = function() { /* "this" will be a timer */
var i,
nominal = this.value,
slider = this.Slider,
labels = this.outputLabels,
output = this.Output,
step;
if (typeof slider !== "undefined" && nominal != slider.value)
slider.value = nominal;
if (typeof output !== "undefined") {
step = Math.round((nominal - output.sliderMin)/output.sliderStep);
if (labels !== null) {
output.innerHTML = labels[step];
} else {
step = step*output.sliderStep + output.sliderMin;
output.innerHTML = step.toPrecision(output.outputPrecision);
}
}
for (i=0; i < this.actions.length; i++) {
this.actions[i].value = nominal;
}
self.applyControls(el, this.actions, false);
self.drawScene();
},
OnSliderInput = function() { /* "this" will be the slider */
this.rgltimer.value = Number(this.value);
this.rgltimer.Tick();
},
addSlider = function(min, max, step, value) {
var slider = document.createElement("input");
slider.type = "range";
slider.min = min;
slider.max = max;
slider.step = step;
slider.value = value;
slider.oninput = OnSliderInput;
slider.sliderActions = control.actions;
slider.sliderScene = this;
slider.className = "rgl-slider";
slider.id = el.id + "-slider";
el.rgltimer.Slider = slider;
slider.rgltimer = el.rgltimer;
el.appendChild(slider);
},
addLabel = function(labels, min, step, precision) {
var output = document.createElement("output");
output.sliderMin = min;
output.sliderStep = step;
output.outputPrecision = precision;
output.className = "rgl-label";
output.id = el.id + "-label";
el.rgltimer.Output = output;
el.rgltimer.outputLabels = labels;
el.appendChild(output);
},
addButton = function(which, label, active) {
var button = document.createElement("input"),
onclicks = {Reverse: function() { this.rgltimer.reverse();},
Play: function() { this.rgltimer.play();
this.value = this.rgltimer.enabled ? this.inactiveValue : this.activeValue; },
Slower: function() { this.rgltimer.slower(); },
Faster: function() { this.rgltimer.faster(); },
Reset: function() { this.rgltimer.reset(); },
Step:  function() { this.rgltimer.step(); }
};
button.rgltimer = el.rgltimer;
button.type = "button";
button.value = label;
button.activeValue = label;
button.inactiveValue = active;
if (which === "Play")
button.rgltimer.PlayButton = button;
button.onclick = onclicks[which];
button.className = "rgl-button";
button.id = el.id + "-" + which;
el.appendChild(button);
};
if (typeof control.reinit !== "undefined" && control.reinit !== null) {
control.actions.reinit = control.reinit;
}
el.rgltimer = new rgltimerClass(Tick, control.start, control.interval, control.stop,
control.step, control.value, control.rate, control.loop, control.actions);
for (var i=0; i < components.length; i++) {
switch(components[i]) {
case "Slider": addSlider(control.start, control.stop,
control.step, control.value);
break;
case "Label": addLabel(control.labels, control.start,
control.step, control.precision);
break;
default:
addButton(components[i], buttonLabels[i], control.pause);
}
}
el.rgltimer.Tick();
};
/**
* Apply all registered controls
* @param { Object } el - DOM element of the control
* @param { Object } x - List of actions to apply
* @param { boolean } [draw=true] - Whether to redraw after applying
*/
rglwidgetClass.prototype.applyControls = function(el, x, draw) {
var self = this, reinit = x.reinit, i, control, type;
for (i = 0; i < x.length; i++) {
control = x[i];
type = control.type;
self[type](el, control);
}
if (typeof reinit !== "undefined" && reinit !== null) {
reinit = [].concat(reinit);
for (i = 0; i < reinit.length; i++)
self.getObj(reinit[i]).initialized = false;
}
if (typeof draw === "undefined" || draw)
self.drawScene();
};
/**
* Handler for scene change
* @param { Object } message - What sort of scene change to do?
*/
rglwidgetClass.prototype.sceneChangeHandler = function(message) {
var self = document.getElementById(message.elementId).rglinstance,
objs = message.objects, mat = message.material,
root = message.rootSubscene,
initSubs = message.initSubscenes,
redraw = message.redrawScene,
skipRedraw = message.skipRedraw,
deletes, subs, allsubs = [], i,j;
if (typeof message.delete !== "undefined") {
deletes = [].concat(message.delete);
if (typeof message.delfromSubscenes !== "undefined")
subs = [].concat(message.delfromSubscenes);
else
subs = [];
for (i = 0; i < deletes.length; i++) {
for (j = 0; j < subs.length; j++) {
self.delFromSubscene(deletes[i], subs[j]);
}
delete self.scene.objects[deletes[i]];
}
}
if (typeof objs !== "undefined") {
Object.keys(objs).forEach(function(key){
key = parseInt(key, 10);
self.scene.objects[key] = objs[key];
self.initObj(key);
var obj = self.getObj(key),
subs = [].concat(obj.inSubscenes), k;
allsubs = allsubs.concat(subs);
for (k = 0; k < subs.length; k++)
self.addToSubscene(key, subs[k]);
});
}
if (typeof mat !== "undefined") {
self.scene.material = mat;
}
if (typeof root !== "undefined") {
self.scene.rootSubscene = root;
}
if (typeof initSubs !== "undefined")
allsubs = allsubs.concat(initSubs);
allsubs = self.unique(allsubs);
for (i = 0; i < allsubs.length; i++) {
self.initSubscene(allsubs[i]);
}
if (typeof skipRedraw !== "undefined") {
root = self.getObj(self.scene.rootSubscene);
root.par3d.skipRedraw = skipRedraw;
}
if (redraw)
self.drawScene();
};
/**
* Set mouse mode for a subscene
* @param { string } mode - name of mode
* @param { number } button - button number (1 to 3)
* @param { number } subscene - subscene id number
* @param { number } stayActive - if truthy, don't clear brush
*/
rglwidgetClass.prototype.setMouseMode = function(mode, button, subscene, stayActive) {
var sub = this.getObj(subscene),
which = ["left", "right", "middle"][button - 1];
if (!stayActive && sub.par3d.mouseMode[which] === "selecting")
this.clearBrush(null);
sub.par3d.mouseMode[which] = mode;
};
/**
* The class of an rgl timer object
* @class
*/
/**
* Construct an rgltimerClass object
* @constructor
* @param { function } Tick - action when timer fires
* @param { number } startTime - nominal start time in seconds
* @param { number } interval - seconds between updates
* @param { number } stopTime - nominal stop time in seconds
* @param { number } stepSize - nominal step size
* @param { number } value - current nominal time
* @param { number } rate - nominal units per second
* @param { string } loop - "none", "cycle" or "oscillate"
* @param { Object } actions - list of actions
*/
rgltimerClass = function(Tick, startTime, interval, stopTime, stepSize, value, rate, loop, actions) {
this.enabled = false;
this.timerId = 0;
/** nominal start time in seconds */
this.startTime = startTime;   
/** current nominal time */      
this.value = value;
/** seconds between updates */                 
this.interval = interval;
/** nominal stop time */           
this.stopTime = stopTime;
/** nominal step size */           
this.stepSize = stepSize;
/** nominal units per second */           
this.rate = rate;
/** "none", "cycle", or "oscillate" */                   
this.loop = loop;
/** real world start time */                   
this.realStart = undefined;
/** multiplier for fast-forward or reverse */         
this.multiplier = 1;                
this.actions = actions;
this.Tick = Tick;
};
/**
* Start playing timer object
*/
rgltimerClass.prototype.play = function() {
if (this.enabled) {
this.enabled = false;
window.clearInterval(this.timerId);
this.timerId = 0;
return;
}
var tick = function(self) {
var now = new Date();
self.value = self.multiplier*self.rate*(now - self.realStart)/1000 + self.startTime;
self.forceToRange();
if (typeof self.Tick !== "undefined") {
self.Tick(self.value);
}
};
this.realStart = new Date() - 1000*(this.value - this.startTime)/this.rate/this.multiplier;
this.timerId = window.setInterval(tick, 1000*this.interval, this);
this.enabled = true;
};
/**
* Force value into legal range
*/
rgltimerClass.prototype.forceToRange = function() {
if (this.value > this.stopTime + this.stepSize/2 || this.value < this.startTime - this.stepSize/2) {
if (!this.loop) {
this.reset();
} else {
var cycle = this.stopTime - this.startTime + this.stepSize,
newval = (this.value - this.startTime) % cycle + this.startTime;
if (newval < this.startTime) {
newval += cycle;
}
this.realStart += (this.value - newval)*1000/this.multiplier/this.rate;
this.value = newval;
}
}
};
/**
* Reset to start values
*/
rgltimerClass.prototype.reset = function() {
this.value = this.startTime;
this.newmultiplier(1);
if (typeof this.Tick !== "undefined") {
this.Tick(this.value);
}
if (this.enabled)
this.play();  /* really pause... */
if (typeof this.PlayButton !== "undefined")
this.PlayButton.value = "Play";
};
/**
* Increase the multiplier to play faster
*/
rgltimerClass.prototype.faster = function() {
this.newmultiplier(Math.SQRT2*this.multiplier);
};
/**
* Decrease the multiplier to play slower
*/
rgltimerClass.prototype.slower = function() {
this.newmultiplier(this.multiplier/Math.SQRT2);
};
/**
* Change sign of multiplier to reverse direction
*/
rgltimerClass.prototype.reverse = function() {
this.newmultiplier(-this.multiplier);
};
/**
* Set multiplier for play speed
* @param { number } newmult - new value
*/
rgltimerClass.prototype.newmultiplier = function(newmult) {
if (newmult != this.multiplier) {
this.realStart += 1000*(this.value - this.startTime)/this.rate*(1/this.multiplier - 1/newmult);
this.multiplier = newmult;
}
};
/**
* Take one step
*/
rgltimerClass.prototype.step = function() {
this.value += this.rate*this.multiplier;
this.forceToRange();
if (typeof this.Tick !== "undefined")
this.Tick(this.value);
};</script>

<div id="unnamed_chunk_5div" class="rglWebGL">

</div>

<script type="text/javascript">
var unnamed_chunk_5div = document.getElementById("unnamed_chunk_5div"),
unnamed_chunk_5rgl = new rglwidgetClass();
unnamed_chunk_5div.width = 673;
unnamed_chunk_5div.height = 481;
unnamed_chunk_5rgl.initialize(unnamed_chunk_5div,
{"material":{"color":"#000000","alpha":1,"lit":true,"ambient":"#000000","specular":"#FFFFFF","emission":"#000000","shininess":50,"smooth":true,"front":"filled","back":"filled","size":3,"lwd":1,"fog":false,"point_antialias":false,"line_antialias":false,"texture":null,"textype":"rgb","texmipmap":false,"texminfilter":"linear","texmagfilter":"linear","texenvmap":false,"depth_mask":true,"depth_test":"less","isTransparent":false,"polygon_offset":[0,0]},"rootSubscene":6,"objects":{"12":{"id":12,"type":"triangles","material":{},"vertices":[[1.54247,1.6458,-1.2074],[1.56725,1.7664,-1.296],[1.54254,1.65902,-1.20187],[1.54254,1.65902,-1.20187],[1.56725,1.7664,-1.296],[1.54254,1.87378,-1.20187],[1.54254,1.87378,-1.20187],[1.56725,1.7664,-1.296],[1.54247,1.887,-1.2074],[1.54247,1.887,-1.2074],[1.56725,1.7664,-1.296],[1.54249,1.89224,-1.22148],[1.54249,1.89224,-1.22148],[1.56725,1.7664,-1.296],[1.54249,1.89224,-1.37052],[1.54249,1.89224,-1.37052],[1.56725,1.7664,-1.296],[1.54247,1.887,-1.3846],[1.54247,1.887,-1.3846],[1.56725,1.7664,-1.296],[1.54254,1.87378,-1.39013],[1.54254,1.87378,-1.39013],[1.56725,1.7664,-1.296],[1.54254,1.65902,-1.39013],[1.54254,1.65902,-1.39013],[1.56725,1.7664,-1.296],[1.54247,1.6458,-1.3846],[1.54247,1.6458,-1.3846],[1.56725,1.7664,-1.296],[1.54249,1.64056,-1.37052],[1.54249,1.64056,-1.37052],[1.56725,1.7664,-1.296],[1.54249,1.64056,-1.22148],[1.54249,1.64056,-1.22148],[1.56725,1.7664,-1.296],[1.54247,1.6458,-1.2074],[1.54247,1.6458,1.2074],[1.54254,1.65902,1.20187],[1.56725,1.7664,1.296],[1.54254,1.65902,1.20187],[1.54254,1.87378,1.20187],[1.56725,1.7664,1.296],[1.54254,1.87378,1.20187],[1.54247,1.887,1.2074],[1.56725,1.7664,1.296],[1.54247,1.887,1.2074],[1.54249,1.89224,1.22148],[1.56725,1.7664,1.296],[1.54249,1.89224,1.22148],[1.54249,1.89224,1.37052],[1.56725,1.7664,1.296],[1.54249,1.89224,1.37052],[1.54247,1.887,1.3846],[1.56725,1.7664,1.296],[1.54247,1.887,1.3846],[1.54254,1.87378,1.39013],[1.56725,1.7664,1.296],[1.54254,1.87378,1.39013],[1.54254,1.65902,1.39013],[1.56725,1.7664,1.296],[1.54254,1.65902,1.39013],[1.54247,1.6458,1.3846],[1.56725,1.7664,1.296],[1.54247,1.6458,1.3846],[1.54249,1.64056,1.37052],[1.56725,1.7664,1.296],[1.54249,1.64056,1.37052],[1.54249,1.64056,1.22148],[1.56725,1.7664,1.296],[1.54249,1.64056,1.22148],[1.54247,1.6458,1.2074],[1.56725,1.7664,1.296],[-1.64774,0.807097,1.09],[-1.68734,0.753372,1.09],[-1.64774,0.807097,1.052],[-1.68734,0.753372,1.052],[-1.64774,0.807097,1.052],[-1.68734,0.753372,1.09],[-1.0236,0.753372,-1.1],[-1.01076,0.69,-1.1],[-1.0236,0.753372,-1.062],[-1.01076,0.69,-1.062],[-1.0236,0.753372,-1.062],[-1.01076,0.69,-1.1],[1.94156,0.56163,-1.1],[1.95888,0.414,-1.1],[1.94156,0.56163,-1.062],[1.95888,0.414,-1.062],[1.94156,0.56163,-1.062],[1.95888,0.414,-1.1],[-0.850175,0.414,1.09],[-1.01076,0.69,1.09],[-0.850175,0.414,1.052],[-1.01076,0.69,1.052],[-0.850175,0.414,1.052],[-1.01076,0.69,1.09],[2.29264,0.414,-1.062],[2.29264,0.414,1.052],[2.29264,0.4692,-1.062],[2.29264,0.4692,1.052],[2.29264,0.4692,-1.062],[2.29264,0.414,1.052],[-0.830148,0.4692,1.052],[0.824928,0.4692,1.052],[-0.850175,0.414,1.052],[0.857956,0.414,1.052],[-0.850175,0.414,1.052],[0.824928,0.4692,1.052],[-2.28344,0.4692,1.052],[-1.87155,0.4692,1.052],[-2.28344,0.414,1.052],[-1.85152,0.414,1.052],[-2.28344,0.414,1.052],[-1.87155,0.4692,1.052],[-1.01076,0.69,-1.062],[-0.997156,0.7452,-1.062],[-1.0236,0.753372,-1.062],[-1.01051,0.808572,-1.062],[-1.0236,0.753372,-1.062],[-0.997156,0.7452,-1.062],[1.97355,0.61683,-1.062],[1.99191,0.4692,-1.062],[1.97355,0.61683,-0.55094],[1.99191,0.4692,-0.55094],[1.97355,0.61683,-0.55094],[1.99191,0.4692,-1.062],[1.14505,0.966418,1.052],[1.03507,0.8772,1.052],[1.14505,0.966418,0.54094],[1.03507,0.8772,0.54094],[1.14505,0.966418,0.54094],[1.03507,0.8772,1.052],[-1.7008,0.808572,1.052],[-1.71527,0.7452,1.052],[-1.7008,0.808572,0.54094],[-1.71527,0.7452,0.54094],[-1.7008,0.808572,0.54094],[-1.71527,0.7452,1.052],[-1.10547,0.898194,-1.062],[-1.04855,0.862297,-1.062],[-1.10547,0.898194,-0.55094],[-1.04855,0.862297,-0.55094],[-1.10547,0.898194,-0.55094],[-1.04855,0.862297,-1.062],[-1.87155,0.4692,-0.55094],[-0.830148,0.4692,-0.55094],[-1.87155,0.4692,0.54094],[-0.830148,0.4692,0.54094],[-1.87155,0.4692,0.54094],[-0.830148,0.4692,-0.55094],[-2.20303,0.3634,-1.1054],[-2.18632,0.371258,-1.12347],[-2.18632,0.371258,1.12347],[-2.18632,0.371258,1.12347],[-2.20303,0.3634,1.1054],[-2.20303,0.3634,-1.1054],[-2.42292,0.422587,-1.11797],[-2.42329,0.417169,-1.11838],[-2.42329,0.417169,1.11838],[-2.42329,0.417169,1.11838],[-2.42292,0.422587,1.11797],[-2.42292,0.422587,-1.11797],[-2.19088,0.422587,1.11797],[-2.19051,0.417169,1.11838],[-2.19051,0.417169,-1.11838],[-2.19051,0.417169,-1.11838],[-2.19088,0.422587,-1.11797],[-2.19088,0.422587,1.11797],[-2.18289,0.49916,-1.12699],[-2.18289,0.49916,1.12699],[-2.18944,0.475849,-1.1196],[-2.18944,0.475849,1.1196],[-2.18944,0.475849,-1.1196],[-2.18289,0.49916,1.12699],[2.20202,0.40434,1.131],[2.20202,0.4899,1.131],[2.18933,0.40434,1.12529],[2.18933,0.4899,1.12529],[2.18933,0.40434,1.12529],[2.20202,0.4899,1.131],[2.18408,0.40434,-1.1115],[2.18408,0.4899,-1.1115],[2.18933,0.40434,-1.12529],[2.18933,0.4899,-1.12529],[2.18933,0.40434,-1.12529],[2.18408,0.4899,-1.1115],[2.33897,0.401345,-1.12529],[2.42102,0.494263,-1.12529],[2.34485,0.399977,-1.1115],[2.42627,0.494263,-1.1115],[2.34485,0.399977,-1.1115],[2.42102,0.494263,-1.12529],[2.18408,0.66838,-1.1115],[2.18933,0.681065,-1.1115],[2.18933,0.66838,-1.12529],[2.19109,0.679314,-1.12338],[2.18933,0.66838,-1.12529],[2.18933,0.681065,-1.1115],[2.20202,0.40434,1.131],[2.18933,0.40434,1.12529],[2.20202,0.391655,1.12529],[2.19109,0.393406,1.12338],[2.20202,0.391655,1.12529],[2.18933,0.40434,1.12529],[2.43497,0.508966,1.1115],[2.42972,0.508966,1.12529],[2.42627,0.494263,1.1115],[2.42102,0.494263,1.12529],[2.42627,0.494263,1.1115],[2.42972,0.508966,1.12529],[2.31518,0.3864,1.1115],[2.20202,0.3864,1.1115],[2.31518,0.3864,-1.1115],[2.20202,0.3864,-1.1115],[2.31518,0.3864,-1.1115],[2.20202,0.3864,1.1115],[2.3276,0.897,0.505],[2.3276,0.7406,0.505],[2.27884,0.897,0.505],[2.27884,0.7406,0.505],[2.27884,0.897,0.505],[2.3276,0.7406,0.505],[2.3276,0.7682,0.92],[2.3276,0.7682,0.645],[2.31564,0.7682,0.92],[2.31564,0.7682,0.645],[2.31564,0.7682,0.92],[2.3276,0.7682,0.645],[1.3271,0.4646,0.8625],[1.5341,0.4646,0.8625],[1.3271,0.4646,0.408501],[1.5341,0.4646,0.408501],[1.3271,0.4646,0.408501],[1.5341,0.4646,0.8625],[1.36074,0.2668,0.8625],[1.3271,0.3358,0.8625],[1.36074,0.2668,0.408501],[1.3271,0.3358,0.408501],[1.36074,0.2668,0.408501],[1.3271,0.3358,0.8625],[1.50046,0.306469,-0.253456],[1.36074,0.306469,-0.253456],[1.50046,0.313808,-0.287501],[1.36074,0.313808,-0.287501],[1.50046,0.313808,-0.287501],[1.36074,0.306469,-0.253456],[1.5341,0.504269,0.253456],[1.3271,0.504269,0.253456],[1.5341,0.511608,0.287501],[1.3271,0.511608,0.287501],[1.5341,0.511608,0.287501],[1.3271,0.504269,0.253456],[1.3271,0.4646,-0.408501],[1.3271,0.504269,-0.321545],[1.5341,0.4646,-0.408501],[1.5341,0.504269,-0.321545],[1.5341,0.4646,-0.408501],[1.3271,0.504269,-0.321545],[1.3271,0.4646,0.166501],[1.3271,0.504269,0.253456],[1.5341,0.4646,0.166501],[1.5341,0.504269,0.253456],[1.5341,0.4646,0.166501],[1.3271,0.504269,0.253456],[-1.4559,0.3358,0.166501],[-1.4559,0.4646,0.166501],[-1.4559,0.3358,-0.1665],[-1.4559,0.4646,-0.1665],[-1.4559,0.3358,-0.1665],[-1.4559,0.4646,0.166501],[-1.28254,0.2668,0.166501],[-1.42226,0.2668,0.166501],[-1.28254,0.2668,-0.1665],[-1.42226,0.2668,-0.1665],[-1.28254,0.2668,-0.1665],[-1.42226,0.2668,0.166501],[-1.4559,0.375469,-0.321545],[-1.4559,0.382808,-0.287501],[-1.4559,0.504269,-0.321545],[-1.4559,0.511608,-0.287501],[-1.4559,0.504269,-0.321545],[-1.4559,0.382808,-0.287501],[-1.2489,0.375469,0.253456],[-1.2489,0.504269,0.253456],[-1.2489,0.382808,0.287501],[-1.2489,0.511608,0.287501],[-1.2489,0.382808,0.287501],[-1.2489,0.504269,0.253456],[-1.2489,0.375469,-0.321545],[-1.28254,0.306469,-0.321545],[-1.2489,0.3358,-0.4085],[-1.28254,0.2668,-0.4085],[-1.2489,0.3358,-0.4085],[-1.28254,0.306469,-0.321545],[-1.28254,0.2668,-0.4085],[-1.28254,0.306469,-0.321545],[-1.42226,0.2668,-0.4085],[-1.42226,0.306469,-0.321545],[-1.42226,0.2668,-0.4085],[-1.28254,0.306469,-0.321545],[-2.3184,0.5244,0.3146],[-2.3184,0.5244,-0.3146],[-2.3184,0.7176,0.3146],[-2.3184,0.7176,-0.3146],[-2.3184,0.7176,0.3146],[-2.3184,0.5244,-0.3146],[-2.33266,0.5244,-0.3146],[-2.33266,0.5244,0.3146],[-2.33266,0.7176,-0.3146],[-2.33266,0.7176,0.3146],[-2.33266,0.7176,-0.3146],[-2.33266,0.5244,0.3146],[0.19872,2.4012,1.09],[0.175483,2.39658,1.09],[0.19872,2.4012,1.055],[0.175483,2.39658,1.055],[0.19872,2.4012,1.055],[0.175483,2.39658,1.09],[1.0028,1.2696,-1.055],[1.0028,1.2696,-1.1],[1.0028,1.242,-1.055],[1.0028,1.242,-1.1],[1.0028,1.242,-1.055],[1.0028,1.2696,-1.1],[0.965317,2.39658,-1.1],[0.965316,2.39658,-1.055],[0.94208,2.4012,-1.1],[0.94208,2.4012,-1.055],[0.94208,2.4012,-1.1],[0.965316,2.39658,-1.055],[1.53896,1.87858,-1.19766],[1.53896,1.65421,-1.19766],[1.54254,1.87378,-1.20187],[1.54254,1.65902,-1.20187],[1.54254,1.87378,-1.20187],[1.53896,1.65421,-1.19766],[1.53894,1.89212,-1.38836],[1.53895,1.89766,-1.37373],[1.54247,1.887,-1.3846],[1.54249,1.89224,-1.37052],[1.54247,1.887,-1.3846],[1.53895,1.89766,-1.37373],[1.5226,1.88048,-1.196],[1.53226,1.88048,-1.196],[1.5226,1.89414,-1.20215],[1.53226,1.89414,-1.20215],[1.5226,1.89414,-1.20215],[1.53226,1.88048,-1.196],[1.55986,1.60503,-1.2524],[1.55986,1.61123,-1.26866],[1.55986,1.57504,-1.285],[1.55986,1.58124,-1.30126],[1.55986,1.57504,-1.285],[1.55986,1.61123,-1.26866],[1.57421,1.79032,-1.2754],[1.57421,1.79032,-1.308],[1.55986,1.79032,-1.2754],[1.55986,1.79032,-1.308],[1.55986,1.79032,-1.2754],[1.57421,1.79032,-1.308],[1.5976,1.6192,-1.08653],[1.61,1.6192,-1.054],[2.29039,1.28184,-1.08653],[2.30279,1.28169,-1.054],[2.29039,1.28184,-1.08653],[1.61,1.6192,-1.054],[2.32089,1.2654,1.044],[2.3276,1.242,1.044],[2.32089,1.2654,-1.054],[2.3276,1.242,-1.054],[2.32089,1.2654,-1.054],[2.3276,1.242,1.044],[-2.3184,1.2696,1.044],[-2.3184,1.6192,1.044],[-2.3184,1.2696,0.6075],[-2.3184,1.6192,0.6075],[-2.3184,1.2696,0.6075],[-2.3184,1.6192,1.044],[-2.27608,1.6192,1.09],[-2.27608,2.53519,1.09],[-2.306,1.6192,1.07653],[-2.306,2.53519,1.07653],[-2.306,1.6192,1.07653],[-2.27608,2.53519,1.09],[-2.3184,0.7176,0.3146],[-2.33266,0.7176,0.3146],[-2.3184,0.7176,0.6075],[-2.33266,0.7176,0.6075],[-2.3184,0.7176,0.6075],[-2.33266,0.7176,0.3146],[-2.30271,2.57293,1.04572],[-2.29541,2.57293,1.06494],[-2.28569,2.58437,1.04203],[-2.28238,2.58437,1.05074],[-2.28569,2.58437,1.04203],[-2.29541,2.57293,1.06494],[-2.28569,2.58437,-1.05203],[-2.30271,2.57293,-1.05572],[-2.30271,2.57293,1.04572],[-2.30271,2.57293,1.04572],[-2.28569,2.58437,1.04203],[-2.28569,2.58437,-1.05203],[-2.27396,2.58437,-1.05157],[-2.28569,2.58437,-1.05203],[-2.28569,2.58437,1.04203],[-2.28569,2.58437,1.04203],[-2.27396,2.58437,1.04157],[-2.27396,2.58437,-1.05157],[1.61,1.6192,1.044],[1.334,2.53519,1.044],[1.5976,1.6192,1.07653],[1.3216,2.53519,1.07653],[1.5976,1.6192,1.07653],[1.334,2.53519,1.044],[1.53226,1.89414,1.20215],[1.53226,1.8998,1.217],[1.53894,1.89212,1.20364],[1.53895,1.89766,1.21827],[1.53894,1.89212,1.20364],[1.53226,1.8998,1.217],[1.53226,1.89414,1.38985],[1.53226,1.88048,1.396],[1.53894,1.89212,1.38836],[1.53896,1.87858,1.39434],[1.53894,1.89212,1.38836],[1.53226,1.88048,1.396],[1.53226,1.633,1.217],[1.53226,1.633,1.375],[1.5226,1.633,1.217],[1.5226,1.633,1.375],[1.5226,1.633,1.217],[1.53226,1.633,1.375],[1.55986,1.60503,1.09],[1.55986,1.57504,1.09],[1.55986,1.60503,1.2524],[1.55986,1.57504,1.285],[1.55986,1.60503,1.2524],[1.55986,1.57504,1.09],[-1.51855,0.8556,1.09],[-1.58847,0.842994,1.09],[-1.51855,0.8556,1.052],[-1.58847,0.842994,1.052],[-1.51855,0.8556,1.052],[-1.58847,0.842994,1.09],[-1.11491,0.842994,-1.1],[-1.06017,0.807097,-1.1],[-1.11491,0.842994,-1.062],[-1.06017,0.807097,-1.062],[-1.11491,0.842994,-1.062],[-1.06017,0.807097,-1.1],[1.80345,0.817334,-1.1],[1.88825,0.6992,-1.1],[1.80345,0.817334,-1.062],[1.88825,0.6992,-1.062],[1.80345,0.817334,-1.062],[1.88825,0.6992,-1.1],[1.0562,0.822,1.09],[0.857956,0.414,1.09],[1.0562,0.822,1.052],[0.857956,0.414,1.052],[1.0562,0.822,1.052],[0.857956,0.414,1.09],[1.95888,0.414,-1.062],[1.99191,0.4692,-1.062],[1.94156,0.56163,-1.062],[1.97355,0.61683,-1.062],[1.94156,0.56163,-1.062],[1.99191,0.4692,-1.062],[1.0562,0.822,1.052],[1.03507,0.8772,1.052],[1.15996,0.911218,1.052],[1.14505,0.966418,1.052],[1.15996,0.911218,1.052],[1.03507,0.8772,1.052],[-1.70125,0.69,1.052],[-1.71527,0.7452,1.052],[-1.68734,0.753372,1.052],[-1.7008,0.808572,1.052],[-1.68734,0.753372,1.052],[-1.71527,0.7452,1.052],[-1.10547,0.898194,-1.062],[-1.11491,0.842994,-1.062],[-1.04855,0.862297,-1.062],[-1.06017,0.807097,-1.062],[-1.04855,0.862297,-1.062],[-1.11491,0.842994,-1.062],[1.82716,0.872534,-1.062],[1.91704,0.7544,-1.062],[1.82716,0.872534,-0.55094],[1.91704,0.7544,-0.55094],[1.82716,0.872534,-0.55094],[1.91704,0.7544,-1.062],[1.42721,1.0396,1.052],[1.28081,1.02016,1.052],[1.42721,1.0396,0.54094],[1.28081,1.02016,0.54094],[1.42721,1.0396,0.54094],[1.28081,1.02016,1.052],[-1.59797,0.898194,1.052],[-1.65961,0.862297,1.052],[-1.59797,0.898194,0.54094],[-1.65961,0.862297,0.54094],[-1.59797,0.898194,0.54094],[-1.65961,0.862297,1.052],[-1.52526,0.9108,-1.062],[-1.17262,0.9108,-1.062],[-1.52526,0.9108,-0.55094],[-1.17262,0.9108,-0.55094],[-1.52526,0.9108,-0.55094],[-1.17262,0.9108,-1.062],[0.824928,0.4692,0.54094],[-0.830148,0.4692,0.54094],[0.824928,0.4692,-0.55094],[-0.830148,0.4692,-0.55094],[0.824928,0.4692,-0.55094],[-0.830148,0.4692,0.54094],[-2.18632,0.516542,-1.12347],[-2.42748,0.516542,-1.12347],[-2.20303,0.5244,-1.1054],[-2.41077,0.5244,-1.1054],[-2.20303,0.5244,-1.1054],[-2.42748,0.516542,-1.12347],[-2.41077,0.3634,1.1054],[-2.41077,0.3634,-1.1054],[-2.20303,0.3634,1.1054],[-2.20303,0.3634,-1.1054],[-2.20303,0.3634,1.1054],[-2.41077,0.3634,-1.1054],[-2.19051,0.470631,-1.11838],[-2.18944,0.475849,-1.1196],[-2.18944,0.475849,1.1196],[-2.18944,0.475849,1.1196],[-2.19051,0.470631,1.11838],[-2.19051,0.470631,-1.11838],[-2.43091,0.49916,1.12699],[-2.43091,0.49916,-1.12699],[-2.42436,0.475849,1.1196],[-2.42436,0.475849,-1.1196],[-2.42436,0.475849,1.1196],[-2.43091,0.49916,-1.12699],[2.20202,0.52578,1.131],[2.20202,0.66838,1.131],[2.18933,0.52578,1.12529],[2.18933,0.66838,1.12529],[2.18933,0.52578,1.12529],[2.20202,0.66838,1.131],[2.43275,0.681065,-1.1115],[2.42006,0.68632,-1.1115],[2.42006,0.68632,1.1115],[2.42006,0.68632,1.1115],[2.43275,0.681065,1.1115],[2.43275,0.681065,-1.1115],[2.20202,0.4899,-1.131],[2.20202,0.50784,-1.131],[2.40833,0.494263,-1.131],[2.41703,0.508966,-1.131],[2.40833,0.494263,-1.131],[2.20202,0.50784,-1.131],[2.20202,0.40434,-1.131],[2.20202,0.391655,-1.12529],[2.18933,0.40434,-1.12529],[2.19109,0.393406,-1.12338],[2.18933,0.40434,-1.12529],[2.20202,0.391655,-1.12529],[2.34485,0.399977,-1.1115],[2.33149,0.389956,-1.1115],[2.33897,0.401345,-1.12529],[2.33066,0.391774,-1.12382],[2.33897,0.401345,-1.12529],[2.33149,0.389956,-1.1115],[2.41703,0.508966,1.131],[2.42972,0.508966,1.12529],[2.42006,0.52578,1.131],[2.43275,0.52578,1.12529],[2.42006,0.52578,1.131],[2.42972,0.508966,1.12529],[2.18933,0.50784,1.12529],[2.20202,0.50784,1.131],[2.18933,0.52578,1.12529],[2.20202,0.52578,1.131],[2.18933,0.52578,1.12529],[2.20202,0.50784,1.131],[2.3276,0.7406,-0.505],[2.3276,0.897,-0.505],[2.27884,0.7406,-0.505],[2.27884,0.897,-0.505],[2.27884,0.7406,-0.505],[2.3276,0.897,-0.505],[2.3276,0.7682,-0.92],[2.3276,1.1914,-0.92],[2.31564,0.7682,-0.92],[2.31564,1.1914,-0.92],[2.31564,0.7682,-0.92],[2.3276,1.1914,-0.92],[1.3271,0.4646,-0.408501],[1.5341,0.4646,-0.408501],[1.3271,0.4646,-0.8625],[1.5341,0.4646,-0.8625],[1.3271,0.4646,-0.8625],[1.5341,0.4646,-0.408501],[1.36074,0.2668,-0.408501],[1.3271,0.3358,-0.408501],[1.36074,0.2668,-0.8625],[1.3271,0.3358,-0.8625],[1.36074,0.2668,-0.8625],[1.3271,0.3358,-0.408501],[1.36074,0.306469,-0.253456],[1.3271,0.375469,-0.253456],[1.36074,0.313808,-0.287501],[1.3271,0.382808,-0.287501],[1.36074,0.313808,-0.287501],[1.3271,0.375469,-0.253456],[1.3271,0.375469,0.253456],[1.3271,0.382808,0.287501],[1.3271,0.504269,0.253456],[1.3271,0.511608,0.287501],[1.3271,0.504269,0.253456],[1.3271,0.382808,0.287501],[1.5341,0.504269,0.253456],[1.5341,0.375469,0.253456],[1.5341,0.4646,0.166501],[1.5341,0.3358,0.166501],[1.5341,0.4646,0.166501],[1.5341,0.375469,0.253456],[1.3271,0.3358,-0.408501],[1.3271,0.375469,-0.321545],[1.3271,0.4646,-0.408501],[1.3271,0.504269,-0.321545],[1.3271,0.4646,-0.408501],[1.3271,0.375469,-0.321545],[-1.4559,0.4646,0.8625],[-1.2489,0.4646,0.8625],[-1.4559,0.4646,0.408501],[-1.2489,0.4646,0.408501],[-1.4559,0.4646,0.408501],[-1.2489,0.4646,0.8625],[-1.42226,0.2668,0.8625],[-1.4559,0.3358,0.8625],[-1.42226,0.2668,0.408501],[-1.4559,0.3358,0.408501],[-1.42226,0.2668,0.408501],[-1.4559,0.3358,0.8625],[-1.28254,0.306469,-0.253456],[-1.42226,0.306469,-0.253456],[-1.28254,0.313808,-0.287501],[-1.42226,0.313808,-0.287501],[-1.28254,0.313808,-0.287501],[-1.42226,0.306469,-0.253456],[-1.2489,0.504269,0.253456],[-1.4559,0.504269,0.253456],[-1.2489,0.511608,0.287501],[-1.4559,0.511608,0.287501],[-1.2489,0.511608,0.287501],[-1.4559,0.504269,0.253456],[-1.4559,0.4646,-0.4085],[-1.4559,0.504269,-0.321545],[-1.2489,0.4646,-0.4085],[-1.2489,0.504269,-0.321545],[-1.2489,0.4646,-0.4085],[-1.4559,0.504269,-0.321545],[-1.4559,0.4646,0.166501],[-1.4559,0.504269,0.253456],[-1.2489,0.4646,0.166501],[-1.2489,0.504269,0.253456],[-1.2489,0.4646,0.166501],[-1.4559,0.504269,0.253456],[-2.34692,0.5244,-0.6075],[-2.34692,0.5244,-0.3146],[-2.34692,0.7176,-0.6075],[-2.34692,0.7176,-0.3146],[-2.34692,0.7176,-0.6075],[-2.34692,0.5244,-0.3146],[-2.3184,0.7176,-0.3146],[-2.33266,0.7176,-0.3146],[-2.3184,0.7176,0.3146],[-2.33266,0.7176,0.3146],[-2.3184,0.7176,0.3146],[-2.33266,0.7176,-0.3146],[0.965317,2.39658,1.09],[0.94208,2.4012,1.09],[0.965316,2.39658,1.055],[0.94208,2.4012,1.055],[0.965316,2.39658,1.055],[0.94208,2.4012,1.09],[1.0028,1.6192,-1.055],[1.0028,1.6192,-1.1],[1.0028,1.2972,-1.055],[1.0028,1.2972,-1.1],[1.0028,1.2972,-1.055],[1.0028,1.6192,-1.1],[0.175483,2.39658,-1.055],[0.175483,2.39658,-1.1],[0.19872,2.4012,-1.055],[0.19872,2.4012,-1.1],[0.19872,2.4012,-1.055],[0.175483,2.39658,-1.1],[1.54254,1.65902,-1.39013],[1.54247,1.6458,-1.3846],[1.53896,1.65421,-1.39434],[1.53894,1.64068,-1.38836],[1.53896,1.65421,-1.39434],[1.54247,1.6458,-1.3846],[1.53894,1.64068,-1.20364],[1.53895,1.63514,-1.21827],[1.54247,1.6458,-1.2074],[1.54249,1.64056,-1.22148],[1.54247,1.6458,-1.2074],[1.53895,1.63514,-1.21827],[1.5226,1.8998,-1.375],[1.53226,1.8998,-1.375],[1.5226,1.89414,-1.38985],[1.53226,1.89414,-1.38985],[1.5226,1.89414,-1.38985],[1.53226,1.8998,-1.375],[1.58746,1.62619,-1.2754],[1.58746,1.61123,-1.26866],[1.58746,1.5962,-1.308],[1.58746,1.58124,-1.30126],[1.58746,1.5962,-1.308],[1.58746,1.61123,-1.26866],[1.58746,1.77707,-1.2754],[1.58746,1.62619,-1.2754],[1.58746,1.77707,-1.308],[1.58746,1.5962,-1.308],[1.58746,1.77707,-1.308],[1.58746,1.62619,-1.2754],[1.0028,1.2972,1.09],[1.0028,1.2696,1.09],[2.26047,1.28219,1.09],[2.27857,1.26554,1.09],[2.26047,1.28219,1.09],[1.0028,1.2696,1.09],[1.0028,1.6192,-1.1],[1.56768,1.6192,-1.1],[1.0028,1.2972,-1.1],[2.26047,1.28219,-1.1],[1.0028,1.2972,-1.1],[1.56768,1.6192,-1.1],[-2.306,1.6192,-1.08653],[-2.306,1.2696,-1.08653],[-2.3184,1.6192,-1.054],[-2.3184,1.2696,-1.054],[-2.3184,1.6192,-1.054],[-2.306,1.2696,-1.08653],[-2.28344,0.414,-1.062],[-2.28344,0.414,1.052],[-2.3184,0.414,-1.054],[-2.3184,0.414,1.044],[-2.3184,0.414,-1.054],[-2.28344,0.414,1.052],[-2.3184,0.7176,0.6075],[-2.33266,0.7176,0.6075],[-2.3184,0.5244,0.6075],[-2.33266,0.5244,0.6075],[-2.3184,0.5244,0.6075],[-2.33266,0.7176,0.6075],[-2.27438,2.58437,1.05435],[-2.27777,2.57293,1.0729],[1.29337,2.57293,1.0729],[1.29337,2.57293,1.0729],[1.28998,2.58437,1.05435],[-2.27438,2.58437,1.05435],[-2.30271,2.57293,-1.05572],[-2.31419,2.55597,-1.0582],[-2.31419,2.55597,1.0482],[-2.31419,2.55597,1.0482],[-2.30271,2.57293,1.04572],[-2.30271,2.57293,-1.05572],[-2.31419,2.55597,-1.0582],[-2.30271,2.57293,-1.05572],[-2.30419,2.55597,-1.08451],[-2.29541,2.57293,-1.07494],[-2.30419,2.55597,-1.08451],[-2.30271,2.57293,-1.05572],[1.5976,1.6192,1.07653],[1.3216,2.53519,1.07653],[1.56768,1.6192,1.09],[1.29168,2.53519,1.09],[1.56768,1.6192,1.09],[1.3216,2.53519,1.07653],[1.53226,1.89414,1.20215],[1.53894,1.89212,1.20364],[1.53226,1.88048,1.196],[1.53896,1.87858,1.19766],[1.53226,1.88048,1.196],[1.53894,1.89212,1.20364],[1.53226,1.63866,1.20215],[1.53226,1.65232,1.196],[1.53894,1.64068,1.20364],[1.53896,1.65421,1.19766],[1.53894,1.64068,1.20364],[1.53226,1.65232,1.196],[1.53226,1.8998,1.375],[1.53226,1.8998,1.217],[1.5226,1.8998,1.375],[1.5226,1.8998,1.217],[1.5226,1.8998,1.375],[1.53226,1.8998,1.217],[1.55986,1.57504,1.09],[1.58746,1.57504,1.09],[1.55986,1.57504,1.285],[1.58746,1.57504,1.285],[1.55986,1.57504,1.285],[1.58746,1.57504,1.09],[-1.58847,0.842994,1.09],[-1.64774,0.807097,1.09],[-1.58847,0.842994,1.052],[-1.64774,0.807097,1.052],[-1.58847,0.842994,1.052],[-1.64774,0.807097,1.09],[-1.06017,0.807097,-1.1],[-1.0236,0.753372,-1.1],[-1.06017,0.807097,-1.062],[-1.0236,0.753372,-1.062],[-1.06017,0.807097,-1.062],[-1.0236,0.753372,-1.1],[1.88825,0.6992,-1.1],[1.94156,0.56163,-1.1],[1.88825,0.6992,-1.062],[1.94156,0.56163,-1.062],[1.88825,0.6992,-1.062],[1.94156,0.56163,-1.1],[0.857956,0.414,1.09],[-0.850175,0.414,1.09],[0.857956,0.414,1.052],[-0.850175,0.414,1.052],[0.857956,0.414,1.052],[-0.850175,0.414,1.09],[2.29264,0.4692,-1.062],[1.99191,0.4692,-1.062],[2.29264,0.414,-1.062],[1.95888,0.414,-1.062],[2.29264,0.414,-1.062],[1.99191,0.4692,-1.062],[0.857956,0.414,1.052],[0.824928,0.4692,1.052],[1.0562,0.822,1.052],[1.03507,0.8772,1.052],[1.0562,0.822,1.052],[0.824928,0.4692,1.052],[-1.85152,0.414,1.052],[-1.87155,0.4692,1.052],[-1.70125,0.69,1.052],[-1.71527,0.7452,1.052],[-1.70125,0.69,1.052],[-1.87155,0.4692,1.052],[-1.04855,0.862297,-1.062],[-1.06017,0.807097,-1.062],[-1.01051,0.808572,-1.062],[-1.0236,0.753372,-1.062],[-1.01051,0.808572,-1.062],[-1.06017,0.807097,-1.062],[1.91704,0.7544,-1.062],[1.97355,0.61683,-1.062],[1.91704,0.7544,-0.55094],[1.97355,0.61683,-0.55094],[1.91704,0.7544,-0.55094],[1.97355,0.61683,-1.062],[1.28081,1.02016,1.052],[1.14505,0.966418,1.052],[1.28081,1.02016,0.54094],[1.14505,0.966418,0.54094],[1.28081,1.02016,0.54094],[1.14505,0.966418,1.052],[-1.65961,0.862297,1.052],[-1.7008,0.808572,1.052],[-1.65961,0.862297,0.54094],[-1.7008,0.808572,0.54094],[-1.65961,0.862297,0.54094],[-1.7008,0.808572,1.052],[-1.17262,0.9108,-1.062],[-1.10547,0.898194,-1.062],[-1.17262,0.9108,-0.55094],[-1.10547,0.898194,-0.55094],[-1.17262,0.9108,-0.55094],[-1.10547,0.898194,-1.062],[-2.24112,0.4692,-0.55094],[-1.87155,0.4692,-0.55094],[-2.24112,0.4692,0.54094],[-1.87155,0.4692,0.54094],[-2.24112,0.4692,0.54094],[-1.87155,0.4692,-0.55094],[-2.18289,0.49916,-1.12699],[-2.43091,0.49916,-1.12699],[-2.18632,0.516542,-1.12347],[-2.42748,0.516542,-1.12347],[-2.18632,0.516542,-1.12347],[-2.43091,0.49916,-1.12699],[-2.20303,0.5244,1.1054],[-2.20303,0.5244,-1.1054],[-2.41077,0.5244,1.1054],[-2.41077,0.5244,-1.1054],[-2.41077,0.5244,1.1054],[-2.20303,0.5244,-1.1054],[-2.19088,0.465213,-1.11797],[-2.19051,0.470631,-1.11838],[-2.19051,0.470631,1.11838],[-2.19051,0.470631,1.11838],[-2.19088,0.465213,1.11797],[-2.19088,0.465213,-1.11797],[-2.43091,0.38864,1.12699],[-2.42436,0.41195,1.1196],[-2.43091,0.38864,-1.12699],[-2.42436,0.41195,-1.1196],[-2.43091,0.38864,-1.12699],[-2.42436,0.41195,1.1196],[2.18933,0.52578,1.12529],[2.18933,0.66838,1.12529],[2.18408,0.52578,1.1115],[2.18408,0.66838,1.1115],[2.18408,0.52578,1.1115],[2.18933,0.66838,1.12529],[2.438,0.66838,-1.1115],[2.43275,0.681065,-1.1115],[2.43275,0.681065,1.1115],[2.43275,0.681065,1.1115],[2.438,0.66838,1.1115],[2.438,0.66838,-1.1115],[2.20202,0.50784,-1.131],[2.20202,0.52578,-1.131],[2.41703,0.508966,-1.131],[2.42006,0.52578,-1.131],[2.41703,0.508966,-1.131],[2.20202,0.52578,-1.131],[2.20202,0.66838,-1.131],[2.18933,0.66838,-1.12529],[2.20202,0.681065,-1.12529],[2.19109,0.679314,-1.12338],[2.20202,0.681065,-1.12529],[2.18933,0.66838,-1.12529],[2.31518,0.3864,-1.1115],[2.31799,0.391744,-1.12529],[2.33149,0.389956,-1.1115],[2.33066,0.391774,-1.12382],[2.33149,0.389956,-1.1115],[2.31799,0.391744,-1.12529],[2.40833,0.494263,1.131],[2.42102,0.494263,1.12529],[2.41703,0.508966,1.131],[2.42972,0.508966,1.12529],[2.41703,0.508966,1.131],[2.42102,0.494263,1.12529],[2.18933,0.50784,1.12529],[2.18933,0.52578,1.12529],[2.18408,0.50784,1.1115],[2.18408,0.52578,1.1115],[2.18408,0.50784,1.1115],[2.18933,0.52578,1.12529],[2.3276,0.897,-0.505],[2.3276,0.897,0.505],[2.27884,0.897,-0.505],[2.27884,0.897,0.505],[2.27884,0.897,-0.505],[2.3276,0.897,0.505],[2.3276,1.1914,-0.92],[2.3276,1.1914,-0.645],[2.31564,1.1914,-0.92],[2.31564,1.1914,-0.645],[2.31564,1.1914,-0.92],[2.3276,1.1914,-0.645],[1.3271,0.4646,0.166501],[1.5341,0.4646,0.166501],[1.3271,0.4646,-0.166501],[1.5341,0.4646,-0.166501],[1.3271,0.4646,-0.166501],[1.5341,0.4646,0.166501],[1.36074,0.2668,0.166501],[1.3271,0.3358,0.166501],[1.36074,0.2668,-0.166501],[1.3271,0.3358,-0.166501],[1.36074,0.2668,-0.166501],[1.3271,0.3358,0.166501],[1.3271,0.382808,-0.287501],[1.3271,0.375469,-0.321545],[1.36074,0.313808,-0.287501],[1.36074,0.306469,-0.321545],[1.36074,0.313808,-0.287501],[1.3271,0.375469,-0.321545],[1.3271,0.511608,0.287501],[1.3271,0.382808,0.287501],[1.3271,0.504269,0.321545],[1.3271,0.375469,0.321545],[1.3271,0.504269,0.321545],[1.3271,0.382808,0.287501],[1.5341,0.3358,0.408501],[1.5341,0.375469,0.321545],[1.5341,0.4646,0.408501],[1.5341,0.504269,0.321545],[1.5341,0.4646,0.408501],[1.5341,0.375469,0.321545],[1.3271,0.3358,-0.166501],[1.3271,0.4646,-0.166501],[1.3271,0.375469,-0.253456],[1.3271,0.504269,-0.253456],[1.3271,0.375469,-0.253456],[1.3271,0.4646,-0.166501],[-1.4559,0.3358,-0.4085],[-1.4559,0.4646,-0.4085],[-1.4559,0.3358,-0.8625],[-1.4559,0.4646,-0.8625],[-1.4559,0.3358,-0.8625],[-1.4559,0.4646,-0.4085],[-1.28254,0.2668,-0.4085],[-1.42226,0.2668,-0.4085],[-1.28254,0.2668,-0.8625],[-1.42226,0.2668,-0.8625],[-1.28254,0.2668,-0.8625],[-1.42226,0.2668,-0.4085],[-1.28254,0.313808,-0.287501],[-1.42226,0.313808,-0.287501],[-1.28254,0.306469,-0.321545],[-1.42226,0.306469,-0.321545],[-1.28254,0.306469,-0.321545],[-1.42226,0.313808,-0.287501],[-1.2489,0.511608,0.287501],[-1.4559,0.511608,0.287501],[-1.2489,0.504269,0.321545],[-1.4559,0.504269,0.321545],[-1.2489,0.504269,0.321545],[-1.4559,0.511608,0.287501],[-1.4559,0.504269,-0.253456],[-1.4559,0.4646,-0.1665],[-1.2489,0.504269,-0.253456],[-1.2489,0.4646,-0.1665],[-1.2489,0.504269,-0.253456],[-1.4559,0.4646,-0.1665],[-1.4559,0.504269,0.321545],[-1.4559,0.4646,0.408501],[-1.2489,0.504269,0.321545],[-1.2489,0.4646,0.408501],[-1.2489,0.504269,0.321545],[-1.4559,0.4646,0.408501],[-2.34692,0.5244,0.3146],[-2.34692,0.5244,0.6075],[-2.34692,0.7176,0.3146],[-2.34692,0.7176,0.6075],[-2.34692,0.7176,0.3146],[-2.34692,0.5244,0.6075],[-2.33266,0.5244,0.3146],[-2.34692,0.5244,0.3146],[-2.33266,0.7176,0.3146],[-2.34692,0.7176,0.3146],[-2.33266,0.7176,0.3146],[-2.34692,0.5244,0.3146],[0.94208,2.4012,1.09],[0.19872,2.4012,1.09],[0.94208,2.4012,1.055],[0.19872,2.4012,1.055],[0.94208,2.4012,1.055],[0.19872,2.4012,1.09],[1.0028,1.2972,-1.055],[1.0028,1.2972,-1.1],[1.0028,1.2696,-1.055],[1.0028,1.2696,-1.1],[1.0028,1.2696,-1.055],[1.0028,1.2972,-1.1],[0.94208,2.4012,-1.1],[0.94208,2.4012,-1.055],[0.19872,2.4012,-1.1],[0.19872,2.4012,-1.055],[0.19872,2.4012,-1.1],[0.94208,2.4012,-1.055],[1.53226,1.63866,-1.38985],[1.53226,1.65232,-1.396],[1.53894,1.64068,-1.38836],[1.53896,1.65421,-1.39434],[1.53894,1.64068,-1.38836],[1.53226,1.65232,-1.396],[1.53226,1.63866,-1.20215],[1.53226,1.633,-1.217],[1.53894,1.64068,-1.20364],[1.53895,1.63514,-1.21827],[1.53894,1.64068,-1.20364],[1.53226,1.633,-1.217],[1.5226,1.89414,-1.38985],[1.53226,1.89414,-1.38985],[1.5226,1.88048,-1.396],[1.53226,1.88048,-1.396],[1.5226,1.88048,-1.396],[1.53226,1.89414,-1.38985],[1.58746,1.60503,-1.2524],[1.58746,1.57504,-1.285],[1.58746,1.61123,-1.26866],[1.58746,1.58124,-1.30126],[1.58746,1.61123,-1.26866],[1.58746,1.57504,-1.285],[1.58358,1.78644,-1.2754],[1.57421,1.79032,-1.2754],[1.58746,1.77707,-1.2754],[1.58746,1.62619,-1.2754],[1.58746,1.77707,-1.2754],[1.57421,1.79032,-1.2754],[1.57421,1.79032,-1.2754],[1.55986,1.79032,-1.2754],[1.58746,1.62619,-1.2754],[1.55986,1.62619,-1.2754],[1.58746,1.62619,-1.2754],[1.55986,1.79032,-1.2754],[1.0028,1.6192,1.09],[1.0028,1.2972,1.09],[1.56768,1.6192,1.09],[2.26047,1.28219,1.09],[1.56768,1.6192,1.09],[1.0028,1.2972,1.09],[2.3276,1.1914,-0.645],[2.3276,1.15,-0.505],[2.3276,0.9936,-0.505],[2.3276,0.9936,0.505],[2.3276,1.15,0.505],[2.3276,1.1914,0.645],[2.3276,0.7682,-0.645],[2.3276,1.1914,-0.645],[2.3276,0.9936,-0.505],[2.3276,0.9936,0.505],[2.3276,1.1914,0.645],[2.3276,0.7682,0.645],[2.3276,0.7682,0.92],[2.3276,1.1914,0.92],[2.3276,1.242,1.044],[2.3276,0.897,0.505],[2.3276,0.9936,0.505],[2.3276,0.7682,0.645],[2.3276,0.7406,0.505],[2.3276,0.897,0.505],[2.3276,0.7682,0.645],[2.3276,0.7682,0.92],[2.3276,1.242,1.044],[2.3276,0.414,1.044],[2.3276,0.7682,0.645],[2.3276,0.7682,0.92],[2.3276,0.414,1.044],[2.3276,0.7682,0.645],[2.3276,0.414,1.044],[2.3276,0.414,-1.054],[2.3276,0.9936,-0.505],[2.3276,0.9936,0.505],[2.3276,0.897,0.505],[2.3276,0.9936,-0.505],[2.3276,0.897,0.505],[2.3276,0.897,-0.505],[2.3276,0.7682,-0.645],[2.3276,0.9936,-0.505],[2.3276,0.897,-0.505],[2.3276,0.7682,-0.645],[2.3276,0.897,-0.505],[2.3276,0.7406,-0.505],[2.3276,1.242,1.044],[2.3276,1.1914,0.92],[2.3276,1.1914,0.645],[2.3276,0.7682,-0.92],[2.3276,0.7682,-0.645],[2.3276,0.7406,-0.505],[2.3276,1.1914,0.645],[2.3276,1.15,0.505],[2.3276,1.15,-0.505],[2.3276,1.1914,0.645],[2.3276,1.15,-0.505],[2.3276,1.1914,-0.645],[2.3276,1.242,-1.054],[2.3276,1.242,1.044],[2.3276,1.1914,0.645],[2.3276,1.242,-1.054],[2.3276,1.1914,0.645],[2.3276,1.1914,-0.645],[2.3276,1.242,-1.054],[2.3276,1.1914,-0.645],[2.3276,1.1914,-0.92],[2.3276,1.242,-1.054],[2.3276,1.1914,-0.92],[2.3276,0.7682,-0.92],[2.3276,0.414,-1.054],[2.3276,1.242,-1.054],[2.3276,0.7682,-0.92],[2.3276,0.414,-1.054],[2.3276,0.7682,-0.92],[2.3276,0.7406,-0.505],[2.3276,0.414,-1.054],[2.3276,0.7406,-0.505],[2.3276,0.7406,0.505],[2.3276,0.414,-1.054],[2.3276,0.7406,0.505],[2.3276,0.7682,0.645],[0.138,1.6192,1.09],[-2.27608,1.6192,1.09],[0.138,1.2696,1.09],[-2.27608,1.2696,1.09],[0.138,1.2696,1.09],[-2.27608,1.6192,1.09],[-2.306,1.6192,1.07653],[-2.306,2.53519,1.07653],[-2.3184,1.6192,1.044],[-2.3184,2.53519,1.044],[-2.3184,1.6192,1.044],[-2.306,2.53519,1.07653],[-2.3184,0.7176,-0.6075],[-2.33266,0.7176,-0.6075],[-2.3184,0.7176,-0.3146],[-2.33266,0.7176,-0.3146],[-2.3184,0.7176,-0.3146],[-2.33266,0.7176,-0.6075],[-2.27777,2.57293,1.0729],[-2.27438,2.58437,1.05435],[-2.29541,2.57293,1.06494],[-2.28238,2.58437,1.05074],[-2.29541,2.57293,1.06494],[-2.27438,2.58437,1.05435],[-2.3184,2.53519,1.044],[-2.306,2.53519,1.07653],[-2.31419,2.55597,1.0482],[-2.30419,2.55597,1.07451],[-2.31419,2.55597,1.0482],[-2.306,2.53519,1.07653],[1.30129,2.58437,1.04203],[1.31831,2.57293,1.04572],[1.31831,2.57293,-1.05572],[1.31831,2.57293,-1.05572],[1.30129,2.58437,-1.05203],[1.30129,2.58437,1.04203],[1.29168,2.53519,1.09],[1.3216,2.53519,1.07653],[1.29565,2.55597,1.08541],[1.31979,2.55597,1.07451],[1.29565,2.55597,1.08541],[1.3216,2.53519,1.07653],[1.53894,1.89212,1.20364],[1.53895,1.89766,1.21827],[1.54247,1.887,1.2074],[1.54249,1.89224,1.22148],[1.54247,1.887,1.2074],[1.53895,1.89766,1.21827],[1.54254,1.87378,1.39013],[1.54247,1.887,1.3846],[1.53896,1.87858,1.39434],[1.53894,1.89212,1.38836],[1.53896,1.87858,1.39434],[1.54247,1.887,1.3846],[1.53226,1.65232,1.396],[1.53226,1.88048,1.396],[1.5226,1.65232,1.396],[1.5226,1.88048,1.396],[1.5226,1.65232,1.396],[1.53226,1.88048,1.396],[1.58746,1.60503,1.09],[1.58746,1.60503,1.2524],[1.58746,1.57504,1.09],[1.58746,1.57504,1.285],[1.58746,1.57504,1.09],[1.58746,1.60503,1.2524],[-1.58847,0.842994,-1.1],[-1.51855,0.8556,-1.1],[-1.58847,0.842994,-1.062],[-1.51855,0.8556,-1.062],[-1.58847,0.842994,-1.062],[-1.51855,0.8556,-1.1],[1.42614,0.9844,-1.1],[1.56425,0.964964,-1.1],[1.42614,0.9844,-1.062],[1.56425,0.964964,-1.062],[1.42614,0.9844,-1.062],[1.56425,0.964964,-1.1],[1.42614,0.9844,1.09],[1.28804,0.964964,1.09],[1.42614,0.9844,1.052],[1.28804,0.964964,1.052],[1.42614,0.9844,1.052],[1.28804,0.964964,1.09],[1.80345,0.817334,-1.062],[1.82716,0.872534,-1.062],[1.69294,0.907981,-1.062],[1.71001,0.963181,-1.062],[1.69294,0.907981,-1.062],[1.82716,0.872534,-1.062],[1.42721,1.0396,1.052],[1.5736,1.02016,1.052],[1.42614,0.9844,1.052],[1.56425,0.964964,1.052],[1.42614,0.9844,1.052],[1.5736,1.02016,1.052],[-1.52526,0.9108,1.052],[-1.51855,0.8556,1.052],[-1.59797,0.898194,1.052],[-1.58847,0.842994,1.052],[-1.59797,0.898194,1.052],[-1.51855,0.8556,1.052],[-1.52526,0.9108,-1.062],[-1.59797,0.898194,-1.062],[-1.51855,0.8556,-1.062],[-1.58847,0.842994,-1.062],[-1.51855,0.8556,-1.062],[-1.59797,0.898194,-1.062],[1.42721,1.0396,-1.062],[1.5736,1.02016,-1.062],[1.42721,1.0396,-0.55094],[1.5736,1.02016,-0.55094],[1.42721,1.0396,-0.55094],[1.5736,1.02016,-1.062],[1.82716,0.872534,1.052],[1.71001,0.963181,1.052],[1.82716,0.872534,0.54094],[1.71001,0.963181,0.54094],[1.82716,0.872534,0.54094],[1.71001,0.963181,1.052],[-1.10547,0.898194,1.052],[-1.17262,0.9108,1.052],[-1.10547,0.898194,0.54094],[-1.17262,0.9108,0.54094],[-1.10547,0.898194,0.54094],[-1.17262,0.9108,1.052],[-1.7008,0.808572,-1.062],[-1.65961,0.862297,-1.062],[-1.7008,0.808572,-0.55094],[-1.65961,0.862297,-0.55094],[-1.7008,0.808572,-0.55094],[-1.65961,0.862297,-1.062],[1.14505,0.966418,-1.062],[1.28081,1.02016,-1.062],[1.14505,0.966418,-0.55094],[1.28081,1.02016,-0.55094],[1.14505,0.966418,-0.55094],[1.28081,1.02016,-1.062],[-2.41077,0.5244,-1.1054],[-2.42748,0.516542,-1.12347],[-2.42748,0.516542,1.12347],[-2.42748,0.516542,1.12347],[-2.41077,0.5244,1.1054],[-2.41077,0.5244,-1.1054],[-2.18632,0.516542,1.12347],[-2.18289,0.49916,1.12699],[-2.18289,0.49916,-1.12699],[-2.18289,0.49916,-1.12699],[-2.18632,0.516542,-1.12347],[-2.18632,0.516542,1.12347],[-2.42436,0.41195,1.1196],[-2.18944,0.41195,1.1196],[-2.42329,0.417169,1.11838],[-2.19051,0.417169,1.11838],[-2.42329,0.417169,1.11838],[-2.18944,0.41195,1.1196],[-2.42292,0.422587,1.11797],[-2.19088,0.422587,1.11797],[-2.42292,0.465213,1.11797],[-2.19088,0.465213,1.11797],[-2.42292,0.465213,1.11797],[-2.19088,0.422587,1.11797],[2.20202,0.50784,1.131],[2.41703,0.508966,1.131],[2.20202,0.52578,1.131],[2.42006,0.52578,1.131],[2.20202,0.52578,1.131],[2.41703,0.508966,1.131],[2.18408,0.40434,-1.1115],[2.18933,0.391655,-1.1115],[2.18933,0.391655,1.1115],[2.18933,0.391655,1.1115],[2.18408,0.40434,1.1115],[2.18408,0.40434,-1.1115],[2.18933,0.52578,-1.12529],[2.18933,0.66838,-1.12529],[2.20202,0.52578,-1.131],[2.20202,0.66838,-1.131],[2.20202,0.52578,-1.131],[2.18933,0.66838,-1.12529],[2.31518,0.3864,-1.1115],[2.33149,0.389956,-1.1115],[2.33149,0.389956,1.1115],[2.33149,0.389956,1.1115],[2.31518,0.3864,1.1115],[2.31518,0.3864,-1.1115],[2.438,0.52578,-1.1115],[2.43497,0.508966,-1.1115],[2.43275,0.52578,-1.12529],[2.42972,0.508966,-1.12529],[2.43275,0.52578,-1.12529],[2.43497,0.508966,-1.1115],[2.42006,0.66838,1.131],[2.43275,0.66838,1.12529],[2.42006,0.681065,1.12529],[2.43099,0.679314,1.12338],[2.42006,0.681065,1.12529],[2.43275,0.66838,1.12529],[2.18408,0.50784,-1.1115],[2.18933,0.50784,-1.12529],[2.18408,0.4899,-1.1115],[2.18933,0.4899,-1.12529],[2.18408,0.4899,-1.1115],[2.18933,0.50784,-1.12529],[2.20202,0.40434,1.131],[2.32565,0.402713,1.131],[2.20202,0.4899,1.131],[2.40833,0.494263,1.131],[2.20202,0.4899,1.131],[2.32565,0.402713,1.131],[2.27884,1.15,0.505],[2.27884,0.9936,0.505],[2.27884,1.15,-0.505],[2.27884,0.9936,-0.505],[2.27884,1.15,-0.505],[2.27884,0.9936,0.505],[1.5341,0.3358,0.8625],[1.5341,0.4646,0.8625],[1.3271,0.3358,0.8625],[1.3271,0.4646,0.8625],[1.3271,0.3358,0.8625],[1.5341,0.4646,0.8625],[1.5341,0.3358,-0.408501],[1.50046,0.2668,-0.408501],[1.5341,0.3358,-0.8625],[1.50046,0.2668,-0.8625],[1.5341,0.3358,-0.8625],[1.50046,0.2668,-0.408501],[1.5341,0.382808,0.287501],[1.50046,0.313808,0.287501],[1.5341,0.375469,0.253456],[1.50046,0.306469,0.253456],[1.5341,0.375469,0.253456],[1.50046,0.313808,0.287501],[1.5341,0.511608,-0.287501],[1.5341,0.504269,-0.253456],[1.5341,0.382808,-0.287501],[1.5341,0.375469,-0.253456],[1.5341,0.382808,-0.287501],[1.5341,0.504269,-0.253456],[1.5341,0.3358,-0.166501],[1.5341,0.375469,-0.253456],[1.5341,0.4646,-0.166501],[1.5341,0.504269,-0.253456],[1.5341,0.4646,-0.166501],[1.5341,0.375469,-0.253456],[1.5341,0.375469,0.321545],[1.5341,0.3358,0.408501],[1.50046,0.306469,0.321545],[1.50046,0.2668,0.408501],[1.50046,0.306469,0.321545],[1.5341,0.3358,0.408501],[-1.4559,0.3358,-0.8625],[-1.2489,0.3358,-0.8625],[-1.42226,0.2668,-0.8625],[-1.28254,0.2668,-0.8625],[-1.42226,0.2668,-0.8625],[-1.2489,0.3358,-0.8625],[-1.2489,0.3358,0.8625],[-1.28254,0.2668,0.8625],[-1.2489,0.3358,0.408501],[-1.28254,0.2668,0.408501],[-1.2489,0.3358,0.408501],[-1.28254,0.2668,0.8625],[-1.42226,0.313808,0.287501],[-1.28254,0.313808,0.287501],[-1.42226,0.306469,0.321545],[-1.28254,0.306469,0.321545],[-1.42226,0.306469,0.321545],[-1.28254,0.313808,0.287501],[-1.2489,0.511608,-0.287501],[-1.4559,0.511608,-0.287501],[-1.2489,0.504269,-0.253456],[-1.4559,0.504269,-0.253456],[-1.2489,0.504269,-0.253456],[-1.4559,0.511608,-0.287501],[-1.28254,0.306469,0.321545],[-1.28254,0.2668,0.408501],[-1.42226,0.306469,0.321545],[-1.42226,0.2668,0.408501],[-1.42226,0.306469,0.321545],[-1.28254,0.2668,0.408501],[-1.42226,0.2668,0.408501],[-1.4559,0.3358,0.408501],[-1.42226,0.306469,0.321545],[-1.4559,0.375469,0.321545],[-1.42226,0.306469,0.321545],[-1.4559,0.3358,0.408501],[-2.2954,2.4104,0.6075],[-2.2954,2.4104,-0.6075],[-2.2954,1.6192,0.6075],[-2.2954,1.6192,-0.6075],[-2.2954,1.6192,0.6075],[-2.2954,2.4104,-0.6075],[-2.3184,0.5244,0.6075],[-2.33266,0.5244,0.6075],[-2.3184,0.5244,0.3146],[-2.33266,0.5244,0.3146],[-2.3184,0.5244,0.3146],[-2.33266,0.5244,0.6075],[1.0028,2.34048,1.09],[0.998178,2.36372,1.09],[1.0028,2.34048,1.055],[0.998178,2.36372,1.055],[1.0028,2.34048,1.055],[0.998178,2.36372,1.09],[0.138,0.4784,1.055],[0.138,0.4784,1.09],[0.7912,0.4784,1.055],[0.7912,0.4784,1.09],[0.7912,0.4784,1.055],[0.138,0.4784,1.09],[0.138,2.34048,-1.055],[0.138,2.34048,-1.1],[0.142622,2.36372,-1.055],[0.142622,2.36372,-1.1],[0.142622,2.36372,-1.055],[0.138,2.34048,-1.1],[1.53226,1.65232,-1.396],[1.53226,1.88048,-1.396],[1.53896,1.65421,-1.39434],[1.53896,1.87858,-1.39434],[1.53896,1.65421,-1.39434],[1.53226,1.88048,-1.396],[1.53226,1.633,-1.217],[1.53226,1.633,-1.375],[1.53895,1.63514,-1.21827],[1.53895,1.63514,-1.37373],[1.53895,1.63514,-1.21827],[1.53226,1.633,-1.375],[1.5226,1.63866,-1.38985],[1.53226,1.63866,-1.38985],[1.5226,1.633,-1.375],[1.53226,1.633,-1.375],[1.5226,1.633,-1.375],[1.53226,1.63866,-1.38985],[1.55986,1.58124,-1.30126],[1.55986,1.5962,-1.308],[1.58746,1.58124,-1.30126],[1.58746,1.5962,-1.308],[1.58746,1.58124,-1.30126],[1.55986,1.5962,-1.308],[1.58746,1.77707,-1.308],[1.58358,1.78644,-1.308],[1.58746,1.77707,-1.2754],[1.58358,1.78644,-1.2754],[1.58746,1.77707,-1.2754],[1.58358,1.78644,-1.308],[2.3152,1.242,1.07653],[2.30849,1.26544,1.07653],[2.28528,1.242,1.09],[2.27857,1.26554,1.09],[2.28528,1.242,1.09],[2.30849,1.26544,1.07653],[2.26047,1.28219,-1.1],[2.29039,1.28184,-1.08653],[2.27857,1.26554,-1.1],[2.30849,1.26544,-1.08653],[2.27857,1.26554,-1.1],[2.29039,1.28184,-1.08653],[-2.27608,1.2696,1.09],[-2.27608,1.6192,1.09],[-2.306,1.2696,1.07653],[-2.306,1.6192,1.07653],[-2.306,1.2696,1.07653],[-2.27608,1.6192,1.09],[2.29264,0.414,1.052],[2.29264,0.414,-1.062],[2.3276,0.414,1.044],[2.3276,0.414,-1.054],[2.3276,0.414,1.044],[2.29264,0.414,-1.062],[-2.3184,0.7176,-0.6075],[-2.3184,0.7176,-0.3146],[-2.3184,0.5244,-0.6075],[-2.3184,0.5244,-0.3146],[-2.3184,0.5244,-0.6075],[-2.3184,0.7176,-0.3146],[-2.27777,2.57293,1.0729],[-2.28005,2.55597,1.08541],[1.29565,2.55597,1.08541],[1.29565,2.55597,1.08541],[1.29337,2.57293,1.0729],[-2.27777,2.57293,1.0729],[-2.31419,2.55597,1.0482],[-2.31419,2.55597,-1.0582],[-2.3184,2.53519,1.044],[-2.3184,2.53519,-1.054],[-2.3184,2.53519,1.044],[-2.31419,2.55597,-1.0582],[-2.27608,2.53519,-1.1],[-2.306,2.53519,-1.08653],[-2.28005,2.55597,-1.09541],[-2.30419,2.55597,-1.08451],[-2.28005,2.55597,-1.09541],[-2.306,2.53519,-1.08653],[1.28956,2.58437,1.04157],[1.28998,2.58437,1.05435],[1.30129,2.58437,1.04203],[1.29798,2.58437,1.05074],[1.30129,2.58437,1.04203],[1.28998,2.58437,1.05435],[1.29168,2.53519,-1.1],[1.3216,2.53519,-1.08653],[1.56768,1.6192,-1.1],[1.5976,1.6192,-1.08653],[1.56768,1.6192,-1.1],[1.3216,2.53519,-1.08653],[1.53896,1.87858,1.19766],[1.54254,1.87378,1.20187],[1.53896,1.65421,1.19766],[1.54254,1.65902,1.20187],[1.53896,1.65421,1.19766],[1.54254,1.87378,1.20187],[1.53894,1.89212,1.38836],[1.54247,1.887,1.3846],[1.53895,1.89766,1.37373],[1.54249,1.89224,1.37052],[1.53895,1.89766,1.37373],[1.54247,1.887,1.3846],[1.5226,1.88048,1.196],[1.5226,1.89414,1.20215],[1.53226,1.88048,1.196],[1.53226,1.89414,1.20215],[1.53226,1.88048,1.196],[1.5226,1.89414,1.20215],[1.55986,1.60503,1.2524],[1.55986,1.57504,1.285],[1.55986,1.61123,1.26866],[1.55986,1.58124,1.30126],[1.55986,1.61123,1.26866],[1.55986,1.57504,1.285],[1.57421,1.79032,1.2754],[1.55986,1.79032,1.2754],[1.57421,1.79032,1.308],[1.55986,1.79032,1.308],[1.57421,1.79032,1.308],[1.55986,1.79032,1.2754],[-1.51855,0.8556,-1.1],[-1.17947,0.8556,-1.1],[-1.51855,0.8556,-1.062],[-1.17947,0.8556,-1.062],[-1.51855,0.8556,-1.062],[-1.17947,0.8556,-1.1],[1.56425,0.964964,-1.1],[1.69294,0.907981,-1.1],[1.56425,0.964964,-1.062],[1.69294,0.907981,-1.062],[1.56425,0.964964,-1.062],[1.69294,0.907981,-1.1],[1.28804,0.964964,1.09],[1.15996,0.911218,1.09],[1.28804,0.964964,1.052],[1.15996,0.911218,1.052],[1.28804,0.964964,1.052],[1.15996,0.911218,1.09],[1.88825,0.6992,-1.062],[1.91704,0.7544,-1.062],[1.80345,0.817334,-1.062],[1.82716,0.872534,-1.062],[1.80345,0.817334,-1.062],[1.91704,0.7544,-1.062],[1.42721,1.0396,1.052],[1.42614,0.9844,1.052],[1.28081,1.02016,1.052],[1.28804,0.964964,1.052],[1.28081,1.02016,1.052],[1.42614,0.9844,1.052],[-1.59797,0.898194,1.052],[-1.58847,0.842994,1.052],[-1.65961,0.862297,1.052],[-1.64774,0.807097,1.052],[-1.65961,0.862297,1.052],[-1.58847,0.842994,1.052],[-1.17262,0.9108,-1.062],[-1.52526,0.9108,-1.062],[-1.17947,0.8556,-1.062],[-1.51855,0.8556,-1.062],[-1.17947,0.8556,-1.062],[-1.52526,0.9108,-1.062],[1.5736,1.02016,-1.062],[1.71001,0.963181,-1.062],[1.5736,1.02016,-0.55094],[1.71001,0.963181,-0.55094],[1.5736,1.02016,-0.55094],[1.71001,0.963181,-1.062],[1.71001,0.963181,1.052],[1.5736,1.02016,1.052],[1.71001,0.963181,0.54094],[1.5736,1.02016,0.54094],[1.71001,0.963181,0.54094],[1.5736,1.02016,1.052],[-1.17262,0.9108,1.052],[-1.52526,0.9108,1.052],[-1.17262,0.9108,0.54094],[-1.52526,0.9108,0.54094],[-1.17262,0.9108,0.54094],[-1.52526,0.9108,1.052],[-1.65961,0.862297,-1.062],[-1.59797,0.898194,-1.062],[-1.65961,0.862297,-0.55094],[-1.59797,0.898194,-0.55094],[-1.65961,0.862297,-0.55094],[-1.59797,0.898194,-1.062],[1.28081,1.02016,-1.062],[1.42721,1.0396,-1.062],[1.28081,1.02016,-0.55094],[1.42721,1.0396,-0.55094],[1.28081,1.02016,-0.55094],[1.42721,1.0396,-1.062],[-2.41077,0.3634,1.1054],[-2.42748,0.371258,1.12347],[-2.42748,0.371258,-1.12347],[-2.42748,0.371258,-1.12347],[-2.41077,0.3634,-1.1054],[-2.41077,0.3634,1.1054],[-2.18632,0.371258,1.12347],[-2.42748,0.371258,1.12347],[-2.20303,0.3634,1.1054],[-2.41077,0.3634,1.1054],[-2.20303,0.3634,1.1054],[-2.42748,0.371258,1.12347],[-2.18944,0.475849,1.1196],[-2.42436,0.475849,1.1196],[-2.19051,0.470631,1.11838],[-2.42329,0.470631,1.11838],[-2.19051,0.470631,1.11838],[-2.42436,0.475849,1.1196],[-2.43091,0.49916,-1.12699],[-2.18289,0.49916,-1.12699],[-2.42436,0.475849,-1.1196],[-2.18944,0.475849,-1.1196],[-2.42436,0.475849,-1.1196],[-2.18289,0.49916,-1.12699],[2.438,0.66838,-1.1115],[2.438,0.52578,-1.1115],[2.43275,0.66838,-1.12529],[2.43275,0.52578,-1.12529],[2.43275,0.66838,-1.12529],[2.438,0.52578,-1.1115],[2.43497,0.508966,-1.1115],[2.438,0.52578,-1.1115],[2.438,0.52578,1.1115],[2.438,0.52578,1.1115],[2.43497,0.508966,1.1115],[2.43497,0.508966,-1.1115],[2.438,0.52578,1.1115],[2.438,0.66838,1.1115],[2.43275,0.52578,1.12529],[2.43275,0.66838,1.12529],[2.43275,0.52578,1.12529],[2.438,0.66838,1.1115],[2.20202,0.3864,-1.1115],[2.18933,0.391655,-1.1115],[2.20202,0.391655,-1.12529],[2.19109,0.393406,-1.12338],[2.20202,0.391655,-1.12529],[2.18933,0.391655,-1.1115],[2.42102,0.494263,-1.12529],[2.42972,0.508966,-1.12529],[2.42627,0.494263,-1.1115],[2.43497,0.508966,-1.1115],[2.42627,0.494263,-1.1115],[2.42972,0.508966,-1.12529],[2.438,0.66838,1.1115],[2.43275,0.681065,1.1115],[2.43275,0.66838,1.12529],[2.43099,0.679314,1.12338],[2.43275,0.66838,1.12529],[2.43275,0.681065,1.1115],[2.18408,0.50784,1.1115],[2.18408,0.4899,1.1115],[2.18933,0.50784,1.12529],[2.18933,0.4899,1.12529],[2.18933,0.50784,1.12529],[2.18408,0.4899,1.1115],[2.42006,0.66838,-1.131],[2.42006,0.52578,-1.131],[2.20202,0.66838,-1.131],[2.20202,0.52578,-1.131],[2.20202,0.66838,-1.131],[2.42006,0.52578,-1.131],[2.3276,1.1914,-0.645],[2.3276,0.7682,-0.645],[2.31564,1.1914,-0.645],[2.31564,0.7682,-0.645],[2.31564,1.1914,-0.645],[2.3276,0.7682,-0.645],[1.5341,0.3358,0.8625],[1.3271,0.3358,0.8625],[1.50046,0.2668,0.8625],[1.36074,0.2668,0.8625],[1.50046,0.2668,0.8625],[1.3271,0.3358,0.8625],[1.5341,0.3358,0.166501],[1.50046,0.2668,0.166501],[1.5341,0.3358,-0.166501],[1.50046,0.2668,-0.166501],[1.5341,0.3358,-0.166501],[1.50046,0.2668,0.166501],[1.36074,0.306469,0.253456],[1.50046,0.306469,0.253456],[1.36074,0.313808,0.287501],[1.50046,0.313808,0.287501],[1.36074,0.313808,0.287501],[1.50046,0.306469,0.253456],[1.5341,0.504269,-0.321545],[1.3271,0.504269,-0.321545],[1.5341,0.511608,-0.287501],[1.3271,0.511608,-0.287501],[1.5341,0.511608,-0.287501],[1.3271,0.504269,-0.321545],[1.50046,0.2668,0.166501],[1.50046,0.306469,0.253456],[1.36074,0.2668,0.166501],[1.36074,0.306469,0.253456],[1.36074,0.2668,0.166501],[1.50046,0.306469,0.253456],[1.36074,0.2668,0.166501],[1.36074,0.306469,0.253456],[1.3271,0.3358,0.166501],[1.3271,0.375469,0.253456],[1.3271,0.3358,0.166501],[1.36074,0.306469,0.253456],[-1.4559,0.4646,-0.4085],[-1.2489,0.4646,-0.4085],[-1.4559,0.4646,-0.8625],[-1.2489,0.4646,-0.8625],[-1.4559,0.4646,-0.8625],[-1.2489,0.4646,-0.4085],[-1.42226,0.2668,-0.4085],[-1.4559,0.3358,-0.4085],[-1.42226,0.2668,-0.8625],[-1.4559,0.3358,-0.8625],[-1.42226,0.2668,-0.8625],[-1.4559,0.3358,-0.4085],[-1.42226,0.306469,-0.253456],[-1.4559,0.375469,-0.253456],[-1.42226,0.313808,-0.287501],[-1.4559,0.382808,-0.287501],[-1.42226,0.313808,-0.287501],[-1.4559,0.375469,-0.253456],[-1.4559,0.375469,0.253456],[-1.4559,0.382808,0.287501],[-1.4559,0.504269,0.253456],[-1.4559,0.511608,0.287501],[-1.4559,0.504269,0.253456],[-1.4559,0.382808,0.287501],[-1.2489,0.504269,0.253456],[-1.2489,0.375469,0.253456],[-1.2489,0.4646,0.166501],[-1.2489,0.3358,0.166501],[-1.2489,0.4646,0.166501],[-1.2489,0.375469,0.253456],[-1.4559,0.3358,-0.4085],[-1.4559,0.375469,-0.321545],[-1.4559,0.4646,-0.4085],[-1.4559,0.504269,-0.321545],[-1.4559,0.4646,-0.4085],[-1.4559,0.375469,-0.321545],[-2.2954,1.6192,0.6075],[-2.2954,1.6192,-0.6075],[-2.2954,1.2696,0.6075],[-2.2954,1.2696,-0.6075],[-2.2954,1.2696,0.6075],[-2.2954,1.6192,-0.6075],[-2.33266,0.5244,0.6075],[-2.34692,0.5244,0.6075],[-2.33266,0.5244,0.3146],[-2.34692,0.5244,0.3146],[-2.33266,0.5244,0.3146],[-2.34692,0.5244,0.6075],[0.998178,2.36372,1.09],[0.985016,2.38342,1.09],[0.998178,2.36372,1.055],[0.985015,2.38342,1.055],[0.998178,2.36372,1.055],[0.985016,2.38342,1.09],[0.7912,0.4784,1.055],[0.7912,0.4784,1.09],[1.0028,0.8924,1.055],[1.0028,0.8924,1.09],[1.0028,0.8924,1.055],[0.7912,0.4784,1.09],[0.142622,2.36372,-1.055],[0.142622,2.36372,-1.1],[0.155784,2.38342,-1.055],[0.155784,2.38342,-1.1],[0.155784,2.38342,-1.055],[0.142622,2.36372,-1.1],[1.53894,1.64068,-1.38836],[1.54247,1.6458,-1.3846],[1.53895,1.63514,-1.37373],[1.54249,1.64056,-1.37052],[1.53895,1.63514,-1.37373],[1.54247,1.6458,-1.3846],[1.53895,1.89766,-1.37373],[1.53895,1.89766,-1.21827],[1.54249,1.89224,-1.37052],[1.54249,1.89224,-1.22148],[1.54249,1.89224,-1.37052],[1.53895,1.89766,-1.21827],[1.53226,1.65232,-1.196],[1.5226,1.65232,-1.196],[1.53226,1.63866,-1.20215],[1.5226,1.63866,-1.20215],[1.53226,1.63866,-1.20215],[1.5226,1.65232,-1.196],[1.55986,1.62619,-1.2754],[1.55986,1.61123,-1.26866],[1.58746,1.62619,-1.2754],[1.58746,1.61123,-1.26866],[1.58746,1.62619,-1.2754],[1.55986,1.61123,-1.26866],[1.58358,1.78644,-1.308],[1.57421,1.79032,-1.308],[1.58358,1.78644,-1.2754],[1.57421,1.79032,-1.2754],[1.58358,1.78644,-1.2754],[1.57421,1.79032,-1.308],[2.30849,1.26544,1.07653],[2.3152,1.242,1.07653],[2.32089,1.2654,1.044],[2.3276,1.242,1.044],[2.32089,1.2654,1.044],[2.3152,1.242,1.07653],[1.0028,1.2696,-1.1],[2.27857,1.26554,-1.1],[1.0028,1.242,-1.1],[2.28528,1.242,-1.1],[1.0028,1.242,-1.1],[2.27857,1.26554,-1.1],[-2.306,1.2696,1.07653],[-2.306,1.6192,1.07653],[-2.3184,1.2696,1.044],[-2.3184,1.6192,1.044],[-2.3184,1.2696,1.044],[-2.306,1.6192,1.07653],[-1.85152,0.414,1.09],[-2.27608,0.414,1.09],[-1.85152,0.414,1.052],[-2.28344,0.414,1.052],[-1.85152,0.414,1.052],[-2.27608,0.414,1.09],[-2.3184,0.7176,0.3146],[-2.3184,0.7176,0.6075],[-2.3184,0.5244,0.3146],[-2.3184,0.5244,0.6075],[-2.3184,0.5244,0.3146],[-2.3184,0.7176,0.6075],[1.29337,2.57293,-1.0829],[1.29565,2.55597,-1.09541],[-2.28005,2.55597,-1.09541],[-2.28005,2.55597,-1.09541],[-2.27777,2.57293,-1.0829],[1.29337,2.57293,-1.0829],[-2.27608,2.53519,1.09],[-2.28005,2.55597,1.08541],[-2.306,2.53519,1.07653],[-2.30419,2.55597,1.07451],[-2.306,2.53519,1.07653],[-2.28005,2.55597,1.08541],[-2.27777,2.57293,-1.0829],[-2.29541,2.57293,-1.07494],[-2.27438,2.58437,-1.06435],[-2.28238,2.58437,-1.06074],[-2.27438,2.58437,-1.06435],[-2.29541,2.57293,-1.07494],[0.998178,2.36372,1.09],[1.0028,2.34048,1.09],[1.0442,2.3483,1.09],[-2.27608,2.53519,1.09],[-2.27608,1.6192,1.09],[0.138,1.6192,1.09],[0.998178,2.36372,1.09],[1.0442,2.3483,1.09],[1.0879,2.392,1.09],[-2.27608,2.53519,1.09],[0.138,1.6192,1.09],[0.138,2.34048,1.09],[0.985016,2.38342,1.09],[0.998178,2.36372,1.09],[1.0879,2.392,1.09],[0.965317,2.39658,1.09],[0.985016,2.38342,1.09],[1.0879,2.392,1.09],[1.0442,1.6767,1.09],[1.0442,2.3483,1.09],[1.0028,2.34048,1.09],[1.0442,1.6767,1.09],[1.0028,2.34048,1.09],[1.0028,1.6192,1.09],[1.0879,1.633,1.09],[1.0442,1.6767,1.09],[1.0028,1.6192,1.09],[1.0879,1.633,1.09],[1.0028,1.6192,1.09],[1.56768,1.6192,1.09],[1.46502,1.633,1.09],[1.0879,1.633,1.09],[1.56768,1.6192,1.09],[1.53021,1.6919,1.09],[1.46502,1.633,1.09],[1.56768,1.6192,1.09],[0.94208,2.4012,1.09],[0.965317,2.39658,1.09],[1.0879,2.392,1.09],[0.94208,2.4012,1.09],[1.0879,2.392,1.09],[1.2489,2.392,1.09],[1.53021,1.6919,1.09],[1.56768,1.6192,1.09],[1.29168,2.53519,1.09],[1.30664,2.35062,1.09],[1.53021,1.6919,1.09],[1.29168,2.53519,1.09],[1.2489,2.392,1.09],[1.30664,2.35062,1.09],[1.29168,2.53519,1.09],[0.94208,2.4012,1.09],[1.2489,2.392,1.09],[1.29168,2.53519,1.09],[0.19872,2.4012,1.09],[0.94208,2.4012,1.09],[1.29168,2.53519,1.09],[0.19872,2.4012,1.09],[1.29168,2.53519,1.09],[-2.27608,2.53519,1.09],[-2.27608,2.53519,1.09],[0.138,2.34048,1.09],[0.142622,2.36372,1.09],[-2.27608,2.53519,1.09],[0.142622,2.36372,1.09],[0.155784,2.38342,1.09],[-2.27608,2.53519,1.09],[0.155784,2.38342,1.09],[0.175483,2.39658,1.09],[-2.27608,2.53519,1.09],[0.175483,2.39658,1.09],[0.19872,2.4012,1.09],[1.29168,2.53519,-1.1],[1.29565,2.55597,-1.09541],[1.3216,2.53519,-1.08653],[1.31979,2.55597,-1.08451],[1.3216,2.53519,-1.08653],[1.29565,2.55597,-1.09541],[1.53226,1.88048,1.196],[1.53896,1.87858,1.19766],[1.53226,1.65232,1.196],[1.53896,1.65421,1.19766],[1.53226,1.65232,1.196],[1.53896,1.87858,1.19766],[1.53226,1.89414,1.38985],[1.53894,1.89212,1.38836],[1.53226,1.8998,1.375],[1.53895,1.89766,1.37373],[1.53226,1.8998,1.375],[1.53894,1.89212,1.38836],[1.5226,1.89414,1.20215],[1.5226,1.8998,1.217],[1.53226,1.89414,1.20215],[1.53226,1.8998,1.217],[1.53226,1.89414,1.20215],[1.5226,1.8998,1.217],[1.55986,1.58124,1.30126],[1.55986,1.5962,1.308],[1.55986,1.61123,1.26866],[1.55986,1.62619,1.2754],[1.55986,1.61123,1.26866],[1.55986,1.5962,1.308],[-1.64774,0.807097,-1.1],[-1.58847,0.842994,-1.1],[-1.64774,0.807097,-1.062],[-1.58847,0.842994,-1.062],[-1.64774,0.807097,-1.062],[-1.58847,0.842994,-1.1],[1.28804,0.964964,-1.1],[1.42614,0.9844,-1.1],[1.28804,0.964964,-1.062],[1.42614,0.9844,-1.062],[1.28804,0.964964,-1.062],[1.42614,0.9844,-1.1],[1.56425,0.964964,1.09],[1.42614,0.9844,1.09],[1.56425,0.964964,1.052],[1.42614,0.9844,1.052],[1.56425,0.964964,1.052],[1.42614,0.9844,1.09],[1.69294,0.907981,-1.062],[1.71001,0.963181,-1.062],[1.56425,0.964964,-1.062],[1.5736,1.02016,-1.062],[1.56425,0.964964,-1.062],[1.71001,0.963181,-1.062],[1.69294,0.907981,1.052],[1.56425,0.964964,1.052],[1.71001,0.963181,1.052],[1.5736,1.02016,1.052],[1.71001,0.963181,1.052],[1.56425,0.964964,1.052],[-1.52526,0.9108,1.052],[-1.17262,0.9108,1.052],[-1.51855,0.8556,1.052],[-1.17947,0.8556,1.052],[-1.51855,0.8556,1.052],[-1.17262,0.9108,1.052],[-1.64774,0.807097,-1.062],[-1.58847,0.842994,-1.062],[-1.65961,0.862297,-1.062],[-1.59797,0.898194,-1.062],[-1.65961,0.862297,-1.062],[-1.58847,0.842994,-1.062],[1.42721,1.0396,-1.062],[1.28081,1.02016,-1.062],[1.42614,0.9844,-1.062],[1.28804,0.964964,-1.062],[1.42614,0.9844,-1.062],[1.28081,1.02016,-1.062],[1.91704,0.7544,1.052],[1.82716,0.872534,1.052],[1.91704,0.7544,0.54094],[1.82716,0.872534,0.54094],[1.91704,0.7544,0.54094],[1.82716,0.872534,1.052],[-1.04855,0.862297,1.052],[-1.10547,0.898194,1.052],[-1.04855,0.862297,0.54094],[-1.10547,0.898194,0.54094],[-1.04855,0.862297,0.54094],[-1.10547,0.898194,1.052],[-1.71527,0.7452,-1.062],[-1.7008,0.808572,-1.062],[-1.71527,0.7452,-0.55094],[-1.7008,0.808572,-0.55094],[-1.71527,0.7452,-0.55094],[-1.7008,0.808572,-1.062],[1.03507,0.8772,-1.062],[1.14505,0.966418,-1.062],[1.03507,0.8772,-0.55094],[1.14505,0.966418,-0.55094],[1.03507,0.8772,-0.55094],[1.14505,0.966418,-1.062],[-2.42748,0.516542,-1.12347],[-2.43091,0.49916,-1.12699],[-2.43091,0.49916,1.12699],[-2.43091,0.49916,1.12699],[-2.42748,0.516542,1.12347],[-2.42748,0.516542,-1.12347],[-2.20303,0.5244,1.1054],[-2.18632,0.516542,1.12347],[-2.18632,0.516542,-1.12347],[-2.18632,0.516542,-1.12347],[-2.20303,0.5244,-1.1054],[-2.20303,0.5244,1.1054],[-2.42329,0.417169,1.11838],[-2.19051,0.417169,1.11838],[-2.42292,0.422587,1.11797],[-2.19088,0.422587,1.11797],[-2.42292,0.422587,1.11797],[-2.19051,0.417169,1.11838],[-2.19088,0.422587,1.11797],[-2.19088,0.422587,-1.11797],[-2.19088,0.465213,1.11797],[-2.19088,0.465213,-1.11797],[-2.19088,0.465213,1.11797],[-2.19088,0.422587,-1.11797],[2.20202,0.4899,1.131],[2.40833,0.494263,1.131],[2.20202,0.50784,1.131],[2.41703,0.508966,1.131],[2.20202,0.50784,1.131],[2.40833,0.494263,1.131],[2.18933,0.391655,-1.1115],[2.20202,0.3864,-1.1115],[2.20202,0.3864,1.1115],[2.20202,0.3864,1.1115],[2.18933,0.391655,1.1115],[2.18933,0.391655,-1.1115],[2.18408,0.52578,-1.1115],[2.18408,0.66838,-1.1115],[2.18933,0.52578,-1.12529],[2.18933,0.66838,-1.12529],[2.18933,0.52578,-1.12529],[2.18408,0.66838,-1.1115],[2.33149,0.389956,-1.1115],[2.34485,0.399977,-1.1115],[2.34485,0.399977,1.1115],[2.34485,0.399977,1.1115],[2.33149,0.389956,1.1115],[2.33149,0.389956,-1.1115],[2.43275,0.52578,-1.12529],[2.42972,0.508966,-1.12529],[2.42006,0.52578,-1.131],[2.41703,0.508966,-1.131],[2.42006,0.52578,-1.131],[2.42972,0.508966,-1.12529],[2.42006,0.68632,1.1115],[2.42006,0.681065,1.12529],[2.43275,0.681065,1.1115],[2.43099,0.679314,1.12338],[2.43275,0.681065,1.1115],[2.42006,0.681065,1.12529],[2.20202,0.50784,-1.131],[2.20202,0.4899,-1.131],[2.18933,0.50784,-1.12529],[2.18933,0.4899,-1.12529],[2.18933,0.50784,-1.12529],[2.20202,0.4899,-1.131],[2.20202,0.52578,1.131],[2.42006,0.52578,1.131],[2.20202,0.66838,1.131],[2.42006,0.66838,1.131],[2.20202,0.66838,1.131],[2.42006,0.52578,1.131],[2.27884,0.897,0.505],[2.27884,0.7406,0.505],[2.27884,0.897,-0.505],[2.27884,0.7406,-0.505],[2.27884,0.897,-0.505],[2.27884,0.7406,0.505],[1.5341,0.4646,-0.8625],[1.5341,0.3358,-0.8625],[1.3271,0.4646,-0.8625],[1.3271,0.3358,-0.8625],[1.3271,0.4646,-0.8625],[1.5341,0.3358,-0.8625],[1.5341,0.4646,0.8625],[1.5341,0.3358,0.8625],[1.5341,0.4646,0.408501],[1.5341,0.3358,0.408501],[1.5341,0.4646,0.408501],[1.5341,0.3358,0.8625],[1.50046,0.306469,0.321545],[1.50046,0.313808,0.287501],[1.5341,0.375469,0.321545],[1.5341,0.382808,0.287501],[1.5341,0.375469,0.321545],[1.50046,0.313808,0.287501],[1.5341,0.375469,-0.321545],[1.5341,0.504269,-0.321545],[1.5341,0.382808,-0.287501],[1.5341,0.511608,-0.287501],[1.5341,0.382808,-0.287501],[1.5341,0.504269,-0.321545],[1.5341,0.504269,-0.321545],[1.5341,0.375469,-0.321545],[1.5341,0.4646,-0.408501],[1.5341,0.3358,-0.408501],[1.5341,0.4646,-0.408501],[1.5341,0.375469,-0.321545],[1.5341,0.375469,0.253456],[1.50046,0.306469,0.253456],[1.5341,0.3358,0.166501],[1.50046,0.2668,0.166501],[1.5341,0.3358,0.166501],[1.50046,0.306469,0.253456],[-1.2489,0.3358,0.8625],[-1.4559,0.3358,0.8625],[-1.28254,0.2668,0.8625],[-1.42226,0.2668,0.8625],[-1.28254,0.2668,0.8625],[-1.4559,0.3358,0.8625],[-1.2489,0.3358,0.166501],[-1.28254,0.2668,0.166501],[-1.2489,0.3358,-0.1665],[-1.28254,0.2668,-0.1665],[-1.2489,0.3358,-0.1665],[-1.28254,0.2668,0.166501],[-1.42226,0.306469,0.253456],[-1.28254,0.306469,0.253456],[-1.42226,0.313808,0.287501],[-1.28254,0.313808,0.287501],[-1.42226,0.313808,0.287501],[-1.28254,0.306469,0.253456],[-1.2489,0.504269,-0.321545],[-1.4559,0.504269,-0.321545],[-1.2489,0.511608,-0.287501],[-1.4559,0.511608,-0.287501],[-1.2489,0.511608,-0.287501],[-1.4559,0.504269,-0.321545],[-1.28254,0.2668,0.166501],[-1.28254,0.306469,0.253456],[-1.42226,0.2668,0.166501],[-1.42226,0.306469,0.253456],[-1.42226,0.2668,0.166501],[-1.28254,0.306469,0.253456],[-1.42226,0.2668,0.166501],[-1.42226,0.306469,0.253456],[-1.4559,0.3358,0.166501],[-1.4559,0.375469,0.253456],[-1.4559,0.3358,0.166501],[-1.42226,0.306469,0.253456],[-2.3184,1.6192,-0.6075],[-2.3184,1.2696,-0.6075],[-2.2954,1.6192,-0.6075],[-2.2954,1.2696,-0.6075],[-2.2954,1.6192,-0.6075],[-2.3184,1.2696,-0.6075],[-2.33266,0.5244,-0.3146],[-2.34692,0.5244,-0.3146],[-2.33266,0.5244,-0.6075],[-2.34692,0.5244,-0.6075],[-2.33266,0.5244,-0.6075],[-2.34692,0.5244,-0.3146],[1.0028,1.6192,1.055],[1.0028,1.6192,1.09],[1.0028,2.34048,1.055],[1.0028,2.34048,1.09],[1.0028,2.34048,1.055],[1.0028,1.6192,1.09],[0.138,1.2696,1.055],[0.138,1.2696,1.09],[0.138,0.4784,1.055],[0.138,0.4784,1.09],[0.138,0.4784,1.055],[0.138,1.2696,1.09],[0.138,1.6192,-1.055],[0.138,1.6192,-1.1],[0.138,2.34048,-1.055],[0.138,2.34048,-1.1],[0.138,2.34048,-1.055],[0.138,1.6192,-1.1],[1.53896,1.65421,-1.39434],[1.53896,1.87858,-1.39434],[1.54254,1.65902,-1.39013],[1.54254,1.87378,-1.39013],[1.54254,1.65902,-1.39013],[1.53896,1.87858,-1.39434],[1.53895,1.63514,-1.21827],[1.53895,1.63514,-1.37373],[1.54249,1.64056,-1.22148],[1.54249,1.64056,-1.37052],[1.54249,1.64056,-1.22148],[1.53895,1.63514,-1.37373],[1.5226,1.65232,-1.396],[1.53226,1.65232,-1.396],[1.5226,1.63866,-1.38985],[1.53226,1.63866,-1.38985],[1.5226,1.63866,-1.38985],[1.53226,1.65232,-1.396],[1.55986,1.57504,-1.285],[1.55986,1.58124,-1.30126],[1.58746,1.57504,-1.285],[1.58746,1.58124,-1.30126],[1.58746,1.57504,-1.285],[1.55986,1.58124,-1.30126],[1.58746,1.60503,-1.09],[1.58746,1.60503,-1.2524],[1.55986,1.60503,-1.09],[1.55986,1.60503,-1.2524],[1.55986,1.60503,-1.09],[1.58746,1.60503,-1.2524],[1.5976,1.6192,1.07653],[2.29039,1.28184,1.07653],[1.61,1.6192,1.044],[2.30279,1.28169,1.044],[1.61,1.6192,1.044],[2.29039,1.28184,1.07653],[2.32089,1.2654,-1.054],[2.30849,1.26544,-1.08653],[2.30279,1.28169,-1.054],[2.29039,1.28184,-1.08653],[2.30279,1.28169,-1.054],[2.30849,1.26544,-1.08653],[-2.306,0.414,1.07653],[-2.306,1.2696,1.07653],[-2.3184,0.414,1.044],[-2.3184,1.2696,1.044],[-2.3184,0.414,1.044],[-2.306,1.2696,1.07653],[1.95888,0.414,-1.1],[2.28528,0.414,-1.1],[1.95888,0.414,-1.062],[2.29264,0.414,-1.062],[1.95888,0.414,-1.062],[2.28528,0.414,-1.1],[-2.3184,0.7176,0.6075],[-2.3184,1.2696,0.6075],[-2.2954,0.7176,0.6075],[-2.2954,1.2696,0.6075],[-2.2954,0.7176,0.6075],[-2.3184,1.2696,0.6075],[-2.27396,2.58437,1.04157],[1.28956,2.58437,1.04157],[-2.27396,2.58437,-1.05157],[1.28956,2.58437,-1.05157],[-2.27396,2.58437,-1.05157],[1.28956,2.58437,1.04157],[1.31831,2.57293,1.04572],[1.30129,2.58437,1.04203],[1.31101,2.57293,1.06494],[1.29798,2.58437,1.05074],[1.31101,2.57293,1.06494],[1.30129,2.58437,1.04203],[1.31979,2.55597,-1.08451],[1.29565,2.55597,-1.09541],[1.31101,2.57293,-1.07494],[1.29337,2.57293,-1.0829],[1.31101,2.57293,-1.07494],[1.29565,2.55597,-1.09541],[-2.27396,2.58437,-1.05157],[-2.27438,2.58437,-1.06435],[-2.28569,2.58437,-1.05203],[-2.28238,2.58437,-1.06074],[-2.28569,2.58437,-1.05203],[-2.27438,2.58437,-1.06435],[-2.27608,2.53519,-1.1],[-2.28005,2.55597,-1.09541],[1.29565,2.55597,-1.09541],[1.29565,2.55597,-1.09541],[1.29168,2.53519,-1.1],[-2.27608,2.53519,-1.1],[1.53226,1.63866,1.38985],[1.53894,1.64068,1.38836],[1.53226,1.65232,1.396],[1.53896,1.65421,1.39434],[1.53226,1.65232,1.396],[1.53894,1.64068,1.38836],[1.53226,1.63866,1.20215],[1.53894,1.64068,1.20364],[1.53226,1.633,1.217],[1.53895,1.63514,1.21827],[1.53226,1.633,1.217],[1.53894,1.64068,1.20364],[1.5226,1.89414,1.38985],[1.5226,1.88048,1.396],[1.53226,1.89414,1.38985],[1.53226,1.88048,1.396],[1.53226,1.89414,1.38985],[1.5226,1.88048,1.396],[1.58746,1.58124,1.30126],[1.58746,1.57504,1.285],[1.58746,1.61123,1.26866],[1.58746,1.60503,1.2524],[1.58746,1.61123,1.26866],[1.58746,1.57504,1.285],[1.58358,1.78644,1.2754],[1.58746,1.77707,1.2754],[1.57421,1.79032,1.2754],[1.55986,1.79032,1.2754],[1.57421,1.79032,1.2754],[1.58746,1.77707,1.2754],[1.58746,1.77707,1.2754],[1.58746,1.62619,1.2754],[1.55986,1.79032,1.2754],[1.55986,1.62619,1.2754],[1.55986,1.79032,1.2754],[1.58746,1.62619,1.2754],[-1.70125,0.69,1.09],[-1.85152,0.414,1.09],[-1.70125,0.69,1.052],[-1.85152,0.414,1.052],[-1.70125,0.69,1.052],[-1.85152,0.414,1.09],[-0.850175,0.414,-1.1],[0.857956,0.414,-1.1],[-0.850175,0.414,-1.062],[0.857956,0.414,-1.062],[-0.850175,0.414,-1.062],[0.857956,0.414,-1.1],[1.94156,0.56163,1.09],[1.88825,0.6992,1.09],[1.94156,0.56163,1.052],[1.88825,0.6992,1.052],[1.94156,0.56163,1.052],[1.88825,0.6992,1.09],[-1.0236,0.753372,1.09],[-1.06017,0.807097,1.09],[-1.0236,0.753372,1.052],[-1.06017,0.807097,1.052],[-1.0236,0.753372,1.052],[-1.06017,0.807097,1.09],[1.95888,0.414,1.052],[1.94156,0.56163,1.052],[1.99191,0.4692,1.052],[1.97355,0.61683,1.052],[1.99191,0.4692,1.052],[1.94156,0.56163,1.052],[-1.01076,0.69,1.052],[-1.0236,0.753372,1.052],[-0.997156,0.7452,1.052],[-1.01051,0.808572,1.052],[-0.997156,0.7452,1.052],[-1.0236,0.753372,1.052],[-2.28344,0.4692,-1.062],[-2.28344,0.414,-1.062],[-1.87155,0.4692,-1.062],[-1.85152,0.414,-1.062],[-1.87155,0.4692,-1.062],[-2.28344,0.414,-1.062],[-0.830148,0.4692,-1.062],[-0.850175,0.414,-1.062],[0.824928,0.4692,-1.062],[0.857956,0.414,-1.062],[0.824928,0.4692,-1.062],[-0.850175,0.414,-1.062],[2.29264,0.4692,-1.062],[2.29264,0.4692,1.052],[2.25032,0.4692,-0.55094],[2.25032,0.4692,0.54094],[2.25032,0.4692,-0.55094],[2.29264,0.4692,1.052],[0.824928,0.4692,1.052],[-0.830148,0.4692,1.052],[0.824928,0.4692,0.54094],[-0.830148,0.4692,0.54094],[0.824928,0.4692,0.54094],[-0.830148,0.4692,1.052],[-2.28344,0.4692,1.052],[-2.24112,0.4692,0.54094],[-1.87155,0.4692,1.052],[-1.87155,0.4692,0.54094],[-1.87155,0.4692,1.052],[-2.24112,0.4692,0.54094],[-1.01051,0.808572,-1.062],[-0.997156,0.7452,-1.062],[-1.01051,0.808572,-0.55094],[-0.997156,0.7452,-0.55094],[-1.01051,0.808572,-0.55094],[-0.997156,0.7452,-1.062],[-1.71527,0.7452,0.54094],[-1.87155,0.4692,0.54094],[-0.830148,0.4692,0.54094],[-1.71527,0.7452,0.54094],[-0.830148,0.4692,0.54094],[-0.997156,0.7452,0.54094],[-1.65961,0.862297,0.54094],[-1.7008,0.808572,0.54094],[-1.71527,0.7452,0.54094],[-1.65961,0.862297,0.54094],[-1.71527,0.7452,0.54094],[-0.997156,0.7452,0.54094],[-0.997156,0.7452,0.54094],[-1.01051,0.808572,0.54094],[-1.04855,0.862297,0.54094],[-1.65961,0.862297,0.54094],[-0.997156,0.7452,0.54094],[-1.04855,0.862297,0.54094],[-1.04855,0.862297,0.54094],[-1.10547,0.898194,0.54094],[-1.17262,0.9108,0.54094],[-1.52526,0.9108,0.54094],[-1.59797,0.898194,0.54094],[-1.65961,0.862297,0.54094],[-1.65961,0.862297,0.54094],[-1.04855,0.862297,0.54094],[-1.17262,0.9108,0.54094],[-1.17262,0.9108,0.54094],[-1.52526,0.9108,0.54094],[-1.65961,0.862297,0.54094],[-2.42748,0.371258,-1.12347],[-2.18632,0.371258,-1.12347],[-2.41077,0.3634,-1.1054],[-2.20303,0.3634,-1.1054],[-2.41077,0.3634,-1.1054],[-2.18632,0.371258,-1.12347],[-2.42436,0.475849,-1.1196],[-2.18944,0.475849,-1.1196],[-2.42329,0.470631,-1.11838],[-2.19051,0.470631,-1.11838],[-2.42329,0.470631,-1.11838],[-2.18944,0.475849,-1.1196],[-2.42329,0.470631,1.11838],[-2.42436,0.475849,1.1196],[-2.42436,0.475849,-1.1196],[-2.42436,0.475849,-1.1196],[-2.42329,0.470631,-1.11838],[-2.42329,0.470631,1.11838],[-2.18289,0.49916,1.12699],[-2.43091,0.49916,1.12699],[-2.18944,0.475849,1.1196],[-2.42436,0.475849,1.1196],[-2.18944,0.475849,1.1196],[-2.43091,0.49916,1.12699],[2.18408,0.50784,-1.1115],[2.18408,0.4899,-1.1115],[2.18408,0.4899,1.1115],[2.18408,0.4899,1.1115],[2.18408,0.50784,1.1115],[2.18408,0.50784,-1.1115],[2.20202,0.3864,-1.1115],[2.20202,0.391655,-1.12529],[2.31518,0.3864,-1.1115],[2.31799,0.391744,-1.12529],[2.31518,0.3864,-1.1115],[2.20202,0.391655,-1.12529],[2.20202,0.3864,1.1115],[2.31518,0.3864,1.1115],[2.20202,0.391655,1.12529],[2.31799,0.391744,1.12529],[2.20202,0.391655,1.12529],[2.31518,0.3864,1.1115],[2.42006,0.66838,-1.131],[2.42006,0.681065,-1.12529],[2.43275,0.66838,-1.12529],[2.43099,0.679314,-1.12338],[2.43275,0.66838,-1.12529],[2.42006,0.681065,-1.12529],[2.20202,0.3864,1.1115],[2.20202,0.391655,1.12529],[2.18933,0.391655,1.1115],[2.19109,0.393406,1.12338],[2.18933,0.391655,1.1115],[2.20202,0.391655,1.12529],[2.31799,0.391744,1.12529],[2.33066,0.391774,1.12382],[2.32565,0.402713,1.131],[2.33897,0.401345,1.12529],[2.32565,0.402713,1.131],[2.33066,0.391774,1.12382],[2.438,0.66838,1.1115],[2.438,0.52578,1.1115],[2.438,0.66838,-1.1115],[2.438,0.52578,-1.1115],[2.438,0.66838,-1.1115],[2.438,0.52578,1.1115],[2.3276,0.9936,-0.505],[2.3276,1.15,-0.505],[2.27884,0.9936,-0.505],[2.27884,1.15,-0.505],[2.27884,0.9936,-0.505],[2.3276,1.15,-0.505],[2.3276,1.1914,0.645],[2.3276,1.1914,0.92],[2.31564,1.1914,0.645],[2.31564,1.1914,0.92],[2.31564,1.1914,0.645],[2.3276,1.1914,0.92],[1.3271,0.3358,0.166501],[1.3271,0.4646,0.166501],[1.3271,0.3358,-0.166501],[1.3271,0.4646,-0.166501],[1.3271,0.3358,-0.166501],[1.3271,0.4646,0.166501],[1.50046,0.2668,0.166501],[1.36074,0.2668,0.166501],[1.50046,0.2668,-0.166501],[1.36074,0.2668,-0.166501],[1.50046,0.2668,-0.166501],[1.36074,0.2668,0.166501],[1.3271,0.375469,-0.321545],[1.3271,0.382808,-0.287501],[1.3271,0.504269,-0.321545],[1.3271,0.511608,-0.287501],[1.3271,0.504269,-0.321545],[1.3271,0.382808,-0.287501],[1.5341,0.375469,0.253456],[1.5341,0.504269,0.253456],[1.5341,0.382808,0.287501],[1.5341,0.511608,0.287501],[1.5341,0.382808,0.287501],[1.5341,0.504269,0.253456],[1.5341,0.375469,-0.321545],[1.50046,0.306469,-0.321545],[1.5341,0.3358,-0.408501],[1.50046,0.2668,-0.408501],[1.5341,0.3358,-0.408501],[1.50046,0.306469,-0.321545],[1.50046,0.2668,-0.408501],[1.50046,0.306469,-0.321545],[1.36074,0.2668,-0.408501],[1.36074,0.306469,-0.321545],[1.36074,0.2668,-0.408501],[1.50046,0.306469,-0.321545],[-1.2489,0.4646,-0.4085],[-1.2489,0.3358,-0.4085],[-1.2489,0.4646,-0.8625],[-1.2489,0.3358,-0.8625],[-1.2489,0.4646,-0.8625],[-1.2489,0.3358,-0.4085],[-1.42226,0.306469,0.321545],[-1.4559,0.375469,0.321545],[-1.42226,0.313808,0.287501],[-1.4559,0.382808,0.287501],[-1.42226,0.313808,0.287501],[-1.4559,0.375469,0.321545],[-1.28254,0.306469,-0.253456],[-1.28254,0.313808,-0.287501],[-1.2489,0.375469,-0.253456],[-1.2489,0.382808,-0.287501],[-1.2489,0.375469,-0.253456],[-1.28254,0.313808,-0.287501],[-1.4559,0.3358,0.166501],[-1.4559,0.375469,0.253456],[-1.4559,0.4646,0.166501],[-1.4559,0.504269,0.253456],[-1.4559,0.4646,0.166501],[-1.4559,0.375469,0.253456],[-1.42226,0.2668,-0.4085],[-1.42226,0.306469,-0.321545],[-1.4559,0.3358,-0.4085],[-1.4559,0.375469,-0.321545],[-1.4559,0.3358,-0.4085],[-1.42226,0.306469,-0.321545],[-2.3184,1.2696,0.6075],[-2.3184,1.6192,0.6075],[-2.2954,1.2696,0.6075],[-2.2954,1.6192,0.6075],[-2.2954,1.2696,0.6075],[-2.3184,1.6192,0.6075],[-2.33266,0.7176,0.6075],[-2.34692,0.7176,0.6075],[-2.33266,0.5244,0.6075],[-2.34692,0.5244,0.6075],[-2.33266,0.5244,0.6075],[-2.34692,0.7176,0.6075],[1.0028,0.8924,1.055],[1.0028,0.8924,1.09],[1.0028,1.242,1.055],[1.0028,1.242,1.09],[1.0028,1.242,1.055],[1.0028,0.8924,1.09],[0.155784,2.38342,1.09],[0.142622,2.36372,1.09],[0.155784,2.38342,1.055],[0.142622,2.36372,1.055],[0.155784,2.38342,1.055],[0.142622,2.36372,1.09],[1.0028,0.8924,-1.1],[0.7912,0.4784,-1.1],[1.0028,0.8924,-1.055],[0.7912,0.4784,-1.055],[1.0028,0.8924,-1.055],[0.7912,0.4784,-1.1],[0.985015,2.38342,-1.055],[0.985016,2.38342,-1.1],[0.998178,2.36372,-1.055],[0.998178,2.36372,-1.1],[0.998178,2.36372,-1.055],[0.985016,2.38342,-1.1],[1.54254,1.87378,-1.20187],[1.54247,1.887,-1.2074],[1.53896,1.87858,-1.19766],[1.53894,1.89212,-1.20364],[1.53896,1.87858,-1.19766],[1.54247,1.887,-1.2074],[1.54254,1.65902,-1.20187],[1.53896,1.65421,-1.19766],[1.54247,1.6458,-1.2074],[1.53894,1.64068,-1.20364],[1.54247,1.6458,-1.2074],[1.53896,1.65421,-1.19766],[1.53226,1.88048,-1.196],[1.5226,1.88048,-1.196],[1.53226,1.65232,-1.196],[1.5226,1.65232,-1.196],[1.53226,1.65232,-1.196],[1.5226,1.88048,-1.196],[1.55986,1.79032,-1.2754],[1.55986,1.79032,-1.308],[1.55986,1.62619,-1.2754],[1.55986,1.5962,-1.308],[1.55986,1.62619,-1.2754],[1.55986,1.79032,-1.308],[2.32089,1.2654,1.044],[2.30279,1.28169,1.044],[2.30849,1.26544,1.07653],[2.29039,1.28184,1.07653],[2.30849,1.26544,1.07653],[2.30279,1.28169,1.044],[2.3276,1.242,-1.054],[2.3276,0.414,-1.054],[2.3152,1.242,-1.08653],[2.3152,0.414,-1.08653],[2.3152,1.242,-1.08653],[2.3276,0.414,-1.054],[1.61,1.6192,1.044],[2.30279,1.28169,1.044],[1.61,1.6192,-1.054],[2.30279,1.28169,-1.054],[1.61,1.6192,-1.054],[2.30279,1.28169,1.044],[-2.3184,1.2696,-0.6075],[-2.3184,1.6192,-0.6075],[-2.3184,1.2696,-1.054],[-2.3184,1.6192,-1.054],[-2.3184,1.2696,-1.054],[-2.3184,1.6192,-0.6075],[-2.27608,1.6192,-1.1],[-2.306,1.6192,-1.08653],[-2.27608,2.53519,-1.1],[-2.306,2.53519,-1.08653],[-2.27608,2.53519,-1.1],[-2.306,1.6192,-1.08653],[-2.27608,0.414,1.09],[-2.306,0.414,1.07653],[-2.28344,0.414,1.052],[-2.3184,0.414,1.044],[-2.28344,0.414,1.052],[-2.306,0.414,1.07653],[1.31979,2.55597,1.07451],[1.32979,2.55597,1.0482],[1.31101,2.57293,1.06494],[1.31831,2.57293,1.04572],[1.31101,2.57293,1.06494],[1.32979,2.55597,1.0482],[1.31831,2.57293,-1.05572],[1.31101,2.57293,-1.07494],[1.30129,2.58437,-1.05203],[1.29798,2.58437,-1.06074],[1.30129,2.58437,-1.05203],[1.31101,2.57293,-1.07494],[-2.27396,2.58437,1.04157],[-2.28569,2.58437,1.04203],[-2.27438,2.58437,1.05435],[-2.28238,2.58437,1.05074],[-2.27438,2.58437,1.05435],[-2.28569,2.58437,1.04203],[1.32979,2.55597,-1.0582],[1.32979,2.55597,1.0482],[1.334,2.53519,-1.054],[1.334,2.53519,1.044],[1.334,2.53519,-1.054],[1.32979,2.55597,1.0482],[1.53226,1.65232,1.396],[1.53896,1.65421,1.39434],[1.53226,1.88048,1.396],[1.53896,1.87858,1.39434],[1.53226,1.88048,1.396],[1.53896,1.65421,1.39434],[1.53226,1.633,1.217],[1.53895,1.63514,1.21827],[1.53226,1.633,1.375],[1.53895,1.63514,1.37373],[1.53226,1.633,1.375],[1.53895,1.63514,1.21827],[1.5226,1.63866,1.38985],[1.5226,1.633,1.375],[1.53226,1.63866,1.38985],[1.53226,1.633,1.375],[1.53226,1.63866,1.38985],[1.5226,1.633,1.375],[1.55986,1.58124,1.30126],[1.58746,1.58124,1.30126],[1.55986,1.5962,1.308],[1.58746,1.5962,1.308],[1.55986,1.5962,1.308],[1.58746,1.58124,1.30126],[1.58746,1.77707,1.308],[1.58746,1.77707,1.2754],[1.58358,1.78644,1.308],[1.58358,1.78644,1.2754],[1.58358,1.78644,1.308],[1.58746,1.77707,1.2754],[-1.17947,0.8556,1.09],[-1.51855,0.8556,1.09],[-1.17947,0.8556,1.052],[-1.51855,0.8556,1.052],[-1.17947,0.8556,1.052],[-1.51855,0.8556,1.09],[-1.17947,0.8556,-1.1],[-1.11491,0.842994,-1.1],[-1.17947,0.8556,-1.062],[-1.11491,0.842994,-1.062],[-1.17947,0.8556,-1.062],[-1.11491,0.842994,-1.1],[1.69294,0.907981,-1.1],[1.80345,0.817334,-1.1],[1.69294,0.907981,-1.062],[1.80345,0.817334,-1.062],[1.69294,0.907981,-1.062],[1.80345,0.817334,-1.1],[1.15996,0.911218,1.09],[1.0562,0.822,1.09],[1.15996,0.911218,1.052],[1.0562,0.822,1.052],[1.15996,0.911218,1.052],[1.0562,0.822,1.09],[1.94156,0.56163,-1.062],[1.97355,0.61683,-1.062],[1.88825,0.6992,-1.062],[1.91704,0.7544,-1.062],[1.88825,0.6992,-1.062],[1.97355,0.61683,-1.062],[1.15996,0.911218,1.052],[1.14505,0.966418,1.052],[1.28804,0.964964,1.052],[1.28081,1.02016,1.052],[1.28804,0.964964,1.052],[1.14505,0.966418,1.052],[-1.65961,0.862297,1.052],[-1.64774,0.807097,1.052],[-1.7008,0.808572,1.052],[-1.68734,0.753372,1.052],[-1.7008,0.808572,1.052],[-1.64774,0.807097,1.052],[-1.17262,0.9108,-1.062],[-1.17947,0.8556,-1.062],[-1.10547,0.898194,-1.062],[-1.11491,0.842994,-1.062],[-1.10547,0.898194,-1.062],[-1.17947,0.8556,-1.062],[1.71001,0.963181,-1.062],[1.82716,0.872534,-1.062],[1.71001,0.963181,-0.55094],[1.82716,0.872534,-0.55094],[1.71001,0.963181,-0.55094],[1.82716,0.872534,-1.062],[1.5736,1.02016,1.052],[1.42721,1.0396,1.052],[1.5736,1.02016,0.54094],[1.42721,1.0396,0.54094],[1.5736,1.02016,0.54094],[1.42721,1.0396,1.052],[-1.52526,0.9108,1.052],[-1.59797,0.898194,1.052],[-1.52526,0.9108,0.54094],[-1.59797,0.898194,0.54094],[-1.52526,0.9108,0.54094],[-1.59797,0.898194,1.052],[-1.59797,0.898194,-1.062],[-1.52526,0.9108,-1.062],[-1.59797,0.898194,-0.55094],[-1.52526,0.9108,-0.55094],[-1.59797,0.898194,-0.55094],[-1.52526,0.9108,-1.062],[2.25032,0.4692,0.54094],[1.99191,0.4692,0.54094],[2.25032,0.4692,-0.55094],[1.99191,0.4692,-0.55094],[2.25032,0.4692,-0.55094],[1.99191,0.4692,0.54094],[-2.42748,0.371258,1.12347],[-2.43091,0.38864,1.12699],[-2.43091,0.38864,-1.12699],[-2.43091,0.38864,-1.12699],[-2.42748,0.371258,-1.12347],[-2.42748,0.371258,1.12347],[-2.18289,0.38864,1.12699],[-2.43091,0.38864,1.12699],[-2.18632,0.371258,1.12347],[-2.42748,0.371258,1.12347],[-2.18632,0.371258,1.12347],[-2.43091,0.38864,1.12699],[-2.19051,0.470631,1.11838],[-2.42329,0.470631,1.11838],[-2.19088,0.465213,1.11797],[-2.42292,0.465213,1.11797],[-2.19088,0.465213,1.11797],[-2.42329,0.470631,1.11838],[-2.43091,0.38864,-1.12699],[-2.42436,0.41195,-1.1196],[-2.18289,0.38864,-1.12699],[-2.18944,0.41195,-1.1196],[-2.18289,0.38864,-1.12699],[-2.42436,0.41195,-1.1196],[2.43275,0.66838,-1.12529],[2.43275,0.52578,-1.12529],[2.42006,0.66838,-1.131],[2.42006,0.52578,-1.131],[2.42006,0.66838,-1.131],[2.43275,0.52578,-1.12529],[2.42627,0.494263,-1.1115],[2.43497,0.508966,-1.1115],[2.43497,0.508966,1.1115],[2.43497,0.508966,1.1115],[2.42627,0.494263,1.1115],[2.42627,0.494263,-1.1115],[2.43275,0.52578,1.12529],[2.43275,0.66838,1.12529],[2.42006,0.52578,1.131],[2.42006,0.66838,1.131],[2.42006,0.52578,1.131],[2.43275,0.66838,1.12529],[2.18408,0.40434,-1.1115],[2.18933,0.40434,-1.12529],[2.18933,0.391655,-1.1115],[2.19109,0.393406,-1.12338],[2.18933,0.391655,-1.1115],[2.18933,0.40434,-1.12529],[2.33897,0.401345,-1.12529],[2.33066,0.391774,-1.12382],[2.32565,0.402713,-1.131],[2.31799,0.391744,-1.12529],[2.32565,0.402713,-1.131],[2.33066,0.391774,-1.12382],[2.438,0.52578,1.1115],[2.43275,0.52578,1.12529],[2.43497,0.508966,1.1115],[2.42972,0.508966,1.12529],[2.43497,0.508966,1.1115],[2.43275,0.52578,1.12529],[2.20202,0.50784,1.131],[2.18933,0.50784,1.12529],[2.20202,0.4899,1.131],[2.18933,0.4899,1.12529],[2.20202,0.4899,1.131],[2.18933,0.50784,1.12529],[2.20202,0.40434,-1.131],[2.20202,0.4899,-1.131],[2.32565,0.402713,-1.131],[2.40833,0.494263,-1.131],[2.32565,0.402713,-1.131],[2.20202,0.4899,-1.131],[2.3276,0.7682,-0.645],[2.3276,0.7682,-0.92],[2.31564,0.7682,-0.645],[2.31564,0.7682,-0.92],[2.31564,0.7682,-0.645],[2.3276,0.7682,-0.92],[1.3271,0.3358,-0.8625],[1.5341,0.3358,-0.8625],[1.36074,0.2668,-0.8625],[1.50046,0.2668,-0.8625],[1.36074,0.2668,-0.8625],[1.5341,0.3358,-0.8625],[1.5341,0.3358,0.8625],[1.50046,0.2668,0.8625],[1.5341,0.3358,0.408501],[1.50046,0.2668,0.408501],[1.5341,0.3358,0.408501],[1.50046,0.2668,0.8625],[1.36074,0.313808,0.287501],[1.50046,0.313808,0.287501],[1.36074,0.306469,0.321545],[1.50046,0.306469,0.321545],[1.36074,0.306469,0.321545],[1.50046,0.313808,0.287501],[1.5341,0.511608,-0.287501],[1.3271,0.511608,-0.287501],[1.5341,0.504269,-0.253456],[1.3271,0.504269,-0.253456],[1.5341,0.504269,-0.253456],[1.3271,0.511608,-0.287501],[1.50046,0.306469,0.321545],[1.50046,0.2668,0.408501],[1.36074,0.306469,0.321545],[1.36074,0.2668,0.408501],[1.36074,0.306469,0.321545],[1.50046,0.2668,0.408501],[1.36074,0.2668,0.408501],[1.3271,0.3358,0.408501],[1.36074,0.306469,0.321545],[1.3271,0.375469,0.321545],[1.36074,0.306469,0.321545],[1.3271,0.3358,0.408501],[-1.4559,0.4646,0.166501],[-1.2489,0.4646,0.166501],[-1.4559,0.4646,-0.1665],[-1.2489,0.4646,-0.1665],[-1.4559,0.4646,-0.1665],[-1.2489,0.4646,0.166501],[-1.42226,0.2668,0.166501],[-1.4559,0.3358,0.166501],[-1.42226,0.2668,-0.1665],[-1.4559,0.3358,-0.1665],[-1.42226,0.2668,-0.1665],[-1.4559,0.3358,0.166501],[-1.4559,0.382808,-0.287501],[-1.4559,0.375469,-0.321545],[-1.42226,0.313808,-0.287501],[-1.42226,0.306469,-0.321545],[-1.42226,0.313808,-0.287501],[-1.4559,0.375469,-0.321545],[-1.4559,0.511608,0.287501],[-1.4559,0.382808,0.287501],[-1.4559,0.504269,0.321545],[-1.4559,0.375469,0.321545],[-1.4559,0.504269,0.321545],[-1.4559,0.382808,0.287501],[-1.2489,0.3358,0.408501],[-1.2489,0.375469,0.321545],[-1.2489,0.4646,0.408501],[-1.2489,0.504269,0.321545],[-1.2489,0.4646,0.408501],[-1.2489,0.375469,0.321545],[-1.4559,0.3358,-0.1665],[-1.4559,0.4646,-0.1665],[-1.4559,0.375469,-0.253456],[-1.4559,0.504269,-0.253456],[-1.4559,0.375469,-0.253456],[-1.4559,0.4646,-0.1665],[-2.2954,1.2696,0.6075],[-2.2954,1.2696,-0.6075],[-2.2954,0.7176,0.6075],[-2.2954,0.7176,-0.6075],[-2.2954,0.7176,0.6075],[-2.2954,1.2696,-0.6075],[-2.3184,0.5244,0.3146],[-2.33266,0.5244,0.3146],[-2.3184,0.5244,-0.3146],[-2.33266,0.5244,-0.3146],[-2.3184,0.5244,-0.3146],[-2.33266,0.5244,0.3146],[0.985016,2.38342,1.09],[0.965317,2.39658,1.09],[0.985015,2.38342,1.055],[0.965316,2.39658,1.055],[0.985015,2.38342,1.055],[0.965317,2.39658,1.09],[0.7912,0.4784,-1.055],[0.138,0.4784,-1.055],[0.138,1.2696,-1.055],[0.19872,2.4012,-1.055],[0.27968,2.3092,-1.055],[0.1932,2.22272,-1.055],[1.0028,0.8924,-1.055],[0.7912,0.4784,-1.055],[0.138,1.2696,-1.055],[1.0028,1.242,-1.055],[1.0028,0.8924,-1.055],[0.138,1.2696,-1.055],[1.0028,1.242,-1.055],[0.138,1.2696,-1.055],[0.138,1.6192,-1.055],[0.175483,2.39658,-1.055],[0.19872,2.4012,-1.055],[0.1932,2.22272,-1.055],[0.81972,2.3092,-1.055],[0.27968,2.3092,-1.055],[0.19872,2.4012,-1.055],[0.81972,2.3092,-1.055],[0.19872,2.4012,-1.055],[0.94208,2.4012,-1.055],[0.9062,2.22272,-1.055],[0.81972,2.3092,-1.055],[0.94208,2.4012,-1.055],[0.9062,2.22272,-1.055],[0.94208,2.4012,-1.055],[0.965316,2.39658,-1.055],[0.9062,2.22272,-1.055],[0.965316,2.39658,-1.055],[0.985015,2.38342,-1.055],[0.155784,2.38342,-1.055],[0.175483,2.39658,-1.055],[0.1932,2.22272,-1.055],[0.9062,2.22272,-1.055],[0.985015,2.38342,-1.055],[0.998178,2.36372,-1.055],[0.9062,1.68268,-1.055],[0.9062,2.22272,-1.055],[0.998178,2.36372,-1.055],[0.142622,2.36372,-1.055],[0.155784,2.38342,-1.055],[0.1932,2.22272,-1.055],[0.138,2.34048,-1.055],[0.142622,2.36372,-1.055],[0.1932,2.22272,-1.055],[0.138,1.6192,-1.055],[0.138,2.34048,-1.055],[0.1932,2.22272,-1.055],[0.138,1.6192,-1.055],[0.1932,2.22272,-1.055],[0.1932,1.68268,-1.055],[0.138,1.6192,-1.055],[0.1932,1.68268,-1.055],[0.27968,1.5962,-1.055],[1.0028,1.242,-1.055],[0.138,1.6192,-1.055],[0.27968,1.5962,-1.055],[1.0028,1.242,-1.055],[0.27968,1.5962,-1.055],[0.81972,1.5962,-1.055],[1.0028,1.242,-1.055],[0.81972,1.5962,-1.055],[0.9062,1.68268,-1.055],[1.0028,1.242,-1.055],[0.9062,1.68268,-1.055],[0.998178,2.36372,-1.055],[0.998178,2.36372,-1.055],[1.0028,2.34048,-1.055],[1.0028,1.6192,-1.055],[0.998178,2.36372,-1.055],[1.0028,1.6192,-1.055],[1.0028,1.2972,-1.055],[0.998178,2.36372,-1.055],[1.0028,1.2972,-1.055],[1.0028,1.2696,-1.055],[0.998178,2.36372,-1.055],[1.0028,1.2696,-1.055],[1.0028,1.242,-1.055],[0.155784,2.38342,-1.055],[0.155784,2.38342,-1.1],[0.175483,2.39658,-1.055],[0.175483,2.39658,-1.1],[0.175483,2.39658,-1.055],[0.155784,2.38342,-1.1],[1.53226,1.63866,-1.38985],[1.53894,1.64068,-1.38836],[1.53226,1.633,-1.375],[1.53895,1.63514,-1.37373],[1.53226,1.633,-1.375],[1.53894,1.64068,-1.38836],[1.53226,1.8998,-1.375],[1.53226,1.8998,-1.217],[1.53895,1.89766,-1.37373],[1.53895,1.89766,-1.21827],[1.53895,1.89766,-1.37373],[1.53226,1.8998,-1.217],[1.53226,1.63866,-1.20215],[1.5226,1.63866,-1.20215],[1.53226,1.633,-1.217],[1.5226,1.633,-1.217],[1.53226,1.633,-1.217],[1.5226,1.63866,-1.20215],[1.55986,1.61123,-1.26866],[1.55986,1.60503,-1.2524],[1.58746,1.61123,-1.26866],[1.58746,1.60503,-1.2524],[1.58746,1.61123,-1.26866],[1.55986,1.60503,-1.2524],[1.58358,1.78644,-1.308],[1.58746,1.77707,-1.308],[1.57421,1.79032,-1.308],[1.55986,1.79032,-1.308],[1.57421,1.79032,-1.308],[1.58746,1.77707,-1.308],[1.58746,1.77707,-1.308],[1.58746,1.5962,-1.308],[1.55986,1.79032,-1.308],[1.55986,1.5962,-1.308],[1.55986,1.79032,-1.308],[1.58746,1.5962,-1.308],[1.0028,1.2696,1.09],[1.0028,1.242,1.09],[2.27857,1.26554,1.09],[2.28528,1.242,1.09],[2.27857,1.26554,1.09],[1.0028,1.242,1.09],[1.0028,1.2972,-1.1],[2.26047,1.28219,-1.1],[1.0028,1.2696,-1.1],[2.27857,1.26554,-1.1],[1.0028,1.2696,-1.1],[2.26047,1.28219,-1.1],[-2.27608,1.6192,-1.1],[-2.27608,1.2696,-1.1],[-2.306,1.6192,-1.08653],[-2.306,1.2696,-1.08653],[-2.306,1.6192,-1.08653],[-2.27608,1.2696,-1.1],[-1.85152,0.414,-1.1],[-1.85152,0.414,-1.062],[-2.27608,0.414,-1.1],[-2.28344,0.414,-1.062],[-2.27608,0.414,-1.1],[-1.85152,0.414,-1.062],[-2.3184,0.5244,-0.6075],[-2.33266,0.5244,-0.6075],[-2.3184,0.7176,-0.6075],[-2.33266,0.7176,-0.6075],[-2.3184,0.7176,-0.6075],[-2.33266,0.5244,-0.6075],[1.28998,2.58437,-1.06435],[1.29337,2.57293,-1.0829],[-2.27777,2.57293,-1.0829],[-2.27777,2.57293,-1.0829],[-2.27438,2.58437,-1.06435],[1.28998,2.58437,-1.06435],[1.29337,2.57293,-1.0829],[1.28998,2.58437,-1.06435],[1.31101,2.57293,-1.07494],[1.29798,2.58437,-1.06074],[1.31101,2.57293,-1.07494],[1.28998,2.58437,-1.06435],[1.29565,2.55597,1.08541],[1.31979,2.55597,1.07451],[1.29337,2.57293,1.0729],[1.31101,2.57293,1.06494],[1.29337,2.57293,1.0729],[1.31979,2.55597,1.07451],[1.29168,2.53519,1.09],[1.29565,2.55597,1.08541],[-2.28005,2.55597,1.08541],[-2.28005,2.55597,1.08541],[-2.27608,2.53519,1.09],[1.29168,2.53519,1.09],[1.54254,1.87378,1.20187],[1.53896,1.87858,1.19766],[1.54247,1.887,1.2074],[1.53894,1.89212,1.20364],[1.54247,1.887,1.2074],[1.53896,1.87858,1.19766],[1.54254,1.65902,1.20187],[1.54247,1.6458,1.2074],[1.53896,1.65421,1.19766],[1.53894,1.64068,1.20364],[1.53896,1.65421,1.19766],[1.54247,1.6458,1.2074],[1.53226,1.88048,1.196],[1.53226,1.65232,1.196],[1.5226,1.88048,1.196],[1.5226,1.65232,1.196],[1.5226,1.88048,1.196],[1.53226,1.65232,1.196],[1.55986,1.79032,1.2754],[1.55986,1.62619,1.2754],[1.55986,1.79032,1.308],[1.55986,1.5962,1.308],[1.55986,1.79032,1.308],[1.55986,1.62619,1.2754],[-1.68734,0.753372,-1.1],[-1.64774,0.807097,-1.1],[-1.68734,0.753372,-1.062],[-1.64774,0.807097,-1.062],[-1.68734,0.753372,-1.062],[-1.64774,0.807097,-1.1],[1.15996,0.911218,-1.1],[1.28804,0.964964,-1.1],[1.15996,0.911218,-1.062],[1.28804,0.964964,-1.062],[1.15996,0.911218,-1.062],[1.28804,0.964964,-1.1],[1.69294,0.907981,1.09],[1.56425,0.964964,1.09],[1.69294,0.907981,1.052],[1.56425,0.964964,1.052],[1.69294,0.907981,1.052],[1.56425,0.964964,1.09],[1.42721,1.0396,-1.062],[1.42614,0.9844,-1.062],[1.5736,1.02016,-1.062],[1.56425,0.964964,-1.062],[1.5736,1.02016,-1.062],[1.42614,0.9844,-1.062],[1.80345,0.817334,1.052],[1.69294,0.907981,1.052],[1.82716,0.872534,1.052],[1.71001,0.963181,1.052],[1.82716,0.872534,1.052],[1.69294,0.907981,1.052],[-1.11491,0.842994,1.052],[-1.17947,0.8556,1.052],[-1.10547,0.898194,1.052],[-1.17262,0.9108,1.052],[-1.10547,0.898194,1.052],[-1.17947,0.8556,1.052],[-1.68734,0.753372,-1.062],[-1.64774,0.807097,-1.062],[-1.7008,0.808572,-1.062],[-1.65961,0.862297,-1.062],[-1.7008,0.808572,-1.062],[-1.64774,0.807097,-1.062],[1.15996,0.911218,-1.062],[1.28804,0.964964,-1.062],[1.14505,0.966418,-1.062],[1.28081,1.02016,-1.062],[1.14505,0.966418,-1.062],[1.28804,0.964964,-1.062],[1.97355,0.61683,1.052],[1.91704,0.7544,1.052],[1.97355,0.61683,0.54094],[1.91704,0.7544,0.54094],[1.97355,0.61683,0.54094],[1.91704,0.7544,1.052],[-1.01051,0.808572,1.052],[-1.04855,0.862297,1.052],[-1.01051,0.808572,0.54094],[-1.04855,0.862297,0.54094],[-1.01051,0.808572,0.54094],[-1.04855,0.862297,1.052],[-1.87155,0.4692,-1.062],[-1.71527,0.7452,-1.062],[-1.87155,0.4692,-0.55094],[-1.71527,0.7452,-0.55094],[-1.87155,0.4692,-0.55094],[-1.71527,0.7452,-1.062],[0.824928,0.4692,-1.062],[1.03507,0.8772,-1.062],[0.824928,0.4692,-0.55094],[1.03507,0.8772,-0.55094],[0.824928,0.4692,-0.55094],[1.03507,0.8772,-1.062],[-0.997156,0.7452,-0.55094],[-0.830148,0.4692,-0.55094],[-1.87155,0.4692,-0.55094],[-0.997156,0.7452,-0.55094],[-1.87155,0.4692,-0.55094],[-1.71527,0.7452,-0.55094],[-1.71527,0.7452,-0.55094],[-1.7008,0.808572,-0.55094],[-1.65961,0.862297,-0.55094],[-0.997156,0.7452,-0.55094],[-1.71527,0.7452,-0.55094],[-1.65961,0.862297,-0.55094],[-1.04855,0.862297,-0.55094],[-1.01051,0.808572,-0.55094],[-0.997156,0.7452,-0.55094],[-1.04855,0.862297,-0.55094],[-0.997156,0.7452,-0.55094],[-1.65961,0.862297,-0.55094],[-1.17262,0.9108,-0.55094],[-1.10547,0.898194,-0.55094],[-1.04855,0.862297,-0.55094],[-1.65961,0.862297,-0.55094],[-1.59797,0.898194,-0.55094],[-1.52526,0.9108,-0.55094],[-1.04855,0.862297,-0.55094],[-1.65961,0.862297,-0.55094],[-1.52526,0.9108,-0.55094],[-1.52526,0.9108,-0.55094],[-1.17262,0.9108,-0.55094],[-1.04855,0.862297,-0.55094],[-2.43091,0.49916,1.12699],[-2.18289,0.49916,1.12699],[-2.42748,0.516542,1.12347],[-2.18632,0.516542,1.12347],[-2.42748,0.516542,1.12347],[-2.18289,0.49916,1.12699],[-2.18944,0.41195,-1.1196],[-2.42436,0.41195,-1.1196],[-2.19051,0.417169,-1.11838],[-2.42329,0.417169,-1.11838],[-2.19051,0.417169,-1.11838],[-2.42436,0.41195,-1.1196],[-2.42292,0.422587,-1.11797],[-2.42292,0.422587,1.11797],[-2.42292,0.465213,-1.11797],[-2.42292,0.465213,1.11797],[-2.42292,0.465213,-1.11797],[-2.42292,0.422587,1.11797],[2.20202,0.681065,-1.12529],[2.42006,0.681065,-1.12529],[2.20202,0.66838,-1.131],[2.42006,0.66838,-1.131],[2.20202,0.66838,-1.131],[2.42006,0.681065,-1.12529],[2.18408,0.66838,1.1115],[2.18933,0.681065,1.1115],[2.18933,0.681065,-1.1115],[2.18933,0.681065,-1.1115],[2.18408,0.66838,-1.1115],[2.18408,0.66838,1.1115],[2.32565,0.402713,1.131],[2.33897,0.401345,1.12529],[2.40833,0.494263,1.131],[2.42102,0.494263,1.12529],[2.40833,0.494263,1.131],[2.33897,0.401345,1.12529],[2.42006,0.681065,1.12529],[2.20202,0.681065,1.12529],[2.42006,0.66838,1.131],[2.20202,0.66838,1.131],[2.42006,0.66838,1.131],[2.20202,0.681065,1.12529],[2.40833,0.494263,-1.131],[2.41703,0.508966,-1.131],[2.42102,0.494263,-1.12529],[2.42972,0.508966,-1.12529],[2.42102,0.494263,-1.12529],[2.41703,0.508966,-1.131],[2.20202,0.68632,1.1115],[2.18933,0.681065,1.1115],[2.20202,0.681065,1.12529],[2.19109,0.679314,1.12338],[2.20202,0.681065,1.12529],[2.18933,0.681065,1.1115],[2.18933,0.50784,-1.12529],[2.18933,0.52578,-1.12529],[2.20202,0.50784,-1.131],[2.20202,0.52578,-1.131],[2.20202,0.50784,-1.131],[2.18933,0.52578,-1.12529],[2.18408,0.40434,-1.1115],[2.18408,0.40434,1.1115],[2.18408,0.4899,-1.1115],[2.18408,0.4899,1.1115],[2.18408,0.4899,-1.1115],[2.18408,0.40434,1.1115],[2.3276,0.9936,0.505],[2.3276,0.9936,-0.505],[2.27884,0.9936,0.505],[2.27884,0.9936,-0.505],[2.27884,0.9936,0.505],[2.3276,0.9936,-0.505],[2.31564,1.1914,-0.645],[2.31564,0.7682,-0.645],[2.31564,1.1914,-0.92],[2.31564,0.7682,-0.92],[2.31564,1.1914,-0.92],[2.31564,0.7682,-0.645],[1.5341,0.4646,0.166501],[1.5341,0.3358,0.166501],[1.5341,0.4646,-0.166501],[1.5341,0.3358,-0.166501],[1.5341,0.4646,-0.166501],[1.5341,0.3358,0.166501],[1.3271,0.382808,0.287501],[1.3271,0.375469,0.253456],[1.36074,0.313808,0.287501],[1.36074,0.306469,0.253456],[1.36074,0.313808,0.287501],[1.3271,0.375469,0.253456],[1.5341,0.382808,-0.287501],[1.50046,0.313808,-0.287501],[1.5341,0.375469,-0.321545],[1.50046,0.306469,-0.321545],[1.5341,0.375469,-0.321545],[1.50046,0.313808,-0.287501],[1.3271,0.3358,0.408501],[1.3271,0.4646,0.408501],[1.3271,0.375469,0.321545],[1.3271,0.504269,0.321545],[1.3271,0.375469,0.321545],[1.3271,0.4646,0.408501],[1.36074,0.2668,-0.166501],[1.3271,0.3358,-0.166501],[1.36074,0.306469,-0.253456],[1.3271,0.375469,-0.253456],[1.36074,0.306469,-0.253456],[1.3271,0.3358,-0.166501],[-1.2489,0.3358,0.8625],[-1.2489,0.4646,0.8625],[-1.4559,0.3358,0.8625],[-1.4559,0.4646,0.8625],[-1.4559,0.3358,0.8625],[-1.2489,0.4646,0.8625],[-1.2489,0.3358,-0.4085],[-1.28254,0.2668,-0.4085],[-1.2489,0.3358,-0.8625],[-1.28254,0.2668,-0.8625],[-1.2489,0.3358,-0.8625],[-1.28254,0.2668,-0.4085],[-1.2489,0.382808,0.287501],[-1.28254,0.313808,0.287501],[-1.2489,0.375469,0.253456],[-1.28254,0.306469,0.253456],[-1.2489,0.375469,0.253456],[-1.28254,0.313808,0.287501],[-1.2489,0.511608,-0.287501],[-1.2489,0.504269,-0.253456],[-1.2489,0.382808,-0.287501],[-1.2489,0.375469,-0.253456],[-1.2489,0.382808,-0.287501],[-1.2489,0.504269,-0.253456],[-1.2489,0.3358,-0.1665],[-1.2489,0.375469,-0.253456],[-1.2489,0.4646,-0.1665],[-1.2489,0.504269,-0.253456],[-1.2489,0.4646,-0.1665],[-1.2489,0.375469,-0.253456],[-1.2489,0.375469,0.321545],[-1.2489,0.3358,0.408501],[-1.28254,0.306469,0.321545],[-1.28254,0.2668,0.408501],[-1.28254,0.306469,0.321545],[-1.2489,0.3358,0.408501],[-2.3184,2.4104,-0.6075],[-2.3184,1.6192,-0.6075],[-2.2954,2.4104,-0.6075],[-2.2954,1.6192,-0.6075],[-2.2954,2.4104,-0.6075],[-2.3184,1.6192,-0.6075],[-2.3184,0.5244,-0.3146],[-2.33266,0.5244,-0.3146],[-2.3184,0.5244,-0.6075],[-2.33266,0.5244,-0.6075],[-2.3184,0.5244,-0.6075],[-2.33266,0.5244,-0.3146],[1.0028,1.2972,1.055],[1.0028,1.2972,1.09],[1.0028,1.6192,1.055],[1.0028,1.6192,1.09],[1.0028,1.6192,1.055],[1.0028,1.2972,1.09],[0.138,1.6192,1.055],[0.138,1.6192,1.09],[0.138,1.2696,1.055],[0.138,1.2696,1.09],[0.138,1.2696,1.055],[0.138,1.6192,1.09],[0.138,1.2696,-1.055],[0.138,1.2696,-1.1],[0.138,1.6192,-1.055],[0.138,1.6192,-1.1],[0.138,1.6192,-1.055],[0.138,1.2696,-1.1],[1.53226,1.89414,-1.20215],[1.53894,1.89212,-1.20364],[1.53226,1.8998,-1.217],[1.53895,1.89766,-1.21827],[1.53226,1.8998,-1.217],[1.53894,1.89212,-1.20364],[1.53226,1.89414,-1.38985],[1.53894,1.89212,-1.38836],[1.53226,1.88048,-1.396],[1.53896,1.87858,-1.39434],[1.53226,1.88048,-1.396],[1.53894,1.89212,-1.38836],[1.53226,1.633,-1.217],[1.5226,1.633,-1.217],[1.53226,1.633,-1.375],[1.5226,1.633,-1.375],[1.53226,1.633,-1.375],[1.5226,1.633,-1.217],[1.55986,1.60503,-1.09],[1.55986,1.60503,-1.2524],[1.55986,1.57504,-1.09],[1.55986,1.57504,-1.285],[1.55986,1.57504,-1.09],[1.55986,1.60503,-1.2524],[1.56768,1.6192,1.09],[2.26047,1.28219,1.09],[1.5976,1.6192,1.07653],[2.29039,1.28184,1.07653],[1.5976,1.6192,1.07653],[2.26047,1.28219,1.09],[2.27857,1.26554,-1.1],[2.30849,1.26544,-1.08653],[2.28528,1.242,-1.1],[2.3152,1.242,-1.08653],[2.28528,1.242,-1.1],[2.30849,1.26544,-1.08653],[-2.27608,0.414,1.09],[-2.27608,1.2696,1.09],[-2.306,0.414,1.07653],[-2.306,1.2696,1.07653],[-2.306,0.414,1.07653],[-2.27608,1.2696,1.09],[1.95888,0.414,1.09],[1.95888,0.414,1.052],[2.28528,0.414,1.09],[2.29264,0.414,1.052],[2.28528,0.414,1.09],[1.95888,0.414,1.052],[-2.3184,0.7176,-0.6075],[-2.3184,0.7176,0.6075],[-2.2954,0.7176,-0.6075],[-2.2954,0.7176,0.6075],[-2.2954,0.7176,-0.6075],[-2.3184,0.7176,0.6075],[2.3276,0.414,1.044],[2.3152,0.414,1.07653],[2.29264,0.414,1.052],[2.28528,0.414,1.09],[2.29264,0.414,1.052],[2.3152,0.414,1.07653],[1.29337,2.57293,1.0729],[1.31101,2.57293,1.06494],[1.28998,2.58437,1.05435],[1.29798,2.58437,1.05074],[1.28998,2.58437,1.05435],[1.31101,2.57293,1.06494],[-2.30419,2.55597,1.07451],[-2.28005,2.55597,1.08541],[-2.29541,2.57293,1.06494],[-2.27777,2.57293,1.0729],[-2.29541,2.57293,1.06494],[-2.28005,2.55597,1.08541],[-2.27396,2.58437,1.04157],[-2.27438,2.58437,1.05435],[1.28998,2.58437,1.05435],[1.28998,2.58437,1.05435],[1.28956,2.58437,1.04157],[-2.27396,2.58437,1.04157],[1.0442,2.3483,-1.1],[1.0028,2.34048,-1.1],[0.998178,2.36372,-1.1],[0.138,1.6192,-1.1],[-2.27608,1.6192,-1.1],[-2.27608,2.53519,-1.1],[1.0879,2.392,-1.1],[1.0442,2.3483,-1.1],[0.998178,2.36372,-1.1],[0.138,2.34048,-1.1],[0.138,1.6192,-1.1],[-2.27608,2.53519,-1.1],[1.0879,2.392,-1.1],[0.998178,2.36372,-1.1],[0.985016,2.38342,-1.1],[1.0879,2.392,-1.1],[0.985016,2.38342,-1.1],[0.965317,2.39658,-1.1],[1.0028,2.34048,-1.1],[1.0442,2.3483,-1.1],[1.0442,1.6767,-1.1],[1.0028,1.6192,-1.1],[1.0028,2.34048,-1.1],[1.0442,1.6767,-1.1],[1.0028,1.6192,-1.1],[1.0442,1.6767,-1.1],[1.0879,1.633,-1.1],[1.56768,1.6192,-1.1],[1.0028,1.6192,-1.1],[1.0879,1.633,-1.1],[1.56768,1.6192,-1.1],[1.0879,1.633,-1.1],[1.46502,1.633,-1.1],[1.56768,1.6192,-1.1],[1.46502,1.633,-1.1],[1.53021,1.6919,-1.1],[1.0879,2.392,-1.1],[0.965317,2.39658,-1.1],[0.94208,2.4012,-1.1],[1.2489,2.392,-1.1],[1.0879,2.392,-1.1],[0.94208,2.4012,-1.1],[1.29168,2.53519,-1.1],[1.56768,1.6192,-1.1],[1.53021,1.6919,-1.1],[1.29168,2.53519,-1.1],[1.53021,1.6919,-1.1],[1.30664,2.35062,-1.1],[1.29168,2.53519,-1.1],[1.30664,2.35062,-1.1],[1.2489,2.392,-1.1],[1.29168,2.53519,-1.1],[1.2489,2.392,-1.1],[0.94208,2.4012,-1.1],[1.29168,2.53519,-1.1],[0.94208,2.4012,-1.1],[0.19872,2.4012,-1.1],[-2.27608,2.53519,-1.1],[1.29168,2.53519,-1.1],[0.19872,2.4012,-1.1],[0.142622,2.36372,-1.1],[0.138,2.34048,-1.1],[-2.27608,2.53519,-1.1],[0.155784,2.38342,-1.1],[0.142622,2.36372,-1.1],[-2.27608,2.53519,-1.1],[0.175483,2.39658,-1.1],[0.155784,2.38342,-1.1],[-2.27608,2.53519,-1.1],[-2.27608,2.53519,-1.1],[0.19872,2.4012,-1.1],[0.175483,2.39658,-1.1],[1.54254,1.65902,1.39013],[1.53896,1.65421,1.39434],[1.54247,1.6458,1.3846],[1.53894,1.64068,1.38836],[1.54247,1.6458,1.3846],[1.53896,1.65421,1.39434],[1.53894,1.64068,1.20364],[1.54247,1.6458,1.2074],[1.53895,1.63514,1.21827],[1.54249,1.64056,1.22148],[1.53895,1.63514,1.21827],[1.54247,1.6458,1.2074],[1.5226,1.8998,1.375],[1.5226,1.89414,1.38985],[1.53226,1.8998,1.375],[1.53226,1.89414,1.38985],[1.53226,1.8998,1.375],[1.5226,1.89414,1.38985],[1.58746,1.62619,1.2754],[1.58746,1.5962,1.308],[1.58746,1.61123,1.26866],[1.58746,1.58124,1.30126],[1.58746,1.61123,1.26866],[1.58746,1.5962,1.308],[1.58746,1.77707,1.2754],[1.58746,1.77707,1.308],[1.58746,1.62619,1.2754],[1.58746,1.5962,1.308],[1.58746,1.62619,1.2754],[1.58746,1.77707,1.308],[-1.70125,0.69,-1.1],[-1.68734,0.753372,-1.1],[-1.70125,0.69,-1.062],[-1.68734,0.753372,-1.062],[-1.70125,0.69,-1.062],[-1.68734,0.753372,-1.1],[1.0562,0.822,-1.1],[1.15996,0.911218,-1.1],[1.0562,0.822,-1.062],[1.15996,0.911218,-1.062],[1.0562,0.822,-1.062],[1.15996,0.911218,-1.1],[1.80345,0.817334,1.09],[1.69294,0.907981,1.09],[1.80345,0.817334,1.052],[1.69294,0.907981,1.052],[1.80345,0.817334,1.052],[1.69294,0.907981,1.09],[-1.11491,0.842994,1.09],[-1.17947,0.8556,1.09],[-1.11491,0.842994,1.052],[-1.17947,0.8556,1.052],[-1.11491,0.842994,1.052],[-1.17947,0.8556,1.09],[1.88825,0.6992,1.052],[1.80345,0.817334,1.052],[1.91704,0.7544,1.052],[1.82716,0.872534,1.052],[1.91704,0.7544,1.052],[1.80345,0.817334,1.052],[-1.06017,0.807097,1.052],[-1.11491,0.842994,1.052],[-1.04855,0.862297,1.052],[-1.10547,0.898194,1.052],[-1.04855,0.862297,1.052],[-1.11491,0.842994,1.052],[-1.70125,0.69,-1.062],[-1.68734,0.753372,-1.062],[-1.71527,0.7452,-1.062],[-1.7008,0.808572,-1.062],[-1.71527,0.7452,-1.062],[-1.68734,0.753372,-1.062],[1.0562,0.822,-1.062],[1.15996,0.911218,-1.062],[1.03507,0.8772,-1.062],[1.14505,0.966418,-1.062],[1.03507,0.8772,-1.062],[1.15996,0.911218,-1.062],[1.99191,0.4692,1.052],[1.97355,0.61683,1.052],[1.99191,0.4692,0.54094],[1.97355,0.61683,0.54094],[1.99191,0.4692,0.54094],[1.97355,0.61683,1.052],[-0.997156,0.7452,1.052],[-1.01051,0.808572,1.052],[-0.997156,0.7452,0.54094],[-1.01051,0.808572,0.54094],[-0.997156,0.7452,0.54094],[-1.01051,0.808572,1.052],[-2.28344,0.4692,-1.062],[-1.87155,0.4692,-1.062],[-2.24112,0.4692,-0.55094],[-1.87155,0.4692,-0.55094],[-2.24112,0.4692,-0.55094],[-1.87155,0.4692,-1.062],[-0.830148,0.4692,-1.062],[0.824928,0.4692,-1.062],[-0.830148,0.4692,-0.55094],[0.824928,0.4692,-0.55094],[-0.830148,0.4692,-0.55094],[0.824928,0.4692,-1.062],[0.824928,0.4692,-0.55094],[1.03507,0.8772,-0.55094],[1.99191,0.4692,-0.55094],[1.97355,0.61683,-0.55094],[1.99191,0.4692,-0.55094],[1.03507,0.8772,-0.55094],[1.03507,0.8772,-0.55094],[1.14505,0.966418,-0.55094],[1.97355,0.61683,-0.55094],[1.91704,0.7544,-0.55094],[1.97355,0.61683,-0.55094],[1.14505,0.966418,-0.55094],[1.14505,0.966418,-0.55094],[1.28081,1.02016,-0.55094],[1.91704,0.7544,-0.55094],[1.82716,0.872534,-0.55094],[1.91704,0.7544,-0.55094],[1.28081,1.02016,-0.55094],[1.28081,1.02016,-0.55094],[1.42721,1.0396,-0.55094],[1.82716,0.872534,-0.55094],[1.71001,0.963181,-0.55094],[1.82716,0.872534,-0.55094],[1.42721,1.0396,-0.55094],[1.42721,1.0396,-0.55094],[1.5736,1.02016,-0.55094],[1.71001,0.963181,-0.55094],[-2.42748,0.516542,1.12347],[-2.18632,0.516542,1.12347],[-2.41077,0.5244,1.1054],[-2.20303,0.5244,1.1054],[-2.41077,0.5244,1.1054],[-2.18632,0.516542,1.12347],[-2.19051,0.417169,-1.11838],[-2.42329,0.417169,-1.11838],[-2.19088,0.422587,-1.11797],[-2.42292,0.422587,-1.11797],[-2.19088,0.422587,-1.11797],[-2.42329,0.417169,-1.11838],[-2.19088,0.422587,-1.11797],[-2.42292,0.422587,-1.11797],[-2.19088,0.465213,-1.11797],[-2.42292,0.465213,-1.11797],[-2.19088,0.465213,-1.11797],[-2.42292,0.422587,-1.11797],[2.20202,0.68632,-1.1115],[2.42006,0.68632,-1.1115],[2.20202,0.681065,-1.12529],[2.42006,0.681065,-1.12529],[2.20202,0.681065,-1.12529],[2.42006,0.68632,-1.1115],[2.18933,0.681065,1.1115],[2.20202,0.68632,1.1115],[2.20202,0.68632,-1.1115],[2.20202,0.68632,-1.1115],[2.18933,0.681065,-1.1115],[2.18933,0.681065,1.1115],[2.33897,0.401345,1.12529],[2.34485,0.399977,1.1115],[2.42102,0.494263,1.12529],[2.42627,0.494263,1.1115],[2.42102,0.494263,1.12529],[2.34485,0.399977,1.1115],[2.42006,0.68632,1.1115],[2.20202,0.68632,1.1115],[2.42006,0.681065,1.12529],[2.20202,0.681065,1.12529],[2.42006,0.681065,1.12529],[2.20202,0.68632,1.1115],[2.438,0.66838,-1.1115],[2.43275,0.66838,-1.12529],[2.43275,0.681065,-1.1115],[2.43099,0.679314,-1.12338],[2.43275,0.681065,-1.1115],[2.43275,0.66838,-1.12529],[2.20202,0.66838,1.131],[2.20202,0.681065,1.12529],[2.18933,0.66838,1.12529],[2.19109,0.679314,1.12338],[2.18933,0.66838,1.12529],[2.20202,0.681065,1.12529],[2.18933,0.50784,-1.12529],[2.18408,0.50784,-1.1115],[2.18933,0.52578,-1.12529],[2.18408,0.52578,-1.1115],[2.18933,0.52578,-1.12529],[2.18408,0.50784,-1.1115],[2.18408,0.66838,1.1115],[2.18408,0.66838,-1.1115],[2.18408,0.52578,1.1115],[2.18408,0.52578,-1.1115],[2.18408,0.52578,1.1115],[2.18408,0.66838,-1.1115],[2.3276,1.15,0.505],[2.3276,0.9936,0.505],[2.27884,1.15,0.505],[2.27884,0.9936,0.505],[2.27884,1.15,0.505],[2.3276,0.9936,0.505],[2.31564,1.1914,0.645],[2.31564,1.1914,0.92],[2.31564,0.7682,0.645],[2.31564,0.7682,0.92],[2.31564,0.7682,0.645],[2.31564,1.1914,0.92],[1.5341,0.4646,-0.408501],[1.5341,0.3358,-0.408501],[1.5341,0.4646,-0.8625],[1.5341,0.3358,-0.8625],[1.5341,0.4646,-0.8625],[1.5341,0.3358,-0.408501],[1.36074,0.306469,0.321545],[1.3271,0.375469,0.321545],[1.36074,0.313808,0.287501],[1.3271,0.382808,0.287501],[1.36074,0.313808,0.287501],[1.3271,0.375469,0.321545],[1.50046,0.306469,-0.253456],[1.50046,0.313808,-0.287501],[1.5341,0.375469,-0.253456],[1.5341,0.382808,-0.287501],[1.5341,0.375469,-0.253456],[1.50046,0.313808,-0.287501],[1.3271,0.3358,0.166501],[1.3271,0.375469,0.253456],[1.3271,0.4646,0.166501],[1.3271,0.504269,0.253456],[1.3271,0.4646,0.166501],[1.3271,0.375469,0.253456],[1.36074,0.2668,-0.408501],[1.36074,0.306469,-0.321545],[1.3271,0.3358,-0.408501],[1.3271,0.375469,-0.321545],[1.3271,0.3358,-0.408501],[1.36074,0.306469,-0.321545],[-1.2489,0.4646,-0.8625],[-1.2489,0.3358,-0.8625],[-1.4559,0.4646,-0.8625],[-1.4559,0.3358,-0.8625],[-1.4559,0.4646,-0.8625],[-1.2489,0.3358,-0.8625],[-1.2489,0.4646,0.8625],[-1.2489,0.3358,0.8625],[-1.2489,0.4646,0.408501],[-1.2489,0.3358,0.408501],[-1.2489,0.4646,0.408501],[-1.2489,0.3358,0.8625],[-1.28254,0.306469,0.321545],[-1.28254,0.313808,0.287501],[-1.2489,0.375469,0.321545],[-1.2489,0.382808,0.287501],[-1.2489,0.375469,0.321545],[-1.28254,0.313808,0.287501],[-1.2489,0.375469,-0.321545],[-1.2489,0.504269,-0.321545],[-1.2489,0.382808,-0.287501],[-1.2489,0.511608,-0.287501],[-1.2489,0.382808,-0.287501],[-1.2489,0.504269,-0.321545],[-1.2489,0.504269,-0.321545],[-1.2489,0.375469,-0.321545],[-1.2489,0.4646,-0.4085],[-1.2489,0.3358,-0.4085],[-1.2489,0.4646,-0.4085],[-1.2489,0.375469,-0.321545],[-1.2489,0.375469,0.253456],[-1.28254,0.306469,0.253456],[-1.2489,0.3358,0.166501],[-1.28254,0.2668,0.166501],[-1.2489,0.3358,0.166501],[-1.28254,0.306469,0.253456],[-2.3184,2.4104,0.6075],[-2.3184,2.4104,-0.6075],[-2.2954,2.4104,0.6075],[-2.2954,2.4104,-0.6075],[-2.2954,2.4104,0.6075],[-2.3184,2.4104,-0.6075],[-2.33266,0.7176,0.3146],[-2.34692,0.7176,0.3146],[-2.33266,0.7176,0.6075],[-2.34692,0.7176,0.6075],[-2.33266,0.7176,0.6075],[-2.34692,0.7176,0.3146],[1.0028,1.2696,1.055],[1.0028,1.2696,1.09],[1.0028,1.2972,1.055],[1.0028,1.2972,1.09],[1.0028,1.2972,1.055],[1.0028,1.2696,1.09],[0.138,2.34048,1.055],[0.138,2.34048,1.09],[0.138,1.6192,1.055],[0.138,1.6192,1.09],[0.138,1.6192,1.055],[0.138,2.34048,1.09],[0.138,0.4784,-1.055],[0.138,0.4784,-1.1],[0.138,1.2696,-1.055],[0.138,1.2696,-1.1],[0.138,1.2696,-1.055],[0.138,0.4784,-1.1],[1.0028,2.34048,-1.055],[1.0028,2.34048,-1.1],[1.0028,1.6192,-1.055],[1.0028,1.6192,-1.1],[1.0028,1.6192,-1.055],[1.0028,2.34048,-1.1],[1.53894,1.89212,-1.20364],[1.54247,1.887,-1.2074],[1.53895,1.89766,-1.21827],[1.54249,1.89224,-1.22148],[1.53895,1.89766,-1.21827],[1.54247,1.887,-1.2074],[1.54254,1.87378,-1.39013],[1.53896,1.87858,-1.39434],[1.54247,1.887,-1.3846],[1.53894,1.89212,-1.38836],[1.54247,1.887,-1.3846],[1.53896,1.87858,-1.39434],[1.53226,1.65232,-1.396],[1.5226,1.65232,-1.396],[1.53226,1.88048,-1.396],[1.5226,1.88048,-1.396],[1.53226,1.88048,-1.396],[1.5226,1.65232,-1.396],[1.58746,1.60503,-1.09],[1.58746,1.57504,-1.09],[1.58746,1.60503,-1.2524],[1.58746,1.57504,-1.285],[1.58746,1.60503,-1.2524],[1.58746,1.57504,-1.09],[2.3152,1.242,1.07653],[2.3152,0.414,1.07653],[2.3276,1.242,1.044],[2.3276,0.414,1.044],[2.3276,1.242,1.044],[2.3152,0.414,1.07653],[2.30849,1.26544,-1.08653],[2.32089,1.2654,-1.054],[2.3152,1.242,-1.08653],[2.3276,1.242,-1.054],[2.3152,1.242,-1.08653],[2.32089,1.2654,-1.054],[-2.306,0.414,-1.08653],[-2.306,1.2696,-1.08653],[-2.27608,0.414,-1.1],[-2.27608,1.2696,-1.1],[-2.27608,0.414,-1.1],[-2.306,1.2696,-1.08653],[-1.85152,0.414,-1.1],[-2.27608,0.414,-1.1],[-2.27608,1.2696,-1.1],[2.28528,0.414,-1.1],[1.95888,0.414,-1.1],[1.94156,0.56163,-1.1],[2.28528,1.242,-1.1],[2.28528,0.414,-1.1],[1.94156,0.56163,-1.1],[-2.27608,1.2696,-1.1],[0.138,1.2696,-1.1],[0.138,0.4784,-1.1],[-1.70125,0.69,-1.1],[-1.85152,0.414,-1.1],[-2.27608,1.2696,-1.1],[1.0028,0.8924,-1.1],[1.0028,1.242,-1.1],[2.28528,1.242,-1.1],[2.28528,1.242,-1.1],[1.94156,0.56163,-1.1],[1.88825,0.6992,-1.1],[2.28528,1.242,-1.1],[1.88825,0.6992,-1.1],[1.80345,0.817334,-1.1],[2.28528,1.242,-1.1],[1.80345,0.817334,-1.1],[1.69294,0.907981,-1.1],[2.28528,1.242,-1.1],[1.69294,0.907981,-1.1],[1.56425,0.964964,-1.1],[-1.68734,0.753372,-1.1],[-1.70125,0.69,-1.1],[-2.27608,1.2696,-1.1],[-1.64774,0.807097,-1.1],[-1.68734,0.753372,-1.1],[-2.27608,1.2696,-1.1],[-1.58847,0.842994,-1.1],[-1.64774,0.807097,-1.1],[-2.27608,1.2696,-1.1],[2.28528,1.242,-1.1],[1.56425,0.964964,-1.1],[1.42614,0.9844,-1.1],[-1.51855,0.8556,-1.1],[-1.58847,0.842994,-1.1],[-2.27608,1.2696,-1.1],[-1.17947,0.8556,-1.1],[-1.51855,0.8556,-1.1],[-2.27608,1.2696,-1.1],[-0.850175,0.414,-1.1],[-1.01076,0.69,-1.1],[-1.0236,0.753372,-1.1],[2.28528,1.242,-1.1],[1.42614,0.9844,-1.1],[1.28804,0.964964,-1.1],[-1.11491,0.842994,-1.1],[-1.17947,0.8556,-1.1],[-2.27608,1.2696,-1.1],[-1.11491,0.842994,-1.1],[-2.27608,1.2696,-1.1],[0.138,0.4784,-1.1],[-1.06017,0.807097,-1.1],[-1.11491,0.842994,-1.1],[0.138,0.4784,-1.1],[-1.0236,0.753372,-1.1],[-1.06017,0.807097,-1.1],[0.138,0.4784,-1.1],[-0.850175,0.414,-1.1],[-1.0236,0.753372,-1.1],[0.138,0.4784,-1.1],[0.857956,0.414,-1.1],[-0.850175,0.414,-1.1],[0.138,0.4784,-1.1],[0.857956,0.414,-1.1],[0.138,0.4784,-1.1],[0.7912,0.4784,-1.1],[1.0562,0.822,-1.1],[0.857956,0.414,-1.1],[0.7912,0.4784,-1.1],[1.0562,0.822,-1.1],[0.7912,0.4784,-1.1],[1.0028,0.8924,-1.1],[1.15996,0.911218,-1.1],[1.0562,0.822,-1.1],[1.0028,0.8924,-1.1],[1.28804,0.964964,-1.1],[1.15996,0.911218,-1.1],[1.0028,0.8924,-1.1],[1.0028,0.8924,-1.1],[2.28528,1.242,-1.1],[1.28804,0.964964,-1.1],[-2.3184,1.2696,-0.6075],[-2.3184,0.7176,-0.6075],[-2.2954,1.2696,-0.6075],[-2.2954,0.7176,-0.6075],[-2.2954,1.2696,-0.6075],[-2.3184,0.7176,-0.6075],[2.28528,0.414,-1.1],[2.3152,0.414,-1.08653],[2.29264,0.414,-1.062],[2.3276,0.414,-1.054],[2.29264,0.414,-1.062],[2.3152,0.414,-1.08653],[-2.31419,2.55597,1.0482],[-2.30419,2.55597,1.07451],[-2.30271,2.57293,1.04572],[-2.29541,2.57293,1.06494],[-2.30271,2.57293,1.04572],[-2.30419,2.55597,1.07451],[-2.3184,2.53519,-1.054],[-2.31419,2.55597,-1.0582],[-2.306,2.53519,-1.08653],[-2.30419,2.55597,-1.08451],[-2.306,2.53519,-1.08653],[-2.31419,2.55597,-1.0582],[1.28956,2.58437,1.04157],[1.30129,2.58437,1.04203],[1.30129,2.58437,-1.05203],[1.30129,2.58437,-1.05203],[1.28956,2.58437,-1.05157],[1.28956,2.58437,1.04157],[1.334,2.53519,-1.054],[1.3216,2.53519,-1.08653],[1.32979,2.55597,-1.0582],[1.31979,2.55597,-1.08451],[1.32979,2.55597,-1.0582],[1.3216,2.53519,-1.08653],[1.53226,1.63866,1.38985],[1.53226,1.633,1.375],[1.53894,1.64068,1.38836],[1.53895,1.63514,1.37373],[1.53894,1.64068,1.38836],[1.53226,1.633,1.375],[1.53226,1.8998,1.375],[1.53895,1.89766,1.37373],[1.53226,1.8998,1.217],[1.53895,1.89766,1.21827],[1.53226,1.8998,1.217],[1.53895,1.89766,1.37373],[1.53226,1.63866,1.20215],[1.53226,1.633,1.217],[1.5226,1.63866,1.20215],[1.5226,1.633,1.217],[1.5226,1.63866,1.20215],[1.53226,1.633,1.217],[1.55986,1.61123,1.26866],[1.58746,1.61123,1.26866],[1.55986,1.60503,1.2524],[1.58746,1.60503,1.2524],[1.55986,1.60503,1.2524],[1.58746,1.61123,1.26866],[1.58358,1.78644,1.308],[1.57421,1.79032,1.308],[1.58746,1.77707,1.308],[1.58746,1.5962,1.308],[1.58746,1.77707,1.308],[1.57421,1.79032,1.308],[1.57421,1.79032,1.308],[1.55986,1.79032,1.308],[1.58746,1.5962,1.308],[1.55986,1.5962,1.308],[1.58746,1.5962,1.308],[1.55986,1.79032,1.308],[-1.85152,0.414,-1.1],[-1.70125,0.69,-1.1],[-1.85152,0.414,-1.062],[-1.70125,0.69,-1.062],[-1.85152,0.414,-1.062],[-1.70125,0.69,-1.1],[0.857956,0.414,-1.1],[1.0562,0.822,-1.1],[0.857956,0.414,-1.062],[1.0562,0.822,-1.062],[0.857956,0.414,-1.062],[1.0562,0.822,-1.1],[1.88825,0.6992,1.09],[1.80345,0.817334,1.09],[1.88825,0.6992,1.052],[1.80345,0.817334,1.052],[1.88825,0.6992,1.052],[1.80345,0.817334,1.09],[-1.06017,0.807097,1.09],[-1.11491,0.842994,1.09],[-1.06017,0.807097,1.052],[-1.11491,0.842994,1.052],[-1.06017,0.807097,1.052],[-1.11491,0.842994,1.09],[1.94156,0.56163,1.052],[1.88825,0.6992,1.052],[1.97355,0.61683,1.052],[1.91704,0.7544,1.052],[1.97355,0.61683,1.052],[1.88825,0.6992,1.052],[-1.0236,0.753372,1.052],[-1.06017,0.807097,1.052],[-1.01051,0.808572,1.052],[-1.04855,0.862297,1.052],[-1.01051,0.808572,1.052],[-1.06017,0.807097,1.052],[-1.85152,0.414,-1.062],[-1.70125,0.69,-1.062],[-1.87155,0.4692,-1.062],[-1.71527,0.7452,-1.062],[-1.87155,0.4692,-1.062],[-1.70125,0.69,-1.062],[0.857956,0.414,-1.062],[1.0562,0.822,-1.062],[0.824928,0.4692,-1.062],[1.03507,0.8772,-1.062],[0.824928,0.4692,-1.062],[1.0562,0.822,-1.062],[2.29264,0.4692,1.052],[1.99191,0.4692,1.052],[2.25032,0.4692,0.54094],[1.99191,0.4692,0.54094],[2.25032,0.4692,0.54094],[1.99191,0.4692,1.052],[-0.830148,0.4692,1.052],[-0.997156,0.7452,1.052],[-0.830148,0.4692,0.54094],[-0.997156,0.7452,0.54094],[-0.830148,0.4692,0.54094],[-0.997156,0.7452,1.052],[-2.28344,0.4692,1.052],[-2.28344,0.4692,-1.062],[-2.24112,0.4692,0.54094],[-2.24112,0.4692,-0.55094],[-2.24112,0.4692,0.54094],[-2.28344,0.4692,-1.062],[-0.997156,0.7452,-1.062],[-0.830148,0.4692,-1.062],[-0.997156,0.7452,-0.55094],[-0.830148,0.4692,-0.55094],[-0.997156,0.7452,-0.55094],[-0.830148,0.4692,-1.062],[0.824928,0.4692,0.54094],[1.99191,0.4692,0.54094],[1.03507,0.8772,0.54094],[1.14505,0.966418,0.54094],[1.03507,0.8772,0.54094],[1.99191,0.4692,0.54094],[1.99191,0.4692,0.54094],[1.97355,0.61683,0.54094],[1.14505,0.966418,0.54094],[1.28081,1.02016,0.54094],[1.14505,0.966418,0.54094],[1.97355,0.61683,0.54094],[1.97355,0.61683,0.54094],[1.91704,0.7544,0.54094],[1.28081,1.02016,0.54094],[1.42721,1.0396,0.54094],[1.28081,1.02016,0.54094],[1.91704,0.7544,0.54094],[1.91704,0.7544,0.54094],[1.82716,0.872534,0.54094],[1.42721,1.0396,0.54094],[1.5736,1.02016,0.54094],[1.42721,1.0396,0.54094],[1.82716,0.872534,0.54094],[1.82716,0.872534,0.54094],[1.71001,0.963181,0.54094],[1.5736,1.02016,0.54094],[-2.43091,0.38864,-1.12699],[-2.18289,0.38864,-1.12699],[-2.42748,0.371258,-1.12347],[-2.18632,0.371258,-1.12347],[-2.42748,0.371258,-1.12347],[-2.18289,0.38864,-1.12699],[-2.42329,0.470631,-1.11838],[-2.19051,0.470631,-1.11838],[-2.42292,0.465213,-1.11797],[-2.19088,0.465213,-1.11797],[-2.42292,0.465213,-1.11797],[-2.19051,0.470631,-1.11838],[-2.42292,0.465213,1.11797],[-2.42329,0.470631,1.11838],[-2.42329,0.470631,-1.11838],[-2.42329,0.470631,-1.11838],[-2.42292,0.465213,-1.11797],[-2.42292,0.465213,1.11797],[-2.18289,0.38864,1.12699],[-2.18944,0.41195,1.1196],[-2.43091,0.38864,1.12699],[-2.42436,0.41195,1.1196],[-2.43091,0.38864,1.12699],[-2.18944,0.41195,1.1196],[2.18408,0.52578,-1.1115],[2.18408,0.50784,-1.1115],[2.18408,0.50784,1.1115],[2.18408,0.50784,1.1115],[2.18408,0.52578,1.1115],[2.18408,0.52578,-1.1115],[2.20202,0.391655,-1.12529],[2.20202,0.40434,-1.131],[2.31799,0.391744,-1.12529],[2.32565,0.402713,-1.131],[2.31799,0.391744,-1.12529],[2.20202,0.40434,-1.131],[2.20202,0.391655,1.12529],[2.31799,0.391744,1.12529],[2.20202,0.40434,1.131],[2.32565,0.402713,1.131],[2.20202,0.40434,1.131],[2.31799,0.391744,1.12529],[2.42006,0.68632,-1.1115],[2.43275,0.681065,-1.1115],[2.42006,0.681065,-1.12529],[2.43099,0.679314,-1.12338],[2.42006,0.681065,-1.12529],[2.43275,0.681065,-1.1115],[2.18408,0.66838,1.1115],[2.18933,0.66838,1.12529],[2.18933,0.681065,1.1115],[2.19109,0.679314,1.12338],[2.18933,0.681065,1.1115],[2.18933,0.66838,1.12529],[2.31518,0.3864,1.1115],[2.33149,0.389956,1.1115],[2.31799,0.391744,1.12529],[2.33066,0.391774,1.12382],[2.31799,0.391744,1.12529],[2.33149,0.389956,1.1115],[2.20202,0.68632,1.1115],[2.42006,0.68632,1.1115],[2.20202,0.68632,-1.1115],[2.42006,0.68632,-1.1115],[2.20202,0.68632,-1.1115],[2.42006,0.68632,1.1115],[2.3276,1.15,-0.505],[2.3276,1.15,0.505],[2.27884,1.15,-0.505],[2.27884,1.15,0.505],[2.27884,1.15,-0.505],[2.3276,1.15,0.505],[2.3276,1.1914,0.92],[2.3276,0.7682,0.92],[2.31564,1.1914,0.92],[2.31564,0.7682,0.92],[2.31564,1.1914,0.92],[2.3276,0.7682,0.92],[1.3271,0.3358,0.8625],[1.3271,0.4646,0.8625],[1.3271,0.3358,0.408501],[1.3271,0.4646,0.408501],[1.3271,0.3358,0.408501],[1.3271,0.4646,0.8625],[1.50046,0.2668,0.8625],[1.36074,0.2668,0.8625],[1.50046,0.2668,0.408501],[1.36074,0.2668,0.408501],[1.50046,0.2668,0.408501],[1.36074,0.2668,0.8625],[1.3271,0.511608,-0.287501],[1.3271,0.382808,-0.287501],[1.3271,0.504269,-0.253456],[1.3271,0.375469,-0.253456],[1.3271,0.504269,-0.253456],[1.3271,0.382808,-0.287501],[1.5341,0.511608,0.287501],[1.5341,0.504269,0.321545],[1.5341,0.382808,0.287501],[1.5341,0.375469,0.321545],[1.5341,0.382808,0.287501],[1.5341,0.504269,0.321545],[1.5341,0.375469,-0.253456],[1.5341,0.3358,-0.166501],[1.50046,0.306469,-0.253456],[1.50046,0.2668,-0.166501],[1.50046,0.306469,-0.253456],[1.5341,0.3358,-0.166501],[1.50046,0.306469,-0.253456],[1.50046,0.2668,-0.166501],[1.36074,0.306469,-0.253456],[1.36074,0.2668,-0.166501],[1.36074,0.306469,-0.253456],[1.50046,0.2668,-0.166501],[-1.2489,0.4646,0.166501],[-1.2489,0.3358,0.166501],[-1.2489,0.4646,-0.1665],[-1.2489,0.3358,-0.1665],[-1.2489,0.4646,-0.1665],[-1.2489,0.3358,0.166501],[-1.4559,0.382808,0.287501],[-1.4559,0.375469,0.253456],[-1.42226,0.313808,0.287501],[-1.42226,0.306469,0.253456],[-1.42226,0.313808,0.287501],[-1.4559,0.375469,0.253456],[-1.2489,0.382808,-0.287501],[-1.28254,0.313808,-0.287501],[-1.2489,0.375469,-0.321545],[-1.28254,0.306469,-0.321545],[-1.2489,0.375469,-0.321545],[-1.28254,0.313808,-0.287501],[-1.4559,0.3358,0.408501],[-1.4559,0.4646,0.408501],[-1.4559,0.375469,0.321545],[-1.4559,0.504269,0.321545],[-1.4559,0.375469,0.321545],[-1.4559,0.4646,0.408501],[-1.42226,0.2668,-0.1665],[-1.4559,0.3358,-0.1665],[-1.42226,0.306469,-0.253456],[-1.4559,0.375469,-0.253456],[-1.42226,0.306469,-0.253456],[-1.4559,0.3358,-0.1665],[-2.3184,1.6192,0.6075],[-2.3184,2.4104,0.6075],[-2.2954,1.6192,0.6075],[-2.2954,2.4104,0.6075],[-2.2954,1.6192,0.6075],[-2.3184,2.4104,0.6075],[-2.33266,0.7176,-0.6075],[-2.34692,0.7176,-0.6075],[-2.33266,0.7176,-0.3146],[-2.34692,0.7176,-0.3146],[-2.33266,0.7176,-0.3146],[-2.34692,0.7176,-0.6075],[1.0028,1.242,1.055],[1.0028,1.242,1.09],[1.0028,1.2696,1.055],[1.0028,1.2696,1.09],[1.0028,1.2696,1.055],[1.0028,1.242,1.09],[0.142622,2.36372,1.09],[0.138,2.34048,1.09],[0.142622,2.36372,1.055],[0.138,2.34048,1.055],[0.142622,2.36372,1.055],[0.138,2.34048,1.09],[0.138,0.4784,-1.055],[0.7912,0.4784,-1.055],[0.138,0.4784,-1.1],[0.7912,0.4784,-1.1],[0.138,0.4784,-1.1],[0.7912,0.4784,-1.055],[0.998178,2.36372,-1.055],[0.998178,2.36372,-1.1],[1.0028,2.34048,-1.055],[1.0028,2.34048,-1.1],[1.0028,2.34048,-1.055],[0.998178,2.36372,-1.1],[1.53226,1.89414,-1.20215],[1.53226,1.88048,-1.196],[1.53894,1.89212,-1.20364],[1.53896,1.87858,-1.19766],[1.53894,1.89212,-1.20364],[1.53226,1.88048,-1.196],[1.53226,1.63866,-1.20215],[1.53894,1.64068,-1.20364],[1.53226,1.65232,-1.196],[1.53896,1.65421,-1.19766],[1.53226,1.65232,-1.196],[1.53894,1.64068,-1.20364],[1.53226,1.8998,-1.375],[1.5226,1.8998,-1.375],[1.53226,1.8998,-1.217],[1.5226,1.8998,-1.217],[1.53226,1.8998,-1.217],[1.5226,1.8998,-1.375],[1.55986,1.57504,-1.09],[1.55986,1.57504,-1.285],[1.58746,1.57504,-1.09],[1.58746,1.57504,-1.285],[1.58746,1.57504,-1.09],[1.55986,1.57504,-1.285],[2.28528,1.242,1.09],[2.28528,0.414,1.09],[2.3152,1.242,1.07653],[2.3152,0.414,1.07653],[2.3152,1.242,1.07653],[2.28528,0.414,1.09],[2.3152,1.242,-1.08653],[2.3152,0.414,-1.08653],[2.28528,1.242,-1.1],[2.28528,0.414,-1.1],[2.28528,1.242,-1.1],[2.3152,0.414,-1.08653],[-2.3184,0.414,-1.054],[-2.3184,1.2696,-1.054],[-2.306,0.414,-1.08653],[-2.306,1.2696,-1.08653],[-2.306,0.414,-1.08653],[-2.3184,1.2696,-1.054],[-2.27608,1.2696,1.09],[-2.27608,0.414,1.09],[-1.85152,0.414,1.09],[1.94156,0.56163,1.09],[1.95888,0.414,1.09],[2.28528,0.414,1.09],[1.94156,0.56163,1.09],[2.28528,0.414,1.09],[2.28528,1.242,1.09],[0.138,0.4784,1.09],[0.138,1.2696,1.09],[-2.27608,1.2696,1.09],[-2.27608,1.2696,1.09],[-1.85152,0.414,1.09],[-1.70125,0.69,1.09],[2.28528,1.242,1.09],[1.0028,1.242,1.09],[1.0028,0.8924,1.09],[1.88825,0.6992,1.09],[1.94156,0.56163,1.09],[2.28528,1.242,1.09],[1.80345,0.817334,1.09],[1.88825,0.6992,1.09],[2.28528,1.242,1.09],[1.69294,0.907981,1.09],[1.80345,0.817334,1.09],[2.28528,1.242,1.09],[1.56425,0.964964,1.09],[1.69294,0.907981,1.09],[2.28528,1.242,1.09],[-2.27608,1.2696,1.09],[-1.70125,0.69,1.09],[-1.68734,0.753372,1.09],[-2.27608,1.2696,1.09],[-1.68734,0.753372,1.09],[-1.64774,0.807097,1.09],[-2.27608,1.2696,1.09],[-1.64774,0.807097,1.09],[-1.58847,0.842994,1.09],[1.42614,0.9844,1.09],[1.56425,0.964964,1.09],[2.28528,1.242,1.09],[-2.27608,1.2696,1.09],[-1.58847,0.842994,1.09],[-1.51855,0.8556,1.09],[-2.27608,1.2696,1.09],[-1.51855,0.8556,1.09],[-1.17947,0.8556,1.09],[-1.0236,0.753372,1.09],[-1.01076,0.69,1.09],[-0.850175,0.414,1.09],[1.28804,0.964964,1.09],[1.42614,0.9844,1.09],[2.28528,1.242,1.09],[-2.27608,1.2696,1.09],[-1.17947,0.8556,1.09],[-1.11491,0.842994,1.09],[0.138,0.4784,1.09],[-2.27608,1.2696,1.09],[-1.11491,0.842994,1.09],[0.138,0.4784,1.09],[-1.11491,0.842994,1.09],[-1.06017,0.807097,1.09],[0.138,0.4784,1.09],[-1.06017,0.807097,1.09],[-1.0236,0.753372,1.09],[0.138,0.4784,1.09],[-1.0236,0.753372,1.09],[-0.850175,0.414,1.09],[0.138,0.4784,1.09],[-0.850175,0.414,1.09],[0.857956,0.414,1.09],[0.7912,0.4784,1.09],[0.138,0.4784,1.09],[0.857956,0.414,1.09],[0.7912,0.4784,1.09],[0.857956,0.414,1.09],[1.0562,0.822,1.09],[1.0028,0.8924,1.09],[0.7912,0.4784,1.09],[1.0562,0.822,1.09],[1.0028,0.8924,1.09],[1.0562,0.822,1.09],[1.15996,0.911218,1.09],[1.0028,0.8924,1.09],[1.15996,0.911218,1.09],[1.28804,0.964964,1.09],[1.28804,0.964964,1.09],[2.28528,1.242,1.09],[1.0028,0.8924,1.09],[-2.306,1.6192,-1.08653],[-2.3184,1.6192,-1.054],[-2.306,2.53519,-1.08653],[-2.3184,2.53519,-1.054],[-2.306,2.53519,-1.08653],[-2.3184,1.6192,-1.054],[-2.3184,0.414,-1.054],[-2.306,0.414,-1.08653],[-2.28344,0.414,-1.062],[-2.27608,0.414,-1.1],[-2.28344,0.414,-1.062],[-2.306,0.414,-1.08653],[-2.30419,2.55597,-1.08451],[-2.29541,2.57293,-1.07494],[-2.28005,2.55597,-1.09541],[-2.27777,2.57293,-1.0829],[-2.28005,2.55597,-1.09541],[-2.29541,2.57293,-1.07494],[1.31979,2.55597,-1.08451],[1.31101,2.57293,-1.07494],[1.32979,2.55597,-1.0582],[1.31831,2.57293,-1.05572],[1.32979,2.55597,-1.0582],[1.31101,2.57293,-1.07494],[1.28956,2.58437,-1.05157],[1.30129,2.58437,-1.05203],[1.28998,2.58437,-1.06435],[1.29798,2.58437,-1.06074],[1.28998,2.58437,-1.06435],[1.30129,2.58437,-1.05203],[1.61,1.6192,-1.054],[1.5976,1.6192,-1.08653],[1.334,2.53519,-1.054],[1.3216,2.53519,-1.08653],[1.334,2.53519,-1.054],[1.5976,1.6192,-1.08653],[1.60861,1.6238,0.9605],[1.58887,1.68934,1.035],[1.61,1.6192,1.044],[1.3813,2.3782,0.0795],[1.40105,2.31266,0.005],[1.40105,2.31266,-0.005],[1.3813,2.3782,0.0795],[1.40105,2.31266,-0.005],[1.3813,2.3782,-0.0795],[1.58887,1.68934,-0.005],[1.40105,2.31266,-0.005],[1.40105,2.31266,0.005],[1.58887,1.68934,-0.005],[1.40105,2.31266,0.005],[1.58887,1.68934,0.005],[1.60861,1.6238,-0.0795],[1.58887,1.68934,-0.005],[1.58887,1.68934,0.005],[1.60861,1.6238,-0.0795],[1.58887,1.68934,0.005],[1.60861,1.6238,0.0795],[1.61,1.6192,1.044],[1.58887,1.68934,1.035],[1.40105,2.31266,1.035],[1.334,2.53519,1.044],[1.61,1.6192,1.044],[1.40105,2.31266,1.035],[1.334,2.53519,1.044],[1.40105,2.31266,1.035],[1.3813,2.3782,0.9605],[1.334,2.53519,1.044],[1.3813,2.3782,0.9605],[1.3813,2.3782,0.0795],[1.334,2.53519,-1.054],[1.334,2.53519,1.044],[1.3813,2.3782,0.0795],[1.334,2.53519,-1.054],[1.3813,2.3782,0.0795],[1.3813,2.3782,-0.0795],[1.334,2.53519,-1.054],[1.3813,2.3782,-0.0795],[1.3813,2.3782,-0.9605],[1.334,2.53519,-1.054],[1.3813,2.3782,-0.9605],[1.40105,2.31266,-1.035],[1.61,1.6192,-1.054],[1.334,2.53519,-1.054],[1.40105,2.31266,-1.035],[1.61,1.6192,-1.054],[1.40105,2.31266,-1.035],[1.58887,1.68934,-1.035],[1.61,1.6192,-1.054],[1.58887,1.68934,-1.035],[1.60861,1.6238,-0.9605],[1.60861,1.6238,0.0795],[1.60861,1.6238,0.9605],[1.61,1.6192,1.044],[1.61,1.6192,-1.054],[1.60861,1.6238,-0.9605],[1.60861,1.6238,-0.0795],[1.61,1.6192,1.044],[1.61,1.6192,-1.054],[1.60861,1.6238,-0.0795],[1.60861,1.6238,-0.0795],[1.60861,1.6238,0.0795],[1.61,1.6192,1.044],[1.53894,1.64068,1.38836],[1.53895,1.63514,1.37373],[1.54247,1.6458,1.3846],[1.54249,1.64056,1.37052],[1.54247,1.6458,1.3846],[1.53895,1.63514,1.37373],[1.53895,1.89766,1.37373],[1.54249,1.89224,1.37052],[1.53895,1.89766,1.21827],[1.54249,1.89224,1.22148],[1.53895,1.89766,1.21827],[1.54249,1.89224,1.37052],[1.53226,1.65232,1.196],[1.53226,1.63866,1.20215],[1.5226,1.65232,1.196],[1.5226,1.63866,1.20215],[1.5226,1.65232,1.196],[1.53226,1.63866,1.20215],[1.55986,1.62619,1.2754],[1.58746,1.62619,1.2754],[1.55986,1.61123,1.26866],[1.58746,1.61123,1.26866],[1.55986,1.61123,1.26866],[1.58746,1.62619,1.2754],[1.58358,1.78644,1.308],[1.58358,1.78644,1.2754],[1.57421,1.79032,1.308],[1.57421,1.79032,1.2754],[1.57421,1.79032,1.308],[1.58358,1.78644,1.2754],[0.138,1.2696,1.055],[0.138,0.4784,1.055],[0.7912,0.4784,1.055],[0.1932,2.22272,1.055],[0.27968,2.3092,1.055],[0.19872,2.4012,1.055],[0.138,1.2696,1.055],[0.7912,0.4784,1.055],[1.0028,0.8924,1.055],[0.138,1.2696,1.055],[1.0028,0.8924,1.055],[1.0028,1.242,1.055],[0.138,1.6192,1.055],[0.138,1.2696,1.055],[1.0028,1.242,1.055],[0.1932,2.22272,1.055],[0.19872,2.4012,1.055],[0.175483,2.39658,1.055],[0.19872,2.4012,1.055],[0.27968,2.3092,1.055],[0.81972,2.3092,1.055],[0.94208,2.4012,1.055],[0.19872,2.4012,1.055],[0.81972,2.3092,1.055],[0.94208,2.4012,1.055],[0.81972,2.3092,1.055],[0.9062,2.22272,1.055],[0.965316,2.39658,1.055],[0.94208,2.4012,1.055],[0.9062,2.22272,1.055],[0.985015,2.38342,1.055],[0.965316,2.39658,1.055],[0.9062,2.22272,1.055],[0.1932,2.22272,1.055],[0.175483,2.39658,1.055],[0.155784,2.38342,1.055],[0.998178,2.36372,1.055],[0.985015,2.38342,1.055],[0.9062,2.22272,1.055],[0.998178,2.36372,1.055],[0.9062,2.22272,1.055],[0.9062,1.68268,1.055],[0.1932,2.22272,1.055],[0.155784,2.38342,1.055],[0.142622,2.36372,1.055],[0.1932,2.22272,1.055],[0.142622,2.36372,1.055],[0.138,2.34048,1.055],[0.1932,2.22272,1.055],[0.138,2.34048,1.055],[0.138,1.6192,1.055],[0.1932,1.68268,1.055],[0.1932,2.22272,1.055],[0.138,1.6192,1.055],[0.27968,1.5962,1.055],[0.1932,1.68268,1.055],[0.138,1.6192,1.055],[0.27968,1.5962,1.055],[0.138,1.6192,1.055],[1.0028,1.242,1.055],[0.81972,1.5962,1.055],[0.27968,1.5962,1.055],[1.0028,1.242,1.055],[0.9062,1.68268,1.055],[0.81972,1.5962,1.055],[1.0028,1.242,1.055],[0.998178,2.36372,1.055],[0.9062,1.68268,1.055],[1.0028,1.242,1.055],[1.0028,1.6192,1.055],[1.0028,2.34048,1.055],[0.998178,2.36372,1.055],[1.0028,1.2972,1.055],[1.0028,1.6192,1.055],[0.998178,2.36372,1.055],[1.0028,1.2696,1.055],[1.0028,1.2972,1.055],[0.998178,2.36372,1.055],[1.0028,1.242,1.055],[1.0028,1.2696,1.055],[0.998178,2.36372,1.055],[-1.68734,0.753372,1.09],[-1.70125,0.69,1.09],[-1.68734,0.753372,1.052],[-1.70125,0.69,1.052],[-1.68734,0.753372,1.052],[-1.70125,0.69,1.09],[-1.01076,0.69,-1.1],[-0.850175,0.414,-1.1],[-1.01076,0.69,-1.062],[-0.850175,0.414,-1.062],[-1.01076,0.69,-1.062],[-0.850175,0.414,-1.1],[1.95888,0.414,1.09],[1.94156,0.56163,1.09],[1.95888,0.414,1.052],[1.94156,0.56163,1.052],[1.95888,0.414,1.052],[1.94156,0.56163,1.09],[-1.01076,0.69,1.09],[-1.0236,0.753372,1.09],[-1.01076,0.69,1.052],[-1.0236,0.753372,1.052],[-1.01076,0.69,1.052],[-1.0236,0.753372,1.09],[2.29264,0.4692,1.052],[2.29264,0.414,1.052],[1.99191,0.4692,1.052],[1.95888,0.414,1.052],[1.99191,0.4692,1.052],[2.29264,0.414,1.052],[-0.850175,0.414,1.052],[-1.01076,0.69,1.052],[-0.830148,0.4692,1.052],[-0.997156,0.7452,1.052],[-0.830148,0.4692,1.052],[-1.01076,0.69,1.052],[-2.28344,0.414,1.052],[-2.28344,0.414,-1.062],[-2.28344,0.4692,1.052],[-2.28344,0.4692,-1.062],[-2.28344,0.4692,1.052],[-2.28344,0.414,-1.062],[-0.850175,0.414,-1.062],[-0.830148,0.4692,-1.062],[-1.01076,0.69,-1.062],[-0.997156,0.7452,-1.062],[-1.01076,0.69,-1.062],[-0.830148,0.4692,-1.062],[2.29264,0.4692,-1.062],[2.25032,0.4692,-0.55094],[1.99191,0.4692,-1.062],[1.99191,0.4692,-0.55094],[1.99191,0.4692,-1.062],[2.25032,0.4692,-0.55094],[1.03507,0.8772,1.052],[0.824928,0.4692,1.052],[1.03507,0.8772,0.54094],[0.824928,0.4692,0.54094],[1.03507,0.8772,0.54094],[0.824928,0.4692,1.052],[-1.71527,0.7452,1.052],[-1.87155,0.4692,1.052],[-1.71527,0.7452,0.54094],[-1.87155,0.4692,0.54094],[-1.71527,0.7452,0.54094],[-1.87155,0.4692,1.052],[-1.04855,0.862297,-1.062],[-1.01051,0.808572,-1.062],[-1.04855,0.862297,-0.55094],[-1.01051,0.808572,-0.55094],[-1.04855,0.862297,-0.55094],[-1.01051,0.808572,-1.062],[0.824928,0.4692,-0.55094],[1.99191,0.4692,-0.55094],[0.824928,0.4692,0.54094],[1.99191,0.4692,0.54094],[0.824928,0.4692,0.54094],[1.99191,0.4692,-0.55094],[-2.18632,0.371258,-1.12347],[-2.18289,0.38864,-1.12699],[-2.18289,0.38864,1.12699],[-2.18289,0.38864,1.12699],[-2.18632,0.371258,1.12347],[-2.18632,0.371258,-1.12347],[-2.42329,0.417169,-1.11838],[-2.42436,0.41195,-1.1196],[-2.42436,0.41195,1.1196],[-2.42436,0.41195,1.1196],[-2.42329,0.417169,1.11838],[-2.42329,0.417169,-1.11838],[-2.19051,0.417169,1.11838],[-2.18944,0.41195,1.1196],[-2.18944,0.41195,-1.1196],[-2.18944,0.41195,-1.1196],[-2.19051,0.417169,-1.11838],[-2.19051,0.417169,1.11838],[-2.18289,0.38864,-1.12699],[-2.18944,0.41195,-1.1196],[-2.18289,0.38864,1.12699],[-2.18944,0.41195,1.1196],[-2.18289,0.38864,1.12699],[-2.18944,0.41195,-1.1196],[2.18933,0.40434,1.12529],[2.18933,0.4899,1.12529],[2.18408,0.40434,1.1115],[2.18408,0.4899,1.1115],[2.18408,0.40434,1.1115],[2.18933,0.4899,1.12529],[2.18933,0.40434,-1.12529],[2.18933,0.4899,-1.12529],[2.20202,0.40434,-1.131],[2.20202,0.4899,-1.131],[2.20202,0.40434,-1.131],[2.18933,0.4899,-1.12529],[2.32565,0.402713,-1.131],[2.40833,0.494263,-1.131],[2.33897,0.401345,-1.12529],[2.42102,0.494263,-1.12529],[2.33897,0.401345,-1.12529],[2.40833,0.494263,-1.131],[2.20202,0.68632,-1.1115],[2.20202,0.681065,-1.12529],[2.18933,0.681065,-1.1115],[2.19109,0.679314,-1.12338],[2.18933,0.681065,-1.1115],[2.20202,0.681065,-1.12529],[2.18408,0.40434,1.1115],[2.18933,0.391655,1.1115],[2.18933,0.40434,1.12529],[2.19109,0.393406,1.12338],[2.18933,0.40434,1.12529],[2.18933,0.391655,1.1115],[2.34485,0.399977,1.1115],[2.33897,0.401345,1.12529],[2.33149,0.389956,1.1115],[2.33066,0.391774,1.12382],[2.33149,0.389956,1.1115],[2.33897,0.401345,1.12529],[2.42627,0.494263,1.1115],[2.34485,0.399977,1.1115],[2.42627,0.494263,-1.1115],[2.34485,0.399977,-1.1115],[2.42627,0.494263,-1.1115],[2.34485,0.399977,1.1115],[2.3276,0.7406,0.505],[2.3276,0.7406,-0.505],[2.27884,0.7406,0.505],[2.27884,0.7406,-0.505],[2.27884,0.7406,0.505],[2.3276,0.7406,-0.505],[2.3276,0.7682,0.645],[2.3276,1.1914,0.645],[2.31564,0.7682,0.645],[2.31564,1.1914,0.645],[2.31564,0.7682,0.645],[2.3276,1.1914,0.645],[1.3271,0.3358,-0.408501],[1.3271,0.4646,-0.408501],[1.3271,0.3358,-0.8625],[1.3271,0.4646,-0.8625],[1.3271,0.3358,-0.8625],[1.3271,0.4646,-0.408501],[1.50046,0.2668,-0.408501],[1.36074,0.2668,-0.408501],[1.50046,0.2668,-0.8625],[1.36074,0.2668,-0.8625],[1.50046,0.2668,-0.8625],[1.36074,0.2668,-0.408501],[1.50046,0.313808,-0.287501],[1.36074,0.313808,-0.287501],[1.50046,0.306469,-0.321545],[1.36074,0.306469,-0.321545],[1.50046,0.306469,-0.321545],[1.36074,0.313808,-0.287501],[1.5341,0.511608,0.287501],[1.3271,0.511608,0.287501],[1.5341,0.504269,0.321545],[1.3271,0.504269,0.321545],[1.5341,0.504269,0.321545],[1.3271,0.511608,0.287501],[1.3271,0.504269,-0.253456],[1.3271,0.4646,-0.166501],[1.5341,0.504269,-0.253456],[1.5341,0.4646,-0.166501],[1.5341,0.504269,-0.253456],[1.3271,0.4646,-0.166501],[1.3271,0.504269,0.321545],[1.3271,0.4646,0.408501],[1.5341,0.504269,0.321545],[1.5341,0.4646,0.408501],[1.5341,0.504269,0.321545],[1.3271,0.4646,0.408501],[-1.4559,0.3358,0.8625],[-1.4559,0.4646,0.8625],[-1.4559,0.3358,0.408501],[-1.4559,0.4646,0.408501],[-1.4559,0.3358,0.408501],[-1.4559,0.4646,0.8625],[-1.28254,0.2668,0.8625],[-1.42226,0.2668,0.8625],[-1.28254,0.2668,0.408501],[-1.42226,0.2668,0.408501],[-1.28254,0.2668,0.408501],[-1.42226,0.2668,0.8625],[-1.4559,0.511608,-0.287501],[-1.4559,0.382808,-0.287501],[-1.4559,0.504269,-0.253456],[-1.4559,0.375469,-0.253456],[-1.4559,0.504269,-0.253456],[-1.4559,0.382808,-0.287501],[-1.2489,0.511608,0.287501],[-1.2489,0.504269,0.321545],[-1.2489,0.382808,0.287501],[-1.2489,0.375469,0.321545],[-1.2489,0.382808,0.287501],[-1.2489,0.504269,0.321545],[-1.2489,0.375469,-0.253456],[-1.2489,0.3358,-0.1665],[-1.28254,0.306469,-0.253456],[-1.28254,0.2668,-0.1665],[-1.28254,0.306469,-0.253456],[-1.2489,0.3358,-0.1665],[-1.28254,0.306469,-0.253456],[-1.28254,0.2668,-0.1665],[-1.42226,0.306469,-0.253456],[-1.42226,0.2668,-0.1665],[-1.42226,0.306469,-0.253456],[-1.28254,0.2668,-0.1665],[-2.33266,0.5244,-0.6075],[-2.34692,0.5244,-0.6075],[-2.33266,0.7176,-0.6075],[-2.34692,0.7176,-0.6075],[-2.33266,0.7176,-0.6075],[-2.34692,0.5244,-0.6075],[-2.34692,0.5244,-0.3146],[-2.33266,0.5244,-0.3146],[-2.34692,0.7176,-0.3146],[-2.33266,0.7176,-0.3146],[-2.34692,0.7176,-0.3146],[-2.33266,0.5244,-0.3146],[0.175483,2.39658,1.09],[0.155784,2.38342,1.09],[0.175483,2.39658,1.055],[0.155784,2.38342,1.055],[0.175483,2.39658,1.055],[0.155784,2.38342,1.09],[1.0028,1.242,-1.055],[1.0028,1.242,-1.1],[1.0028,0.8924,-1.055],[1.0028,0.8924,-1.1],[1.0028,0.8924,-1.055],[1.0028,1.242,-1.1],[0.965316,2.39658,-1.055],[0.965317,2.39658,-1.1],[0.985015,2.38342,-1.055],[0.985016,2.38342,-1.1],[0.985015,2.38342,-1.055],[0.965317,2.39658,-1.1],[1.53226,1.88048,-1.196],[1.53226,1.65232,-1.196],[1.53896,1.87858,-1.19766],[1.53896,1.65421,-1.19766],[1.53896,1.87858,-1.19766],[1.53226,1.65232,-1.196],[1.53226,1.89414,-1.38985],[1.53226,1.8998,-1.375],[1.53894,1.89212,-1.38836],[1.53895,1.89766,-1.37373],[1.53894,1.89212,-1.38836],[1.53226,1.8998,-1.375],[1.5226,1.89414,-1.20215],[1.53226,1.89414,-1.20215],[1.5226,1.8998,-1.217],[1.53226,1.8998,-1.217],[1.5226,1.8998,-1.217],[1.53226,1.89414,-1.20215],[1.55986,1.62619,-1.2754],[1.55986,1.5962,-1.308],[1.55986,1.61123,-1.26866],[1.55986,1.58124,-1.30126],[1.55986,1.61123,-1.26866],[1.55986,1.5962,-1.308],[2.30849,1.26544,1.07653],[2.29039,1.28184,1.07653],[2.27857,1.26554,1.09],[2.26047,1.28219,1.09],[2.27857,1.26554,1.09],[2.29039,1.28184,1.07653],[1.56768,1.6192,-1.1],[1.5976,1.6192,-1.08653],[2.26047,1.28219,-1.1],[2.29039,1.28184,-1.08653],[2.26047,1.28219,-1.1],[1.5976,1.6192,-1.08653],[2.30279,1.28169,1.044],[2.32089,1.2654,1.044],[2.30279,1.28169,-1.054],[2.32089,1.2654,-1.054],[2.30279,1.28169,-1.054],[2.32089,1.2654,1.044],[0.138,1.6192,-1.1],[0.138,1.2696,-1.1],[-2.27608,1.6192,-1.1],[-2.27608,1.2696,-1.1],[-2.27608,1.6192,-1.1],[0.138,1.2696,-1.1],[-2.3184,1.6192,-1.054],[-2.3184,1.6192,-0.6075],[-2.3184,2.4104,-0.6075],[-2.3184,2.53519,-1.054],[-2.3184,1.6192,-1.054],[-2.3184,2.4104,-0.6075],[-2.3184,2.4104,0.6075],[-2.3184,1.6192,0.6075],[-2.3184,1.6192,1.044],[-2.3184,2.4104,0.6075],[-2.3184,1.6192,1.044],[-2.3184,2.53519,1.044],[-2.3184,2.4104,0.6075],[-2.3184,2.53519,1.044],[-2.3184,2.53519,-1.054],[-2.3184,2.53519,-1.054],[-2.3184,2.4104,-0.6075],[-2.3184,2.4104,0.6075],[-2.3184,0.7176,-0.6075],[-2.3184,1.2696,-0.6075],[-2.3184,1.2696,-1.054],[-2.3184,0.7176,-0.6075],[-2.3184,1.2696,-1.054],[-2.3184,0.414,-1.054],[-2.3184,1.2696,1.044],[-2.3184,1.2696,0.6075],[-2.3184,0.7176,0.6075],[-2.3184,0.414,1.044],[-2.3184,1.2696,1.044],[-2.3184,0.7176,0.6075],[-2.3184,0.414,-1.054],[-2.3184,0.414,1.044],[-2.3184,0.7176,0.6075],[-2.3184,0.7176,0.6075],[-2.3184,0.7176,-0.6075],[-2.3184,0.414,-1.054],[-2.30271,2.57293,-1.05572],[-2.28569,2.58437,-1.05203],[-2.29541,2.57293,-1.07494],[-2.28238,2.58437,-1.06074],[-2.29541,2.57293,-1.07494],[-2.28569,2.58437,-1.05203],[1.31831,2.57293,1.04572],[1.32979,2.55597,1.0482],[1.32979,2.55597,-1.0582],[1.32979,2.55597,-1.0582],[1.31831,2.57293,-1.05572],[1.31831,2.57293,1.04572],[1.28956,2.58437,-1.05157],[1.28998,2.58437,-1.06435],[-2.27438,2.58437,-1.06435],[-2.27438,2.58437,-1.06435],[-2.27396,2.58437,-1.05157],[1.28956,2.58437,-1.05157],[1.334,2.53519,1.044],[1.32979,2.55597,1.0482],[1.3216,2.53519,1.07653],[1.31979,2.55597,1.07451],[1.3216,2.53519,1.07653],[1.32979,2.55597,1.0482],[1.53896,1.65421,1.39434],[1.54254,1.65902,1.39013],[1.53896,1.87858,1.39434],[1.54254,1.87378,1.39013],[1.53896,1.87858,1.39434],[1.54254,1.65902,1.39013],[1.53895,1.63514,1.21827],[1.54249,1.64056,1.22148],[1.53895,1.63514,1.37373],[1.54249,1.64056,1.37052],[1.53895,1.63514,1.37373],[1.54249,1.64056,1.22148],[1.5226,1.65232,1.396],[1.5226,1.63866,1.38985],[1.53226,1.65232,1.396],[1.53226,1.63866,1.38985],[1.53226,1.65232,1.396],[1.5226,1.63866,1.38985],[1.55986,1.57504,1.285],[1.58746,1.57504,1.285],[1.55986,1.58124,1.30126],[1.58746,1.58124,1.30126],[1.55986,1.58124,1.30126],[1.58746,1.57504,1.285],[1.58746,1.60503,1.09],[1.55986,1.60503,1.09],[1.58746,1.60503,1.2524],[1.55986,1.60503,1.2524],[1.58746,1.60503,1.2524],[1.55986,1.60503,1.09]],"colors":[[1,0,0,1]],"centers":[[1.550753,1.690407,-1.23509],[1.550777,1.7664,-1.233247],[1.550753,1.842393,-1.23509],[1.550737,1.848547,-1.241627],[1.550743,1.850293,-1.296],[1.550737,1.848547,-1.350374],[1.550753,1.842393,-1.35691],[1.550777,1.7664,-1.358753],[1.550753,1.690407,-1.35691],[1.550737,1.684253,-1.350374],[1.550743,1.682507,-1.296],[1.550737,1.684253,-1.241627],[1.550753,1.690407,1.23509],[1.550777,1.7664,1.233247],[1.550753,1.842393,1.23509],[1.550737,1.848547,1.241627],[1.550743,1.850293,1.296],[1.550737,1.848547,1.350373],[1.550753,1.842393,1.35691],[1.550777,1.7664,1.358753],[1.550753,1.690407,1.35691],[1.550737,1.684253,1.350373],[1.550743,1.682507,1.296],[1.550737,1.684253,1.241627],[-1.66094,0.7891887,1.077333],[-1.67414,0.7712804,1.064667],[-1.01932,0.732248,-1.087333],[-1.01504,0.7111241,-1.074667],[1.947333,0.5124201,-1.087333],[1.953107,0.46321,-1.074667],[-0.9037034,0.506,1.077333],[-0.9572317,0.598,1.064667],[2.29264,0.4324,-0.3573334],[2.29264,0.4508,0.3473334],[-0.2851317,0.4508,1.052],[0.2775697,0.4324,1.052],[-2.146143,0.4508,1.052],[-2.00217,0.4324,1.052],[-1.010505,0.729524,-1.062],[-1.010422,0.769048,-1.062],[1.97967,0.56762,-0.8916467],[1.98579,0.51841,-0.7212933],[1.10839,0.9366788,0.8816468],[1.07173,0.9069394,0.7112934],[-1.705623,0.787448,0.8816468],[-1.710447,0.766324,0.7112934],[-1.086497,0.8862283,-0.8916467],[-1.067523,0.8742627,-0.7212933],[-1.524416,0.4692,-0.18698],[-1.177282,0.4692,0.17698],[-2.19189,0.3686387,-0.3684667],[-2.19746,0.3660193,0.37449],[-2.423167,0.418975,-0.3726567],[-2.423043,0.420781,0.3727934],[-2.190634,0.418975,0.3726567],[-2.190757,0.420781,-0.3727934],[-2.185073,0.4913897,-0.3732],[-2.187257,0.4836193,0.3756633],[2.19779,0.43286,1.129097],[2.19356,0.46138,1.127193],[2.18583,0.43286,-1.116097],[2.18758,0.46138,-1.120693],[2.36828,0.4318617,-1.120693],[2.39738,0.4628344,-1.116097],[2.18758,0.6726084,-1.116097],[2.189917,0.676253,-1.120057],[2.19779,0.4001117,1.127193],[2.194147,0.396467,1.124653],[2.43032,0.504065,1.116097],[2.42567,0.499164,1.120693],[2.27746,0.3864,0.3705],[2.23974,0.3864,-0.3705],[2.311347,0.8448667,0.505],[2.295094,0.7927334,0.505],[2.323613,0.7682,0.8283334],[2.319627,0.7682,0.7366667],[1.3961,0.4646,0.711167],[1.4651,0.4646,0.559834],[1.349527,0.2898,0.711167],[1.338313,0.3128,0.559834],[1.453887,0.3089153,-0.2648043],[1.407313,0.3113617,-0.2761527],[1.4651,0.5067154,0.2648043],[1.3961,0.5091617,0.2761527],[1.3961,0.477823,-0.3795157],[1.4651,0.491046,-0.3505304],[1.3961,0.477823,0.195486],[1.4651,0.491046,0.224471],[-1.4559,0.3787333,0.05550067],[-1.4559,0.4216667,-0.05549967],[-1.329113,0.2668,0.05550067],[-1.375687,0.2668,-0.05549967],[-1.4559,0.4208487,-0.310197],[-1.4559,0.4662283,-0.298849],[-1.2489,0.4208487,0.2648043],[-1.2489,0.4662283,0.2761527],[-1.260113,0.339246,-0.35053],[-1.271327,0.303023,-0.379515],[-1.329113,0.280023,-0.379515],[-1.375687,0.293246,-0.35053],[-2.3184,0.5888,0.1048667],[-2.3184,0.6532,-0.1048667],[-2.33266,0.5888,-0.1048667],[-2.33266,0.6532,0.1048667],[0.1909743,2.39966,1.078333],[0.1832287,2.39812,1.066667],[1.0028,1.2604,-1.07],[1.0028,1.2512,-1.085],[0.957571,2.39812,-1.085],[0.9498254,2.39966,-1.07],[1.540153,1.80219,-1.199063],[1.541347,1.729003,-1.200467],[1.54012,1.89226,-1.38223],[1.541303,1.8923,-1.376283],[1.52582,1.885033,-1.19805],[1.52904,1.889587,-1.2001],[1.55986,1.5971,-1.268687],[1.55986,1.58917,-1.284973],[1.569427,1.79032,-1.286267],[1.564643,1.79032,-1.297133],[1.832663,1.506747,-1.075687],[2.067727,1.394243,-1.064843],[2.323127,1.2576,0.3446667],[2.325363,1.2498,-0.3546667],[-2.3184,1.386133,0.8985001],[-2.3184,1.502667,0.753],[-2.286053,1.92453,1.08551],[-2.296027,2.22986,1.08102],[-2.323153,0.7176,0.4122334],[-2.327907,0.7176,0.5098667],[-2.294603,2.576744,1.050897],[-2.287827,2.580557,1.05257],[-2.297037,2.576744,-0.35401],[-2.291363,2.580557,0.34524],[-2.28178,2.58437,-0.3538567],[-2.27787,2.58437,0.34401],[1.513867,1.92453,1.054843],[1.417733,2.22986,1.065687],[1.534487,1.895353,1.207597],[1.536717,1.896527,1.21297],[1.534487,1.888913,1.391403],[1.53672,1.883727,1.3929],[1.52904,1.633,1.269667],[1.52582,1.633,1.322333],[1.55986,1.595033,1.144133],[1.55986,1.585037,1.209133],[-1.541857,0.851398,1.077333],[-1.565163,0.847196,1.064667],[-1.096663,0.8310283,-1.087333],[-1.078417,0.8190627,-1.074667],[1.831717,0.777956,-1.087333],[1.859983,0.738578,-1.074667],[0.9901187,0.686,1.077333],[0.9240373,0.5500001,1.064667],[1.964117,0.48161,-1.062],[1.969007,0.54922,-1.062],[1.083743,0.8701394,1.052],[1.11336,0.9182787,1.052],[-1.701287,0.729524,1.052],[-1.701137,0.769048,1.052],[-1.089643,0.8678284,-1.062],[-1.074543,0.8374627,-1.062],[1.85712,0.833156,-0.8916467],[1.88708,0.7937781,-0.7212933],[1.37841,1.03312,0.8816468],[1.32961,1.02664,0.7112934],[-1.618517,0.8862283,0.8816468],[-1.639063,0.8742627,0.7112934],[-1.407713,0.9108,-0.8916467],[-1.290167,0.9108,-0.7212933],[0.273236,0.4692,0.17698],[-0.278456,0.4692,-0.18698],[-2.272277,0.5191613,-1.117447],[-2.347093,0.5217807,-1.111423],[-2.341523,0.3634,0.3684667],[-2.272277,0.3634,-0.3684667],[-2.189797,0.4741097,-0.3727933],[-2.190154,0.4723704,0.3732],[-2.428727,0.4913897,0.3732],[-2.426543,0.4836193,-0.3756633],[2.19779,0.5733134,1.129097],[2.19356,0.6208467,1.127193],[2.42429,0.6845683,-0.3705],[2.42852,0.6828167,0.3705],[2.27079,0.4973343,-1.131],[2.34246,0.5036897,-1.131],[2.19779,0.4001117,-1.127193],[2.194147,0.396467,-1.124653],[2.338437,0.3970927,-1.116097],[2.333707,0.3943583,-1.120203],[2.42227,0.5145707,1.129097],[2.42751,0.5201754,1.127193],[2.19356,0.5138201,1.127193],[2.19779,0.5198001,1.129097],[2.311347,0.7927334,-0.505],[2.295094,0.8448667,-0.505],[2.323613,0.9092667,-0.92],[2.319627,1.050333,-0.92],[1.3961,0.4646,-0.559834],[1.4651,0.4646,-0.711167],[1.349527,0.2898,-0.559834],[1.338313,0.3128,-0.711167],[1.349527,0.3319153,-0.2648043],[1.338313,0.3573617,-0.2761527],[1.3271,0.4208487,0.2648043],[1.3271,0.4662283,0.2761527],[1.5341,0.4481127,0.224471],[1.5341,0.3919563,0.195486],[1.3271,0.3919563,-0.3795157],[1.3271,0.4481127,-0.3505304],[-1.3869,0.4646,0.711167],[-1.3179,0.4646,0.559834],[-1.433473,0.2898,0.711167],[-1.444687,0.3128,0.559834],[-1.329113,0.3089153,-0.2648043],[-1.375687,0.3113617,-0.2761527],[-1.3179,0.5067154,0.2648043],[-1.3869,0.5091617,0.2761527],[-1.3869,0.477823,-0.379515],[-1.3179,0.491046,-0.35053],[-1.3869,0.477823,0.195486],[-1.3179,0.491046,0.224471],[-2.34692,0.5888,-0.5098667],[-2.34692,0.6532,-0.4122334],[-2.323153,0.7176,-0.1048667],[-2.327907,0.7176,0.1048667],[0.957571,2.39812,1.078333],[0.9498254,2.39966,1.066667],[1.0028,1.511867,-1.07],[1.0028,1.404533,-1.085],[0.1832287,2.39812,-1.07],[0.1909743,2.39966,-1.085],[1.541323,1.65301,-1.38969],[1.540123,1.646897,-1.3891],[1.54012,1.64054,-1.20977],[1.541303,1.6405,-1.215717],[1.52582,1.897913,-1.37995],[1.52904,1.896027,-1.3849],[1.58746,1.611207,-1.28402],[1.58746,1.596223,-1.29264],[1.58746,1.726777,-1.286267],[1.58746,1.666487,-1.297133],[1.422023,1.282997,1.09],[1.84728,1.272443,1.09],[1.191093,1.511867,-1.1],[1.610317,1.39953,-1.1],[-2.310133,1.502667,-1.075687],[-2.314267,1.386133,-1.064843],[-2.295094,0.414,-0.3546667],[-2.306747,0.414,0.3473334],[-2.323153,0.6532,0.6075],[-2.327907,0.5888,0.6075],[-1.08626,2.576744,1.066717],[0.1029901,2.580557,1.060533],[-2.310363,2.561623,-0.35524],[-2.306537,2.567277,0.3460667],[-2.30703,2.561623,-1.066143],[-2.30077,2.567277,-1.071723],[1.495627,1.92453,1.08102],[1.393653,2.22986,1.08551],[1.534487,1.888914,1.200597],[1.53672,1.883727,1.1991],[1.534487,1.643887,1.200597],[1.53672,1.64907,1.1991],[1.52904,1.8998,1.322333],[1.52582,1.8998,1.269667],[1.56906,1.57504,1.155],[1.57826,1.57504,1.22],[-1.608227,0.8310283,1.077333],[-1.627983,0.8190627,1.064667],[-1.04798,0.7891887,-1.087333],[-1.03579,0.7712804,-1.074667],[1.90602,0.6533433,-1.087333],[1.92379,0.6074867,-1.074667],[0.288579,0.414,1.077333],[-0.280798,0.414,1.064667],[2.192397,0.4508,-1.062],[2.081143,0.4324,-1.062],[0.913028,0.5684,1.052],[0.972066,0.7228,1.052],[-1.808107,0.5244,1.052],[-1.76269,0.6348,1.052],[-1.039743,0.8259887,-1.062],[-1.031427,0.7896804,-1.062],[1.935877,0.7085433,-0.8916467],[1.954713,0.6626867,-0.7212933],[1.235557,1.002246,0.8816468],[1.190303,0.984332,0.7112934],[-1.67334,0.8443887,0.8816468],[-1.68707,0.8264803,0.7112934],[-1.150237,0.906598,-0.8916467],[-1.127853,0.9023961,-0.7212933],[-2.11793,0.4692,-0.18698],[-1.99474,0.4692,0.17698],[-2.266707,0.504954,-1.125817],[-2.348237,0.510748,-1.124643],[-2.272277,0.5244,0.3684667],[-2.341523,0.5244,-0.3684667],[-2.190634,0.468825,-0.3726567],[-2.190757,0.467019,0.3727934],[-2.428727,0.39641,0.3732001],[-2.426544,0.40418,-0.3756634],[2.18758,0.5733134,1.120693],[2.18583,0.6208467,1.116097],[2.4345,0.6768367,-0.3705],[2.43625,0.6726084,0.3705],[2.27369,0.5141954,-1.131],[2.34637,0.5201754,-1.131],[2.19779,0.6726084,-1.127193],[2.194147,0.676253,-1.124653],[2.321553,0.3893667,-1.116097],[2.326714,0.391158,-1.120203],[2.41546,0.499164,1.129097],[2.42259,0.504065,1.127193],[2.18758,0.5138201,1.120693],[2.18583,0.5198001,1.116097],[2.311347,0.897,-0.1683333],[2.295094,0.897,0.1683333],[2.323613,1.1914,-0.8283334],[2.319627,1.1914,-0.7366667],[1.3961,0.4646,0.05550034],[1.4651,0.4646,-0.05550034],[1.349527,0.2898,0.05550034],[1.338313,0.3128,-0.05550034],[1.338313,0.3573617,-0.298849],[1.349527,0.3319153,-0.310197],[1.3271,0.4662283,0.298849],[1.3271,0.4208487,0.310197],[1.5341,0.3919563,0.3795157],[1.5341,0.4481127,0.3505304],[1.3271,0.3919563,-0.195486],[1.3271,0.4481127,-0.224471],[-1.4559,0.3787333,-0.5598333],[-1.4559,0.4216667,-0.7111667],[-1.329113,0.2668,-0.5598333],[-1.375687,0.2668,-0.7111667],[-1.329113,0.3113617,-0.298849],[-1.375687,0.3089153,-0.310197],[-1.3179,0.5091617,0.298849],[-1.3869,0.5067154,0.310197],[-1.3869,0.491046,-0.2244707],[-1.3179,0.477823,-0.1954853],[-1.3869,0.491046,0.3505304],[-1.3179,0.477823,0.3795157],[-2.34692,0.5888,0.4122334],[-2.34692,0.6532,0.5098667],[-2.337413,0.5888,0.3146],[-2.342167,0.6532,0.3146],[0.6942934,2.4012,1.078333],[0.4465067,2.4012,1.066667],[1.0028,1.288,-1.07],[1.0028,1.2788,-1.085],[0.6942934,2.4012,-1.085],[0.4465067,2.4012,-1.07],[1.534487,1.643887,-1.391403],[1.53672,1.64907,-1.3929],[1.534487,1.637447,-1.207597],[1.536717,1.636273,-1.21297],[1.52582,1.889587,-1.3919],[1.52904,1.885033,-1.39395],[1.58746,1.5971,-1.268687],[1.58746,1.58917,-1.284973],[1.58175,1.78461,-1.2754],[1.583043,1.731193,-1.2754],[1.573843,1.73561,-1.2754],[1.56906,1.6809,-1.2754],[1.191093,1.511867,1.09],[1.610317,1.39953,1.09],[2.3276,1.111667,-0.5516667],[2.3276,1.111667,0.5516667],[2.3276,0.9844,-0.5983334],[2.3276,0.9844,0.5983334],[2.3276,1.0672,0.9613334],[2.3276,0.8862666,0.5516667],[2.3276,0.8019333,0.5516667],[2.3276,0.8080667,1.002667],[2.3276,0.6501334,0.8696668],[2.3276,0.5320667,0.2116667],[2.3276,0.9614,0.1683333],[2.3276,0.9292001,-0.1683333],[2.3276,0.8862667,-0.5516667],[2.3276,0.8019333,-0.5516667],[2.3276,1.208267,0.8696667],[2.3276,0.759,-0.6900001],[2.3276,1.1638,0.215],[2.3276,1.1776,-0.1683333],[2.3276,1.225133,0.2116667],[2.3276,1.208267,-0.3513333],[2.3276,1.208267,-0.873],[2.3276,1.0672,-0.9646667],[2.3276,0.8080667,-1.009333],[2.3276,0.6409333,-0.8263334],[2.3276,0.6317334,-0.3513333],[2.3276,0.6409334,0.03199999],[-0.6666933,1.502667,1.09],[-1.471387,1.386133,1.09],[-2.310133,1.92453,1.065687],[-2.314267,2.22986,1.054843],[-2.323153,0.7176,-0.5098667],[-2.327907,0.7176,-0.4122334],[-2.28252,2.576744,1.064063],[-2.284057,2.580557,1.056677],[-2.312863,2.542117,1.056243],[-2.308126,2.549043,1.066413],[1.312637,2.576744,0.34401],[1.306963,2.580557,-0.35524],[1.302977,2.542117,1.08398],[1.312347,2.549043,1.078817],[1.54012,1.89226,1.20977],[1.541303,1.8923,1.215717],[1.541323,1.879787,1.38969],[1.540123,1.8859,1.3891],[1.52904,1.728373,1.396],[1.52582,1.804427,1.396],[1.58746,1.595033,1.144133],[1.58746,1.585037,1.209133],[-1.565163,0.847196,-1.087333],[-1.541857,0.851398,-1.074667],[1.472177,0.9779214,-1.087333],[1.518213,0.9714426,-1.074667],[1.380107,0.9779214,1.077333],[1.334073,0.9714426,1.064667],[1.774517,0.8659496,-1.062],[1.74337,0.9145653,-1.062],[1.47565,1.01472,1.052],[1.52133,0.9898413,1.052],[-1.54726,0.888198,1.052],[-1.56833,0.865596,1.052],[-1.54726,0.8881981,-1.062],[-1.56833,0.865596,-1.062],[1.476007,1.03312,-0.8916467],[1.524804,1.02664,-0.7212933],[1.78811,0.9027497,0.8816468],[1.74906,0.9329654,0.7112934],[-1.127853,0.9023961,0.8816468],[-1.150237,0.906598,0.7112934],[-1.68707,0.8264803,-0.8916467],[-1.67334,0.8443887,-0.7212933],[1.190303,0.984332,-0.8916467],[1.235557,1.002246,-0.7212933],[-2.42191,0.5191613,-0.3684667],[-2.41634,0.5217807,0.37449],[-2.184033,0.504954,0.37449],[-2.185177,0.510748,-0.3756633],[-2.345697,0.4136897,1.119193],[-2.267747,0.4154294,1.118787],[-2.345573,0.4367957,1.11797],[-2.268227,0.4510044,1.11797],[2.27369,0.5141954,1.131],[2.34637,0.5201754,1.131],[2.18758,0.3958833,-0.3705],[2.18583,0.4001117,0.3705],[2.19356,0.5733134,-1.127193],[2.19779,0.6208467,-1.129097],[2.326054,0.3887707,-0.3705],[2.320617,0.3875853,0.3705],[2.43524,0.5201754,-1.116097],[2.43248,0.5145707,-1.120693],[2.42429,0.6726084,1.127193],[2.427933,0.676253,1.124653],[2.18583,0.50186,-1.116097],[2.18758,0.49588,-1.120693],[2.24323,0.4323177,1.131],[2.312,0.462292,1.131],[2.27884,1.097867,0.1683333],[2.27884,1.045733,-0.1683333],[1.4651,0.3787333,0.8625001],[1.3961,0.4216667,0.8625001],[1.522887,0.3128,-0.559834],[1.511673,0.2898,-0.711167],[1.522887,0.3573617,0.2761527],[1.511673,0.3319153,0.2648043],[1.5341,0.4662283,-0.2761527],[1.5341,0.4208487,-0.2648043],[1.5341,0.3919563,-0.195486],[1.5341,0.4481127,-0.224471],[1.522887,0.339246,0.3505304],[1.511673,0.303023,0.3795157],[-1.375687,0.3128,-0.8625001],[-1.3179,0.2898,-0.8625001],[-1.260113,0.3128,0.711167],[-1.271327,0.2898,0.559834],[-1.375687,0.3113617,0.298849],[-1.329113,0.3089153,0.310197],[-1.3179,0.5091617,-0.2761527],[-1.3869,0.5067154,-0.2648043],[-1.329113,0.293246,0.3505304],[-1.375687,0.280023,0.3795157],[-1.433473,0.303023,0.3795156],[-1.444687,0.339246,0.3505304],[-2.2954,2.146667,0.2025],[-2.2954,1.882933,-0.2025],[-2.323153,0.5244,0.5098667],[-2.327907,0.5244,0.4122334],[1.001259,2.348227,1.078333],[0.9997187,2.355973,1.066667],[0.3557333,0.4784,1.066667],[0.5734667,0.4784,1.078333],[0.1395407,2.348227,-1.07],[0.1410813,2.355973,-1.085],[1.534493,1.729003,-1.395447],[1.536727,1.804423,-1.394893],[1.53449,1.633713,-1.27009],[1.53672,1.634427,-1.322333],[1.52582,1.636773,-1.3849],[1.52904,1.634887,-1.37995],[1.56906,1.586227,-1.303507],[1.57826,1.591213,-1.305753],[1.586167,1.780193,-1.297133],[1.584873,1.783317,-1.286267],[2.30299,1.249813,1.08102],[2.29078,1.25766,1.08551],[2.276477,1.276523,-1.09551],[2.292483,1.27094,-1.09102],[-2.286053,1.386133,1.08551],[-2.296027,1.502667,1.08102],[2.304293,0.414,0.3446667],[2.315947,0.414,-0.3573334],[-2.3184,0.6532,-0.5098667],[-2.3184,0.5888,-0.4122334],[-1.08739,2.561623,1.08124],[0.10375,2.567277,1.07707],[-2.315593,2.549043,0.3446667],[-2.316997,2.542117,-0.3560667],[-2.287377,2.542117,-1.09398],[-2.296747,2.549043,-1.088817],[1.29361,2.58437,1.045983],[1.296417,2.58437,1.04904],[1.393653,2.22986,-1.09551],[1.495627,1.92453,-1.09102],[1.540153,1.80219,1.199063],[1.541347,1.729003,1.200467],[1.54012,1.89226,1.38223],[1.541303,1.8923,1.376283],[1.52582,1.885033,1.19805],[1.52904,1.889587,1.2001],[1.55986,1.5971,1.268687],[1.55986,1.58917,1.284973],[1.569427,1.79032,1.286267],[1.564643,1.79032,1.297133],[-1.405523,0.8556001,-1.087333],[-1.292497,0.8556001,-1.074667],[1.607147,0.9459697,-1.087333],[1.650043,0.9269754,-1.074667],[1.245347,0.9470487,1.077333],[1.202653,0.9291334,1.064667],[1.86958,0.756978,-1.062],[1.849217,0.814756,-1.062],[1.378053,1.01472,1.052],[1.331663,0.9898413,1.052],[-1.61535,0.8678284,1.052],[-1.63194,0.8374627,1.052],[-1.29245,0.8924,-1.062],[-1.40776,0.874,-1.062],[1.61907,1.001167,-0.8916467],[1.66454,0.982174,-0.7212933],[1.66454,0.982174,0.8816468],[1.61907,1.001167,0.7112934],[-1.290167,0.9108,0.8816468],[-1.407713,0.9108,0.7112934],[-1.639063,0.8742627,-0.8916467],[-1.618517,0.8862283,-0.7212933],[1.32961,1.02664,-0.8916467],[1.37841,1.03312,-0.7212933],[-2.42191,0.3686387,0.3684667],[-2.41634,0.3660193,-0.37449],[-2.272277,0.3686387,1.117447],[-2.347093,0.3660193,1.111423],[-2.268103,0.4741096,1.119193],[-2.346054,0.4723704,1.118787],[-2.346054,0.4913897,-1.124527],[-2.265563,0.4836193,-1.122063],[2.43625,0.6208467,-1.116097],[2.4345,0.5733134,-1.120693],[2.43699,0.5201754,-0.3705],[2.43598,0.5145707,0.3705],[2.43625,0.5733134,1.116097],[2.4345,0.6208467,1.120693],[2.19779,0.3899034,-1.116097],[2.194147,0.3922387,-1.120057],[2.42567,0.499164,-1.120693],[2.43032,0.504065,-1.116097],[2.4345,0.6726084,1.116097],[2.432163,0.676253,1.120057],[2.18583,0.50186,1.116097],[2.18758,0.49588,1.120693],[2.34738,0.6208467,-1.131],[2.2747,0.5733134,-1.131],[2.323613,1.050333,-0.645],[2.319627,0.9092667,-0.645],[1.453887,0.3128,0.8625001],[1.3961,0.2898,0.8625001],[1.522887,0.3128,0.05550034],[1.511673,0.2898,-0.05550034],[1.407313,0.3089153,0.2648043],[1.453887,0.3113617,0.2761527],[1.4651,0.5067154,-0.310197],[1.3961,0.5091617,-0.298849],[1.453887,0.280023,0.195486],[1.407313,0.293246,0.224471],[1.349527,0.303023,0.195486],[1.338313,0.339246,0.224471],[-1.3869,0.4646,-0.5598333],[-1.3179,0.4646,-0.7111667],[-1.433473,0.2898,-0.5598333],[-1.444687,0.3128,-0.7111667],[-1.433473,0.3319153,-0.2648043],[-1.444687,0.3573617,-0.2761527],[-1.4559,0.4208487,0.2648043],[-1.4559,0.4662283,0.2761527],[-1.2489,0.4481127,0.224471],[-1.2489,0.3919563,0.195486],[-1.4559,0.3919563,-0.379515],[-1.4559,0.4481127,-0.35053],[-2.2954,1.502667,0.2025],[-2.2954,1.386133,-0.2025],[-2.337413,0.5244,0.5098667],[-2.342167,0.5244,0.4122334],[0.9937906,2.370287,1.078333],[0.9894029,2.376853,1.066667],[0.8617333,0.6164,1.066667],[0.9322667,0.7544,1.078333],[0.1470093,2.370287,-1.07],[0.1513967,2.376853,-1.085],[1.54012,1.64054,-1.38223],[1.541303,1.6405,-1.376283],[1.54013,1.895853,-1.32084],[1.54131,1.894047,-1.27009],[1.52904,1.647767,-1.19805],[1.52582,1.643213,-1.2001],[1.56906,1.621203,-1.273153],[1.57826,1.616217,-1.270907],[1.580457,1.787733,-1.297133],[1.577333,1.789027,-1.286267],[2.31486,1.257613,1.065687],[2.32123,1.2498,1.054843],[1.428057,1.259047,-1.1],[1.85555,1.249847,-1.1],[-2.310133,1.386133,1.065687],[-2.314267,1.502667,1.054843],[-1.99304,0.414,1.077333],[-2.137013,0.414,1.064667],[-2.3184,0.6532,0.4122334],[-2.3184,0.5888,0.5098667],[0.10299,2.561623,-1.09124],[-1.08815,2.567277,-1.08707],[-2.287377,2.542117,1.08398],[-2.296747,2.549043,1.078817],[-2.28252,2.576744,-1.074063],[-2.284057,2.580557,-1.066677],[1.015059,2.350833,1.09],[-1.471387,1.92453,1.09],[1.043426,2.368007,1.09],[-0.6666933,2.164957,1.09],[1.023698,2.379714,1.09],[1.012744,2.390667,1.09],[1.0304,2.121827,1.09],[1.0166,1.878793,1.09],[1.044967,1.642967,1.09],[1.21946,1.6238,1.09],[1.373533,1.6284,1.09],[1.52097,1.648033,1.09],[0.9984323,2.396594,1.09],[1.09296,2.395067,1.09],[1.46319,1.948763,1.09],[1.376177,2.19257,1.09],[1.282407,2.425937,1.09],[1.160887,2.442797,1.09],[0.8108267,2.445863,1.09],[-0.2618933,2.490527,1.09],[-0.6651527,2.41313,1.09],[-0.6592246,2.427444,1.09],[-0.648271,2.438397,1.09],[-0.633959,2.444324,1.09],[1.302977,2.542117,-1.09398],[1.312347,2.549043,-1.088817],[1.534493,1.803793,1.196553],[1.536727,1.72837,1.197107],[1.534487,1.895353,1.384403],[1.536717,1.896527,1.37903],[1.52582,1.896027,1.2071],[1.52904,1.897913,1.21205],[1.55986,1.596223,1.29264],[1.55986,1.611207,1.28402],[-1.627983,0.8190627,-1.087333],[-1.608227,0.8310283,-1.074667],[1.334073,0.9714426,-1.087333],[1.380107,0.9779214,-1.074667],[1.518213,0.9714426,1.077333],[1.472177,0.9779214,1.064667],[1.655733,0.9453753,-1.062],[1.615953,0.9827683,-1.062],[1.655733,0.9453753,1.052],[1.615953,0.9827683,1.052],[-1.405477,0.8924,1.052],[-1.290213,0.874,1.052],[-1.63194,0.8374627,-1.062],[-1.61535,0.8678284,-1.062],[1.378053,1.01472,-1.062],[1.331663,0.9898413,-1.062],[1.88708,0.7937781,0.8816468],[1.85712,0.833156,0.7112934],[-1.067523,0.8742627,0.8816468],[-1.086497,0.8862283,0.7112934],[-1.710447,0.766324,-0.8916467],[-1.705623,0.787448,-0.7212933],[1.07173,0.9069394,-0.8916467],[1.10839,0.9366788,-0.7212933],[-2.429767,0.504954,-0.37449],[-2.428623,0.510748,0.3756633],[-2.19189,0.5191613,0.3684667],[-2.19746,0.5217807,-0.37449],[-2.345573,0.418975,1.118243],[-2.268103,0.420781,1.118107],[-2.19088,0.4367957,0.3726567],[-2.19088,0.4510044,-0.3726567],[2.27079,0.4973343,1.131],[2.34246,0.5036897,1.131],[2.19779,0.3881517,-0.3705],[2.19356,0.3899034,0.3705],[2.18583,0.5733134,-1.116097],[2.18758,0.6208467,-1.120693],[2.340397,0.3966367,-0.3705],[2.335943,0.3932963,0.3705],[2.42751,0.5201754,-1.127193],[2.42227,0.5145707,-1.129097],[2.42429,0.6828167,1.116097],[2.427933,0.6804814,1.120057],[2.19779,0.50186,-1.129097],[2.19356,0.49588,-1.127193],[2.2747,0.5733134,1.131],[2.34738,0.6208467,1.131],[2.27884,0.8448667,0.1683333],[2.27884,0.7927334,-0.1683333],[1.4651,0.4216667,-0.8625001],[1.3961,0.3787333,-0.8625001],[1.5341,0.4216667,0.711167],[1.5341,0.3787333,0.559834],[1.511673,0.3319153,0.310197],[1.522887,0.3573617,0.298849],[1.5341,0.4208487,-0.310197],[1.5341,0.4662283,-0.298849],[1.5341,0.4481127,-0.3505304],[1.5341,0.3919563,-0.3795156],[1.522887,0.339246,0.224471],[1.511673,0.303023,0.195486],[-1.329113,0.3128,0.8625001],[-1.3869,0.2898,0.8625001],[-1.260113,0.3128,0.05550067],[-1.271327,0.2898,-0.05549967],[-1.375687,0.3089153,0.2648043],[-1.329113,0.3113617,0.2761527],[-1.3179,0.5067154,-0.310197],[-1.3869,0.5091617,-0.298849],[-1.329113,0.280023,0.195486],[-1.375687,0.293246,0.224471],[-1.433473,0.303023,0.195486],[-1.444687,0.339246,0.224471],[-2.310733,1.502667,-0.6075],[-2.303067,1.386133,-0.6075],[-2.337413,0.5244,-0.4122334],[-2.342167,0.5244,-0.5098667],[1.0028,1.859627,1.066667],[1.0028,2.100054,1.078333],[0.138,1.005867,1.066667],[0.138,0.7421333,1.078333],[0.138,1.859627,-1.07],[0.138,2.100054,-1.085],[1.540153,1.730603,-1.392937],[1.541347,1.803793,-1.391533],[1.54013,1.636947,-1.27116],[1.54131,1.638753,-1.32191],[1.52582,1.647767,-1.39395],[1.52904,1.643213,-1.3919],[1.56906,1.577107,-1.29042],[1.57826,1.579173,-1.29584],[1.57826,1.60503,-1.144133],[1.56906,1.60503,-1.198267],[1.832663,1.506747,1.065687],[2.067727,1.394243,1.054843],[2.310723,1.270843,-1.064843],[2.300557,1.276323,-1.075687],[-2.310133,0.6992,1.065687],[-2.314267,0.9844,1.054843],[2.06768,0.414,-1.087333],[2.178934,0.414,-1.074667],[-2.310733,0.9016001,0.6075],[-2.303067,1.0856,0.6075],[-1.08612,2.58437,0.3438566],[0.1017199,2.58437,-0.3538567],[1.310203,2.576744,1.050897],[1.303427,2.580557,1.05257],[1.308817,2.561623,-1.084953],[1.30001,2.567277,-1.084417],[-2.27801,2.58437,-1.055983],[-2.280817,2.58437,-1.05904],[-1.086827,2.549043,-1.09694],[0.10375,2.542117,-1.09847],[1.534487,1.643887,1.391403],[1.53672,1.64907,1.3929],[1.534487,1.637447,1.207597],[1.536717,1.636273,1.21297],[1.52582,1.889587,1.3919],[1.52904,1.885033,1.39395],[1.58746,1.58917,1.284973],[1.58746,1.5971,1.268687],[1.58175,1.78461,1.2754],[1.573843,1.785903,1.2754],[1.57826,1.731193,1.2754],[1.56906,1.6809,1.2754],[-1.75134,0.598,1.077333],[-1.80143,0.506,1.064667],[-0.280798,0.414,-1.087333],[0.288579,0.414,-1.074667],[1.92379,0.6074867,1.077333],[1.90602,0.6533433,1.064667],[-1.03579,0.7712804,1.077333],[-1.04798,0.7891887,1.064667],[1.964117,0.48161,1.052],[1.969007,0.54922,1.052],[-1.010505,0.729524,1.052],[-1.010422,0.769048,1.052],[-2.146143,0.4508,-1.062],[-2.00217,0.4324,-1.062],[-0.2851317,0.4508,-1.062],[0.2775697,0.4324,-1.062],[2.278533,0.4692,-0.18698],[2.264427,0.4692,0.3473334],[0.273236,0.4692,0.8816468],[-0.278456,0.4692,0.7112934],[-2.132037,0.4692,0.8816468],[-1.99474,0.4692,0.7112934],[-1.006059,0.787448,-0.8916467],[-1.001607,0.766324,-0.7212933],[-1.472323,0.5612,0.54094],[-1.180858,0.6532,0.54094],[-1.691893,0.8053564,0.54094],[-1.457345,0.7842323,0.54094],[-1.018739,0.8053564,0.54094],[-1.235105,0.8232647,0.54094],[-1.10888,0.8904303,0.54094],[-1.59428,0.8904304,0.54094],[-1.293593,0.8784647,0.54094],[-1.452497,0.8946323,0.54094],[-2.341523,0.3686387,-1.117447],[-2.266707,0.3660193,-1.111423],[-2.345697,0.4741096,-1.119193],[-2.267747,0.4723704,-1.118787],[-2.424003,0.4741097,0.3727933],[-2.423647,0.4723704,-0.3732],[-2.267747,0.4913897,1.124527],[-2.348237,0.4836193,1.122063],[2.18408,0.49588,-0.3705],[2.18408,0.50186,0.3705],[2.23974,0.3881517,-1.116097],[2.278397,0.389933,-1.120693],[2.23974,0.3881517,1.116097],[2.278397,0.389933,1.120693],[2.42429,0.6726084,-1.127193],[2.427933,0.676253,-1.124653],[2.19779,0.3899034,1.116097],[2.194147,0.3922387,1.120057],[2.324767,0.3954103,1.126703],[2.33176,0.3986107,1.126703],[2.438,0.6208467,0.3705],[2.438,0.5733134,-0.3705],[2.311347,1.045733,-0.505],[2.295094,1.097867,-0.505],[2.323613,1.1914,0.7366667],[2.319627,1.1914,0.8283334],[1.3271,0.3787333,0.05550034],[1.3271,0.4216667,-0.05550034],[1.453887,0.2668,0.05550034],[1.407313,0.2668,-0.05550034],[1.3271,0.4208487,-0.310197],[1.3271,0.4662283,-0.298849],[1.5341,0.4208487,0.2648043],[1.5341,0.4662283,0.2761527],[1.522887,0.339246,-0.3505304],[1.511673,0.303023,-0.3795156],[1.453887,0.280023,-0.3795157],[1.407313,0.293246,-0.3505304],[-1.2489,0.4216667,-0.5598333],[-1.2489,0.3787333,-0.7111667],[-1.433473,0.3319153,0.310197],[-1.444687,0.3573617,0.298849],[-1.271327,0.3319153,-0.2648043],[-1.260113,0.3573617,-0.2761527],[-1.4559,0.3919563,0.195486],[-1.4559,0.4481127,0.224471],[-1.433473,0.303023,-0.379515],[-1.444687,0.339246,-0.35053],[-2.310733,1.386133,0.6075],[-2.303067,1.502667,0.6075],[-2.337413,0.6532,0.6075],[-2.342167,0.5888,0.6075],[1.0028,1.008933,1.066667],[1.0028,1.125467,1.078333],[0.1513967,2.376853,1.078333],[0.1470093,2.370287,1.066667],[0.9322667,0.7544,-1.085],[0.8617333,0.6164,-1.07],[0.989403,2.376853,-1.07],[0.9937906,2.370287,-1.085],[1.541323,1.879787,-1.20231],[1.540123,1.8859,-1.2029],[1.541323,1.65301,-1.20231],[1.540123,1.646897,-1.2029],[1.52904,1.804427,-1.196],[1.52582,1.728373,-1.196],[1.55986,1.73561,-1.286267],[1.55986,1.670903,-1.297133],[2.310723,1.270843,1.054843],[2.300557,1.276323,1.065687],[2.323467,0.966,-1.064843],[2.319334,0.69,-1.075687],[1.84093,1.506697,0.3446667],[2.07186,1.394193,-0.3546667],[-2.3184,1.386133,-0.7563334],[-2.3184,1.502667,-0.9051667],[-2.286053,1.92453,-1.09551],[-2.296027,2.22986,-1.09102],[-2.288507,0.414,1.072843],[-2.302613,0.414,1.05751],[1.320197,2.561623,1.06255],[1.319703,2.567277,1.052953],[1.310203,2.576744,-1.060897],[1.303427,2.580557,-1.06257],[-2.27801,2.58437,1.045983],[-2.280817,2.58437,1.04904],[1.331193,2.549043,-0.3546667],[1.332597,2.542117,0.3460667],[1.534493,1.729003,1.395447],[1.536727,1.804423,1.394893],[1.53449,1.633713,1.27009],[1.53672,1.634427,1.322333],[1.52582,1.636773,1.3849],[1.52904,1.634887,1.37995],[1.56906,1.586227,1.303507],[1.57826,1.591213,1.305753],[1.586167,1.780193,1.297133],[1.584873,1.783317,1.286267],[-1.292497,0.8556001,1.077333],[-1.405523,0.8556001,1.064667],[-1.15795,0.851398,-1.087333],[-1.13643,0.847196,-1.074667],[1.729777,0.8777654,-1.087333],[1.766613,0.8475497,-1.074667],[1.125373,0.8814787,1.077333],[1.090787,0.8517393,1.064667],[1.934453,0.6258867,-1.062],[1.92628,0.6901433,-1.062],[1.197683,0.9475333,1.052],[1.237967,0.9838473,1.052],[-1.669383,0.8259887,1.052],[-1.678627,0.7896804,1.052],[-1.15252,0.888198,-1.062],[-1.133283,0.865596,-1.062],[1.74906,0.9329654,-0.8916467],[1.78811,0.9027497,-0.7212933],[1.524804,1.02664,0.8816468],[1.476007,1.03312,0.7112934],[-1.549497,0.906598,0.8816468],[-1.573733,0.9023961,0.7112934],[-1.573733,0.9023961,-0.8916467],[-1.549497,0.906598,-0.7212933],[2.164183,0.4692,0.17698],[2.078047,0.4692,-0.18698],[-2.429767,0.382846,0.37449],[-2.428623,0.377052,-0.3756633],[-2.266707,0.382846,1.125817],[-2.348237,0.377052,1.124643],[-2.268227,0.468825,1.118243],[-2.345697,0.467019,1.118107],[-2.346054,0.39641,-1.124527],[-2.265563,0.40418,-1.122063],[2.42852,0.6208467,-1.127193],[2.42429,0.5733134,-1.129097],[2.43207,0.504065,-0.3705],[2.42917,0.499164,0.3705],[2.42852,0.5733134,1.127193],[2.42429,0.6208467,1.129097],[2.18758,0.4001117,-1.116097],[2.189917,0.396467,-1.120057],[2.33176,0.3986107,-1.126703],[2.324767,0.3954103,-1.126703],[2.43524,0.5201754,1.116097],[2.43248,0.5145707,1.120693],[2.19779,0.50186,1.129097],[2.19356,0.49588,1.127193],[2.24323,0.4323177,-1.131],[2.312,0.462292,-1.131],[2.323613,0.7682,-0.7366667],[2.319627,0.7682,-0.8283334],[1.407313,0.3128,-0.8625001],[1.4651,0.2898,-0.8625001],[1.522887,0.3128,0.711167],[1.511673,0.2898,0.559834],[1.407313,0.3113617,0.298849],[1.453887,0.3089153,0.310197],[1.4651,0.5091617,-0.2761527],[1.3961,0.5067154,-0.2648043],[1.453887,0.293246,0.3505304],[1.407313,0.280023,0.3795157],[1.349527,0.303023,0.3795156],[1.338313,0.339246,0.3505304],[-1.3869,0.4646,0.05550067],[-1.3179,0.4646,-0.05549967],[-1.433473,0.2898,0.05550067],[-1.444687,0.3128,-0.05549967],[-1.444687,0.3573617,-0.298849],[-1.433473,0.3319153,-0.310197],[-1.4559,0.4662283,0.298849],[-1.4559,0.4208487,0.310197],[-1.2489,0.3919563,0.3795157],[-1.2489,0.4481127,0.3505304],[-1.4559,0.3919563,-0.1954853],[-1.4559,0.4481127,-0.2244707],[-2.2954,1.0856,0.2025],[-2.2954,0.9016001,-0.2025],[-2.323153,0.5244,0.1048667],[-2.327907,0.5244,-0.1048667],[0.9784493,2.387807,1.078333],[0.9718827,2.392193,1.066667],[0.3557333,0.7421333,-1.055],[0.2238667,2.31104,-1.055],[0.644,0.8801333,-1.055],[0.7145333,1.134667,-1.055],[0.4262667,1.376933,-1.055],[0.1891343,2.340167,-1.055],[0.4327067,2.339867,-1.055],[0.6535067,2.370533,-1.055],[0.8893334,2.31104,-1.055],[0.9378654,2.340167,-1.055],[0.952177,2.33424,-1.055],[0.1748223,2.33424,-1.055],[0.963131,2.323287,-1.055],[0.9368593,2.089707,-1.055],[0.1638687,2.323287,-1.055],[0.1579407,2.308973,-1.055],[0.1564,2.0608,-1.055],[0.1748,1.841533,-1.055],[0.2036267,1.632693,-1.055],[0.4734933,1.4858,-1.055],[0.7007334,1.478133,-1.055],[0.9095733,1.50696,-1.055],[0.9690593,1.7628,-1.055],[1.001259,2.1078,-1.055],[1.001259,1.76004,-1.055],[1.001259,1.643507,-1.055],[1.001259,1.625107,-1.055],[0.1623503,2.387807,-1.07],[0.1689167,2.392193,-1.085],[1.534487,1.637447,-1.384403],[1.536717,1.636273,-1.37903],[1.53449,1.899087,-1.32191],[1.53672,1.898373,-1.269667],[1.52904,1.636773,-1.2071],[1.52582,1.634887,-1.21205],[1.56906,1.609163,-1.26324],[1.57826,1.607097,-1.25782],[1.58175,1.78461,-1.308],[1.573843,1.785903,-1.308],[1.57826,1.721197,-1.308],[1.56906,1.660907,-1.308],[1.428057,1.259047,1.09],[1.85555,1.249847,1.09],[1.422023,1.282997,-1.1],[1.84728,1.272443,-1.1],[-2.286053,1.502667,-1.09551],[-2.296027,1.386133,-1.09102],[-1.99304,0.414,-1.087333],[-2.137013,0.414,-1.074667],[-2.323153,0.5888,-0.6075],[-2.327907,0.6532,-0.6075],[0.10186,2.576744,-1.076717],[-1.08739,2.580557,-1.070534],[1.29812,2.576744,-1.074063],[1.299657,2.580557,-1.066677],[1.302937,2.561623,1.077607],[1.308057,2.567277,1.070783],[0.1024266,2.549043,1.08694],[-1.08815,2.542117,1.08847],[1.541323,1.879787,1.20231],[1.540123,1.8859,1.2029],[1.541323,1.65301,1.20231],[1.540123,1.646897,1.2029],[1.52904,1.804427,1.196],[1.52582,1.728373,1.196],[1.55986,1.73561,1.286267],[1.55986,1.670903,1.297133],[-1.67414,0.7712804,-1.087333],[-1.66094,0.7891887,-1.074667],[1.202653,0.9291334,-1.087333],[1.245347,0.9470487,-1.074667],[1.650043,0.9269754,1.077333],[1.607147,0.9459697,1.064667],[1.47565,1.01472,-1.062],[1.52133,0.9898413,-1.062],[1.774517,0.8659496,1.052],[1.74337,0.9145653,1.052],[-1.133283,0.865596,1.052],[-1.15252,0.8881981,1.052],[-1.678627,0.7896804,-1.062],[-1.669383,0.8259887,-1.062],[1.197683,0.9475333,-1.062],[1.237967,0.9838473,-1.062],[1.954713,0.6626867,0.8816468],[1.935877,0.7085433,0.7112934],[-1.02319,0.8264803,0.8816468],[-1.03587,0.8443887,0.7112934],[-1.819457,0.5612,-0.8916467],[-1.767363,0.6532,-0.7212933],[0.8949754,0.6052001,-0.8916467],[0.9650227,0.7412,-0.7212933],[-1.232951,0.5612,-0.55094],[-1.527992,0.6532,-0.55094],[-1.691893,0.8053564,-0.55094],[-1.457345,0.7842323,-0.55094],[-1.018739,0.8053564,-0.55094],[-1.235105,0.8232647,-0.55094],[-1.10888,0.8904304,-0.55094],[-1.59428,0.8904303,-0.55094],[-1.41114,0.8784647,-0.55094],[-1.24881,0.8946323,-0.55094],[-2.347093,0.504954,1.125817],[-2.265563,0.510748,1.124643],[-2.268103,0.4136897,-1.119193],[-2.346054,0.4154294,-1.118787],[-2.42292,0.4367957,-0.3726567],[-2.42292,0.4510044,0.3726567],[2.2747,0.6768367,-1.127193],[2.34738,0.6726084,-1.129097],[2.18758,0.6768367,0.3705],[2.18583,0.6726084,-0.3705],[2.35765,0.4327737,1.129097],[2.38944,0.4632903,1.127193],[2.34738,0.6768367,1.127193],[2.2747,0.6726084,1.129097],[2.41546,0.499164,-1.129097],[2.42259,0.504065,-1.127193],[2.19779,0.6828167,1.116097],[2.194147,0.6804814,1.120057],[2.19356,0.5138201,-1.127193],[2.19779,0.5198001,-1.129097],[2.18408,0.43286,-0.3705],[2.18408,0.46138,0.3705],[2.311347,0.9936001,0.1683333],[2.295094,0.9936001,-0.1683333],[2.31564,1.050333,-0.7366667],[2.31564,0.9092667,-0.8283334],[1.5341,0.4216667,0.05550034],[1.5341,0.3787333,-0.05550034],[1.338313,0.3573617,0.2761527],[1.349527,0.3319153,0.2648043],[1.522887,0.3573617,-0.298849],[1.511673,0.3319153,-0.310197],[1.3271,0.3919563,0.3795156],[1.3271,0.4481127,0.3505304],[1.349527,0.303023,-0.195486],[1.338313,0.339246,-0.224471],[-1.3179,0.3787333,0.8625001],[-1.3869,0.4216667,0.8625001],[-1.260113,0.3128,-0.5598333],[-1.271327,0.2898,-0.7111667],[-1.260113,0.3573617,0.2761527],[-1.271327,0.3319153,0.2648043],[-1.2489,0.4662283,-0.2761527],[-1.2489,0.4208487,-0.2648043],[-1.2489,0.3919563,-0.1954853],[-1.2489,0.4481127,-0.2244707],[-1.260113,0.339246,0.3505304],[-1.271327,0.303023,0.3795157],[-2.310733,2.146667,-0.6075],[-2.303067,1.882933,-0.6075],[-2.323153,0.5244,-0.4122334],[-2.327907,0.5244,-0.5098667],[1.0028,1.404533,1.066667],[1.0028,1.511867,1.078333],[0.138,1.502667,1.066667],[0.138,1.386133,1.078333],[0.138,1.386133,-1.07],[0.138,1.502667,-1.085],[1.534487,1.895353,-1.207597],[1.536717,1.896527,-1.21297],[1.534487,1.888914,-1.391403],[1.53672,1.883727,-1.3929],[1.52904,1.633,-1.269667],[1.52582,1.633,-1.322333],[1.55986,1.595033,-1.144133],[1.55986,1.585037,-1.209133],[1.808583,1.506863,1.08551],[2.049487,1.39441,1.08102],[2.29078,1.25766,-1.09551],[2.30299,1.249813,-1.09102],[-2.286053,0.6992,1.08551],[-2.296027,0.9844,1.08102],[2.06768,0.414,1.077333],[2.178933,0.414,1.064667],[-2.310733,0.7176,-0.2025],[-2.303067,0.7176,0.2025],[2.311813,0.414,1.05751],[2.297707,0.414,1.072843],[1.29812,2.576744,1.064063],[1.299657,2.580557,1.056677],[-2.293217,2.561623,1.074953],[-2.28441,2.567277,1.074417],[-1.08612,2.58437,1.05009],[0.10186,2.58437,1.04583],[1.015059,2.350833,-1.1],[-1.471387,1.92453,-1.1],[1.043426,2.368007,-1.1],[-0.6666933,2.164957,-1.1],[1.023698,2.379714,-1.1],[1.012744,2.390667,-1.1],[1.0304,2.121827,-1.1],[1.0166,1.878793,-1.1],[1.044967,1.642967,-1.1],[1.21946,1.6238,-1.1],[1.373533,1.6284,-1.1],[1.52097,1.648033,-1.1],[0.9984323,2.396594,-1.1],[1.09296,2.395067,-1.1],[1.46319,1.948764,-1.1],[1.376177,2.19257,-1.1],[1.282407,2.425937,-1.1],[1.160887,2.442797,-1.1],[0.8108267,2.445863,-1.1],[-0.2618933,2.490527,-1.1],[-0.6651527,2.41313,-1.1],[-0.6592246,2.427444,-1.1],[-0.648271,2.438397,-1.1],[-0.633959,2.444323,-1.1],[1.541323,1.65301,1.38969],[1.540123,1.646897,1.3891],[1.54012,1.64054,1.20977],[1.541303,1.6405,1.215717],[1.52582,1.897913,1.37995],[1.52904,1.896027,1.3849],[1.58746,1.611207,1.28402],[1.58746,1.596223,1.29264],[1.58746,1.726777,1.286267],[1.58746,1.666487,1.297133],[-1.696613,0.7111241,-1.087333],[-1.691977,0.732248,-1.074667],[1.090787,0.8517393,-1.087333],[1.125373,0.8814787,-1.074667],[1.766613,0.8475497,1.077333],[1.729777,0.8777654,1.064667],[-1.13643,0.847196,1.077333],[-1.15795,0.851398,1.064667],[1.86958,0.756978,1.052],[1.849217,0.814756,1.052],[-1.074543,0.8374627,1.052],[-1.089643,0.8678284,1.052],[-1.701287,0.729524,-1.062],[-1.701137,0.769048,-1.062],[1.083743,0.8701393,-1.062],[1.11336,0.9182787,-1.062],[1.98579,0.51841,0.8816468],[1.97967,0.56762,0.7112934],[-1.001607,0.766324,0.8816468],[-1.006059,0.787448,0.7112934],[-2.132037,0.4692,-0.8916467],[-1.99474,0.4692,-0.7212933],[-0.278456,0.4692,-0.8916467],[0.273236,0.4692,-0.7212933],[1.283969,0.6052001,-0.55094],[1.666843,0.65441,-0.55094],[1.384557,0.8201494,-0.55094],[1.678547,0.779216,-0.55094],[1.447633,0.9136593,-0.55094],[1.675003,0.8823647,-0.55094],[1.511727,0.9774314,-0.55094],[1.654793,0.9584384,-0.55094],[1.570273,1.007647,-0.55094],[-2.341523,0.5191613,1.117447],[-2.266707,0.5217807,1.111423],[-2.268227,0.418975,-1.118243],[-2.345697,0.420781,-1.118107],[-2.268227,0.4367957,-1.11797],[-2.345573,0.4510044,-1.11797],[2.2747,0.6845683,-1.116097],[2.34738,0.6828167,-1.120693],[2.19779,0.6845683,0.3705],[2.19356,0.6828167,-0.3705],[2.36828,0.4318616,1.120693],[2.39738,0.4628343,1.116097],[2.34738,0.6845683,1.116097],[2.2747,0.6828167,1.120693],[2.4345,0.6726084,-1.116097],[2.432163,0.676253,-1.120057],[2.19779,0.6726084,1.127193],[2.194147,0.676253,1.124653],[2.18758,0.5138201,-1.120693],[2.18583,0.5198001,-1.116097],[2.18408,0.6208467,0.3705],[2.18408,0.5733134,-0.3705],[2.311347,1.097867,0.505],[2.295094,1.045733,0.505],[2.31564,1.050333,0.7366667],[2.31564,0.9092667,0.8283334],[1.5341,0.4216667,-0.559834],[1.5341,0.3787333,-0.711167],[1.349527,0.3319153,0.310197],[1.338313,0.3573617,0.298849],[1.511673,0.3319153,-0.2648043],[1.522887,0.3573617,-0.2761527],[1.3271,0.3919563,0.195486],[1.3271,0.4481127,0.224471],[1.349527,0.303023,-0.3795157],[1.338313,0.339246,-0.3505304],[-1.3179,0.4216667,-0.8625001],[-1.3869,0.3787333,-0.8625001],[-1.2489,0.4216667,0.711167],[-1.2489,0.3787333,0.559834],[-1.271327,0.3319153,0.310197],[-1.260113,0.3573617,0.298849],[-1.2489,0.4208487,-0.310197],[-1.2489,0.4662283,-0.298849],[-1.2489,0.4481127,-0.35053],[-1.2489,0.3919563,-0.379515],[-1.260113,0.339246,0.224471],[-1.271327,0.303023,0.195486],[-2.310733,2.4104,0.2025],[-2.303067,2.4104,-0.2025],[-2.337413,0.7176,0.4122334],[-2.342167,0.7176,0.5098667],[1.0028,1.2788,1.066667],[1.0028,1.288,1.078333],[0.138,2.100054,1.066667],[0.138,1.859627,1.078333],[0.138,0.7421333,-1.07],[0.138,1.005867,-1.085],[1.0028,2.100054,-1.07],[1.0028,1.859627,-1.085],[1.54012,1.89226,-1.20977],[1.541303,1.8923,-1.215717],[1.541323,1.879787,-1.38969],[1.540123,1.8859,-1.3891],[1.52904,1.728373,-1.396],[1.52582,1.804427,-1.396],[1.58746,1.595033,-1.144133],[1.58746,1.585037,-1.209133],[2.319334,0.966,1.065687],[2.323467,0.69,1.054843],[2.31486,1.257613,-1.075687],[2.32123,1.2498,-1.064843],[-2.296027,0.6992,-1.09102],[-2.286053,0.9844,-1.09551],[-2.13456,0.6992,-1.1],[2.061907,0.46321,-1.1],[2.170707,0.73921,-1.1],[-0.6666933,1.005867,-1.1],[-1.94295,0.7912,-1.1],[1.430293,1.125467,-1.1],[2.038363,0.8342767,-1.1],[1.992327,0.9195113,-1.1],[1.927223,0.989105,-1.1],[1.84749,1.038315,-1.1],[-1.888223,0.9043241,-1.1],[-1.870387,0.9433564,-1.1],[-1.83743,0.9732304,-1.1],[1.758557,1.063788,-1.1],[-1.794367,0.989398,-1.1],[-1.658033,0.9936001,-1.1],[-0.9615116,0.6191241,-1.1],[1.666487,1.063788,-1.1],[-1.523487,0.989398,-1.1],[-1.08433,0.8636647,-1.1],[-0.6790267,0.709497,-1.1],[-0.64859,0.679623,-1.1],[-0.5785917,0.5485907,-1.1],[0.04859366,0.4354667,-1.1],[0.5957187,0.4569333,-1.1],[0.9017853,0.5714667,-1.1],[0.9500667,0.7309334,-1.1],[1.072987,0.875206,-1.1],[1.150267,0.9228606,-1.1],[1.525373,1.033121,-1.1],[-2.310733,1.0856,-0.6075],[-2.303067,0.9016001,-0.6075],[2.297707,0.414,-1.082843],[2.311813,0.414,-1.06751],[-2.30703,2.561623,1.056143],[-2.30077,2.567277,1.061723],[-2.312863,2.542117,-1.066243],[-2.308127,2.549043,-1.076413],[1.29738,2.58437,0.3438567],[1.29347,2.58437,-0.35401],[1.328463,2.542117,-1.066243],[1.323727,2.549043,-1.076413],[1.534487,1.637447,1.384403],[1.536717,1.636273,1.37903],[1.53449,1.899087,1.32191],[1.53672,1.898373,1.269667],[1.52904,1.636773,1.2071],[1.52582,1.634887,1.21205],[1.56906,1.609163,1.26324],[1.57826,1.607097,1.25782],[1.58175,1.78461,1.308],[1.583043,1.721197,1.308],[1.573843,1.725613,1.308],[1.56906,1.660907,1.308],[-1.80143,0.506,-1.087333],[-1.75134,0.598,-1.074667],[0.9240373,0.5500001,-1.087333],[0.9901187,0.686,-1.074667],[1.859983,0.738578,1.077333],[1.831717,0.777956,1.064667],[-1.078417,0.8190627,1.077333],[-1.096663,0.8310283,1.064667],[1.934453,0.6258867,1.052],[1.92628,0.6901433,1.052],[-1.031427,0.7896804,1.052],[-1.039743,0.8259887,1.052],[-1.808107,0.5244,-1.062],[-1.76269,0.6348,-1.062],[0.913028,0.5684,-1.062],[0.972066,0.7228,-1.062],[2.17829,0.4692,0.8816468],[2.078047,0.4692,0.7112934],[-0.8858174,0.5612,0.8816468],[-0.9414867,0.6532,0.7112934],[-2.269334,0.4692,0.17698],[-2.255227,0.4692,-0.3573334],[-0.9414867,0.6532,-0.8916467],[-0.8858174,0.5612,-0.7212933],[1.283969,0.6052,0.54094],[1.390677,0.7709394,0.54094],[1.703503,0.6841494,0.54094],[1.46647,0.8678026,0.54094],[1.7238,0.79713,0.54094],[1.541687,0.9380534,0.54094],[1.723803,0.8888447,0.54094],[1.609323,0.9774314,0.54094],[1.70359,0.9519584,0.54094],[-2.347093,0.382846,-1.125817],[-2.265563,0.377052,-1.124643],[-2.345573,0.468825,-1.118243],[-2.268103,0.467019,-1.118107],[-2.423167,0.468825,0.3726567],[-2.423043,0.467019,-0.3727934],[-2.267747,0.39641,1.124527],[-2.348237,0.40418,1.122063],[2.18408,0.5138201,-0.3705],[2.18408,0.5198001,0.3705],[2.240677,0.395913,-1.127193],[2.281887,0.399599,-1.129097],[2.240677,0.395913,1.127193],[2.281887,0.399599,1.129097],[2.42429,0.6828167,-1.116097],[2.427933,0.6804814,-1.120057],[2.18758,0.6726084,1.116097],[2.189917,0.676253,1.120057],[2.321554,0.3893667,1.116097],[2.326714,0.391158,1.120203],[2.2747,0.68632,0.3705],[2.34738,0.68632,-0.3705],[2.311347,1.15,-0.1683333],[2.295094,1.15,0.1683333],[2.323613,1.050333,0.92],[2.319627,0.9092667,0.92],[1.3271,0.3787333,0.711167],[1.3271,0.4216667,0.559834],[1.453887,0.2668,0.711167],[1.407313,0.2668,0.559834],[1.3271,0.4662283,-0.2761527],[1.3271,0.4208487,-0.2648043],[1.5341,0.4662283,0.298849],[1.5341,0.4208487,0.310197],[1.522887,0.339246,-0.224471],[1.511673,0.303023,-0.195486],[1.453887,0.293246,-0.224471],[1.407313,0.280023,-0.195486],[-1.2489,0.4216667,0.05550067],[-1.2489,0.3787333,-0.05549967],[-1.444687,0.3573617,0.2761527],[-1.433473,0.3319153,0.2648043],[-1.260113,0.3573617,-0.298849],[-1.271327,0.3319153,-0.310197],[-1.4559,0.3919563,0.3795156],[-1.4559,0.4481127,0.3505304],[-1.433473,0.303023,-0.1954853],[-1.444687,0.339246,-0.2244707],[-2.310733,1.882933,0.6075],[-2.303067,2.146667,0.6075],[-2.337413,0.7176,-0.5098667],[-2.342167,0.7176,-0.4122334],[1.0028,1.2512,1.066667],[1.0028,1.2604,1.078333],[0.1410813,2.355973,1.078333],[0.1395407,2.348227,1.066667],[0.3557333,0.4784,-1.07],[0.5734667,0.4784,-1.085],[0.9997187,2.355973,-1.07],[1.001259,2.348227,-1.085],[1.534487,1.888913,-1.200597],[1.53672,1.883727,-1.1991],[1.534487,1.643887,-1.200597],[1.53672,1.64907,-1.1991],[1.52904,1.8998,-1.322333],[1.52582,1.8998,-1.269667],[1.56906,1.57504,-1.155],[1.57826,1.57504,-1.22],[2.295254,0.966,1.08551],[2.305227,0.69,1.08102],[2.305227,0.966,-1.09102],[2.295254,0.69,-1.09551],[-2.314267,0.6992,-1.064843],[-2.310133,0.9844,-1.075687],[-2.13456,0.6992,1.09],[2.061907,0.46321,1.09],[2.170707,0.73921,1.09],[-0.6666933,1.005867,1.09],[-1.94295,0.7912,1.09],[1.430293,1.125467,1.09],[2.038363,0.8342767,1.09],[1.992327,0.9195113,1.09],[1.927223,0.989105,1.09],[1.84749,1.038315,1.09],[-1.888223,0.904324,1.09],[-1.870387,0.9433564,1.09],[-1.83743,0.9732304,1.09],[1.758557,1.063788,1.09],[-1.794367,0.989398,1.09],[-1.658033,0.9936001,1.09],[-0.9615116,0.6191241,1.09],[1.666487,1.063788,1.09],[-1.523487,0.989398,1.09],[-1.08433,0.8636646,1.09],[-0.6790267,0.709497,1.09],[-0.64859,0.679623,1.09],[-0.5785917,0.5485907,1.09],[0.04859366,0.4354667,1.09],[0.5957187,0.4569333,1.09],[0.9017854,0.5714667,1.09],[0.9500667,0.7309334,1.09],[1.072987,0.875206,1.09],[1.150267,0.9228606,1.09],[1.525373,1.033121,1.09],[-2.310133,1.92453,-1.075687],[-2.314267,2.22986,-1.064843],[-2.302613,0.414,-1.06751],[-2.288507,0.414,-1.082843],[-2.293217,2.561623,-1.084953],[-2.28441,2.567277,-1.084417],[1.320197,2.561623,-1.07255],[1.319703,2.567277,-1.062953],[1.29361,2.58437,-1.055983],[1.296417,2.58437,-1.05904],[1.513867,1.92453,-1.064843],[1.417733,2.22986,-1.075687],[1.602494,1.644113,1.013167],[1.394467,2.334507,0.0265],[1.387883,2.356353,-0.001666668],[1.463657,2.104887,-0.001666667],[1.526263,1.897113,0.001666667],[1.59545,1.667493,-0.0265],[1.60203,1.645647,0.001666668],[1.533307,1.873733,1.038],[1.44835,2.155684,1.041],[1.372117,2.408683,1.013167],[1.365533,2.43053,0.6946666],[1.349767,2.48286,0.02316667],[1.365533,2.43053,-0.3513333],[1.365533,2.43053,-0.698],[1.372117,2.408684,-1.0165],[1.44835,2.155684,-1.047667],[1.533307,1.873733,-1.041333],[1.602494,1.644113,-1.0165],[1.609073,1.622267,0.6946667],[1.609073,1.622267,-0.6980001],[1.609537,1.620733,-0.02983333],[1.609073,1.622267,0.348],[1.54012,1.64054,1.38223],[1.541303,1.6405,1.376283],[1.54013,1.895853,1.32084],[1.54131,1.894047,1.27009],[1.52904,1.647767,1.19805],[1.52582,1.643213,1.2001],[1.56906,1.621203,1.273153],[1.57826,1.616217,1.270907],[1.580457,1.787733,1.297133],[1.577333,1.789027,1.286267],[0.3557333,0.7421333,1.055],[0.2238667,2.31104,1.055],[0.644,0.8801333,1.055],[0.7145333,1.134667,1.055],[0.4262667,1.376933,1.055],[0.1891343,2.340167,1.055],[0.4327067,2.339867,1.055],[0.6535067,2.370533,1.055],[0.8893334,2.31104,1.055],[0.9378654,2.340167,1.055],[0.952177,2.33424,1.055],[0.1748223,2.33424,1.055],[0.963131,2.323287,1.055],[0.9368593,2.089707,1.055],[0.1638687,2.323287,1.055],[0.1579407,2.308973,1.055],[0.1564,2.0608,1.055],[0.1748,1.841533,1.055],[0.2036267,1.632693,1.055],[0.4734933,1.4858,1.055],[0.7007334,1.478133,1.055],[0.9095733,1.50696,1.055],[0.9690593,1.7628,1.055],[1.001259,2.1078,1.055],[1.001259,1.76004,1.055],[1.001259,1.643507,1.055],[1.001259,1.625107,1.055],[-1.691977,0.732248,1.077333],[-1.696613,0.7111241,1.064667],[-0.9572317,0.598,-1.087333],[-0.9037034,0.506,-1.074667],[1.953107,0.46321,1.077333],[1.947333,0.5124201,1.064667],[-1.01504,0.7111241,1.077333],[-1.01932,0.732248,1.064667],[2.192397,0.4508,1.052],[2.081143,0.4324,1.052],[-0.8970277,0.5244,1.052],[-0.9460213,0.6348,1.052],[-2.28344,0.4324,0.3473334],[-2.28344,0.4508,-0.3573334],[-0.8970277,0.5244,-1.062],[-0.9460213,0.6348,-1.062],[2.17829,0.4692,-0.8916467],[2.078047,0.4692,-0.7212934],[0.9650227,0.7412,0.8816468],[0.8949754,0.6052001,0.7112934],[-1.767363,0.6532,0.8816468],[-1.819457,0.5612,0.7112934],[-1.03587,0.8443887,-0.8916467],[-1.02319,0.8264803,-0.7212933],[1.213922,0.4692,-0.18698],[1.602916,0.4692,0.17698],[-2.184033,0.382846,-0.37449],[-2.185177,0.377052,0.3756633],[-2.424003,0.4136897,-0.3727933],[-2.423647,0.4154293,0.3732],[-2.189797,0.4136897,0.3727933],[-2.190154,0.4154293,-0.3732],[-2.185073,0.39641,-0.3732001],[-2.187257,0.40418,0.3756634],[2.18758,0.43286,1.120693],[2.18583,0.46138,1.116097],[2.19356,0.43286,-1.127193],[2.19779,0.46138,-1.129097],[2.35765,0.4327737,-1.129097],[2.38944,0.4632903,-1.127193],[2.19779,0.6828167,-1.116097],[2.194147,0.6804814,-1.120057],[2.18758,0.4001117,1.116097],[2.189917,0.396467,1.120057],[2.338437,0.3970927,1.116097],[2.333707,0.3943583,1.120203],[2.39913,0.4628344,0.3705],[2.37199,0.4314057,-0.3705],[2.311347,0.7406,0.1683333],[2.295094,0.7406,-0.1683333],[2.323613,0.9092667,0.645],[2.319627,1.050333,0.645],[1.3271,0.3787333,-0.559834],[1.3271,0.4216667,-0.711167],[1.453887,0.2668,-0.559834],[1.407313,0.2668,-0.711167],[1.453887,0.3113617,-0.298849],[1.407313,0.3089153,-0.310197],[1.4651,0.5091617,0.298849],[1.3961,0.5067154,0.310197],[1.3961,0.491046,-0.224471],[1.4651,0.477823,-0.195486],[1.3961,0.491046,0.3505304],[1.4651,0.477823,0.3795157],[-1.4559,0.3787333,0.711167],[-1.4559,0.4216667,0.559834],[-1.329113,0.2668,0.711167],[-1.375687,0.2668,0.559834],[-1.4559,0.4662283,-0.2761527],[-1.4559,0.4208487,-0.2648043],[-1.2489,0.4662283,0.298849],[-1.2489,0.4208487,0.310197],[-1.260113,0.339246,-0.2244707],[-1.271327,0.303023,-0.1954853],[-1.329113,0.293246,-0.2244707],[-1.375687,0.280023,-0.1954853],[-2.337413,0.5888,-0.6075],[-2.342167,0.6532,-0.6075],[-2.342167,0.5888,-0.3146],[-2.337413,0.6532,-0.3146],[0.1689167,2.392193,1.078333],[0.1623503,2.387807,1.066667],[1.0028,1.125467,-1.07],[1.0028,1.008933,-1.085],[0.9718827,2.392193,-1.07],[0.9784493,2.387807,-1.085],[1.534493,1.803793,-1.196553],[1.536727,1.72837,-1.197107],[1.534487,1.895353,-1.384403],[1.536717,1.896527,-1.37903],[1.52582,1.896027,-1.2071],[1.52904,1.897913,-1.21205],[1.55986,1.611207,-1.28402],[1.55986,1.596223,-1.29264],[2.292483,1.27094,1.08102],[2.276477,1.276523,1.08551],[1.808583,1.506863,-1.09551],[2.049487,1.39441,-1.09102],[2.308823,1.27626,0.3446667],[2.314857,1.27083,-0.3546667],[-0.6666933,1.502667,-1.1],[-1.471387,1.386133,-1.1],[-2.3184,1.882933,-0.7563334],[-2.3184,2.188263,-0.9051667],[-2.3184,1.882933,0.753],[-2.3184,2.188263,0.8985],[-2.3184,2.493593,0.1991667],[-2.3184,2.451997,-0.3513333],[-2.3184,1.0856,-0.7563334],[-2.3184,0.8004,-0.9051666],[-2.3184,1.0856,0.753],[-2.3184,0.8004,0.8985001],[-2.3184,0.5152,0.1991667],[-2.3184,0.6164,-0.3513333],[-2.294603,2.576744,-1.060897],[-2.287827,2.580557,-1.06257],[1.325963,2.561623,0.34524],[1.322137,2.567277,-0.3560667],[0.10172,2.58437,-1.06009],[-1.08626,2.58437,-1.05583],[1.328463,2.542117,1.056243],[1.323727,2.549043,1.066413],[1.540153,1.730603,1.392937],[1.541347,1.803793,1.391533],[1.54013,1.636947,1.27116],[1.54131,1.638753,1.32191],[1.52582,1.647767,1.39395],[1.52904,1.643213,1.3919],[1.56906,1.577107,1.29042],[1.57826,1.579173,1.29584],[1.57826,1.60503,1.144133],[1.56906,1.60503,1.198267]],"normals":[[0.9823695,-0.07654645,0.1705602],[0.9823695,-0.07654645,0.1705602],[0.9823695,-0.07654645,0.1705602],[0.9672287,0,0.253907],[0.9672287,0,0.253907],[0.9672287,0,0.253907],[0.9823695,0.07654645,0.1705602],[0.9823695,0.07654645,0.1705602],[0.9823695,0.07654645,0.1705602],[0.9855622,0.1582176,0.06028473],[0.9855622,0.1582176,0.06028473],[0.9855622,0.1582176,0.06028473],[0.9811877,0.1930563,0],[0.9811877,0.1930563,0],[0.9811877,0.1930563,0],[0.9855622,0.1582176,-0.06028473],[0.9855622,0.1582176,-0.06028473],[0.9855622,0.1582176,-0.06028473],[0.9823695,0.07654645,-0.1705602],[0.9823695,0.07654645,-0.1705602],[0.9823695,0.07654645,-0.1705602],[0.9672287,0,-0.253907],[0.9672287,0,-0.253907],[0.9672287,0,-0.253907],[0.9823695,-0.07654645,-0.1705602],[0.9823695,-0.07654645,-0.1705602],[0.9823695,-0.07654645,-0.1705602],[0.9855622,-0.1582184,-0.06028366],[0.9855622,-0.1582184,-0.06028366],[0.9855622,-0.1582184,-0.06028366],[0.9811876,-0.1930564,0],[0.9811876,-0.1930564,0],[0.9811876,-0.1930564,0],[0.9855622,-0.1582184,0.06028366],[0.9855622,-0.1582184,0.06028366],[0.9855622,-0.1582184,0.06028366],[0.9823695,-0.07654645,-0.1705603],[0.9823695,-0.07654645,-0.1705603],[0.9823695,-0.07654645,-0.1705603],[0.9672287,0,-0.253907],[0.9672287,0,-0.253907],[0.9672287,0,-0.253907],[0.9823695,0.07654645,-0.1705603],[0.9823695,0.07654645,-0.1705603],[0.9823695,0.07654645,-0.1705603],[0.9855623,0.1582176,-0.0602847],[0.9855623,0.1582176,-0.0602847],[0.9855623,0.1582176,-0.0602847],[0.9811877,0.1930563,0],[0.9811877,0.1930563,0],[0.9811877,0.1930563,0],[0.9855622,0.1582176,0.06028468],[0.9855622,0.1582176,0.06028468],[0.9855622,0.1582176,0.06028468],[0.9823695,0.07654645,0.1705603],[0.9823695,0.07654645,0.1705603],[0.9823695,0.07654645,0.1705603],[0.9672287,-0,0.253907],[0.9672287,-0,0.253907],[0.9672287,-0,0.253907],[0.9823695,-0.07654645,0.1705603],[0.9823695,-0.07654645,0.1705603],[0.9823695,-0.07654645,0.1705603],[0.9855622,-0.1582184,0.06028364],[0.9855622,-0.1582184,0.06028364],[0.9855622,-0.1582184,0.06028364],[0.9811876,-0.1930564,0],[0.9811876,-0.1930564,0],[0.9811876,-0.1930564,0],[0.9855622,-0.1582184,-0.06028363],[0.9855622,-0.1582184,-0.06028363],[0.9855622,-0.1582184,-0.06028363],[0.8049616,-0.5933269,0],[0.8049616,-0.5933269,0],[0.8049616,-0.5933269,0],[0.8049616,-0.5933269,0],[0.8049616,-0.5933269,0],[0.8049616,-0.5933269,0],[-0.980085,-0.1985785,0],[-0.980085,-0.1985785,0],[-0.980085,-0.1985785,0],[-0.980085,-0.1985785,0],[-0.980085,-0.1985785,0],[-0.980085,-0.1985785,0],[-0.9931883,-0.1165206,0],[-0.9931883,-0.1165206,0],[-0.9931883,-0.1165206,0],[-0.9931883,-0.1165206,0],[-0.9931883,-0.1165206,0],[-0.9931883,-0.1165206,0],[-0.8643441,-0.5029009,0],[-0.8643441,-0.5029009,0],[-0.8643441,-0.5029009,0],[-0.8643441,-0.5029009,0],[-0.8643441,-0.5029009,0],[-0.8643441,-0.5029009,0],[-1,0,0],[-1,0,0],[-1,0,0],[-1,0,0],[-1,0,0],[-1,0,0],[-0,0,-1],[-0,0,-1],[-0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[-0,0,-1],[-0,0,-1],[-0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,1],[0,0,1],[0,0,1],[-0,0,0.9999999],[-0,0,0.9999999],[-0,0,0.9999999],[-0.9923552,-0.1234143,0],[-0.9923552,-0.1234143,0],[-0.9923552,-0.1234143,0],[-0.9923552,-0.1234143,0],[-0.9923552,-0.1234143,0],[-0.9923552,-0.1234143,0],[0.6299935,-0.7766004,0],[0.6299935,-0.7766004,0],[0.6299935,-0.7766004,0],[0.6299935,-0.7766004,0],[0.6299935,-0.7766004,0],[0.6299935,-0.7766004,0],[0.9749084,-0.2226065,0],[0.9749084,-0.2226065,0],[0.9749084,-0.2226065,0],[0.9749084,-0.2226065,0],[0.9749084,-0.2226065,0],[0.9749084,-0.2226065,0],[-0.5334361,-0.8458403,0],[-0.5334361,-0.8458403,0],[-0.5334361,-0.8458403,0],[-0.5334361,-0.8458403,0],[-0.5334361,-0.8458403,0],[-0.5334361,-0.8458403,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0.42555,-0.9049349,0],[0.42555,-0.9049349,0],[0.42555,-0.9049349,0],[0.42555,-0.9049349,0],[0.42555,-0.9049349,0],[0.42555,-0.9049349,0],[-0.9976759,0.06813686,0],[-0.9976759,0.06813686,0],[-0.9976759,0.06813686,0],[-0.9976759,0.06813686,0],[-0.9976759,0.06813686,0],[-0.9976759,0.06813686,0],[0.9976759,0.06813686,0],[0.9976759,0.06813686,0],[0.9976759,0.06813686,0],[0.9976759,0.06813686,-0],[0.9976759,0.06813686,-0],[0.9976759,0.06813686,-0],[0.962717,-0.2705105,0],[0.962717,-0.2705105,0],[0.962717,-0.2705105,0],[0.962717,-0.2705105,0],[0.962717,-0.2705105,0],[0.962717,-0.2705105,0],[-0.4103398,0,0.9119328],[-0.4103398,0,0.9119328],[-0.4103398,0,0.9119328],[-0.4103398,0,0.9119328],[-0.4103398,0,0.9119328],[-0.4103398,0,0.9119328],[-0.9345582,-0,-0.3558106],[-0.9345582,-0,-0.3558106],[-0.9345582,-0,-0.3558106],[-0.9345582,0,-0.3558106],[-0.9345582,0,-0.3558106],[-0.9345582,0,-0.3558106],[0.6994621,-0.6176515,-0.3595267],[0.6994621,-0.6176515,-0.3595267],[0.6994621,-0.6176515,-0.3595267],[0.727269,-0.6280276,-0.2768777],[0.727269,-0.6280276,-0.2768777],[0.727269,-0.6280276,-0.2768777],[-0.8716245,0.360758,-0.3318501],[-0.8716245,0.360758,-0.3318501],[-0.8716245,0.360758,-0.3318501],[-0.9676845,0.1855876,-0.1707163],[-0.9676845,0.1855876,-0.1707163],[-0.9676845,0.1855876,-0.1707163],[-0.3796015,-0.379746,0.8436206],[-0.3796015,-0.379746,0.8436206],[-0.3796015,-0.379746,0.8436206],[-0.1996425,-0.1997184,0.9592994],[-0.1996425,-0.1997184,0.9592994],[-0.1996425,-0.1997184,0.9592994],[0.8178453,-0.4839251,0.3113608],[0.8178453,-0.4839251,0.3113608],[0.8178453,-0.4839251,0.3113608],[0.8178453,-0.4839251,0.3113608],[0.8178453,-0.4839251,0.3113608],[0.8178453,-0.4839251,0.3113608],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,1,-0],[0,1,-0],[0,1,-0],[0,1,0],[0,1,0],[0,1,0],[0,1,0],[0,1,0],[0,1,0],[0,1,0],[0,1,0],[0,1,0],[-0.8988641,-0.4382276,0],[-0.8988641,-0.4382276,0],[-0.8988641,-0.4382276,0],[-0.8988641,-0.4382276,0],[-0.8988641,-0.4382276,0],[-0.8988641,-0.4382276,0],[0,-0.9775449,-0.210727],[0,-0.9775449,-0.210727],[0,-0.9775449,-0.210727],[-0,-0.9775449,-0.210727],[-0,-0.9775449,-0.210727],[-0,-0.9775449,-0.210727],[0,0.9775449,-0.210727],[0,0.9775449,-0.210727],[0,0.9775449,-0.210727],[0,0.9775449,-0.210727],[0,0.9775449,-0.210727],[0,0.9775449,-0.210727],[0,0.9097998,-0.4150473],[0,0.9097998,-0.4150473],[0,0.9097998,-0.4150473],[0,0.9097998,-0.4150473],[0,0.9097998,-0.4150473],[0,0.9097998,-0.4150473],[0,0.909798,-0.4150513],[0,0.909798,-0.4150513],[0,0.909798,-0.4150513],[0,0.909798,-0.4150513],[0,0.909798,-0.4150513],[0,0.909798,-0.4150513],[-1,-0,0],[-1,-0,0],[-1,-0,0],[-1,0,0],[-1,0,0],[-1,0,0],[0,-0.9999999,0],[0,-0.9999999,0],[0,-0.9999999,0],[0,-0.9999999,0],[0,-0.9999999,0],[0,-0.9999999,0],[-1,0,-0],[-1,0,-0],[-1,0,-0],[-1,0,0],[-1,0,0],[-1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,-0,0],[1,-0,0],[1,-0,0],[0.8814223,-0.4297241,0.1960408],[0.8814223,-0.4297241,0.1960408],[0.8814223,-0.4297241,0.1960408],[0.8814223,-0.4297241,0.1960408],[0.8814223,-0.4297241,0.1960408],[0.8814223,-0.4297241,0.1960408],[0,-0.909798,0.4150513],[0,-0.909798,0.4150513],[0,-0.909798,0.4150513],[0,-0.909798,0.4150513],[0,-0.909798,0.4150513],[0,-0.909798,0.4150513],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[-1,0,0],[-1,0,0],[-1,0,0],[-1,0,0],[-1,0,0],[-1,0,0],[0.1950071,-0.9808018,0],[0.1950071,-0.9808018,0],[0.1950071,-0.9808018,0],[0.1950071,-0.9808018,0],[0.1950071,-0.9808018,0],[0.1950071,-0.9808018,0],[-1,0,0],[-1,0,0],[-1,0,0],[-1,0,0],[-1,0,0],[-1,0,0],[-0.1950071,-0.9808018,-4.391028e-06],[-0.1950071,-0.9808018,-4.391028e-06],[-0.1950071,-0.9808018,-4.391028e-06],[-0.1950153,-0.9808002,0],[-0.1950153,-0.9808002,0],[-0.1950153,-0.9808002,0],[0.7618072,-0,0.6478036],[0.7618072,-0,0.6478036],[0.7618072,-0,0.6478036],[0.7618073,0,0.6478037],[0.7618073,0,0.6478037],[0.7618073,0,0.6478037],[0.8662896,0.4669726,-0.1774228],[0.8662896,0.4669726,-0.1774228],[0.8662896,0.4669726,-0.1774228],[0.868475,0.4641962,-0.1739923],[0.868475,0.4641962,-0.1739923],[0.868475,0.4641962,-0.1739923],[0,0.4105328,0.9118459],[0,0.4105328,0.9118459],[0,0.4105328,0.9118459],[-0,0.4105328,0.9118459],[-0,0.4105328,0.9118459],[-0,0.4105328,0.9118459],[-0.9999999,-0,0],[-0.9999999,-0,0],[-0.9999999,-0,0],[-1,0,0],[-1,0,0],[-1,0,0],[0,1,-0],[0,1,-0],[0,1,-0],[0,1,0],[0,1,0],[0,1,0],[0.4318366,0.8868036,-0.164611],[0.4318366,0.8868036,-0.164611],[0.4318366,0.8868036,-0.164611],[0.4322751,0.8873095,-0.160686],[0.4322751,0.8873095,-0.160686],[0.4322751,0.8873095,-0.160686],[0.9612595,0.2756446,0],[0.9612595,0.2756446,0],[0.9612595,0.2756446,0],[0.9612595,0.2756446,0],[0.9612595,0.2756446,0],[0.9612595,0.2756446,0],[-1,-0,0],[-1,-0,0],[-1,-0,0],[-1,0,0],[-1,0,0],[-1,0,0],[-0.4105169,0,0.911853],[-0.4105169,0,0.911853],[-0.4105169,0,0.911853],[-0.4105169,0,0.911853],[-0.4105169,0,0.911853],[-0.4105169,0,0.911853],[0,1,0],[0,1,0],[0,1,0],[0,1,0],[0,1,0],[0,1,0],[-0.5172721,0.8329639,0.1964703],[-0.5172721,0.8329639,0.1964703],[-0.5172721,0.8329639,0.1964703],[-0.5171735,0.83301,0.1965351],[-0.5171735,0.83301,0.1965351],[-0.5171735,0.83301,0.1965351],[-0.5578406,0.8299481,0],[-0.5578406,0.8299481,0],[-0.5578406,0.8299481,0],[-0.5578406,0.8299481,0],[-0.5578406,0.8299481,0],[-0.5578406,0.8299481,0],[0,1,0],[0,1,0],[0,1,0],[0,1,0],[0,1,0],[0,1,0],[0.8994444,0.2710146,0.3428571],[0.8994444,0.2710146,0.3428571],[0.8994444,0.2710146,0.3428571],[0.8994444,0.2710146,0.3428571],[0.8994444,0.2710146,0.3428571],[0.8994444,0.2710146,0.3428571],[0.3403893,0.8786296,-0.3348809],[0.3403893,0.8786296,-0.3348809],[0.3403893,0.8786296,-0.3348809],[0.3440155,0.8780355,-0.3327264],[0.3440155,0.8780355,-0.3327264],[0.3440155,0.8780355,-0.3327264],[0.3112625,0.3901391,0.866549],[0.3112625,0.3901391,0.866549],[0.3112625,0.3901391,0.866549],[0.3229925,0.3827518,0.86555],[0.3229925,0.3827518,0.86555],[0.3229925,0.3827518,0.86555],[0,-1,-0],[0,-1,-0],[0,-1,-0],[0,-1,0],[0,-1,0],[0,-1,0],[-1,0,0],[-1,0,0],[-1,0,0],[-1,-0,0],[-1,-0,0],[-1,-0,0],[0.1774316,-0.9841331,0],[0.1774316,-0.9841331,0],[0.1774316,-0.9841331,0],[0.1774316,-0.9841331,0],[0.1774316,-0.9841331,0],[0.1774316,-0.9841331,0],[-0.5483773,-0.836231,0],[-0.5483773,-0.836231,0],[-0.5483773,-0.836231,0],[-0.5483773,-0.836231,0],[-0.5483773,-0.836231,0],[-0.5483773,-0.836231,0],[-0.8123698,-0.5831425,0],[-0.8123698,-0.5831425,0],[-0.8123698,-0.5831425,0],[-0.8123698,-0.5831425,0],[-0.8123698,-0.5831425,0],[-0.8123698,-0.5831425,0],[0.8994453,-0.4370334,0],[0.8994453,-0.4370334,0],[0.8994453,-0.4370334,0],[0.8994453,-0.4370334,0],[0.8994453,-0.4370334,0],[0.8994453,-0.4370334,0],[0,0,1],[0,0,1],[0,0,1],[-0,0,1],[-0,0,1],[-0,0,1],[0,0,-1],[0,0,-1],[0,0,-1],[-0,0,-0.9999999],[-0,0,-0.9999999],[-0,0,-0.9999999],[0,0,-1],[0,0,-1],[0,0,-1],[-0,0,-0.9999999],[-0,0,-0.9999999],[-0,0,-0.9999999],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[-0.7958434,-0.6055024,0],[-0.7958434,-0.6055024,0],[-0.7958434,-0.6055024,0],[-0.7958434,-0.6055024,0],[-0.7958434,-0.6055024,0],[-0.7958434,-0.6055024,0],[0.1316319,-0.9912987,0],[0.1316319,-0.9912987,0],[0.1316319,-0.9912987,0],[0.1316319,-0.9912987,0],[0.1316319,-0.9912987,0],[0.1316319,-0.9912987,0],[0.5032468,-0.8641427,0],[0.5032468,-0.8641427,0],[0.5032468,-0.8641427,0],[0.5032468,-0.8641427,0],[0.5032468,-0.8641427,0],[0.5032468,-0.8641427,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,0.9170429,-0.3987886],[0,0.9170429,-0.3987886],[0,0.9170429,-0.3987886],[0,0.9170429,-0.3987886],[0,0.9170429,-0.3987886],[0,0.9170429,-0.3987886],[0,-0.9999999,0],[0,-0.9999999,0],[0,-0.9999999,0],[0,-0.9999999,-0],[0,-0.9999999,-0],[0,-0.9999999,-0],[0.9796152,-0.2008836,0],[0.9796152,-0.2008836,0],[0.9796152,-0.2008836,0],[0.9796151,-0.2008836,0],[0.9796151,-0.2008836,0],[0.9796151,-0.2008836,0],[-0.962717,-0.2705105,0],[-0.962717,-0.2705105,0],[-0.962717,-0.2705105,0],[-0.962717,-0.2705105,-0],[-0.962717,-0.2705105,-0],[-0.962717,-0.2705105,-0],[-0.4103397,0,0.9119327],[-0.4103397,0,0.9119327],[-0.4103397,0,0.9119327],[-0.4103397,0,0.9119327],[-0.4103397,0,0.9119327],[-0.4103397,0,0.9119327],[0.3825956,0.9239159,-0],[0.3825956,0.9239159,-0],[0.3825956,0.9239159,-0],[0.3825956,0.9239159,0],[0.3825956,0.9239159,0],[0.3825956,0.9239159,0],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[-0.3796015,-0.379746,-0.8436206],[-0.3796015,-0.379746,-0.8436206],[-0.3796015,-0.379746,-0.8436206],[-0.1996425,-0.1997184,-0.9592994],[-0.1996425,-0.1997184,-0.9592994],[-0.1996425,-0.1997184,-0.9592994],[0.5689232,-0.7584897,-0.3178359],[0.5689232,-0.7584897,-0.3178359],[0.5689232,-0.7584897,-0.3178359],[0.7355877,-0.6612584,-0.1471323],[0.7355877,-0.6612584,-0.1471323],[0.7355877,-0.6612584,-0.1471323],[0.4092226,-0.07374049,0.9094499],[0.4092226,-0.07374049,0.9094499],[0.4092226,-0.07374049,0.9094499],[0.409216,-0.07374512,0.9094525],[0.409216,-0.07374512,0.9094525],[0.409216,-0.07374512,0.9094525],[-0.4103397,0,0.9119327],[-0.4103397,0,0.9119327],[-0.4103397,0,0.9119327],[-0.4103397,0,0.9119327],[-0.4103397,0,0.9119327],[-0.4103397,0,0.9119327],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,1,0],[0,1,0],[0,1,0],[0,1,0],[0,1,0],[0,1,0],[-0.8988641,-0.4382276,0],[-0.8988641,-0.4382276,0],[-0.8988641,-0.4382276,0],[-0.8988641,-0.4382276,0],[-0.8988641,-0.4382276,0],[-0.8988641,-0.4382276,0],[-0.8948799,-0.4362851,-0.09404895],[-0.8948799,-0.4362851,-0.09404895],[-0.8948799,-0.4362851,-0.09404895],[-0.8948799,-0.4362851,-0.09404895],[-0.8948799,-0.4362851,-0.09404895],[-0.8948799,-0.4362851,-0.09404895],[-1,0,-0],[-1,0,-0],[-1,0,-0],[-1,0,0],[-1,0,0],[-1,0,0],[1,-0,0],[1,-0,0],[1,-0,0],[1,0,0],[1,0,0],[1,0,0],[-1,0,-0],[-1,0,-0],[-1,0,-0],[-1,0,0],[-1,0,0],[-1,0,0],[0,1,0],[0,1,0],[0,1,0],[0,1,0],[0,1,0],[0,1,0],[-0.8988641,-0.4382276,0],[-0.8988641,-0.4382276,0],[-0.8988641,-0.4382276,0],[-0.8988641,-0.4382276,0],[-0.8988641,-0.4382276,0],[-0.8988641,-0.4382276,0],[0,-0.9775449,-0.210727],[0,-0.9775449,-0.210727],[0,-0.9775449,-0.210727],[-0,-0.9775449,-0.210727],[-0,-0.9775449,-0.210727],[-0,-0.9775449,-0.210727],[0,0.977545,-0.210727],[0,0.977545,-0.210727],[0,0.977545,-0.210727],[0,0.977545,-0.210727],[0,0.977545,-0.210727],[0,0.977545,-0.210727],[0,0.909798,-0.4150514],[0,0.909798,-0.4150514],[0,0.909798,-0.4150514],[0,0.909798,-0.4150514],[0,0.909798,-0.4150514],[0,0.909798,-0.4150514],[0,0.909798,-0.4150513],[0,0.909798,-0.4150513],[0,0.909798,-0.4150513],[0,0.909798,-0.4150513],[0,0.909798,-0.4150513],[0,0.909798,-0.4150513],[-1,0,0],[-1,0,0],[-1,0,0],[-1,0,0],[-1,0,0],[-1,0,0],[0,1,0],[0,1,0],[0,1,0],[0,1,0],[0,1,0],[0,1,0],[-0.1950071,-0.9808018,5.642004e-06],[-0.1950071,-0.9808018,5.642004e-06],[-0.1950071,-0.9808018,5.642004e-06],[-0.1950153,-0.9808003,0],[-0.1950153,-0.9808003,0],[-0.1950153,-0.9808003,0],[-0.9999999,0,0],[-0.9999999,0,0],[-0.9999999,0,0],[-0.9999999,0,0],[-0.9999999,0,0],[-0.9999999,0,0],[0.1950071,-0.9808019,0],[0.1950071,-0.9808019,0],[0.1950071,-0.9808019,0],[0.1950071,-0.9808019,0],[0.1950071,-0.9808019,0],[0.1950071,-0.9808019,0],[0.8494051,-0.207478,-0.4852462],[0.8494051,-0.207478,-0.4852462],[0.8494051,-0.207478,-0.4852462],[0.8422657,-0.218962,-0.4925892],[0.8422657,-0.218962,-0.4925892],[0.8422657,-0.218962,-0.4925892],[0.8662897,-0.4669719,0.177424],[0.8662897,-0.4669719,0.177424],[0.8662897,-0.4669719,0.177424],[0.8684799,-0.4641892,0.1739858],[0.8684799,-0.4641892,0.1739858],[0.8684799,-0.4641892,0.1739858],[0,0.9344294,-0.3561485],[0,0.9344294,-0.3561485],[0,0.9344294,-0.3561485],[0,0.9344294,-0.3561485],[0,0.9344294,-0.3561485],[0,0.9344294,-0.3561485],[1,0,0],[1,0,0],[1,0,0],[1,0,-0],[1,0,-0],[1,0,-0],[1,-0,0],[1,-0,0],[1,-0,0],[1,0,0],[1,0,0],[1,0,0],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[-0,0,-1],[-0,0,-1],[-0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[-0.9344156,0,-0.3561847],[-0.9344156,0,-0.3561847],[-0.9344156,0,-0.3561847],[-0.9344156,-0,-0.3561847],[-0.9344156,-0,-0.3561847],[-0.9344156,-0,-0.3561847],[0,-1,-0],[0,-1,-0],[0,-1,-0],[0,-1,0],[0,-1,0],[0,-1,0],[-0,0,1],[-0,0,1],[-0,0,1],[0,0,1],[0,0,1],[0,0,1],[-0,0.8511584,0.524909],[-0,0.8511584,0.524909],[-0,0.8511584,0.524909],[0,0.8511584,0.524909],[0,0.8511584,0.524909],[0,0.8511584,0.524909],[-0.8281289,0.5605377,0],[-0.8281289,0.5605377,0],[-0.8281289,0.5605377,0],[-0.8281289,0.5605377,0],[-0.8281289,0.5605377,0],[-0.8281289,0.5605377,0],[-0.7712992,0.5649389,-0.2931583],[-0.7712992,0.5649389,-0.2931583],[-0.7712992,0.5649389,-0.2931583],[-0.7714956,0.5647374,-0.2930295],[-0.7714956,0.5647374,-0.2930295],[-0.7714956,0.5647374,-0.2930295],[0.4074135,0.1227591,0.9049555],[0.4074135,0.1227591,0.9049555],[0.4074135,0.1227591,0.9049555],[0.4074135,0.1227591,0.9049555],[0.4074135,0.1227591,0.9049555],[0.4074135,0.1227591,0.9049555],[0.3112625,0.3901391,-0.866549],[0.3112625,0.3901391,-0.866549],[0.3112625,0.3901391,-0.866549],[0.3229924,0.3827519,-0.8655502],[0.3229924,0.3827519,-0.8655502],[0.3229924,0.3827519,-0.8655502],[0.311262,-0.3901364,-0.8665504],[0.311262,-0.3901364,-0.8665504],[0.311262,-0.3901364,-0.8665504],[0.3225104,-0.3830546,-0.8655959],[0.3225104,-0.3830546,-0.8655959],[0.3225104,-0.3830546,-0.8655959],[0,1,-0],[0,1,-0],[0,1,-0],[0,1,0],[0,1,0],[0,1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0.5180455,-0.855353,0],[0.5180455,-0.855353,0],[0.5180455,-0.855353,0],[0.5180455,-0.855353,0],[0.5180455,-0.855353,0],[0.5180455,-0.855353,0],[-0.8266612,-0.5627,0],[-0.8266612,-0.5627,0],[-0.8266612,-0.5627,0],[-0.8266612,-0.5627,0],[-0.8266612,-0.5627,0],[-0.8266612,-0.5627,0],[-0.9324376,-0.3613309,0],[-0.9324376,-0.3613309,0],[-0.9324376,-0.3613309,0],[-0.9324376,-0.3613309,0],[-0.9324376,-0.3613309,0],[-0.9324376,-0.3613309,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[-0,0,1],[-0,0,1],[-0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,-1],[0,0,-1],[0,0,-1],[-0,0,-1],[-0,0,-1],[-0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[-0,0,-0.9999999],[-0,0,-0.9999999],[-0,0,-0.9999999],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[-0.925001,-0.3799649,0],[-0.925001,-0.3799649,0],[-0.925001,-0.3799649,0],[-0.925001,-0.3799649,0],[-0.925001,-0.3799649,0],[-0.925001,-0.3799649,0],[0.3680699,-0.9297982,0],[0.3680699,-0.9297982,0],[0.3680699,-0.9297982,0],[0.3680699,-0.9297982,0],[0.3680699,-0.9297982,0],[0.3680699,-0.9297982,0],[0.7936012,-0.6084384,0],[0.7936012,-0.6084384,0],[0.7936012,-0.6084384,0],[0.7936012,-0.6084384,0],[0.7936012,-0.6084384,0],[0.7936012,-0.6084384,0],[-0.1845051,-0.9828316,0],[-0.1845051,-0.9828316,0],[-0.1845051,-0.9828316,0],[-0.1845051,-0.9828316,0],[-0.1845051,-0.9828316,0],[-0.1845051,-0.9828316,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,0.1984798,-0.9801049],[0,0.1984798,-0.9801049],[0,0.1984798,-0.9801049],[0,0.1984798,-0.980105],[0,0.1984798,-0.980105],[0,0.1984798,-0.980105],[0,0.9999999,-0],[0,0.9999999,-0],[0,0.9999999,-0],[0,0.9999999,0],[0,0.9999999,0],[0,0.9999999,0],[0.9976759,-0.06813686,0],[0.9976759,-0.06813686,0],[0.9976759,-0.06813686,0],[0.9976759,-0.06813686,0],[0.9976759,-0.06813686,0],[0.9976759,-0.06813686,0],[-0.962714,0.2705211,0],[-0.962714,0.2705211,0],[-0.962714,0.2705211,0],[-0.9627141,0.2705211,0],[-0.9627141,0.2705211,0],[-0.9627141,0.2705211,0],[-0.9345581,0,0.3558105],[-0.9345581,0,0.3558105],[-0.9345581,0,0.3558105],[-0.9345581,0,0.3558105],[-0.9345581,0,0.3558105],[-0.9345581,0,0.3558105],[0.9239908,0.3824147,-0],[0.9239908,0.3824147,-0],[0.9239908,0.3824147,-0],[0.9239908,0.3824147,0],[0.9239908,0.3824147,0],[0.9239908,0.3824147,0],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[-0.3796015,0.379746,-0.8436206],[-0.3796015,0.379746,-0.8436206],[-0.3796015,0.379746,-0.8436206],[-0.1996425,0.1997184,-0.9592994],[-0.1996425,0.1997184,-0.9592994],[-0.1996425,0.1997184,-0.9592994],[0.2019745,-0.9263822,-0.3178397],[0.2019745,-0.9263822,-0.3178397],[0.2019745,-0.9263822,-0.3178397],[0.01942625,-0.9889103,-0.1472382],[0.01942625,-0.9889103,-0.1472382],[0.01942625,-0.9889103,-0.1472382],[0.3987478,-0.2359485,0.8861876],[0.3987478,-0.2359485,0.8861876],[0.3987478,-0.2359485,0.8861876],[0.3987543,-0.2359459,0.8861853],[0.3987543,-0.2359459,0.8861853],[0.3987543,-0.2359459,0.8861853],[-0.9345581,0,0.3558105],[-0.9345581,0,0.3558105],[-0.9345581,0,0.3558105],[-0.9345581,0,0.3558105],[-0.9345581,0,0.3558105],[-0.9345581,0,0.3558105],[0,-1,-0],[0,-1,-0],[0,-1,-0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,-0],[0,-1,-0],[0,-1,-0],[0,-1,0],[0,-1,0],[0,-1,0],[0,1,0],[0,1,0],[0,1,0],[0,1,0],[0,1,0],[0,1,0],[-0.8988641,-0.4382276,0],[-0.8988641,-0.4382276,0],[-0.8988641,-0.4382276,0],[-0.8988641,-0.4382276,0],[-0.8988641,-0.4382276,0],[-0.8988641,-0.4382276,0],[-0.8948797,-0.436285,0.09405172],[-0.8948797,-0.436285,0.09405172],[-0.8948797,-0.436285,0.09405172],[-0.8948797,-0.436285,0.09405172],[-0.8948797,-0.436285,0.09405172],[-0.8948797,-0.436285,0.09405172],[-1,0,0],[-1,0,0],[-1,0,0],[-1,-0,0],[-1,-0,0],[-1,-0,0],[1,0,-0],[1,0,-0],[1,0,-0],[1,0,0],[1,0,0],[1,0,0],[-1,-0,0],[-1,-0,0],[-1,-0,0],[-1,0,0],[-1,0,0],[-1,0,0],[-0.9999999,-0,0],[-0.9999999,-0,0],[-0.9999999,-0,0],[-0.9999999,0,0],[-0.9999999,0,0],[-0.9999999,0,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-0.9775435,0.210733],[0,-0.9775435,0.210733],[0,-0.9775435,0.210733],[0,-0.9775435,0.210733],[0,-0.9775435,0.210733],[0,-0.9775435,0.210733],[-0,0.9775437,0.210733],[-0,0.9775437,0.210733],[-0,0.9775437,0.210733],[0,0.9775437,0.210733],[0,0.9775437,0.210733],[0,0.9775437,0.210733],[0,0.9097998,0.4150473],[0,0.9097998,0.4150473],[0,0.9097998,0.4150473],[0,0.9097998,0.4150473],[0,0.9097998,0.4150473],[0,0.9097998,0.4150473],[0,0.9097998,0.4150473],[0,0.9097998,0.4150473],[0,0.9097998,0.4150473],[0,0.9097998,0.4150473],[0,0.9097998,0.4150473],[0,0.9097998,0.4150473],[-1,0,0],[-1,0,0],[-1,0,0],[-1,0,0],[-1,0,0],[-1,0,0],[0,0,-1],[0,0,-1],[0,0,-1],[-0,0,-1],[-0,0,-1],[-0,0,-1],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[-1,0,0],[-1,0,0],[-1,0,0],[-1,0,0],[-1,0,0],[-1,0,0],[0,-1,-0],[0,-1,-0],[0,-1,-0],[0,-1,0],[0,-1,0],[0,-1,0],[0.311262,-0.3901364,-0.8665504],[0.311262,-0.3901364,-0.8665504],[0.311262,-0.3901364,-0.8665504],[0.3225104,-0.3830546,-0.8655959],[0.3225104,-0.3830546,-0.8655959],[0.3225104,-0.3830546,-0.8655959],[0.3403893,-0.8786296,0.3348809],[0.3403893,-0.8786296,0.3348809],[0.3403893,-0.8786296,0.3348809],[0.3440155,-0.8780355,0.3327264],[0.3440155,-0.8780355,0.3327264],[0.3440155,-0.8780355,0.3327264],[0,0.4105328,-0.9118459],[0,0.4105328,-0.9118459],[0,0.4105328,-0.9118459],[0,0.4105328,-0.9118459],[0,0.4105328,-0.9118459],[0,0.4105328,-0.9118459],[0.9999999,0,0],[0.9999999,0,0],[0.9999999,0,0],[1,0,0],[1,0,0],[1,0,0],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[-0,0,1],[-0,0,1],[-0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[1,0,0],[1,0,0],[1,0,0],[1,0,-0],[1,0,-0],[1,0,-0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,-0],[1,0,-0],[1,0,-0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[0.9999999,-0,0],[0.9999999,-0,0],[0.9999999,-0,0],[1,0,0],[1,0,0],[1,0,0],[0.9999999,-0,0],[0.9999999,-0,0],[0.9999999,-0,0],[0.9999999,-0,0],[0.9999999,-0,0],[0.9999999,-0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[0.9999999,0,0],[0.9999999,0,0],[0.9999999,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,-0],[1,0,-0],[1,0,-0],[1,0,-0],[1,0,-0],[1,0,-0],[-0,0,1],[-0,0,1],[-0,0,1],[0,0,1],[0,0,1],[0,0,1],[-0.9344155,0,0.3561847],[-0.9344155,0,0.3561847],[-0.9344155,0,0.3561847],[-0.9344155,0,0.3561847],[-0.9344155,0,0.3561847],[-0.9344155,0,0.3561847],[0,1,0],[0,1,0],[0,1,0],[0,1,0],[0,1,0],[0,1,0],[-0.218006,0.8479846,0.4831104],[-0.218006,0.8479846,0.4831104],[-0.218006,0.8479846,0.4831104],[-0.2180028,0.8479827,0.4831149],[-0.2180028,0.8479827,0.4831149],[-0.2180028,0.8479827,0.4831149],[-0.9280505,0.116522,0.3537584],[-0.9280505,0.116522,0.3537584],[-0.9280505,0.116522,0.3537584],[-0.9285354,0.1151894,0.3529212],[-0.9285354,0.1151894,0.3529212],[-0.9285354,0.1151894,0.3529212],[0.5578406,0.8299481,0],[0.5578406,0.8299481,0],[0.5578406,0.8299481,0],[0.5578406,0.8299481,-0],[0.5578406,0.8299481,-0],[0.5578406,0.8299481,-0],[0.4074484,0.1220679,0.9050333],[0.4074484,0.1220679,0.9050333],[0.4074484,0.1220679,0.9050333],[0.4083769,0.1234887,0.9044219],[0.4083769,0.1234887,0.9044219],[0.4083769,0.1234887,0.9044219],[0.8662897,0.4669719,-0.177424],[0.8662897,0.4669719,-0.177424],[0.8662897,0.4669719,-0.177424],[0.8684763,0.4641938,-0.1739914],[0.8684763,0.4641938,-0.1739914],[0.8684763,0.4641938,-0.1739914],[0.8492427,0.2075781,0.4854875],[0.8492427,0.2075781,0.4854875],[0.8492427,0.2075781,0.4854875],[0.8422273,0.2188509,0.4927044],[0.8422273,0.2188509,0.4927044],[0.8422273,0.2188509,0.4927044],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[0.1774316,-0.9841331,0],[0.1774316,-0.9841331,0],[0.1774316,-0.9841331,0],[0.1774316,-0.9841331,0],[0.1774316,-0.9841331,0],[0.1774316,-0.9841331,0],[-0.1393552,-0.9902424,0],[-0.1393552,-0.9902424,0],[-0.1393552,-0.9902424,0],[-0.1393552,-0.9902424,0],[-0.1393552,-0.9902424,0],[-0.1393552,-0.9902424,0],[0.1393652,-0.9902411,0],[0.1393652,-0.9902411,0],[0.1393652,-0.9902411,0],[0.1393652,-0.9902411,0],[0.1393652,-0.9902411,0],[0.1393652,-0.9902411,0],[0,0,1],[0,0,1],[0,0,1],[-0,0,1],[-0,0,1],[-0,0,1],[-0,0,-1],[-0,0,-1],[-0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[-0,0,1],[-0,0,1],[-0,0,1],[0,0,1],[0,0,1],[0,0,1],[-0.1316406,-0.9912975,0],[-0.1316406,-0.9912975,0],[-0.1316406,-0.9912975,0],[-0.1316406,-0.9912975,0],[-0.1316406,-0.9912975,0],[-0.1316406,-0.9912975,0],[-0.6119635,-0.7908861,0],[-0.6119635,-0.7908861,0],[-0.6119635,-0.7908861,0],[-0.6119635,-0.7908861,0],[-0.6119635,-0.7908861,0],[-0.6119635,-0.7908861,0],[-0.1845051,-0.9828316,0],[-0.1845051,-0.9828316,0],[-0.1845051,-0.9828316,0],[-0.1845051,-0.9828316,0],[-0.1845051,-0.9828316,0],[-0.1845051,-0.9828316,0],[0.7936012,-0.6084384,0],[0.7936012,-0.6084384,0],[0.7936012,-0.6084384,0],[0.7936012,-0.6084384,0],[0.7936012,-0.6084384,0],[0.7936012,-0.6084384,0],[0.3680699,-0.9297982,0],[0.3680699,-0.9297982,0],[0.3680699,-0.9297982,0],[0.3680699,-0.9297982,0],[0.3680699,-0.9297982,0],[0.3680699,-0.9297982,0],[-0.42555,0.9049349,0],[-0.42555,0.9049349,0],[-0.42555,0.9049349,0],[-0.42555,0.9049349,0],[-0.42555,0.9049349,0],[-0.42555,0.9049349,0],[0.9810798,0.193604,0],[0.9810798,0.193604,0],[0.9810798,0.193604,0],[0.9810798,0.193604,-0],[0.9810798,0.193604,-0],[0.9810798,0.193604,-0],[0,0.2276432,0.9737446],[0,0.2276432,0.9737446],[0,0.2276432,0.9737446],[-0,0.2276432,0.9737446],[-0,0.2276432,0.9737446],[-0,0.2276432,0.9737446],[0,0,1],[0,0,1],[0,0,1],[-0,0,1],[-0,0,1],[-0,0,1],[0,0,1],[0,0,1],[0,0,1],[-0,0,1],[-0,0,1],[-0,0,1],[-0.9239847,-0.3824295,0],[-0.9239847,-0.3824295,0],[-0.9239847,-0.3824295,0],[-0.9239847,-0.3824295,-0],[-0.9239847,-0.3824295,-0],[-0.9239847,-0.3824295,-0],[-0.4103397,-0,-0.9119327],[-0.4103397,-0,-0.9119327],[-0.4103397,-0,-0.9119327],[-0.4103397,0,-0.9119327],[-0.4103397,0,-0.9119327],[-0.4103397,0,-0.9119327],[0.2130209,-0.9770476,0],[0.2130209,-0.9770476,0],[0.2130209,-0.9770476,0],[0.2130209,-0.9770476,0],[0.2130209,-0.9770476,0],[0.2130209,-0.9770476,0],[0.9215846,-0.1660794,-0.3508553],[0.9215846,-0.1660794,-0.3508553],[0.9215846,-0.1660794,-0.3508553],[0.9215846,-0.1660794,-0.3508553],[0.9215846,-0.1660794,-0.3508553],[0.9215846,-0.1660794,-0.3508553],[0.3795954,0.3797471,0.8436229],[0.3795954,0.3797471,0.8436229],[0.3795954,0.3797471,0.8436229],[0.1996382,0.1997179,0.9593004],[0.1996382,0.1997179,0.9593004],[0.1996382,0.1997179,0.9593004],[-0.9345582,0,-0.3558106],[-0.9345582,0,-0.3558106],[-0.9345582,0,-0.3558106],[-0.9345582,0,-0.3558106],[-0.9345582,0,-0.3558106],[-0.9345582,0,-0.3558106],[0,0,1],[0,0,1],[0,0,1],[-0,0,1],[-0,0,1],[-0,0,1],[1,-0,0],[1,-0,0],[1,-0,0],[1,0,0],[1,0,0],[1,0,0],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0.8988635,-0.4382288,0],[0.8988635,-0.4382288,0],[0.8988635,-0.4382288,0],[0.8988635,-0.4382288,0],[0.8988635,-0.4382288,0],[0.8988635,-0.4382288,0],[0.8948793,-0.4362864,0.09404916],[0.8948793,-0.4362864,0.09404916],[0.8948793,-0.4362864,0.09404916],[0.8948793,-0.4362864,0.09404916],[0.8948793,-0.4362864,0.09404916],[0.8948793,-0.4362864,0.09404916],[1,0,0],[1,0,0],[1,0,0],[1,0,-0],[1,0,-0],[1,0,-0],[1,0,-0],[1,0,-0],[1,0,-0],[1,0,0],[1,0,0],[1,0,0],[0.8814219,-0.4297254,-0.1960392],[0.8814219,-0.4297254,-0.1960392],[0.8814219,-0.4297254,-0.1960392],[0.8814219,-0.4297254,-0.1960392],[0.8814219,-0.4297254,-0.1960392],[0.8814219,-0.4297254,-0.1960392],[-0,0,-0.9999999],[-0,0,-0.9999999],[-0,0,-0.9999999],[0,0,-1],[0,0,-1],[0,0,-1],[0.8988641,-0.4382276,0],[0.8988641,-0.4382276,0],[0.8988641,-0.4382276,0],[0.8988641,-0.4382276,0],[0.8988641,-0.4382276,0],[0.8988641,-0.4382276,0],[-0,-0.9775435,-0.210733],[-0,-0.9775435,-0.210733],[-0,-0.9775435,-0.210733],[0,-0.9775435,-0.210733],[0,-0.9775435,-0.210733],[0,-0.9775435,-0.210733],[-0,0.977545,0.210727],[-0,0.977545,0.210727],[-0,0.977545,0.210727],[0,0.977545,0.210727],[0,0.977545,0.210727],[0,0.977545,0.210727],[0,-0.9097999,-0.4150473],[0,-0.9097999,-0.4150473],[0,-0.9097999,-0.4150473],[0,-0.9097999,-0.4150473],[0,-0.9097999,-0.4150473],[0,-0.9097999,-0.4150473],[-0.8814226,-0.4297242,-0.1960386],[-0.8814226,-0.4297242,-0.1960386],[-0.8814226,-0.4297242,-0.1960386],[-0.8814226,-0.4297242,-0.1960386],[-0.8814226,-0.4297242,-0.1960386],[-0.8814226,-0.4297242,-0.1960386],[-0.9999999,0,0],[-0.9999999,0,0],[-0.9999999,0,0],[-0.9999999,0,0],[-0.9999999,0,0],[-0.9999999,0,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[-0.9807911,-0.1950615,0],[-0.9807911,-0.1950615,0],[-0.9807911,-0.1950615,0],[-0.9807911,-0.1950615,0],[-0.9807911,-0.1950615,0],[-0.9807911,-0.1950615,0],[0,1,0],[0,1,0],[0,1,0],[0,1,-0],[0,1,-0],[0,1,-0],[0.9807909,-0.195062,0],[0.9807909,-0.195062,0],[0.9807909,-0.195062,0],[0.9807909,-0.195062,0],[0.9807909,-0.195062,0],[0.9807909,-0.195062,0],[0.2404871,0,-0.9706523],[0.2404871,0,-0.9706523],[0.2404871,0,-0.9706523],[0.2404871,0,-0.9706523],[0.2404871,0,-0.9706523],[0.2404871,0,-0.9706523],[0.3046618,-0.9524606,0],[0.3046618,-0.9524606,0],[0.3046618,-0.9524606,0],[0.3046618,-0.9524605,0],[0.3046618,-0.9524605,0],[0.3046618,-0.9524605,0],[-0,-0.9344294,-0.3561485],[-0,-0.9344294,-0.3561485],[-0,-0.9344294,-0.3561485],[0,-0.9344294,-0.3561485],[0,-0.9344294,-0.3561485],[0,-0.9344294,-0.3561485],[0,-0.4107704,-0.9117389],[0,-0.4107704,-0.9117389],[0,-0.4107704,-0.9117389],[0,-0.4107704,-0.9117389],[0,-0.4107704,-0.9117389],[0,-0.4107704,-0.9117389],[0.9239194,0.3825871,0],[0.9239194,0.3825871,0],[0.9239194,0.3825871,0],[0.9239194,0.3825871,0],[0.9239194,0.3825871,0],[0.9239194,0.3825871,0],[0.4077114,0.1167135,0.9056211],[0.4077114,0.1167135,0.9056211],[0.4077114,0.1167135,0.9056211],[0.4080538,0.1163152,0.905518],[0.4080538,0.1163152,0.905518],[0.4080538,0.1163152,0.905518],[0.3782091,0.4111477,-0.8294067],[0.3782091,0.4111477,-0.8294067],[0.3782091,0.4111477,-0.8294067],[0.3748817,0.4137425,-0.8296269],[0.3748817,0.4137425,-0.8296269],[0.3748817,0.4137425,-0.8296269],[-0.410517,0,0.911853],[-0.410517,0,0.911853],[-0.410517,0,0.911853],[-0.410517,0,0.911853],[-0.410517,0,0.911853],[-0.410517,0,0.911853],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,-0],[0,-1,-0],[0,-1,-0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[-0,0.5935985,0.8047613],[-0,0.5935985,0.8047613],[-0,0.5935985,0.8047613],[0,0.5935985,0.8047613],[0,0.5935985,0.8047613],[0,0.5935985,0.8047613],[-0.9800876,0.1985657,0],[-0.9800876,0.1985657,0],[-0.9800876,0.1985657,0],[-0.9800875,0.1985656,0],[-0.9800875,0.1985656,0],[-0.9800875,0.1985656,0],[-0.4074471,0.122066,-0.9050341],[-0.4074471,0.122066,-0.9050341],[-0.4074471,0.122066,-0.9050341],[-0.4083784,0.1234911,-0.9044208],[-0.4083784,0.1234911,-0.9044208],[-0.4083784,0.1234911,-0.9044208],[0,0.9999999,0],[0,0.9999999,0],[0,0.9999999,0],[0,1,0],[0,1,0],[0,1,0],[0.4074133,0.1227591,-0.9049555],[0.4074133,0.1227591,-0.9049555],[0.4074133,0.1227591,-0.9049555],[0.4074133,0.1227591,-0.9049555],[0.4074133,0.1227591,-0.9049555],[0.4074133,0.1227591,-0.9049555],[0.7618073,0,-0.6478037],[0.7618073,0,-0.6478037],[0.7618073,0,-0.6478037],[0.7618073,0,-0.6478037],[0.7618073,0,-0.6478037],[0.7618073,0,-0.6478037],[0.8662896,0.4669726,0.1774228],[0.8662896,0.4669726,0.1774228],[0.8662896,0.4669726,0.1774228],[0.868475,0.4641961,0.1739923],[0.868475,0.4641961,0.1739923],[0.868475,0.4641961,0.1739923],[0,0.4105328,-0.9118459],[0,0.4105328,-0.9118459],[0,0.4105328,-0.9118459],[0,0.4105328,-0.9118459],[0,0.4105328,-0.9118459],[0,0.4105328,-0.9118459],[-0.9999999,0,0],[-0.9999999,0,0],[-0.9999999,0,0],[-1,0,0],[-1,0,0],[-1,0,0],[0,1,0],[0,1,0],[0,1,0],[0,1,0],[0,1,0],[0,1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[-0.4048769,-0.9143712,0],[-0.4048769,-0.9143712,0],[-0.4048769,-0.9143712,0],[-0.4048769,-0.9143712,0],[-0.4048769,-0.9143712,0],[-0.4048769,-0.9143712,0],[0.386941,-0.9221044,0],[0.386941,-0.9221044,0],[0.386941,-0.9221044,0],[0.386941,-0.9221044,0],[0.386941,-0.9221044,0],[0.386941,-0.9221044,0],[0,0,1],[0,0,1],[0,0,1],[-0,0,1],[-0,0,1],[-0,0,1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-0.9999999],[0,0,-0.9999999],[0,0,-0.9999999],[0,0,-1],[0,0,-1],[0,0,-1],[-0,0,1],[-0,0,1],[-0,0,1],[0,0,1],[0,0,1],[0,0,1],[-0.3854305,-0.9227369,0],[-0.3854305,-0.9227369,0],[-0.3854305,-0.9227369,0],[-0.3854305,-0.9227369,0],[-0.3854305,-0.9227369,0],[-0.3854305,-0.9227369,0],[-0.3854305,-0.9227369,0],[-0.3854305,-0.9227369,0],[-0.3854305,-0.9227369,0],[-0.3854305,-0.9227369,0],[-0.3854305,-0.9227369,0],[-0.3854305,-0.9227369,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0.5032468,-0.8641427,0],[0.5032468,-0.8641427,0],[0.5032468,-0.8641427,0],[0.5032468,-0.8641427,0],[0.5032468,-0.8641427,0],[0.5032468,-0.8641427,0],[0.1316319,-0.9912987,0],[0.1316319,-0.9912987,0],[0.1316319,-0.9912987,0],[0.1316319,-0.9912987,0],[0.1316319,-0.9912987,0],[0.1316319,-0.9912987,0],[-0.42555,-0.9049349,-0],[-0.42555,-0.9049349,-0],[-0.42555,-0.9049349,-0],[-0.42555,-0.9049349,0],[-0.42555,-0.9049349,0],[-0.42555,-0.9049349,0],[0,-0.9170429,0.3987886],[0,-0.9170429,0.3987886],[0,-0.9170429,0.3987886],[0,-0.9170429,0.3987886],[0,-0.9170429,0.3987886],[0,-0.9170429,0.3987886],[0,-0.2276851,0.9737348],[0,-0.2276851,0.9737348],[0,-0.2276851,0.9737348],[0,-0.2276851,0.9737348],[0,-0.2276851,0.9737348],[0,-0.2276851,0.9737348],[-0,-0.3021924,-0.953247],[-0,-0.3021924,-0.953247],[-0,-0.3021924,-0.953247],[0,-0.3021924,-0.953247],[0,-0.3021924,-0.953247],[0,-0.3021924,-0.953247],[0.9345635,0,-0.3557964],[0.9345635,0,-0.3557964],[0.9345635,0,-0.3557964],[0.9345635,0,-0.3557964],[0.9345635,0,-0.3557964],[0.9345635,0,-0.3557964],[0.9841471,-0.1773539,0],[0.9841471,-0.1773539,0],[0.9841471,-0.1773539,0],[0.9841471,-0.1773539,0],[0.9841471,-0.1773539,0],[0.9841471,-0.1773539,0],[0.9345635,0,0.3557964],[0.9345635,0,0.3557964],[0.9345635,0,0.3557964],[0.9345635,-0,0.3557964],[0.9345635,-0,0.3557964],[0.9345635,-0,0.3557964],[-0.3608875,-0.8714774,-0.3320955],[-0.3608875,-0.8714774,-0.3320955],[-0.3608875,-0.8714774,-0.3320955],[-0.1847878,-0.9679556,-0.1700452],[-0.1847878,-0.9679556,-0.1700452],[-0.1847878,-0.9679556,-0.1700452],[0.8178453,-0.4839251,-0.3113608],[0.8178453,-0.4839251,-0.3113608],[0.8178453,-0.4839251,-0.3113608],[0.8178453,-0.4839251,-0.3113608],[0.8178453,-0.4839251,-0.3113608],[0.8178453,-0.4839251,-0.3113608],[0.8716339,0.3607456,0.3318386],[0.8716339,0.3607456,0.3318386],[0.8716339,0.3607456,0.3318386],[0.9676845,0.1855876,0.1707163],[0.9676845,0.1855876,0.1707163],[0.9676845,0.1855876,0.1707163],[-0.9345582,0,0.3558106],[-0.9345582,0,0.3558106],[-0.9345582,0,0.3558106],[-0.9345582,0,0.3558106],[-0.9345582,0,0.3558106],[-0.9345582,0,0.3558106],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[-0,0,1],[-0,0,1],[-0,0,1],[0,0,1],[0,0,1],[0,0,1],[0.8988634,-0.4382288,0],[0.8988634,-0.4382288,0],[0.8988634,-0.4382288,0],[0.8988634,-0.4382288,0],[0.8988634,-0.4382288,0],[0.8988634,-0.4382288,0],[0,-0.9775449,0.210727],[0,-0.9775449,0.210727],[0,-0.9775449,0.210727],[0,-0.9775449,0.210727],[0,-0.9775449,0.210727],[0,-0.9775449,0.210727],[0,0.9775437,-0.210733],[0,0.9775437,-0.210733],[0,0.9775437,-0.210733],[0,0.9775437,-0.210733],[0,0.9775437,-0.210733],[0,0.9775437,-0.210733],[0,-0.909798,0.4150513],[0,-0.909798,0.4150513],[0,-0.909798,0.4150513],[0,-0.909798,0.4150513],[0,-0.909798,0.4150513],[0,-0.909798,0.4150513],[-0.8814223,-0.4297241,0.1960408],[-0.8814223,-0.4297241,0.1960408],[-0.8814223,-0.4297241,0.1960408],[-0.8814223,-0.4297241,0.1960408],[-0.8814223,-0.4297241,0.1960408],[-0.8814223,-0.4297241,0.1960408],[0,1,0],[0,1,0],[0,1,0],[0,1,0],[0,1,0],[0,1,0],[-0.898864,-0.4382276,0],[-0.898864,-0.4382276,0],[-0.898864,-0.4382276,0],[-0.898864,-0.4382276,0],[-0.898864,-0.4382276,0],[-0.898864,-0.4382276,0],[-0.8948799,-0.4362851,-0.09404895],[-0.8948799,-0.4362851,-0.09404895],[-0.8948799,-0.4362851,-0.09404895],[-0.8948799,-0.4362851,-0.09404895],[-0.8948799,-0.4362851,-0.09404895],[-0.8948799,-0.4362851,-0.09404895],[-1,0,-0],[-1,0,-0],[-1,0,-0],[-1,0,0],[-1,0,0],[-1,0,0],[1,-0,0],[1,-0,0],[1,-0,0],[1,0,0],[1,0,0],[1,0,0],[-1,0,-0],[-1,0,-0],[-1,0,-0],[-1,0,0],[-1,0,0],[-1,0,0],[-1,0,0],[-1,0,0],[-1,0,0],[-1,0,0],[-1,0,0],[-1,0,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[-0.8314918,-0.5555371,0],[-0.8314918,-0.5555371,0],[-0.8314918,-0.5555371,0],[-0.831472,-0.5555666,2.407606e-05],[-0.831472,-0.5555666,2.407606e-05],[-0.831472,-0.5555666,2.407606e-05],[-0.8904347,0.455111,0],[-0.8904347,0.455111,0],[-0.8904347,0.455111,0],[-0.8904347,0.455111,0],[-0.8904347,0.455111,0],[-0.8904347,0.455111,0],[0.8314921,-0.5555366,0],[0.8314921,-0.5555366,0],[0.8314921,-0.5555366,0],[0.8314921,-0.5555366,0],[0.8314921,-0.5555366,0],[0.8314921,-0.5555366,0],[0.8662896,-0.4669726,-0.1774228],[0.8662896,-0.4669726,-0.1774228],[0.8662896,-0.4669726,-0.1774228],[0.8684785,-0.4641916,-0.1739866],[0.8684785,-0.4641916,-0.1739866],[0.8684785,-0.4641916,-0.1739866],[0.8372375,0.5468395,0],[0.8372375,0.5468395,0],[0.8372375,0.5468395,0],[0.8372375,0.5468395,-0],[0.8372375,0.5468395,-0],[0.8372375,0.5468395,-0],[0,-0.4105298,0.9118472],[0,-0.4105298,0.9118472],[0,-0.4105298,0.9118472],[0,-0.4105298,0.9118472],[0,-0.4105298,0.9118472],[0,-0.4105298,0.9118472],[0,0.4107764,0.9117361],[0,0.4107764,0.9117361],[0,0.4107764,0.9117361],[0,0.4107764,0.9117361],[0,0.4107764,0.9117361],[0,0.4107764,0.9117361],[0.3825871,0.9239194,0],[0.3825871,0.9239194,0],[0.3825871,0.9239194,0],[0.3825871,0.9239194,0],[0.3825871,0.9239194,0],[0.3825871,0.9239194,0],[0.9027781,0.2584338,0.3438078],[0.9027781,0.2584338,0.3438078],[0.9027781,0.2584338,0.3438078],[0.9025766,0.258817,0.3440482],[0.9025766,0.258817,0.3440482],[0.9025766,0.258817,0.3440482],[-0,0,-1],[-0,0,-1],[-0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[-0.9344156,0,0.3561847],[-0.9344156,0,0.3561847],[-0.9344156,0,0.3561847],[-0.9344156,0,0.3561847],[-0.9344156,0,0.3561847],[-0.9344156,0,0.3561847],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[0,0.5935985,-0.8047613],[0,0.5935985,-0.8047613],[0,0.5935985,-0.8047613],[0,0.5935985,-0.8047613],[0,0.5935985,-0.8047613],[0,0.5935985,-0.8047613],[-0.4074471,0.122066,0.9050341],[-0.4074471,0.122066,0.9050341],[-0.4074471,0.122066,0.9050341],[-0.4083783,0.1234911,0.9044208],[-0.4083783,0.1234911,0.9044208],[-0.4083783,0.1234911,0.9044208],[-0.218006,0.8479845,-0.4831104],[-0.218006,0.8479845,-0.4831104],[-0.218006,0.8479845,-0.4831104],[-0.2180028,0.8479828,-0.4831149],[-0.2180028,0.8479828,-0.4831149],[-0.2180028,0.8479828,-0.4831149],[0,-0,1],[0,-0,1],[0,-0,1],[0,0,1],[0,0,1],[0,0,1],[0,-0,1],[0,-0,1],[0,-0,1],[0,-0,0.9999999],[0,-0,0.9999999],[0,-0,0.9999999],[0,-0,1],[0,-0,1],[0,-0,1],[0,-0,0.9999999],[0,-0,0.9999999],[0,-0,0.9999999],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,0.9999999],[0,0,0.9999999],[0,0,0.9999999],[-0,0,1],[-0,0,1],[-0,0,1],[-0,0,0.9999999],[-0,0,0.9999999],[-0,0,0.9999999],[-0,-0,1],[-0,-0,1],[-0,-0,1],[0,-0,1],[0,-0,1],[0,-0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,-0,1],[0,-0,1],[0,-0,1],[0,-0,0.9999999],[0,-0,0.9999999],[0,-0,0.9999999],[0,0,1],[0,0,1],[0,0,1],[0,-0,0.9999999],[0,-0,0.9999999],[0,-0,0.9999999],[0,-0,1],[0,-0,1],[0,-0,1],[0,-0,1],[0,-0,1],[0,-0,1],[0,-0,1],[0,-0,1],[0,-0,1],[0.4074484,0.1220679,-0.9050333],[0.4074484,0.1220679,-0.9050333],[0.4074484,0.1220679,-0.9050333],[0.4083768,0.1234887,-0.9044219],[0.4083768,0.1234887,-0.9044219],[0.4083768,0.1234887,-0.9044219],[0.240487,0,-0.9706523],[0.240487,0,-0.9706523],[0.240487,0,-0.9706523],[0.2404871,0,-0.9706523],[0.2404871,0,-0.9706523],[0.2404871,0,-0.9706523],[0.3403893,0.8786296,0.3348809],[0.3403893,0.8786296,0.3348809],[0.3403893,0.8786296,0.3348809],[0.3440206,0.8780348,0.3327234],[0.3440206,0.8780348,0.3327234],[0.3440206,0.8780348,0.3327234],[0,0.9344294,-0.3561485],[0,0.9344294,-0.3561485],[0,0.9344294,-0.3561485],[0,0.9344294,-0.3561485],[0,0.9344294,-0.3561485],[0,0.9344294,-0.3561485],[-1,0,-0],[-1,0,-0],[-1,0,-0],[-1,0,0],[-1,0,0],[-1,0,0],[0.5180455,-0.855353,0],[0.5180455,-0.855353,0],[0.5180455,-0.855353,0],[0.5180455,-0.855353,0],[0.5180455,-0.855353,0],[0.5180455,-0.855353,0],[0.1393652,-0.9902411,0],[0.1393652,-0.9902411,0],[0.1393652,-0.9902411,0],[0.1393652,-0.9902411,0],[0.1393652,-0.9902411,0],[0.1393652,-0.9902411,0],[-0.1393552,-0.9902424,0],[-0.1393552,-0.9902424,0],[-0.1393552,-0.9902424,0],[-0.1393552,-0.9902424,0],[-0.1393552,-0.9902424,0],[-0.1393552,-0.9902424,0],[0,0,1],[0,0,1],[0,0,1],[-0,0,1],[-0,0,1],[-0,0,1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[-0,0,-1],[-0,0,-1],[-0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,1],[0,0,1],[0,0,1],[-0,0,0.9999999],[-0,0,0.9999999],[-0,0,0.9999999],[-0,0,1],[-0,0,1],[-0,0,1],[0,0,1],[0,0,1],[0,0,1],[-0.7958434,-0.6055024,0],[-0.7958434,-0.6055024,0],[-0.7958434,-0.6055024,0],[-0.7958434,-0.6055024,0],[-0.7958434,-0.6055024,0],[-0.7958434,-0.6055024,0],[-0.5334361,-0.8458403,0],[-0.5334361,-0.8458403,0],[-0.5334361,-0.8458403,0],[-0.5334361,-0.8458403,0],[-0.5334361,-0.8458403,0],[-0.5334361,-0.8458403,0],[0.9749084,-0.2226065,0],[0.9749084,-0.2226065,0],[0.9749084,-0.2226065,0],[0.9749084,-0.2226065,0],[0.9749084,-0.2226065,0],[0.9749084,-0.2226065,0],[0.6299935,-0.7766004,0],[0.6299935,-0.7766004,0],[0.6299935,-0.7766004,0],[0.6299935,-0.7766004,0],[0.6299935,-0.7766004,0],[0.6299935,-0.7766004,0],[-0.9810798,0.193604,0],[-0.9810798,0.193604,0],[-0.9810798,0.193604,0],[-0.9810798,0.193604,0],[-0.9810798,0.193604,0],[-0.9810798,0.193604,0],[0.42555,0.9049349,0],[0.42555,0.9049349,0],[0.42555,0.9049349,0],[0.42555,0.9049349,-0],[0.42555,0.9049349,-0],[0.42555,0.9049349,-0],[0,0.07545072,0.9971495],[0,0.07545072,0.9971495],[0,0.07545072,0.9971495],[-0,0.07545071,0.9971495],[-0,0.07545071,0.9971495],[-0,0.07545071,0.9971495],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[0,0,1],[0,0,1],[0,0,1],[-0,0,1],[-0,0,1],[-0,0,1],[-0.3826017,-0.9239134,0],[-0.3826017,-0.9239134,0],[-0.3826017,-0.9239134,0],[-0.3826017,-0.9239134,-0],[-0.3826017,-0.9239134,-0],[-0.3826017,-0.9239134,-0],[-0.9345581,-0,-0.3558105],[-0.9345581,-0,-0.3558105],[-0.9345581,-0,-0.3558105],[-0.9345581,0,-0.3558105],[-0.9345581,0,-0.3558105],[-0.9345581,0,-0.3558105],[0.6000378,-0.7999718,0],[0.6000378,-0.7999718,0],[0.6000378,-0.7999718,0],[0.6000378,-0.7999718,0],[0.6000378,-0.7999718,0],[0.6000378,-0.7999718,0],[0.409216,-0.07374512,-0.9094525],[0.409216,-0.07374512,-0.9094525],[0.409216,-0.07374512,-0.9094525],[0.4092226,-0.07374048,-0.9094499],[0.4092226,-0.07374048,-0.9094499],[0.4092226,-0.07374048,-0.9094499],[0.3608816,0.8714796,0.3320963],[0.3608816,0.8714796,0.3320963],[0.3608816,0.8714796,0.3320963],[0.1847838,0.9679565,0.1700447],[0.1847838,0.9679565,0.1700447],[0.1847838,0.9679565,0.1700447],[-0.4103397,0,-0.9119328],[-0.4103397,0,-0.9119328],[-0.4103397,0,-0.9119328],[-0.4103397,-0,-0.9119328],[-0.4103397,-0,-0.9119328],[-0.4103397,-0,-0.9119328],[0,0,1],[0,0,1],[0,0,1],[-0,0,1],[-0,0,1],[-0,0,1],[1,-0,0],[1,-0,0],[1,-0,0],[1,0,0],[1,0,0],[1,0,0],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0.9999999,-0,0],[0.9999999,-0,0],[0.9999999,-0,0],[0.9999999,0,0],[0.9999999,0,0],[0.9999999,0,0],[0.8948791,-0.4362863,-0.09405198],[0.8948791,-0.4362863,-0.09405198],[0.8948791,-0.4362863,-0.09405198],[0.8948791,-0.4362863,-0.09405198],[0.8948791,-0.4362863,-0.09405198],[0.8948791,-0.4362863,-0.09405198],[1,0,0],[1,0,0],[1,0,0],[1,-0,0],[1,-0,0],[1,-0,0],[1,-0,0],[1,-0,0],[1,-0,0],[1,0,0],[1,0,0],[1,0,0],[0.8814216,-0.4297253,0.1960413],[0.8814216,-0.4297253,0.1960413],[0.8814216,-0.4297253,0.1960413],[0.8814216,-0.4297253,0.1960413],[0.8814216,-0.4297253,0.1960413],[0.8814216,-0.4297253,0.1960413],[-0,0,0.9999999],[-0,0,0.9999999],[-0,0,0.9999999],[0,0,1],[0,0,1],[0,0,1],[0.8988641,-0.4382276,0],[0.8988641,-0.4382276,0],[0.8988641,-0.4382276,0],[0.8988641,-0.4382276,0],[0.8988641,-0.4382276,0],[0.8988641,-0.4382276,0],[0,-0.9775449,0.210727],[0,-0.9775449,0.210727],[0,-0.9775449,0.210727],[0,-0.9775449,0.210727],[0,-0.9775449,0.210727],[0,-0.9775449,0.210727],[0,0.9775437,-0.210733],[0,0.9775437,-0.210733],[0,0.9775437,-0.210733],[0,0.9775437,-0.210733],[0,0.9775437,-0.210733],[0,0.9775437,-0.210733],[0,-0.909798,0.4150513],[0,-0.909798,0.4150513],[0,-0.909798,0.4150513],[0,-0.909798,0.4150513],[0,-0.909798,0.4150513],[0,-0.909798,0.4150513],[-0.8814223,-0.4297241,0.1960408],[-0.8814223,-0.4297241,0.1960408],[-0.8814223,-0.4297241,0.1960408],[-0.8814223,-0.4297241,0.1960408],[-0.8814223,-0.4297241,0.1960408],[-0.8814223,-0.4297241,0.1960408],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[-1,0,0],[-1,0,0],[-1,0,0],[-1,0,0],[-1,0,0],[-1,0,0],[0.9999999,0,0],[0.9999999,0,0],[0.9999999,0,0],[0.9999999,0,0],[0.9999999,0,0],[0.9999999,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[0.7618072,0,-0.6478036],[0.7618072,0,-0.6478036],[0.7618072,0,-0.6478036],[0.7618073,0,-0.6478037],[0.7618073,0,-0.6478037],[0.7618073,0,-0.6478037],[0.837243,-0.5468311,0],[0.837243,-0.5468311,0],[0.837243,-0.5468311,0],[0.837243,-0.5468311,0],[0.837243,-0.5468311,0],[0.837243,-0.5468311,0],[-0,-0.4105298,-0.9118472],[-0,-0.4105298,-0.9118472],[-0,-0.4105298,-0.9118472],[0,-0.4105298,-0.9118472],[0,-0.4105298,-0.9118472],[0,-0.4105298,-0.9118472],[0,-0.9343773,-0.3562853],[0,-0.9343773,-0.3562853],[0,-0.9343773,-0.3562853],[0,-0.9343773,-0.3562853],[0,-0.9343773,-0.3562853],[0,-0.9343773,-0.3562853],[0,1,-0],[0,1,-0],[0,1,-0],[0,1,0],[0,1,0],[0,1,0],[0.4318366,0.8868036,0.1646111],[0.4318366,0.8868036,0.1646111],[0.4318366,0.8868036,0.1646111],[0.4322751,0.8873095,0.1606856],[0.4322751,0.8873095,0.1606856],[0.4322751,0.8873095,0.1606856],[0.6483614,0.7204047,-0.2462609],[0.6483614,0.7204047,-0.2462609],[0.6483614,0.7204047,-0.2462609],[0.6510112,0.7184958,-0.244843],[0.6510112,0.7184958,-0.244843],[0.6510112,0.7184958,-0.244843],[-0.9344155,0,0.3561847],[-0.9344155,0,0.3561847],[-0.9344155,0,0.3561847],[-0.9344155,0,0.3561847],[-0.9344155,0,0.3561847],[-0.9344155,0,0.3561847],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,1,0],[0,1,0],[0,1,0],[0,1,0],[0,1,0],[0,1,0],[0.5172729,0.8329642,0.1964674],[0.5172729,0.8329642,0.1964674],[0.5172729,0.8329642,0.1964674],[0.5171643,0.8330148,0.1965387],[0.5171643,0.8330148,0.1965387],[0.5171643,0.8330148,0.1965387],[0.3328105,0.5881914,-0.737067],[0.3328105,0.5881914,-0.737067],[0.3328105,0.5881914,-0.737067],[0.3325891,0.5883551,-0.7370365],[0.3325891,0.5883551,-0.7370365],[0.3325891,0.5883551,-0.7370365],[0,1,-0],[0,1,-0],[0,1,-0],[0,0.9999999,0],[0,0.9999999,0],[0,0.9999999,0],[0,0.2156895,-0.976462],[0,0.2156895,-0.976462],[0,0.2156895,-0.976462],[0,0.2156895,-0.9764621],[0,0.2156895,-0.9764621],[0,0.2156895,-0.9764621],[0.311262,-0.3901364,0.8665504],[0.311262,-0.3901364,0.8665504],[0.311262,-0.3901364,0.8665504],[0.3225103,-0.3830546,0.8655959],[0.3225103,-0.3830546,0.8655959],[0.3225103,-0.3830546,0.8655959],[0.3403893,-0.8786296,-0.3348809],[0.3403893,-0.8786296,-0.3348809],[0.3403893,-0.8786296,-0.3348809],[0.3440155,-0.8780355,-0.3327264],[0.3440155,-0.8780355,-0.3327264],[0.3440155,-0.8780355,-0.3327264],[0,0.4105328,0.9118459],[0,0.4105328,0.9118459],[0,0.4105328,0.9118459],[0,0.4105328,0.9118459],[0,0.4105328,0.9118459],[0,0.4105328,0.9118459],[1,-0,0],[1,-0,0],[1,-0,0],[0.9999999,0,0],[0.9999999,0,0],[0.9999999,0,0],[0,0,-1],[0,0,-1],[0,0,-1],[-0,-0,-1],[-0,-0,-1],[-0,-0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0.8782638,-0.4781764,0],[0.8782638,-0.4781764,0],[0.8782638,-0.4781764,0],[0.8782638,-0.4781764,0],[0.8782638,-0.4781764,0],[0.8782638,-0.4781764,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[-0.9324376,-0.3613309,0],[-0.9324376,-0.3613309,0],[-0.9324376,-0.3613309,0],[-0.9324376,-0.3613309,0],[-0.9324376,-0.3613309,0],[-0.9324376,-0.3613309,0],[-0.8266612,-0.5627,0],[-0.8266612,-0.5627,0],[-0.8266612,-0.5627,0],[-0.8266612,-0.5627,0],[-0.8266612,-0.5627,0],[-0.8266612,-0.5627,0],[0,0,-0.9999999],[0,0,-0.9999999],[0,0,-0.9999999],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-0.9999999],[0,0,-0.9999999],[0,0,-0.9999999],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,-1,-0],[0,-1,-0],[0,-1,-0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,-0],[0,-1,-0],[0,-1,-0],[-0.978511,-0.2061948,0],[-0.978511,-0.2061948,0],[-0.978511,-0.2061948,0],[-0.978511,-0.2061948,0],[-0.978511,-0.2061948,0],[-0.978511,-0.2061948,0],[0,0,1],[0,0,1],[0,0,1],[0,0,0.9999999],[0,0,0.9999999],[0,0,0.9999999],[-0,0,0.9999999],[-0,0,0.9999999],[-0,0,0.9999999],[0,0,1],[0,0,1],[0,0,1],[0,0,0.9999999],[0,0,0.9999999],[0,0,0.9999999],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[-0,0,0.9999999],[-0,0,0.9999999],[-0,0,0.9999999],[0,0,1],[0,0,1],[0,0,1],[-0,0,1],[-0,0,1],[-0,0,1],[-0,-0.9170429,-0.3987886],[-0,-0.9170429,-0.3987886],[-0,-0.9170429,-0.3987886],[0,-0.9170429,-0.3987886],[0,-0.9170429,-0.3987886],[0,-0.9170429,-0.3987886],[-0,-0.2276851,-0.9737348],[-0,-0.2276851,-0.9737348],[-0,-0.2276851,-0.9737348],[0,-0.2276851,-0.9737348],[0,-0.2276851,-0.9737348],[0,-0.2276851,-0.9737348],[-0.9796152,-0.2008836,-0],[-0.9796152,-0.2008836,-0],[-0.9796152,-0.2008836,-0],[-0.9796151,-0.2008836,0],[-0.9796151,-0.2008836,0],[-0.9796151,-0.2008836,0],[0,-0.3021924,0.953247],[0,-0.3021924,0.953247],[0,-0.3021924,0.953247],[0,-0.3021924,0.953247],[0,-0.3021924,0.953247],[0,-0.3021924,0.953247],[-1,0,0],[-1,0,0],[-1,0,0],[-1,-0,-0],[-1,-0,-0],[-1,-0,-0],[0,-0.9344506,-0.3560927],[0,-0.9344506,-0.3560927],[0,-0.9344506,-0.3560927],[0.0007155413,-0.9324825,-0.3612144],[0.0007155413,-0.9324825,-0.3612144],[0.0007155413,-0.9324825,-0.3612144],[0,-0.9344506,0.3560927],[0,-0.9344506,0.3560927],[0,-0.9344506,0.3560927],[0.0007155408,-0.9324825,0.3612144],[0.0007155408,-0.9324825,0.3612144],[0.0007155408,-0.9324825,0.3612144],[0.3795954,0.3797471,-0.8436229],[0.3795954,0.3797471,-0.8436229],[0.3795954,0.3797471,-0.8436229],[0.1996382,0.1997179,-0.9593004],[0.1996382,0.1997179,-0.9593004],[0.1996382,0.1997179,-0.9593004],[-0.3608875,-0.8714774,0.3320955],[-0.3608875,-0.8714774,0.3320955],[-0.3608875,-0.8714774,0.3320955],[-0.1847878,-0.9679556,0.1700452],[-0.1847878,-0.9679556,0.1700452],[-0.1847878,-0.9679556,0.1700452],[0.1000978,-0.5135276,0.8522147],[0.1000978,-0.5135276,0.8522147],[0.1000978,-0.5135276,0.8522147],[0.3230919,-0.4114231,0.8522574],[0.3230919,-0.4114231,0.8522574],[0.3230919,-0.4114231,0.8522574],[1,-0,0],[1,-0,0],[1,-0,0],[1,0,0],[1,0,0],[1,0,0],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,-1,-0],[0,-1,-0],[0,-1,-0],[0,-1,0],[0,-1,0],[0,-1,0],[-1,-0,0],[-1,-0,0],[-1,-0,0],[-1,0,0],[-1,0,0],[-1,0,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[-1,0,-0],[-1,0,-0],[-1,0,-0],[-1,0,0],[-1,0,0],[-1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,-0,0],[1,-0,0],[1,-0,0],[0.8814219,-0.4297254,0.1960392],[0.8814219,-0.4297254,0.1960392],[0.8814219,-0.4297254,0.1960392],[0.8814219,-0.4297254,0.1960392],[0.8814219,-0.4297254,0.1960392],[0.8814219,-0.4297254,0.1960392],[0,-0.9097999,0.4150473],[0,-0.9097999,0.4150473],[0,-0.9097999,0.4150473],[0,-0.9097999,0.4150473],[0,-0.9097999,0.4150473],[0,-0.9097999,0.4150473],[0.9999999,-0,0],[0.9999999,-0,0],[0.9999999,-0,0],[0.9999999,0,0],[0.9999999,0,0],[0.9999999,0,0],[-0.8948797,-0.436285,-0.09405173],[-0.8948797,-0.436285,-0.09405173],[-0.8948797,-0.436285,-0.09405173],[-0.8948797,-0.436285,-0.09405173],[-0.8948797,-0.436285,-0.09405173],[-0.8948797,-0.436285,-0.09405173],[0.8948799,-0.4362852,-0.09404894],[0.8948799,-0.4362852,-0.09404894],[0.8948799,-0.4362852,-0.09404894],[0.8948799,-0.4362852,-0.09404894],[0.8948799,-0.4362852,-0.09404894],[0.8948799,-0.4362852,-0.09404894],[-1,0,-0],[-1,0,-0],[-1,0,-0],[-1,0,0],[-1,0,0],[-1,0,0],[-0.8814222,-0.4297241,0.1960408],[-0.8814222,-0.4297241,0.1960408],[-0.8814222,-0.4297241,0.1960408],[-0.8814222,-0.4297241,0.1960408],[-0.8814222,-0.4297241,0.1960408],[-0.8814222,-0.4297241,0.1960408],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[-0,0,1],[-0,0,1],[-0,0,1],[0,0,1],[0,0,1],[0,0,1],[-1,0,0],[-1,0,0],[-1,0,0],[-1,0,0],[-1,0,0],[-1,0,0],[0.8314921,-0.5555367,0],[0.8314921,-0.5555367,0],[0.8314921,-0.5555367,0],[0.8314921,-0.5555367,0],[0.8314921,-0.5555367,0],[0.8314921,-0.5555367,0],[-0.8904347,0.455111,0],[-0.8904347,0.455111,0],[-0.8904347,0.455111,0],[-0.8904347,0.455111,0],[-0.8904347,0.455111,0],[-0.8904347,0.455111,0],[-0.831472,-0.5555666,-1.872248e-05],[-0.831472,-0.5555666,-1.872248e-05],[-0.831472,-0.5555666,-1.872248e-05],[-0.8314918,-0.5555371,-0],[-0.8314918,-0.5555371,-0],[-0.8314918,-0.5555371,-0],[0.8492427,0.2075781,0.4854875],[0.8492427,0.2075781,0.4854875],[0.8492427,0.2075781,0.4854875],[0.8422273,0.2188509,0.4927044],[0.8422273,0.2188509,0.4927044],[0.8422273,0.2188509,0.4927044],[0.8494051,-0.207478,0.4852462],[0.8494051,-0.207478,0.4852462],[0.8494051,-0.207478,0.4852462],[0.8422657,-0.218962,0.4925892],[0.8422657,-0.218962,0.4925892],[0.8422657,-0.218962,0.4925892],[-0,0,1],[-0,0,1],[-0,0,1],[0,0,1],[0,0,1],[0,0,1],[-1,0,0],[-1,0,0],[-1,0,0],[-1,0,-0],[-1,0,-0],[-1,0,-0],[0.6483614,0.7204047,0.2462609],[0.6483614,0.7204047,0.2462609],[0.6483614,0.7204047,0.2462609],[0.6510111,0.7184958,0.244843],[0.6510111,0.7184958,0.244843],[0.6510111,0.7184958,0.244843],[0.9344155,0,-0.3561847],[0.9344155,0,-0.3561847],[0.9344155,0,-0.3561847],[0.9344155,0,-0.3561847],[0.9344155,0,-0.3561847],[0.9344155,0,-0.3561847],[0.4379662,0.8989915,0],[0.4379662,0.8989915,0],[0.4379662,0.8989915,0],[0.4379662,0.8989915,0],[0.4379662,0.8989915,0],[0.4379662,0.8989915,0],[-1,-0,0],[-1,-0,0],[-1,-0,0],[-1,0,0],[-1,0,0],[-1,0,0],[-0.4105169,0,-0.911853],[-0.4105169,0,-0.911853],[-0.4105169,0,-0.911853],[-0.4105169,0,-0.911853],[-0.4105169,0,-0.911853],[-0.4105169,0,-0.911853],[0,-1,0],[0,-1,0],[0,-1,0],[-0,-1,0],[-0,-1,0],[-0,-1,0],[0.7713996,0.5647819,0.2931964],[0.7713996,0.5647819,0.2931964],[0.7713996,0.5647819,0.2931964],[0.7713572,0.5649562,0.2929721],[0.7713572,0.5649562,0.2929721],[0.7713572,0.5649562,0.2929721],[0.5172729,0.8329642,-0.1964674],[0.5172729,0.8329642,-0.1964674],[0.5172729,0.8329642,-0.1964674],[0.5171643,0.8330148,-0.1965387],[0.5171643,0.8330148,-0.1965387],[0.5171643,0.8330148,-0.1965387],[-0,1,0],[-0,1,0],[-0,1,0],[0,1,0],[0,1,0],[0,1,0],[0.9800876,0.1985657,0],[0.9800876,0.1985657,0],[0.9800876,0.1985657,0],[0.9800875,0.1985656,-0],[0.9800875,0.1985656,-0],[0.9800875,0.1985656,-0],[0.240487,0,0.9706523],[0.240487,0,0.9706523],[0.240487,0,0.9706523],[0.2404871,0,0.9706523],[0.2404871,0,0.9706523],[0.2404871,0,0.9706523],[0.3046618,-0.9524606,0],[0.3046618,-0.9524606,0],[0.3046618,-0.9524606,0],[0.3046618,-0.9524605,0],[0.3046618,-0.9524605,0],[0.3046618,-0.9524605,0],[0,-0.9344294,0.3561485],[0,-0.9344294,0.3561485],[0,-0.9344294,0.3561485],[0,-0.9344294,0.3561485],[0,-0.9344294,0.3561485],[0,-0.9344294,0.3561485],[0,-0.4107704,0.9117389],[0,-0.4107704,0.9117389],[0,-0.4107704,0.9117389],[0,-0.4107704,0.9117389],[0,-0.4107704,0.9117389],[0,-0.4107704,0.9117389],[0.9239194,0.3825871,-0],[0.9239194,0.3825871,-0],[0.9239194,0.3825871,-0],[0.9239194,0.3825871,0],[0.9239194,0.3825871,0],[0.9239194,0.3825871,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[-0.1916416,-0.981465,0],[-0.1916416,-0.981465,0],[-0.1916416,-0.981465,0],[-0.1916416,-0.981465,0],[-0.1916416,-0.981465,0],[-0.1916416,-0.981465,0],[-0.6342,-0.773169,0],[-0.6342,-0.773169,0],[-0.6342,-0.773169,0],[-0.6342,-0.773169,0],[-0.6342,-0.773169,0],[-0.6342,-0.773169,0],[0.6519735,-0.7582417,0],[0.6519735,-0.7582417,0],[0.6519735,-0.7582417,0],[0.6519735,-0.7582417,0],[0.6519735,-0.7582417,0],[0.6519735,-0.7582417,0],[0,0,1],[0,0,1],[0,0,1],[-0,0,1],[-0,0,1],[-0,0,1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,0.9999999],[0,0,0.9999999],[0,0,0.9999999],[0,0,1],[0,0,1],[0,0,1],[-0.6119635,-0.7908861,0],[-0.6119635,-0.7908861,0],[-0.6119635,-0.7908861,0],[-0.6119635,-0.7908861,0],[-0.6119635,-0.7908861,0],[-0.6119635,-0.7908861,0],[-0.1316406,-0.9912975,0],[-0.1316406,-0.9912975,0],[-0.1316406,-0.9912975,0],[-0.1316406,-0.9912975,0],[-0.1316406,-0.9912975,0],[-0.1316406,-0.9912975,0],[0.1708248,-0.9853013,0],[0.1708248,-0.9853013,0],[0.1708248,-0.9853013,0],[0.1708248,-0.9853013,0],[0.1708248,-0.9853013,0],[0.1708248,-0.9853013,0],[0.1708248,-0.9853013,0],[0.1708248,-0.9853013,0],[0.1708248,-0.9853013,0],[0.1708248,-0.9853013,0],[0.1708248,-0.9853013,0],[0.1708248,-0.9853013,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[-0.9810797,-0.1936043,-0],[-0.9810797,-0.1936043,-0],[-0.9810797,-0.1936043,-0],[-0.9810798,-0.1936043,0],[-0.9810798,-0.1936043,0],[-0.9810798,-0.1936043,0],[0,-0.1984801,0.980105],[0,-0.1984801,0.980105],[0,-0.1984801,0.980105],[0,-0.1984801,0.9801049],[0,-0.1984801,0.9801049],[0,-0.1984801,0.9801049],[0,-0.07545072,0.9971495],[0,-0.07545072,0.9971495],[0,-0.07545072,0.9971495],[0,-0.07545071,0.9971495],[0,-0.07545071,0.9971495],[0,-0.07545071,0.9971495],[0,0.302204,-0.9532433],[0,0.302204,-0.9532433],[0,0.302204,-0.9532433],[0,0.302204,-0.9532433],[0,0.302204,-0.9532433],[0,0.302204,-0.9532433],[0.4103333,0,-0.9119356],[0.4103333,0,-0.9119356],[0.4103333,0,-0.9119356],[0.4103333,0,-0.9119356],[0.4103333,0,-0.9119356],[0.4103333,0,-0.9119356],[0.8606254,-0.5092384,0],[0.8606254,-0.5092384,0],[0.8606254,-0.5092384,0],[0.8606254,-0.5092384,0],[0.8606254,-0.5092384,0],[0.8606254,-0.5092384,0],[0.4103333,0,0.9119356],[0.4103333,0,0.9119356],[0.4103333,0,0.9119356],[0.4103333,-0,0.9119356],[0.4103333,-0,0.9119356],[0.4103333,-0,0.9119356],[-0.8716245,-0.360758,-0.3318501],[-0.8716245,-0.360758,-0.3318501],[-0.8716245,-0.360758,-0.3318501],[-0.9676845,-0.1855876,-0.1707162],[-0.9676845,-0.1855876,-0.1707162],[-0.9676845,-0.1855876,-0.1707162],[0.3230919,-0.411423,-0.8522574],[0.3230919,-0.411423,-0.8522574],[0.3230919,-0.411423,-0.8522574],[0.1000978,-0.5135276,-0.8522146],[0.1000978,-0.5135276,-0.8522146],[0.1000978,-0.5135276,-0.8522146],[0.9215846,-0.1660794,0.3508553],[0.9215846,-0.1660794,0.3508553],[0.9215846,-0.1660794,0.3508553],[0.9215846,-0.1660794,0.3508553],[0.9215846,-0.1660794,0.3508553],[0.9215846,-0.1660794,0.3508553],[-0.4103397,0,0.9119328],[-0.4103397,0,0.9119328],[-0.4103397,0,0.9119328],[-0.4103397,0,0.9119328],[-0.4103397,0,0.9119328],[-0.4103397,0,0.9119328],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,1,-0],[0,1,-0],[0,1,-0],[0,1,0],[0,1,0],[0,1,0],[-0,0,-1],[-0,0,-1],[-0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0.8988635,-0.4382288,0],[0.8988635,-0.4382288,0],[0.8988635,-0.4382288,0],[0.8988635,-0.4382288,0],[0.8988635,-0.4382288,0],[0.8988635,-0.4382288,0],[-0,-0.9775435,-0.210733],[-0,-0.9775435,-0.210733],[-0,-0.9775435,-0.210733],[0,-0.9775435,-0.210733],[0,-0.9775435,-0.210733],[0,-0.9775435,-0.210733],[-0,0.9775449,0.210727],[-0,0.9775449,0.210727],[-0,0.9775449,0.210727],[0,0.9775449,0.210727],[0,0.9775449,0.210727],[0,0.9775449,0.210727],[0,-0.9097999,-0.4150473],[0,-0.9097999,-0.4150473],[0,-0.9097999,-0.4150473],[0,-0.9097999,-0.4150473],[0,-0.9097999,-0.4150473],[0,-0.9097999,-0.4150473],[-0.8814226,-0.4297242,-0.1960386],[-0.8814226,-0.4297242,-0.1960386],[-0.8814226,-0.4297242,-0.1960386],[-0.8814226,-0.4297242,-0.1960386],[-0.8814226,-0.4297242,-0.1960386],[-0.8814226,-0.4297242,-0.1960386],[0,0.9999999,0],[0,0.9999999,0],[0,0.9999999,0],[0,0.9999999,0],[0,0.9999999,0],[0,0.9999999,0],[-0.8988641,-0.4382276,0],[-0.8988641,-0.4382276,0],[-0.8988641,-0.4382276,0],[-0.8988641,-0.4382276,0],[-0.8988641,-0.4382276,0],[-0.8988641,-0.4382276,0],[-0.8948797,-0.436285,0.09405172],[-0.8948797,-0.436285,0.09405172],[-0.8948797,-0.436285,0.09405172],[-0.8948797,-0.436285,0.09405172],[-0.8948797,-0.436285,0.09405172],[-0.8948797,-0.436285,0.09405172],[-1,0,0],[-1,0,0],[-1,0,0],[-1,-0,0],[-1,-0,0],[-1,-0,0],[1,0,-0],[1,0,-0],[1,0,-0],[1,0,0],[1,0,0],[1,0,0],[-1,-0,0],[-1,-0,0],[-1,-0,0],[-1,0,0],[-1,0,0],[-1,0,0],[-1,0,0],[-1,0,0],[-1,0,0],[-1,0,0],[-1,0,0],[-1,0,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[-0.5554989,-0.8315173,1.607586e-05],[-0.5554989,-0.8315173,1.607586e-05],[-0.5554989,-0.8315173,1.607586e-05],[-0.5554989,-0.8315173,1.607586e-05],[-0.5554989,-0.8315173,1.607586e-05],[-0.5554989,-0.8315173,1.607586e-05],[0,0,-1],[0,0,-1],[0,0,-1],[-0,0,-1],[-0,0,-1],[-0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,-0,-1],[0,-0,-1],[0,-0,-1],[0,-0,-1],[0,-0,-1],[0,-0,-1],[0,-0,-1],[0,-0,-1],[0,-0,-1],[0,-0,-1],[0,-0,-1],[0,-0,-1],[0,0,-0.9999999],[0,0,-0.9999999],[0,0,-0.9999999],[0,-0,-1],[0,-0,-1],[0,-0,-1],[0,-0,-1],[0,-0,-1],[0,-0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,-0,-1],[0,-0,-1],[0,-0,-1],[0,-0,-1],[0,-0,-1],[0,-0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-0.9999999],[0,0,-0.9999999],[0,0,-0.9999999],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[-0,-0,-1],[-0,-0,-1],[-0,-0,-1],[-0,-0,-1],[-0,-0,-1],[-0,-0,-1],[-0,-0,-1],[-0,-0,-1],[-0,-0,-1],[-0,-0,-1],[-0,-0,-1],[-0,-0,-1],[0.5554983,-0.8315177,0],[0.5554983,-0.8315177,0],[0.5554983,-0.8315177,0],[0.5554983,-0.8315177,0],[0.5554983,-0.8315177,0],[0.5554983,-0.8315177,0],[0.3403893,-0.8786296,-0.3348809],[0.3403893,-0.8786296,-0.3348809],[0.3403893,-0.8786296,-0.3348809],[0.3440206,-0.8780348,-0.3327234],[0.3440206,-0.8780348,-0.3327234],[0.3440206,-0.8780348,-0.3327234],[0.3046618,0.9524606,0],[0.3046618,0.9524606,0],[0.3046618,0.9524606,0],[0.3046618,0.9524605,-0],[0.3046618,0.9524605,-0],[0.3046618,0.9524605,-0],[0,-0.9344294,0.3561485],[0,-0.9344294,0.3561485],[0,-0.9344294,0.3561485],[0,-0.9344294,0.3561485],[0,-0.9344294,0.3561485],[0,-0.9344294,0.3561485],[0,0.9343763,0.3562876],[0,0.9343763,0.3562876],[0,0.9343763,0.3562876],[0,0.9343763,0.3562876],[0,0.9343763,0.3562876],[0,0.9343763,0.3562876],[0,0,-1],[0,0,-1],[0,0,-1],[-0,-0,-1],[-0,-0,-1],[-0,-0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[-0,0,-1],[-0,0,-1],[-0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[-0.410517,0,-0.911853],[-0.410517,0,-0.911853],[-0.410517,0,-0.911853],[-0.410517,-0,-0.911853],[-0.410517,-0,-0.911853],[-0.410517,-0,-0.911853],[0,-1,-0],[0,-1,-0],[0,-1,-0],[0,-0.9999999,0],[0,-0.9999999,0],[0,-0.9999999,0],[0,0,-1],[0,0,-1],[0,0,-1],[-0,0,-1],[-0,0,-1],[-0,0,-1],[0,0.8511584,-0.524909],[0,0.8511584,-0.524909],[0,0.8511584,-0.524909],[0,0.8511584,-0.524909],[0,0.8511584,-0.524909],[0,0.8511584,-0.524909],[0.2180051,0.847984,-0.4831117],[0.2180051,0.847984,-0.4831117],[0.2180051,0.847984,-0.4831117],[0.2180068,0.847985,-0.4831094],[0.2180068,0.847985,-0.4831094],[0.2180068,0.847985,-0.4831094],[0.3327683,0.5883326,0.7369735],[0.3327683,0.5883326,0.7369735],[0.3327683,0.5883326,0.7369735],[0.3326469,0.5881618,0.7371646],[0.3326469,0.5881618,0.7371646],[0.3326469,0.5881618,0.7371646],[0,0.2156895,0.976462],[0,0.2156895,0.976462],[0,0.2156895,0.976462],[-0,0.2156895,0.9764621],[-0,0.2156895,0.9764621],[-0,0.2156895,0.9764621],[0.8492427,0.2075781,-0.4854875],[0.8492427,0.2075781,-0.4854875],[0.8492427,0.2075781,-0.4854875],[0.8422273,0.2188509,-0.4927044],[0.8422273,0.2188509,-0.4927044],[0.8422273,0.2188509,-0.4927044],[0.8494051,-0.207478,-0.4852462],[0.8494051,-0.207478,-0.4852462],[0.8494051,-0.207478,-0.4852462],[0.8422657,-0.218962,-0.4925892],[0.8422657,-0.218962,-0.4925892],[0.8422657,-0.218962,-0.4925892],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[-1,0,0],[-1,0,0],[-1,0,0],[-1,-0,0],[-1,-0,0],[-1,-0,0],[0.8049616,-0.5933269,0],[0.8049616,-0.5933269,0],[0.8049616,-0.5933269,0],[0.8049616,-0.5933269,0],[0.8049616,-0.5933269,0],[0.8049616,-0.5933269,0],[0.386941,-0.9221044,0],[0.386941,-0.9221044,0],[0.386941,-0.9221044,0],[0.386941,-0.9221044,0],[0.386941,-0.9221044,0],[0.386941,-0.9221044,0],[-0.4048769,-0.9143712,0],[-0.4048769,-0.9143712,0],[-0.4048769,-0.9143712,0],[-0.4048769,-0.9143712,0],[-0.4048769,-0.9143712,0],[-0.4048769,-0.9143712,0],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[-0,0,-1],[-0,0,-1],[-0,0,-1],[0,0,1],[0,0,1],[0,0,1],[-0,0,1],[-0,0,1],[-0,0,1],[0,0,1],[0,0,1],[0,0,1],[-0,0,1],[-0,0,1],[-0,0,1],[-0.925001,-0.3799649,0],[-0.925001,-0.3799649,0],[-0.925001,-0.3799649,0],[-0.925001,-0.3799649,0],[-0.925001,-0.3799649,0],[-0.925001,-0.3799649,0],[-0.8161331,-0.5778638,0],[-0.8161331,-0.5778638,0],[-0.8161331,-0.5778638,0],[-0.8161331,-0.5778638,0],[-0.8161331,-0.5778638,0],[-0.8161331,-0.5778638,0],[0.8701845,-0.492726,0],[0.8701845,-0.492726,0],[0.8701845,-0.492726,0],[0.8701845,-0.492726,0],[0.8701845,-0.492726,0],[0.8701845,-0.492726,0],[0.88901,-0.4578879,0],[0.88901,-0.4578879,0],[0.88901,-0.4578879,0],[0.88901,-0.4578879,0],[0.88901,-0.4578879,0],[0.88901,-0.4578879,0],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-0.9999999],[0,0,-0.9999999],[0,0,-0.9999999],[0,-0,-0.9999999],[0,-0,-0.9999999],[0,-0,-0.9999999],[0,0,-1],[0,0,-1],[0,0,-1],[-0,-0,-0.9999999],[-0,-0,-0.9999999],[-0,-0,-0.9999999],[0,0,-1],[0,0,-1],[0,0,-1],[-0,-0,-1],[-0,-0,-1],[-0,-0,-1],[0,-0,-0.9999999],[0,-0,-0.9999999],[0,-0,-0.9999999],[0,0,-1],[0,0,-1],[0,0,-1],[-0,-0,-1],[-0,-0,-1],[-0,-0,-1],[0,0.1984798,0.9801049],[0,0.1984798,0.9801049],[0,0.1984798,0.9801049],[-0,0.1984798,0.980105],[-0,0.1984798,0.980105],[-0,0.1984798,0.980105],[0,0.2276432,-0.9737446],[0,0.2276432,-0.9737446],[0,0.2276432,-0.9737446],[0,0.2276432,-0.9737446],[0,0.2276432,-0.9737446],[0,0.2276432,-0.9737446],[-1,0,0],[-1,0,0],[-1,0,0],[-1,0,0],[-1,0,0],[-1,0,0],[0,0.4104696,-0.9118742],[0,0.4104696,-0.9118742],[0,0.4104696,-0.9118742],[0,0.4104696,-0.9118742],[0,0.4104696,-0.9118742],[0,0.4104696,-0.9118742],[-0.9239847,0.3824295,0],[-0.9239847,0.3824295,0],[-0.9239847,0.3824295,0],[-0.9239847,0.3824295,0],[-0.9239847,0.3824295,0],[-0.9239847,0.3824295,0],[0.3468219,-0.3132194,0.8840861],[0.3468219,-0.3132194,0.8840861],[0.3468219,-0.3132194,0.8840861],[0.3857889,-0.3406662,0.8573875],[0.3857889,-0.3406662,0.8573875],[0.3857889,-0.3406662,0.8573875],[-0,0.4104696,0.9118742],[-0,0.4104696,0.9118742],[-0,0.4104696,0.9118742],[0,0.4104696,0.9118742],[0,0.4104696,0.9118742],[0,0.4104696,0.9118742],[0.3987478,-0.2359485,-0.8861876],[0.3987478,-0.2359485,-0.8861876],[0.3987478,-0.2359485,-0.8861876],[0.3987543,-0.2359458,-0.8861853],[0.3987543,-0.2359458,-0.8861853],[0.3987543,-0.2359458,-0.8861853],[-0.3608875,0.8714774,0.3320955],[-0.3608875,0.8714774,0.3320955],[-0.3608875,0.8714774,0.3320955],[-0.1847878,0.9679556,0.1700452],[-0.1847878,0.9679556,0.1700452],[-0.1847878,0.9679556,0.1700452],[-0.4103397,-0,-0.9119327],[-0.4103397,-0,-0.9119327],[-0.4103397,-0,-0.9119327],[-0.4103397,0,-0.9119327],[-0.4103397,0,-0.9119327],[-0.4103397,0,-0.9119327],[-0.9999999,0,0],[-0.9999999,0,0],[-0.9999999,0,0],[-0.9999999,0,0],[-0.9999999,0,0],[-0.9999999,0,0],[0,1,-0],[0,1,-0],[0,1,-0],[0,1,0],[0,1,0],[0,1,0],[1,-0,0],[1,-0,0],[1,-0,0],[1,0,0],[1,0,0],[1,0,0],[1,-0,0],[1,-0,0],[1,-0,0],[1,0,0],[1,0,0],[1,0,0],[-0.8948799,-0.4362852,0.09404894],[-0.8948799,-0.4362852,0.09404894],[-0.8948799,-0.4362852,0.09404894],[-0.8948799,-0.4362852,0.09404894],[-0.8948799,-0.4362852,0.09404894],[-0.8948799,-0.4362852,0.09404894],[0.8948791,-0.4362863,0.09405193],[0.8948791,-0.4362863,0.09405193],[0.8948791,-0.4362863,0.09405193],[0.8948791,-0.4362863,0.09405193],[0.8948791,-0.4362863,0.09405193],[0.8948791,-0.4362863,0.09405193],[-1,-0,0],[-1,-0,0],[-1,-0,0],[-1,0,0],[-1,0,0],[-1,0,0],[-0.8814223,-0.4297241,-0.1960408],[-0.8814223,-0.4297241,-0.1960408],[-0.8814223,-0.4297241,-0.1960408],[-0.8814223,-0.4297241,-0.1960408],[-0.8814223,-0.4297241,-0.1960408],[-0.8814223,-0.4297241,-0.1960408],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0.898864,-0.4382276,0],[0.898864,-0.4382276,0],[0.898864,-0.4382276,0],[0.898864,-0.4382276,0],[0.898864,-0.4382276,0],[0.898864,-0.4382276,0],[0.8948799,-0.4362851,0.09404895],[0.8948799,-0.4362851,0.09404895],[0.8948799,-0.4362851,0.09404895],[0.8948799,-0.4362851,0.09404895],[0.8948799,-0.4362851,0.09404895],[0.8948799,-0.4362851,0.09404895],[1,0,0],[1,0,0],[1,0,0],[1,0,-0],[1,0,-0],[1,0,-0],[1,0,-0],[1,0,-0],[1,0,-0],[1,0,0],[1,0,0],[1,0,0],[0.8814226,-0.4297242,-0.1960386],[0.8814226,-0.4297242,-0.1960386],[0.8814226,-0.4297242,-0.1960386],[0.8814226,-0.4297242,-0.1960386],[0.8814226,-0.4297242,-0.1960386],[0.8814226,-0.4297242,-0.1960386],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[-1,0,0],[-1,0,0],[-1,0,0],[-1,0,0],[-1,0,0],[-1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[0.3403893,0.8786296,0.3348809],[0.3403893,0.8786296,0.3348809],[0.3403893,0.8786296,0.3348809],[0.3440155,0.8780355,0.3327264],[0.3440155,0.8780355,0.3327264],[0.3440155,0.8780355,0.3327264],[0.3112625,0.3901391,-0.866549],[0.3112625,0.3901391,-0.866549],[0.3112625,0.3901391,-0.866549],[0.3229924,0.3827519,-0.8655502],[0.3229924,0.3827519,-0.8655502],[0.3229924,0.3827519,-0.8655502],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[-1,0,0],[-1,0,0],[-1,0,0],[-1,0,0],[-1,0,0],[-1,0,0],[0.3137324,0.6449382,0.6968692],[0.3137324,0.6449382,0.6968692],[0.3137324,0.6449382,0.6968692],[0.3175515,0.6521121,0.688412],[0.3175515,0.6521121,0.688412],[0.3175515,0.6521121,0.688412],[0.4080538,0.1163152,-0.905518],[0.4080538,0.1163152,-0.905518],[0.4080538,0.1163152,-0.905518],[0.4077113,0.1167135,-0.9056211],[0.4077113,0.1167135,-0.9056211],[0.4077113,0.1167135,-0.9056211],[-0.4105169,0,0.9118529],[-0.4105169,0,0.9118529],[-0.4105169,0,0.9118529],[-0.4105169,0,0.9118529],[-0.4105169,0,0.9118529],[-0.4105169,0,0.9118529],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,-0],[0,-1,-0],[0,-1,-0],[0,1,0],[0,1,0],[0,1,0],[0,1,-0],[0,1,-0],[0,1,-0],[0,-1,-0],[0,-1,-0],[0,-1,-0],[0,-1,0],[0,-1,0],[0,-1,0],[0.2180051,0.847984,0.4831117],[0.2180051,0.847984,0.4831117],[0.2180051,0.847984,0.4831117],[0.2180068,0.847985,0.4831093],[0.2180068,0.847985,0.4831093],[0.2180068,0.847985,0.4831093],[-0.3328119,0.5881917,0.7370664],[-0.3328119,0.5881917,0.7370664],[-0.3328119,0.5881917,0.7370664],[-0.3325911,0.5883548,0.7370358],[-0.3325911,0.5883548,0.7370358],[-0.3325911,0.5883548,0.7370358],[-0,0.9999999,0],[-0,0.9999999,0],[-0,0.9999999,0],[0,1,-0],[0,1,-0],[0,1,-0],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-0.9999999],[0,0,-0.9999999],[0,0,-0.9999999],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-0.9999999],[0,0,-0.9999999],[0,0,-0.9999999],[0,-0,-1],[0,-0,-1],[0,-0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,-0,-1],[0,-0,-1],[0,-0,-1],[0,0,-0.9999999],[0,0,-0.9999999],[0,0,-0.9999999],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-0.9999999],[0,0,-0.9999999],[0,0,-0.9999999],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-0.9999999],[0,0,-0.9999999],[0,0,-0.9999999],[-0,0,-1],[-0,0,-1],[-0,0,-1],[0,0,-0.9999999],[0,0,-0.9999999],[0,0,-0.9999999],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[-0,0,-1],[-0,0,-1],[-0,0,-1],[0.8494051,-0.207478,0.4852462],[0.8494051,-0.207478,0.4852462],[0.8494051,-0.207478,0.4852462],[0.8422657,-0.218962,0.4925892],[0.8422657,-0.218962,0.4925892],[0.8422657,-0.218962,0.4925892],[0.8662898,-0.4669719,-0.177424],[0.8662898,-0.4669719,-0.177424],[0.8662898,-0.4669719,-0.177424],[0.8684798,-0.4641891,-0.1739857],[0.8684798,-0.4641891,-0.1739857],[0.8684798,-0.4641891,-0.1739857],[0,0.9344294,0.3561485],[0,0.9344294,0.3561485],[0,0.9344294,0.3561485],[0,0.9344294,0.3561485],[0,0.9344294,0.3561485],[0,0.9344294,0.3561485],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,-0],[1,0,-0],[1,0,-0],[0.9767476,-0.2143927,0],[0.9767476,-0.2143927,0],[0.9767476,-0.2143927,0],[0.9767476,-0.2143927,0],[0.9767476,-0.2143927,0],[0.9767476,-0.2143927,0],[0.6519735,-0.7582417,0],[0.6519735,-0.7582417,0],[0.6519735,-0.7582417,0],[0.6519735,-0.7582417,0],[0.6519735,-0.7582417,0],[0.6519735,-0.7582417,0],[-0.6342,-0.773169,0],[-0.6342,-0.773169,0],[-0.6342,-0.773169,0],[-0.6342,-0.773169,0],[-0.6342,-0.773169,0],[-0.6342,-0.773169,0],[-0.1916416,-0.981465,0],[-0.1916416,-0.981465,0],[-0.1916416,-0.981465,0],[-0.1916416,-0.981465,0],[-0.1916416,-0.981465,0],[-0.1916416,-0.981465,0],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[-0,0,-1],[-0,0,-1],[-0,0,-1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,0.9999999],[0,0,0.9999999],[0,0,0.9999999],[-0.9923552,-0.1234143,0],[-0.9923552,-0.1234143,0],[-0.9923552,-0.1234143,0],[-0.9923552,-0.1234143,0],[-0.9923552,-0.1234143,0],[-0.9923552,-0.1234143,0],[-0.978511,-0.2061948,0],[-0.978511,-0.2061948,0],[-0.978511,-0.2061948,0],[-0.978511,-0.2061948,0],[-0.978511,-0.2061948,0],[-0.978511,-0.2061948,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-0,-1],[0,-0,-1],[0,-0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,-0,-1],[0,-0,-1],[0,-0,-1],[0,0,-0.9999999],[0,0,-0.9999999],[0,0,-0.9999999],[0,-0,-1],[0,-0,-1],[0,-0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,-0,-1],[0,-0,-1],[0,-0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[-0,-0,-0.9999999],[-0,-0,-0.9999999],[-0,-0,-0.9999999],[0,0.9170429,0.3987886],[0,0.9170429,0.3987886],[0,0.9170429,0.3987886],[-0,0.9170429,0.3987886],[-0,0.9170429,0.3987886],[-0,0.9170429,0.3987886],[0,0.07545072,-0.9971495],[0,0.07545072,-0.9971495],[0,0.07545072,-0.9971495],[0,0.07545071,-0.9971495],[0,0.07545071,-0.9971495],[0,0.07545071,-0.9971495],[0,0,-1],[0,0,-1],[0,0,-1],[-0,0,-1],[-0,0,-1],[-0,0,-1],[0,0.9344506,-0.3560927],[0,0.9344506,-0.3560927],[0,0.9344506,-0.3560927],[0,0.9344506,-0.3560927],[0,0.9344506,-0.3560927],[0,0.9344506,-0.3560927],[-0.3826017,0.9239134,0],[-0.3826017,0.9239134,0],[-0.3826017,0.9239134,0],[-0.3826017,0.9239134,0],[-0.3826017,0.9239134,0],[-0.3826017,0.9239134,0],[0.6994621,-0.6176515,0.3595266],[0.6994621,-0.6176515,0.3595266],[0.6994621,-0.6176515,0.3595266],[0.727269,-0.6280276,0.2768776],[0.727269,-0.6280276,0.2768776],[0.727269,-0.6280276,0.2768776],[-0,0.9344506,0.3560927],[-0,0.9344506,0.3560927],[-0,0.9344506,0.3560927],[0,0.9344506,0.3560927],[0,0.9344506,0.3560927],[0,0.9344506,0.3560927],[0.8716339,0.3607456,-0.3318386],[0.8716339,0.3607456,-0.3318386],[0.8716339,0.3607456,-0.3318386],[0.9676845,0.1855876,-0.1707162],[0.9676845,0.1855876,-0.1707162],[0.9676845,0.1855876,-0.1707162],[-0.3796015,0.379746,0.8436206],[-0.3796015,0.379746,0.8436206],[-0.3796015,0.379746,0.8436206],[-0.1996425,0.1997184,0.9592994],[-0.1996425,0.1997184,0.9592994],[-0.1996425,0.1997184,0.9592994],[-0.9345581,0,-0.3558105],[-0.9345581,0,-0.3558105],[-0.9345581,0,-0.3558105],[-0.9345581,0,-0.3558105],[-0.9345581,0,-0.3558105],[-0.9345581,0,-0.3558105],[-1,0,0],[-1,0,0],[-1,0,0],[-1,0,0],[-1,0,0],[-1,0,0],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[0.9999999,-0,0],[0.9999999,-0,0],[0.9999999,-0,0],[0.9999999,0,0],[0.9999999,0,0],[0.9999999,0,0],[-0.8948797,-0.436285,-0.09405173],[-0.8948797,-0.436285,-0.09405173],[-0.8948797,-0.436285,-0.09405173],[-0.8948797,-0.436285,-0.09405173],[-0.8948797,-0.436285,-0.09405173],[-0.8948797,-0.436285,-0.09405173],[0.8948793,-0.4362864,-0.09404921],[0.8948793,-0.4362864,-0.09404921],[0.8948793,-0.4362864,-0.09404921],[0.8948793,-0.4362864,-0.09404921],[0.8948793,-0.4362864,-0.09404921],[0.8948793,-0.4362864,-0.09404921],[-1,0,-0],[-1,0,-0],[-1,0,-0],[-1,0,0],[-1,0,0],[-1,0,0],[-0.8814226,-0.4297242,0.1960386],[-0.8814226,-0.4297242,0.1960386],[-0.8814226,-0.4297242,0.1960386],[-0.8814226,-0.4297242,0.1960386],[-0.8814226,-0.4297242,0.1960386],[-0.8814226,-0.4297242,0.1960386],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0.9999999,-0,0],[0.9999999,-0,0],[0.9999999,-0,0],[0.9999999,0,0],[0.9999999,0,0],[0.9999999,0,0],[0.8948797,-0.436285,-0.09405172],[0.8948797,-0.436285,-0.09405172],[0.8948797,-0.436285,-0.09405172],[0.8948797,-0.436285,-0.09405172],[0.8948797,-0.436285,-0.09405172],[0.8948797,-0.436285,-0.09405172],[1,0,0],[1,0,0],[1,0,0],[1,-0,0],[1,-0,0],[1,-0,0],[1,-0,0],[1,-0,0],[1,-0,0],[1,0,0],[1,0,0],[1,0,0],[0.8814223,-0.4297241,0.1960408],[0.8814223,-0.4297241,0.1960408],[0.8814223,-0.4297241,0.1960408],[0.8814223,-0.4297241,0.1960408],[0.8814223,-0.4297241,0.1960408],[0.8814223,-0.4297241,0.1960408],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,-0],[0,-1,-0],[0,-1,-0],[0,1,0],[0,1,0],[0,1,0],[0,1,0],[0,1,0],[0,1,0],[-1,0,0],[-1,0,0],[-1,0,0],[-1,0,0],[-1,0,0],[-1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[-1,0,0],[-1,0,0],[-1,0,0],[-1,0,0],[-1,0,0],[-1,0,0],[0.8662898,0.4669719,0.177424],[0.8662898,0.4669719,0.177424],[0.8662898,0.4669719,0.177424],[0.8684763,0.4641937,0.1739914],[0.8684763,0.4641937,0.1739914],[0.8684763,0.4641937,0.1739914],[0.8492427,0.2075781,-0.4854875],[0.8492427,0.2075781,-0.4854875],[0.8492427,0.2075781,-0.4854875],[0.8422273,0.2188509,-0.4927044],[0.8422273,0.2188509,-0.4927044],[0.8422273,0.2188509,-0.4927044],[0,0,-1],[0,0,-1],[0,0,-1],[-0,0,-1],[-0,0,-1],[-0,0,-1],[1,-0,0],[1,-0,0],[1,-0,0],[1,0,0],[1,0,0],[1,0,0],[0.9344155,-0,0.3561847],[0.9344155,-0,0.3561847],[0.9344155,-0,0.3561847],[0.9344155,0,0.3561847],[0.9344155,0,0.3561847],[0.9344155,0,0.3561847],[0.9027781,0.2584338,-0.3438077],[0.9027781,0.2584338,-0.3438077],[0.9027781,0.2584338,-0.3438077],[0.9025766,0.258817,-0.3440482],[0.9025766,0.258817,-0.3440482],[0.9025766,0.258817,-0.3440482],[-0.4105169,-0,-0.9118529],[-0.4105169,-0,-0.9118529],[-0.4105169,-0,-0.9118529],[-0.4105169,0,-0.9118529],[-0.4105169,0,-0.9118529],[-0.4105169,0,-0.9118529],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-0.9999999],[0,0,-0.9999999],[0,0,-0.9999999],[0,0,-1],[0,0,-1],[0,0,-1],[-0,-0,-1],[-0,-0,-1],[-0,-0,-1],[0,0,-0.9999999],[0,0,-0.9999999],[0,0,-0.9999999],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-0.9999999],[0,0,-0.9999999],[0,0,-0.9999999],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-0.9999999],[0,0,-0.9999999],[0,0,-0.9999999],[0,0,-0.9999999],[0,0,-0.9999999],[0,0,-0.9999999],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-0.9999999],[0,0,-0.9999999],[0,0,-0.9999999],[-0,0,-1],[-0,0,-1],[-0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-0.9999999],[0,0,-0.9999999],[0,0,-0.9999999],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-0.9999999],[0,0,-0.9999999],[0,0,-0.9999999],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[-0,0,-1],[-0,0,-1],[-0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[-0,-1,0],[-0,-1,0],[-0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[-0.7712991,0.5649387,0.2931583],[-0.7712991,0.5649387,0.2931583],[-0.7712991,0.5649387,0.2931583],[-0.7714956,0.5647374,0.2930295],[-0.7714956,0.5647374,0.2930295],[-0.7714956,0.5647374,0.2930295],[-0.9280505,0.116522,-0.3537584],[-0.9280505,0.116522,-0.3537584],[-0.9280505,0.116522,-0.3537584],[-0.9285353,0.1151894,-0.3529212],[-0.9285353,0.1151894,-0.3529212],[-0.9285353,0.1151894,-0.3529212],[0,0.9999999,0],[0,0.9999999,0],[0,0.9999999,0],[-0,1,0],[-0,1,0],[-0,1,0],[0.9280494,0.1165211,-0.3537614],[0.9280494,0.1165211,-0.3537614],[0.9280494,0.1165211,-0.3537614],[0.9285359,0.1151842,-0.3529214],[0.9285359,0.1151842,-0.3529214],[0.9285359,0.1151842,-0.3529214],[0.3403893,-0.8786296,0.3348809],[0.3403893,-0.8786296,0.3348809],[0.3403893,-0.8786296,0.3348809],[0.3440206,-0.8780348,0.3327233],[0.3440206,-0.8780348,0.3327233],[0.3440206,-0.8780348,0.3327233],[0.3046618,0.9524606,0],[0.3046618,0.9524606,0],[0.3046618,0.9524606,0],[0.3046618,0.9524605,0],[0.3046618,0.9524605,0],[0.3046618,0.9524605,0],[0,-0.9344294,-0.3561485],[0,-0.9344294,-0.3561485],[0,-0.9344294,-0.3561485],[0,-0.9344294,-0.3561485],[0,-0.9344294,-0.3561485],[0,-0.9344294,-0.3561485],[0,0.9343763,-0.3562876],[0,0.9343763,-0.3562876],[0,0.9343763,-0.3562876],[0,0.9343763,-0.3562876],[0,0.9343763,-0.3562876],[0,0.9343763,-0.3562876],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[-0,0,1],[-0,0,1],[-0,0,1],[0,0,1],[0,0,1],[0,0,1],[0.8782638,-0.4781764,0],[0.8782638,-0.4781764,0],[0.8782638,-0.4781764,0],[0.8782638,-0.4781764,0],[0.8782638,-0.4781764,0],[0.8782638,-0.4781764,0],[0.8994453,-0.4370334,0],[0.8994453,-0.4370334,0],[0.8994453,-0.4370334,0],[0.8994453,-0.4370334,0],[0.8994453,-0.4370334,0],[0.8994453,-0.4370334,0],[-0.8123698,-0.5831425,0],[-0.8123698,-0.5831425,0],[-0.8123698,-0.5831425,0],[-0.8123698,-0.5831425,0],[-0.8123698,-0.5831425,0],[-0.8123698,-0.5831425,0],[-0.5483773,-0.836231,0],[-0.5483773,-0.836231,0],[-0.5483773,-0.836231,0],[-0.5483773,-0.836231,0],[-0.5483773,-0.836231,0],[-0.5483773,-0.836231,0],[0,0,-0.9999999],[0,0,-0.9999999],[0,0,-0.9999999],[0,0,-0.9999999],[0,0,-0.9999999],[0,0,-0.9999999],[0,0,-1],[0,0,-1],[0,0,-1],[-0,0,-1],[-0,0,-1],[-0,0,-1],[0,0,0.9999999],[0,0,0.9999999],[0,0,0.9999999],[0,0,0.9999999],[0,0,0.9999999],[0,0,0.9999999],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[-0.8555613,-0.5177015,0],[-0.8555613,-0.5177015,0],[-0.8555613,-0.5177015,0],[-0.8555613,-0.5177015,0],[-0.8555613,-0.5177015,0],[-0.8555613,-0.5177015,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,-0],[0,-1,-0],[0,-1,-0],[-0.8555613,-0.5177015,0],[-0.8555613,-0.5177015,0],[-0.8555613,-0.5177015,0],[-0.8555613,-0.5177015,0],[-0.8555613,-0.5177015,0],[-0.8555613,-0.5177015,0],[0,0,1],[0,0,1],[0,0,1],[-0,0,1],[-0,0,1],[-0,0,1],[0,0,1],[0,0,1],[0,0,1],[-0,0,1],[-0,0,1],[-0,0,1],[0,0,1],[0,0,1],[0,0,1],[-0,0,1],[-0,0,1],[-0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[-0,-0.1984801,-0.980105],[-0,-0.1984801,-0.980105],[-0,-0.1984801,-0.980105],[0,-0.1984801,-0.9801049],[0,-0.1984801,-0.9801049],[0,-0.1984801,-0.9801049],[-0,-0.07545072,-0.9971495],[-0,-0.07545072,-0.9971495],[-0,-0.07545072,-0.9971495],[0,-0.07545071,-0.9971495],[0,-0.07545071,-0.9971495],[0,-0.07545071,-0.9971495],[-0.9976759,-0.06813686,-0],[-0.9976759,-0.06813686,-0],[-0.9976759,-0.06813686,-0],[-0.9976759,-0.06813686,0],[-0.9976759,-0.06813686,0],[-0.9976759,-0.06813686,0],[0,0.302204,0.9532433],[0,0.302204,0.9532433],[0,0.302204,0.9532433],[0,0.302204,0.9532433],[0,0.302204,0.9532433],[0,0.302204,0.9532433],[-1,0,0],[-1,0,0],[-1,0,0],[-1,-0,-0],[-1,-0,-0],[-1,-0,-0],[0.0003149735,-0.4104696,-0.9118742],[0.0003149735,-0.4104696,-0.9118742],[0.0003149735,-0.4104696,-0.9118742],[-0.006032857,-0.4584159,-0.8887173],[-0.006032857,-0.4584159,-0.8887173],[-0.006032857,-0.4584159,-0.8887173],[0.0003149741,-0.4104696,0.9118742],[0.0003149741,-0.4104696,0.9118742],[0.0003149741,-0.4104696,0.9118742],[-0.006032852,-0.4584158,0.8887173],[-0.006032852,-0.4584158,0.8887173],[-0.006032852,-0.4584158,0.8887173],[0.3608816,0.8714796,-0.3320963],[0.3608816,0.8714796,-0.3320963],[0.3608816,0.8714796,-0.3320963],[0.1847838,0.9679565,-0.1700447],[0.1847838,0.9679565,-0.1700447],[0.1847838,0.9679565,-0.1700447],[-0.8716245,0.360758,0.3318501],[-0.8716245,0.360758,0.3318501],[-0.8716245,0.360758,0.3318501],[-0.9676845,0.1855876,0.1707162],[-0.9676845,0.1855876,0.1707162],[-0.9676845,0.1855876,0.1707162],[0.2019745,-0.9263823,0.3178397],[0.2019745,-0.9263823,0.3178397],[0.2019745,-0.9263823,0.3178397],[0.01942625,-0.9889103,0.1472382],[0.01942625,-0.9889103,0.1472382],[0.01942625,-0.9889103,0.1472382],[0,1,0],[0,1,0],[0,1,0],[0,1,0],[0,1,0],[0,1,0],[0,-1,-0],[0,-1,-0],[0,-1,-0],[0,-1,0],[0,-1,0],[0,-1,0],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[-0.9999999,-0,0],[-0.9999999,-0,0],[-0.9999999,-0,0],[-0.9999999,0,0],[-0.9999999,0,0],[-0.9999999,0,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[-1,0,0],[-1,0,0],[-1,0,0],[-1,-0,0],[-1,-0,0],[-1,-0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,-0],[1,0,-0],[1,0,-0],[0.8814216,-0.4297253,-0.1960413],[0.8814216,-0.4297253,-0.1960413],[0.8814216,-0.4297253,-0.1960413],[0.8814216,-0.4297253,-0.1960413],[0.8814216,-0.4297253,-0.1960413],[0.8814216,-0.4297253,-0.1960413],[0,-0.909798,-0.4150513],[0,-0.909798,-0.4150513],[0,-0.909798,-0.4150513],[0,-0.909798,-0.4150513],[0,-0.909798,-0.4150513],[0,-0.909798,-0.4150513],[1,-0,0],[1,-0,0],[1,-0,0],[1,0,0],[1,0,0],[1,0,0],[-0.8948799,-0.4362852,0.09404894],[-0.8948799,-0.4362852,0.09404894],[-0.8948799,-0.4362852,0.09404894],[-0.8948799,-0.4362852,0.09404894],[-0.8948799,-0.4362852,0.09404894],[-0.8948799,-0.4362852,0.09404894],[0.8948797,-0.436285,0.09405173],[0.8948797,-0.436285,0.09405173],[0.8948797,-0.436285,0.09405173],[0.8948797,-0.436285,0.09405173],[0.8948797,-0.436285,0.09405173],[0.8948797,-0.436285,0.09405173],[-1,-0,0],[-1,-0,0],[-1,-0,0],[-1,0,0],[-1,0,0],[-1,0,0],[-0.8814226,-0.4297242,-0.1960386],[-0.8814226,-0.4297242,-0.1960386],[-0.8814226,-0.4297242,-0.1960386],[-0.8814226,-0.4297242,-0.1960386],[-0.8814226,-0.4297242,-0.1960386],[-0.8814226,-0.4297242,-0.1960386],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,1,0],[0,1,0],[0,1,0],[0,1,0],[0,1,0],[0,1,0],[-0.9999999,0,0],[-0.9999999,0,0],[-0.9999999,0,0],[-0.9999999,0,0],[-0.9999999,0,0],[-0.9999999,0,0],[0.9807909,-0.1950621,0],[0.9807909,-0.1950621,0],[0.9807909,-0.1950621,0],[0.9807909,-0.1950621,0],[0.9807909,-0.1950621,0],[0.9807909,-0.1950621,0],[0,1,0],[0,1,0],[0,1,0],[0,1,0],[0,1,0],[0,1,0],[-0.980791,-0.1950614,0],[-0.980791,-0.1950614,0],[-0.980791,-0.1950614,0],[-0.980791,-0.1950614,-0],[-0.980791,-0.1950614,-0],[-0.980791,-0.1950614,-0],[0.3112625,0.3901391,0.866549],[0.3112625,0.3901391,0.866549],[0.3112625,0.3901391,0.866549],[0.3229925,0.3827518,0.86555],[0.3229925,0.3827518,0.86555],[0.3229925,0.3827518,0.86555],[0.311262,-0.3901364,0.8665504],[0.311262,-0.3901364,0.8665504],[0.311262,-0.3901364,0.8665504],[0.3225103,-0.3830546,0.8655959],[0.3225103,-0.3830546,0.8655959],[0.3225103,-0.3830546,0.8655959],[0,1,0],[0,1,0],[0,1,0],[0,1,0],[0,1,0],[0,1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,-0],[0,-1,-0],[0,-1,-0],[0.4105169,-0,0.911853],[0.4105169,-0,0.911853],[0.4105169,-0,0.911853],[0.4105169,0,0.911853],[0.4105169,0,0.911853],[0.4105169,0,0.911853],[0.4105169,0,-0.911853],[0.4105169,0,-0.911853],[0.4105169,0,-0.911853],[0.4105169,0,-0.911853],[0.4105169,0,-0.911853],[0.4105169,0,-0.911853],[-0.9344155,-0,-0.3561847],[-0.9344155,-0,-0.3561847],[-0.9344155,-0,-0.3561847],[-0.9344155,0,-0.3561847],[-0.9344155,0,-0.3561847],[-0.9344155,0,-0.3561847],[0,0,1],[0,0,1],[0,0,1],[0,-0,0.9999999],[0,-0,0.9999999],[0,-0,0.9999999],[0,-0,1],[0,-0,1],[0,-0,1],[0,0,1],[0,0,1],[0,0,1],[0,-0,0.9999999],[0,-0,0.9999999],[0,-0,0.9999999],[-0,0,1],[-0,0,1],[-0,0,1],[0,-0,1],[0,-0,1],[0,-0,1],[0,-0,1],[0,-0,1],[0,-0,1],[0,-0,0.9999999],[0,-0,0.9999999],[0,-0,0.9999999],[0,-0,1],[0,-0,1],[0,-0,1],[0,-0,1],[0,-0,1],[0,-0,1],[0,-0,0.9999999],[0,-0,0.9999999],[0,-0,0.9999999],[0,-0,0.9999999],[0,-0,0.9999999],[0,-0,0.9999999],[0,-0,1],[0,-0,1],[0,-0,1],[0,-0,1],[0,-0,1],[0,-0,1],[0,-0,1],[0,-0,1],[0,-0,1],[-0,-0,0.9999999],[-0,-0,0.9999999],[-0,-0,0.9999999],[0,-0,1],[0,-0,1],[0,-0,1],[-0,-0,1],[-0,-0,1],[-0,-0,1],[0,0,0.9999999],[0,0,0.9999999],[0,0,0.9999999],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,0.9999999],[0,0,0.9999999],[0,0,0.9999999],[-0,0,1],[-0,0,1],[-0,0,1],[0,-0,1],[0,-0,1],[0,-0,1],[0,0,1],[0,0,1],[0,0,1],[0,-0,1],[0,-0,1],[0,-0,1],[0,-0,1],[0,-0,1],[0,-0,1],[0,0,1],[0,0,1],[0,0,1],[-0.9344155,0,-0.3561847],[-0.9344155,0,-0.3561847],[-0.9344155,0,-0.3561847],[-0.9344155,0,-0.3561847],[-0.9344155,0,-0.3561847],[-0.9344155,0,-0.3561847],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,-0],[0,-1,-0],[0,-1,-0],[-0.3328118,0.5881917,-0.7370664],[-0.3328118,0.5881917,-0.7370664],[-0.3328118,0.5881917,-0.7370664],[-0.332591,0.5883548,-0.7370358],[-0.332591,0.5883548,-0.7370358],[-0.332591,0.5883548,-0.7370358],[0.7713996,0.5647818,-0.2931964],[0.7713996,0.5647818,-0.2931964],[0.7713996,0.5647818,-0.2931964],[0.7713572,0.5649562,-0.2929721],[0.7713572,0.5649562,-0.2929721],[0.7713572,0.5649562,-0.2929721],[0,1,0],[0,1,0],[0,1,0],[0,1,0],[0,1,0],[0,1,0],[0.8994442,0.2710145,-0.3428575],[0.8994442,0.2710145,-0.3428575],[0.8994442,0.2710145,-0.3428575],[0.8994442,0.2710145,-0.3428575],[0.8994442,0.2710145,-0.3428575],[0.8994442,0.2710145,-0.3428575],[0.9574969,0.2884437,-4.845662e-05],[0.9574969,0.2884437,-4.845662e-05],[0.9574969,0.2884437,-4.845662e-05],[0.9574718,0.2885268,0],[0.9574718,0.2885268,0],[0.9574718,0.2885268,0],[0.9574718,0.2885268,0],[0.9574718,0.2885268,0],[0.9574718,0.2885268,0],[0.9574771,0.288509,-0],[0.9574771,0.288509,-0],[0.9574771,0.288509,-0],[0.9574771,0.288509,0],[0.9574771,0.288509,0],[0.9574771,0.288509,0],[0.9575121,0.2883932,-0],[0.9575121,0.2883932,-0],[0.9575121,0.2883932,-0],[0.9575121,0.2883931,0],[0.9575121,0.2883931,0],[0.9575121,0.2883931,0],[0.957477,0.288509,0.0005070638],[0.957477,0.288509,0.0005070638],[0.957477,0.288509,0.0005070638],[0.9574795,0.2885013,-0.0001401839],[0.9574795,0.2885013,-0.0001401839],[0.9574795,0.2885013,-0.0001401839],[0.9574807,0.2884972,-2.838519e-05],[0.9574807,0.2884972,-2.838519e-05],[0.9574807,0.2884972,-2.838519e-05],[0.9574849,0.2884834,0],[0.9574849,0.2884834,0],[0.9574849,0.2884834,0],[0.9574849,0.2884834,0],[0.9574849,0.2884834,0],[0.9574849,0.2884834,0],[0.9574848,0.2884834,0],[0.9574848,0.2884834,0],[0.9574848,0.2884834,0],[0.9574849,0.2884834,0],[0.9574849,0.2884834,0],[0.9574849,0.2884834,0],[0.9574804,0.2884983,2.734142e-05],[0.9574804,0.2884983,2.734142e-05],[0.9574804,0.2884983,2.734142e-05],[0.9574795,0.2885013,6.619826e-05],[0.9574795,0.2885013,6.619826e-05],[0.9574795,0.2885013,6.619826e-05],[0.9574771,0.288509,-0.0002409421],[0.9574771,0.288509,-0.0002409421],[0.9574771,0.288509,-0.0002409421],[0.9574984,0.2884386,4.354864e-05],[0.9574984,0.2884386,4.354864e-05],[0.9574984,0.2884386,4.354864e-05],[0.9572536,0.28925,0],[0.9572536,0.28925,0],[0.9572536,0.28925,0],[0.9572536,0.28925,-0],[0.9572536,0.28925,-0],[0.9572536,0.28925,-0],[0.9572537,0.28925,-0],[0.9572537,0.28925,-0],[0.9572537,0.28925,-0],[0.9572536,0.28925,0],[0.9572536,0.28925,0],[0.9572536,0.28925,0],[0.8662896,-0.4669726,0.1774228],[0.8662896,-0.4669726,0.1774228],[0.8662896,-0.4669726,0.1774228],[0.8684785,-0.4641916,0.1739866],[0.8684785,-0.4641916,0.1739866],[0.8684785,-0.4641916,0.1739866],[0.8372375,0.5468395,0],[0.8372375,0.5468395,0],[0.8372375,0.5468395,0],[0.8372375,0.5468395,0],[0.8372375,0.5468395,0],[0.8372375,0.5468395,0],[0,-0.4105298,-0.9118472],[0,-0.4105298,-0.9118472],[0,-0.4105298,-0.9118472],[0,-0.4105298,-0.9118472],[0,-0.4105298,-0.9118472],[0,-0.4105298,-0.9118472],[0,0.4107764,-0.9117361],[0,0.4107764,-0.9117361],[0,0.4107764,-0.9117361],[0,0.4107764,-0.9117361],[0,0.4107764,-0.9117361],[0,0.4107764,-0.9117361],[0.3825871,0.9239194,-0],[0.3825871,0.9239194,-0],[0.3825871,0.9239194,-0],[0.3825871,0.9239194,0],[0.3825871,0.9239194,0],[0.3825871,0.9239194,0],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,-0,1],[0,-0,1],[0,-0,1],[0,-0,1],[0,-0,1],[0,-0,1],[-0,0,1],[-0,0,1],[-0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,-0,1],[0,-0,1],[0,-0,1],[-0,0,1],[-0,0,1],[-0,0,1],[-0,0,1],[-0,0,1],[-0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[-0,0,0.9999999],[-0,0,0.9999999],[-0,0,0.9999999],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[-0,0,0.9999999],[-0,0,0.9999999],[-0,0,0.9999999],[-0,0,1],[-0,0,1],[-0,0,1],[-0,0,1],[-0,0,1],[-0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0.9767476,-0.2143927,0],[0.9767476,-0.2143927,0],[0.9767476,-0.2143927,0],[0.9767476,-0.2143927,0],[0.9767476,-0.2143927,0],[0.9767476,-0.2143927,0],[-0.8643441,-0.5029009,0],[-0.8643441,-0.5029009,0],[-0.8643441,-0.5029009,0],[-0.8643441,-0.5029009,0],[-0.8643441,-0.5029009,0],[-0.8643441,-0.5029009,0],[-0.9931883,-0.1165206,0],[-0.9931883,-0.1165206,0],[-0.9931883,-0.1165206,0],[-0.9931883,-0.1165206,0],[-0.9931883,-0.1165206,0],[-0.9931883,-0.1165206,0],[-0.980085,-0.1985785,0],[-0.980085,-0.1985785,0],[-0.980085,-0.1985785,0],[-0.980085,-0.1985785,0],[-0.980085,-0.1985785,0],[-0.980085,-0.1985785,0],[0,0,-1],[0,0,-1],[0,0,-1],[0,-0,-1],[0,-0,-1],[0,-0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[0,0,1],[0,0,1],[0,0,1],[-0,0,1],[-0,0,1],[-0,0,1],[0,-1,-0],[0,-1,-0],[0,-1,-0],[0,-1,0],[0,-1,0],[0,-1,0],[0.88901,-0.4578879,0],[0.88901,-0.4578879,0],[0.88901,-0.4578879,0],[0.88901,-0.4578879,0],[0.88901,-0.4578879,0],[0.88901,-0.4578879,0],[0.8701845,-0.492726,0],[0.8701845,-0.492726,0],[0.8701845,-0.492726,0],[0.8701845,-0.492726,0],[0.8701845,-0.492726,0],[0.8701845,-0.492726,0],[-0.8161331,-0.5778638,0],[-0.8161331,-0.5778638,0],[-0.8161331,-0.5778638,0],[-0.8161331,-0.5778638,0],[-0.8161331,-0.5778638,0],[-0.8161331,-0.5778638,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0.9810797,-0.1936043,0],[0.9810797,-0.1936043,0],[0.9810797,-0.1936043,0],[0.9810798,-0.1936043,0],[0.9810798,-0.1936043,0],[0.9810798,-0.1936043,0],[-0.9796228,0.2008462,0],[-0.9796228,0.2008462,0],[-0.9796228,0.2008462,0],[-0.9796228,0.2008461,0],[-0.9796228,0.2008461,0],[-0.9796228,0.2008461,0],[0.9796228,0.2008462,0],[0.9796228,0.2008462,0],[0.9796228,0.2008462,0],[0.9796228,0.2008461,-0],[0.9796228,0.2008461,-0],[0.9796228,0.2008461,-0],[0.962714,0.2705211,0],[0.962714,0.2705211,0],[0.962714,0.2705211,0],[0.9627141,0.2705211,0],[0.9627141,0.2705211,0],[0.9627141,0.2705211,0],[-0.9345582,0,0.3558106],[-0.9345582,0,0.3558106],[-0.9345582,0,0.3558106],[-0.9345582,0,0.3558106],[-0.9345582,0,0.3558106],[-0.9345582,0,0.3558106],[-0.4103398,-0,-0.9119328],[-0.4103398,-0,-0.9119328],[-0.4103398,-0,-0.9119328],[-0.4103398,0,-0.9119328],[-0.4103398,0,-0.9119328],[-0.4103398,0,-0.9119328],[0.3468219,-0.3132194,-0.8840861],[0.3468219,-0.3132194,-0.8840861],[0.3468219,-0.3132194,-0.8840861],[0.3857889,-0.3406662,-0.8573875],[0.3857889,-0.3406662,-0.8573875],[0.3857889,-0.3406662,-0.8573875],[-0.3608875,0.8714774,-0.3320955],[-0.3608875,0.8714774,-0.3320955],[-0.3608875,0.8714774,-0.3320955],[-0.1847878,0.9679556,-0.1700452],[-0.1847878,0.9679556,-0.1700452],[-0.1847878,0.9679556,-0.1700452],[-0.8716245,-0.360758,0.3318501],[-0.8716245,-0.360758,0.3318501],[-0.8716245,-0.360758,0.3318501],[-0.9676845,-0.1855876,0.1707163],[-0.9676845,-0.1855876,0.1707163],[-0.9676845,-0.1855876,0.1707163],[0.5689232,-0.7584897,0.3178359],[0.5689232,-0.7584897,0.3178359],[0.5689232,-0.7584897,0.3178359],[0.7355877,-0.6612584,0.1471323],[0.7355877,-0.6612584,0.1471323],[0.7355877,-0.6612584,0.1471323],[0.7568582,-0.6535791,0],[0.7568582,-0.6535791,0],[0.7568582,-0.6535791,0],[0.7568582,-0.6535791,0],[0.7568582,-0.6535791,0],[0.7568582,-0.6535791,0],[0,1,-0],[0,1,-0],[0,1,-0],[0,1,0],[0,1,0],[0,1,0],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[-0.9999999,-0,0],[-0.9999999,-0,0],[-0.9999999,-0,0],[-0.9999999,0,0],[-0.9999999,0,0],[-0.9999999,0,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-0.9775435,0.210733],[0,-0.9775435,0.210733],[0,-0.9775435,0.210733],[0,-0.9775435,0.210733],[0,-0.9775435,0.210733],[0,-0.9775435,0.210733],[-0,0.9775437,0.210733],[-0,0.9775437,0.210733],[-0,0.9775437,0.210733],[0,0.9775437,0.210733],[0,0.9775437,0.210733],[0,0.9775437,0.210733],[0,0.909798,0.4150513],[0,0.909798,0.4150513],[0,0.909798,0.4150513],[0,0.909798,0.4150513],[0,0.909798,0.4150513],[0,0.909798,0.4150513],[0,0.9097998,0.4150473],[0,0.9097998,0.4150473],[0,0.9097998,0.4150473],[0,0.9097998,0.4150473],[0,0.9097998,0.4150473],[0,0.9097998,0.4150473],[-0.9999999,-0,0],[-0.9999999,-0,0],[-0.9999999,-0,0],[-0.9999999,0,0],[-0.9999999,0,0],[-0.9999999,0,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[0,-1,0],[-1,0,0],[-1,0,0],[-1,0,0],[-1,-0,0],[-1,-0,0],[-1,-0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,-0],[1,0,-0],[1,0,-0],[0.8814226,-0.4297242,-0.1960386],[0.8814226,-0.4297242,-0.1960386],[0.8814226,-0.4297242,-0.1960386],[0.8814226,-0.4297242,-0.1960386],[0.8814226,-0.4297242,-0.1960386],[0.8814226,-0.4297242,-0.1960386],[0,-0.9097999,-0.4150473],[0,-0.9097999,-0.4150473],[0,-0.9097999,-0.4150473],[0,-0.9097999,-0.4150473],[0,-0.9097999,-0.4150473],[0,-0.9097999,-0.4150473],[0,0,-1],[0,0,-1],[0,0,-1],[-0,0,-1],[-0,0,-1],[-0,0,-1],[0,0,1],[0,0,1],[0,0,1],[-0,0,1],[-0,0,1],[-0,0,1],[0.5554983,-0.8315177,0],[0.5554983,-0.8315177,0],[0.5554983,-0.8315177,0],[0.5554983,-0.8315177,0],[0.5554983,-0.8315177,0],[0.5554983,-0.8315177,0],[-1,0,0],[-1,0,0],[-1,0,0],[-1,0,0],[-1,0,0],[-1,0,0],[-0.5554988,-0.8315173,-1.250832e-05],[-0.5554988,-0.8315173,-1.250832e-05],[-0.5554988,-0.8315173,-1.250832e-05],[-0.5554988,-0.8315173,-1.250832e-05],[-0.5554988,-0.8315173,-1.250832e-05],[-0.5554988,-0.8315173,-1.250832e-05],[0.2404871,-0,0.9706523],[0.2404871,-0,0.9706523],[0.2404871,-0,0.9706523],[0.2404871,0,0.9706523],[0.2404871,0,0.9706523],[0.2404871,0,0.9706523],[0.3403893,0.8786296,-0.3348809],[0.3403893,0.8786296,-0.3348809],[0.3403893,0.8786296,-0.3348809],[0.3440206,0.8780348,-0.3327233],[0.3440206,0.8780348,-0.3327233],[0.3440206,0.8780348,-0.3327233],[0,0.9344294,0.3561485],[0,0.9344294,0.3561485],[0,0.9344294,0.3561485],[-0,0.9344294,0.3561485],[-0,0.9344294,0.3561485],[-0,0.9344294,0.3561485],[-1,0,0],[-1,0,0],[-1,0,0],[-1,0,0],[-1,0,0],[-1,0,0],[0.3748817,0.4137425,0.8296269],[0.3748817,0.4137425,0.8296269],[0.3748817,0.4137425,0.8296269],[0.3782091,0.4111477,0.8294067],[0.3782091,0.4111477,0.8294067],[0.3782091,0.4111477,0.8294067],[0.3137324,0.6449381,-0.6968693],[0.3137324,0.6449381,-0.6968693],[0.3137324,0.6449381,-0.6968693],[0.3175517,0.6521124,-0.6884115],[0.3175517,0.6521124,-0.6884115],[0.3175517,0.6521124,-0.6884115],[0.6689631,0.7432957,0],[0.6689631,0.7432957,0],[0.6689631,0.7432957,0],[0.6689631,0.7432957,0],[0.6689631,0.7432957,0],[0.6689631,0.7432957,0],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[0,0,-1],[-1,0,0],[-1,0,0],[-1,0,0],[-0.9999999,0,0],[-0.9999999,0,0],[-0.9999999,0,0],[-1,0,0],[-1,0,0],[-1,0,0],[-1,0,0],[-1,0,0],[-1,0,0],[-1,0,-0],[-1,0,-0],[-1,0,-0],[-0.9999999,0,0],[-0.9999999,0,0],[-0.9999999,0,0],[-1,-0,-0],[-1,-0,-0],[-1,-0,-0],[-1,0,0],[-1,0,0],[-1,0,0],[-1,0,0],[-1,0,0],[-1,0,0],[-1,-0,0],[-1,-0,0],[-1,-0,0],[-1,0,0],[-1,0,0],[-1,0,0],[-1,-0,0],[-1,-0,0],[-1,-0,0],[-0.5172721,0.8329639,-0.1964703],[-0.5172721,0.8329639,-0.1964703],[-0.5172721,0.8329639,-0.1964703],[-0.5171735,0.8330099,-0.1965351],[-0.5171735,0.8330099,-0.1965351],[-0.5171735,0.8330099,-0.1965351],[0.8281263,0.5605416,0],[0.8281263,0.5605416,0],[0.8281263,0.5605416,0],[0.8281263,0.5605417,-0],[0.8281263,0.5605417,-0],[0.8281263,0.5605417,-0],[0,1,0],[0,1,0],[0,1,0],[-0,1,0],[-0,1,0],[-0,1,0],[0.9280494,0.1165212,0.3537614],[0.9280494,0.1165212,0.3537614],[0.9280494,0.1165212,0.3537614],[0.9285359,0.1151842,0.3529214],[0.9285359,0.1151842,0.3529214],[0.9285359,0.1151842,0.3529214],[0.7618073,0,0.6478037],[0.7618073,0,0.6478037],[0.7618073,0,0.6478037],[0.7618073,0,0.6478037],[0.7618073,0,0.6478037],[0.7618073,0,0.6478037],[0.837243,-0.5468311,0],[0.837243,-0.5468311,0],[0.837243,-0.5468311,0],[0.837243,-0.5468311,0],[0.837243,-0.5468311,0],[0.837243,-0.5468311,0],[0,-0.4105298,0.9118472],[0,-0.4105298,0.9118472],[0,-0.4105298,0.9118472],[0,-0.4105298,0.9118472],[0,-0.4105298,0.9118472],[0,-0.4105298,0.9118472],[0,-0.9343773,0.3562853],[0,-0.9343773,0.3562853],[0,-0.9343773,0.3562853],[0,-0.9343773,0.3562853],[0,-0.9343773,0.3562853],[0,-0.9343773,0.3562853],[0,1,0],[0,1,0],[0,1,0],[0,1,0],[0,1,0],[0,1,0]],"ignoreExtent":false,"flags":3},"10":{"id":10,"type":"light","vertices":[[0,0,1]],"colors":[[1,1,1,1],[1,1,1,1],[1,1,1,1]],"viewpoint":true,"finite":false},"9":{"id":9,"type":"background","material":{"fog":true},"colors":[[0.2980392,0.2980392,0.2980392,1]],"centers":[[0,0,0]],"sphere":false,"fogtype":"none","flags":0},"11":{"id":11,"type":"background","material":{"lit":false,"back":"lines"},"colors":[[1,1,1,1]],"centers":[[0,0,0]],"sphere":false,"fogtype":"none","flags":0},"6":{"id":6,"type":"subscene","par3d":{"antialias":8,"FOV":0,"ignoreExtent":false,"listeners":6,"mouseMode":{"left":"trackball","right":"zoom","middle":"fov","wheel":"pull"},"observer":[0,0,6.072288],"modelMatrix":[[0.5,0,-0.8660254,-0.001772463],[-0.2241439,0.9659258,-0.1294095,-1.376215],[0.8365163,0.258819,0.4829629,-6.444221],[0,0,0,1]],"projMatrix":[[0.3360869,0,0,0],[0,0.4705217,0,0],[0,0,-0.3293652,-2],[0,0,0,1]],"skipRedraw":false,"userMatrix":[[0.5,0,-0.8660254,0],[-0.2241439,0.9659258,-0.1294095,0],[0.8365163,0.258819,0.4829629,0],[0,0,0,1]],"userProjection":[[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]],"scale":[1,1,1],"viewport":{"x":0,"y":0,"width":1,"height":1},"zoom":0.7,"bbox":[-2.43091,2.438,0.2668,2.58437,-1.396,1.396],"windowRect":[60,105,732,585],"family":"sans","font":1,"cex":1,"useFreeType":true,"fontname":"/Library/Frameworks/R.framework/Versions/4.0/Resources/library/rgl/fonts/FreeSans.ttf","maxClipPlanes":6,"glVersion":2.1,"activeSubscene":0},"embeddings":{"viewport":"replace","projection":"replace","model":"replace","mouse":"replace"},"objects":[11,12,10],"subscenes":[],"flags":515}},"snapshot":"data:image/png;base64,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","width":673,"height":481,"sphereVerts":{"vb":[[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.07465783,0.1464466,0.2126075,0.2705981,0.3181896,0.3535534,0.3753303,0.3826834,0.3753303,0.3535534,0.3181896,0.2705981,0.2126075,0.1464466,0.07465783,0,0,0.1379497,0.2705981,0.3928475,0.5,0.5879378,0.6532815,0.6935199,0.7071068,0.6935199,0.6532815,0.5879378,0.5,0.3928475,0.2705981,0.1379497,0,0,0.18024,0.3535534,0.51328,0.6532815,0.7681778,0.8535534,0.9061274,0.9238795,0.9061274,0.8535534,0.7681778,0.6532815,0.51328,0.3535534,0.18024,0,0,0.1950903,0.3826834,0.5555702,0.7071068,0.8314696,0.9238795,0.9807853,1,0.9807853,0.9238795,0.8314696,0.7071068,0.5555702,0.3826834,0.1950903,0,0,0.18024,0.3535534,0.51328,0.6532815,0.7681778,0.8535534,0.9061274,0.9238795,0.9061274,0.8535534,0.7681778,0.6532815,0.51328,0.3535534,0.18024,0,0,0.1379497,0.2705981,0.3928475,0.5,0.5879378,0.6532815,0.6935199,0.7071068,0.6935199,0.6532815,0.5879378,0.5,0.3928475,0.2705981,0.1379497,0,0,0.07465783,0.1464466,0.2126075,0.2705981,0.3181896,0.3535534,0.3753303,0.3826834,0.3753303,0.3535534,0.3181896,0.2705981,0.2126075,0.1464466,0.07465783,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0,-0.07465783,-0.1464466,-0.2126075,-0.2705981,-0.3181896,-0.3535534,-0.3753303,-0.3826834,-0.3753303,-0.3535534,-0.3181896,-0.2705981,-0.2126075,-0.1464466,-0.07465783,-0,-0,-0.1379497,-0.2705981,-0.3928475,-0.5,-0.5879378,-0.6532815,-0.6935199,-0.7071068,-0.6935199,-0.6532815,-0.5879378,-0.5,-0.3928475,-0.2705981,-0.1379497,-0,-0,-0.18024,-0.3535534,-0.51328,-0.6532815,-0.7681778,-0.8535534,-0.9061274,-0.9238795,-0.9061274,-0.8535534,-0.7681778,-0.6532815,-0.51328,-0.3535534,-0.18024,-0,-0,-0.1950903,-0.3826834,-0.5555702,-0.7071068,-0.8314696,-0.9238795,-0.9807853,-1,-0.9807853,-0.9238795,-0.8314696,-0.7071068,-0.5555702,-0.3826834,-0.1950903,-0,-0,-0.18024,-0.3535534,-0.51328,-0.6532815,-0.7681778,-0.8535534,-0.9061274,-0.9238795,-0.9061274,-0.8535534,-0.7681778,-0.6532815,-0.51328,-0.3535534,-0.18024,-0,-0,-0.1379497,-0.2705981,-0.3928475,-0.5,-0.5879378,-0.6532815,-0.6935199,-0.7071068,-0.6935199,-0.6532815,-0.5879378,-0.5,-0.3928475,-0.2705981,-0.1379497,-0,-0,-0.07465783,-0.1464466,-0.2126075,-0.2705981,-0.3181896,-0.3535534,-0.3753303,-0.3826834,-0.3753303,-0.3535534,-0.3181896,-0.2705981,-0.2126075,-0.1464466,-0.07465783,-0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[-1,-0.9807853,-0.9238795,-0.8314696,-0.7071068,-0.5555702,-0.3826834,-0.1950903,0,0.1950903,0.3826834,0.5555702,0.7071068,0.8314696,0.9238795,0.9807853,1,-1,-0.9807853,-0.9238795,-0.8314696,-0.7071068,-0.5555702,-0.3826834,-0.1950903,0,0.1950903,0.3826834,0.5555702,0.7071068,0.8314696,0.9238795,0.9807853,1,-1,-0.9807853,-0.9238795,-0.8314696,-0.7071068,-0.5555702,-0.3826834,-0.1950903,0,0.1950903,0.3826834,0.5555702,0.7071068,0.8314696,0.9238795,0.9807853,1,-1,-0.9807853,-0.9238795,-0.8314696,-0.7071068,-0.5555702,-0.3826834,-0.1950903,0,0.1950903,0.3826834,0.5555702,0.7071068,0.8314696,0.9238795,0.9807853,1,-1,-0.9807853,-0.9238795,-0.8314696,-0.7071068,-0.5555702,-0.3826834,-0.1950903,0,0.1950903,0.3826834,0.5555702,0.7071068,0.8314696,0.9238795,0.9807853,1,-1,-0.9807853,-0.9238795,-0.8314696,-0.7071068,-0.5555702,-0.3826834,-0.1950903,0,0.1950903,0.3826834,0.5555702,0.7071068,0.8314696,0.9238795,0.9807853,1,-1,-0.9807853,-0.9238795,-0.8314696,-0.7071068,-0.5555702,-0.3826834,-0.1950903,0,0.1950903,0.3826834,0.5555702,0.7071068,0.8314696,0.9238795,0.9807853,1,-1,-0.9807853,-0.9238795,-0.8314696,-0.7071068,-0.5555702,-0.3826834,-0.1950903,0,0.1950903,0.3826834,0.5555702,0.7071068,0.8314696,0.9238795,0.9807853,1,-1,-0.9807853,-0.9238795,-0.8314696,-0.7071068,-0.5555702,-0.3826834,-0.1950903,0,0.1950903,0.3826834,0.5555702,0.7071068,0.8314696,0.9238795,0.9807853,1,-1,-0.9807853,-0.9238795,-0.8314696,-0.7071068,-0.5555702,-0.3826834,-0.1950903,0,0.1950903,0.3826834,0.5555702,0.7071068,0.8314696,0.9238795,0.9807853,1,-1,-0.9807853,-0.9238795,-0.8314696,-0.7071068,-0.5555702,-0.3826834,-0.1950903,0,0.1950903,0.3826834,0.5555702,0.7071068,0.8314696,0.9238795,0.9807853,1,-1,-0.9807853,-0.9238795,-0.8314696,-0.7071068,-0.5555702,-0.3826834,-0.1950903,0,0.1950903,0.3826834,0.5555702,0.7071068,0.8314696,0.9238795,0.9807853,1,-1,-0.9807853,-0.9238795,-0.8314696,-0.7071068,-0.5555702,-0.3826834,-0.1950903,0,0.1950903,0.3826834,0.5555702,0.7071068,0.8314696,0.9238795,0.9807853,1,-1,-0.9807853,-0.9238795,-0.8314696,-0.7071068,-0.5555702,-0.3826834,-0.1950903,0,0.1950903,0.3826834,0.5555702,0.7071068,0.8314696,0.9238795,0.9807853,1,-1,-0.9807853,-0.9238795,-0.8314696,-0.7071068,-0.5555702,-0.3826834,-0.1950903,0,0.1950903,0.3826834,0.5555702,0.7071068,0.8314696,0.9238795,0.9807853,1,-1,-0.9807853,-0.9238795,-0.8314696,-0.7071068,-0.5555702,-0.3826834,-0.1950903,0,0.1950903,0.3826834,0.5555702,0.7071068,0.8314696,0.9238795,0.9807853,1,-1,-0.9807853,-0.9238795,-0.8314696,-0.7071068,-0.5555702,-0.3826834,-0.1950903,0,0.1950903,0.3826834,0.5555702,0.7071068,0.8314696,0.9238795,0.9807853,1],[0,0.1950903,0.3826834,0.5555702,0.7071068,0.8314696,0.9238795,0.9807853,1,0.9807853,0.9238795,0.8314696,0.7071068,0.5555702,0.3826834,0.1950903,0,0,0.18024,0.3535534,0.51328,0.6532815,0.7681778,0.8535534,0.9061274,0.9238795,0.9061274,0.8535534,0.7681778,0.6532815,0.51328,0.3535534,0.18024,0,0,0.1379497,0.2705981,0.3928475,0.5,0.5879378,0.6532815,0.6935199,0.7071068,0.6935199,0.6532815,0.5879378,0.5,0.3928475,0.2705981,0.1379497,0,0,0.07465783,0.1464466,0.2126075,0.2705981,0.3181896,0.3535534,0.3753303,0.3826834,0.3753303,0.3535534,0.3181896,0.2705981,0.2126075,0.1464466,0.07465783,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0,-0.07465783,-0.1464466,-0.2126075,-0.2705981,-0.3181896,-0.3535534,-0.3753303,-0.3826834,-0.3753303,-0.3535534,-0.3181896,-0.2705981,-0.2126075,-0.1464466,-0.07465783,-0,-0,-0.1379497,-0.2705981,-0.3928475,-0.5,-0.5879378,-0.6532815,-0.6935199,-0.7071068,-0.6935199,-0.6532815,-0.5879378,-0.5,-0.3928475,-0.2705981,-0.1379497,-0,-0,-0.18024,-0.3535534,-0.51328,-0.6532815,-0.7681778,-0.8535534,-0.9061274,-0.9238795,-0.9061274,-0.8535534,-0.7681778,-0.6532815,-0.51328,-0.3535534,-0.18024,-0,-0,-0.1950903,-0.3826834,-0.5555702,-0.7071068,-0.8314696,-0.9238795,-0.9807853,-1,-0.9807853,-0.9238795,-0.8314696,-0.7071068,-0.5555702,-0.3826834,-0.1950903,-0,-0,-0.18024,-0.3535534,-0.51328,-0.6532815,-0.7681778,-0.8535534,-0.9061274,-0.9238795,-0.9061274,-0.8535534,-0.7681778,-0.6532815,-0.51328,-0.3535534,-0.18024,-0,-0,-0.1379497,-0.2705981,-0.3928475,-0.5,-0.5879378,-0.6532815,-0.6935199,-0.7071068,-0.6935199,-0.6532815,-0.5879378,-0.5,-0.3928475,-0.2705981,-0.1379497,-0,-0,-0.07465783,-0.1464466,-0.2126075,-0.2705981,-0.3181896,-0.3535534,-0.3753303,-0.3826834,-0.3753303,-0.3535534,-0.3181896,-0.2705981,-0.2126075,-0.1464466,-0.07465783,-0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.07465783,0.1464466,0.2126075,0.2705981,0.3181896,0.3535534,0.3753303,0.3826834,0.3753303,0.3535534,0.3181896,0.2705981,0.2126075,0.1464466,0.07465783,0,0,0.1379497,0.2705981,0.3928475,0.5,0.5879378,0.6532815,0.6935199,0.7071068,0.6935199,0.6532815,0.5879378,0.5,0.3928475,0.2705981,0.1379497,0,0,0.18024,0.3535534,0.51328,0.6532815,0.7681778,0.8535534,0.9061274,0.9238795,0.9061274,0.8535534,0.7681778,0.6532815,0.51328,0.3535534,0.18024,0,0,0.1950903,0.3826834,0.5555702,0.7071068,0.8314696,0.9238795,0.9807853,1,0.9807853,0.9238795,0.8314696,0.7071068,0.5555702,0.3826834,0.1950903,0]],"it":[[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,238,239,240,241,242,243,244,245,246,247,248,249,250,251,252,253,255,256,257,258,259,260,261,262,263,264,265,266,267,268,269,270,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,238,239,240,241,242,243,244,245,246,247,248,249,250,251,252,253,255,256,257,258,259,260,261,262,263,264,265,266,267,268,269,270],[17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,238,239,240,241,242,243,244,245,246,247,248,249,250,251,252,253,255,256,257,258,259,260,261,262,263,264,265,266,267,268,269,270,272,273,274,275,276,277,278,279,280,281,282,283,284,285,286,287,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,239,240,241,242,243,244,245,246,247,248,249,250,251,252,253,254,256,257,258,259,260,261,262,263,264,265,266,267,268,269,270,271,273,274,275,276,277,278,279,280,281,282,283,284,285,286,287,288],[18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,239,240,241,242,243,244,245,246,247,248,249,250,251,252,253,254,256,257,258,259,260,261,262,263,264,265,266,267,268,269,270,271,273,274,275,276,277,278,279,280,281,282,283,284,285,286,287,288,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,239,240,241,242,243,244,245,246,247,248,249,250,251,252,253,254,256,257,258,259,260,261,262,263,264,265,266,267,268,269,270,271]],"material":[],"normals":null,"texcoords":[[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.1875,0.1875,0.1875,0.1875,0.1875,0.1875,0.1875,0.1875,0.1875,0.1875,0.1875,0.1875,0.1875,0.1875,0.1875,0.1875,0.1875,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.3125,0.3125,0.3125,0.3125,0.3125,0.3125,0.3125,0.3125,0.3125,0.3125,0.3125,0.3125,0.3125,0.3125,0.3125,0.3125,0.3125,0.375,0.375,0.375,0.375,0.375,0.375,0.375,0.375,0.375,0.375,0.375,0.375,0.375,0.375,0.375,0.375,0.375,0.4375,0.4375,0.4375,0.4375,0.4375,0.4375,0.4375,0.4375,0.4375,0.4375,0.4375,0.4375,0.4375,0.4375,0.4375,0.4375,0.4375,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5625,0.5625,0.5625,0.5625,0.5625,0.5625,0.5625,0.5625,0.5625,0.5625,0.5625,0.5625,0.5625,0.5625,0.5625,0.5625,0.5625,0.625,0.625,0.625,0.625,0.625,0.625,0.625,0.625,0.625,0.625,0.625,0.625,0.625,0.625,0.625,0.625,0.625,0.6875,0.6875,0.6875,0.6875,0.6875,0.6875,0.6875,0.6875,0.6875,0.6875,0.6875,0.6875,0.6875,0.6875,0.6875,0.6875,0.6875,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.8125,0.8125,0.8125,0.8125,0.8125,0.8125,0.8125,0.8125,0.8125,0.8125,0.8125,0.8125,0.8125,0.8125,0.8125,0.8125,0.8125,0.875,0.875,0.875,0.875,0.875,0.875,0.875,0.875,0.875,0.875,0.875,0.875,0.875,0.875,0.875,0.875,0.875,0.9375,0.9375,0.9375,0.9375,0.9375,0.9375,0.9375,0.9375,0.9375,0.9375,0.9375,0.9375,0.9375,0.9375,0.9375,0.9375,0.9375,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1]],"meshColor":"vertices"},"context":{"shiny":false,"rmarkdown":"github_document"},"crosstalk":{"key":[],"group":[],"id":[],"options":[]}});
unnamed_chunk_5rgl.prefix = "unnamed_chunk_5";
</script>

<p id="unnamed_chunk_5debug">

You must enable Javascript to view this page properly.

</p>

<script>unnamed_chunk_5rgl.start();</script>

## Some details

**dracor** is deliberately intended as a minimal decoder package without
any dependencies besides the Rcpp package. It accepts raw bytes, a file
or a URL as input and can produce either an
[rgl](https://cran.r-project.org/package=rgl) `mesh3d` object as output
or a list containing points and 0-indexed faces. It essentially
replicates the most basic decoding ability of the `draco_decoder`
command line tool.

If you just want a result as close as possible to what the Draco library
would give then set `mesh3d=FALSE`

``` r
car.m=dracor::draco_decode(carurl, mesh3d=FALSE)
str(car.m)
#> List of 2
#>  $ points: num [1:3, 1:1856] 1.54 1.65 -1.21 1.57 1.77 ...
#>  $ faces : int [1:3, 1:1744] 0 1 2 2 1 3 3 1 4 4 ...
```

## Acknowledgements

Many thanks to the authors of:

  - [Draco library](https://github.com/google/draco)
  - [Rcpp package](https://cran.r-project.org/package=Rcpp)
