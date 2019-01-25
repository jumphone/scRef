<a href="/">
<img src="/Reference/images/NewLogo.png" width="280">
</a>

* A reference-based toolkit for single-cell analysis <a href='https://github.com/jumphone/scRef/wiki#table-of-content'>Visit Wiki</a>

* A web resource of references. <a href='/Reference/'>Visit dbRef</a>

# Table of content

* [Citation](#Citation)
* [Download](#Download)
* [Requirement](#Requirement)
* [Input format](#Input-format)
* [Usage (Wiki)](https://github.com/jumphone/scRef/wiki#table-of-content)
* [License](#License)


# Citation

Feng Zhang, Yaguang Dou, Rohit R Rao, Q. Richard Lu, Weidong Tian; SCREF: a reference-based solution for single-cell analysis, Coming Soon

# Download

### Whole repository:

Download: https://github.com/jumphone/scRef/archive/master.zip

### Only "scRef.R":

curl https://raw.githubusercontent.com/jumphone/scRef/master/scRef.R > scRef.R


# Requirement

R: 3.5.0

Denpendent R packages: pcaPP, MASS, ggplot2, igraph, pastecs

# Input format

Single-cell and reference expression matrix must be tab-delimited with header and row names.

For reference expression matrix, we recommend using TPM, FPKM, RPKM, RPK, or UMI matrix. 

For single-cell expression matrix, we recommend using UMI matrix.

# Usage
   
See [Wiki](https://github.com/jumphone/scRef/wiki#table-of-content)

   
# License

    Copyright (c) 2019 Zhang, Feng

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.





