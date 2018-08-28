## SFAmix

SFAmix is a sparse matrix decomposition tool. Given a matrix **Y** with dimension of **P** by **N**, BicMix decompose it into the product of two sparse matrices **LAM** and dense matrix **X**

## How to install

git clone https://github.com/chuangao/SFAmixC.git <br/>
cd SFAmixC <br/>
make <br/>

This will generate the executible **SFAmix** (The compilation may take a few minutes as it also compiles a local copy of gsl) <br/>

For linux, the default compilation have openmp turned on <br/>
For macOS, openmp is turned off because it is not supported by default (yet) <br/>

## Usage
### To run without using openmp
./SFAmix --y input_file --out output_directory <br/>

**Please no headers in the input matrix, no missing values, just pure numbers, ideally quantile normalized** <br/>
**Also no corrections of confounding beforehand, BicMix will handle that in the dense components** <br/>
**For a gene expression matrix, it is prefered that each gene is a row and each sample is a column** <br/> 

### To run using openmp <br/>
OMP_NUM_THREADS=10 ./SFAmix --y input_file --out output_directory <br/>

### Parameters accepted
**--y** The matrix that is to be decomposed <br/>
**--out** The output directory to dump the results <br/>
**--nf** The initial number of factors that the algorithm should start with, BicMix will shrink this number to the true number (assume such a number can be found) <br/>
**--itr** The maximum number of iterations to run, default to 5000 <br/>
**--interval** How many iterations that the number of nonzero elements remain unchanged (to claim that the algorithm converges), default to 500 iterations <br/>
**--write_itr** Every how many iterations that the result will be written, default to 50 <br/>

### Output files
**LAM** The loading matrix <br/>
**EX** The factor matrix <br/>
**Z** A vector indicating whether the component (column) in **LAM** is sparse <br/>
**count_lam** the number of nonzero elements in each component(column) of **LAM** <br/>
**PSI** Variance for the error matrix <br/>
**EXX** E(XX^T) <br/>

#### The many iterations of output is for temporary access to intermediate results when the algorithm runs too long 


