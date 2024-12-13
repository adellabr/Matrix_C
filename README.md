# Matrix_C
Implementation of the matrix.h library

## Matrix

A matrix is a collection of numbers arranged into a fixed number of rows and columns.

Matrix A is a rectangular table of numbers arranged in m rows and n columns

### Matrix structure in C language:

```c
typedef struct matrix_struct {
    double** matrix;
    int rows;
    int columns;
} matrix_t;
```
## Matrix operations

All operations (except matrix comparison) return the resulting code:
- 0 - OK
- 1 - Error, incorrect matrix
- 2 - Calculation error (mismatched matrix sizes; matrix for which calculations cannot be performed, etc.)

### Creating matrices (create_matrix)

```c
int s21_create_matrix(int rows, int columns, matrix_t *result);
```

### Cleaning of matrices (remove_matrix)

```c
void s21_remove_matrix(matrix_t *A);
```

### Matrix comparison (eq_matrix)

```c
#define SUCCESS 1
#define FAILURE 0

int s21_eq_matrix(matrix_t *A, matrix_t *B);
```

The matrices A, B are equal |A = B| if they have the same dimensions and the corresponding elements are identical, thus for all i and j: A(i,j) = B(i,j)

The comparison is including 7 decimal places.

### Adding (sum_matrix) and subtracting matrices (sub_matrix)

```c
int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
```

The sum of two matrices A = m × n and B = m × n of the same size is a matrix C = m × n = A + B of the same size whose elements are defined by the equations C(i,j) = A(i,j) + B(i,j).

The difference of two matrices A = m × n and B = m × n of the same size is a matrix C = m × n = A - B of the same size whose elements are defined by the equations C(i,j) = A(i,j) - B(i,j).


```
            1 2 3   1 0 0   2 2 3
С = A + B = 0 4 5 + 2 0 0 = 2 4 5
            0 0 6   3 4 1   3 4 7
```

### Matrix multiplication by scalar (mult_number). Multiplication of two matrices (mult_matrix)

```c
int s21_mult_number(matrix_t *A, double number, matrix_t *result);
int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
```

The product of the matrix A = m × n by the number λ is the matrix B = m × n = λ × A whose elements are defined by the equations B = λ × A(i,j).

```
                1 2 3   2 4 6   
B = 2 × A = 2 × 0 4 2 = 0 8 4 
                2 3 4   4 6 8   
```

The product of A = m × k by B = k × n is a matrix C = m × n = A × B of size m × n whose elements are defined by the equation C(i,j) = A(i,1) × B(1,j) + A(i,2) × B(2,j) + ... + A(i,k) × B(k,j).

```
            1 4    1 -1  1    9 11 17   
C = A × B = 2 5  × 2  3  4 = 12 13 22
            3 6              15 15 27
```

### Matrix transpose (transpose)

```c
int s21_transpose(matrix_t *A, matrix_t *result);
```

The transpose of matrix A is in switching its rows with its columns with their numbers retained

### Minor of matrix and matrix of algebraic complements (calc_complements)
```c
int s21_calc_complements(matrix_t *A, matrix_t *result);
```

Minor M(i,j) is a (n-1)-order determinant obtained by deleting out the i-th row and the j-th column from the matrix A.

The algebraic complement of a matrix element is the value of the minor multiplied by -1^(i+j).


### Matrix determinant

```c
int s21_determinant(matrix_t *A, double *result);
```

The determinant is a number that is associated to each square matrix and calculated from the elements using special formulas. \
Tip: The determinant can only be calculated for a square matrix.

The determinant of a matrix equals the sum of the products of elements of the row (column) and the corresponding algebraic complements.

### Inverse of the matrix (inverse_matrix)

```c
int s21_inverse_matrix(matrix_t *A, matrix_t *result);
```

A matrix A to the power of -1 is called the inverse of a square matrix A if the product of these matrices equals the identity matrix.

If the determinant of the matrix is zero, then it does not have an inverse.

The formula to calculate the inverse of matrix is $`A^{-1}=\frac{1} {|A|} × A_*^T`$
