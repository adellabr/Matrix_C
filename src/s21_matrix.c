#include "s21_matrix.h"

int s21_create_matrix(int rows, int columns, matrix_t *result) {
  if (result == NULL || rows <= 0 || columns <= 0) {
    return 1;
  }
  int error = 0;
  if ((result->matrix = (double **)calloc(rows, sizeof(double *)))) {
    for (int i = 0; i < rows && !error; i++) {
      result->matrix[i] = (double *)calloc(columns, sizeof(double));
      if (result->matrix[i] == NULL) error = 1;
    }
    if (!error) {
      result->rows = rows;
      result->columns = columns;
    }
  } else
    error = 1;
  return error;
}

void s21_remove_matrix(matrix_t *A) {
  if (A == NULL) return;
  if (A != NULL && A->rows != 0) {
    for (int i = 0; i < A->rows; i++) {
      if (A->matrix[i] != NULL) free(A->matrix[i]);
    }
    if (A->matrix != NULL) free(A->matrix);
    A->matrix = NULL;
    A->rows = 0;
    A->columns = 0;
  }
}

int s21_eq_matrix(matrix_t *A, matrix_t *B) {
  if (A == NULL || B == NULL || A->rows <= 0 || A->columns <= 0 ||
      B->rows <= 0 || B->columns <= 0 || A->rows != B->rows ||
      A->columns != B->columns)
    return 0;
  int diff = 1;
  for (int i = 0; i < A->rows && diff; i++) {
    for (int j = 0; j < A->columns && diff; j++) {
      if (fabs(A->matrix[i][j] - B->matrix[i][j]) > s21_EPS) diff = 0;
    }
  }
  return diff;
}

int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  if (B == NULL || result == NULL || A == NULL || A->rows <= 0 ||
      A->columns <= 0 || B->rows <= 0 || B->columns <= 0)
    return 1;
  if (A->rows != B->rows || A->columns != B->columns) return 2;
  int error = 0;
  error = s21_create_matrix(A->rows, A->columns, result);
  for (int i = 0; i < A->rows; i++) {
    for (int j = 0; j < A->columns; j++) {
      result->matrix[i][j] = A->matrix[i][j] + B->matrix[i][j];
    }
  }
  return error;
}

int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  if (A == NULL || B == NULL || result == NULL || A->rows <= 0 ||
      A->columns <= 0 || B->rows <= 0 || B->columns <= 0)
    return 1;
  if (A->rows != B->rows || A->columns != B->columns) return 2;
  int error = 0;
  error = s21_create_matrix(A->rows, A->columns, result);
  for (int i = 0; i < A->rows; i++) {
    for (int j = 0; j < A->columns; j++) {
      result->matrix[i][j] = A->matrix[i][j] - B->matrix[i][j];
    }
  }
  return error;
}

int s21_mult_number(matrix_t *A, double number, matrix_t *result) {
  if (A == NULL || result == NULL || A->rows <= 0 || A->columns <= 0) return 1;
  int error = 0;
  error = s21_create_matrix(A->rows, A->columns, result);
  for (int i = 0; i < A->rows; i++) {
    for (int j = 0; j < A->columns; j++) {
      result->matrix[i][j] = A->matrix[i][j] * number;
    }
  }
  return error;
}

int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  if (A == NULL || B == NULL || result == NULL || A->rows <= 0 ||
      A->columns <= 0 || B->rows <= 0 || B->columns <= 0)
    return 1;
  if (A->rows != B->columns || A->columns != B->rows) return 2;
  int error = 0;
  error = s21_create_matrix(A->rows, B->columns, result);
  if (!error) {
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < B->columns; j++) {
        for (int k = 0; k < A->columns; k++) {
          result->matrix[i][j] += A->matrix[i][k] * B->matrix[k][j];
        }
      }
    }
  }
  return error;
}

int s21_transpose(matrix_t *A, matrix_t *result) {
  if (A == NULL || result == NULL || A->rows <= 0 || A->columns <= 0) return 1;
  int error = 0;
  error = s21_create_matrix(A->columns, A->rows, result);
  for (int i = 0; i < A->rows; i++) {
    for (int j = 0; j < A->columns; j++) {
      result->matrix[j][i] = A->matrix[i][j];
    }
  }
  return error;
}

void calc_minor(matrix_t *A, matrix_t * compl, int size, int i, int j) {
  for (int row = 0, row1 = 0; row < size; row++) {
    for (int col = 0, col1 = 0; col < size; col++) {
      if (row != i && col != j) {
        compl ->matrix[row1][col1++] = A->matrix[row][col];
        if (col1 == size - 1) {
          col1 = 0;
          row1++;
        }
      }
    }
  }
}

double determinant(matrix_t *A, int size) {
  double result = 0;
  if (size == 1)
    return A->matrix[0][0];
  else {
    matrix_t compl = {0};
    s21_create_matrix(size, size, &compl );
    for (int j = 0; j < size; j++) {
      calc_minor(A, &compl, size, 0, j);
      result += pow(-1, j) * A->matrix[0][j] * determinant(&compl, size - 1);
    }

    s21_remove_matrix(&compl );
  }
  return result;
}

int s21_determinant(matrix_t *A, double *result) {
  if (A == NULL || A->rows <= 0 || A->columns <= 0)
    return 1;
  else if (A->rows != A->columns)
    return 2;
  else {
    *result = determinant(A, A->rows);
  }
  return 0;
}

int s21_calc_complements(matrix_t *A, matrix_t *result) {
  if (A == NULL || result == NULL || A->rows <= 0 || A->columns <= 0) return 1;
  if (A->rows != A->columns) return 2;
  double res = 0;
  int error = 0;
  error = s21_create_matrix(A->rows, A->columns, result);
  if (A->rows == 1)
    result->matrix[0][0] = 1;
  else if (!error) {
    matrix_t compl = {0};
    error = s21_create_matrix(A->rows - 1, A->rows - 1, &compl );
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        calc_minor(A, &compl, A->rows, i, j);
        s21_determinant(&compl, &res);
        result->matrix[i][j] = pow(-1, (i + j)) * res;
      }
    }
    s21_remove_matrix(&compl );
  }
  return error;
}

int s21_inverse_matrix(matrix_t *A, matrix_t *result) {
  if (A == NULL || result == NULL || A->rows <= 0 || A->columns <= 0) return 1;
  double det = 0;
  s21_determinant(A, &det);
  if (A->rows != A->columns || det == 0) return 2;
  s21_create_matrix(A->rows, A->columns, result);
  int error = 0;
  if (A->rows == 1)
    result->matrix[0][0] = 1 / A->matrix[0][0];
  else {
    matrix_t trans = {0};
    matrix_t minors = {0};
    error = s21_create_matrix(A->columns, A->rows, &minors);
    if (!error) error = s21_create_matrix(A->columns, A->rows, &trans);
    if (!error) error = s21_calc_complements(A, &minors);
    if (!error) error = s21_transpose(&minors, &trans);
    if (!error) error = s21_mult_number(&trans, (double)1 / det, result);
    if (!error) s21_remove_matrix(&trans);
  }
  return error;
}