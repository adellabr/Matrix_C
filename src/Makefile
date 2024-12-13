CC = gcc
CFLAGS = -Wall -Werror -Wextra -std=c11
GCOVFLAGS = -Wall -Werror -Wextra --coverage
SRC = s21_matrix.c
OBJS = $(SRC:.c=.o)

OS := $(shell uname -s)

ifeq ($(OS),Linux)
  LIBS =  -pthread -lcheck -lm -lrt -lsubunit
else
  LIBS = -lcheck -lm
endif


all: gcov_report s21_matrix.a test

gcov_report: clean s21_matrix_test.c $(SRC)
	$(CC) $(GCOVFLAGS) -c $(SRC) 
	ar rcs s21_matrix.a $(OBJS)
	$(CC) $(GCOVFLAGS) -c s21_matrix_test.c -o s21_matrix_test.o
	$(CC) $(GCOVFLAGS) s21_matrix_test.o s21_matrix.a -o gcov_unit_test $(LIBS)
	./gcov_unit_test
	lcov -t "test" -o test.info --no-external -c -d .
	genhtml -o report/ test.info
	open ./report/index.html
	make clean_without_report

test: s21_matrix.a
	$(CC) $(CFLAGS) s21_matrix_test.c $< -o unit_test $(LIBS)
	./unit_test

matrix_tests.o: s21_matrix_test.c
	$(CC) $(CFLAGS) -c s21_matrix_test.c -o $@

s21_matrix.a: $(OBJS)
	ar rcs s21_matrix.a $(OBJS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clang:
	clang-format -i *.c *.h 

clean_without_report: 
	rm -f  *.o *.gcno *.gcov *.gcda *.info *.a main_test unit_test gcov_unit_test

clean: clean_without_report
	rm -rf report test_log/*.log check