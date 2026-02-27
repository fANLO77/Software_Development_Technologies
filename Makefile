CC = gcc
CFLAGS = -O2 -Wall
LDFLAGS = -lm

# Указываем пути поиска заголовочных файлов
# Добавляем корень (.), папку OpenBLAS и папку fakeblas
INCLUDES = -I. -I./OpenBLAS -I./fakeblas

OPENBLAS_LIB = ./OpenBLAS/libopenblas.a
TEST_SRC = tests/test_level2_full.c

all: test_openblas test_fakeblas

# ================= OpenBLAS =================

test_openblas: $(TEST_SRC)
	$(CC) $(CFLAGS) $(INCLUDES) $(TEST_SRC) $(OPENBLAS_LIB) $(LDFLAGS) -o test_openblas

run_openblas: test_openblas
	./test_openblas

# ================= FakeBLAS =================

fakeblas/libfakeblas.a: fakeblas/fake_cblas.c
	$(CC) $(CFLAGS) $(INCLUDES) -c fakeblas/fake_cblas.c -o fakeblas/fake_cblas.o
	ar rcs fakeblas/libfakeblas.a fakeblas/fake_cblas.o

test_fakeblas: $(TEST_SRC) fakeblas/libfakeblas.a
	$(CC) $(CFLAGS) $(INCLUDES) $(TEST_SRC) fakeblas/libfakeblas.a $(LDFLAGS) -o test_fakeblas

run_fakeblas: test_fakeblas
	./test_fakeblas

clean:
	rm -f test_openblas test_fakeblas fakeblas/*.o fakeblas/*.a