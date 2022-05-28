package com.zivlazarov.linearalgebra;

import java.util.Random;

public class TwoDimMatrice {

    private final int numRows;
    private final int numCols;

    private final int[] shape;

    private final double[][] matrice;

    private static final Random random = new Random();

    public TwoDimMatrice(int numRows, int numCols) {
        this.numRows = numRows;
        this.numCols = numCols;
        shape = new int[] {numRows, numCols};

        matrice = new double[numRows][numCols];
    }

    public TwoDimMatrice(double[][] newMatrice) {
        matrice = newMatrice;
        numRows = newMatrice.length;
        numCols = newMatrice[0].length;
        shape = new int[] {numRows, numCols};
    }

    public boolean forwardElimination(double[][] mat, int row, int col) {
        int M = mat.length;
        int N = mat[0].length;

        // matrice isn't invertible
        if (M > N || M < N) return false;

        if (col >= N && row >= M) {
            return backElimination(mat, row - 1, col - 1);
        }

        if (mat[row][col] == 0) {
            for (int r = row + 1; r < M; r++) {
                if (mat[r][col] != 0) {
                    int counter = 1;
                    for (int c = 0; c < N; c++) {
                        if (mat[r][c] == 0) {
                            counter++;
                            continue;
                        }
                        double temp = mat[row][c];
                        mat[row][c] = mat[r][c];
                        mat[r][c] = temp;
                    }
                    if (counter == N) break;
                }
            }
        }

        if (mat[row][col] != 1 && mat[row][col] != 0) {
            for (int c = 0; c < N; c++) {
                mat[row][c] *= 1.0 / mat[row][col];
            }
        }

        for (int r = row + 1; r < M; r++) {
            double fact = mat[r][col];
            for (int c = 0; c < N; c++) {
                mat[r][c] -= mat[row][c] * fact;
            }
        }
        return forwardElimination(mat, row + 1, col + 1);
    }

    public boolean backElimination(double[][] mat, int row, int col) {
        int M = mat.length;
        int N = mat[0].length;

        if (col < 0 || row < 0) return true;

        if (mat[row][col] == 0) {
            for (int r = row + 1; r < M; r++) {
                if (mat[r][col] != 0) {
                    int counter = 1;
                    for (int c = 0; c < N; c++) {
                        if (mat[r][c] == 0) {
                            counter++;
                            continue;
                        }
                        double temp = mat[row][c];
                        mat[row][c] = mat[r][c];
                        mat[r][c] = temp;
                    }
                    if (counter == N) break;
                }
            }
        }

        if (mat[row][col] != 1 && mat[row][col] != 0) {
            for (int c = 0; c < N; c++) {
                mat[row][c] *= 1.0 / mat[row][col];
            }
        }

        for (int r = row - 1; r >= 0; r--) {
            double fact = mat[r][col];
            for (int c = 0; c < N; c++) {
                mat[r][c] -= mat[row][c] * fact;
            }
        }

        return backElimination(mat, row - 1, col - 1);
    }

    public TwoDimMatrice add(TwoDimMatrice other) {
        if (shape[0] != other.shape[0] || shape[1] != other.shape[1]) {
            throw new RuntimeException("Cannot add matrices with different shapes!");
        }
        double[][] newMatrice = new double[numRows][numCols];
        for (int r = 0; r < numRows; r++) {
            for (int c = 0; c < numCols; c++) {
                newMatrice[r][c] = matrice[r][c] + other.matrice[r][c];
            }
        }
        return new TwoDimMatrice(newMatrice);
    }

    public TwoDimMatrice sub(TwoDimMatrice other) {
        if (shape[0] != other.shape[0] || shape[1] != other.shape[1]) {
            throw new RuntimeException("Cannot add matrices with different shapes!");
        }
        double[][] newMatrice = new double[numRows][numCols];
        for (int r = 0; r < numRows; r++) {
            for (int c = 0; c < numCols; c++) {
                newMatrice[r][c] = matrice[r][c] - other.matrice[r][c];
            }
        }
        return new TwoDimMatrice(newMatrice);
    }

    public TwoDimMatrice multiply(TwoDimMatrice other) {
        if (shape[1] != other.shape[0]) {
            throw new RuntimeException("Cannot multiply matrice with " + shape[1]
                    + " columns and a matrice with " + other.shape[0] + " rows!");
        }
        double[][] newMatrice = new double[numRows][other.shape[1]];
        for (int r = 0; r < numRows; r++) {
            for (int c = 0; c < other.shape[1]; c++) {
                for (int col = 0; col < numCols; col++) {
                    newMatrice[r][c] += matrice[r][col] * other.matrice[col][c];
                }
            }
        }
        return new TwoDimMatrice(newMatrice);
    }

    public TwoDimMatrice power(int n) {
        TwoDimMatrice newMat = new TwoDimMatrice(matrice);
        for (int i = 0; i < n-1; i++) {
            newMat = newMat.multiply(this);
        }
        return newMat;
    }

    public TwoDimMatrice scalarMultiply(double x) {
        TwoDimMatrice newMat = new TwoDimMatrice(shape[0], shape[1]);
        for (int r = 0; r < shape[0]; r++) {
            for (int c = 0; c < shape[1]; c++) {
                newMat.matrice[r][c] = matrice[r][c] * x;
            }
        }
        return newMat;
    }

    public double determinant() {
        return determinant(matrice, 0);
    }

    public double determinant(double[][] mat, double sum) {
        if (mat.length != mat[0].length) throw new RuntimeException("Determinant isn't defined for non squared matrices!");
        int n = mat.length;
        if (n <= 2) {
            sum = mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];
            return sum;
        }

        for (int i = 0; i < n; i++) {
            double[][] tempMinor = new double[n-1][n-1];
            for (int r = 1; r < n; r++) {
                for (int c = 0; c < n; c++) {
                    if (c == i) continue;
                    if (c < i) tempMinor[r-1][c] = mat[r][c];
                    else tempMinor[r-1][c-1] = mat[r][c];
                }
            }
            sum += Math.pow(-1, i+2) * mat[0][i] * determinant(tempMinor, 0);
        }
        return sum;
    }

    public TwoDimMatrice transpose() {
        double[][] newMatrice = new double[numCols][numRows];
        for (int r = 0; r < numRows; r++) {
            for (int c = 0; c < numCols; c++) {
                newMatrice[c][r] = matrice[r][c];
            }
        }
        return new TwoDimMatrice(newMatrice);
    }

    public void print() {
        for (int i = 0; i < numRows; i++) {
            System.out.println();
            for (int j = 0; j < numCols; j++) {
                System.out.print(matrice[i][j] + " ");
            }
        }
    }

    public static TwoDimMatrice I(int n) {
        double[][] m = new double[n][n];

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i == j) m[i][j] = 1;
                else m[i][j] = 0;
            }
        }
        return new TwoDimMatrice(m);
    }

    public static TwoDimMatrice generate(int rows, int cols, int max) {
        double[][] m = new double[rows][cols];
        for (int r = 0; r < rows; r++) {
            for (int c = 0; c < cols; c++) {
                m[r][c] = Math.round(random.nextDouble() * max);
            }
        }
        return new TwoDimMatrice(m);
    }

    public int[] getShape() {
        return shape;
    }

    public double[][] getMatrice() {
        return matrice;
    }
}
