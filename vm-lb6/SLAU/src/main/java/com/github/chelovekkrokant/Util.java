package com.github.chelovekkrokant;

import java.util.Random;

public class Util {

    // Исходная матрица А
    public static double[][] A = {
            {12.60920, 13.56706, 15.94644, 10.36991},
            {11.75549, 12.01742, 6.66735, 7.32291},
            {9.56322, 13.66453, 11.77007, 15.87784},
            {14.12295, 14.83661, 10.16751, 14.30176}
    };

    // Матрица А с ошибкой 0.05
    public static double[][] A_WITH_ERROR = {
            {12.64626, 13.59580, 15.90055, 10.41349},
            {11.77364, 12.06634, 6.64271, 7.37197},
            {9.60998, 13.63690, 11.76722, 15.85680},
            {14.10258, 14.83700, 10.13298, 14.31976}
    };

    // Исходная матрица А c элементами Матрицы Гильберта
    public static double[][] A_HILBERT = {
            {1.00000,   13.5670, 15.94644, 10.36991},
            {11.75549, 0.33333, 6.66735, 7.32291},
            {9.56322, 13.66453, 0.20000, 15.87784},
            {14.12295, 14.83661, 10.16751, 0.14286}
    };

    // Исходная матрица А c элементами Матрицы Гильберта с ошибкой 0.05
    public static double[][] A_HILBERT_WITH_ERROR = {
            {1.01895, 13.59372, 15.90824, 10.37957},
            {11.75064, 0.29938, 6.64799, 7.30878},
            {9.52879, 13.68825, 0.21474, 15.83027},
            {14.11733, 14.81182, 10.16147, 0.15056}
    };

    // Исходный вектор B
    public static double[] B = {
            9.25920, 13.98666, 10.58777, 10.27069
    };

    // Вектор B с ошибкой 0.05
    public static double[] B_WITH_ERROR = {
            9.23842, 14.0607, 10.55909, 10.30941
    };

    public static void main(String[] args){
        addErrorHilbertMatrix();
    }

    static void findDet(){
        System.out.println("Определитель сгенерированной матрицы = " + calculateDetFull(A));
    }

    static void generate(){
        System.out.println("Сгенерированная матрица:");
        printMatrix(generateMatrix());

        System.out.println("Сгенерированный вектор:");
        printVector(generateVector());
    }

    static void addErrorVector(){
        System.out.println("Базовый вектора:");
        printVector(B);
        System.out.println("Добавлена ошибка в вектор свободных членов:");
        printVector(addErrorToVector(B));
    }

    static void addErrorMatrix(){
        System.out.println("Базовая матрица:");
        printMatrix(A);
        System.out.println("Добавлена ошибка в значения элементов матрицы:");
        printMatrix(addErrorToMatrix(A));
    }

    static void addErrorHilbertMatrix(){
        System.out.println("Базовая матрица:");
        printMatrix(A_HILBERT);
        System.out.println("Добавлена ошибка в значения элементов матрицы с элементами матрицы Гильберта:");
        printMatrix(addErrorToMatrix(A_HILBERT));
    }

    static void replaceDiagonalValues(){
        System.out.println("Базовая матрица:");
        printMatrix(A);
        System.out.println("Элементы главной диагонали заменены на элементы матрицы Гильберта:");
        printMatrix(replaceDiagonalValuesToHilberts(A));
    }

    static double[][] replaceDiagonalValuesToHilberts(double[][] matrix) {
        int n = matrix.length; double[][] result = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i == j) {
                    result[i][j] = 1.0 / (i + j + 1);
                } else {
                    result[i][j] = matrix[i][j];
                }
            }
        }
        return result;
    }

    static double[][] addErrorToMatrix(double[][] matrix) {
        Random rand = new Random();
        double[][] result = deepCopy(matrix);
        for (int i = 0; i < result.length; i++) {
            for (int j = 0; j < result[0].length; j++) {
                result[i][j] += 0.1 * (rand.nextDouble() - 0.5);
            }
        }
        return result;
    }

    static double[] addErrorToVector(double[] vector) {
        Random rand = new Random();
        double[] result = vector.clone();
        for (int i = 0; i < result.length; i++) {
            result[i] += 0.1 * (rand.nextDouble() - 0.5);
        }
        return result;
    }

    static double[][] generateMatrix() {
        Random rand = new Random();
        double[][] matrix = new double[4][4];
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                matrix[i][j] = rand.nextDouble() * 10 + 6;
            }
        }
        return matrix;
    }

    static double[] generateVector() {
        Random rand = new Random();
        double[] vector = new double[4];
        for (int i = 0; i < 4; i++) {
            vector[i] = rand.nextDouble() * 10 + 6;
        }
        return vector;
    }

    static void printMatrix(double[][] matrix) {
        for (double[] row : matrix) {
            for (double value : row) {
                System.out.printf("%10.5f ", value);
            }
            System.out.println();
        }
    }

    static void printVector(double[] vector) {
        for (double v : vector) {
            System.out.printf("%10.5f ", v);
        }
        System.out.println();
    }

    static double[][] deepCopy(double[][] matrix) {
        double[][] copy = new double[matrix.length][matrix[0].length];
        for (int i = 0; i < matrix.length; i++) {
            System.arraycopy(matrix[i], 0, copy[i], 0, matrix[0].length);
        }
        return copy;
    }

    static double calculateDetFull(double[][] matrix) {
        if (matrix.length == 2) {
            return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
        }
        if (matrix.length == 4) {
            double det = 0;
            for (int i = 0; i < 4; i++) {
                double[][] minor = new double[3][3];
                for (int j = 1; j < 4; j++) {
                    int minorCol = 0;
                    for (int k = 0; k < 4; k++) {
                        if (k == i) continue;
                        minor[j-1][minorCol++] = matrix[j][k];
                    }
                }
                det += (i % 2 == 0 ? 1 : -1) * matrix[0][i] * calculateDetMinor(minor);
            }
            return det;
        }
        return 0;
    }

    static double calculateDetMinor(double[][] minor) {
        return minor[0][0] * (minor[1][1]*minor[2][2] - minor[1][2]*minor[2][1])
                - minor[0][1] * (minor[1][0]*minor[2][2] - minor[1][2]*minor[2][0])
                + minor[0][2] * (minor[1][0]*minor[2][1] - minor[1][1]*minor[2][0]);
    }

}
