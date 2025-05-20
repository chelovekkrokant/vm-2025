package com.github.chelovekkrokant;

import static com.github.chelovekkrokant.Util.*;

public class SLAU {
    
    public static void main(String[] args) {
        System.out.println("\n\nЗАДАНИЕ 1.\nРАБОТА С БАЗОВЫМИ ЗНАЧЕНИЯМИ");
        printMatrix(A);
        printVector(B);
        // Проверяем невырожденность
        System.out.printf("Определитель матрицы A: %f\n", calculateDetFull(A));
        if (Math.abs(calculateDetFull(A)) < 1e-8) {
            System.out.println("Матрица вырождена. Решение невозможно.");
            return;
        }
        // Решение методом обратной матрицы
        double[] X = solveCramer(A, B);
        System.out.println("Решение X1:");
        Util.printVector(X);
        // Оцениваем числа обусловленности
        evaluateConditionNumbers(A, B, X);
        /// ///////////////////////////////////////////////////////////////////


        System.out.println("\n\nЗАДАНИЕ 2.\nРАБОТА С ВЕКТОРОМ В С ОШИБКАМИ");
        printMatrix(A);
        printVector(B_WITH_ERROR);
        double[] X_star = solveCramer(A, B_WITH_ERROR);
        System.out.println("Решение X2:");
        Util.printVector(X_star);
        evaluateConditionNumbers(A, B_WITH_ERROR, X_star);
        calculateAllErrorsAndCondition(
                A, A,
                B, B_WITH_ERROR,
                X, X_star, 1);
        /// ///////////////////////////////////////////////////////////////////


        System.out.println("\n\nЗАДАНИЕ 3.\nРАБОТА С МАТРИЦЕЙ А С ОШИБКАМИ");
        printMatrix(A_WITH_ERROR);
        printVector(B);
        X_star = solveCramer(A_WITH_ERROR, B);
        System.out.println("Решение X3:");
        Util.printVector(X_star);
        evaluateConditionNumbers(A_WITH_ERROR, B, X_star);
        calculateAllErrorsAndCondition(
                A, A_WITH_ERROR,
                B, B,
                X, X_star, 2);
        /// ///////////////////////////////////////////////////////////////////


        System.out.println("\n\nЗАДАНИЕ 4.\nРАБОТА С МАТРИЦЕЙ А И ВЕКТОРОМ В С ОШИБКАМИ");
        printMatrix(A_WITH_ERROR);
        printVector(B_WITH_ERROR);
        X_star = solveCramer(A_WITH_ERROR, B_WITH_ERROR);
        System.out.println("Решение X4:");
        Util.printVector(X_star);
        evaluateConditionNumbers(A_WITH_ERROR, B_WITH_ERROR, X_star);
        calculateAllErrorsAndCondition(
                A, A_WITH_ERROR,
                B, B_WITH_ERROR,
                X, X_star, 3);
        /// ///////////////////////////////////////////////////////////////////


        System.out.println("\n\nЗАДАНИЕ 5.\nРАБОТА С МАТРИЦЕЙ (ГИЛЬБЕРТ)");
        printMatrix(A_HILBERT);
        printVector(B);
        System.out.printf("Определитель матрицы Mod: %f\n", calculateDetFull(A_HILBERT));
        if (Math.abs(calculateDetFull(A_HILBERT)) < 1e-8) {
            System.out.println("Модифицированная матрица вырождена. Решение невозможно.");
            return;
        }
        X = solveCramer(A_HILBERT, B);
        System.out.println("Решение X5:");
        Util.printVector(X);
        evaluateConditionNumbers(A_HILBERT, B, X);
        /// ///////////////////////////////////////////////////////////////////


        System.out.println("\n\nЗАДАНИЕ 6.\nРАБОТА С МАТРИЦЕЙ (ГИЛЬБЕРТ) И ОШИБКАМИ ВЕКТОРА В");
        printMatrix(A_HILBERT);
        printVector(B_WITH_ERROR);
        X_star = solveCramer(A_HILBERT, B_WITH_ERROR);
        System.out.println("Решение X6:");
        Util.printVector(X_star);
        evaluateConditionNumbers(A_HILBERT, B_WITH_ERROR, X_star);
        calculateAllErrorsAndCondition(
                A_HILBERT, A_HILBERT,
                B, B_WITH_ERROR,
                X, X_star, 1);
        /// ///////////////////////////////////////////////////////////////////


        System.out.println("\n\nЗАДАНИЕ 7.\nРАБОТА С ОШИБКАМИ В МАТРИЦЕ (ГИЛЬБЕРТ)");
        printMatrix(A_HILBERT_WITH_ERROR);
        printVector(B);
        X_star = solveCramer(A_HILBERT_WITH_ERROR, B);
        System.out.println("Решение X7:");
        Util.printVector(X_star);
        evaluateConditionNumbers(A_HILBERT_WITH_ERROR, B, X_star);
        calculateAllErrorsAndCondition(
                A_HILBERT, A_HILBERT_WITH_ERROR,
                B, B,
                X, X_star, 2);
        /// ///////////////////////////////////////////////////////////////////


        System.out.println("\n\nЗАДАНИЕ 8.\nРАБОТА С ОШИБКАМИ В ВЕКТОРЕ В И МАТРИЦЕ (ГИЛЬБЕРТ)");
        printMatrix(A_HILBERT_WITH_ERROR);
        printVector(B_WITH_ERROR);
        X_star = solveCramer(A_HILBERT_WITH_ERROR, B_WITH_ERROR);
        System.out.println("Решение X8:");
        Util.printVector(X_star);
        evaluateConditionNumbers(A_HILBERT_WITH_ERROR, B_WITH_ERROR, X_star);
        calculateAllErrorsAndCondition(
                A_HILBERT, A_HILBERT_WITH_ERROR,
                B, B_WITH_ERROR,
                X, X_star, 3);
    }

    /**
     * Решает систему линейных уравнений методом Крамера
     *
     * @param coefficientMatrix матрица коэффициентов системы (n x n)
     * @param rightHandSide вектор правых частей системы (длина n)
     * @return вектор решения системы
     * @throws IllegalArgumentException если матрица вырождена (определитель = 0)
     */
    static double[] solveCramer(double[][] coefficientMatrix, double[] rightHandSide) {
        // Проверка размеров входных данных
        int systemSize = coefficientMatrix.length;
        if (systemSize != rightHandSide.length) {
            throw new IllegalArgumentException(
                    "Размер матрицы коэффициентов (" + systemSize + "x" + systemSize +
                            ") не соответствует размеру вектора правых частей (" + rightHandSide.length + ")"
            );
        }

        // Вычисляем определитель основной матрицы
        double mainDeterminant = calculateDetFull(coefficientMatrix);

        // Проверка на вырожденность матрицы
        if (Math.abs(mainDeterminant) < 1e-10) {
            throw new IllegalArgumentException(
                    "Определитель матрицы равен " + mainDeterminant +
                            ". Система не имеет единственного решения."
            );
        }

        // Инициализация вектора результатов
        double[] solution = new double[systemSize];

        // Вычисление каждого компонента решения по правилу Крамера
        for (int variableIndex = 0; variableIndex < systemSize; variableIndex++) {
            // Создаем модифицированную матрицу, заменяя columnIndex-ый столбец на вектор правых частей
            double[][] A_HILBERT = createModifiedMatrix(coefficientMatrix, rightHandSide, variableIndex);

            // Вычисляем определитель модифицированной матрицы
            double modifiedDeterminant = calculateDetFull(A_HILBERT);

            // Вычисляем значение текущей переменной
            solution[variableIndex] = modifiedDeterminant / mainDeterminant;
        }

        return solution;
    }

    /**
     * Создает модифицированную матрицу для метода Крамера, заменяя указанный столбец на вектор правых частей
     *
     * @param originalMatrix исходная матрица коэффициентов
     * @param rightHandSide вектор правых частей
     * @param columnToReplace индекс столбца для замены
     * @return модифицированная матрица
     */
    private static double[][] createModifiedMatrix(double[][] originalMatrix, double[] rightHandSide, int columnToReplace) {
        int size = originalMatrix.length;
        double[][] A_HILBERT = new double[size][size];

        // Копируем исходную матрицу, заменяя указанный столбец
        for (int row = 0; row < size; row++) {
            System.arraycopy(originalMatrix[row], 0, A_HILBERT[row], 0, size);
            A_HILBERT[row][columnToReplace] = rightHandSide[row];
        }

        return A_HILBERT;
    }
    

    static double[][] invertMatrix(double[][] matrix) {
        int n = matrix.length;
        double[][] a = Util.deepCopy(matrix);
        double[][] inv = new double[n][n];

        for (int i = 0; i < n; i++) inv[i][i] = 1.0;

        for (int i = 0; i < n; i++) {
            double diag = a[i][i];
            for (int j = 0; j < n; j++) {
                a[i][j] /= diag;
                inv[i][j] /= diag;
            }
            for (int k = 0; k < n; k++) {
                if (k != i) {
                    double factor = a[k][i];
                    for (int j = 0; j < n; j++) {
                        a[k][j] -= factor * a[i][j];
                        inv[k][j] -= factor * inv[i][j];
                    }
                }
            }
        }
        return inv;
    }

    static double calculateColumnNorm(double[][] matrix) {
        double max = 0;
        for (int j = 0; j < matrix[0].length; j++) {
            double sum = 0;
            for (int i = 0; i < matrix.length; i++) {
                sum += Math.abs(matrix[i][j]);
            }
            if (sum > max) max = sum;
        }
        return max;
    }

    static double calculateVectorNorm(double[] vector) {
        double sum = 0;
        for (double v : vector) {
            sum += Math.abs(v);
        }
        return sum;
    }

    static double matrixDifferenceNorm(double[][] A, double[][] A_star) {
        double max = 0;
        for (int j = 0; j < A[0].length; j++) { // идём по столбцам
            double sum = 0;
            for (int i = 0; i < A.length; i++) { // суммируем по строкам
                sum += Math.abs(A[i][j] - A_star[i][j]);
            }
            if (sum > max) max = sum;
        }
        return max;
    }

    static double vectorDifferenceNorm(double[] B, double[] B_star) {
        double sum = 0;
        for (int i = 0; i < B.length; i++) {
            sum += Math.abs(B[i] - B_star[i]);
        }
        return sum;
    }

    static void evaluateConditionNumbers(double[][] A, double[] B, double[] X) {
        double[][] invertedA = invertMatrix(A);
        double columnNormA = calculateColumnNorm(A);
        double columnNormInvertedA = calculateColumnNorm(invertedA);
        double conditionNumber = columnNormA * columnNormInvertedA;
        System.out.printf("Абсолютное число обусловленности ν_Δ       : %.6f\n", columnNormInvertedA);
        System.out.printf("Естественное число обусловленности ν_δ     : %.6f\n", columnNormInvertedA * (calculateVectorNorm(B)/calculateVectorNorm(X)));
        System.out.printf("Стандартное число обусловленности          : %.6f\n", conditionNumber);
    }

    public static double evaluateConditionNumber(double deltaX, double deltaA, double deltaB, int type) {
        if (type == 3){
            return deltaX / (deltaA + deltaB);
        }
        else if (type == 2) {
            return deltaX / (deltaA);
        } else {
            return deltaX / (deltaB);
        }
    }

    public static void calculateAllErrorsAndCondition(
            double[][] A, double[][] A_star,
            double[] B, double[] B_star,
            double[] X, double[] X_star, int type) {

        double[][] invertedA = invertMatrix(A_star);
        double columnNormA = calculateColumnNorm(A_star);
        double columnNormInvertedA = calculateColumnNorm(invertedA);
        double conditionNumber = columnNormA * columnNormInvertedA;
        double normB = calculateVectorNorm(B);
        double normX = calculateVectorNorm(X);

        double matrixAbsoluteError = matrixDifferenceNorm(A, A_star);
        double vectorAbsoluteError = vectorDifferenceNorm(B, B_star);
        double solutionAbsoluteError = vectorDifferenceNorm(X, X_star);

        double matrixRelativeError = matrixAbsoluteError / columnNormA;
        double vectorRelativeError = vectorAbsoluteError / normB;
        double solutionRelativeError = solutionAbsoluteError / normX;

        double conditionEstimate = evaluateConditionNumber(solutionRelativeError, matrixRelativeError, vectorRelativeError, type);

        System.out.println("\nАнализ погрешностей:");
        System.out.printf("Абсолютная погрешность матрицы       ∆(A*) : %.6f\n", matrixAbsoluteError);
        System.out.printf("Относительная погрешность матрицы    δ(A*) : %.6f\n", matrixRelativeError);
        System.out.printf("Абсолютная погрешность вектора В     ∆(B*) : %.6f\n", vectorAbsoluteError);
        System.out.printf("Относительная погрешность вектора В  δ(B*) : %.6f\n", vectorRelativeError);
        System.out.printf("Абсолютная погрешность вектора Х     ∆(X*) : %.6f\n", solutionAbsoluteError);
        System.out.printf("Относительная погрешность вектора Х  δ(X*) : %.6f\n", solutionRelativeError);
        System.out.printf("Оценка числа обусловленности               : cond(A) > %.8f\n", conditionEstimate);
        System.out.printf("Фактическое число обусловленности          : %.8f\n", conditionNumber);
    }

}