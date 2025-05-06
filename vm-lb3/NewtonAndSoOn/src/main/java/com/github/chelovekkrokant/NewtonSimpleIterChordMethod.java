package com.github.chelovekkrokant;

public class NewtonSimpleIterChordMethod {

    static double function(double x){
        return (Math.log(x*x+4) / (x+6) - 9);
    }

    static double derivativeFunction(double x){
        return ((2 * x) / ((x + 6) * (x * x + 4))) - ((Math.log(x * x + 4)) / ((x + 6) * (x + 6)));
    }

    static double doubleDerivativeFunction(double x) {
        return (2 * Math.log(x * x + 4)) / Math.pow(x + 6, 3)
                - (2 * x) / (Math.pow(x + 6, 2) * (x * x + 4))
                + (2 * ((x + 6) * (x * x + 4) - x * (x * x + 2 * x * (x + 6) + 4)))
                / (Math.pow(x + 6, 2) * Math.pow(x * x + 4, 2));
    }

    static void getV (double x1, double x2, double eps) {
        double delta = Math.abs(x1 - x2);
        System.out.println("Basic eps   :   " + eps);
        System.out.println("Delta       :   x1 - x2         = " + delta);
        System.out.println("V           :   eps / delta     = " + eps / delta);
    }

    static void getVDelta (double x) {
        System.out.println("V(delta)    :   1 / deriv(f)    = " + 1 / Math.abs(derivativeFunction(x)));
    }

    static double newtonMethod(double initialGuess, double epsilon) {
        System.out.println("Метод Ньютона:");
        System.out.println("Приближение :   " + initialGuess);
        // Константы метода
        final double LOWER_M = -18.238091;
        final double UPPER_M = 172.072830;

        // Вычисление адаптированной точности
        final double adaptedepsilon = Math.sqrt(2 * LOWER_M * epsilon / UPPER_M);

        double currentX = initialGuess;
        double previousX = 0;
        int iterationCount = 0;

        do {
            iterationCount++;

            double functionValue = function(currentX);

            // Проверка на точное решение
            if (functionValue == 0.0) {
                processResult(currentX, previousX, epsilon, iterationCount);
                return currentX;
            }

            double derivativeValue = derivativeFunction(currentX);

            // Проверка на невозможность продолжения
            if (derivativeValue == 0.0) {
                return Double.NaN;
            }

            // Вычисление и применение шага Ньютона
            double newtonStep = functionValue / derivativeValue;
            previousX = currentX;
            currentX -= newtonStep;

        } while (Math.abs(currentX - previousX) >= adaptedepsilon);

        processResult(currentX, previousX, epsilon, iterationCount);
        return currentX;
    }

    // Вспомогательный метод для обработки результатов
    private static void processResult(double currentX, double previousX,
                                      double epsilon, int iterations) {
        getV(currentX, previousX, epsilon);
        getVDelta(currentX);
        System.out.println("Number of iterations : " + iterations);
    }

    static double phi (double x, double m, double M) {
        return x - 2 * function(x) / (m + M) ;
    }

    static double phiDerivative(double x, double m, double M) {
        return 1 - 2 / (m + M) * derivativeFunction(x);
    }

    static void getSimpleIterationVDelta(double x, double m, double M) {
        System.out.println("V delta = " + 1 / (1 - Math.abs(phiDerivative(x, m, M))));
    }

    static boolean checkPhi(double m, double M, double step) {
        for (double x = -5.65; x <= -5.55; x += step) {
            if (phiDerivative(x, m, M) >= 1) {
                return false;
            }
        }
        return true;
    }

    static double chordMethod(double leftBound, double rightBound, double epsilon) {
        System.out.println("Метод хорд:");

        // Вычисляем значения функции на границах интервала
        double leftValue = function(leftBound);
        double rightValue = function(rightBound);

        // Проверка корректности интервала
        if (leftValue * rightValue > 0.0) {
            System.out.println("Интервал задан неверно");
            return Double.NaN;
        }

        // Проверка на точное решение на границах
        if (leftValue == 0.0) return leftBound;
        if (rightValue == 0.0) return rightBound;

        // Инициализация переменных
        double currentApproximation = 0;
        double previousApproximation = 0;
        int iterationCount = 0;
        final double convergenceConstant = 1.8;

        double currentValue;
        do {
            iterationCount++;

            // Вычисление нового приближения
            currentApproximation = leftBound - (rightBound - leftBound) * leftValue
                    / (rightValue - leftValue);
            currentValue = function(currentApproximation);

            // Проверка на точное решение
            if (currentValue == 0.0) {
                processResult(currentApproximation, previousApproximation, epsilon, iterationCount);
                return currentApproximation;
            }

            // Обновление границ интервала
            if (currentValue * leftValue < 0.0) {
                previousApproximation = rightBound;
                rightBound = currentApproximation;
                rightValue = currentValue;
            } else {
                previousApproximation = leftBound;
                leftBound = currentApproximation;
                leftValue = currentValue;
            }

        } while (Math.abs(currentValue) >= epsilon);

        processResult(currentApproximation, previousApproximation, epsilon, iterationCount);
        return currentApproximation;
    }

    static double chordMethodGar(double leftBound, double rightBound, int n) {
        System.out.println("Метод хорд (Гарвик):");

        // Вычисляем значения функции на границах интервала
        double leftValue = function(leftBound);
        double rightValue = function(rightBound);

        // Проверка корректности интервала
        if (leftValue * rightValue > 0.0) {
            System.out.println("Интервал задан неверно");
            return Double.NaN;
        }

        // Проверка на точное решение на границах
        if (leftValue == 0.0) return leftBound;
        if (rightValue == 0.0) return rightBound;

        // Инициализация переменных
        double currentApproximation = 0;
        double previousApproximation = 0;
        int iterationCount = 0;
        final double convergenceConstant = 1.8;

        double currentValue;
        do {
            iterationCount++;

            // Вычисление нового приближения
            currentApproximation = leftBound - (rightBound - leftBound) * leftValue
                    / (rightValue - leftValue);
            currentValue = function(currentApproximation);

            // Проверка на точное решение
            if (currentValue == 0.0) {
                processResult(currentApproximation, previousApproximation, n, iterationCount);
                return currentApproximation;
            }

            // Обновление границ интервала
            if (currentValue * leftValue < 0.0) {
                previousApproximation = rightBound;
                rightBound = currentApproximation;
                rightValue = currentValue;
            } else {
                previousApproximation = leftBound;
                leftBound = currentApproximation;
                leftValue = currentValue;
            }

        } while (Math.abs(currentValue) >= n);

        processResult(currentApproximation, previousApproximation, n, iterationCount);
        return currentApproximation;
    }

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////
    static double simpleIterationsMethod(double initialGuess, double epsilon) {
        System.out.println("Метод простых итераций:");

        // Параметры метода
        final double UPPER_M = -18.238091;
        final double LOWER_M = -30.134383;

        // Проверка условий сходимости
        if (!checkConvergenceConditions(LOWER_M, UPPER_M, 0.001)) {
            return Double.NaN;
        }

        // Инициализация приближений
        double currentApproximation = phi(initialGuess, LOWER_M, UPPER_M);
        double nextApproximation = phi(currentApproximation, LOWER_M, UPPER_M);
        int iterationCount = 0;

        double q = (UPPER_M - LOWER_M) / (UPPER_M + LOWER_M);

        // Основной итерационный процесс
        while (Math.abs(currentApproximation - nextApproximation) > (1 - q) / q * epsilon) {
            iterationCount++;
            double previousApproximation = currentApproximation;
            currentApproximation = nextApproximation;
            nextApproximation = phi(currentApproximation, LOWER_M, UPPER_M);

        }

        // Обработка и вывод результатов
        processIterationResult(nextApproximation, currentApproximation, epsilon, iterationCount);
        getSimpleIterationVDelta(nextApproximation, LOWER_M, UPPER_M);

        return nextApproximation;
    }
    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////
    static double simpleIterationsMethodSave(double initialGuess, double epsilon) {
        System.out.println("Метод простых итераций:");

        // Параметры метода
        final double UPPER_M = -18.238091;
        final double LOWER_M = -30.134383;

        // Проверка условий сходимости
        if (!checkConvergenceConditions(LOWER_M, UPPER_M, 0.001)) {
            return Double.NaN;
        }

        // Инициализация приближений
        double currentApproximation = phi(initialGuess, LOWER_M, UPPER_M);
        double nextApproximation = phi(currentApproximation, LOWER_M, UPPER_M);
        int iterationCount = 0;
        double q = phiDerivative(currentApproximation, LOWER_M, UPPER_M);

        // Основной итерационный процесс
        while (Math.abs(currentApproximation - nextApproximation) > (1 - q) / q * epsilon) {
            iterationCount++;

            // Обновление приближений
            double previousApproximation = currentApproximation;
            currentApproximation = nextApproximation;
            nextApproximation = phi(currentApproximation, LOWER_M, UPPER_M);

            // Обновление коэффициента сходимости
            q = phiDerivative(nextApproximation, LOWER_M, UPPER_M);
        }

        // Обработка и вывод результатов
        processIterationResult(nextApproximation, currentApproximation, epsilon, iterationCount);
        getSimpleIterationVDelta(nextApproximation, LOWER_M, UPPER_M);

        return nextApproximation;
    }



    // Вспомогательные методы
    private static boolean checkConvergenceConditions(double lowerBound, double upperBound, double tolerance) {
        return checkPhi(lowerBound, upperBound, tolerance);
    }

    private static void processIterationResult(double current, double previous,
                                               double epsilon, int iterations) {
        getV(current, previous, epsilon);
        System.out.println("Number of iterations : " + (iterations + 1));
    }


    static void launchAlgsEps(){
        double left = -5.65, right = -5.55;
        double eps = 0.00001;

        double x0 = -5.64;
        double result = newtonMethod(x0, eps);
        System.out.println("Root x : " + result);

        System.out.println("################");

        result = simpleIterationsMethod(x0, eps);
        System.out.println("Root x : " + result);

        System.out.println("################");

        result = chordMethod(left, right, eps);
        System.out.println("Root x : " + result);

    }

    public static void main(String[] args) {
        launchAlgsEps();
    }

    ///////////////////////////////////////////////

    public static void findMinMaxOfDerivFunctions(double leftBorder, double righBorder){
        System.out.println("Первая и вторая производные монотонные на участке локализации.\n" +
                "Поэтому максимальное и минимальное значения будут находится на концах промежутка локализации.");
        double derLeftVal = derivativeFunction(leftBorder);
        double derRightVal = derivativeFunction(righBorder);
        double doubleDerLeftVal = doubleDerivativeFunction(leftBorder);
        double doubleDerRightVal = doubleDerivativeFunction(righBorder);
        if ( derRightVal > derLeftVal) {
            System.out.printf("Минимальное значение первой производной на [%f; %f] : %f\n",leftBorder,righBorder, derLeftVal);
            System.out.printf("Максимальное значение первой производной на [%f; %f] : %f\n",leftBorder,righBorder, derRightVal);
        } else {
            System.out.printf("Минимальное значение первой производной на [%f; %f] : %f\n",leftBorder,righBorder, derRightVal);
            System.out.printf("Максимальное значение первой производной на [%f; %f] : %f\n",leftBorder,righBorder, derLeftVal);
        }
        if ( doubleDerRightVal > doubleDerLeftVal) {
            System.out.printf("Минимальное значение второй производной на [%f; %f] : %f\n",leftBorder,righBorder, doubleDerLeftVal);
            System.out.printf("Максимальное значение второй производной на [%f; %f] : %f\n",leftBorder,righBorder, doubleDerRightVal);
        } else {
            System.out.printf("Минимальное значение второй производной на [%f; %f] : %f\n",leftBorder,righBorder, doubleDerRightVal);
            System.out.printf("Максимальное значение второй производной на [%f; %f] : %f\n",leftBorder,righBorder, doubleDerLeftVal);
        }
    }

}