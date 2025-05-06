package com.github.chelovekkrokant;

public class GarwickNewtonSimpleIterChordMethod {

    static double function(double x){
        return (Math.log(x*x+4) / (x+6) - 9);
    }

    static double derivativeFunction(double x){
        return ((2 * x) / ((x + 6) * (x * x + 4))) - ((Math.log(x * x + 4)) / ((x + 6) * (x + 6)));
    }

    static double phi (double x, double m, double M) {
        return x - 2 / (m + M) * function(x);
    }

    static double chordMethod(double a, double b, int maxIterations) {
        System.out.println("Метод хорд:");

        double leftValue = function(a);
        double rightValue = function(b);
        double currentX = 0, currentY, prevX = a, prevPrevX = b;
        int iterationCount = 0;

        if (leftValue * rightValue > 0.0) {
            System.out.println("Интервал задан неверно");
            return Double.NaN;
        }
        if (leftValue == 0.0) return leftValue;
        if (rightValue == 0.0) return rightValue;

        while (iterationCount < maxIterations) {
            currentX = a - (b - a) * leftValue / (rightValue - leftValue);
            currentY = function(currentX);

            if (currentY == 0.0) {
                return currentX;
            }

            if (iterationCount < 2) {
                System.out.println("Итерация: " + (iterationCount + 1));
            }

            if (iterationCount >= 2) {
                double q = Math.abs(currentX - prevX) / Math.abs(prevPrevX - prevX);
                System.out.println("Итерация: " + (iterationCount + 1) + ", q = " + q);
                if (q >= 1 || currentX == prevX) return prevX;
            }

            prevPrevX = prevX;
            if (currentY * leftValue < 0.0) {
                prevX = b;
                b = currentX;
                rightValue = currentY;
            } else {
                prevX = a;
                a = currentX;
                leftValue = currentY;
            }
            iterationCount++;
        }

        return currentX;
    }

    static double newtonMethod(double initialGuess, int maxIterations) {
        System.out.println("Метод Ньютона:");

        double functionValue, derivativeValue, delta;
        double prevX = 0, prevPrevX = 0;
        int iterationCount = 0;
        double x = initialGuess;

        while (iterationCount < maxIterations) {
            functionValue = function(x);

            if (functionValue == 0.0) {
                return x;
            }

            if (iterationCount < 2) {
                double q = Math.abs(x - prevX) / Math.abs(prevPrevX - prevX);
                System.out.println("Итерация: " + (iterationCount + 1));
            }

            if (iterationCount >= 2) {
                double q = Math.abs(x - prevX) / Math.abs(prevPrevX - prevX);
                System.out.println("Итерация: " + (iterationCount + 1) + ", q = " + q);
                if (q >= 1 || x == prevX) return prevX;
            }

            derivativeValue = derivativeFunction(x);
            if (derivativeValue == 0.0) return Double.NaN;

            delta = functionValue / derivativeValue;
            prevPrevX = prevX;
            prevX = x;
            x -= delta;

            iterationCount++;
        }

        return x;
    }

    static double simpleIterationsMethod(double initialX, int maxIterations) {
        System.out.println("Метод простых итераций:");

        double upperBound = -18.238091;
        double lowerBound = -30.134383;
        double x0 = initialX;
        double x1 = phi(x0, lowerBound, upperBound);
        double x2 = phi(x1, lowerBound, upperBound);
        int iterationCount = 0;

        while (iterationCount < maxIterations) {
            if (iterationCount < 2) {
                double q = Math.abs(x2 - x1) / Math.abs(x0 - x1);
                System.out.println("Итерация: " + (iterationCount + 1));
            }

            if (iterationCount >= 2) {
                double q = Math.abs(x2 - x1) / Math.abs(x0 - x1);
                System.out.println("Итерация: " + (iterationCount + 1) + ", q = " + q);
                if (q >= 1 || x2 == x1) return x1;
            }

            x0 = x1;
            x1 = x2;
            x2 = phi(x1, lowerBound, upperBound);

            iterationCount++;
        }

        return x2;
    }

    static void launchAlgsGarwick(){
        double left = -5.65, right = -5.55;
        int n = 100;
        double x0 = -5.64;
        double result = newtonMethod(x0, n);
        System.out.println("Root x : " + result);

        System.out.println("################");

        result = simpleIterationsMethod(x0, n);
        System.out.println("Root x : " + result);

        System.out.println("################");

        result = chordMethod(left, right, n);
        System.out.println("x = " + result);
    }

    public static void main(String[] args) {
        launchAlgsGarwick();
    }
}