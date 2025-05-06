package com.github.chelovekkrokant;

import java.util.ArrayList;
import java.util.List;

public class Bisection {
    static int iterationCounter = 0;

    static double function(double x){
        return (Math.log(x*x+4) / (x+6) - 9);
    }

    static double derivativeFunction(double x){
        return ((2 * x) / ((x + 6) * (x * x + 4))) - ((Math.log(x * x + 4)) / ((x + 6) * (x + 6)));
    }

    static int calculateTheoreticalNumberOfIterations(double a, double b, double eps){
        double absDiff = Math.abs(b - a); // Вычисляем |b - a|
        int n = 0;
        // Находим минимальное n, при котором |b - a| / (2^(n+1)) < eps
        while (absDiff / Math.pow(2, n + 1) >= eps) {
            n++;
        }
        return n;
    }

    static double findRootUsingBisection(double leftBorder, double rightBorder, double eps){
        double middleValue = (leftBorder + rightBorder) / 2;
        double differenceBetweenFunctionValues = 0;
        double previousFunctionValue = 0;
        double currentFunctionValue = 0;
        System.out.printf("Запуск алгоритма на промежутке [%.5f;%.5f] при eps = %.5f\n", leftBorder, rightBorder, eps);
        System.out.printf("Значение функции на левом конце промежутка локализации : %.5f\n", function(leftBorder));
        System.out.printf("Значение функции на правом конце промежутка локализации : %.5f\n", function(rightBorder));
        System.out.println("Ожидаемое количество итераций = " + calculateTheoreticalNumberOfIterations(leftBorder, rightBorder, eps));

        while(true) {
            iterationCounter += 1;
            if(Math.abs(rightBorder - leftBorder) >= eps){
                middleValue = (leftBorder + rightBorder) / 2;

                if (function(leftBorder) * function(middleValue) < 0){
                    rightBorder = middleValue;
                } else {
                    if (function(middleValue) * function(rightBorder) < 0) {
                        leftBorder = middleValue;
                    }
                }


                previousFunctionValue = currentFunctionValue;
                currentFunctionValue = function(middleValue);
                differenceBetweenFunctionValues = Math.abs(currentFunctionValue - previousFunctionValue);
                System.out.printf("Итерация %d; Границы : (%f, %f); Разница значений функции на границах : %f\n", iterationCounter, leftBorder, rightBorder, differenceBetweenFunctionValues);
            } else {
                break;
            }
        }

        System.out.printf("Найденный корень = %.5f\n", middleValue);
        System.out.printf("Значение производной в найденном корне = %.5f\n", derivativeFunction(middleValue));
        System.out.printf("При eps = %.5f: V(deriv) = %.5f, V(delta) = %.5f\n", eps, 1 / Math.abs(derivativeFunction(middleValue)), eps / differenceBetweenFunctionValues);
        iterationCounter = 0;
        return middleValue;
    }

    static double findGarwickInterval(double leftBorder, double rightBorder, double eps){
        double middleValue = (leftBorder + rightBorder) / 2;
        iterationCounter = 0;

        List<Double> midValues = new ArrayList<>();
        double garWickInterval = 1;
        while(true) {
            middleValue = (leftBorder + rightBorder) / 2;
            midValues.add(middleValue);
            iterationCounter += 1;

            if (function(leftBorder) * function(middleValue) < 0){
                rightBorder = middleValue;
            } else {
                if (function(middleValue) * function(rightBorder) < 0) {
                    leftBorder = middleValue;
                }
            }
            System.out.println("Итерация " + iterationCounter + ". Интервал между последними значениями x = " + garWickInterval);
            if (midValues.size() >= 3){
                double prevValue = midValues.getLast();
                double prevPrevValue = midValues.get(midValues.size() - 2);
                double prevPrevPrevValue = midValues.get(midValues.size() - 3);
                garWickInterval = Math.abs(prevValue - prevPrevValue);
                if (Math.abs(prevValue - prevPrevValue) / Math.abs(prevPrevValue - prevPrevPrevValue) >= 1) {
                    break;
                }
            }
        }
        System.out.println("Найден интервал разболтки " + garWickInterval + " после " + iterationCounter +" итераций");
        iterationCounter = 0;
        return middleValue;
    }



}
