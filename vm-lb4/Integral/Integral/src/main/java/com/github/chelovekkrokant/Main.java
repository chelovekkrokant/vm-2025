package com.github.chelovekkrokant;

import java.util.Locale;

public class Main {
    static class Function {
        public double calculate(double x) {
            return Math.sqrt(1 - 0.49 * Math.sin(x) * Math.sin(x));
        }

        public double calculateFirstDerivative(double x) {
            double sinX = Math.sin(x);
            double cosX = Math.cos(x);
            double denominator = Math.sqrt(1 - 0.49 * sinX * sinX);
            return (-0.49 * sinX * cosX) / denominator;
        }

        public double calculateSecondDerivative(double x) {
            double sinX = Math.sin(x);
            double cosX = Math.cos(x);
            double term = 1 - 0.49 * sinX * sinX;
            double sqrtTerm = Math.sqrt(term);

            double numerator = -0.49 * (cosX * cosX - sinX * sinX) * sqrtTerm;
            numerator += (0.49 * 0.49 * sinX * sinX * cosX * cosX) / sqrtTerm;

            return numerator / term;
        }

        public double calculateThirdDerivative(double x) {
            double sin = Math.sin(x);
            double cos = Math.cos(x);
            double sin2 = sin * sin;
            double term = 1 - 0.49 * sin2;

            double numerator = 0.49 * sin * (1.47 - 4.41 * sin2 + 2.94 * Math.pow(sin, 4)
                    - 0.7203 * sin2 * (1 - sin2));
            double denominator = Math.pow(term, 2.5);

            return numerator / denominator;
        }

        public double calculateFourthDerivative(double x) {
            double sin = Math.sin(x);
            double cos = Math.cos(x);
            double sin2 = sin * sin;
            double term = 1 - 0.49 * sin2;

            double numerator = 0.49 * (2.205 - 13.23 * sin2 + 19.845 * Math.pow(sin, 4)
                    - 9.9225 * Math.pow(sin, 6) + 0.7203 * sin2 * (1 - sin2)
                    * (7 - 6 * sin2));
            double denominator = Math.pow(term, 3.5);

            return numerator / denominator;
        }
    }

    // Метод средних прямоугольников, основанный на суммировании площадей прямоугольников
    // С высотой значений функции в средних точках интервалов
    static class RECT {
        public double calculate(double x, int n, double h) {
            Function f = new Function();
            double integralSum = 0;
            for (int i = 0; i < n; ++i) {
                integralSum += f.calculate(x + i * h + h / 2);
            }
            return integralSum * h;
        }
    }
    // Метод трапеций, осованный на суммировании площадей трапеций
    // С перепадом высот, соответствующим значениям в крайних точкках
    static class TRAP {
        public double calculate(double x, int n, double h) {
            Function f = new Function();
            double integralSum = (f.calculate(x) + f.calculate(x + n * h)) / 2;
            for (int i = 1; i < n; ++i) {
                integralSum += f.calculate(x + i * h);
            }
            return integralSum * h;
        }
    }
    // Метод Симпсона, основанный на суммировании площадей параболических трапеций
    // Со "сглаженным" ребром, связывающим значения в крайних точках
    static class SIMPS {
        public double calculate(double x, int n, double h) {
            Function f = new Function();
            // Вычисление четного количества интервалов
            int m = n / 2;
            double integralSum = f.calculate(x) + f.calculate(x + h * n);
            // Умножение четных узлов на 2
            for (int i = 1; i < m; ++i) {
                integralSum += 2 * f.calculate(x + 2 * i * h);
            }
            // Умножение нечетных узлов на 4
            for (int i = 0; i < m; ++i) {
                integralSum += 4 * f.calculate(x + (2 * i + 1) * h);
            }
            return integralSum * h / 3; // Формула Симпсона
        }
    }

    public static void main(String[] args) {
        // Настройка локали для корректного вывода чисел с точкой
        Locale.setDefault(Locale.US);

        // Инициализация переменных
        int n = 0;
        // Текущее количество интервалов
        int next_n = 1;
        // Следующее количество интервалов
        int iter = 1;
        // Счетчик итераций
        double a = 0;
        // Нижний предел интегрирования
        double b = 1;
        // Верхний предел интегрирования
        double eps = 0.001;
        // Требуемая точность
        double h;
        // Шаг интегрирования

        // Переменные для хранения результатов и погрешностей
        double INT_VAL_PR = 0, INT_VAL_TR = 0, INT_VAL_SI = 0;
        double RES_INT_VAL_PR = 0, RES_INT_VAL_TR = 0, RES_INT_VAL_SI = 0;
        double H_VAL_PR = 0, H_VALUE_TR = 0, H_VAL_SI = 0;

        // Флаги для контроля точности методов
        boolean flagP = true, flagTR = true, flagS = true;

        // Переменные для хранения оптимальных параметров
        int N_VAL_PR = 0, N_VALUE_TR = 0, N_VAL_SI = 0;

        // Создание экземпляров классов методов интегрирования
        RECT IntP = new RECT();
        TRAP IntTR = new TRAP();
        SIMPS IntS = new SIMPS();
        Function ff = new Function();
        // Основной цикл
        while (next_n != n) {
            n = next_n;
            h = (b - a) / n; // Вычисление текущего шага

            // Вывод информации о текущей итерации
            System.out.printf("Итерация %d:\tКоличество шагов n = %d; Длина шага h = %.3f\n", iter, n, h);

            // Блок обработки метода прямоугольников
            if (flagP) {
                H_VAL_PR = h;    // Сохраняем текущий шаг
                N_VAL_PR = n;    // Сохраняем текущее n
                // Уточнение результата по Рунге
                INT_VAL_PR = IntP.calculate(a, 2 * n, h / 2) +
                        (IntP.calculate(a, 2 * n, h / 2) - IntP.calculate(a, n, h)) / 3;
                next_n = 2 * n; // Удваиваем количество интервалов
                // Вычисление текущей погрешности
                System.out.printf("Значение интеграла по методу прямоугольников  = %.6f ± %.6f\n",
                        INT_VAL_PR, Math.abs(IntP.calculate(a, n, h) - IntP.calculate(a, 2 * n, h / 2)) / 3);
            } else {
                // Вывод финального результата для метода
                RES_INT_VAL_PR = Math.abs(IntP.calculate(a, N_VAL_PR, H_VAL_PR) - IntP.calculate(a, 2 * N_VAL_PR, H_VAL_PR / 2)) / 3;
                System.out.printf("Значение интеграла по методу прямоугольников  = %.6f ± %.6f При n = %d, h = %.3f\n",
                        INT_VAL_PR, RES_INT_VAL_PR, N_VAL_PR, H_VAL_PR);
            }
            // Проверка достижения требуемой точности
            if (Math.abs(IntP.calculate(a, n, h) - IntP.calculate(a, 2 * n, h / 2)) / 3  < eps)
                flagP = false;

            // Блок для метода трапеций
            if (flagTR) {
                H_VALUE_TR = h;
                N_VALUE_TR = n;
                INT_VAL_TR = IntTR.calculate(a, 2 * n, h / 2) +
                        (IntTR.calculate(a, 2 * n, h / 2) - IntTR.calculate(a, n, h)) / 3;
                next_n = 2 * n;
                System.out.printf("Значение интеграла по методу трапеций  = %.6f ± %.6f\n",
                        INT_VAL_TR, Math.abs(IntTR.calculate(a, n, h) - IntTR.calculate(a, 2 * n, h / 2)) / 3);
            } else {
                RES_INT_VAL_TR = Math.abs(IntTR.calculate(a, N_VALUE_TR, H_VALUE_TR) - IntTR.calculate(a, 2 * N_VALUE_TR, H_VALUE_TR / 2)) / 3;
                System.out.printf("Значение интеграла по методу трапеций  = %.6f ± %.6f При n = %d, h = %.3f\n",
                        INT_VAL_TR, RES_INT_VAL_TR, N_VALUE_TR, H_VALUE_TR);
            }
            // /(2^k -1 ) рассматриваем порядок точности = 2
            if (Math.abs(IntTR.calculate(a, n, h) - IntTR.calculate(a, 2 * n, h / 2)) / 3 < eps)
                flagTR = false;

            // Блок для метода Симпсона
            if (flagS) {
                H_VAL_SI = h;
                N_VAL_SI = n;
                INT_VAL_SI = IntS.calculate(a, 2 * n, h / 2) +
                        (IntS.calculate(a, 2 * n, h / 2) - IntS.calculate(a, n, h)) / 15;
                next_n = 2 * n;
                System.out.printf("Значение интеграла по методу Симпсона  = %.6f ± %.6f\n",
                        INT_VAL_SI, Math.abs(IntS.calculate(a, n, h) - IntS.calculate(a, 2 * n, h / 2)) / 15);
            } else {
                // Вывод финального результата для метода
                RES_INT_VAL_SI = Math.abs(IntS.calculate(a, N_VAL_SI, H_VAL_SI) - IntS.calculate(a, 2 * N_VAL_SI, H_VAL_SI / 2)) / 15;
                System.out.printf("Значение интеграла по методу Симпсона  = %.6f ± %.6f При n = %d, h = %.3f\n",
                        INT_VAL_SI, RES_INT_VAL_SI, N_VAL_SI, H_VAL_SI);
            }
            // порядок К для метода симпсона = 4
            if (Math.abs(IntS.calculate(a, n, h) - IntS.calculate(a, 2 * n, h / 2)) / 15 < eps)
                flagS = false;

            System.out.println();
            iter++;
        }
    }
}