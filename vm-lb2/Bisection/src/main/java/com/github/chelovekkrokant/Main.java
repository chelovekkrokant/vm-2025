package com.github.chelovekkrokant;

import static com.github.chelovekkrokant.Bisection.findGarwickInterval;

public class Main{
    public static void main(String[] args){

//        Bisection.findRootUsingBisection(-5.65d, -5.55d, 0.00001);
//        Bisection.findRootUsingBisection(-5.99d, 0.d, 0.0001);
//        Bisection.findRootUsingBisection(-5.99d, 0.d, 0.001);
//        Bisection.findRootUsingBisection(-5.99d, 0.d, 0.01);
//        Bisection.findRootUsingBisection(-5.99d, 0.d, 0.1);

        Bisection.findGarwickInterval(-5.65d, 5.55d, 0.0001);
    }

}
