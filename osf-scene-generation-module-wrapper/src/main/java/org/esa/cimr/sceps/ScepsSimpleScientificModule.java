package org.esa.cimr.sceps;

/**
 * Class containing just a main method for first simple tests.
 *
 */
public class ScepsSimpleScientificModule {

    /**
     * The main method.
     *
     * @param args - Program arguments
     */
    public static void main(String[] args){
        if (args.length == 0) {
            System.out.println("No arguments given");
        } else {
            for (String arg : args) {
                System.out.println("arg = " + arg);
            }
        }
    }
}
