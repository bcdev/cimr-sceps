package org.esa.cimr.sceps;

/**
 * Class which could serve as Java 'Scientific Module' for simple tests.
 *
 */
public class ScepsSimpleJavaScientificModule {

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
