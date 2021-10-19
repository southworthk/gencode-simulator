/*
 * Copyright 2021 Kerry Mark Southworth
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated
 * documentation files (the "Software"), to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software,
 * and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all copies or substantial portions
 * of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED
 * TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
 * CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

package com.sqnsolutions.simulator;

import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.Random;
import java.util.Scanner;
import java.util.concurrent.ThreadLocalRandom;

/**
 * This application simulates genetic code evolution by introducing random changes to a given sequence,
 * translating it into its amino acid equivalent, and comparing the mutated source to a known target.
 * As the simulation concludes, it reports the number of successful trials for a given set of parameters.
 */
public class Simulator {

    private static final Random rand = new Random();
    private static final DateFormat dateFormat = new SimpleDateFormat("MM/dd/yyyy HH:mm:ss");
    private static Long totalTrials = 100000000L; // the total number of trials to complete, 100 million
   //private static Long totalTrials = 10000000L;
    private static int abandonmentThreshold = 100; // the point at which to start a new trial if there is no match
    private static int trialStatusThreshold = 100000000; // the display frequency, default to one hundred million
    //private static int trialStatusThreshold = 1000000;

    private static final String VERSION = "1.12";
    private static ArrayList<Integer> aaDiff;

    /**
     * The entry point of the Simulator application.
     *
     * @param args the input arguments, which include such things as the desired sequence length and the
     *             percentage difference between source and target
     */
    public static void main( String[] args )
    {
        aaDiff = new ArrayList<>();
        Scanner scanner = new Scanner(System.in);

        System.out.println("Version " + VERSION);
        System.out.println("Default settings:");
        System.out.println("Total trials: " + totalTrials);
        System.out.println("Abandonment threshold: " + abandonmentThreshold);
        System.out.println("Trial status threshold: " + trialStatusThreshold);

        System.out.println("How many nucleotides should be in the simulated string (must be a multiple of 3 "+
                "greater than 50)?");
        int nucleotideLen = scanner.nextInt();
        if ((nucleotideLen < 51) || ((nucleotideLen % 3) > 0)) {
            System.out.println("Number of simulated nucleotide string must be a multiple of 3 and greater than 50");
            System.exit(1);
        }
        System.out.println("By what percentage should the target and source differ (must be whole number greater "+
                "than 0 and less than 100)?");
        int pctDiff = scanner.nextInt();
        if ((pctDiff > 99) || (pctDiff < 1)) {
            System.out.println("Percent difference must be 1 or greater and ");
            System.exit(1);
        }

        System.out.println("Do you want to change the existing defaults for total trials, abandonment "+
                        "threshold and so on (y/n)?");
        String changeDefaults = scanner.next().toUpperCase();

        if(changeDefaults.substring(0,1).equals("Y")) {
            System.out.println("How many hundreds of millions of trials should be carried out for the specified strings? "+
                    "Default is 1 billion.");
            totalTrials = scanner.nextInt() * 100000000L; // multiply by a hundred million

            System.out.println("After how many attempts should a new trial start if there is not match? "+
                    "Default is 100.");
            abandonmentThreshold = scanner.nextInt();

            System.out.println("After how many millions of trials should the status be reported? "+
                            "Default is 100,000,000.");
            trialStatusThreshold = scanner.nextInt() * 1000000; // multiply by a million
        }
        // create a nucleotide sequence that once translated, will serve as the target
        String targetSeq = getRandomNSequence(nucleotideLen).toString();
        System.out.println( "\nTarget nucleotide sequence:          "+targetSeq );
        String targetSeqTranslated = translateNSeq( new StringBuilder(targetSeq));

        // create a copy of the target and make random changes, this copy will be where mutations occur
        String copySeq = createDifferentSequence(new StringBuilder(targetSeq), targetSeqTranslated, pctDiff).toString();
        System.out.println( "Source nucleotide sequence:          "+getMarkedDifferences(targetSeq,
                copySeq,false));

        System.out.println("\nNucleotide sequence length: "+nucleotideLen);
        System.out.println( "Specified difference: "+pctDiff+"%");
        System.out.println( "Actual similarity of nucleotide sequences: "+calculateNucleotideSimilarity(targetSeq,
                copySeq)+"%");

        String copySeqTranslated = translateNSeq( new StringBuilder(copySeq));

        System.out.println( "\nTarget amino acid sequence:          "+targetSeqTranslated);
        System.out.println( "Source amino acid sequence:          "+getMarkedDifferences(targetSeqTranslated,
                copySeqTranslated,true));

        System.out.println("\nAmino acid sequence length: "+(nucleotideLen/3));
        System.out.println( "Similarity of amino acid sequences: "+calculateAASimilarity(targetSeqTranslated,
                copySeqTranslated)+"%");
        System.out.println( "Position of amino acid difference(s): "+getAADifferenceStr());
        simulate(totalTrials, copySeq, targetSeqTranslated);

    }

    /**
     * Method creates a string that is appended to "Position of amino acid difference(s): " in the preliminary
     * display.
     * @return a string with the position of each amino acid difference
     */
    private static String getAADifferenceStr(){
        StringBuilder sb = new StringBuilder();
        int aaSize = aaDiff.size();
        int index = 1;

        for(Integer pos : aaDiff){
            if (index == aaSize) {
                sb.append((pos + 1) + " ");
            } else {
                sb.append((pos + 1) + ", ");
            }
            index++;
        }

        return sb.toString();
    }

    /**
     * Method compares two strings to determine exactly where they differ and updates the ArrayList aaDiff accordingly.
     * @param one a sequence of genetic code
     * @param two another sequence of genetic code
     * @param aaString boolean indicating whether the strings are comprised of amino acid residues (true) or
     *                 nucleotides (false)
     * @return a string with known difference(s) as the source string
     */
    private static String getMarkedDifferences(String one, String two, boolean aaString){
        int oneLen = one.length();
        int twoLen = two.length();

        StringBuilder stringBuilder = new StringBuilder();

        if (oneLen == twoLen){
            char[] oneArr = one.toCharArray();
            char[] twoArr = two.toCharArray();
            for(int i=0; i < oneLen; i++){
                if(oneArr[i] == twoArr[i]){
                    stringBuilder.append(twoArr[i]);
                }else{
                    if (aaString) aaDiff.add(i);
                    stringBuilder.append(twoArr[i]);
                }
            }
        }

        return stringBuilder.toString();
    }

    /**
     * Method for carrying out the mutation of the source string and its comparison with a target for a set
     * number of trials.
     * @param trials the number of trials the simulation with run
     * @param copySeq a copy of the nucleotide source string that will undergo random changes
     * @param targetSeqTranslated the translated target string that will be used as a basis of comparison
     */
    private static void simulate(long trials, String copySeq, String targetSeqTranslated){
        int matches = 0;
        long startTime = System.currentTimeMillis();
        int copySeqLen = copySeq.length();
        StringBuilder mutantCopy = new StringBuilder(copySeq);
        double totalAvgSimilarity = 0.0;
        double maxSimilarity = 0.0;
        long reportCounter = 0;
        long totalTests = 0;

        for (long i = 0; i < trials; i++)
        {
            long generation = 1;
            while (true)
            {
                mutate(mutantCopy, copySeqLen);
                String mutantCopyTranslated = translateNSeq(mutantCopy);

                if(mutantCopyTranslated.equals(targetSeqTranslated))
                {
                    reportCounter++;
                    matches++;
                    break;
                }
                else if (generation > abandonmentThreshold) // abandon trial if it hasn't found a match by a certain point
                {
                    /*
                     the user can look at the reported maxSimilarity to decide whether the abandonment
                      threshold needs to be adjusted
                    */
                    double similarity = calculateAASimilarity(targetSeqTranslated, mutantCopyTranslated);
                    maxSimilarity = Math.max(similarity, maxSimilarity);
                    totalAvgSimilarity = (totalAvgSimilarity + similarity);
                    mutantCopy = new StringBuilder(copySeq);
                    generation = 0;
                    reportCounter++;
                    break;
                }
                else
                {
                    generation++;
                }
                totalTests++;
            }
            if(reportCounter > trialStatusThreshold) // print status update
            {
                printStatus(startTime, matches, totalTests, totalAvgSimilarity, maxSimilarity);
                reportCounter = 0;
            }

        }
        // print final status
        printStatus(startTime, matches, totalTests, totalAvgSimilarity, maxSimilarity);

    }

    /**
     * Method displays a status report to the console.
     * @param startTime when the simulation started
     * @param matches how many times the source string evolved to match the target string
     * @param totalTests the total number of tests run to this point
     * @param totalAvgSimilarity the average similarity between source and target at the point of abandonment
     * @param maxSimilarity the maximum or closest similarity between source and target at abandonment
     */
    private static void printStatus(long startTime, int matches, long totalTests, double totalAvgSimilarity, double maxSimilarity){
        long endTime = System.currentTimeMillis();
        float duration = (endTime - startTime) / 1000; // elapsed seconds

        int trials = Math.round(totalTests/ abandonmentThreshold) + matches;
        System.out.println("*********************************************************************************");

        Date date = new Date();
        System.out.println("Completed: "+dateFormat.format(date));
        System.out.println("Trials: "+trials+", total matches: "+matches+", elapsed time: "+duration+" seconds");
        System.out.printf("Similarity of AA seq at abandonment, avg: %.2f%%, max: %.2f%% %n",((totalAvgSimilarity/(totalTests))*100),maxSimilarity);

    }

    /**
     * Method creates a string that differs from the target by a specified percentage.
     * @param target the target nucleotide sequence
     * @param targetAA the target amino acid sequence
     * @param pctDiff the desired percentage which may differ from the actual percentage because of the need
     *                for whole number differences
     * @return a source string that differs from the target string by a whole number
     */
    private static StringBuilder createDifferentSequence(StringBuilder target, String targetAA, int pctDiff)
    {
        StringBuilder sequence;

        while(true)
        {
            sequence = target;
            int seqLen = sequence.length();
            if(seqLen < 51) break;

            int numberOfChanges = Math.round((pctDiff / (float) 100) * seqLen);

            for (int i = 0; i < numberOfChanges; i++)
            {
                sequence = mutate(sequence, seqLen);
            }

            String sequenceAA = translateNSeq(sequence);

            // Exit loop if the target and the copy don't have the same translation (AA sequence)
            if(!targetAA.equals(sequenceAA))
            {
                break;
            }
        }

        return sequence;
    }

    /**
     * Method calculates the actual percentage similarity between nucleotide strings.
     * @param target the target sequence
     * @param copy the source sequence
     * @return a number representing the calculated percentage similarity between the two strings
     */
    private static double calculateNucleotideSimilarity(String target, String copy){
        double pctSimilarity;

        int len = target.length();

        char[] targetArr = target.toCharArray();
        char[] copyArr = copy.toCharArray();

        int matchCount = 0;

        for (int i = 0; i < len; i++ ){
            if( targetArr[i] == copyArr[i]){
                matchCount++;
            }
        }

        pctSimilarity = (matchCount / (float) len) * 100;

        return round(pctSimilarity, 3);
    }

    /**
     * Method calculates the actual percentage similarity between amino acid strings.
     * @param target the target sequence
     * @param copy the source sequence
     * @return a number representing the calculated percentage similarity between the two strings
     */
    private static double calculateAASimilarity(String target, String copy){
        double pctSimilarity;

        int len = target.length();

        char[] targetArr = target.toCharArray();
        char[] copyArr = copy.toCharArray();

        int matchCount = 0;

        for (int i = 0; i < len; i++ ){
            if( targetArr[i] == copyArr[i]){
                matchCount++;
            }
        }

        pctSimilarity = (matchCount / (float) len) * 100;

        return round(pctSimilarity, 3);
    }

    /**
     * Method rounds doubles to specified precision.
     *
     * @param value     the value
     * @param precision the precision
     * @return the rounded double
     */
    public static double round (double value, int precision) {
        int scale = (int) Math.pow(10, precision);
        return (double) Math.round(value * scale) / scale;
    }

    /**
     * Method introduces a random mutation into a string of a given length.
     *
     * @param original the string to be mutated
     * @param seqLen the length of the string
     * @return
     */
    private static StringBuilder mutate(StringBuilder original, int seqLen)
    {
        StringBuilder mutant = original;

        // get random position
        int randPos = rand.nextInt(seqLen);

        // get string at that position
        String currentValue = String.valueOf(mutant.charAt(randPos));

        // get random nucleotide, ensure that it is not the same as the current value
        String randomNucleotide;
        do {
            randomNucleotide = getRandomNucleotide();
        } while (currentValue.equals(randomNucleotide));

        // make substitution
        mutant.setCharAt(randPos, randomNucleotide.charAt(0));

        return mutant;

    }

    /**
     * Method returns a random nucleotide.
     *
     * @return the random value
     */
    private static String getRandomNucleotide(){

        String randomNucleotide = null;

        int randomNum = ThreadLocalRandom.current().nextInt(1, 4 + 1);

        switch(randomNum)
        {
            case 1:
                randomNucleotide = "a"; // Adenine
                break;
            case 2:
                randomNucleotide = "c"; // Cytosine
                break;
            case 3:
                randomNucleotide = "g"; // Guanine
                break;
            case 4:
                randomNucleotide = "u"; // Uracil
                break;
        }

        return randomNucleotide;
    }

    /**
     * Method returns a StringBuilder object that corresponds to the translation
     * of the source argument.
     *
     * @param source a string that corresponds to a ribonucleic acid (RNA) sequence
     * @return       a string that corresponds to the translation of the source
     *               to a sequence of amino acid residues
     */
    private static String translateNSeq(StringBuilder source)
    {

        // instantiate new StringBuilder object to hold the translated string
        StringBuilder translation = new StringBuilder();

        // calculate the length of the source sequence
        int len = source.length();

        // translate the source sequence three characters at a time.
        for( int i = 0; i < len; i+=3 )
        {
            // get the next three nucleotide symbols in the sequence
            String codon = source.substring(i,i+3);

            // if the string is exactly three symbols in length
            // translate it and add to the translation StringBuilder object
            if (codon.length() == 3)
            {
                translation.append(translateNtoA(codon));
            }
        }

        return translation.toString();
    }

    /**
     * Method translates a nucleotide sequence into its amino acid equivalent.
     *
     * @param source nucleotide sequence
     * @return equivalent amino acid sequence
     */
    private static String translateNtoA(String source)
    {
        String translation = null;

        switch(source.toUpperCase())
        {
            case "UUU":
            case "UUC":
                translation = "F"; // phenylalanine
                break;
            case "UUA":
            case "UUG":
            case "CUU":
            case "CUC":
            case "CUA":
            case "CUG":
                translation = "L"; // leucine
                break;
            case "AUU":
            case "AUC":
            case "AUA":
                translation = "I"; // isoleucine
                break;
            case "AUG":
                translation = "M"; // methionine
                break;
            case "GUU":
            case "GUC":
            case "GUA":
            case "GUG":
                translation = "V"; // valine
                break;
            case "UCU":
            case "UCC":
            case "UCA":
            case "UCG":
            case "AGU":
            case "AGC":
                translation = "S"; // serine
                break;
            case "CCU":
            case "CCC":
            case "CCA":
            case "CCG":
                translation = "P"; // proline
                break;
            case "ACU":
            case "ACC":
            case "ACA":
            case "ACG":
                translation = "T"; // threonine
                break;
            case "GCU":
            case "GCC":
            case "GCA":
            case "GCG":
                translation = "A"; // alanine
                break;
            case "UAU":
            case "UAC":
                translation = "Y"; // tyrosine
                break;
            case "UAA":
            case "UAG":
            case "UGA":
                translation = "*"; // stop signal, no amino acid equivalent
                break;
            case "CAU":
            case "CAC":
                translation = "H"; // histidine
                break;
            case "CAA":
            case "CAG":
                translation = "Q"; // glutamine
                break;
            case "AAU":
            case "AAC":
                translation = "N"; // asparagine
                break;
            case "AAA":
            case "AAG":
                translation = "K"; // lysine
                break;
            case "GAU":
            case "GAC":
                translation = "D"; // aspartic acid
                break;
            case "GAA":
            case "GAG":
                translation = "E"; // glutamic acid
                break;
            case "UGU":
            case "UGC":
                translation = "C"; // cysteine
                break;
            case "UGG":
                translation = "W"; // tryptophan
                break;
            case "CGU":
            case "CGC":
            case "CGA":
            case "CGG":
            case "AGA":
            case "AGG":
                translation = "R"; // arginine
                break;
            case "GGU":
            case "GGC":
            case "GGA":
            case "GGG":
                translation = "G"; // glycine
                break;
        }

        return translation;
    }

    /**
     * Method generates a random nucleotide sequence of a given length.
     *
     * @param length of the desired sequence
     * @return a string of nucleotides for a specified length
     */
    private static StringBuilder getRandomNSequence(int length)
    {
        StringBuilder nSeq = new StringBuilder();

        for( int i = 0; i < length; i++)
        {
            nSeq.append(getRandomNucleotide());
        }

        return nSeq;
    }

}
