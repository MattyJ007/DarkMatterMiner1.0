package sample;
import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Random;

import static jdk.nashorn.internal.objects.NativeString.indexOf;
class AnalyseSequence {
    static void analyseSequence(Sequence newSeq, boolean secureRandom){
        System.out.println(newSeq.getName());
        setCountForObservedORFs();
        setTotalExpectedFrequencies();
        newSeq.setGcContent(gCcount((newSeq.getRawSeq())));
        generatePermutations(newSeq.getRawSeq(), secureRandom);

        chiSquare(allPermutationFreqArrays, newSeq);

    }

    private static Integer[] totalExpectedFrequenciesTri = new Integer[64];
    private static Integer[] totalExpectedFrequenciesDi = new Integer[16];
    private static ArrayList<ArrayList<Integer>> allPermutationFreqArrays = new ArrayList<>();
    private static ArrayList<Integer> trinucFreqFrame1Temp = new ArrayList<>();
    private static ArrayList<Integer> trinucFreqFrame2Temp = new ArrayList<>();
    private static ArrayList<Integer> trinucFreqFrame3Temp = new ArrayList<>();
    private static ArrayList<Integer> dinucFreqFrame1Temp = new ArrayList<>();
    private static ArrayList<Integer> dinucFreqFrame2Temp = new ArrayList<>();
    private static ArrayList<Integer> longestORFs = new ArrayList<>();
    private static ArrayList<Integer> observedLongestORFs = new ArrayList<>();
    private static int countForObservedORFs;
    private static void setCountForObservedORFs(){
        countForObservedORFs = 0;
    }
    private static void setTotalExpectedFrequencies(){
        for (int j = 0; j<64;j++){
            totalExpectedFrequenciesTri[j] = 0;
        }
        for (int j = 0; j<16;j++){
            totalExpectedFrequenciesDi[j] = 0;
        }
    }
    private static void clearFrameArrays(){
        trinucFreqFrame1Temp.clear();
        trinucFreqFrame2Temp.clear();
        trinucFreqFrame3Temp.clear();
        dinucFreqFrame1Temp.clear();
        dinucFreqFrame2Temp.clear();
    }

    private static double gCcount(String seq){
        //** gets GC content of each sequence
        int countGC = countMatches(seq, "G") + countMatches(seq, "C");
        return (double)countGC/seq.length();
    }

    private static int countMatches(String str, String sub) {
        if (str.isEmpty() || sub.isEmpty()) {return -1;}
        int count = 0;
        int idx = 0;
        while ((idx = indexOf(str, sub, idx)) != -1) {
            count++;
            idx += sub.length();
        }
        return count;
    }

    private static void generatePermutations(String seqLine, boolean secureRandom){
        int permutations = DarkMatterMinerUI.getPermutations();
        SecureRandom sRand;
        Random nRand;
//        OpenReadingFrameDistances.stopper(seqLine,0);
        int randomNum;
        int len;
        String temp;
        String randSeq;
        String[] seq;
        getFrequency(seqLine);
        if (secureRandom) {
            sRand = new SecureRandom();
            for (int j = 0; j< permutations; j++){
                len = seqLine.length();
                seq = seqLine.split("");
                for (int i=0; i<len; i++){
                    randomNum = sRand.nextInt(len);
                    temp = seq[i];
                    seq[i] = seq[randomNum];
                    seq[randomNum] = temp;
                }
                randSeq = String.join("", (CharSequence[]) seq);
                getFrequency(randSeq);
//                OpenReadingFrameDistances.stopper(randSeq, 1);
            }
        }
        else {
            nRand = new Random();
            for (int j = 0; j< permutations; j++){
                len = seqLine.length();
                seq = seqLine.split("");
                for (int i=0; i<len; i++){
                    randomNum = nRand.nextInt(len);
                    temp = seq[i];
                    seq[i] = seq[randomNum];
                    seq[randomNum] = temp;
                }
                randSeq = String.join("", (CharSequence[]) seq);
                getFrequency(randSeq);
//                OpenReadingFrameDistances.stopper(randSeq, 1);
            }
        }
    }

    private static void getFrequency(String seq){
        clearFrameArrays();
        ArrayList<ArrayList<Integer>> orfindeces = new ArrayList<>();
        ArrayList<ArrayList<Integer>> unsortedMotifFrequencies = new ArrayList<>();
        String motifString = "AAA,AAT,AAG,AAC,ATA,ATT,ATG,ATC,AGA,AGT,AGG,AGC,ACA,ACT,ACG,ACC,TAA,TAT,TAG,TAC,TTA,TTT,TTG,TTC,TGA,TGT,TGG,TGC,TCA,TCT,TCG,TCC,GAA,GAT,GAG,GAC,GTA,GTT,GTG,GTC,GGA,GGT,GGG,GGC,GCA,GCT,GCG,GCC,CAA,CAT,CAG,CAC,CTA,CTT,CTG,CTC,CGA,CGT,CGG,CGC,CCA,CCT,CCG,CCC,AA,AT,AG,AC,TA,TT,TG,TC,GA,GT,GG,GC,CA,CT,CG,CC";
        Integer[] stopCodons = {16,18,20,24,28,52};
        String[] motifList = motifString.split(",");
        for (String motif: motifList){
            ArrayList<Integer> motifIndeces = new ArrayList<>();
            int index = seq.indexOf(motif);
            while(index >= 0) {
                motifIndeces.add(index);
                index = seq.indexOf(motif, index+1);
            }
            unsortedMotifFrequencies.add(motifIndeces);
        }
        for (int d: stopCodons){
            orfindeces.add(unsortedMotifFrequencies.get(d));
        }
        getORFLoci(orfindeces, seq.length());
        sortFrequenciesIntoFrames(unsortedMotifFrequencies);
    }

    private static void sortFrequenciesIntoFrames(ArrayList<ArrayList<Integer>> unsortedMotifFrequencies){
        int frame;
        for (int j = 0; j<80; j++){
            int motifFrameT1 = 0;
            int motifFrameT2 = 0;
            int motifFrameT3 = 0;
            int motifFrameD1 = 0;
            int motifFrameD2 = 0;
            if(j<64){
                for(int motifIndex: unsortedMotifFrequencies.get(j)) {
                    frame = motifIndex % 3;
                    if (frame == 0) {
                        motifFrameT1++;
                    }
                    else if (frame == 1) {
                        motifFrameT2++;
                    }
                    else {
                        motifFrameT3++;
                    }
                }
                trinucFreqFrame1Temp.add(motifFrameT1);
                trinucFreqFrame2Temp.add(motifFrameT2);
                trinucFreqFrame3Temp.add(motifFrameT3);
            }
            else{
                for(int motifIndex: unsortedMotifFrequencies.get(j)) {
                    frame = motifIndex % 2;
                    if (frame == 0) {
                        motifFrameD1++;
                    }
                    else {
                        motifFrameD2++;
                    }
                }
                dinucFreqFrame1Temp.add(motifFrameD1);
                dinucFreqFrame2Temp.add(motifFrameD2);
            }
        }
        for(int i = 0; i<64; i++) {
            totalExpectedFrequenciesTri[i] += trinucFreqFrame1Temp.get(i) + trinucFreqFrame2Temp.get(i) + trinucFreqFrame3Temp.get(i);
        }
        for(int i = 0; i<16; i++) {
            totalExpectedFrequenciesDi[i] += dinucFreqFrame1Temp.get(i) + dinucFreqFrame2Temp.get(i);
        }
        ArrayList<Integer> trinucFreqFrame1 = new ArrayList<>(trinucFreqFrame1Temp);
        ArrayList<Integer> trinucFreqFrame2 = new ArrayList<>(trinucFreqFrame2Temp);
        ArrayList<Integer> trinucFreqFrame3 = new ArrayList<>(trinucFreqFrame3Temp);
        ArrayList<Integer> dinucFreqFrame1 = new ArrayList<>(dinucFreqFrame1Temp);
        ArrayList<Integer> dinucFreqFrame2 = new ArrayList<>(dinucFreqFrame2Temp);
        allPermutationFreqArrays.add(trinucFreqFrame1);
        allPermutationFreqArrays.add(trinucFreqFrame2);
        allPermutationFreqArrays.add(trinucFreqFrame3);
        allPermutationFreqArrays.add(dinucFreqFrame1);
        allPermutationFreqArrays.add(dinucFreqFrame2);
    }

    private static void chiSquare( ArrayList<ArrayList<Integer>> observedArrays, Sequence newSeq){
        double exp;
        int count = 0;
        ArrayList<Double> observedChiValuesT = new ArrayList<>();
        ArrayList<Double> expectedChiValuesT = new ArrayList<>();
        ArrayList<Double> observedChiValuesD = new ArrayList<>();
        ArrayList<Double> expectedChiValuesD = new ArrayList<>();
        for (ArrayList<Integer> observedArray :observedArrays) {
            double chiT = 0;
            if (observedArray.size() == 64){
                for(int n = 0; n<64;n++){
                    exp = (totalExpectedFrequenciesTri[n]-observedArray.get(n))/((double) (DarkMatterMinerUI.getPermutations()*3)+2);
//                    System.out.println(exp + " = " + totalExpectedFrequenciesTri[n]+" - "+ observedArray.get(n)+" / "+ DarkMatterMinerUI.getPermutations()+" * 3" + " + 2");
                    chiT += Math.pow(exp - observedArray.get(n),2)/exp;
//                    System.out.println(chi + " +=( " +exp + " - "+observedArray.get(n)+")**2/ "+exp );
                }
//                System.out.println(chi + "  ======");
                if (count <3){
                    observedChiValuesT.add(chiT);
                    expectedChiValuesT.add(chiT);
                }
                else{
                    expectedChiValuesT.add(chiT);
                }
            }
            double chiD = 0;
            if(observedArray.size() == 16){
                for(int n = 0; n<16;n++){
                    exp = (totalExpectedFrequenciesDi[n]-observedArray.get(n))/((double) (DarkMatterMinerUI.getPermutations()*2)+1);
                    chiD += Math.pow(exp - observedArray.get(n),2)/exp;
                }
                if (count>=3 && count<5){
                    observedChiValuesD.add(chiD);
                    expectedChiValuesD.add(chiD);
                }
                else{
                    expectedChiValuesD.add(chiD);
                }
            }
            count++;
        }
        Collections.sort(expectedChiValuesT);
        Collections.sort(expectedChiValuesD);
        getMotifPValues(observedChiValuesD, observedChiValuesT,expectedChiValuesD,expectedChiValuesT, newSeq);
    }

    private static void getMotifPValues(ArrayList<Double> observedChiValuesD, ArrayList<Double> observedChiValuesT, ArrayList<Double> expectedChiValuesD, ArrayList<Double> expectedChiValuesT, Sequence newSeq){
        newSeq.setDinuc1((1-expectedChiValuesD.indexOf(observedChiValuesD.get(0))/(double) expectedChiValuesD.size()));
        newSeq.setDinuc2((1-expectedChiValuesD.indexOf(observedChiValuesD.get(1))/(double) expectedChiValuesD.size()));
        if ((1-expectedChiValuesD.indexOf(observedChiValuesD.get(0))/(double) expectedChiValuesD.size())<(1-expectedChiValuesD.indexOf(observedChiValuesD.get(1))/(double) expectedChiValuesD.size())){
            newSeq.setDinucleotidePvalue((1-expectedChiValuesD.indexOf(observedChiValuesD.get(0))/(double) expectedChiValuesD.size()));
        }
        else {
            newSeq.setDinucleotidePvalue((1-expectedChiValuesD.indexOf(observedChiValuesD.get(1))/(double) expectedChiValuesD.size()));
        }
        newSeq.setTrinucelotidePValue((1-expectedChiValuesT.indexOf(observedChiValuesT.get(0))/(double) expectedChiValuesT.size()));
        newSeq.setTrinuc1((1-expectedChiValuesT.indexOf(observedChiValuesT.get(0))/(double) expectedChiValuesT.size()));

        if ((1-expectedChiValuesT.indexOf(observedChiValuesT.get(1))/(double) expectedChiValuesT.size()) < newSeq.getTrinucelotidePValue()){
            newSeq.setTrinucelotidePValue((1-expectedChiValuesT.indexOf(observedChiValuesT.get(1))/(double) expectedChiValuesT.size()));
        }
        newSeq.setTrinuc2((1-expectedChiValuesT.indexOf(observedChiValuesT.get(1))/(double) expectedChiValuesT.size()));
        if ((1-expectedChiValuesT.indexOf(observedChiValuesT.get(2))/(double) expectedChiValuesT.size()) < newSeq.getTrinucelotidePValue()){
            newSeq.setTrinucelotidePValue((1-expectedChiValuesT.indexOf(observedChiValuesT.get(2))/(double) expectedChiValuesT.size()));
        }
        newSeq.setTrinuc3((1-expectedChiValuesT.indexOf(observedChiValuesT.get(2))/(double) expectedChiValuesT.size()));
//        System.out.println(newSeq.getName());
//        System.out.println(newSeq.getDinucleotidePValue());
//        System.out.println(newSeq.getTrinucelotidePValue());
    }

    private static void getORFLoci(ArrayList<ArrayList<Integer>> orfindeces, int len){
        //** values in arrays indicate the starting loci of stop codons
        int numberOfCodons = len%3;
        ArrayList<Integer> stopCodonsFrame1F;
        ArrayList<Integer> stopCodonsFrame2F;
        ArrayList<Integer> stopCodonsFrame3F;
        ArrayList<Integer> stopCodonsFrame1R;
        ArrayList<Integer> stopCodonsFrame2R;
        ArrayList<Integer> stopCodonsFrame3R;
        if (numberOfCodons == 0) {
            stopCodonsFrame1F = new ArrayList<>(Arrays.asList(-3,len));
            stopCodonsFrame2F = new ArrayList<>(Arrays.asList(-2,len-2));
            stopCodonsFrame3F = new ArrayList<>(Arrays.asList(-1,len-1));
            stopCodonsFrame1R = new ArrayList<>(Arrays.asList(-3,len));
            stopCodonsFrame2R = new ArrayList<>(Arrays.asList(-2,len-2));
            stopCodonsFrame3R = new ArrayList<>(Arrays.asList(-1,len-1));
        }
        else if(numberOfCodons ==1){
            stopCodonsFrame1F = new ArrayList<>(Arrays.asList(-3, len-1));
            stopCodonsFrame2F = new ArrayList<>(Arrays.asList(-2, len));
            stopCodonsFrame3F = new ArrayList<>(Arrays.asList(-1, len-2));
            stopCodonsFrame1R = new ArrayList<>(Arrays.asList(-3, len-1));
            stopCodonsFrame2R = new ArrayList<>(Arrays.asList(-2, len));
            stopCodonsFrame3R = new ArrayList<>(Arrays.asList(-1, len-2));
        }
        else{
            stopCodonsFrame1F = new ArrayList<>(Arrays.asList(-3, len-2));
            stopCodonsFrame2F = new ArrayList<>(Arrays.asList(-2, len-1));
            stopCodonsFrame3F = new ArrayList<>(Arrays.asList(-1, len));
            stopCodonsFrame1R = new ArrayList<>(Arrays.asList(-3, len-2));
            stopCodonsFrame2R = new ArrayList<>(Arrays.asList(-2, len-1));
            stopCodonsFrame3R = new ArrayList<>(Arrays.asList(-1, len));
        }
        int frame;
        for(int c = 0; c<6; c++){
            for (int locus: orfindeces.get(c)){
                frame = locus%3;
                if (c <3 ){
                    if (frame == 0){
                        stopCodonsFrame1F.add(locus);
                    }
                    else if (frame == 1) {
                        stopCodonsFrame2F.add(locus);
                    }
                    else {
                        stopCodonsFrame3F.add(locus);
                    }
                }
                else {
                    if (frame == 0){
                        stopCodonsFrame1R.add(locus);
                    }
                    else if (frame == 1) {
                        stopCodonsFrame2R.add(locus);
                    }
                    else {
                        stopCodonsFrame3R.add(locus);
                    }
                }
            }
        }
        Collections.sort(stopCodonsFrame1F);
        Collections.sort(stopCodonsFrame2F);
        Collections.sort(stopCodonsFrame3F);
        Collections.sort(stopCodonsFrame1R);
        Collections.sort(stopCodonsFrame2R);
        Collections.sort(stopCodonsFrame3R);
        getORFLengths(stopCodonsFrame1F, 1);
        getORFLengths(stopCodonsFrame2F, 2);
        getORFLengths(stopCodonsFrame3F, 3);
        getORFLengths(stopCodonsFrame1R, 4);
        getORFLengths(stopCodonsFrame2R, 5);
        getORFLengths(stopCodonsFrame3R, 6);
    }

    private static void getORFLengths(ArrayList<Integer> stopCodons, int frame){
        int tempLocus1;
        int tempLocus2 = 0;
        int tempLocus3;
        int lengthOfORF = 0;
        int count = 0;
        for (int loci: stopCodons){
            tempLocus1 = loci;
            if(count == 1){
                tempLocus3 = tempLocus2 - tempLocus1;
                if(((tempLocus3-3)/3) > lengthOfORF){
                    lengthOfORF = ((tempLocus3-3)/3);
                }
            }
            else {
                count = 1;
            }
            tempLocus2 = tempLocus1;
        }
        longestORFs.add(lengthOfORF);
        if (countForObservedORFs <6){
            observedLongestORFs.add(lengthOfORF);
        }
        countForObservedORFs ++;
    }

}
