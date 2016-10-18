package sample;
import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Random;
import java.util.stream.Collectors;

import static jdk.nashorn.internal.objects.NativeString.indexOf;
class AnalyseSequence {
    static void analyseSequence(Sequence newSeq, boolean secureRandom){
        resetArrayListsAndCount();
        setTotalExpectedFrequencies();
        newSeq.setGcContent(gCcount((newSeq.getRawSeq())));
        generatePermutations(newSeq.getRawSeq(), secureRandom);
        chiSquare();
        getPvalues(newSeq);
    }
    private static Short[] totalExpectedFrequenciesTri = new Short[64];
    private static Short[] totalExpectedFrequenciesDi = new Short[16];
    private static ArrayList<ArrayList<Short>> allPermutationFreqArrays;
    private static ArrayList<Short> trinucFreqFrame1Temp;
    private static ArrayList<Short> trinucFreqFrame2Temp;
    private static ArrayList<Short> trinucFreqFrame3Temp;
    private static ArrayList<Short> dinucFreqFrame1Temp;
    private static ArrayList<Short> dinucFreqFrame2Temp;
    private static ArrayList<Float> observedChiValuesT;
    private static ArrayList<Float> expectedChiValuesT;
    private static ArrayList<Float> observedChiValuesD;
    private static ArrayList<Float> expectedChiValuesD;
    private static ArrayList<Short> longestORFs;
    private static ArrayList<Short> observedLongestORFs;
    private static byte countForObservedORFs;

    private static void resetArrayListsAndCount(){
        allPermutationFreqArrays = new ArrayList<>();
        observedChiValuesT = new ArrayList<>();
        observedChiValuesD = new ArrayList<>();
        expectedChiValuesT = new ArrayList<>();
        expectedChiValuesD = new ArrayList<>();
        countForObservedORFs = 0;
        observedLongestORFs = new ArrayList<>();
        longestORFs = new ArrayList<>();
    }
    private static void setTotalExpectedFrequencies(){
        for (byte j = 0; j<64;j++){
            totalExpectedFrequenciesTri[j] = 0;
        }
        for (byte j = 0; j<16;j++){
            totalExpectedFrequenciesDi[j] = 0;
        }
    }
    private static void clearFrameArrays(){
        trinucFreqFrame1Temp = new ArrayList<>();
        trinucFreqFrame2Temp = new ArrayList<>();
        trinucFreqFrame3Temp = new ArrayList<>();
        dinucFreqFrame1Temp = new ArrayList<>();
        dinucFreqFrame2Temp = new ArrayList<>();
    }

    private static float gCcount(String seq){
        //** gets GC content of each sequence
        short countGC = (short) (countMatches(seq, "G") + countMatches(seq, "C"));
        return (float)countGC/seq.length();
    }

    private static short countMatches(String str, String sub) {
        if (str.isEmpty() || sub.isEmpty()) {return -1;}
        short count = 0;
        short idx = 0;
        while ((idx = (short) indexOf(str, sub, idx)) != -1) {
            count++;
            idx += sub.length();
        }
        return count;
    }

    private static void generatePermutations(String newSeq, boolean secureRandom){
        short permutations = (short) DMMController.getPermutations();
        SecureRandom sRand;
        Random nRand;
        short randomNum;
        short len;
        String temp;
        String randSeq;
        String[] seq;
        getFrequency(newSeq);
        if (secureRandom) {
            sRand = new SecureRandom();
            for (short j = 0; j< permutations; j++){
                len = (short) newSeq.length();
                seq = newSeq.split("");
                for (short i=0; i<len; i++){
                    randomNum = (short) sRand.nextInt(len);
                    temp = seq[i];
                    seq[i] = seq[randomNum];
                    seq[randomNum] = temp;
                }
                randSeq = String.join("", (CharSequence[]) seq);
                getFrequency(randSeq);
            }
        }
        else {
            nRand = new Random();
            for (short j = 0; j< permutations; j++){
                len = (short) newSeq.length();
                seq = newSeq.split("");
                for (short i=0; i<len; i++){
                    randomNum = (short) nRand.nextInt(len);
                    temp = seq[i];
                    seq[i] = seq[randomNum];
                    seq[randomNum] = temp;
                }
                randSeq = String.join("", (CharSequence[]) seq);
                getFrequency(randSeq);
            }
        }
    }

    private static void getFrequency(String seq){
        clearFrameArrays();
        ArrayList<ArrayList<Short>> orfindeces = new ArrayList<>();
        ArrayList<ArrayList<Short>> unsortedMotifFrequencies = new ArrayList<>();
        String motifString = "AAA,AAT,AAG,AAC,ATA,ATT,ATG,ATC,AGA,AGT,AGG,AGC,ACA,ACT,ACG,ACC,TAA,TAT,TAG,TAC,TTA,TTT,TTG,TTC,TGA,TGT,TGG,TGC,TCA,TCT,TCG,TCC,GAA,GAT,GAG,GAC,GTA,GTT,GTG,GTC,GGA,GGT,GGG,GGC,GCA,GCT,GCG,GCC,CAA,CAT,CAG,CAC,CTA,CTT,CTG,CTC,CGA,CGT,CGG,CGC,CCA,CCT,CCG,CCC,AA,AT,AG,AC,TA,TT,TG,TC,GA,GT,GG,GC,CA,CT,CG,CC";
        Short[] stopCodons = {16,18,20,24,28,52};
        String[] motifList = motifString.split(",");
        for (String motif: motifList){
            ArrayList<Short> motifIndeces = new ArrayList<>();
            short index = (short) seq.indexOf(motif);
            while(index >= 0) {
                motifIndeces.add(index);
                index = (short) seq.indexOf(motif, index+1);
            }
            unsortedMotifFrequencies.add(motifIndeces);
        }
        for (short d: stopCodons){
            orfindeces.add(unsortedMotifFrequencies.get(d));
        }
        getORFLoci(orfindeces, (short) seq.length());
        sortFrequenciesIntoFrames(unsortedMotifFrequencies);
    }

    private static void sortFrequenciesIntoFrames(ArrayList<ArrayList<Short>> unsortedMotifFrequencies){
        byte frame;
        for (short j = 0; j<80; j++){
            short motifFrameT1 = 0;
            short motifFrameT2 = 0;
            short motifFrameT3 = 0;
            short motifFrameD1 = 0;
            short motifFrameD2 = 0;
            if(j<64){
                for(short motifIndex: unsortedMotifFrequencies.get(j)) {
                    frame = (byte) (motifIndex % 3);
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
                for(short motifIndex: unsortedMotifFrequencies.get(j)) {
                    frame = (byte) (motifIndex % 2);
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
        for(byte i = 0; i<64; i++) {
             (totalExpectedFrequenciesTri[i]) = (short) ((totalExpectedFrequenciesTri[i]) +trinucFreqFrame1Temp.get(i) + trinucFreqFrame2Temp.get(i) + trinucFreqFrame3Temp.get(i));
        }
        for(byte i = 0; i<16; i++) {
            totalExpectedFrequenciesDi[i] = (short) (totalExpectedFrequenciesDi[i] + dinucFreqFrame1Temp.get(i) + dinucFreqFrame2Temp.get(i));
        }
        ArrayList<Short> trinucFreqFrame1 = new ArrayList<>(trinucFreqFrame1Temp);
        ArrayList<Short> trinucFreqFrame2 = new ArrayList<>(trinucFreqFrame2Temp);
        ArrayList<Short> trinucFreqFrame3 = new ArrayList<>(trinucFreqFrame3Temp);
        ArrayList<Short> dinucFreqFrame1 = new ArrayList<>(dinucFreqFrame1Temp);
        ArrayList<Short> dinucFreqFrame2 = new ArrayList<>(dinucFreqFrame2Temp);
        allPermutationFreqArrays.add(trinucFreqFrame1);
        allPermutationFreqArrays.add(trinucFreqFrame2);
        allPermutationFreqArrays.add(trinucFreqFrame3);
        allPermutationFreqArrays.add(dinucFreqFrame1);
        allPermutationFreqArrays.add(dinucFreqFrame2);
    }

    private static void chiSquare(){
        float exp;
        short count = 0;
        for (ArrayList<Short> observedArray :allPermutationFreqArrays) {
            float chiT = 0;
            if (observedArray.size() == 64){
                for(short n = 0; n<64;n++){
                    exp = (totalExpectedFrequenciesTri[n]-observedArray.get(n))/((float) (DMMController.getPermutations()*3)+2);
//                    System.out.println(exp + " = " + totalExpectedFrequenciesTri[n]+" - "+ observedArray.get(n)+" / "+ DarkMatterMinerMain.getPermutations()+" * 3" + " + 2");
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
            float chiD = 0;
            if(observedArray.size() == 16){
                for(short n = 0; n<16;n++){
                    exp = (totalExpectedFrequenciesDi[n]-observedArray.get(n))/((float) (DMMController.getPermutations()*2)+1);
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
    }

    private static void getPvalues(Sequence newSeq){
        ArrayList<Float> orfPValues = new ArrayList<>();
        Collections.sort(longestORFs);
        orfPValues.addAll(observedLongestORFs.stream().map(h -> (1 - (longestORFs.indexOf(h) / (float) longestORFs.size()))).collect(Collectors.toList()));
        newSeq.setOrfPvalues(orfPValues);
        ArrayList<Float> motifFrequenciesT = new ArrayList<>();
        ArrayList<Float> motifFrequenciesD = new ArrayList<>();
        motifFrequenciesT.addAll(observedChiValuesT.stream().map(h -> (1 - (expectedChiValuesT.indexOf(h) / (float) expectedChiValuesT.size()))).collect(Collectors.toList()));
        motifFrequenciesD.addAll(observedChiValuesD.stream().map(h -> (1 - (expectedChiValuesD.indexOf(h) / (float) expectedChiValuesD.size()))).collect(Collectors.toList()));
        newSeq.setMotifPValues(motifFrequenciesT, motifFrequenciesD);
    }

    private static void getORFLoci(ArrayList<ArrayList<Short>> orfindeces, short len){
        //** values in arrays indicate the starting loci of stop codons
        short len1 = (short) (len-1);
        short len2 = (short) (len-2);
        byte numberOfCodons = (byte) (len%3);
        ArrayList<Short> stopCodonsFrame1F;
        ArrayList<Short> stopCodonsFrame2F;
        ArrayList<Short> stopCodonsFrame3F;
        ArrayList<Short> stopCodonsFrame1R;
        ArrayList<Short> stopCodonsFrame2R;
        ArrayList<Short> stopCodonsFrame3R;
        if (numberOfCodons == 0) {
            stopCodonsFrame1F = new ArrayList<>(Arrays.asList((short)-3,len));
            stopCodonsFrame2F = new ArrayList<>(Arrays.asList((short)-2,len2));
            stopCodonsFrame3F = new ArrayList<>(Arrays.asList((short)-1,len1));
            stopCodonsFrame1R = new ArrayList<>(Arrays.asList((short)-3,len));
            stopCodonsFrame2R = new ArrayList<>(Arrays.asList((short)-2,len2));
            stopCodonsFrame3R = new ArrayList<>(Arrays.asList((short)-1,len1));
        }
        else if(numberOfCodons ==1){
            stopCodonsFrame1F = new ArrayList<>(Arrays.asList((short)-3,len1));
            stopCodonsFrame2F = new ArrayList<>(Arrays.asList((short)-2, len));
            stopCodonsFrame3F = new ArrayList<>(Arrays.asList((short)-1, len2));
            stopCodonsFrame1R = new ArrayList<>(Arrays.asList((short)-3, len1));
            stopCodonsFrame2R = new ArrayList<>(Arrays.asList((short)-2, len));
            stopCodonsFrame3R = new ArrayList<>(Arrays.asList((short)-1, len2));
        }
        else{
            stopCodonsFrame1F = new ArrayList<>(Arrays.asList((short)-3, len2));
            stopCodonsFrame2F = new ArrayList<>(Arrays.asList((short)-2, len1));
            stopCodonsFrame3F = new ArrayList<>(Arrays.asList((short)-1, len));
            stopCodonsFrame1R = new ArrayList<>(Arrays.asList((short)-3, len2));
            stopCodonsFrame2R = new ArrayList<>(Arrays.asList((short)-2, len1));
            stopCodonsFrame3R = new ArrayList<>(Arrays.asList((short)-1, len));
        }
        byte frame;
        for(short c = 0; c<6; c++){
            for (short locus: orfindeces.get(c)){
                frame = (byte) (locus%3);
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
        getORFLengths(stopCodonsFrame1F);
        getORFLengths(stopCodonsFrame2F);
        getORFLengths(stopCodonsFrame3F);
        getORFLengths(stopCodonsFrame1R);
        getORFLengths(stopCodonsFrame2R);
        getORFLengths(stopCodonsFrame3R);
    }

    private static void getORFLengths(ArrayList<Short> stopCodons){
        short tempLocus1;
        short tempLocus2 = 0;
        short tempLocus3;
        short lengthOfORF = 0;
        boolean prepareTempVariables = false;
        for (short loci: stopCodons){
            tempLocus1 = loci;
            if(prepareTempVariables){
                tempLocus3 = (short) (tempLocus1 - tempLocus2);
                if(((tempLocus3-3)/3) > lengthOfORF){
                    lengthOfORF = (short) ((tempLocus3-3)/3);
                }
            }
            else {
                prepareTempVariables = true;
            }
            tempLocus2 = tempLocus1;
        }
        longestORFs.add(lengthOfORF);
        if (countForObservedORFs <6){
            observedLongestORFs.add(lengthOfORF);
            countForObservedORFs ++;
        }
    }
}
