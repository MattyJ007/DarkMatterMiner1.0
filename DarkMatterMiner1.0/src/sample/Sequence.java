package sample;
import java.util.ArrayList;
import java.util.Collections;
class Sequence {
    //** Attributes of sequence
    private String name, rawSeq,translatedAA, transcribedmRNA, microsatelliteMotif;
    private int length, rankTri, rankBestORFframeTri,rankDi, rankOrf, rankTot, numSSRrepeats, sSRstartloci,frameWithLongestORF;
    private double gcContent, dinucleotidePValue, trinucelotidePValue, orfLengthPValue, trinuc1, trinuc2, trinuc3, triNucFreqOfLongestORFframe, dinuc1, dinuc2, orfLenP1, orfLenP2,orfLenP3,orfLenP4,orfLenP5,orfLenP6;
//    private boolean microsatellite;

    Sequence(String name, String raw, int len){
        //** initialises Sequence objects with basic information before statistical analysis
        this.name = name;
        rawSeq = raw;
        length = len;

    }
        //** Inputs missing data but this method is most likely going to become redundant
//    Sequence(String n, int l, double gc, double d1, double d2, double d3, double o1, double o2, double o3, double o4, double o5, double o6, double t1, double t2, double t3){
//        name = n;
//        length = l;
//        gcContent = gc;
//        dinucleotidePValue = d1;
//        diNucP2 = d2;
//        diNucP3 = d3;
//        orfLenP1 = o1;
//        orfLenP2 = o2;
//        orfLenP3 = o3;
//        orfLenP4 = o4;
//        orfLenP5 = o5;
//        orfLenP6 = o6;
//        trinucelotidePValue = t1;
//        triBias2 = t2;
//        triBias3 = t3;
//    }

        //** method used when writing csv file
    String getSequence(){
        return (name + "\t" + length + "\t" + gcContent + "\t" +
                frameWithLongestORF + "\t"+
                orfLengthPValue + "\t"+ rankOrf + "\t"+
                orfLenP1+ "\t" +orfLenP2+ "\t" + orfLenP3+ "\t" +
                orfLenP4+ "\t" +orfLenP5+ "\t" +orfLenP6+ "\t" +
                triNucFreqOfLongestORFframe + "\t"+ rankBestORFframeTri +"\t"+
                trinucelotidePValue+ "\t"  + rankTri + "\t"+
                trinuc1+ "\t"  + trinuc2+ "\t"  +trinuc3+ "\t"  +
                dinucleotidePValue+ "\t" + rankDi + "\t"+
                dinuc1+ "\t"  +dinuc2 + "\t"+ rankTot+"\t"+
                translatedAA+"\t"+transcribedmRNA
//                +"\t"+rawSeq
        );
    }
    String getFasSeqInfo(){
        return (name+"\n"+rawSeq+"\n");
    }


    //** Getters and setters of variables
    String getRawSeq() {
        return rawSeq;
    }

    int getFrameWithLongestORF() {
        return frameWithLongestORF;
    }
    int getLength(){
        return length;
    }

    double getDinucleotidePValue() {
        return dinucleotidePValue;
    }

    double getTrinucelotidePValue() {
        return trinucelotidePValue;
    }

    double getOrfLengthPValue() {
        return orfLengthPValue;
    }

    void setRankTri(int rankTri) {
        this.rankTri = rankTri;
    }

    void setRankBestORFframeTri(int rankBestORFframeTri) {
        this.rankBestORFframeTri = rankBestORFframeTri;
    }

    void setRankDi(int rankDi) {
        this.rankDi = rankDi;
    }

    void setRankOrf(int rankOrf) {
        this.rankOrf = rankOrf;
    }

    void setRankTot() {
        this.rankTot = rankDi + rankTri + rankOrf;
    }

    int getRankTot() {
        return rankTot;
    }

    //    public void setMicrosatellite(boolean microsatellite) {
//        this.microsatellite = microsatellite;
//    }
//
    void setTranslatedAA(String translatedAA) {
        this.translatedAA = translatedAA;
    }

    void setTranscribedmRNA(String transcribedmRNA) {
         this.transcribedmRNA = transcribedmRNA;
    }

    //
//    public void setMicrosatelliteMotif(String microsatelliteMotif) {
//        this.microsatelliteMotif = microsatelliteMotif;
//    }
//
//    public void setNumSSRrepeats(int numSSRrepeats) {
//        this.numSSRrepeats = numSSRrepeats;
//    }
//
//    public void setsSRstartloci(int sSRstartloci) {
//        this.sSRstartloci = sSRstartloci;
//    }
    void setGcContent(double gcContent) {
        this.gcContent = gcContent;
    }
    void setOrfPvalues(ArrayList<Double> pValues){
        frameWithLongestORF = pValues.indexOf((Collections.min(pValues)));
        orfLengthPValue = Collections.min(pValues);
        orfLenP1 = pValues.get(0);
        orfLenP2 = pValues.get(1);
        orfLenP3 = pValues.get(2);
        orfLenP4 = pValues.get(3);
        orfLenP5 = pValues.get(4);
        orfLenP6 = pValues.get(5);
    }
    void setMotifPValues(ArrayList<Double> pValues1, ArrayList<Double> pValues2){
        trinucelotidePValue = Collections.min(pValues1);
        trinuc1 = pValues1.get(0);
        trinuc2 = pValues1.get(1);
        trinuc3 = pValues1.get(2);
        dinucleotidePValue = Collections.min(pValues2);
        dinuc1 = pValues2.get(0);
        dinuc2 = pValues2.get(1);
        if (frameWithLongestORF%3 == 0){
            triNucFreqOfLongestORFframe = trinuc1;
        }
        if (frameWithLongestORF%3 == 1){
            triNucFreqOfLongestORFframe = trinuc2;
        }
        if (frameWithLongestORF%3 == 2){
            triNucFreqOfLongestORFframe = trinuc3;
        }
    }
}
