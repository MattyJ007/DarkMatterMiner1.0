package sample;
import java.util.ArrayList;
import java.util.Collections;
class Sequence {
    //** Attributes of sequence
    private String name, rawSeq,translatedAA, transcribedmRNA;
//    String microsatelliteMotif;
    private int length;
    private byte frameWithLongestORF;
    private float rankTot;
//            , rankTri,rankDi, rankOrf, rankBestORFframeTri;
//    short numSSRrepeats, sSRstartloci,
    private float gcContent, dinucleotidePValue, trinucelotidePValue, orfLengthPValue, trinuc1, trinuc2, trinuc3, triNucFreqOfLongestORFframe, dinuc1, dinuc2, orfLenP1, orfLenP2,orfLenP3,orfLenP4,orfLenP5,orfLenP6;
//    private boolean microsatellite;

    Sequence(String name, String raw, int len){
        //** initialises Sequence objects with basic information before statistical analysis
        this.name = name;
        rawSeq = raw;
        length = len;

    }
    //** method used when writing csv file
    String getSequence(){
        return (name + "\t" + length + "\t" + gcContent + "\t" +
                frameWithLongestORF + "\t"+
                orfLengthPValue + "\t"+triNucFreqOfLongestORFframe + "\t"+
                trinucelotidePValue + "\t"  +
//                rankBestORFframeTri +"\t"+ rankTri + "\t"+ rankOrf + "\t"+rankDi + "\t"+
                dinucleotidePValue + "\t" +
                orfLenP1 + "\t" +orfLenP2+ "\t" + orfLenP3+ "\t" + orfLenP4+ "\t" +orfLenP5+ "\t" +orfLenP6+ "\t" +
                trinuc1 + "\t"  + trinuc2+ "\t"  +trinuc3+ "\t"  +
                dinuc1 + "\t" + dinuc2 + "\t"+ rankTot+"\t"+
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

    short getFrameWithLongestORF() {
        return frameWithLongestORF;
    }

    int getLength(){
        return length;
    }

//    float getDinucleotidePValue() {
//        return dinucleotidePValue;
//    }
//
//    float getTrinucelotidePValue() {
//        return trinucelotidePValue;
//    }
//
//    float getOrfLengthPValue() {
//        return orfLengthPValue;
//    }

//    void setRankTri(int rankTri) {
//        this.rankTri = rankTri;
//    }
//
//    void setRankBestORFframeTri( int rankBestORFframeTri) {
//        this.rankBestORFframeTri = rankBestORFframeTri;
//    }
//
//    void setRankDi(int rankDi) {
//        this.rankDi = rankDi;
//    }
//
//    void setRankOrf(int rankOrf) {
//        this.rankOrf = rankOrf;
//    }
//
    void setRankTot(){
        rankTot = (triNucFreqOfLongestORFframe + dinucleotidePValue+orfLengthPValue);
    }

//    float getRankBestORFframeTri(){
//        return rankBestORFframeTri;
//    }

    float getRankTot() {
        return rankTot;
    }
    //    public void setMicrosatellite(boolean microsatellite) {
//        this.microsatellite = microsatellite;
//    }
    void setTranslatedAA(String translatedAA) {
        this.translatedAA = translatedAA;
    }
    void setTranscribedmRNA(String transcribedmRNA) {
         this.transcribedmRNA = transcribedmRNA;
    }
//    public void setMicrosatelliteMotif(String microsatelliteMotif) {
//        this.microsatelliteMotif = microsatelliteMotif;
//    }
//    public void setNumSSRrepeats(short numSSRrepeats) {
//        this.numSSRrepeats = numSSRrepeats;
//    }
//    public void setsSRstartloci(short sSRstartloci) {
//        this.sSRstartloci = sSRstartloci;
//    }
    void setGcContent(float gcContent) {
        this.gcContent = gcContent;
    }
    void setOrfPvalues(ArrayList<Float> pValues){
        frameWithLongestORF = (byte) pValues.indexOf((Collections.min(pValues)));
        orfLengthPValue = Collections.min(pValues);
        orfLenP1 = pValues.get(0);
        orfLenP2 = pValues.get(1);
        orfLenP3 = pValues.get(2);
        orfLenP4 = pValues.get(3);
        orfLenP5 = pValues.get(4);
        orfLenP6 = pValues.get(5);
    }
    void setMotifPValues(ArrayList<Float> pValues1, ArrayList<Float> pValues2){
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