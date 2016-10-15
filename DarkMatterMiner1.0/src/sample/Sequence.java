package sample;

class Sequence {
    //** Attributes of sequence
    private String name, rawSeq,translatedAA, microsatelliteMotif;
    private int length, rank, numSSRrepeats, sSRstartloci;
    private double gcContent, dinucleotidePValue, triBias1,triBias2,triBias3;
    //**Only best p-values will be kept at the end
    private double orfLenP1 = 1;
    private double orfLenP2 = 1;
    private double orfLenP3 = 1;
    private double orfLenP4 = 1;
    private double orfLenP5 = 1;
    private double orfLenP6 = 1;
    private boolean microsatellite;


        //** initialises Sequence objects with basic information before statistical analysis
    Sequence(String name, String raw, int len){
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
//        triBias1 = t1;
//        triBias2 = t2;
//        triBias3 = t3;
//    }

        //** method used when writing csv file
//    public String getSequence(){
//        return (name + "," + length + "," + gcContent + "," +
//                dinucleotidePValue+ "," + diNucP2+ "," +diNucP3+ "," +
//                orfLenP1+ "," +orfLenP2+ "," + orfLenP3+ "," +
//                orfLenP4+ "," +orfLenP5+ "," +orfLenP6+ "," +
//                triBias1+ "," +triBias2+ "," +triBias3);
//    }


    //** Getters and setters of variables
    public String getName() {
        return name;
    }

    public String getRawSeq() {
        return rawSeq;
    }

    public void setRawSeq(String rawSeq) {
        this.rawSeq = rawSeq;
    }

    public int getLength() {
        return length;
    }

    public void setMicrosatellite(boolean microsatellite) {
        this.microsatellite = microsatellite;
    }

    public void setTranslatedAA(String translatedAA) {
        this.translatedAA = translatedAA;
    }

    public void setMicrosatelliteMotif(String microsatelliteMotif) {
        this.microsatelliteMotif = microsatelliteMotif;
    }

    public void setNumSSRrepeats(int numSSRrepeats) {
        this.numSSRrepeats = numSSRrepeats;
    }

    public void setsSRstartloci(int sSRstartloci) {
        this.sSRstartloci = sSRstartloci;
    }

    public void setGcContent(double gcContent) {
        this.gcContent = gcContent;
    }

    public int getRank() {
        return rank;
    }

    public double getGcContent() {
        return gcContent;
    }

    public double getDinucleotidePValue() {
        return dinucleotidePValue;
    }

    public double getOrfLenP1() {
        return orfLenP1;
    }

    public double getOrfLenP2() {
        return orfLenP2;
    }

    public double getOrfLenP3() {
        return orfLenP3;
    }

    public double getOrfLenP4() {
        return orfLenP4;
    }

    public double getOrfLenP5() {
        return orfLenP5;
    }

    public double getOrfLenP6() {
        return orfLenP6;
    }

    public double getTriBias1() {
        return triBias1;
    }

    public double getTriBias2() {
        return triBias2;
    }

    public double getTriBias3() {
        return triBias3;
    }

    public void setRank(int rank) {
        this.rank = rank;
    }

    public void setTriBias1(double triBias1) {
        this.triBias1 = triBias1;
    }

    public void setTriBias2(double triBias2) {
        this.triBias2 = triBias2;
    }

    public void setTriBias3(double triBias3) {
        this.triBias3 = triBias3;
    }

    public void setDinucleotidePvalue(double dinucleotidePValue) {
        this.dinucleotidePValue = dinucleotidePValue;
    }

    public void setOrfLenP1(double orfLenP1) {
        this.orfLenP1 = orfLenP1;
    }

    public void setOrfLenP2(double orfLenP2) {
        this.orfLenP2 = orfLenP2;
    }

    public void setOrfLenP3(double orfLenP3) {
        this.orfLenP3 = orfLenP3;
    }

    public void setOrfLenP4(double orfLenP4) {
        this.orfLenP4 = orfLenP4;
    }

    public void setOrfLenP5(double orfLenP5) {
        this.orfLenP5 = orfLenP5;
    }

    public void setOrfLenP6(double orfLenP6) {
        this.orfLenP6 = orfLenP6;
    }



}
