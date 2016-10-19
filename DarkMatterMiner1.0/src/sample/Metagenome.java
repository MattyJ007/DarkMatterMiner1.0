package sample;
import java.io.*;
import java.util.*;

class Metagenome {
    //**Stores all data on every sequence
    private static ArrayList<Sequence> sequences = new ArrayList<>();
    static void create(String inOutFolder, boolean secureRandom, final ProgressCheck progressFile, final ProgressCheck progressFolder){
        File folder = new File(inOutFolder);
        File[] listOfFilesTemp = folder.listFiles();
        short folderSize;
        short fileCount = 0;
        assert listOfFilesTemp != null;
        for(File file: listOfFilesTemp) {
            if (file.getName().substring(file.getName().length() - 4).equals(".fas") || file.getName().substring(file.getName().length() - 4).equals("asta")) {
                fileCount++;
            }
        }
        folderSize = fileCount;
        fileCount = 0;

        for(File file: listOfFilesTemp){
            fileCount++;
            progressFolder.setProgressNum(fileCount/((double) folderSize));
            //** inputs fas and fasta files into gmato. prevents other files being submitted accidently
            if (file.getName().substring(file.getName().length()-4).equals(".fas") || file.getName().substring(file.getName().length()-4).equals("asta")) {
                //** java running perl through command line to run GMATo
                try {
                    FileReader reader = new FileReader(file);
                    BufferedReader bread = new BufferedReader(reader);
                    if (bread.readLine() == null){
//                        progressString.setProgressString(file.getName() + " is Empty and has been skipped.");
                        continue;
                    }
                    //                progressString.setProgressString("Finding SSRs...");
                    runGMATo(file);
//                progressString.setProgressString("Labelling SSRs...");
                    removeSSRs(file);
                    getData(secureRandom, progressFile);
                    try{
                        getRankings();
                    }
                    catch (Exception e){
                        System.out.println(e.getMessage() + " ------ getRankings");
                    }
                    outputData(file);
                    sequences.clear();
                    //**clears sequences of previous file - therefore new file starts with unassigned variable.
                }
                catch (Exception e){
                    System.out.println(e.getMessage()+ "----- List object issue " + e.getCause());
                }
            }
        }
//        progressString.setProgressString("-------------------------------\n          Complete\n------------------------------");
        System.out.println("Complete!!!");
    }
    private static void runGMATo(File input){
        String i = input.getAbsolutePath();
        //** GMATo settings
        String que = "Running statement: \n -r " + DMMController.getMotifRepeats() + " -m " + DMMController.getMinMotifLen() + " -x " + DMMController.getMaxMotifLen() + " -s " + 0 + " -i " + i;
        que = que.replaceAll("Running statement: \n","");
        BufferedReader bufreader;
        Properties pp = System.getProperties();
        String os = pp.getProperty("os.name");
        if ((os.startsWith("win")) || (os.startsWith("Win")))
        {
//            String winfilepath = System.getProperty("user.dir");
//            System.out.println("\nFile path:\n" + winfilepath + "\n");
            try
            {
                String wincmd = "perl gmat.pl " + que;
                //** runs GMATo through command line
                Process winproc = Runtime.getRuntime().exec("cmd /c" + wincmd);
                if (winproc != null)
                {
//                    System.out.println("\nProgram starts running...\n\n");
                    winproc.waitFor();
                    String winline;
                    //** Streams command line output to java console
                    bufreader = new BufferedReader(new InputStreamReader(winproc.getInputStream(), "GB2312"));
                    while (((winline = bufreader.readLine()) != null))
                    {
                        System.out.println(winline + "\n");
                    }
                }
                else
                {
                    System.out.println("Program is finished?! ------ Something is funny in GMATo");
                }
            }
            catch (Exception e)
            {
                System.out.println(e.getMessage() + "runGMATo");
            }
        }
    }
    private static void removeSSRs(File input){
        //** SSRs found by GMATo are removed
        try(
                FileReader fRead = new FileReader(input+".ssr");
                BufferedReader bRead = new BufferedReader(fRead);
                FileReader fRead2 = new FileReader(input);
                BufferedReader bRead2 = new BufferedReader(fRead2)
        ){
            ArrayList<String> noSSRFilenames = new ArrayList<>();
            while (true){
                String line = bRead.readLine();
                if (line == null){
                    break;
                }
                else{
                    String [] nLine = line.split("\t");
                    if (!nLine[0].equals("Name")){
                        if (!noSSRFilenames.contains(nLine[0])){
                            noSSRFilenames.add(nLine[0]);
                        }
                    }
                }
            }
            boolean match = false;
            Sequence raw;
            String seqName = "";
            while (true) {
                String seqLine = bRead2.readLine();
                if (seqLine == null) {
                    break;
                }
                else {
                    if (seqLine.charAt(0) == '>') {
                        match = false;
                        for (String name : noSSRFilenames) {
                            if (seqLine.equals(name)) {
                                match = true;
                                break;
                            }
                        }
                        if (!match) {
                            seqName = seqLine;
                        }
                    }
                    //** Sequence objects are created and added to Metagenome
                    else if(!match && seqLine.length() > DMMController.getIgnoreShortSeq() ){
                        raw  = new Sequence(seqName,seqLine,seqLine.length());
                        sequences.add(raw);
                    }
                }
            }
        }
        catch (Exception e){
            System.out.println(e.getMessage() + "----------Remove SSRs issue");
        }
    }
    private static void getData(boolean secureRandom, final ProgressCheck progressFile){
        try {
            //** Counts number of seqs in file for progress output
            int lineCount = 0;
            int totSeqNum = sequences.size();
            //** Not necessary - just useful to see progress
            //** Iterates through all valid sequences
            for (Sequence newSequence:sequences) {
                AnalyseSequence.analyseSequence(newSequence, secureRandom);
                lineCount++;
                progressFile.setProgressNum(lineCount/((double) totSeqNum));
            }
        }
        catch (Exception e){
            System.out.println(e.getMessage()+"\n^^^^^^ Metagenome getData");
        }
    }
    private static void getRankings()throws Exception{
//        int metagenomeLength = sequences.size();
//        Collections.sort(sequences, new ItemComparator(ItemComparator.Field.ORFLEN));
//        Collections.reverse(sequences);
//        for(int weight = 1; weight<= metagenomeLength; weight++){
//            sequences.get(weight-1).setRankOrf(weight);
//        }
//        Collections.sort(sequences, new ItemComparator(ItemComparator.Field.TRINUC));
//        Collections.reverse(sequences);
//        for(int weight = 1; weight<= metagenomeLength; weight++){
//            sequences.get(weight-1).setRankTri(weight);
//        }
//        Collections.sort(sequences, new ItemComparator(ItemComparator.Field.DINUC));
//        Collections.reverse(sequences);
//        for(int weight = 1; weight<= metagenomeLength; weight++){
//            sequences.get(weight-1).setRankDi(weight);
//        }
//        Collections.sort(sequences, new ItemComparator(ItemComparator.Field.BEST_ORF_TRI_FREQ));
//        Collections.reverse(sequences);
//        for(int weight = 1; weight<= metagenomeLength; weight++){
//            sequences.get(weight-1).setRankBestORFframeTri(weight);
//            sequences.get(weight-1).setRankTot();
//        }
        for (Sequence newSeq: sequences){
            newSeq.setRankTot();
        }
        Collections.sort(sequences, new ItemComparator(ItemComparator.Field.TOTRANK));
//        Collections.reverse(sequences);
    }
    private static final Map<String, Character> codonsMap;
    static {
        codonsMap = new HashMap<>();
        codonsMap.put("UUU", 'F');
        codonsMap.put("UUC", 'F');
        codonsMap.put("UUA", 'L');
        codonsMap.put("TTT", 'F');
        codonsMap.put("TTC", 'F');
        codonsMap.put("TTA", 'L');
        codonsMap.put("TTG", 'L');
        codonsMap.put("TCT", 'S');
        codonsMap.put("TCC", 'S');
        codonsMap.put("TCA", 'S');
        codonsMap.put("TCG", 'S');
        codonsMap.put("TAT", 'Y');
        codonsMap.put("TAC", 'Y');
        codonsMap.put("TAA", '_');
        codonsMap.put("TAG", '_');
        codonsMap.put("TGT", 'C');
        codonsMap.put("TGC", 'C');
        codonsMap.put("TGA", '_');
        codonsMap.put("TGG", 'W');
        codonsMap.put("CTT", 'L');
        codonsMap.put("CTC", 'L');
        codonsMap.put("CTA", 'L');
        codonsMap.put("CTG", 'L');
        codonsMap.put("CCT", 'P');
        codonsMap.put("CCC", 'P');
        codonsMap.put("CCA", 'P');
        codonsMap.put("CCG", 'P');
        codonsMap.put("CAT", 'H');
        codonsMap.put("CAC", 'H');
        codonsMap.put("CAA", 'Q');
        codonsMap.put("CAG", 'Q');
        codonsMap.put("CGT", 'R');
        codonsMap.put("CGC", 'R');
        codonsMap.put("CGA", 'R');
        codonsMap.put("CGG", 'R');
        codonsMap.put("ATT", 'I');
        codonsMap.put("ATC", 'I');
        codonsMap.put("ATA", 'I');
        codonsMap.put("ATG", 'M');
        codonsMap.put("ACT", 'T');
        codonsMap.put("ACC", 'T');
        codonsMap.put("ACA", 'T');
        codonsMap.put("ACG", 'T');
        codonsMap.put("AAT", 'N');
        codonsMap.put("AAC", 'N');
        codonsMap.put("AAA", 'K');
        codonsMap.put("AAG", 'K');
        codonsMap.put("AGT", 'S');
        codonsMap.put("AGC", 'S');
        codonsMap.put("AGA", 'R');
        codonsMap.put("AGG", 'R');
        codonsMap.put("GTT", 'V');
        codonsMap.put("GTC", 'V');
        codonsMap.put("GTA", 'V');
        codonsMap.put("GTG", 'V');
        codonsMap.put("GCT", 'A');
        codonsMap.put("GCC", 'A');
        codonsMap.put("GCA", 'A');
        codonsMap.put("GCG", 'A');
        codonsMap.put("GAT", 'D');
        codonsMap.put("GAC", 'D');
        codonsMap.put("GAA", 'E');
        codonsMap.put("GAG", 'E');
        codonsMap.put("GGT", 'G');
        codonsMap.put("GGC", 'G');
        codonsMap.put("GGA", 'G');
        codonsMap.put("GGG", 'G');
    }
    private static void translate(Sequence seq) {
        //temp is a multiple of 3
        // every 3 characters correspond to a valid codon?
        int frame;
        if(seq.getFrameWithLongestORF()==2){
            frame =0;
        }
        else if (seq.getFrameWithLongestORF()==0){
            //** TEst annomalies??? this frame produces the best aa
            frame = 0;
        }
        else if (seq.getFrameWithLongestORF()==1){
            frame = 2;
        }
        else if(seq.getFrameWithLongestORF()==3){
            frame = 0;
            compliment(seq);
        }
        else if(seq.getFrameWithLongestORF()==4){
            frame = 1;
            compliment(seq);
        }
        else{
            frame = 2;
            compliment(seq);
        }
        String temp = seq.getRawSeq();
        StringBuilder finalreturn = new StringBuilder();
        String codon;
        for (int i = frame; i < seq.getLength() - 2; i+=3) {
            codon = temp.substring(i, i+3);
            if(codon.contains("N") || codon.contains("n")){
                finalreturn.append("^");
            }
            else {
                finalreturn.append(codonsMap.get(codon));
            }
        }
        seq.setTranslatedAA(finalreturn.toString());
    }
    private static final Map<String, String> mRNA = new HashMap<>();
    static {
        mRNA.put("T", "A");
        mRNA.put("A", "U");
        mRNA.put("G", "C");
        mRNA.put("C", "G");

    }
    private static void transcribe(Sequence seq){
        int frame;
        if(seq.getFrameWithLongestORF()==0){
            frame =0;
        }
        else if (seq.getFrameWithLongestORF()==1){
            frame = 1;
        }
        else if (seq.getFrameWithLongestORF()==2){
            frame = 2;
        }
        else if(seq.getFrameWithLongestORF()==3){
            frame = 0;
            compliment(seq);
        }
        else if(seq.getFrameWithLongestORF()==4){
            frame = 1;
            compliment(seq);
        }
        else{
            frame = 2;
            compliment(seq);
        }
        String temp = seq.getRawSeq();
        StringBuilder finalreturn = new StringBuilder();
        String codon;
        for (int i = frame; i < seq.getLength(); i++) {
            codon = temp.substring(i, i+1);
            if(codon.contains("N")){
                finalreturn.append("NNN");
            }
            else{
                finalreturn.append(mRNA.get(codon));
            }
        }
        seq.setTranscribedmRNA(finalreturn.toString());
    }
    private static final Map<String,String> compliments = new HashMap<>();
    static {
        compliments.put("A","T");
        compliments.put("T","A");
        compliments.put("C","G");
        compliments.put("G","C");
        compliments.put("N","N");
    }
    private static String compliment(Sequence seq){
        String temp = seq.getRawSeq();
        StringBuilder complimentrayStrand = new StringBuilder();
        String nucleotide;
        for (int i = 0; i < temp.length(); i++) {
            nucleotide = temp.substring(i, i+1);
            complimentrayStrand.append(compliments.get(nucleotide));
        }
        return complimentrayStrand.toString();
    }
    private static void outputData(File input){
        String labelString = "Name\t" +
                "Length\t" +
                "GC Content\t"+
                "Frame with Longest ORF\t"+
                "Longest ORF P-value\t"+
//                "Longest ORF Rank\t"+
                "Tri-nuc Freq of Longest ORF Frame P-value\t"+
                "Best Tri-nuc P-value\t"+
//                "Rank Tri-nuc Longest ORF\t"+
//                "Rank Best Tri-nuc\t"+
                "Best Di-nuc P-value\t" +
//                "Rank Best Di-nuc\t"+
                "P-value of ORF lengths frame 1\t" +
                "P-value of ORF lengths frame 2\t" +
                "P-value of ORF lengths frame 3\t" +
                "P-value of ORF lengths frame 4\t" +
                "P-value of ORF lengths frame 5\t" +
                "P-value of ORF lengths frame 6\t" +
                "Tri-nuc P-value Frame 1\t"+
                "Tri-nuc P-value Frame 2\t"+
                "Tri-nuc P-value Frame 3\t"+
                "Di-nuc P-value Frame 1\t"+
                "Di-nuc P-value Frame 2\t"+
                "Rank Total\t"+
                "Amino Acid Seq of Longest ORF Frame\t"+
                "mRNA seq of Longest ORF Frame\t"+
//                "Original seq"+
                "\n";
        try (
                FileWriter writerCSV = new FileWriter(input+"_DMM.csv")
//                FileWriter writerFas = new FileWriter(input+"DMM_BestPotentialSeqs.fas")
        ) {
            writerCSV.write(labelString);
            for(int best = 0; best < ( sequences.size() * DMMController.getTopResults()); best++){
//                writerFas.write(sequences.get(best).getFasSeqInfo());
                translate(sequences.get(best));
                transcribe(sequences.get(best));
            }
            for (Sequence seq: sequences){
                writerCSV.write(seq.getSequence()+"\n");
            }
        }
        catch (Exception e){
            System.out.println(e.getMessage() + " - - - - - outputData");
        }
    }
}
