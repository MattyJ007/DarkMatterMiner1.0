package sample;
import java.io.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Properties;
class Metagenome {
    //**Stores all data on every sequence
    private static ArrayList<Sequence> sequences = new ArrayList<>();
    static void create(String inOutFolder, boolean secureRandom){
        File folder = new File(inOutFolder);
        File[] listOfFilesTemp = folder.listFiles();
        assert listOfFilesTemp != null;
        for(File file: listOfFilesTemp){
            //** inputs fas and fasta files into gmato. prevents other files being submitted accidently
            if (file.getName().substring(file.getName().length()-4).equals(".fas") || file.getName().substring(file.getName().length()-4).equals("asta")) {
                //** java running perl through command line to run GMATo
                runGMATo(file);
                removeSSRs(file);
//                String seqName = "";
                getData(secureRandom);
                try{
                    getRankings();
                }
                catch (Exception e){
                    System.out.println(e.getMessage() + " ------ getRankings");
                }
                outputData(file);
                sequences.clear();
            }
        }
    }
    private static void runGMATo(File input){
        String i = input.getAbsolutePath();
        //** GMATo settings
        String que = "Running statement: \n -r " + DarkMatterMinerUI.getMotifRepeats() + " -m " + DarkMatterMinerUI.getMinMotifLen() + " -x " + DarkMatterMinerUI.getMaxMotifLen() + " -s " + 0 + " -i " + i;
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

            int match = 0;
            Sequence raw;
            String seqName = "";
            while (true) {
                String seqLine = bRead2.readLine();
                if (seqLine == null) {
                    break;
                }
                else {
                    if (seqLine.charAt(0) == '>') {
                        match = 0;
                        for (String name : noSSRFilenames) {
                            if (seqLine.equals(name)) {
                                match = 1;
                                break;
                            }
                        }
                        if (match == 0) {
                            seqName = seqLine;
                        }
                    }
                    //** Sequence objects are created and added to Metagenome
                    else if(match == 0 && seqLine.length() > DarkMatterMinerUI.getIgnoreShortSeq() ){
                        raw  = new Sequence(seqName,seqLine,seqLine.length());
                        sequences.add(raw);
                    }

                }
            }
        }
        catch (Exception e){
            System.out.println(e.getMessage());
        }
    }
    private static void getData(boolean secureRandom){
        try {
            //** Counts number of seqs in file for progress output
//                    int lineCount = 0;
//                    int totSeqNum = sequences.size();
            //** Not necessary - just useful to see progress

//                    String seqTot;
            //** Iterates through all valid sequences
            for (Sequence newSequence:sequences) {
//                        Variables.resetDistancesList();
                AnalyseSequence.analyseSequence(newSequence, secureRandom);
//                        lineCount++;
//                        System.out.println(lineCount+"/"+totSeqNum+"\n------");
//                    }
//                    SortSeqs.bigSort();
//                    OutputPotentialDarkMatterFas.writeOutSeqs(j);
            }
            //** clears sequences of previous file - therefore new file starts with unassigned variable.
        }
        catch (Exception e){
            System.out.println(e.getMessage()+"\n^^^^^^ Metagenome getData");
        }
    }
    private static void getRankings()throws Exception{
        int metagenomeLength = sequences.size();
        Collections.sort(sequences, new ItemComparator(ItemComparator.Field.ORFLEN));
        for(int weight = 1; weight<= metagenomeLength; weight++){
            sequences.get(weight-1).setRankOrf(weight);
        }
        Collections.sort(sequences, new ItemComparator(ItemComparator.Field.TRINUC));
        for(int weight = 1; weight<= metagenomeLength; weight++){
            sequences.get(weight-1).setRankTri(weight);
        }
        Collections.sort(sequences, new ItemComparator(ItemComparator.Field.BEST_ORF_TRI_FREQ));
        for(int weight = 1; weight<= metagenomeLength; weight++){
            sequences.get(weight-1).setRankBestORFframeTri(weight);
        }
        Collections.sort(sequences, new ItemComparator(ItemComparator.Field.DINUC));
        for(int weight = 1; weight<= metagenomeLength; weight++){
            sequences.get(weight-1).setRankDi(weight);
            sequences.get(weight-1).setRankTot();
        }
        Collections.sort(sequences, new ItemComparator(ItemComparator.Field.TOTRANK));
    }
    private static void translateSeq(Sequence bestSeq){
        bestSeq.getFrameWithLongestORF();
        bestSeq.getRawSeq();
    }
    private static void outputData(File input){
        String labelString = "Sequence Name\t" +
                "Length of Sequence\t" +
                "GC content\t"+
                "frame with longest ORF\t"+
                "best ORF p-value\t"+
                "ORF Rank\t"+
                "p-value of ORF lengths frame 1\t" +
                "p-value of ORF lengths frame 2\t" +
                "p-value of ORF lengths frame 3\t" +
                "p-value of ORF lengths frame 4\t" +
                "p-value of ORF lengths frame 5\t" +
                "p-value of ORF lengths frame 6\t" +
                "triNuc Freq of Longest ORF frame\t"+
                "Rank trinuc LongestORF\t"+
                "best Trinuc P-value\t"+
                "Rank best Trinuc\t"+
                "trinucPvalue Frame 1\t"+
                "trinucPvalue Frame 2\t"+
                "trinucPvalue Frame 3\t"+
                "best dinuc P-value\t" +
                "Rank Best dinuc\t"+
                "dinucPvalue Frame 1\t"+
                "dinucPvalue Frame 2"+
                "Rank Total"+
                "\n";
        try (
                FileWriter writerCSV = new FileWriter(input+"_DMM.csv");
                FileWriter writerFas = new FileWriter(input+"DMM_BestPotentialSeqs.fas")
        ) {
            writerCSV.write(labelString);
            for(int best = 0; best < ( sequences.size() * DarkMatterMinerUI.getTopResults()); best++){
                writerFas.write(sequences.get(best).getFasSeqInfo() + "\n");
                translateSeq(sequences.get(best));
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
