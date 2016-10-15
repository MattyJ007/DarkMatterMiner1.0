package sample;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Properties;

class Metagenome {
    //**Stores all data on every sequence
    private static ArrayList<Sequence> sequences = new ArrayList<>();

    static void create(String inOutFolder, double topPercentage, boolean secureRandom){
        File folder = new File(inOutFolder);
        File[] listOfFilesTemp = folder.listFiles();
        assert listOfFilesTemp != null;
        for(File j: listOfFilesTemp){
            //** inputs fas and fasta files into gmato. prevents other files being submitted accidently
            if (j.getName().substring(j.getName().length()-4).equals(".fas") || j.getName().substring(j.getName().length()-4).equals("asta")) {
                //** java running perl through command line to run GMATo
                runGMAT(j);
//                String seqName = "";
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
////                        OrfPValue.pValueMaker(newSequence);
//                        lineCount++;
//                        System.out.println(lineCount+"/"+totSeqNum+"\n------");
//                    }
//                    SortSeqs.bigSort();
//                    OutputPotentialDarkMatterFas.writeOutSeqs(j);
//                    if (Variables.isCsv() == true){
//                        WriteCsv.file(j);
                    }
                    sequences.clear();
                }
                catch (Exception e){
                    System.out.println(e.getMessage()+"\n^^^^^^ Metagenome");
                }
            }
        }
    }

    private static void runGMAT(File input){
        String i = input.getAbsolutePath();
        //** GMATo settings
        String que = "Running statement: \n -r " + DarkMatterMinerUI.getMotifRepeats() + " -m " + DarkMatterMinerUI.getMinMotifLen() + " -x " + DarkMatterMinerUI.getMaxMotifLen() + " -s " + 0 + " -i " + i;
        que = que.replaceAll("Running statement: \n","");
        BufferedReader bufreader;
        Properties pp = System.getProperties();
        String os = pp.getProperty("os.name");
        if ((os.startsWith("win")) || (os.startsWith("Win")))
        {
            String winfilepath = System.getProperty("user.dir");
            System.out.println("\nFile path:\n" + winfilepath + "\n");
            try
            {
                String wincmd = "perl gmat.pl " + que;
                //** runs GMATo through command line
                Process winproc = Runtime.getRuntime().exec("cmd /c" + wincmd);
                if (winproc != null)
                {
                    System.out.println("\nProgram starts running...\n\n");
                    winproc.waitFor();
                    String winline;
                    //** Streams command line output to java console
                    bufreader = new BufferedReader(new InputStreamReader(winproc.getInputStream(), "GB2312"));
                    System.out.println(bufreader.readLine());
                    while (((winline = bufreader.readLine()) != null))
                    {
                        System.out.println(winline + "\n");
                    }
                }
                else
                {
                    System.out.println("Program is finished?!");
                }
            }
            catch (Exception e)
            {
                System.out.println(e.getMessage());
            }
        }
        removeSSRs(input);
    }
    //** SSRs found by GMATo are removed
    private static void removeSSRs(File input){
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
}
