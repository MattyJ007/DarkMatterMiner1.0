package sample;
import javafx.fxml.FXML;
import javafx.scene.control.Button;
import javafx.scene.control.ProgressBar;
import javafx.scene.control.RadioButton;
import javafx.scene.control.TextArea;
import java.io.BufferedReader;
import java.io.File;
import java.io.InputStreamReader;
public class DMMController {
    private static boolean start; //** variable that prevents launching Metagenome if User settings aren't correct
    private static boolean secureRandom; //** Allows user to use either Random() or SecureRandom()
    private static double topResults = 0.005; //** Allows user to choose number of seqs for end fasta file
    private static int motifRepeats = 4;
    private static int permutations = 1000;
    private static int minMotifLen = 2;
    private static int maxMotifLen = 10;
    private static int ignoreShortSeq = 20;
    private static String inputFolder;
    private final ProgressNumber progressFile = new ProgressNumber();
    private final ProgressNumber progressFolder = new ProgressNumber();
    //** Getters and Setters of above variables^^^
    @FXML
    public Button runButt;
    public Button helpButt;
    public RadioButton secureRandomx; //** Allows user to use either Random() or SecureRandom()
    public TextArea topResultsx; //** Allows user to choose number of seqs for end fasta file
    public TextArea motifRepeatsx;
    public TextArea permutationsx;
    public TextArea minMotifLenx;
    public TextArea maxMotifLenx;
    public TextArea ignoreShortSeqx;
    public TextArea inputFolderx;
    public TextArea userOutput;
    public ProgressBar fileProgress;
    public ProgressBar folderProgress;

    public void runDMM(){
        userOutput.setText("Checking entries...");
        secureRandomx.isSelected();
        userOutput.setText(checkUserInput(inputFolderx.getText(), motifRepeatsx.getText(),
                minMotifLenx.getText(), maxMotifLenx.getText(), permutationsx.getText(),
                secureRandomx.selectedProperty().getValue(), ignoreShortSeqx.getText(),
                topResultsx.getText()));
        if (isStart()){
            Thread t1 = new Thread(new RunClass());
            t1.start();
            runButt.setDisable(true);
        }
    }

    private class RunClass implements Runnable{
        @Override
        public void run() {
            try{
                fileProgress.progressProperty().bind(progressFile.progressNumProperty());
                folderProgress.progressProperty().bind(progressFolder.progressNumProperty());
                Metagenome.create(getInputFolder(),isSecureRandom(), progressFile, progressFolder);
            }
            catch (Exception e){
                System.out.println(e.getMessage() + "Run issue");
            }
        }
    }

    public void helpDMM(){
        userOutput.setText("HELP!!");
    }

    private static void setStart(boolean startT) {
        start = startT;
    }

    private static boolean isStart() {
        return start;
    }

    private static boolean isSecureRandom() {
        return secureRandom;
    }

    private static void setSecureRandom(boolean secureRandomT) {
        secureRandom = secureRandomT;
    }

    static double getTopResults() {
        return topResults;
    }

    private static void setTopResults(double topResultsT) {
        topResults = topResultsT;
    }

    static int getMotifRepeats() {
        return motifRepeats;
    }

    private static void setMotifRepeats(int motifRepeatsT) {
        motifRepeats = motifRepeatsT;
    }

    static int getPermutations() {
        return permutations;
    }

    private static void setPermutations(int permutationsT) {
        permutations = permutationsT;
    }

    static int getMinMotifLen() {
        return minMotifLen;
    }

    private static void setMinMotifLen(int minMotifLenT) {
        minMotifLen = minMotifLenT;
    }

    static int getMaxMotifLen() {
        return maxMotifLen;
    }

    private static void setMaxMotifLen(int maxMotifLenT) {
        maxMotifLen = maxMotifLenT;
    }

    static int getIgnoreShortSeq() {
        return ignoreShortSeq;
    }

    private static void setIgnoreShortSeq(int ignoreShortSeqT) {
        ignoreShortSeq = ignoreShortSeqT;
    }

    private static String getInputFolder() {
        return inputFolder;
    }

    private static void setInputFolder(String inputFolderT) {
        inputFolder = inputFolderT;
    }

    //** Checks User Input >>
    private static String checkUserInput(String inFile, String sSRrep, String minSSR, String maxSSR, String permutations, boolean rand,String ignoreSeqs,String topPercentage){
        setStart(false);

        if (inFile.isEmpty()){
            return "No input Folder";
        }
        File folder = new File(inFile);
        File[] listOfFiles = folder.listFiles();
        if(!folder.exists()){
            return "Folder does not exist";
        }
        if(listOfFiles == null || listOfFiles.length == 0){
            return "No files in Folder";
        }

        int fas = 0;
        //** checks file type
        for(File i: listOfFiles){
            if (i.getName().substring(i.getName().length()-4).equals(".fas") || i.getName().substring(i.getName().length()-4).equals("asta")){
                fas = 1;
                break;
            }
        }
        if (fas == 0){
            return "No valid file";
        }
        if (minSSR.isEmpty()){
            return "No Min value entered";
        }
        if (maxSSR.isEmpty()){
            return "No Max value entered";
        }
        if (sSRrep.isEmpty()){
            return "No SSR repeat value entered";
        }
        if (permutations.isEmpty()){
            return "No permutation value entered";
        }
        if (ignoreSeqs.isEmpty()){
            return "No ignore seq value entered";
        }

        try {
            //** ensures user has entered numbers
            int perm = Integer.parseInt(permutations);
            int minSSR1 = Integer.parseInt(minSSR);
            int maxSSR1 = Integer.parseInt(maxSSR);
            int ssrRep = Integer.parseInt(sSRrep);
            int ignored = Integer.parseInt(ignoreSeqs);
            double topResults = Double.parseDouble(topPercentage);
            if (perm < 50){
                return "Permutation value too low";
            }
            if (minSSR1 < 1){
                return "Min SSR length value too low";
            }
            if (maxSSR1 <= minSSR1){
                return "Max SSR length too low";
            }
            if (ssrRep <2){
                return "SSR repeat value too low";
            }
            if (ignored < 20){
                return "Ignored value too low";
            }
            if (topResults >= 1 || topResults <= 0){
                return "Top result must be between 0 and 1";
            }
            //** Stores valid user settings - if GMATo and Perl aren't installed. install them and restart program.
            setInputFolder(inFile);
            setSecureRandom(rand);
            setMaxMotifLen(maxSSR1);
            setMinMotifLen(minSSR1);
            setMotifRepeats(ssrRep);
            setPermutations(perm);
            setIgnoreShortSeq(ignored);
            setTopResults(topResults);
            return perlAndGMAToInstalled();
        }
        catch (Exception e){
            return e.getMessage();
        }
    }
    //** checks for GMATo and Perl
    private static String perlAndGMAToInstalled(){
        try{
            BufferedReader bufreader;
            Process winproc = Runtime.getRuntime().exec("cmd /c perl -version");
            boolean perl = false;
            if (winproc != null)
            {
                winproc.waitFor();
                bufreader = new BufferedReader(new InputStreamReader(winproc.getInputStream()));
                for (int i = 0; i < 2;i++){
                    String winline = bufreader.readLine();
                    if(winline.startsWith("This is perl")){
                        perl = true;
                        break;
                    }
                }
            }
            if (!perl){
                return("Perl might not be installed don this computer");
            }
            String winfilepath = System.getProperty("user.dir");
            String gmat = winfilepath +"\\gmat.pl";
            String format = winfilepath +"\\formatchunk.pl";
            String gssr = winfilepath +"\\gssr.pl";
            String gsts = winfilepath +"\\gsts.pl";
            String gssrtrim = winfilepath +"\\gssrtrim.pl";
            File gmatPresent = new File(gmat);
            File formatchunkPresent = new File(format);
            File gssrPresent = new File(gssr);
            File gstsPresent = new File(gsts);
            File gssrtrimPresent = new File(gssrtrim);
            if (!gmatPresent.exists()){
                return("gmat.pl is not in the this folder");
            }
            if (!formatchunkPresent.exists()){
                return("Formatchunk.pl is missing from folder");
            }
            if (!gssrPresent.exists()){
                return("gssr.pl is missing");
            }
            if (!gssrtrimPresent.exists()){
                return("gssrtrim.pl is missing");
            }
            if (!gstsPresent.exists()){
                return("gsts.pl is missing");
            }

        }
        catch (Exception e){
            return(e.toString());
        }
        setStart(true);
        return "All settings are valid!";
    }
    //** new thread which creates the Metagenome
}
