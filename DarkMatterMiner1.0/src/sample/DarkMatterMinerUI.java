package sample;

import javafx.application.Application;
import javafx.geometry.Insets;
import javafx.geometry.Pos;
import javafx.scene.Scene;
import javafx.scene.control.*;
import javafx.scene.layout.GridPane;
import javafx.scene.layout.HBox;
import javafx.scene.text.Font;
import javafx.scene.text.FontWeight;
import javafx.scene.text.Text;
import javafx.stage.Stage;
import java.io.BufferedReader;
import java.io.File;
import java.io.InputStreamReader;

public class DarkMatterMinerUI extends Application {
    private static boolean start; //** variable that prevents launching Metagenome if User settings aren't correct
    private static boolean secureRandom; //** Allows user to use either Random() or SecureRandom()
    private static double topResults = 0.005; //** Allows user to choose number of seqs for end fasta file
    private static int motifRepeats = 4;
    private static int permutations = 1000;
    private static int minMotifLen = 2;
    private static int maxMotifLen = 10;
    private static int ignoreShortSeq = 20;
    private static String inputFolder;

    //** Getters and Setters of above variables^^^

    private static void setStart(boolean start) {
        DarkMatterMinerUI.start = start;
    }

    private static boolean isStart() {
        return start;
    }

    private static boolean isSecureRandom() {
        return secureRandom;
    }

    private static void setSecureRandom(boolean secureRandom) {
        DarkMatterMinerUI.secureRandom = secureRandom;
    }

    private static double getTopResults() {
        return topResults;
    }

    private static void setTopResults(double topResults) {
        DarkMatterMinerUI.topResults = topResults;
    }

    static int getMotifRepeats() {
        return motifRepeats;
    }

    private static void setMotifRepeats(int motifRepeats) {
        DarkMatterMinerUI.motifRepeats = motifRepeats;
    }

    static int getPermutations() {
        return permutations;
    }

    private static void setPermutations(int permutations) {
        DarkMatterMinerUI.permutations = permutations;
    }

    static int getMinMotifLen() {
        return minMotifLen;
    }

    private static void setMinMotifLen(int minMotifLen) {
        DarkMatterMinerUI.minMotifLen = minMotifLen;
    }

    static int getMaxMotifLen() {
        return maxMotifLen;
    }

    private static void setMaxMotifLen(int maxMotifLen) {
        DarkMatterMinerUI.maxMotifLen = maxMotifLen;
    }

    static int getIgnoreShortSeq() {
        return ignoreShortSeq;
    }

    private static void setIgnoreShortSeq(int ignoreShortSeq) {
        DarkMatterMinerUI.ignoreShortSeq = ignoreShortSeq;
    }

    private static String getInputFolder() {
        return inputFolder;
    }

    private static void setInputFolder(String inputFolder) {
        DarkMatterMinerUI.inputFolder = inputFolder;
    }

    //**Launches DarkMatterMiner GUI

    public void start(Stage primaryStage) throws Exception{
        GridPane grid = new GridPane();
        grid.setAlignment(Pos.CENTER);
        grid.setHgap(10);
        grid.setVgap(10);
        grid.setPadding(new Insets(25, 25, 25, 25));

        Scene scene = new Scene(grid, 650, 500);

        Label folder = new Label("Input/Output folder: ");
        folder.setPadding(new Insets(0, 20, 0, 0));
        TextField folderInputText = new TextField();
        folderInputText.setMinWidth(350);

        //** Horizontal box used to keep the grid neat
        HBox folderIn = new HBox();
        folderIn.setPadding(new Insets(0, 5, 0, 0));
        folderIn.setAlignment(Pos.CENTER_LEFT);
        folderIn.getChildren().add(folder);
        folderIn.getChildren().add(folderInputText);

        Text gMAToParameters = new Text("GMATo parameters");
        gMAToParameters.setFont(Font.font("Tahoma", FontWeight.NORMAL, 15));

        Label minMotifRepeats = new Label("Min motif repeats: ");
        TextField minMotifRepeatstextBox = new TextField();
        minMotifRepeatstextBox.setText("6");
        minMotifRepeatstextBox.setMaxWidth(50);
        minMotifRepeats.setPadding(new Insets(0, 10, 0, 0));

        Label minMotifLength = new Label("Min motif length: ");
        minMotifLength.setPadding(new Insets(0, 10, 0, 30));
        TextField minMotifLengthTextBox = new TextField();
        minMotifLengthTextBox.setText("2");
        minMotifLengthTextBox.setMaxWidth(50);

        Label maxMotifLength = new Label("Max motif length: ");
        maxMotifLength.setPadding(new Insets(0, 10, 0, 30));
        TextField maxMotifLengthTextBox = new TextField();
        maxMotifLengthTextBox.setText("10");
        maxMotifLengthTextBox.setMaxWidth(50);
        //**HBox Saves space
        HBox gmatoParameters = new HBox();
        gmatoParameters.setPadding(new Insets(0, 5, 0, 0));
        gmatoParameters.setAlignment(Pos.CENTER_LEFT);
        gmatoParameters.getChildren().add(minMotifRepeats);
        gmatoParameters.getChildren().add(minMotifRepeatstextBox);
        gmatoParameters.getChildren().add(minMotifLength);
        gmatoParameters.getChildren().add(minMotifLengthTextBox);
        gmatoParameters.getChildren().add(maxMotifLength);
        gmatoParameters.getChildren().add(maxMotifLengthTextBox);

        Text darkMatterParameters = new Text("Dark Matter Parameters");
        darkMatterParameters.setFont(Font.font("Tahoma", FontWeight.NORMAL, 15));

        Label permutations = new Label("Permutations: ");
        TextField permutationstextBox = new TextField();
        permutationstextBox.setMaxWidth(50);
        permutationstextBox.setText("1000");

        CheckBox random = new CheckBox("SecureRandom()");
        Label shortSequenceLength = new Label("Ignore sequences shorter than (bp) : ");
        TextField shortSequenceLengthTextBox = new TextField();
        shortSequenceLengthTextBox.setMaxWidth(50);
        shortSequenceLengthTextBox.setText("30");

        Label topFractionOfResults = new Label("Fraction of Potential Dark Matter to Keep : ");
        TextField topFractionTextBox = new TextField();
        topFractionTextBox.setMaxWidth(50);
        topFractionTextBox.setText("0.005");

        grid.add(folderIn, 0, 0, 2, 1);
        grid.add(gMAToParameters, 0, 1, 2, 1);
        grid.add(gmatoParameters, 0, 2, 6, 1);
        grid.add(darkMatterParameters, 0, 3, 6, 1);
        grid.add(permutations, 0, 4);
        grid.add(permutationstextBox, 1, 4);
        grid.add(random, 0, 5);
        grid.add(shortSequenceLength, 0, 6);
        grid.add(shortSequenceLengthTextBox, 1, 6);
        grid.add(topFractionOfResults, 0, 7);
        grid.add(topFractionTextBox, 1, 7);

        Button runbtn = new Button("Run");
        grid.add(runbtn, 0, 8);
        ScrollPane scroll = new ScrollPane();
        grid.add(scroll, 0, 9, 2, 1);
        scroll.setContent(Controller.getGuiOutputText());

        runbtn.setOnAction(event ->
                {
                    //** If the UserInput check has already been run and successfully
                    //**passed the run button will change function and become an exit button
                    //**this prevents users trying to submit a new folder whilst the program is still running.
                    if (isStart()) {
                        System.exit(0);
                    }
                    //** sets the text in the GUI's output window to indicate validity of user input
                    Controller.getGuiOutputText().setText(checkUserInput(folderInputText.getText(), minMotifRepeatstextBox.getText(),
                            minMotifLengthTextBox.getText(), maxMotifLengthTextBox.getText(),
                            permutationstextBox.getText(), random.selectedProperty().getValue(), shortSequenceLengthTextBox.getText(), topFractionTextBox.getText()));
                    //** if user input is valid the run button changes to an exit button
                    //** a new thread is launched that handles the background processing thus keeping the GUI active
                    if (isStart()) {
                        runbtn.setText("Exit");
                        Thread t1 = new Thread(new RunClass());
                        t1.start();
                    }
                }
        );
        //** Creates the GUI window
        primaryStage.setTitle("DarkMatterMiner1.0");
        primaryStage.setScene(scene);
        primaryStage.show();
    }

    //** Checks User Input >>

    private static String checkUserInput(String inFile, String sSRrep, String minSSR, String maxSSR,
                               String permutations, boolean rand,String ignoreSeqs,String topPercentage){
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
            if (perm <= 50){
                return "Permutation value too low";
            }
            if (minSSR1 <= 1){
                return "Min SSR length value too low";
            }
            if (maxSSR1 <= minSSR1){
                return "Max SSR length too low";
            }
            if (ssrRep <=2){
                return "SSR repeat value too low";
            }
            if (ignored <= 20){
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
    private class RunClass implements Runnable{
        @Override
        public void run() {
            try{
                System.out.println("\nAnalysis starting...");
                Metagenome.create(getInputFolder(),isSecureRandom());
                System.out.println("\n-------------------------------------\n             Complete\n-------------------------------------");
            }
            catch (Exception e){
                System.out.println(e.getMessage());
            }
        }
    }

    public static void main(String[] args) {
        launch(args);
    }
}
