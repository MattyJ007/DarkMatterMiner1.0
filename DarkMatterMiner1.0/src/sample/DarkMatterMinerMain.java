package sample;

import javafx.application.Application;
import javafx.fxml.FXMLLoader;
import javafx.scene.Parent;
import javafx.scene.Scene;
import javafx.stage.Stage;

public class DarkMatterMinerMain extends Application {
    //**Launches DarkMatterMiner GUI
    public void start(Stage primaryStage) throws Exception{
        Parent root = FXMLLoader.load(getClass().getResource("DarkMatterMinerUI.fxml"));
        Scene scene = new Scene(root);
        primaryStage.setTitle("Dark Matter Miner 1.0");
        primaryStage.setScene(scene);
        primaryStage.show();
        primaryStage.setOnCloseRequest(t -> System.exit(0));
    }

    public static void main(String[] args) {
        launch(args);
    }
}
