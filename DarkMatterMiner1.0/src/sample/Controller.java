package sample;

import javafx.scene.control.TextArea;

public class Controller {
    private static TextArea guiOutputText = new TextArea();

    public static TextArea getGuiOutputText(){
        return guiOutputText;
    }
}
