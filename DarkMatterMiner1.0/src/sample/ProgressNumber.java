package sample;

import javafx.beans.property.DoubleProperty;
import javafx.beans.property.SimpleDoubleProperty;

public class ProgressNumber {
    private DoubleProperty progressNum;

    public final double getProgressNum(){
        if(progressNum != null)
            return progressNum.get();
        return 0;
    }
    public final DoubleProperty progressNumProperty(){
        if(progressNum == null){
            progressNum = new SimpleDoubleProperty(0);
        }
        return progressNum;
    }
    public void setProgressNum (double progress) {
        this.progressNumProperty().set(progress);
    }
}
