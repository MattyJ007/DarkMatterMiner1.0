package sample;

import javafx.beans.property.DoubleProperty;
import javafx.beans.property.SimpleDoubleProperty;
//import javafx.beans.property.SimpleStringProperty;
//import javafx.beans.property.StringProperty;

class ProgressCheck {
    private DoubleProperty progressNum;
//    private StringProperty progressString;

//    public final double getProgressNum(){
//        if(progressNum != null)
//            return progressNum.get();
//        return 0;
//    }
//    public String getProgressString() {
//        if(progressString != null){
//            return progressString.get();
//        }
//        return "";
//    }
//
//    final public StringProperty progressStringProperty() {
//        if (progressString == null){
//            progressString = new SimpleStringProperty("");
//        }
//        return progressString;
//    }
//
//    public void setProgressString(String progressString) {
//        this.progressString.set(progressString);
//    }

    final DoubleProperty progressNumProperty(){
        if(progressNum == null){
            progressNum = new SimpleDoubleProperty(0);
        }
        return progressNum;
    }
    void setProgressNum(double progress) {
        this.progressNumProperty().set(progress);
    }
}
