package sample;

import javafx.beans.property.DoubleProperty;
import javafx.beans.property.SimpleDoubleProperty;

class ProgressNumber {
    private DoubleProperty progressNum;

//    public final double getProgressNum(){
//        if(progressNum != null)
//            return progressNum.get();
//        return 0;
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
