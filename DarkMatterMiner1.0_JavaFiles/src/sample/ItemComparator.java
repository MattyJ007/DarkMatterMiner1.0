package sample;

import java.util.Comparator;

class ItemComparator implements Comparator<Sequence> {
        enum Field {
            TOTRANK
//            , BEST_ORF_TRI_FREQ,DINUC, ORFLEN, TRINUC,
        }
        private Field field;
        ItemComparator(Field field) {
        this.field = field;
    }
    @Override
    public int compare(Sequence seq1, Sequence seq2) {
        double comparison = 0;
        switch (field) {
//            case DINUC:
//                comparison = seq1.getDinucleotidePValue() - seq2.getDinucleotidePValue();
//                break;
//            case ORFLEN:
//                comparison = seq1.getOrfLengthPValue() - seq2.getOrfLengthPValue();
//                break;
//            case TRINUC:
//                comparison = seq1.getTrinucelotidePValue() - seq2.getTrinucelotidePValue();
//                break;
//            case BEST_ORF_TRI_FREQ:
//                comparison =  seq1.getRankBestORFframeTri() - seq2.getRankBestORFframeTri();
//                break;
            case TOTRANK:
                comparison =  seq1.getRankTot() - seq2.getRankTot();
                break;
            default:
                System.out.println("\n\nItem Comparator issue!!\n\n-----------------------\n\n");
                break;
        }
        if (comparison > 0){
            return 1;
        }
        else if(comparison < 0){
            return -1;
        }
        else {
            return 0;
        }
    }
}
