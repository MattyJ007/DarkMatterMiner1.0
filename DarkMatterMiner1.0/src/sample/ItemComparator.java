package sample;

import java.util.Comparator;

class ItemComparator implements Comparator<Sequence> {
        enum Field {
            DINUC, ORFLEN, TRINUC, TOTRANK, BEST_ORF_TRI_FREQ
        }
        private Field field;
        ItemComparator(Field field) {
        this.field = field;
    }
    @Override
    public int compare(Sequence seq1, Sequence seq2) {
        double comparison = 0;
        switch (field) {
            case DINUC:
                comparison = seq1.getDinucleotidePValue() - seq2.getDinucleotidePValue();
                break;
            case ORFLEN:
                comparison = seq1.getOrfLengthPValue() - seq2.getOrfLengthPValue();
                break;
            case TRINUC:
                comparison = seq1.getTrinucelotidePValue() - seq2.getTrinucelotidePValue();
                break;
            case TOTRANK:
                comparison = (double) (seq1.getRankTot() - seq2.getRankTot());
                break;
            case BEST_ORF_TRI_FREQ:
                comparison = (double) (seq1.getRankTot() - seq2.getRankTot());
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
//        return (int) (comparison*100000000);
    }
}
