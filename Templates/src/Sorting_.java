import java.util.ArrayList;
import java.util.Collections;

public class Sorting_ {

    //For sorting arrays
    public static void sort(int[] arr)
    {
        //because Arrays.sort() uses quicksort which is dumb -> O(n^2) worstcase when array already sorted.
        //Collections.sort() uses merge sort     -> O(nlogn) always
        ArrayList<Integer> ls = new ArrayList<Integer>();
        for(int x: arr)
            ls.add(x);
        Collections.sort(ls);
        for(int i=0; i < arr.length; i++)
            arr[i] = ls.get(i);
    }
    public static void sortReverse(int[] arr)
    {
        //because Arrays.sort() uses quicksort which is dumb -> O(n^2) worstcase when array already sorted.
        //Collections.sort() uses merge sort     -> O(nlogn) always
        ArrayList<Integer> ls = new ArrayList<Integer>();
        for(int x: arr)
            ls.add(x);
        ls.sort((a, b) -> b - a);
        for(int i=0; i < arr.length; i++)
            arr[i] = ls.get(i);
    }

}
