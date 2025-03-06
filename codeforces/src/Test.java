import java.util.*;

import static java.lang.Math.min;

public class Test {


    public static void main(String[] args) {
        ArrayList<TC> list = new ArrayList<>();
//        list.add(new TC(10,new ArrayProps(5,1,10,true)));
//        list.add(new TC(10,new ArrayProps(5,1,10,false)));
//        list.add(new TC(10,new ArrayProps(10,1,100,true)));
//        list.add(new TC(10,new ArrayProps(1,1,15,true)));
//        list.add(new TC(10,new ArrayProps(2,1,15,true)));
//        list.add(new TC(10,new ArrayProps(6,1,10,true)));

        list.add(new TC(100,new StringProps(10,true,"a-z"),new StringProps(10,true,"a-z")));
        list.add(new TC(100,new StringProps(9,true,"01")));

        for(TC tc:list){
            int noOfTestCases = tc.noOfTestCases;
            for (int i = 0; i < noOfTestCases; i++) {
                // get props
                Props props1 = tc.propsList.get(0);
                Props props2 = tc.propsList.get(1);

                // create from props
                String u = Generators.generateString(props1);
                String v = Generators.generateString(props2);

                
//
//                String bruteForceAnswer = Main.bruteForce(u,v);
//                String optimalAnswer = Main.optimal(u,v);

//                if(bruteForceAnswer!=null && bruteForceAnswer.length()==u.length()+v.length())bruteForceAnswer=null;
//                if(bruteForceAnswer==null || optimalAnswer==null || bruteForceAnswer.length()!=optimalAnswer.length()){
//                    System.out.println("Failed on: "+"u="+u+", v="+v);
//                    System.out.println("Expected: "+bruteForceAnswer);
//                    System.out.println("Got: "+optimalAnswer);
//                }
            }
        }
    }



    private static boolean compareArrays(int[] arr1, int[] arr2) {
        for (int i = 0; i <arr1.length ; i++) {
            if(arr1[i]!=arr2[i])return false;
        }
        return true;
    }
}

class TC{
    int noOfTestCases;
    List<Props> propsList;
    public TC(int noOfTestCases, Props... propsList) {
        this.noOfTestCases = noOfTestCases;
        this.propsList = Arrays.asList(propsList);  // Converts varargs to a List
    }
}


interface Props{}
class ArrayProps implements Props{
    int arraySize;
    int lowerBound;
    int upperBound;
    boolean isDuplicatesAllowed;
    ArrayProps(int arraySize, int lowerBound, int upperBound, boolean isDuplicatesAllowed){
        this.arraySize = arraySize;
        this.lowerBound = lowerBound;
        this.upperBound = upperBound;
        this.isDuplicatesAllowed = isDuplicatesAllowed;
    }
}
class StringProps implements Props{
    int length;
    boolean duplicatesAllowed;
    String charsAllowed;
    StringProps(int length, boolean duplicatesAllowed, String charsAllowed){
        this.length=length;
        this.duplicatesAllowed=duplicatesAllowed;
        this.charsAllowed=charsAllowed;
    }
}




class Generators{

    static Random random = new Random();

    static int[] generateArray(Props props){
        ArrayProps arrayProps = (ArrayProps)props;
        int size=arrayProps.arraySize;
        int lowerLimit = arrayProps.lowerBound;
        int upperLimit = arrayProps.upperBound;
        boolean duplicatesAllowed = arrayProps.isDuplicatesAllowed;

        int[] arr = new int[size];
        if (duplicatesAllowed) {
            for (int i = 0; i < size; i++) {
                arr[i] = random.nextInt(upperLimit - lowerLimit + 1) + lowerLimit;
            }
        } else {
            Set<Integer> usedNumbers = new HashSet<>();
            int i = 0;
            while (i < size) {
                int num = random.nextInt(upperLimit - lowerLimit + 1) + lowerLimit;
                if (!usedNumbers.contains(num)) {
                    arr[i] = num;
                    usedNumbers.add(num);
                    i++;
                }
            }
        }
        return arr;
    }

    private static String expandAllowed(String charsAllowed) {
        if (charsAllowed.length() == 3 && charsAllowed.charAt(1) == '-') {
            char start = charsAllowed.charAt(0);
            char end = charsAllowed.charAt(2);
            if (start <= end) {
                StringBuilder sb = new StringBuilder();
                for (char c = start; c <= end; c++) {
                    sb.append(c);
                }
                return sb.toString();
            }
        }
        return charsAllowed;  // Return as-is if format doesn't match "X-Y"
    }

    // charsAllowed - A string containing all allowed characters, or a shorthand range like "a-z".
    public static String generateString(Props props) {
        StringProps stringProps = (StringProps)props;
        int length = stringProps.length;
        boolean duplicatesAllowed = stringProps.duplicatesAllowed;
        String charsAllowed = stringProps.charsAllowed;

        // Expand the shorthand if needed.
        String allowed = expandAllowed(charsAllowed);

        if (allowed.isEmpty()) {
            throw new IllegalArgumentException("Allowed characters string cannot be empty.");
        }

        if (!duplicatesAllowed && length > allowed.length()) {
            throw new IllegalArgumentException("When duplicates are not allowed, length cannot exceed the number of allowed characters.");
        }

        StringBuilder result = new StringBuilder(length);
        if (duplicatesAllowed) {
            // For each position, choose a random character from allowed.
            for (int i = 0; i < length; i++) {
                int index = random.nextInt(allowed.length());
                result.append(allowed.charAt(index));
            }
        } else {
            // Convert the allowed characters to a list, shuffle it, and then pick the first 'length' characters.
            List<Character> charList = new ArrayList<>();
            for (char ch : allowed.toCharArray()) {
                charList.add(ch);
            }
            Collections.shuffle(charList, random);
            for (int i = 0; i < length; i++) {
                result.append(charList.get(i));
            }
        }
        return result.toString();
    }

}





