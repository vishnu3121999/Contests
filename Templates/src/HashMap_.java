import java.util.HashMap;
import java.util.Map;
import java.util.LinkedHashMap;
import java.util.stream.Collectors;

public class HashMap_ {
    //custom multiset (replace with HashMap if needed)
    public static void push(Map<Integer, Integer> map, int k, int v)
    {
        //map[k] += v;
        if(!map.containsKey(k))
            map.put(k, v);
        else
            map.put(k, map.get(k)+v);
    }
    public static void pull(Map<Integer, Integer> map, int k, int v)
    {
        //assumes map[k] >= v
        //map[k] -= v
        int lol = map.get(k);
        if(lol == v)
            map.remove(k);
        else
            map.put(k, lol-v);
    }
    public static HashMap<Character, Integer> sortByValue(HashMap<Character, Integer> hm)
    {
        HashMap<Character, Integer> temp
                = hm.entrySet()
                .stream()
                .sorted((i1, i2)
                        -> i1.getValue().compareTo(
                        i2.getValue()))
                .collect(Collectors.toMap(
                        Map.Entry::getKey,
                        Map.Entry::getValue,
                        (e1, e2) -> e1,
                        LinkedHashMap::new));

        return temp;
    }

}
