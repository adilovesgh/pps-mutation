package mutation.sim;

import java.util.*;

public class Listagen {
    public Listagen(Mutagen mutagen) {
        List<String> patterns = mutagen.getPatterns();
        List<String> actions = mutagen.getActions();
        data = new HashSet<Long>();
        for (int i = 0; i < patterns.size(); ++ i)
            this.add(patterns.get(i), actions.get(i));
    }

    private void recursive(String[] l, int depth, String action, char[] pattern, int i) {
        if (i >= l.length && i >= depth) {
            // parse result, insert pattern->result
            char[] result = new char[depth];
            for (int j = 0; j < depth; ++ j) {
                if (j >= action.length()) result[j] = pattern[j];
                else {
                    if (action.charAt(j) >= '0' && action.charAt(j) <= '9')
                        result[j] = pattern[action.charAt(j) - '0'];
                    else result[j] = action.charAt(j);
                }
            }
//            Log.record(String.valueOf(pattern) + " => " + String.valueOf(result) + " : " + ((transform(pattern) << 24l) + transform(result)));
            data.add((transform(pattern) << 24l) + transform(result));
            return;
        }
        if (i < l.length) {
            for (int j = 0; j < l[i].length(); ++ j) {
                pattern[i] = l[i].charAt(j);
                recursive(l, depth, action, pattern, i + 1);
            }
        } else {
            for (int j = 0; j < alphabet.length; ++ j) {
                pattern[i] = alphabet[j];
                recursive(l, depth, action, pattern, i + 1);
            }
        }
    }

    private void add(String pattern, String action) {
        String[] l = pattern.split(";");
        int depth = Math.max(action.length(), l.length);
        for (int i = 0; i < action.length(); ++ i)
            if (action.charAt(i) >= '0' && action.charAt(i) <= '9')
                depth = Math.max(depth, 1 + action.charAt(i) - '0');
        char[] pat = new char[depth];
        recursive(l, depth, action, pat, 0);
    }

    private long transform(char[] a) {
        long ret = 0;
        for (int i = 0; i < a.length; ++ i)
            ret = ret * 4l + translate(a[i]);
        return ret + (a.length << 20);
    }

    private int translate(char c) {
        if (c == 'a') return 0;
        if (c == 'c') return 1;
        if (c == 'g') return 2;
        if (c == 't') return 3;
        return -1;
    }

    public boolean equals(Listagen other) {
        return data.equals(other.data);
    }

    Set<Long> data;
    static char[] alphabet = "abcd".toCharArray();
}
