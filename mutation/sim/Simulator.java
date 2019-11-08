package mutation.sim;

import javax.tools.JavaCompiler;
import javax.tools.StandardJavaFileManager;
import javax.tools.ToolProvider;
import java.awt.*;
import java.io.*;
import java.lang.reflect.InvocationTargetException;
import java.net.URI;
import java.net.URISyntaxException;
import java.text.DecimalFormat;
import java.util.*;
import java.util.List;
import java.util.concurrent.TimeoutException;

public class Simulator {
    private static final String root = "mutation";
    private static Random random;
    private static long timeLimit = 60000;
    private static boolean gui = false;
    private static double fps = 5;
    private static long seed = 1;
    private static String cfgPath = "mutation/g5/mutagen_single_rule.cfg";
    private static String name;

    private static int m = 15;
    private static int trials = 100000;
    private static Mutagen target, result;
    private static Player player;

    public static void main(String[] args) throws Exception {
        parseArgs(args);

        Log.record("Starting game with " + trials + " trials");

        HTTPServer server = null;
        if (gui) {
            timeLimit *= 1000;
            server = new HTTPServer();
            Log.record("Hosting HTTP Server on " + server.addr());
            if (!Desktop.isDesktopSupported())
                Log.record("Desktop operations not supported");
            else if (!Desktop.getDesktop().isSupported(Desktop.Action.BROWSE))
                Log.record("Desktop browse operation not supported");
            else {
                try {
                    Desktop.getDesktop().browse(new URI("http://localhost:" + server.port()));
                } catch (URISyntaxException e) {
                    StringWriter errors = new StringWriter();
                    e.printStackTrace(new PrintWriter(errors));
                    System.err.println(errors.toString());
                }
            }
        }

        // Loading mutagen
        target = loadMutagenConfig(cfgPath);
        Console console = new Console(trials, m, target, gui, server, name, 1000.0 / fps);

        Timer thread = new Timer();
        thread.start();

        Log.record("Player " + name + " starts!");
        thread.call_start(() -> {
            console.setEndTime(System.currentTimeMillis() + timeLimit);
            return player.Play(console, m);
        });
        try {
            result = thread.call_wait(timeLimit);
        } catch (TimeoutException e) {
            System.err.println("Player timed out.");
            result = console.getLastGuess();
            // console.reportScore("", "", "Score pending");
            // char[] pool = {'a', 'c', 'g', 't'};
            // char[] data = new char[1000000];
            // random = new Random();
            // for (int i = 0; i < 1000000; ++ i)
            //     data[i] = pool[Math.abs(random.nextInt() % 4)];
            // String testStr = String.valueOf(data);
            // Set<Long> s1 = target.jaccardSet(testStr);
            // Set<Long> s2 = result.jaccardSet(testStr);
            // int intersection = 0;
            // for (Long entry : s2)
            //     if (s1.contains(entry))
            //         ++ intersection;
            // int union = s1.size() + s2.size() - intersection;
            //
            // double score = (double)intersection / (double) union;
            // System.out.println("Jaccard score is: " + intersection + "/" + union + " (" + String.format("%.4f", score * 100.0) + "%).");
            // console.reportScore(String.join("@", result.getPatterns()), String.join("@", result.getActions()), "Jaccard score: " + intersection + "/" + union + " (" + String.format("%.4f", score * 100.0) + "%)");
            // System.exit(-1);
        } catch (Exception e) {
            StringWriter errors = new StringWriter();
            e.printStackTrace(new PrintWriter(errors));
            System.err.println(errors.toString());
            System.exit(-1);
        }
        long elapsedTime = thread.getElapsedTime();
        Log.record("Player finished in " + elapsedTime + "ms.");

        System.out.println("Player " + name + " made " + console.getNumGuesses() + " guesses and " + console.getNumExps() + " experiments");
        if (console.isCorrect() || console.testEquiv(result)) {
            System.out.println("Correct!");
        } else {
            console.reportScore("", "", "Score pending");
            System.out.println("Failed, calculating Jaccard score.");
            List<String> patterns = result.getPatterns();
            List<String> actions = result.getActions();
            System.out.println("Player's guess: ");
            for (int i = 0; i < patterns.size(); ++ i)
                System.out.println(patterns.get(i) + " => " + actions.get(i));

            char[] pool = {'a', 'c', 'g', 't'};
            char[] data = new char[1000000];
            random = new Random();
            for (int i = 0; i < 1000000; ++ i)
                data[i] = pool[Math.abs(random.nextInt() % 4)];
            String testStr = String.valueOf(data);
            Set<Long> s1 = target.jaccardSet(testStr);
            Log.record("Target mutations: " + s1.size());
            Set<Long> s2 = result.jaccardSet(testStr);
            Log.record("Result mutations: " + s2.size());
            int intersection = 0;
            for (Long entry : s2)
                if (s1.contains(entry))
                    ++ intersection;
            int union = s1.size() + s2.size() - intersection;

            double score = (double)intersection / (double) union;
            System.out.println("Jaccard score is: " + intersection + "/" + union + " (" + String.format("%.4f", score * 100.0) + "%).");
            console.reportScore(String.join("@", result.getPatterns()), String.join("@", result.getActions()), "Jaccard score: " + intersection + "/" + union + " (" + String.format("%.4f", score * 100.0) + "%)");
        }

        System.exit(0);
    }

    private static Mutagen loadMutagenConfig(String path) throws IOException {
        Mutagen result = new Mutagen();
        BufferedReader br = new BufferedReader(new FileReader(path));
        String line;
        while ((line = br.readLine()) != null) {
            if (line.length() == 0 || line.charAt(0) == '#') continue;
            String[] strs = line.split("@");
            if (strs.length != 2)
                throw new IOException("Bad config file syntax");
            // if (strs[0].split(";").length != strs[1].length())
            //     throw new IOException("Bad config file syntax");
            result.add(strs[0], strs[1]);
        }
        return result;
    }

    private static void parseArgs(String[] args) throws IOException, ClassNotFoundException, InstantiationException, InvocationTargetException, IllegalAccessException, NoSuchMethodException {
        int i = 0;
        for (; i < args.length; ++i) {
            switch (args[i].charAt(0)) {
                case '-':
                    if (args[i].startsWith("-p") || args[i].equals("--player")) {
                        // TODO
                        if (++i == args.length) {
                            throw new IllegalArgumentException("Missing player name");
                        }
                        name = args[i];
                        player = loadPlayer(name);
                    } else if (args[i].equals("-m")) {
                        if (++i == args.length) {
                            throw new IllegalArgumentException("Missing number of mutations per experiments");
                        }
                        m = Integer.parseInt(args[i]);
                    } else if (args[i].equals("-c") || args[i].equals("--cfgpath")) {
                        if (++i == args.length) {
                            throw new IllegalArgumentException("Missing config file.");
                        }
                        cfgPath = args[i];
                    } else if (args[i].equals("-t") || args[i].equals("--trials")) {
                        if (++i == args.length) {
                            throw new IllegalArgumentException("Missing maximum number of experiments");
                        }
                        trials = Integer.parseInt(args[i]);
                    } else if (args[i].equals("-s") || args[i].equals("--seed")) {
                        if (++i == args.length) {
                            throw new IllegalArgumentException("Missing seed");
                        }
                        seed = Long.parseLong(args[i]);
                    } else if (args[i].equals("-tl") || args[i].equals("--timelimit")) {
                        if (++i == args.length) {
                            throw new IllegalArgumentException("Missing time limit");
                        }
                        timeLimit = Integer.parseInt(args[i]);
                    } else if (args[i].equals("-l") || args[i].equals("--logfile")) {
                        if (++i == args.length) {
                            throw new IllegalArgumentException("Missing logfile name");
                        }
                        Log.setLogFile(args[i]);
                    } else if (args[i].equals("-g") || args[i].equals("--gui")) {
                        gui = true;
                    } else if (args[i].equals("--fps")) {
                        if (++i == args.length) {
                            throw new IllegalArgumentException("Missing fps");
                        }
                        fps = Double.parseDouble(args[i]);
                    } else if (args[i].equals("-v") || args[i].equals("--verbose")) {
                        Log.activate();
                    } else {
                        throw new IllegalArgumentException("Unknown argument \"" + args[i] + "\"");
                    }
                    break;
                default:
                    throw new IllegalArgumentException("Unknown argument \"" + args[i] + "\"");
            }
        }

        Log.record("Time limit: " + timeLimit + "ms");
        Log.record("GUI " + (gui ? "enabled" : "disabled"));
        if (gui)
            Log.record("FPS: " + fps);
    }

    private static Set<File> directory(String path, String extension) {
        Set<File> files = new HashSet<File>();
        Set<File> prev_dirs = new HashSet<File>();
        prev_dirs.add(new File(path));
        do {
            Set<File> next_dirs = new HashSet<File>();
            for (File dir : prev_dirs)
                for (File file : dir.listFiles())
                    if (!file.canRead()) ;
                    else if (file.isDirectory())
                        next_dirs.add(file);
                    else if (file.getPath().endsWith(extension))
                        files.add(file);
            prev_dirs = next_dirs;
        } while (!prev_dirs.isEmpty());
        return files;
    }

    public static Player loadPlayer(String name) throws IOException, ClassNotFoundException, InstantiationException, InvocationTargetException, IllegalAccessException, NoSuchMethodException {
        String sep = File.separator;
        Set<File> player_files = directory(root + sep + name, ".java");
        File class_file = new File(root + sep + name + sep + "Player.class");
        long class_modified = class_file.exists() ? class_file.lastModified() : -1;
        if (class_modified < 0 || class_modified < last_modified(player_files) ||
                class_modified < last_modified(directory(root + sep + "sim", ".java"))) {
            JavaCompiler compiler = ToolProvider.getSystemJavaCompiler();
            if (compiler == null)
                throw new IOException("Cannot find Java compiler");
            StandardJavaFileManager manager = compiler.
                    getStandardFileManager(null, null, null);
//            long files = player_files.size();
            Log.record("Compiling for player " + name);
            if (!compiler.getTask(null, manager, null, Arrays.asList("-g"), null,
                    manager.getJavaFileObjectsFromFiles(player_files)).call())
                throw new IOException("Compilation failed");
            class_file = new File(root + sep + name + sep + "Player.class");
            if (!class_file.exists())
                throw new FileNotFoundException("Missing class file");
        }
        ClassLoader loader = Simulator.class.getClassLoader();
        if (loader == null)
            throw new IOException("Cannot find Java class loader");
        @SuppressWarnings("rawtypes")
        Class<?> raw_class = loader.loadClass(root + "." + name + ".Player");
        return (Player) raw_class.newInstance();
    }

    private static long last_modified(Iterable<File> files) {
        long last_date = 0;
        for (File file : files) {
            long date = file.lastModified();
            if (last_date < date)
                last_date = date;
        }
        return last_date;
    }
}
