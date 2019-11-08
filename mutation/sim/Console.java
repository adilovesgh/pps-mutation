package mutation.sim;

import java.lang.*;
import java.io.*;
import java.net.*;
import java.awt.Desktop;

// This class wraps the Mutagen with limited usage
public class Console {

    public Console(int limit, int m, Mutagen mutagen, boolean gui, HTTPServer server, String playerName, double refresh) {
        this.mutagen = mutagen;
        this.listagen = new Listagen(mutagen);
        this.limit = limit;
        this.m = m;
        this.gui = gui;
        this.server = server;
        this.playerName = playerName;
        this.refresh = refresh;
        this.lastGuess = new Mutagen();

        if (gui) {
            sendGUI(server, state("", "", String.join("@", mutagen.getPatterns()), String.join("@", mutagen.getActions()), "", "", ""));
        }
    }

    public String Mutate(String genome) {
        if (correct || numExps >= limit || genome.length() != 1000)
            return "";
        ++numExps;
        String mutated = mutagen.Mutate(genome, m);
        if (gui) {
            sendGUI(server, state(genome, mutated, "", "", "", "", ""));
        }
        return mutated;
    }

    public boolean Guess(Mutagen other) {
        if (correct) return correct;
        ++ numGuesses;
        correct = mutagen.equals(other);
        lastGuess.getPatterns().clear();
        lastGuess.getActions().clear();
        lastGuess.getPatterns().addAll(other.getPatterns());
        lastGuess.getActions().addAll(other.getActions());
        if (gui) {
            String score = "";
            if (correct) score = "Correct!";
            sendGUI(server, state("", "", "", "", String.join("@", other.getPatterns()), String.join("@", other.getActions()), score));
        }
        return correct;
    }

    public boolean testEquiv(Mutagen other) {
        if (correct) return correct;
        ++ numGuesses;
        correct = listagen.equals(new Listagen(other));
        lastGuess.getPatterns().clear();
        lastGuess.getActions().clear();
        lastGuess.getPatterns().addAll(other.getPatterns());
        lastGuess.getActions().addAll(other.getActions());
        if (gui) {
            String score = "";
            if (correct) score = "Correct!";
            sendGUI(server, state("", "", "", "", String.join("@", other.getPatterns()), String.join("@", other.getActions()), score));
        }
        return correct;
    }

    public int getNumExpsLeft() { return limit - numExps; }

    public long getTimeLeft() {
        return endTime - System.currentTimeMillis();
    }

    public int getNumGuesses() {
        return numGuesses;
    }

    public int getNumExps() {
        return numExps;
    }

    public boolean isCorrect() {
        return correct;
    }

    public int getNumberOfMutations() {
        return mutagen.getNumberOfMutations();
    }

    public Mutagen getLastGuess() {
        return lastGuess;
    }

    public void reportScore(String patterns, String actions, String score) {
        if (gui)
            sendGUI(server, state("", "", "", "", patterns, actions, score));
    }

    public void setEndTime(long time) {
        endTime = time;
    }

    private String state(String genome, String mutated, String tPatterns, String tActions, String patterns, String actions, String score) {
        return playerName + "," + refresh + "," + numExps + "," + numGuesses + "," + score + ","
                + genome + "," + mutated + "," + tPatterns + "," + tActions + "," + patterns + "," + actions;
    }

    private void sendGUI(HTTPServer server, String content) {
        if (server == null) return;
        String path = null;
        for (; ; ) {
            for (; ; ) {
                try {
                    path = server.request();
                    break;
                } catch (IOException e) {
                    Log.record("HTTP request error " + e.getMessage());
                }
            }
            if (path.equals("data.txt")) {
                try {
                    server.reply(content);
                } catch (IOException e) {
                    Log.record("HTTP dynamic reply error " + e.getMessage());
                }
                return;
            }
            if (path.equals("")) path = "index.html";
            else if (!Character.isLetter(path.charAt(0))) {
                Log.record("Potentially malicious HTTP request \"" + path + "\"");
                break;
            }

            File file = new File("statics" + File.separator + path);
            if (file == null) {
                Log.record("Unknown HTTP request \"" + path + "\"");
            } else {
                try {
                    server.reply(file);
                } catch (IOException e) {
                    Log.record("HTTP static reply error " + e.getMessage());
                }
            }
        }
    }

    private int numExps = 0, limit, m = 1, numGuesses = 0;
    private boolean correct = false, gui = false;
    private Mutagen mutagen, lastGuess;
    private HTTPServer server;
    private String playerName;
    private double refresh;
    private Listagen listagen;
    private long endTime;
}
