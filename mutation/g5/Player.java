package mutation.g5;

import java.awt.Point;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.*;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.TimeoutException;

import javafx.util.Pair;
import mutation.sim.Console;
import mutation.sim.Log;
import mutation.sim.Mutagen;

public class Player extends mutation.sim.Player {

	private Random random;
	private Console console;
	private int numMutations;
	private int m;
	private List<Map<Character, Integer>> patternTracker = new ArrayList<>();
	private List<Map<Character, Integer>> actionTracker = new ArrayList<>();
	private static final int DISTANCE_THRESHOLD_SB = 12;
	private static final double PERCENTAGE_DISTANCE_THRESHOLD_MB = 0.15; 
	private static final int NONEXISTENT_BASE_THRESHOLD = 50;
	private static final double PERCENT_NONEXISTENT_BASE_THRESHOLD = 0.05;
	private static final double PERCENTAGE_MUTATION_SIZE_DELTA_THRESHOLD = 0.25;
	private List<Integer> numMatchesForMutationSizes = new ArrayList<>();
	private List<Mutagen> possibleMutagens = new ArrayList<>();
	private List<String> genomes = new ArrayList<>();
	private List<String> mutatedGenomes = new ArrayList<>();
	private List<Integer> numMutationsList = new ArrayList<>();
	private List<Mutagen> mutagensGuessed = new ArrayList<>();
	private List<IntervalArchive> intervalArchives = new ArrayList<>();
	private Map<Integer, List<IntervalArchive>> multipleRuleTracker = new TreeMap<>();
	private Map<Integer, Integer> intervalSizeTracker = new TreeMap<>();

	public Player() { 
		random = new Random();
	}

	private String randomString() {
		char[] pool = {'a', 'c', 'g', 't'};
		String result = "";
		for (int i = 0; i < 1000; ++ i)
			result += pool[Math.abs(random.nextInt() % 4)];
		return result;
	}

	@Override
	public Mutagen Play(Console console, int m) {
		while(true) {
			this.m = m;
			this.console = console;
			setUpStructures();
			int numBasesMutatedPerMutation = 0;
			int totalNumMutations = 0;
			int totalNumBasesMutated = 0;
			List<Mutagen> approach2Mutagens = new ArrayList<>();
			int prevIterationError = 0;
			int errorCount = 0;
			for(int i = 0; i < 100; i++) {
				//				String genome = "tgccacggtctgttttgcatgcacctagttcgttttctggggcgttagacagactatgttgcttcaccccgaaccccaaatttatacatttctgcagttcccgtttgaatgtgctaatcggtaacatctgccttgtgacgcagaaggtatcaagctggaatgccacacggacactcgagcgctgatctcccgggttggatgacagcagtagaaagttaagggattgcgatggggcagcaaagtgtgagcaatggttagattcgatttcgcgtcccatttgaatgcacggtgaacgttttcacttagcaccgttatgggtgggaagtgtaggattttggatgggctgcatagcagtcgacgtcaccgtaaaatggaccgccgctatatatcggagatttaacgagccctagagaatggggattataacctatgtcagtcatccatattagacgactcgacgcgagaacctgtgattataggaatccgaacgagttcaaatccaaaaggcggtccgctcaaggccgcctcggctcccagactgtctaaaacgcctgtccctacagtgtcttattatcgcgacgtcattagggatgggaaggtgtacaggcgaaattaccgtggtacacaataaatcataatcctaggtgagcagcgcattctatgtgagtaactaggtgtctaaggtgacggtattgctacggatagggttggtgcggacgatgtgagctagttagtacgcagattgcgaacatttcccggcttacacggctcgagtcgtctggccggtagggaattcttatgtggataaagccgacagtatgcagaaaggttgcattagaaataattggacgcgggttcggtatcgcctcggcgtaacaagagatttatgataatcgtgggtacaaaagcaggtgtcccgggtctgtattgcatggttattcagtcaccaagggctaatatgaatggttgacacgagcacggaaataccagcacagtttca";
				String genome = randomString();
				String mutatedGenome = console.Mutate(genome);
				genomes.add(genome);
				mutatedGenomes.add(mutatedGenome);
	
				this.numMutations = console.getNumberOfMutations();
				numMutationsList.add(numMutations);
				totalNumMutations += numMutations;
	
				if(errorCount <= 5) {
					try {
						Process process2 = Runtime.getRuntime().exec("python3 mutation/g5/approach.py " + genome + " " + mutatedGenome + " " + numMutations + " " + i + " " + m + " 0");
						if(!process2.waitFor(20, TimeUnit.SECONDS)) {
							prevIterationError = 1;
							errorCount++;
						}
						else {
							prevIterationError = 0;
							InputStream is2 = process2.getInputStream();
							BufferedReader br2 = new BufferedReader(new InputStreamReader(is2));
							String l;
							while((l = br2.readLine()) != null) {
								Mutagen mutagen = new Mutagen();
								mutagen.add(l.split("@")[0], l.split("@")[1]);
								this.possibleMutagens.add(mutagen);
								if(!mutagensGuessed.contains(mutagen) && console.testEquiv(mutagen))
									return mutagen;
								mutagensGuessed.add(mutagen);
							}
							br2.close();
							is2.close();
						}
						process2.destroy();
					} catch (Exception e) {
						e.printStackTrace();
					}
				}
	
				int numBasesMutated = 0;
				for(int j = 0; j < genome.length(); j++)
					if(genome.charAt(j) != mutatedGenome.charAt(j))
						numBasesMutated++;
	
				//				System.out.println("Number of mutations: " + numMutations + ", number of bases mutated: " + numBasesMutated);
				totalNumBasesMutated += numBasesMutated;
				numBasesMutatedPerMutation = (int) Math.max(numBasesMutatedPerMutation, Math.ceil(numBasesMutated * 1.0 / numMutations));
	
				List<Pair<Integer, Integer>> intervals = findIntervals(genome, mutatedGenome, numMutations);
				for(Pair<Integer, Integer> interval : intervals) {
					IntervalArchive intervalObj = new IntervalArchive();
					intervalObj.genome = genome;
					intervalObj.mutatedGenome = mutatedGenome;
					int x = interval.getKey();
					int y = interval.getValue();
					intervalObj.hcIntervals.add(new Point(x, y));
					int intervalSize = y + 1 - x;
					if(!intervalSizeTracker.containsKey(intervalSize))
						intervalSizeTracker.put(intervalSize, 0);
					intervalSizeTracker.put(intervalSize, intervalSizeTracker.get(intervalSize) + 1);
					intervalArchives.add(intervalObj);
					if(!multipleRuleTracker.containsKey(numBasesMutatedPerMutation))
						multipleRuleTracker.put(numBasesMutatedPerMutation, new ArrayList<>());
					multipleRuleTracker.get(numBasesMutatedPerMutation).add(intervalObj);
				}
			}
			//			System.out.println("Mutated Base Score: " + totalNumBasesMutated * 1.0 / totalNumMutations);
			System.out.println();
			for(Integer key : intervalSizeTracker.keySet()) {
				System.out.println("Number of intervals of size " + key + ": " + intervalSizeTracker.get(key));
			}
			System.out.println();
	
			//			for(Integer key : multipleRuleTracker.keySet()) {
			//				System.out.println("Number of experiments of size " + key + ": " + multipleRuleTracker.get(key).size());
			//			}
			//			System.out.println();
	
	
			int maxNumIntervals = -1;
			int secondMaxNumIntervals = -1;
			int thirdMaxNumIntervals = -1;
			int intervalSizeMax = -1;
			int intervalSizeSecondMax = -1;
			int intervalSizeThirdMax = -1;
			for(Integer i : intervalSizeTracker.keySet()) {
				int numIntervals = intervalSizeTracker.get(i);
				if(numIntervals > maxNumIntervals) {
					thirdMaxNumIntervals = secondMaxNumIntervals;
					intervalSizeThirdMax = intervalSizeSecondMax;
					secondMaxNumIntervals = maxNumIntervals;
					intervalSizeSecondMax = intervalSizeMax;
					maxNumIntervals = numIntervals;
					intervalSizeMax = i;
				}
				else if(numIntervals > secondMaxNumIntervals) {
					thirdMaxNumIntervals = secondMaxNumIntervals;
					intervalSizeThirdMax = intervalSizeSecondMax;
					secondMaxNumIntervals = numIntervals;
					intervalSizeSecondMax = i;
				}
				else if(numIntervals > thirdMaxNumIntervals) {
					thirdMaxNumIntervals = numIntervals;
					intervalSizeThirdMax = i;
				}
			}
			System.out.println("Max number of intervals - Size: " + intervalSizeMax + ", number of intervals: " + intervalSizeTracker.get(intervalSizeMax));
			System.out.println("2nd max number of intervals - Size: " + intervalSizeSecondMax + ", number of intervals: " + intervalSizeTracker.get(intervalSizeSecondMax));
			System.out.println("3rd max number of intervals - Size: " + intervalSizeThirdMax + ", number of intervals: " + intervalSizeTracker.get(intervalSizeThirdMax));
	
			List<Mutagen> possMutagensForMax = new ArrayList<>();
			List<Mutagen> possMutagensForSecondMax = new ArrayList<>();
			List<Mutagen> possMutagensForThirdMax = new ArrayList<>();
			if(intervalSizeMax >= 1 && intervalSizeMax <= 10) {				
				setUpTrackers();
				implementTrackers(intervalSizeMax);
				possMutagensForMax = getPossibleMutagens(intervalSizeMax);
				for(Mutagen mutagen : possMutagensForMax) {
					System.out.println("**************************************************");
					for(int i = 0; i < mutagen.getPatterns().size(); i++)
						System.out.println(mutagen.getPatterns().get(i) + "@" + mutagen.getActions().get(i));
					System.out.println("**************************************************");
					if(!mutagensGuessed.contains(mutagen) && console.testEquiv(mutagen))
						return mutagen;
					mutagensGuessed.add(mutagen);
				}
				this.possibleMutagens.addAll(possMutagensForMax);
			}
			if(intervalSizeSecondMax >= 1 && intervalSizeSecondMax <= 10) {
				setUpTrackers();
				implementTrackers(intervalSizeSecondMax);
				possMutagensForSecondMax = getPossibleMutagens(intervalSizeSecondMax);
				for(Mutagen mutagen : possMutagensForSecondMax) {
					System.out.println("**************************************************");
					for(int i = 0; i < mutagen.getPatterns().size(); i++)
						System.out.println(mutagen.getPatterns().get(i) + "@" + mutagen.getActions().get(i));
					System.out.println("**************************************************");
					if(!mutagensGuessed.contains(mutagen) && console.testEquiv(mutagen))
						return mutagen;
					mutagensGuessed.add(mutagen);
				}
				this.possibleMutagens.addAll(possMutagensForSecondMax);
			}
			if(intervalSizeThirdMax >= 1 && intervalSizeThirdMax <= 10) {
				setUpTrackers();
				implementTrackers(intervalSizeThirdMax);
				possMutagensForThirdMax = getPossibleMutagens(intervalSizeThirdMax);
				for(Mutagen mutagen : possMutagensForThirdMax) {
					System.out.println("**************************************************");
					for(int i = 0; i < mutagen.getPatterns().size(); i++)
						System.out.println(mutagen.getPatterns().get(i) + "@" + mutagen.getActions().get(i));
					System.out.println("**************************************************");
					if(!mutagensGuessed.contains(mutagen) && console.testEquiv(mutagen))
						return mutagen;
					mutagensGuessed.add(mutagen);
				}
				this.possibleMutagens.addAll(possMutagensForThirdMax);
			}
	
			List<Mutagen> twoRuleMutagens = new ArrayList<>();
			List<Mutagen> twoRuleMutagens12 = new ArrayList<>();
			for(Mutagen mutagenMax : possMutagensForMax) {
				for(Mutagen mutagenSecondMax : possMutagensForSecondMax) {
					Mutagen mutagen = new Mutagen();
					mutagen.add(mutagenMax.getPatterns().get(0), mutagenMax.getActions().get(0));
					mutagen.add(mutagenSecondMax.getPatterns().get(0), mutagenSecondMax.getActions().get(0));
					twoRuleMutagens12.add(mutagen);
					//					System.out.println("**************************************************");
					//					for(int i = 0; i < mutagen.getPatterns().size(); i++)
					//						System.out.println(mutagen.getPatterns().get(i) + "@" + mutagen.getActions().get(i));
					//					System.out.println("**************************************************");
				}
			}
			twoRuleMutagens.addAll(twoRuleMutagens12);
	
			List<Mutagen> twoRuleMutagens13 = new ArrayList<>();
			for(Mutagen mutagenMax : possMutagensForMax) {
				for(Mutagen mutagenThirdMax : possMutagensForThirdMax) {
					Mutagen mutagen = new Mutagen();
					mutagen.add(mutagenMax.getPatterns().get(0), mutagenMax.getActions().get(0));
					mutagen.add(mutagenThirdMax.getPatterns().get(0), mutagenThirdMax.getActions().get(0));
					twoRuleMutagens12.add(mutagen);
					//					System.out.println("**************************************************");
					//					for(int i = 0; i < mutagen.getPatterns().size(); i++)
					//						System.out.println(mutagen.getPatterns().get(i) + "@" + mutagen.getActions().get(i));
					//					System.out.println("**************************************************");
				}
			}
			twoRuleMutagens.addAll(twoRuleMutagens13);
	
			List<Mutagen> twoRuleMutagens23 = new ArrayList<>();
			for(Mutagen mutagenSecondMax : possMutagensForSecondMax) {
				for(Mutagen mutagenThirdMax : possMutagensForThirdMax) {
					Mutagen mutagen = new Mutagen();
					mutagen.add(mutagenSecondMax.getPatterns().get(0), mutagenSecondMax.getActions().get(0));
					mutagen.add(mutagenThirdMax.getPatterns().get(0), mutagenThirdMax.getActions().get(0));
					twoRuleMutagens23.add(mutagen);
					//					System.out.println("**************************************************");
					//					for(int i = 0; i < mutagen.getPatterns().size(); i++)
					//						System.out.println(mutagen.getPatterns().get(i) + "@" + mutagen.getActions().get(i));
					//					System.out.println("**************************************************");
				}
			}
			twoRuleMutagens.addAll(twoRuleMutagens23);
	
			for(Mutagen mutagen : twoRuleMutagens) {
				if(!mutagensGuessed.contains(mutagen) && console.testEquiv(mutagen))
					return mutagen;
				mutagensGuessed.add(mutagen);
			}
			this.possibleMutagens.addAll(twoRuleMutagens);
	
			List<Mutagen> threeRuleMutagens = new ArrayList<>();
			for(Mutagen mutagenMax : possMutagensForMax) {
				for(Mutagen mutagenSecondMax : possMutagensForSecondMax) {
					for(Mutagen mutagenThirdMax : possMutagensForThirdMax) {
						Mutagen mutagen = new Mutagen();
						mutagen.add(mutagenMax.getPatterns().get(0), mutagenMax.getActions().get(0));
						mutagen.add(mutagenSecondMax.getPatterns().get(0), mutagenSecondMax.getActions().get(0));
						mutagen.add(mutagenThirdMax.getPatterns().get(0), mutagenThirdMax.getActions().get(0));
						threeRuleMutagens.add(mutagen);
						//						System.out.println("**************************************************");
						//						for(int i = 0; i < mutagen.getPatterns().size(); i++)
						//							System.out.println(mutagen.getPatterns().get(i) + "@" + mutagen.getActions().get(i));
						//						System.out.println("**************************************************");
					}
				}
			}
			for(Mutagen mutagen : threeRuleMutagens) {
				if(!mutagensGuessed.contains(mutagen) && console.testEquiv(mutagen))
					return mutagen;
				mutagensGuessed.add(mutagen);
			}
			this.possibleMutagens.addAll(threeRuleMutagens);
	
	//		if(maxNumIntervals >= 100 && secondMaxNumIntervals < 100 && thirdMaxNumIntervals < 100)
	//			if(possMutagensForMax.size() > 0)
	//				return possMutagensForMax.get(0);
	//		if(maxNumIntervals < 100 && secondMaxNumIntervals >= 100 && thirdMaxNumIntervals < 100)
	//			if(possMutagensForSecondMax.size() > 0)
	//				return possMutagensForSecondMax.get(0);
	//		if(maxNumIntervals < 100 && secondMaxNumIntervals < 100 && thirdMaxNumIntervals >= 100)
	//			if(possMutagensForThirdMax.size() > 0)
	//				return possMutagensForThirdMax.get(0);
	//		if(maxNumIntervals >= 100 && secondMaxNumIntervals >= 100 && thirdMaxNumIntervals < 100)
	//			if(twoRuleMutagens12.size() > 0)
	//				return twoRuleMutagens12.get(0);
	//		if(maxNumIntervals >= 100 && secondMaxNumIntervals < 100 && thirdMaxNumIntervals >= 100)
	//			if(twoRuleMutagens13.size() > 0)
	//				return twoRuleMutagens13.get(0);
	//		if(maxNumIntervals < 100 && secondMaxNumIntervals >= 100 && thirdMaxNumIntervals >= 100)
	//			if(twoRuleMutagens23.size() > 0)
	//				return twoRuleMutagens23.get(0);
	//		if(maxNumIntervals >= 100 && secondMaxNumIntervals >= 100 && thirdMaxNumIntervals >= 100)
	//			if(threeRuleMutagens.size() > 0)
	//				return threeRuleMutagens.get(0);

		}
		
//		for(Mutagen mutagen : this.possibleMutagens)
//			return mutagen;
//
//		
//		Mutagen mutagen = new Mutagen();
//		mutagen.add("a;c;c", "att");
//		mutagen.add("g;c;c", "gtt");
//		return mutagen;
	}

	public ArrayList<Pair<Integer, Integer>> findIntervals(String genome, String result, Integer mp){
		ArrayList<Pair<Integer, Integer>> intervals = new ArrayList<>();
		if (mp == 0) return intervals;
		ArrayList<Integer> diffs = new ArrayList<>();
		for (int i=0; i<genome.length(); i++){
			if(genome.charAt(i) != result.charAt(i))
				diffs.add(i);
		}
		Integer startidx = -1;
		Integer curidx = -1;
		for(Integer d: diffs){
			if(startidx == -1){
				startidx = d;
				curidx = d;
			}
			else if (d - curidx < 10){
				curidx = d;
			}
			else{
				intervals.add(new Pair<Integer, Integer>(startidx, curidx));
			}
		}
		if (startidx != -1 && curidx != -1)
			intervals.add(new Pair<Integer, Integer>(startidx, curidx));
		return intervals;
	}

	private List<Mutagen> getPossibleMutagens(int numBasesMutatedPerMutation) {
		if(numBasesMutatedPerMutation == 1)
			return getOneBaseMutationMutagens();
		return getMultipleBaseMutationMutagens(numBasesMutatedPerMutation);
	}

	private List<Mutagen> getOneBaseMutationMutagens() {

		System.out.println("Action tracker: ");
		for(int j = 0; j < actionTracker.size(); j++)
			System.out.println("  " + j + ". " + actionTracker.get(j));
		System.out.println("Pattern tracker: ");
		for(int j = 0; j < patternTracker.size(); j++)
			System.out.println("  " + j + ". " + patternTracker.get(j));
		System.out.println();

		List<Mutagen> mutagens = new ArrayList<>();
		Map<Integer, List<String>> interestingPatternsMap = new TreeMap<>();
		for(int i = 0; i < patternTracker.size(); i++) {
			Map<Character, Integer> patternOccurrences = patternTracker.get(i);
			List<String> interestingPattern = new ArrayList<>();
			int totalCharacters = 0;
			for(Character character : patternOccurrences.keySet())
				totalCharacters += patternOccurrences.get(character);
			for(Character character : patternOccurrences.keySet()) {
				if((patternOccurrences.get(character) * 1.0 / totalCharacters) > PERCENT_NONEXISTENT_BASE_THRESHOLD)
					interestingPattern.add(character + "");
			}
			if(interestingPattern.size() != 4) {
				interestingPatternsMap.put(i, interestingPattern);
			}
		}

		Map<Character, Integer> actionOccurrences = actionTracker.get(0);
		List<Integer> possibleLocationsForAction = new ArrayList<>();
		char letterMutation = getLetterMutation(actionOccurrences);    		
		if(letterMutation != 'x') {
			for(int i = 0; i < patternTracker.size(); i++) {
				Map<Character, Integer> patternOccurrences = patternTracker.get(i);
				boolean isValidLocationForAction = true;
				for(Character character : patternOccurrences.keySet()) {
					if((character != letterMutation && patternOccurrences.get(character) == 0) ||
							(character == letterMutation && patternOccurrences.get(character) != 0)) {
						isValidLocationForAction = false;
						break;
					}
				}
				if(isValidLocationForAction)
					possibleLocationsForAction.add(i);
			}

			/*
			 * Determine pattern
			 */
			String pattern = "";
			int offsetForBaseMutationInAction = 0;
			if(interestingPatternsMap.size() == 0) {
				pattern += "acgt";
			}
			else {
				/*
				 * If we have acgt;acgt;gt;acgt;acgt;act@012c, then
				 * "offsetForBaseMutationInAction" would be 3, since "c"
				 * is in offset 3 in the action.
				 */
				List<Integer> interestingPatternLocations = new ArrayList<>();
				for(Integer interestingPatternLocation : interestingPatternsMap.keySet())
					interestingPatternLocations.add(interestingPatternLocation);
				int currentLocation = interestingPatternLocations.get(0);
				offsetForBaseMutationInAction = 9 - currentLocation;

				for(int i = 0; i < interestingPatternLocations.size(); i++) {
					int nextLocation = interestingPatternLocations.get(i);
					while(currentLocation < nextLocation) {
						pattern += "acgt;";
						currentLocation++;
					}
					if(i == interestingPatternLocations.size() - 1)
						pattern += String.join("", interestingPatternsMap.get(nextLocation));
					else
						pattern += String.join("", interestingPatternsMap.get(nextLocation)) + ";";
					currentLocation++;
				}
			}
			/*
			 * Determine action
			 */
			for(int i = 0; i < 10 - offsetForBaseMutationInAction; i++) {
				String newPattern = pattern;
				if(newPattern.length() != 4 || newPattern.contains(";"))
					for(int j = 0; j < i; j++)
						newPattern = "acgt;" + newPattern;

				String action = "";
				for(int j = 0; j < offsetForBaseMutationInAction + i; j++)
					action += Integer.toString(j);
				action += letterMutation;

				String[] patternElements = newPattern.split(";");
				if(patternElements.length > 10) {
					newPattern = patternElements[patternElements.length - 10];
					for(int abc = patternElements.length - 9; abc <= patternElements.length - 1; abc++) {
						newPattern += ";" + patternElements[abc];
					}
				}
				if(action.length() > 10)
					action = action.substring(action.length() - 10, action.length());

				System.out.println("Pattern: " + newPattern + ", action: " + action);

				//				System.out.println("Predicted pattern: " + newPattern);
				//				System.out.println("Predicted action: " + action);
				//				System.out.println("Predicted mutated base: " + letterMutation);
				//				System.out.println("Offset for base mutation in action: " + offsetForBaseMutationInAction);
				//				System.out.println("Predicted rule: " + newPattern + "@" + action);
				//				System.out.println();

				Mutagen mutagen = new Mutagen();
				mutagen.add(newPattern, action);
				mutagens.add(mutagen);
			}
		}

		// Either the mutation is numerical, or the mutation could be either base or numerical
		for(int i = 0; i < patternTracker.size(); i++) {
			Map<Character, Integer> patternOccurrences = patternTracker.get(i);
			double distance = getDistance(actionOccurrences, patternOccurrences);
			if(distance <= (DISTANCE_THRESHOLD_SB * m * 1.0 / 15))
				possibleLocationsForAction.add(i);
		}

		/*
		 * Determine pattern
		 */
		String pattern = "";
		int locationForFirstPattern = 0;
		boolean firstPatternFound = false;
		if(interestingPatternsMap.size() == 0)
			pattern += "acgt";
		else {
			/*
			 * If we have acgt;acgt;gt;acgt;acgt;act@012c, then
			 * "offsetForBaseMutationInAction" would be 3, since "c"
			 * is in offset 3 in the action.
			 */
			List<Integer> interestingPatternLocations = new ArrayList<>();
			for(Integer interestingPatternLocation : interestingPatternsMap.keySet())
				interestingPatternLocations.add(interestingPatternLocation);
			int currentLocation = interestingPatternLocations.get(0);
			locationForFirstPattern = currentLocation;
			firstPatternFound = true;

			for(int i = 0; i < interestingPatternLocations.size(); i++) {
				int nextLocation = interestingPatternLocations.get(i);
				while(currentLocation < nextLocation) {
					pattern += "acgt;";
					currentLocation++;
				}
				if(i == interestingPatternLocations.size() - 1)
					pattern += String.join("", interestingPatternsMap.get(nextLocation));
				else
					pattern += String.join("", interestingPatternsMap.get(nextLocation)) + ";";
				currentLocation++;
			}
		}

		/*
		 * Determine action
		 */
		for(int i = 0; i < possibleLocationsForAction.size(); i++) {
			String intermediatePattern = pattern;
			int locationForAction = possibleLocationsForAction.get(i);
			if(intermediatePattern.length() != 4 || intermediatePattern.contains(";"))
				if(locationForAction < locationForFirstPattern)
					for(int j = 0; j < locationForFirstPattern - locationForAction; j++)
						intermediatePattern = "acgt;" + intermediatePattern;

			int indexOfFirstInterestingPattern = 0;
			if(!firstPatternFound)
				locationForFirstPattern = locationForAction;
			if(interestingPatternsMap.size() != 0) {
				String firstInterestingPattern = String.join("", interestingPatternsMap.get(locationForFirstPattern));
				indexOfFirstInterestingPattern = Arrays.asList(intermediatePattern.split(";")).indexOf(firstInterestingPattern);					
			}
			int offsetForBaseMutationInAction = 9 + indexOfFirstInterestingPattern - locationForFirstPattern;
			int numberMutation = offsetForBaseMutationInAction + locationForAction - 9;

			for(int j = 0; j < 10 - offsetForBaseMutationInAction; j++) {
				String newPattern = intermediatePattern;

				if(numberMutation + j > 9)
					continue;

				if(newPattern.length() != 4 || newPattern.contains(";"))
					for(int k = 0; k < j; k++)
						newPattern = "acgt;" + newPattern;

				String action = "";
				for(int k = 0; k < offsetForBaseMutationInAction + j; k++) {
					action += Integer.toString(k);
				}					
				action += Integer.toString(numberMutation + j);

				String[] patternElements = newPattern.split(";");
				if(patternElements.length > 10) {
					newPattern = patternElements[patternElements.length - 10];
					for(int abc = patternElements.length - 9; abc <= patternElements.length - 1; abc++) {
						newPattern += ";" + patternElements[abc];
					}
				}
				if(action.length() > 10)
					action = action.substring(action.length() - 10, action.length());

				System.out.println("Pattern: " + newPattern + ", action: " + action);
				System.out.println();

				System.out.println("Predicted pattern: " + newPattern);
				System.out.println("Predicted action: " + action);
				System.out.println("Predicted number mutation: " + numberMutation);
				System.out.println("Location of action: " + locationForAction);
				System.out.println("Index of first interesting pattern: " + indexOfFirstInterestingPattern);
				System.out.println("Location for first interesting pattern: " + locationForFirstPattern);
				System.out.println("Offset for base mutation in action: " + offsetForBaseMutationInAction);
				System.out.println("Predicted rule: " + newPattern + "@" + action);

				Mutagen mutagen = new Mutagen();
				mutagen.add(newPattern, action);
				mutagens.add(mutagen);
			}
		}

		List<Mutagen> newMutagens = new ArrayList<>();
		for(Mutagen mutagen : mutagens) {
			String mutagenPattern = mutagen.getPatterns().get(0);
			List<String> oldPatternElementsAsList = Arrays.asList(mutagenPattern.split(";"));
			List<String> patternElementsAsList = new ArrayList<String>(oldPatternElementsAsList);
			while(patternElementsAsList.size() > 1) {
				if(patternElementsAsList.get(patternElementsAsList.size() - 1).length() >= 3) {
					int lastIndex = patternElementsAsList.size() - 1;
					patternElementsAsList.remove(lastIndex);
					String newMutagenPattern = String.join(";", patternElementsAsList);
					Mutagen newMutagen = new Mutagen();
					newMutagen.add(newMutagenPattern, mutagen.getActions().get(0));
					newMutagens.add(newMutagen);

					//					System.out.println("Predicted pattern: " + newMutagenPattern);
					//					System.out.println("Predicted action: " + mutagen.getActions().get(0));
					//					System.out.println("Predicted rule: " + newMutagenPattern + "@" + mutagen.getActions().get(0));
					//					System.out.println();
				}
				else
					break;
			}
		}
		mutagens.addAll(newMutagens);
		return mutagens;
	}

	private List<Mutagen> getMultipleBaseMutationMutagens(int mutationSize) {
		List<Mutagen> mutagens = new ArrayList<>();
		Map<Integer, List<String>> interestingPatternsMap = new TreeMap<>();
		for(int i = 0; i < patternTracker.size(); i++) {
			Map<Character, Integer> patternOccurrences = patternTracker.get(i);
			List<String> interestingPattern = new ArrayList<>();

			String letter = "";
			if(!(letter = getLetterMutationForMultipleBases(patternOccurrences) + "").equals("x"))
				interestingPattern.add(letter);
			else {
				int totalCharacters = 0;
				for(Character character : patternOccurrences.keySet())
					totalCharacters += patternOccurrences.get(character);
				for(Character character : patternOccurrences.keySet()) {
					if((patternOccurrences.get(character) * 1.0 / totalCharacters) > PERCENT_NONEXISTENT_BASE_THRESHOLD)
						interestingPattern.add(character + "");
				}
			}

			if(interestingPattern.size() != 4)
				interestingPatternsMap.put(i, interestingPattern);
		}

		/*
		 * Determine pattern
		 */
		String pattern = "";
		int offsetForBaseMutationInAction = 0;
		int locationForFirstPattern = 0;
		boolean firstPatternFound = false;
		if(interestingPatternsMap.size() == 0) {
			pattern += "acgt";
		}
		else {
			/*
			 * If we have acgt;acgt;gt;acgt;acgt;act@012c, then
			 * "offsetForBaseMutationInAction" would be 3, since "c"
			 * is in offset 3 in the action.
			 */
			List<Integer> interestingPatternLocations = new ArrayList<>();
			for(Integer interestingPatternLocation : interestingPatternsMap.keySet())
				interestingPatternLocations.add(interestingPatternLocation);
			int currentLocation = interestingPatternLocations.get(0);
			offsetForBaseMutationInAction = 9 + 1 - mutationSize - currentLocation;
			locationForFirstPattern = currentLocation;
			firstPatternFound = true;

			for(int i = 0; i < Math.min(interestingPatternLocations.size(), 10); i++) {
				int nextLocation = interestingPatternLocations.get(i);
				while(currentLocation < nextLocation) {
					pattern += "acgt;";
					currentLocation++;
				}
				if(interestingPatternsMap.get(nextLocation).size() == 0) {
					List<String> interestingPattern = new ArrayList<>();
					interestingPattern.add("a");
					interestingPattern.add("c");
					interestingPattern.add("t");
					interestingPattern.add("g");
					interestingPatternsMap.put(nextLocation, interestingPattern);
				}

				if(i == Math.min(interestingPatternLocations.size(), 10) - 1)
					pattern += String.join("", interestingPatternsMap.get(nextLocation));
				else
					pattern += String.join("", interestingPatternsMap.get(nextLocation)) + ";";
				currentLocation++;
			}
		}

		String[] patternElements = pattern.split(";");
		if(patternElements.length > 10) {
			pattern = patternElements[patternElements.length - 10];
			for(int abc = patternElements.length - 9; abc <= patternElements.length - 1; abc++) {
				pattern += ";" + patternElements[abc];
			}
		}

		System.out.println("Pattern: " + pattern);
		System.out.println();

		List<List<String>> actionElementsToAdd = new ArrayList<>();
		for(int xyz = 0; xyz < mutationSize; xyz++) {
			actionElementsToAdd.add(new ArrayList<>());
			Map<Character, Integer> actionOccurrences = actionTracker.get(xyz);
			List<Integer> possibleLocationsForAction = new ArrayList<>();
			char letterMutation = getLetterMutationForMultipleBases(actionOccurrences);    		
			if(letterMutation != 'x') {
				for(int i = 0; i < patternTracker.size(); i++) {
					Map<Character, Integer> patternOccurrences = patternTracker.get(i);
					boolean isValidLocationForAction = true;
					for(Character character : patternOccurrences.keySet()) {
						if((character != letterMutation && patternOccurrences.get(character) == 0) ||
								(character == letterMutation && patternOccurrences.get(character) != 0)) {
							isValidLocationForAction = false;
							break;
						}
					}
					if(isValidLocationForAction)
						possibleLocationsForAction.add(i);
				}


				/*
				 * Add action element
				 */
				if(xyz == 0) {
					String actionElement = "";
					for(int i = 0; i < offsetForBaseMutationInAction; i++)
						actionElement += Integer.toString(i);
					actionElement += letterMutation + "";
					actionElementsToAdd.get(xyz).add(actionElement);
				}
				else
					actionElementsToAdd.get(xyz).add(letterMutation + "");

				//				System.out.println("Action: " + xyz);
				//				System.out.println("Location for first pattern: " + locationForFirstPattern);
				//				System.out.println("Offset for base mutation in action: " + offsetForBaseMutationInAction);
				//				System.out.println("Letter mutation: " + letterMutation);
				//				System.out.println();

			}

			// Either the mutation is numerical, or the mutation could be either base or numerical
			//			for(int i = 0; i < patternTracker.size(); i++) {
			//				Map<Character, Integer> patternOccurrences = patternTracker.get(i);
			//				double distance = getDistance(actionOccurrences, patternOccurrences);
			//				int totalOccurrences = 0;
			//				for(Character character : patternOccurrences.keySet())
			//					totalOccurrences += patternOccurrences.get(character);
			//				double distancePercentage = distance / totalOccurrences;
			//				if(distancePercentage <= PERCENTAGE_DISTANCE_THRESHOLD_MB)
			//					possibleLocationsForAction.add(i);
			//			}
			double minDistance = Double.MAX_VALUE;
			double secondMinDistance = Double.MAX_VALUE;
			double thirdMinDistance = Double.MAX_VALUE;
			int patternInWindowMin = -1;
			int patternInWindowSecondMin = -1;
			int patternInWindowThirdMin = -1;
			for(int i = 0; i < patternTracker.size(); i++) {
				Map<Character, Integer> patternOccurrences = patternTracker.get(i);
				double distance = getDistance(actionOccurrences, patternOccurrences);
				if(distance < minDistance) {
					thirdMinDistance = secondMinDistance;
					patternInWindowThirdMin = patternInWindowSecondMin;
					secondMinDistance = minDistance;
					patternInWindowSecondMin = patternInWindowMin;
					minDistance = distance;
					patternInWindowMin = i;
				}
				else if(distance < secondMinDistance) {
					thirdMinDistance = secondMinDistance;
					patternInWindowThirdMin = patternInWindowSecondMin;
					secondMinDistance = distance;
					patternInWindowSecondMin = i;
				}
				else if(distance < thirdMinDistance) {
					thirdMinDistance = distance;
					patternInWindowThirdMin = i;
				}
			}
			possibleLocationsForAction.add(patternInWindowMin);
			possibleLocationsForAction.add(patternInWindowSecondMin);
			//			possibleLocationsForAction.add(patternInWindowThirdMin);

			/*
			 * Add action element
			 */
			for(int i = 0; i < possibleLocationsForAction.size(); i++) {
				String intermediatePattern = pattern;
				int locationForAction = possibleLocationsForAction.get(i);
				if(intermediatePattern.length() != 4 || intermediatePattern.contains(";"))
					if(locationForAction < locationForFirstPattern)
						for(int j = 0; j < locationForFirstPattern - locationForAction; j++)
							intermediatePattern = "acgt;" + intermediatePattern;

				int indexOfFirstInterestingPattern = 0;
				if(!firstPatternFound)
					locationForFirstPattern = locationForAction;
				if(interestingPatternsMap.size() != 0) {
					String firstInterestingPattern = String.join("", interestingPatternsMap.get(locationForFirstPattern));
					indexOfFirstInterestingPattern = Arrays.asList(intermediatePattern.split(";")).indexOf(firstInterestingPattern);
				}
				int offsetForNumberMutationInAction = (9 + 1 - mutationSize) + indexOfFirstInterestingPattern - locationForFirstPattern;
				int numberMutation = offsetForNumberMutationInAction + locationForAction - (9 + 1 - mutationSize);

				//				System.out.println("Action: " + xyz);
				//				System.out.println("Possible location for action: " + locationForAction);
				//				System.out.println("Index of first interesting pattern: " + indexOfFirstInterestingPattern);
				//				System.out.println("Location for first pattern: " + locationForFirstPattern);
				//				System.out.println("Offset for number mutation in action: " + offsetForNumberMutationInAction);
				//				System.out.println("Number mutation: " + numberMutation);
				//				System.out.println();

				if(!actionElementsToAdd.get(xyz).contains(Integer.toString(numberMutation)) && numberMutation <= 9 && numberMutation >= 0) {
					if(xyz == 0) {
						String actionElement = "";
						for(int j = 0; j < offsetForNumberMutationInAction; j++)
							actionElement += Integer.toString(j);
						actionElement += Integer.toString(numberMutation);
						if(actionElementsToAdd.get(xyz).size() < 2)
							actionElementsToAdd.get(xyz).add(actionElement);
					}
					else if(actionElementsToAdd.get(xyz).size() < 2)
						actionElementsToAdd.get(xyz).add(Integer.toString(numberMutation));	
				}
			}
		}
		List<String> actionStrings = new ArrayList<>();
		if(actionElementsToAdd.size() != 0)
			getActionStringsForMultipleBases(actionElementsToAdd, actionStrings, "", 0);
		//		System.out.println(actionElementsToAdd);
		for(String action : actionStrings) {
			Mutagen mutagen = new Mutagen();
			mutagen.add(pattern, action);
			mutagens.add(mutagen);
		}

		List<Mutagen> newMutagens = new ArrayList<>();
		for(Mutagen mutagen : mutagens) {
			String mutagenPattern = mutagen.getPatterns().get(0);
			List<String> oldPatternElementsAsList = Arrays.asList(mutagenPattern.split(";"));
			List<String> patternElementsAsList = new ArrayList<String>(oldPatternElementsAsList);
			while(patternElementsAsList.size() > 1) {
				if(patternElementsAsList.get(patternElementsAsList.size() - 1).length() >= 3) {
					int lastIndex = patternElementsAsList.size() - 1;
					patternElementsAsList.remove(lastIndex);
					String newMutagenPattern = String.join(";", patternElementsAsList);
					Mutagen newMutagen = new Mutagen();
					newMutagen.add(newMutagenPattern, mutagen.getActions().get(0));
					newMutagens.add(newMutagen);

					//					System.out.println("Predicted pattern: " + newMutagenPattern);
					//					System.out.println("Predicted action: " + mutagen.getActions().get(0));
					//					System.out.println("Predicted rule: " + newMutagenPattern + "@" + mutagen.getActions().get(0));
					//					System.out.println();
				}
				else
					break;
			}
		}
		mutagens.addAll(newMutagens);
		return mutagens;
	}

	private double getDistance(Map<Character, Integer> set1,
			Map<Character, Integer> set2) {
		double distance = 0.0;
		for(Character character : set1.keySet())
			distance += Math.pow(set1.get(character) - set2.get(character), 2);
		return Math.sqrt(distance);
	}

	private void getActionStringsForMultipleBases(List<List<String>> actionElementsToAdd, List<String> actionStrings, String currString, int listNum) {
		for(int i = 0; i < actionElementsToAdd.get(listNum).size(); i++) {
			if(listNum == actionElementsToAdd.size() - 1)
				actionStrings.add(currString + actionElementsToAdd.get(listNum).get(i));
			else {
				String newString = currString + actionElementsToAdd.get(listNum).get(i);
				int newListNum = listNum + 1;
				getActionStringsForMultipleBases(actionElementsToAdd, actionStrings, newString, newListNum);
			}
		}
	}

	private void setUpStructures() {
		numMatchesForMutationSizes = new ArrayList<>();
		for(int i = 0; i < 10; i++)
			numMatchesForMutationSizes.add(0);
		genomes = new ArrayList<>();
		mutatedGenomes = new ArrayList<>();
		intervalArchives = new ArrayList<>();
		possibleMutagens = new ArrayList<>();
		multipleRuleTracker = new TreeMap<>();
		intervalSizeTracker = new TreeMap<>();
	}

	private void setUpTrackers() {
		patternTracker = new ArrayList<>();
		actionTracker = new ArrayList<>();

		for(int i = 0; i < 19; i++) {
			Map<Character, Integer> occurrences = new HashMap<>();
			occurrences.put('a', 0);
			occurrences.put('c', 0);
			occurrences.put('g', 0);
			occurrences.put('t', 0);
			patternTracker.add(occurrences);			
		}
		for(int i = 0; i < 10; i++) {
			Map<Character, Integer> occurrences = new HashMap<>();
			occurrences.put('a', 0);
			occurrences.put('c', 0);
			occurrences.put('g', 0);
			occurrences.put('t', 0);
			actionTracker.add(occurrences);			
		}
	}

	private Character getLetterMutation(Map<Character, Integer> occurrences) {
		int numAs = occurrences.get('a');
		int numCs = occurrences.get('c');
		int numGs = occurrences.get('g');
		int numTs = occurrences.get('t');
		if(numAs != 0 && numCs == 0 && numGs == 0 && numTs == 0)
			return 'a';
		if(numAs == 0 && numCs != 0 && numGs == 0 && numTs == 0)
			return 'c';
		if(numAs == 0 && numCs == 0 && numGs != 0 && numTs == 0)
			return 'g';
		if(numAs == 0 && numCs == 0 && numGs == 0 && numTs != 0)
			return 't';
		return 'x';
	}

	private Character getLetterMutationForMultipleBases(Map<Character, Integer> occurrences) {
		int numAs = occurrences.get('a');
		int numCs = occurrences.get('c');
		int numGs = occurrences.get('g');
		int numTs = occurrences.get('t');
		if((numAs != 0) && (numCs <= 0.2 * numAs) && (numGs <= 0.2 * numAs) && (numTs <= 0.2 * numAs))
			return 'a';
		if((numAs <= 0.2 * numCs) && (numCs != 0) && (numGs <= 0.2 * numCs) && (numTs <= 0.2 * numCs))
			return 'c';
		if((numAs <= 0.2 * numGs) && (numCs <= 0.2 * numGs) && (numGs != 0) && (numTs <= 0.2 * numGs))
			return 'g';
		if((numAs <= 0.2 * numTs) && (numCs <= 0.2 * numTs) && (numGs <= 0.2 * numTs) && (numTs != 0))
			return 't';
		return 'x';
	}

	private void implementTrackers(int mutationSize) {
		for(IntervalArchive intervalObj : intervalArchives) {
			List<Window> windows = createWindows(intervalObj, mutationSize);
			//			System.out.println("Number of windows: " + windows.size());
			//			for(int i = 0; i < windows.size(); i++) {
			//				System.out.println("Window " + (i + 1));
			//				Window window = windows.get(i);
			//				window.print();
			//			}

			for(Window window : windows) {
				int start = window.windowCoordinates.x;
				int baseLocation = window.baseLocation;
				baseLocation += start > baseLocation ? 1000 : 0;
				String genomeValue = "", mutatedGenomeValue = "";
				for(int i = 0; i < patternTracker.size(); i++) {
					Map<Character, Integer> occurrences = patternTracker.get(i);
					Character genomeBase = intervalObj.genome.charAt((start + i) % intervalObj.genome.length());
					Character mutatedGenomeBase = intervalObj.mutatedGenome.charAt((start + i) % intervalObj.mutatedGenome.length());
					genomeValue += genomeBase;
					mutatedGenomeValue += mutatedGenomeBase;
					Character baseToIncrement;
					if(start + i < baseLocation)
						baseToIncrement = mutatedGenomeBase;
					else {
						if(mutationSize == 1) {
							if(start + i == baseLocation)
								actionTracker.get(0).put(mutatedGenomeBase, actionTracker.get(0).get(mutatedGenomeBase) + 1);
						}
						else if(start + i >= baseLocation && start + i < baseLocation + mutationSize) {
							int index = start + i - baseLocation + mutationSize - window.mutatedBases.size();
							actionTracker.get(index).put(mutatedGenomeBase, actionTracker.get(index).get(mutatedGenomeBase) + 1);
						}
						baseToIncrement = genomeBase;
					}
					occurrences.put(baseToIncrement, occurrences.get(baseToIncrement) + 1);
				}
				//				System.out.println("Genome value: " + genomeValue);
				//				System.out.println("Mutated genome value: " + mutatedGenomeValue);
			}
		}
	}

	private List<Window> createWindows(IntervalArchive intervalObj, int mutationSize) {
		String genome = intervalObj.genome;
		String mutatedGenome = intervalObj.mutatedGenome;
		List<Point> hcIntervals = intervalObj.hcIntervals;

		List<Window> windows = new ArrayList<>();
		if(mutationSize == 1) {
			for(Point hcInterval : hcIntervals) {
				if(hcInterval.x != hcInterval.y)
					continue;
				int i = hcInterval.x;
				if(genome.charAt(i) != mutatedGenome.charAt(i)) {
					Window window = new Window();
					window.originalBase = genome.charAt(i);
					window.mutatedBase = mutatedGenome.charAt(i);
					window.baseLocation = i;
					window.windowCoordinates.x = i < 9 ? genome.length() - (9 - i) : i - 9;
					window.windowCoordinates.y = i > genome.length() - 10 ? 9 - (genome.length() - i) : i + 9;
					windows.add(window);
				}
			}
		}
		else {
			List<Point> newIntervals = new ArrayList<>();
			for(int i = 0; i < hcIntervals.size(); i++)
				if((hcIntervals.get(i).y + 1 - hcIntervals.get(i).x) == mutationSize)
					newIntervals.add(hcIntervals.get(i));

			List<Point> cyclicIntervals = new ArrayList<>(hcIntervals);
			for(int i = 0; i < hcIntervals.size() - 1; i++)
				cyclicIntervals.add(hcIntervals.get(i));

			for(int i = 0; i < hcIntervals.size(); i++) {
				for(int j = 1; j < hcIntervals.size(); j++) {
					int size1 = cyclicIntervals.get(i + j).y + 1 - cyclicIntervals.get(i).x;
					int size2 = cyclicIntervals.get(i + j).y + 1001 - cyclicIntervals.get(i).x;
					Point point = new Point(cyclicIntervals.get(i).x, cyclicIntervals.get(i + j).y);
					if(((size1 == mutationSize)  || (size2 == mutationSize)) && !newIntervals.contains(point))
						newIntervals.add(point);
				}
			}
			//			System.out.println("*******************************************");
			//			System.out.println("HC Intervals: " + hcIntervals);
			//			System.out.println();
			//			System.out.println("New Intervals: " + newIntervals);
			//			System.out.println("*******************************************");

			for(Point newInterval : newIntervals) {
				int i = newInterval.x;
				int j = newInterval.y;
				if(j < i)
					j += 1000;

				String newOriginalGenome = genome + genome;
				String newMutatedGenome = mutatedGenome + mutatedGenome;

				String originalSubstring = newOriginalGenome.substring(i, j + 1);
				List<Character> originalBases = new ArrayList<>();
				for(int k = 0; k < originalSubstring.length(); k++)
					originalBases.add(originalSubstring.charAt(k));
				String mutatedSubstring = newMutatedGenome.substring(i, j + 1);
				List<Character> mutatedBases = new ArrayList<>();
				for(int k = 0; k < mutatedSubstring.length(); k++)
					mutatedBases.add(mutatedSubstring.charAt(k));

				Window window = new Window();
				window.originalBases = originalBases;
				window.mutatedBases = mutatedBases;
				window.baseLocation = i;
				window.windowCoordinates.x = i < (10 - mutationSize) ? 
						genome.length() - ((10 - mutationSize) - i) : 
							i - (10 - mutationSize);
						window.windowCoordinates.y = i > genome.length() - 10 ? 9 - (genome.length() - i) : i + 9;
						windows.add(window);
						//				window.print();
						//				System.out.println();
			}
		}

		/*
		 * Examples of overlapping windows:
		 * 
		 * (1, 19) and (4, 22)
		 * (974, 992) and (986, 4)
		 * (986, 4) and (992, 10)
		 * (992, 10) and (4, 22)
		 * 
		 * Examples of non-overlapping windows:
		 * 
		 * (2, 20) and (34, 52)
		 * (974, 992) and (995, 13)
		 * (995, 13) and (34, 52)
		 */    					
		for(int i = 0; i < windows.size(); i++) {
			for(int j = 0; j < windows.size(); j++) {
				if(windows.get(i).equals(windows.get(j)))
					continue;
				int x_i = windows.get(i).windowCoordinates.x;
				int y_i = windows.get(i).windowCoordinates.y;
				int x_j = windows.get(j).windowCoordinates.x;
				int y_j = windows.get(j).windowCoordinates.y;

				if((x_i > x_j && x_i < y_j) || (y_i > x_j && y_i < y_j) ||
						(x_j > x_i && x_j < y_i) || (y_j > x_i && y_j < y_i))
					windows.get(i).overlappingWindows.add(windows.get(j));
				else {
					y_i += (x_i > y_i) ? 1000 : 0;
					y_j += (x_j > y_j) ? 1000 : 0;    				
					if((x_i > x_j && x_i < y_j) || (y_i > x_j && y_i < y_j) ||
							(x_j > x_i && x_j < y_i) || (y_j > x_i && y_j < y_i))
						windows.get(i).overlappingWindows.add(windows.get(j));
				}	
			}
		}
		return windows;
	}

	class Window {
		public char originalBase;
		public List<Character> originalBases;
		public char mutatedBase;
		public List<Character> mutatedBases;
		public int baseLocation;
		public Point windowCoordinates = new Point();
		public List<Window> overlappingWindows = new ArrayList<>();

		public void print() {
			//			System.out.println("(1) Original base: " + originalBase);
			//			System.out.println("(2) Mutated base: " + mutatedBase);
			System.out.println("(3) Original bases: " + originalBases);
			System.out.println("(4) Mutated bases: " + mutatedBases);
			System.out.println("(5) Base location: " + baseLocation);
			System.out.println("(6) Window coordinates: (" + windowCoordinates.x + ", " + windowCoordinates.y + ")");
			//			System.out.println("(7) Number of overlapping windows: " + overlappingWindows.size());
		}
	}

	class IntervalArchive {
		public String genome;
		public String mutatedGenome;
		public List<Point> hcIntervals = new ArrayList<>();

		public void print() {
			//			System.out.println("(1) Genome: " + genome);
			//			System.out.println("(2) Mutated genome: " + mutatedGenome);
			System.out.println("(3) Intervals: " + hcIntervals);
		}
	}
}
