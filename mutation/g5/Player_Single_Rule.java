package mutation.g5;

import java.awt.Point;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.*;
import mutation.sim.Console;
import mutation.sim.Mutagen;

public class Player_Single_Rule extends mutation.sim.Player {

	private Random random;
	private Console console;
	private int numMutations;
	private int m;
	private List<Map<Character, Integer>> patternTracker = new ArrayList<>();
	private List<Map<Character, Integer>> actionTracker = new ArrayList<>();
	private static final int DISTANCE_THRESHOLD_SB = 12;
	private static final double PERCENTAGE_DISTANCE_THRESHOLD_MB = 0.15; 
	private static final int NONEXISTENT_BASE_THRESHOLD = 20;
	private static final double PERCENTAGE_MUTATION_SIZE_DELTA_THRESHOLD = 0.25;
	private List<Integer> numMatchesForMutationSizes = new ArrayList<>();
	private List<Mutagen> possibleMutagens = new ArrayList<>();
	private List<String> genomes = new ArrayList<>();
	private List<String> mutatedGenomes = new ArrayList<>();
	private List<Integer> numMutationsList = new ArrayList<>();
	private List<Mutagen> mutagensGuessed = new ArrayList<>();
	private List<IntervalArchive> intervalArchives = new ArrayList<>();

	public Player_Single_Rule() { 
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
		this.m = m;
		this.console = console;
		if(mutagensGuessed.size() == 0) {
			setUpStructures();
			int numBasesMutatedPerMutation = 0;
			int totalNumMutations = 0;
			for(int i = 0; i < 100; i++) {
//				String genome = "tgccacggtctgttttgcatgcacctagttcgttttctggggcgttagacagactatgttgcttcaccccgaaccccaaatttatacatttctgcagttcccgtttgaatgtgctaatcggtaacatctgccttgtgacgcagaaggtatcaagctggaatgccacacggacactcgagcgctgatctcccgggttggatgacagcagtagaaagttaagggattgcgatggggcagcaaagtgtgagcaatggttagattcgatttcgcgtcccatttgaatgcacggtgaacgttttcacttagcaccgttatgggtgggaagtgtaggattttggatgggctgcatagcagtcgacgtcaccgtaaaatggaccgccgctatatatcggagatttaacgagccctagagaatggggattataacctatgtcagtcatccatattagacgactcgacgcgagaacctgtgattataggaatccgaacgagttcaaatccaaaaggcggtccgctcaaggccgcctcggctcccagactgtctaaaacgcctgtccctacagtgtcttattatcgcgacgtcattagggatgggaaggtgtacaggcgaaattaccgtggtacacaataaatcataatcctaggtgagcagcgcattctatgtgagtaactaggtgtctaaggtgacggtattgctacggatagggttggtgcggacgatgtgagctagttagtacgcagattgcgaacatttcccggcttacacggctcgagtcgtctggccggtagggaattcttatgtggataaagccgacagtatgcagaaaggttgcattagaaataattggacgcgggttcggtatcgcctcggcgtaacaagagatttatgataatcgtgggtacaaaagcaggtgtcccgggtctgtattgcatggttattcagtcaccaagggctaatatgaatggttgacacgagcacggaaataccagcacagtttca";
				String genome = randomString();
				String mutatedGenome = console.Mutate(genome);
				genomes.add(genome);
				mutatedGenomes.add(mutatedGenome);

				this.numMutations = console.getNumberOfMutations();
				numMutationsList.add(numMutations);
				totalNumMutations += numMutations;

				int numBasesMutated = 0;
				for(int j = 0; j < genome.length(); j++)
					if(genome.charAt(j) != mutatedGenome.charAt(j))
						numBasesMutated++;

				int iterNumBasesMutatedPerMutation = (int) Math.ceil(numBasesMutated * 1.0 / numMutations);
				computeWindowSizes(genome, mutatedGenome, iterNumBasesMutatedPerMutation);

				try {
					Process process = Runtime.getRuntime().exec("python3 mutation/g5/Mutagen.py " + genome + " " + mutatedGenome + " " + numMutations);
					process.waitFor();
					InputStream is = process.getInputStream();
					String line = "";
					BufferedReader br = new BufferedReader(new InputStreamReader(is));
					while((line = br.readLine()) != null) {
//						System.out.println("Interval: " + line);

						IntervalArchive intervalObj = new IntervalArchive();
						intervalObj.genome = genome;
						intervalObj.mutatedGenome = mutatedGenome;					
						String intervalString = line.substring(1, line.length() - 1);
						String[] intervalElements = intervalString.split(", ");
						for(int j = 0; j < intervalElements.length; j += 2) {
							int x = Integer.parseInt(intervalElements[j].substring(1));
							int y = Integer.parseInt(intervalElements[j + 1].substring(0, intervalElements[j + 1].length() - 1));
							intervalObj.hcIntervals.add(new Point(x, y));
						}
						intervalArchives.add(intervalObj);
					}
					br.close();
					is.close();
					process.destroy();
				} catch (IOException | InterruptedException e) {
//					e.printStackTrace();
				}

				numBasesMutatedPerMutation = (int) Math.max(numBasesMutatedPerMutation, iterNumBasesMutatedPerMutation);        	
			}

			List<Mutagen> approach2Mutagens = new ArrayList<>();
			try {			
				String genomesString = genomes.toString();
				genomesString = genomesString.replace(" ", "");
				String mutatedGenomesString = mutatedGenomes.toString();
				mutatedGenomesString = mutatedGenomesString.replace(" ", "");
				String numMutationsListString = numMutationsList.toString();
				numMutationsListString = numMutationsListString.replace(" ", "");
				Process process2 = Runtime.getRuntime().exec("python3 mutation/g5/newapproach.py " + genomesString + " " + mutatedGenomesString + " " + numMutationsListString);
				process2.waitFor();
				InputStream is2 = process2.getInputStream();
				BufferedReader br2 = new BufferedReader(new InputStreamReader(is2));
//				String mutagenSetString = br2.readLine();
//				System.out.println(mutagenSetString);
//				if(mutagenSetString.startsWith(prefix))
				String l;
//				while((l = br2.readLine()) != null)
//					System.out.println(l);
				String line = ((l = br2.readLine()) != null) ? l : "";
				System.out.println("Line, part 1: " + line);
				if(line != "") {
					System.out.println("Line, part 2: " + line);
					if(line.substring(0, 1).equals("1")) {
						System.out.println("Line, part 3: " + line);
						line = line.substring(4, line.length() - 2);
						System.out.println("Line, part 4: " + line);
						Mutagen mutagen = new Mutagen();
						mutagen.add(line.split("@")[0], line.split("@")[1]);
						this.possibleMutagens.add(mutagen);
						if(console.Guess(mutagen))
							return mutagen;
					}
				}
				
				br2.close();
				is2.close();
				process2.destroy();
			} catch (Exception e) {
				e.printStackTrace();
			}

			int mostLikelyMutationSize = numBasesMutatedPerMutation;
			if(mostLikelyMutationSize != 1) {
				System.out.println("Total number of mutations: " + totalNumMutations);
				int prevMutationSize = 0;
				int currMutationSize = numMatchesForMutationSizes.get(0);
				mostLikelyMutationSize = 1;
				System.out.println("Number of matches of size 1: " + numMatchesForMutationSizes.get(0));
				for(int i = 1; i < numMatchesForMutationSizes.size(); i++) {
					prevMutationSize = currMutationSize;
					currMutationSize = numMatchesForMutationSizes.get(i);
					double percentDifference;
					if(prevMutationSize == 0 && currMutationSize != 0)
						percentDifference = Double.MAX_VALUE;
					else if(prevMutationSize == 0 && currMutationSize == 0)
						percentDifference = 0;
					else
						percentDifference = (currMutationSize - prevMutationSize) * 1.0 / prevMutationSize;

					if(percentDifference < 0)
						break;
					if(percentDifference > PERCENTAGE_MUTATION_SIZE_DELTA_THRESHOLD)
						mostLikelyMutationSize = i + 1;
					System.out.println("Number of matches of size " + (i + 1) + ": " + numMatchesForMutationSizes.get(i));
				}
			}
			System.out.println("Most likely mutation size: " + mostLikelyMutationSize);
			
			List<Integer> mutationSizeOrder = new ArrayList<>();
			mutationSizeOrder.add(mostLikelyMutationSize);
			if(mostLikelyMutationSize > 1)
				mutationSizeOrder.add(mostLikelyMutationSize - 1);
			if(mostLikelyMutationSize < 10)
				mutationSizeOrder.add(mostLikelyMutationSize + 1);
			for(int i = 1; i <= 10; i++)
				if(!mutationSizeOrder.contains(i))
					mutationSizeOrder.add(i);
			
			for(int i = 0; i < 10; i++) {				
				int mutationSizeChosen = mutationSizeOrder.get(i);
				System.out.println("Mutation size: " + mutationSizeChosen);
				setUpTrackers();
				implementTrackers(mutationSizeChosen);

				System.out.println("Action tracker: ");
				for(int j = 0; j < actionTracker.size(); j++)
					System.out.println("  " + j + ". " + actionTracker.get(j));
				System.out.println("Pattern tracker: ");
				for(int j = 0; j < patternTracker.size(); j++)
					System.out.println("  " + j + ". " + patternTracker.get(j));
				System.out.println();

				List<Mutagen> possibleMutagens = getPossibleMutagens(mutationSizeChosen);
				this.possibleMutagens.addAll(possibleMutagens);
				for(Mutagen mutagen : possibleMutagens) {
					if(console.Guess(mutagen))
						return mutagen;
				}
			}
		}
		
		for(Mutagen mutagen : this.possibleMutagens) {
			if(!mutagensGuessed.contains(mutagen)) {
				mutagensGuessed.add(mutagen);
				return mutagen;
			}
		}
		
		Mutagen mutagen = new Mutagen();
		mutagen.add("a;c;c", "att");
		mutagen.add("g;c;c", "gtt");
		return mutagen;
	}

	private void computeWindowSizes(String genome, String mutatedGenome, int numBasesMutatedPerMutation) {
		for(int i = 1; i <= 10; i++) {
			int numMatches = 0;
			List<List<Integer>> matchesFound = new ArrayList<>();
			for(int j = 0; j < genome.length(); j++) {
				int numBasesMutated = 0;
				List<Integer> possibleMatch = new ArrayList<>();
				for(int k = 0; k < i; k++)
					if(genome.charAt((j + k) % genome.length()) != mutatedGenome.charAt((j + k) % mutatedGenome.length())) {
						possibleMatch.add((j + k) % genome.length());
						numBasesMutated++;
					}
				if(numBasesMutated == numBasesMutatedPerMutation) {
					boolean matchExistsAlready = false;
					for(int k = 0; k < matchesFound.size(); k++) {
						if(matchesFound.get(k).equals(possibleMatch)) {
							matchExistsAlready = true;
							break;
						}
					}
					if(!matchExistsAlready) {
						matchesFound.add(possibleMatch);
						numMatches++;
					}
				}
			}
			numMatchesForMutationSizes.set(i - 1, numMatchesForMutationSizes.get(i - 1) + numMatches);
		}
	}

	private List<Mutagen> getPossibleMutagens(int numBasesMutatedPerMutation) {
		if(numBasesMutatedPerMutation == 1)
			return getOneBaseMutationMutagens();
		return getMultipleBaseMutationMutagens(numBasesMutatedPerMutation);
	}

	private List<Mutagen> getOneBaseMutationMutagens() {
		List<Mutagen> mutagens = new ArrayList<>();
		Map<Integer, List<String>> interestingPatternsMap = new TreeMap<>();
		for(int i = 0; i < patternTracker.size(); i++) {
			Map<Character, Integer> patternOccurrences = patternTracker.get(i);
			List<String> interestingPattern = new ArrayList<>();
			for(Character character : patternOccurrences.keySet()) {
				if(patternOccurrences.get(character) > (NONEXISTENT_BASE_THRESHOLD * m * 1.0 / 15))
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
				pattern += "actg";
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
						pattern += "actg;";
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
						newPattern = "actg;" + newPattern;

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
			pattern += "actg";
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
					pattern += "actg;";
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
						intermediatePattern = "actg;" + intermediatePattern;
			
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
						newPattern = "actg;" + newPattern;
				
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
				for(Character character : patternOccurrences.keySet()) {
					if(patternOccurrences.get(character) > (NONEXISTENT_BASE_THRESHOLD * m * 1.0 / 15))
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
			pattern += "actg";
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
					pattern += "actg;";
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
							intermediatePattern = "actg;" + intermediatePattern;

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
				System.out.println("Genome value: " + genomeValue);
				System.out.println("Mutated genome value: " + mutatedGenomeValue);
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