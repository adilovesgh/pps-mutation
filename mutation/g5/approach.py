import os
import sys
import glob
import math
import random
import pickle
from typing import *

def intersect(action1, action2):
    return [a1.intersection(a2) for a1, a2 in zip(action1, action2)]

def parse_offset(elements, offset):
    offset = abs(offset)
    return set(e if isinstance(e, str) else e-offset for e in elements if isinstance(e, str) or (e>=offset and e<offset+10))

def action_potential(action):
    prod = 1
    for ch in action: prod *= len(ch)
    return prod

def find_intervals(genome, result, mp):
    if mp == 0: return []

    diffs = [idx for idx, (c1, c2) in enumerate(zip(genome, result)) if c1 != c2]

    intervals = []

    startidx = None
    for d in diffs:
        if startidx is None:
            startidx, curidx = d, d

        elif d - curidx < 10:
            curidx = d

        else:
            intervals.append((startidx, curidx))
            startidx, curidx = d, d

    intervals.append((startidx, curidx))
    return intervals

def hack_action(enumidxes, ch):
    if isinstance(ch, str) and len(ch)==1:
        return set(enumidxes[ch])

    tmp = set()
    for c in ch:
        for t in enumidxes[c]:
            tmp.add(t)
    return tmp

def get_action(before, after):
    enumidxes = {}
    for key in 'acgt':
        enumidxes[key] = [key] + [i for i, ch in enumerate(before) if (isinstance(ch, str) and ch==key) or (isinstance(ch, set) and key in ch)]
#         key: [key] + [i for i, ch in enumerate(before) if ch==key] for key in 'actg'}

    return [hack_action(enumidxes, ch) for ch in after]

def get_offset(idx1, idx2):
    offset = None

    # idx1, idx2 should already be sorted by (exp, idx) format
    for (e1, i1), (e2, i2) in zip(idx1, idx2):
        if offset is None and e1==e2:
            offset = i2 - i1
        elif e1==e2:
            if i2-i1 != offset:
                return None

    return offset

def reduce_hypothesis(actionlist):
    """ Isomorphisms are hard to detect """
    newactionlist = []
    for actions1 in actionlist:
        pattern1, action1, idxes1 = actions1

        for nalidx, actions2 in enumerate(newactionlist):
            pattern2, action2, idxes2 = actions2

            # No offset check
            if (p for _, p in pattern1) == (p for _, p in pattern2) and action1==action2:
                break

            # Pattern 1 offset check
            flag = None
            for offset in range(9):
                if flag is not None: break

                if pattern1[offset][1] != 'acgt' or pattern2[9-offset][1] != 'acgt': break
                if action1[offset] != {offset} or action2[9-offset] != {9-offset}:   break

                for i in range(9-offset):
                    idx = offset + i + 1

                    if pattern2[i][1] != pattern1[idx][1]: break
                    if action2[i]  != parse_offset(action1[idx], offset+1): break
                else:
                    flag = (1, offset+1)

            # Pattern 2 offset check
            for offset in range(9):
                if flag is not None: break

                if pattern2[offset][1] != 'acgt' or pattern1[9-offset][1] != 'acgt': break
                if action2[offset] != {offset} or action1[9-offset] != {9-offset}:   break

                for i in range(9-offset):
                    idx = offset + i + 1

                    if pattern1[i][1] != pattern2[idx][1]: break
                    if action1[i]  != parse_offset(action2[idx], offset+1): break

                else:
                    flag = (2, offset+1)

            if flag is not None:
                mode, offset = flag
                if mode == 2:
#                     newactionlist[nalidx] = actions1
                    tmp_pattern = pattern2[offset:] + get_pattern([{}]*offset)
                    new_pattern = union_pattern(pattern1, tmp_pattern)
                    newactionlist[nalidx] = [new_pattern, action1, idxes1]
                else:
                    tmp_pattern = pattern1[offset:] + get_pattern([{}]*offset)
                    new_pattern = union_pattern(pattern2, tmp_pattern)
                    newactionlist[nalidx] = [new_pattern, action2, idxes2]
                break
        else:
            newactionlist.append(actions1)

    return newactionlist

mask = {'a':0, 'c':1, 'g':2, 't':3}

def get_pattern(before):
    pattern = []
    for ch in before:
        counts = [0., 0., 0., 0.]
        if isinstance(ch, str):
            counts[mask[ch]] += 1.
        else:
            l = len(ch)
            if l != 4:
                for c in ch:
                    counts[mask[c]] += 1./(l+(l!=1))
        pattern.append((counts, 'acgt'))
    return pattern

def parse_pattern(pattern):
     return [parse_count(count) for count, _ in pattern]

def parse_count(count, print_exp=False):
    # Use a simple 5% threshold
    total_count = sum(count)+1
    exp = [(e+1)/total_count for e in count]

    return ''.join("acgt"[i] for i, e in enumerate(exp) if e > 0.03)

def union_pattern(pattern1, pattern2):
    pattern = []
    for (p1, c1), (p2, c2) in zip(pattern1, pattern2):
        counts = [p1[i]+p2[i] for i in range(4)]
        string = parse_count(counts)
        pattern.append((counts, string))
    return pattern

def produce_guess(actionlist):
    results = set()

    for pattern, action, _ in actionlist:
        new_action = [set(a) for a in action]
        for i in range(len(action)):
            for ch in list(action[i]):
                if not isinstance(ch, int): continue

                if len(pattern[ch][1]) == 1:
                    new_action[i].remove(ch)
                    new_action[i].add(pattern[ch][1])

        if action_potential(new_action) != 1: continue

        # Converting from list of sets to just a list of int/char
        flat_action = [pos.pop() for pos in new_action]

        nums_referenced = [ch for idx, ch in enumerate(flat_action) if idx != ch]

        i = 0
        for (_, patpos), actpos in zip(pattern, flat_action):
            if i in nums_referenced: break
            if not (patpos=="acgt" and actpos==i): break
            i += 1

        pattern_string = ";".join(pat[1] for pat in pattern[i:])
        while pattern_string.endswith(";acgt"):
            pattern_string = pattern_string[:-5]

        while flat_action[-1] == len(flat_action)-1:
            flat_action.pop()

        action_string = "".join(
            f"{tmp-i}" if isinstance(tmp, int) else tmp
            for tmp in flat_action[i:])

        results.add(f"{pattern_string}@{action_string}")
        if len(results) > 10: break

    # print("DEBUG : ", len(pattern_strings), len(action_strings))
    return results

def new_approach_one(genome, result, exp, start):
    l = len(genome)

    possible_hypothesis = []

    for idx in range(l - 9):
        action = get_action(genome[idx:idx+10], result[idx:idx+10])

        if action_potential(action) > 0:
            possible_hypothesis.append((get_pattern(genome[idx:idx+10]), action, [(exp, start+idx)]))

    return possible_hypothesis

def new_approach_two(genome, result, exp, start):
    l = len(genome)

    alldiffs = set(idx for idx, (c1, c2) in enumerate(zip(genome, result)) if c1 != c2)

    enumidxes = {}
    for key in 'atcg':
        enumidxes[key] = [key] + [i for i, ch in enumerate(genome) if ch==key]

    separate_actions = []
    for idx1 in range(l-19):
        remdiff1 = set(i for i in alldiffs if not idx1<=i<idx1+10)
        if len(remdiff1) > 10 or max(remdiff1)-min(remdiff1) > 10: continue

        action1  = [parse_offset(enumidxes[ch], idx1) for ch in result[idx1:idx1+10]]

        for idx2 in range(idx1+10, l-9):

            # If the two windows don't cover all differences, it cannot be possible
            remdiff2 = set(i for i in remdiff1 if not idx2<=i<idx2+10)
            if len(remdiff2) > 0: continue

            action2  = [parse_offset(enumidxes[ch], idx2) for ch in result[idx2:idx2+10]]
            action   = intersect(action1, action2)

            # Check if the intersection is possible
            if set() in action: continue

            # Add to possible actions
            separate_actions.append((
                union_pattern(get_pattern(genome[idx1:idx1+10]), get_pattern(genome[idx2:idx2+10])),
                action, [(exp, start+idx1), (exp, start+idx2)]))

    combined_actions = []
    for idx1 in range(l-9):
        remdiff1 = set(i for i in alldiffs if not idx1<=i<idx1+10)
        if len(remdiff1) > 10 or (remdiff1 and max(remdiff1)-min(remdiff1) > 10): continue

        for idx2 in range(max(0, idx1-9), min(l-9, idx1+10)):
            # If the two windows don't cover all differences, it cannot be possible
            remdiff2 = set(i for i in remdiff1 if not idx2<=i<idx2+10)
            if len(remdiff2) > 0: continue

            minidx, maxidx = min(idx1, idx2), max(idx1, idx2)
            intersection_idxes = set(range(maxidx, minidx+10))

            idx1off = idx1 - minidx
            idx2off = idx2 - minidx

            before = [genome[i] for i in range(minidx, maxidx+10)]
            middle = [set('actg') for i in range(minidx, maxidx+10)]
            after  = [result[i] for i in range(minidx, maxidx+10)]

            for i in range(10):
                if i+idx1 not in intersection_idxes:
                    middle[i+idx1-minidx] = after[i+idx1-minidx]

                if i+idx2 not in intersection_idxes:
                    middle[i+idx2-minidx] = before[i+idx2-minidx]

            for i in range(20):
                action1 = get_action(before[idx1-minidx:idx1-minidx+10], middle[idx1-minidx:idx1-minidx+10])
                action2 = get_action(middle[idx2-minidx:idx2-minidx+10], after [idx2-minidx:idx2-minidx+10])

                action   = intersect(action1, action2)

                if action_potential(action) == 0:
                    break

                prev_middle = middle[:]

                # Check overlaps calculated from before
                for i in range(minidx + 10 - maxidx):
                    from_before = set(before[dst+idx1-minidx] for dst in action[i+maxidx-idx1] if isinstance(dst, int))
                    middle[i+maxidx-minidx] = from_before.intersection(middle[i+maxidx-minidx])

                    # Convert set to character if possible
                    if len(middle[i+maxidx-minidx]) == 1:
                        middle[i+maxidx-minidx] = list(middle[i+maxidx-minidx])[0]

                # Check overlaps calculated from after
                for i in range(10):
                    if len(action[i]) != 1: continue

                    src = list(action[i])[0]
                    if isinstance(src, str): continue

                    if after[i + (idx2 - minidx)] not in middle[src+idx2-minidx]: # Handle both str and set()
                        middle[src+idx2-minidx] = set()
                    else:
                        middle[src+idx2-minidx] = after[i + (idx2 - minidx)]

                if prev_middle == middle: break

            if action_potential(action) == 0: continue

            combined_actions.append((
                union_pattern(
                    get_pattern(before[idx1-minidx:idx1+10-minidx]),
                    get_pattern(middle[idx2-minidx:idx2+10-minidx])),
                action, [(exp, start+idx1), (exp, start+idx2)]))

    return separate_actions + combined_actions

def parse_experiment_results(before, after, num_mutations, actionlist, exp, start):
    # print(f"PARSING [{before}=>{after}] {num_mutations}")

    if num_mutations == 1:
        actions = new_approach_one(before, after, exp, start)
    elif num_mutations == 2:
        actions = new_approach_two(before, after, exp, start)

    # print(actions)

    if actionlist is None:
        return reduce_hypothesis(actions)

    possible_actions = []
    for pattern1, a1, idx1 in actionlist:
        for pattern2, a2, idx2 in actions:
            action = intersect(a1, a2)

            # Actions are not possible
            if action_potential(action) == 0: continue

            # Probability of pattern matching is super low
            flag = True
            for (count1, chars1), (count2, chars2) in zip(pattern1, pattern2):
                if num_mutations == 1:
                    for i, k in enumerate("acgt"):
                        if count2[i] == 1:
                            if k not in chars1:
                                # print(f"[{count1} {count2}] {k} was not in {chars1}")
                                flag = False
                                # print(num_mutations, (count1,chars1), (count2, chars2), (i,k))
                                break

                if num_mutations == 2:
                    for i, k in enumerate("acgt"):
                        if count2[i] >= 1:
                            if k not in chars1:
                                flag = False
                                # print(num_mutations, (count1,chars1), (count2, chars2), (i,k))

                if not flag:

                    break

            if not flag: continue
            possible_actions.append((
                union_pattern(pattern1, pattern2), action, idx1+idx2))

    return reduce_hypothesis(possible_actions)

def recursive(genome, result, actionlist, num_left, num_mutations, experiment, intervals):
    if not intervals:
        return actionlist

    # print()
    # print(f"({num_left}, {num_mutations}) {intervals}")
    # if actionlist: print_actionlists([actionlist], "Input actionlist")

    cap_mutations = [(num_mutation, num_mutation + num_left) for num_mutation in num_mutations]
    (minM, maxM), (d1, d2), idx = min(zip(cap_mutations, intervals, range(len(intervals))))

    newintervals = [d for d in intervals if d != (d1, d2)]

    onelist, twolist = [], []
    if minM <= 1 <= maxM:
        before, after = genome[d2-9 : d1+10], result[d2-9 : d1+10]
        one_actionlist = parse_experiment_results(before, after, 1, actionlist, experiment, d2-9)

#         print_actionlists([one_actionlist], f"Calling at {len(intervals)} with [{before}] -> [{after}]")
#         if len(intervals) == 2:
#             print(f"{before}, {after}, 1, {actionlist}, {experiment}, {d2-9}")

        onelist = recursive(genome, result, one_actionlist, num_left,
                [n for i, n in enumerate(num_mutations) if i!= idx], experiment, newintervals)

    if minM <= 2 <= maxM:
        before, after = genome[d1-9 : d2+9], result[d1-9 : d2+9]
        two_actionlist = parse_experiment_results(before, after, 2, actionlist, experiment, d1-9)

        twolist = recursive(genome, result, two_actionlist, num_left - (2-minM),
                [n for i, n in enumerate(num_mutations) if i!= idx], experiment, newintervals)

#     print_actionlists([onelist, twolist], f"RECURSIVE DEBUG {len(intervals)}")
    return onelist + twolist

if __name__ == "__main__":
    genome = sys.argv[1]
    result = sys.argv[2]
    mp = int(sys.argv[3])

    experiment = int(sys.argv[4])
    num_m = int(sys.argv[5])
    prev_error = int(sys.argv[6])

    # At experiment 0 or error, reset everything
    if experiment == 0 or prev_error == 1:
        with open("prev_guesses", "wb") as handle:
            pickle.dump(set(), handle)
        with open("prev_actionlist", "wb") as handle:
            pickle.dump(None, handle)
        with open("memory", "wb") as handle:
            pickle.dump([], handle)

    # Load previous states
    with open("prev_guesses", "rb") as handle:
        prev_guesses = pickle.load(handle)
    with open("prev_actionlist", "rb") as handle:
        actionlist = pickle.load(handle)
    with open("memory", "rb") as handle:
        prev_memory = pickle.load(handle)

    # Rotates if necessary
    diffs = [idx for idx, (c1, c2) in enumerate(zip(genome, result)) if c1 != c2]

    if diffs[-1] > 990: # Possibility of a rollover
        for i in range(len(diffs)-1):
            if diffs[i+1] - diffs[i] > 20:
                genome = genome[diffs[i]+10:] + genome[:diffs[i]+10]
                result = result[diffs[i]+10:] + result[:diffs[i]+10]
                break

    intervals = find_intervals(genome, result, mp)
    num_mutations = [(d2-d1+26)//18 for (d1, d2) in intervals]
    num_left = mp - sum(num_mutations)

    if max(num_mutations) + num_left <= (1 if num_m < 1 else 2):
        actionlist = recursive(genome, result, actionlist, num_left, num_mutations, experiment, intervals)

        if len(actionlist) == 0 or len(actionlist) > 1000:
            actionlist = None

        if actionlist is not None:
            guess = produce_guess(actionlist)

            guess = guess.difference(prev_guesses)
            if len(guess) < 5:
                for g in guess:
                    prev_guesses.add(g)
                    print(g)

    with open("prev_guesses", "wb") as handle:
        pickle.dump(prev_guesses, handle)
    with open("prev_actionlist", "wb") as handle:
        pickle.dump(actionlist, handle)
    with open("memory", "wb") as handle:
        pickle.dump(prev_memory + [(genome, result)], handle)
