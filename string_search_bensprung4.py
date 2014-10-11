#!/usr/bin/python
# Fall 2014 
# Intro to Computational Biology
# HW 2
# Ben Sprung
# Looks for a pattern p against a string s using (1) a naive search, and
# (2) the bad character rule, (3) the good suffix rule
# One string s is read from the first line of a file called "database_ben.txt"
# One or more patterns p are read from a file called "pattern_ben.txt". Each 
# line is used as a separate p

import sys
from time import sleep

# animation of search
animation_speed = 0.15

# open files, read the database string s and search patterns
with open("database_ben.txt") as f:
    s = f.readline().rstrip()
patterns = []
with open("pattern_ben.txt") as g:
    for line in g:
        patterns.append(line.rstrip())   

alphabet = ["A","C","G","T"]
#alphabet = [chr(i) for i in range(32,127)]   #most standard characters

def naive_search(pattern, string):
    p, s = pattern, string
    total_shifts = 0
    for i in range(len(s) - len(p) + 1) :  #need +1 for last position
        #animation:
        sys.stdout.write(i * " " + p + "\r")
        sleep(animation_speed); sys.stdout.flush()

        perfect_match = True
        for j in range(len(p)):
            if p[j] != s[i+j]: 
                perfect_match = False
                total_shifts += 1
                break        #no need to check more after first mismatch
        if perfect_match:
            return i, total_shifts        #we get first match only
    return -1, total_shifts

def bad_character_search(pattern,string):
    """From a given position search right to left. If we hit a mismatched
    character, look up the shift in the precomputed shifts dictionary for the
    bad character rule, which is indexed by the character in the database that
    does not match, and then the position of the mismatch"""
    p, s = pattern, string
    bad_char_shifts_dict = compute_bad_character_shifts(pattern)
    i = 0
    total_shifts = 0
    while i < (len(s) - len(p) + 1):
        #animation:
        sys.stdout.write(i * " " + p + "\r")
        sleep(animation_speed); sys.stdout.flush()
        
        perfect_match = True
        for j in reversed(range(len(p))):
            if p[j] != s[i + j]:
                perfect_match = False
                i += bad_char_shifts_dict[s[i+j]][j]
                total_shifts += 1
                break
        if perfect_match:
            return i, total_shifts
    return -1, total_shifts

def good_suffix_search(pattern,string):
    """From a given position search right to left. If we hit a mismatched
    character, look up the shift in the precomputed shifts dictionary for the
    good suffix rule, which is indexed by the position of the pattern that
    is causing the mismatch"""
    p, s = pattern, string
    good_suffix_shifts_dict = compute_good_suffix_shifts(p)
    i = 0
    total_shifts = 0
    while i < (len(s) - len(p) + 1):
        #animation:
        sys.stdout.write(i * " " + p + "\r")
        sleep(animation_speed); sys.stdout.flush()
        
        perfect_match = True
        for j in reversed(range(len(p))):
            if p[j] != s[i + j]:
                perfect_match = False
                i += good_suffix_shifts_dict[j+1] #it's indexed to follow notes
                total_shifts += 1
                break
        if perfect_match:
            return i, total_shifts
    return -1, total_shifts

def compute_bad_character_shifts(pattern):
    """For every position in the pattern, for every character in the alphabet,
    calculate the distance, from right to left, to the next occurence of that 
    character. If the current position is the character in question, store a 0.
    If the character never occurs, store the number of positions from the
    character to the left end of the pattern."""
    p = pattern
    # initialize dictionary in which results will be stored with all zeros
    bad_char_shifts_dict = {}
    for char in alphabet:
        bad_char_shifts_dict[char] = [0 for i in range(len(p))]

    for i in range(len(p)):
        for char in alphabet:
            # look left from present position until hit first character of 
            # pattern or until hit alphabet character
            j = 0
            while True:
                if p[i-j] == char:
                    bad_char_shifts_dict[char][i] = j
                    break
                if j == i:
                    bad_char_shifts_dict[char][i] = i + 1
                    # plus one b/c if we are in leftmost position we still
                    # want to shift forward one (not zero)
                    break
                j += 1
    #print bad_char_shifts_dict
    return bad_char_shifts_dict

def compute_good_suffix_shifts(pattern):
    """For each substring of a pattern, find the position in which it matches
    the full pattern, starting with right alignment and sliding the substring 
    left by one iteratively. The leftmost character of the substring is treated
    as not-that-character. The max possible offset is the length of the 
    full string minus one. Store the results in a dictionary where the keys are
    the positions of the pattern to database mismatch that imply the substrings
    and the values are the offsets."""
    good_substring_shifts_dict = {}

    # make all the substrings of the given pattern
    substrings = []
    for i in range(len(pattern)):
        substrings.append(pattern[i:])
    
    for substring in substrings:
        offset = 0
        while True:
            # these start positions will handle lining up and overhang
            # worked these out on paper hard to recap here
            start_position_pattern = len(pattern) - len(substring) - offset
            start_position_substring = start_position_pattern
            if start_position_substring < 0: start_position_substring *= -1
            else: start_position_substring = 0
            if start_position_pattern < 0: start_position_pattern = 0

            match = True
            # first, check the first (not-) character which is NOT allowed to 
            # equal the corresponding character on the pattern
            # Note that we don't always check the first position of the 
            # substring: when it is shifted leftwards enough, that position
            # is out of the picture
            if start_position_substring == 0:
                if substring[0] == pattern[start_position_pattern]:
                    match = False
                # since we checked, advance both patterns' start position
                start_position_pattern += 1
                start_position_substring += 1
            # now check the remainder of the substring
            # against the corresponding slice of the pattern
            if substring[start_position_substring:] !=\
               pattern[start_position_pattern:\
               start_position_pattern+len(substring)-start_position_substring]:
                match = False
            if match:
                good_substring_shifts_dict[len(pattern)-len(substring)+1] = \
                                                                        offset
                break
            # if we don't find a match, stop when rightmost substring character
            # is aligned with leftmost pattern character
            if offset == len(pattern) - 1:
                good_substring_shifts_dict[len(pattern)-len(substring)+1] = \
                                                                    offset + 1
                break
            offset += 1
    return good_substring_shifts_dict

if __name__ == "__main__":
    #print compute_good_suffix_shifts("GCGACG")
    for p in patterns:
        print "---\n" + "naive search:\n" + s
        result, shifts = naive_search(p,s)
        print "\nMatch at position:", result, "Total shifts:", shifts

        print "---\n" + "bad character search:\n" + s
        result, shifts = bad_character_search(p,s)
        print "\nMatch at position:", result, "Total shifts:", shifts

        print "---\n" + "good suffix search:\n" + s
        result, shifts = good_suffix_search(p,s)
        print "\nMatch at position:", result, "Total shifts:", shifts
