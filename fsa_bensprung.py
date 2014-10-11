#!/usr/bin/python
# Fall 2014 
# Intro to Computational Biology
# HW 3
# Ben Sprung
# Assignment: Make a finite state automaton that parses the SOX2 binding sites
# CATTGTTG and CAGGCTTG (option 2)

import sys

# For highlighting, cribbed from stackexchange
BOLD, END = '\033[1m','\033[0m'

# Open file and read the database string
with open("fsa_database_ben.txt") as f:
    # Grab the first line only
    database = f.readline().rstrip()

"""
Transition table for CATTGTTG and CAGGCTTG (SOX2 sites):

   | A  | C  | G  | T  |
q0 | q0 | q1 | q0 | q0 | 
q1 | q2 | q0 | q0 | q0 | 
q2 | q0 | q0 | q4 | q3 | 
q3 | q0 | q0 | q0 | q5 | 
q4 | q0 | q0 | q6 | q0 | 
q5 | q0 | q0 | q7 | q0 | 
q6 | q0 | q7 | q0 | q0 | 
q7 | q0 | q0 | q0 | q8 | 
q8 | q0 | q0 | q0 | q9 | 
q9 | q0 | q0 | q10| q0 | 
"""

# Represent transition table as list of lists
# e.g. list[2] is the transitions for q2.
sox2_transition_table = [[0,1,0,0],
                         [2,0,0,0],
                         [0,0,4,3],
                         [0,0,0,5],
                         [0,0,6,0],
                         [0,0,7,0],
                         [0,7,0,0],
                         [0,0,0,8],
                         [0,0,0,9],
                         [0,0,10,0]]

# We will want to be able to convert from characters 
# to transition table column index.
character_lookup_dict = {"A":0, "C":1, "G":2, "T":3}

def fsa_search(database, transition_table):
    db, tt = database, transition_table
    state = 0   # i.e. q0  
    p = 0       # Counter for where we are in the database
    matches = []
    while True:
        try:
            c = character_lookup_dict[db[p]] # Turn A into 0, C into 1, etc.
        except IndexError:
            return matches      # Not sure if bad practice to detect reaching 
                                # the end of the string this way.
        state = tt[state][c]

        if state == 10:
            state = 0           # Reset the state.
            matches.append(p-7) # Subtract length of binding site minus 1
                                # to get start position.
        p += 1

if __name__ == "__main__":
    matches = fsa_search(database, sox2_transition_table)
    # Print results
    if matches == []:
        print "Did not find CAGGCTTG or CAGGCTTG in:"
        print database
    else:
        print "Found a SOX2 binding site at position(s): " + \
               str(matches).strip("[]")
        # Print database string with highlighted results
        i = 0
        for match in matches:
            sys.stdout.write(database[i:match]+BOLD+
                             database[match:match+8]+END)
            i = match + 8
        sys.stdout.write(database[i:]+"\n")
