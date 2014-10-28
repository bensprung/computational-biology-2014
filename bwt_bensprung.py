#!/usr/bin/python
# Fall 2014 
# Intro to Computational Biology
# HW 5
# Ben Sprung
# Assignment: use the BWT string "AG$ATGGACAGCAGCTA" to creage the function
# N[X] and the table B[X,i]. Find the position of the substring GAG in the
# original string if it exists. Otherwise state that it doesn't exist. Use the
# recursive searching scheme shown in class. List all the values of positions
# you find (i.e., minimum and maximum row positions for "G", then "AG", then
# "GAG". Repeat for "GCTT".

def define_characters(bwt_string):
    """  Compute a sorted (alphabetized) list of the characters that make up
    the BWT string, except for "$". (It is not presumed that the
    characters are always going to be A,C,G,T).

    Returns:
    The list of characters as a string, because I need to use "find" on it
    later
    """

    chars = list(set(bwt_string.replace("$","")))
    chars.sort()
    return "".join(c for c in chars)

def compute_N(db_string):
    """ Precompute the function N[X] for Burrows-Wheeler transforms. 
    
    Arguments:
    db_string: the search string or its BWT transform

    Globals:
    chars: sorted list of all chars (except "$") of the string
   
    Returns:
    For each character in the database string, compute the number of symbols in
    the database string that is lexicographically less than the character.
    Returned as a dictionary.

    Examples:
    n("ACAACGT", "A") = 0
    n("ACAACGT", "C") = 3
    n("ACAACGT", "G") = 5
    """

    # strip out the "$" in the BWT if it exists
    db_string = db_string.replace("$","")

    i = 0
    output = {}
    for char in chars:
        output[char] = i
        i += db_string.count(char)

    return output

def compute_B(bwt_string):
    """ Precompute the function B[x,i], for all i, for a BWT transform
    For each character x, compute the number of x's found at or before position
    i in the BWT string.

    Arguments:
    bwt_string: a Blake-Wheeler transform string

    Globals:
    chars: sorted list of all chars (except "$") of the string

    Returns:
    A matrix, implemented as a list of lists, where the rows are the characters,
    and the columns are the positions, and the values are as described above.
    The rows are ordered alphabetically.
    """

    # Initialize a matrix with as many columns as the length of the BWT string
    # and as many rows as unique characters
    output_matrix = [[0 for i in range(len(bwt_string))] for j in \
                                                 range(len(chars))]
    row = 0
    for char in chars:
        i = 0
        number_found = 0
        while i < len(bwt_string):
            if bwt_string[i] == char:
                number_found += 1
            output_matrix[row][i] = number_found
            i += 1
        row += 1

    return output_matrix

def min_bwt(bwt_string, substring):
    """ The min function for BWT substring finding.
    It is defined recursively: min(xW) = N(x) + B[x, min(W)-1] + 1
    In this notation, W is a substring we are looking for and x is a single
    letter that extends W as a prefix

    The boundary condition when W is a single letter is min(W) = N(W) + 1

    Arguments:
    bwt_string: The BWT string
    substring: the substring we are searching for

    Assumed to exist as globals:
    N_dict: a dictionary representing the output of N[x]
    B_matrix: a list of lists representing the output of B[x,i]

    Returns:
    A number that is the minimum row index of the sorted circular permutation
    table.
    """

    if len(substring) == 1:
        return N_dict[substring] + 1
    else:
        x = substring[0]
        W = substring[1:]
        x_row = chars.find(x)
        return N_dict[x] + B_matrix[x_row][min_bwt(bwt_string,W)-1] + 1


def max_bwt(bwt_string, substring):
    """ The max function for BWT substring finding.
    It is defined recursively: max(xW) = N(x) + B[x, max(W)]

    The boundary condition when W is a single letter is max(W) = N(Y)
    where Y is the symbol that lexicographically follows X

    Arguments:
    bwt_string: The BWT string
    substring: the substring we are searching for

    Assumed to exist as globals:
    N_dict: a dictionary representing the output of N[x]
    B_matrix: a list of lists representing the output of B[x,i]

    Returns:
    A number that is the maximum row index of the sorted circular permutation
    table.
    """

    if len(substring) == 1: #boundary condition
        try: 
            next_character = chars[chars.find(substring) + 1]
        except IndexError:
            next_character = chars[0]
        return N_dict[next_character] 
    else:
        x = substring[0]
        W = substring[1:]
        x_row = chars.find(x)
        return N_dict[x] + B_matrix[x_row][max_bwt(bwt_string,W)]

def recover_circular_permutation_table(bwt_string):
    """ Need to recover the CPT to get the numbers that will tell us the
    starting positions of the substrings once we have the min and max. 

    Arguments:
    bwt_string: the BWT string

    Returns: the CPT as a list rows, which are just strings. 
    """
    # Initialize table as list of strings with as many rows as BWT string has
    # characters
    cp_table = []
    for row in bwt_string:
        cp_table.append("")

    # Reconstruct the table
    for n in range(len(bwt_string)):
        i = 0
        for char in bwt_string:
            cp_table[i] = char + cp_table[i] # prepend
            i += 1
        cp_table.sort()
    # print cp_table
    return cp_table

def bwt_find(bwt, substring):
    """A place to execute all the other functions"""

    print "BWT string is %s" % bwt

    global chars; global N_dict; global B_matrix
    chars = define_characters(bwt)
    N_dict = compute_N(bwt)
    B_matrix = compute_B(bwt)
    min_find = min_bwt(bwt, substring)
    max_find = max_bwt(bwt, substring)
    
    # Recover the CPT so we can conver the min and max to positions in the
    # original string
    cpt = recover_circular_permutation_table(bwt)
    # May as well get the original string itself for checking of results. The
    # row that starts with 0 is the original string with "$" at the end
    for row in cpt:
        if row[-1] == "$":
            original_db_string = row.strip("$")
    print "Original database string is", original_db_string

    # Turn the min and max into positions
    # The distance of the "$" from the right end of the row is the starting
    # position on the original string
    found_positions = []
    for i in range(min_find, max_find + 1): # range(x, x+1) is [x]
        found_positions.append(cpt[i][::-1].find("$"))   # [::-1] reverses
        found_positions.sort()
    if found_positions:
        print "Substring %s found at position(s): %s" % \
                            (substring, str(found_positions).strip("[]"))
    else:
        print "Substring %s not found" % substring
    print

        
if __name__ == "__main__":
    #bwt_find("TC$AAACG", "A")
    #bwt_find("TC$AAACG", "AC")
    #bwt_find("TC$AAACG", "ACAA")
    bwt_find("AG$ATGGACAGCAGCTA", "GAG")
    bwt_find("AG$ATGGACAGCAGCTA", "GCTT")
    bwt_find("AG$ATGGACAGCAGCTA", "CG")

