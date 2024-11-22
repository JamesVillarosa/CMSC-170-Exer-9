# Villarosa, James Carl V.
# GH4L

import os                                                                               # Import we need in program

def read_input(input_filename):                                                         # Function for reading the input file
    with open(input_filename, 'r') as file:                                             # Opening input file
        lines = file.readlines()                                                        # Reading file per line

    number_of_string = int(lines[0].strip())                                            # Contain the first line to number of DNA sequence
    
    DNA_sequences = []                                                                  # Initialize array to hold DNA sequences
    for i in range(1, number_of_string + 1):                                            # Read the second line up to number of sequences
        DNA_sequences.append(lines[i].strip())                                          # Contain only the sequence no space or newline

    states = lines[number_of_string + 1].strip().split()                                # Split the line into list of states (ex. A,C,G,T)
    observables = lines[number_of_string + 2].strip().split()                           # Split the line into list of observables (ex. H,L)


    observation_probability = {}                                                        # Initialize dictionary to store observation probability for each state
    for i in range(len(states)):                                                        # We will store observation probabilities for each state

        probs = lines[number_of_string + 3 + i].strip().split()                         # Get the probability for the current state
        probs = [float(p) for p in probs]                                               # Convert probability from string to float

        observation_probability[states[i]] = {}                                         # Create dictionary for current state

        for j in range(len(observables)):                                               # Assign the probability to state for each observable
            observation_probability[states[i]][observables[j]] = probs[j]               # Store probability for each observable

    number_of_cases = int(lines[number_of_string + 3 + len(states)].strip())

    cases = []                                                                          # Initialize list to store cases (ex. G1 given H1)
    for i in range(number_of_cases):                                                    # Loop through number of cases
        cases.append(lines[number_of_string + 4 + len(states) + i].strip())             # Add case to list

    return DNA_sequences, states, observables, observation_probability, cases           # Return the datas we gathered in input file


def initial_probabilities(sequence):                                                    # Function for computing initial probability per sequence
    count = {'A': 0, 'C': 0, 'G': 0, 'T': 0}                                            # Initialize dictionary to store the count of A,C,G,T
    count[sequence[0]] += 1                                                             # We will look at the first nucleotide(A) in sequence to start

    total = 0                                                                           # Initialize total
    for key in count:
        total += count[key]                                                             # Add count with index key in total

    initial_probability = {}                                                            # Create dictionary to store initial probability for each nucleotides
    for key in count:                                                                   # Calculate probability for each nucleotide
            initial_probability[key] = count[key] / total                               # Formula is count of key / total, then contain it to initial probability dictionary

    return initial_probability                                                          # Return initial probability dictionary


def transition_probabilities(sequence):                                                 # Function for computing transition probability
    count = {'A': 0, 'C': 0, 'G': 0, 'T': 0}                                            # Initialize again the dictionary for storing the count of nucleotide
    transitions = {}                                                                    # Dictionary to store transitions of nucleotides (ex. AA, AC, CA, etc)

                                                                                        # We consider all possible transition of nucleotides
    for s1 in count:                                                                    # First nucleotide
        for s2 in count:                                                                # Second nucleotide
            transitions[s1 + s2] = 0                                                    # Initialize the transition count for s1 -> s2 to 0

    previous = sequence[0]                                                              # Start with the first nucleotide in sequence as previouse nucleotide
    count[previous] += 1                                                                # Count first nucleotide occurence

    for i in range(1, len(sequence)):                                                   # Loop through seqence starting from second nucleotide up to the end
        current = sequence[i]                                                           # Make current nucleotide the sequence at index i
        count[current] += 1                                                             # Increment count of current nucleotide
        transitions[previous + current] += 1                                            # Increment transiton count for previous -> current nucleotide
        previous = current                                                              # Update previous nucleotide to be current for next transition

    transition_probability = {}                                                         # Now, we will calculate the transition probabilities based on transition counts we get earlier
    for s1 in count:                                                                    # We loop over each possible starting nucleotide
        transition_probability[s1] = {}                                                 # Initialize dictionary for transition from s1
        total_transitions = 0                                                           # Initialize varabile because we will calculate the total number of transition from s1

        for s2 in count:                                                                # Loop through all possible ending nucleotide
            total_transitions += transitions[s1 + s2]                                   # Add  transition from s1 to any s2 to total transition

        for s2 in count:                                                                # We calculate the transition probability
            if total_transitions == 0:                                                  # Check if total transition is 0 to prevent division by 0
                transition_probability[s1][s2] = 0                                      # If 0, we just equate it by 0
            else:
                transition_probability[s1][s2] = transitions[s1 + s2] / total_transitions   # Formula of transition probability

    return transition_probability                                                       # Return dictionary


def observation_probabilities(gene_probability, observation_probability):               # Function that calculates the overall observation probability based on its parameters
                                                                                        # Finding out how likely H and L, considering all posible gen state (A,C,G,T)
    observation_probabilities = {}                                                      # Creating dictionary to store probabilities of each observation

    for observable in observation_probability[list(observation_probability.keys())[0]]: # Loop to set up initial keys on observation probabilities
        observation_probabilities[observable] = 0                                       # Initialize observables to 0

        for state in gene_probability:                                                  # Traverse every state in gene probability (A,C,G,T)
            observation_probabilities[observable] += gene_probability[state] * observation_probability[state][observable]   # Formula of observation probability. 

    return observation_probabilities                                                    # Return the dictionary with data gathered

def conditional_probability(gene_probability, observation_probabilities, gene_state, observable, observation_probability):  # Function for calculating the probability of gene state (ex. G) given observable state (ex. H) 
    
                                                                                        # We clean first the observable input. From H1 to H because we care only on the letter
    result = ''                                                                         # Initialize result that hold the cleaned observable                              
    for char in observable:                                                             # Traversing in each character of observable
        if not char.isdigit():                                                          # Check if character is not number, if not we keep it
            result += char                                                              # Update the observable to cleaned version (ex. H)
    observable = result

    if gene_state in gene_probability and observable in observation_probabilities:      # Now we check if gene state exist in gene probability and observable exist in observation probability we computed earlier
        numerator = observation_probability[gene_state][observable] * gene_probability[gene_state]  # If yes, we compute for the numeration. this is the formula of numerator in conditional probability
        denominator = observation_probabilities[observable]                             # Then we compute for the denominator. This is the overall probability of observing observable state
        
        if denominator == 0:
            return 0
        else:
            return numerator / denominator                                                  # We calculate and return the conditional probability we compute
    else:
        return None                                                                     # If either gene state or observation dont exist, we return none

def write_output(filename, sequence, initial_probability, transition_probability, observation_probability, cases):
    with open(filename, 'a') as file:                                                   # Open file with the filename
        file.write(sequence + "\n")                                                     # DNA sequence in file

        # file.write("Initial gene probabilities:\n")
        # for state in initial_probability:
        #     prob = initial_probability[state]
        #     file.write("P(" + state + "|Start) = " + str(prob) + "\n")                # For testing only

        # file.write("Transition probabilities:\n")
        # for s1 in transition_probability:
        #     for s2 in transition_probability[s1]:
        #         prob = transition_probability[s1][s2]
        #         file.write("P(" + s2 + "|" + s1 + ") = " + str(prob) + "\n")

        gene_probability = {}                                                           # We create dictionary for us to keep track the probability for each gene (A,C,G,T) at each step
        for key in initial_probability:                                                 # Traversing every gene
            gene_probability[key] = initial_probability[key]                            # Initialize them with initial probability

        for iteration in range(1, len(sequence) + 2):                                   # We traverse through each positition in sequence
            next_gene_probability = {'A': 0, 'C': 0, 'G': 0, 'T': 0}                    # Initialize dictionary to store the next gene probability

            for curr_state in next_gene_probability:                                    # We will calculate the probabilities for the next step using transition probabilities
                for prev_state in gene_probability:                                     # The outer loop is we traverse gene probability on current state and inner loop is for previous state                      
                    next_gene_probability[curr_state] += gene_probability[prev_state] * transition_probability[prev_state][curr_state]  # Formula in getting the gene probability. This gives us total probability of current state

            obs_probs = observation_probabilities(next_gene_probability, observation_probability)   # Call for observation probability

            for case in cases:                                                          # We will compute for the gene observable pairs we want to compute
                parts = case.split()                                                    # Split cases into part (G1 given H1) become 'G1' and 'H1'
                gene_state = parts[0][0]                                                # Get gene state (G from G1)
                observable = parts[2]                                                   # Get observable state (H from H1)

                gene_state_iteration = gene_state + str(iteration)                      # Concat the gene state and iteration for comparing

                if gene_state_iteration == parts[0] and observable == parts[2]:         # Check if current case matches gene and observable
                    probability_given_observable = conditional_probability(next_gene_probability, obs_probs, gene_state, observable, observation_probability) # If yes, we call contionional probability to calculate conditional probability of the gene state given observable

                    if probability_given_observable is not None:                        # If conditional probability we compute earlier is not none, meaning we have a calculated value
                        file.write(f"{gene_state}{iteration} given {observable} = {probability_given_observable}\n")    # We write the calculated value in file together with the gene state and observable

            for key in gene_probability:                                                # Update the gene probabilities for the next iteration
                gene_probability[key] = next_gene_probability[key]

def main():
    print("\nWelcome to HMM for Predicting Region of Coding DNA!")                      # Print a welcome message
    input_filename = input("Enter the filename to read: ")                              # Ask user to input the file to read

    while not os.path.exists(input_filename):                                           # If inputted file name does not exist in current directory, ask again until it became valid
        print("Error: The file does not exist. Please try again.")                      # Prompt user to try again
        input_filename = input("Enter the filename to read: ")                          # As user again

    sequences, states, observables, observation_probability, cases = read_input(input_filename) # Call for read input function and contain the following returning value we need for HMM

    output_filename = input("Enter the filename to write: ")                            # Ask user to enter file name of the output
    with open(output_filename, 'w') as file:                                            # Open output file and ensure the file is properly closed after writing
        pass

    for sequence in sequences:                                                          # We process each sequence by HMM
        initial_probability = initial_probabilities(sequence)                           # Call first initial probability function to compute the initial probability of each sequence
        transition_probability = transition_probabilities(sequence)                     # Second is we calculate the transition probability for each sequnce
        write_output(output_filename, sequence, initial_probability, transition_probability, observation_probability, cases) # We will call the function for computing first the HMM and directly writing the result in output file

    print(f"Results written to {output_filename}")                                      # Prompt user that results are written to given text file name

main()                                                                                  # Call for main function