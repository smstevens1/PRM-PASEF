# PRM-PASEF

This python script requires two input files within same directory as script: 
1) The gene target list in .txt format (e.g., exported from IPA or user-generated)
2) DIA output either from DIA-NN (v. 1.8.1) or TIMS DIA-NN (.tsv format)

Once executed, the user will be asked several prompts (description below).
1) Enter name of IPA file (txt format):
Specify the gene input list (include .txt extension) and ensure appropriate format (example IPA input is provided).
2) Enter name of MS file (tsv format):
Specify the DIA output file (include .tsv extension) from either DIA-NN or TIMS DIA-NN.
3) Enter MS file type (TIMSDIANN/DIANN):
Specify if this output file was generated from TIMSDIANN or DIANN.
4) Enter name of output file (csv format):
Name your output file (include .csv extension)
5) Enter RT range (seconds): [Suggested:240]
Specify the retention time range window for the precursor ions (our suggested is 4 min but could change based on your experiment).
6) Enter IM tolerance (percent): [Suggested:4]
Specify the ion mobility measurement tolerance (our suggested window is 4% of the measured IM value from the DIA output).
7) Enter charge-state relative intensity tolerance (fraction): [Suggested:0.2]
In the case of a precursor ion with multiple charge states, only the predominant charge state will be selected if the lower 
abundance charge state is below 0.2 relative to other charge state. If this criterion is not met, both charge states are kept.
8) Enter undigested sequence relative intensity tolerance (fraction): [Suggested:0.1]
In the case of a peptide with missed cleavages, the missed cleaved peptide will be dropped if below this intensity tolerance. If
the missed cleavage is higher abundance, both peptides will be dropped.
9) Enter a positive integer to print only the top N highest intensity sequences. (Enter -1 for all sequences):
Precursors can be limited to a certain number based on intensity (e.g., 2 would indicate to only include the top 2 most intense 
peptide precursors). If -1 is entered, then all peptides will be included without intensity filtering.
10) Drop seqeunces contaning Methionine (YES/NO):
If yes, then all methionine-containing sequences will be dropped from the PRM list.

The .csv file can be directly imported into Bruker Compass HyStar for PRM-PASEF method generation.
