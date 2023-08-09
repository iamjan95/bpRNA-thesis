The code base for bpRNA consists of 2 projects:-
Both of them are built on visual studio.

1. bpRNA :- In this we convert the RNA sequence to bpRNA sequence. 
There are 3 files in the project : - 
	- FileInputs - where file Handling (reading and writing) is done
	- bpRNASequencer - has the core logic 
	- main 

The format for running is :-
exe-file Input-filename Output-filename

Input file can be bpseq file or dot bracket file.

2. Align :- In this we align two bpRNA sequences and output their alignment and distance.
There are 4 files in this project:-
	- FileHandler - where the bpRNA sequence and dotbracket sequence is read
	- Align - where the alignment logic is written
	- Constants - has a lot of constants that are used in the project
	- main

The format for running is :-
exe-file Input-filename1 input-filename2 <bandwidth>

The input files in this case are the outputs obtained from 1st project