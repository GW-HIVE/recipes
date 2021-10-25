Use-case Examples
-----------------

There are three example data files as well as an example recipe to demonstrate how actions are used on real data. 


**Example Recipes and the actions they contain:** 

example_part1.json
	- union
	- jointables
	- makecombofield
	- splitcol
example_part2.json
	- all actions from part 1
	- filterin
	- filterout
	- expandrows
	- addnormalizedcol
	- addaveragecol

Actions not yet added:
	- transposecol
	- addconstantfield


The example recipe takes all of the example data files as input. 

In the recipe each of the actions are used. To see the outcome of any one action, change the tableid option listed in the outputfiles section. 

Be sure to change the name of the outputfile in the filepath option in the outputfiles section, otherwise each time you run the recipe it will overwrite the previous output file.

**Example Data: 
example_mutations_1.csv
example_mutations_2.csv
example_mapping.csv**

The contents of the data files are the following:

example_mutations_1.csv and example_mutations_2.csv contain SNP data gatehred from UCSC genome browser. 
- There are both unique and overlapping values between these two mutation files 

example_mapping.csv contains uniprot terms and their corresponsing gene symbols, which allows us to show how you can map terms from one table to terms from one table to another
