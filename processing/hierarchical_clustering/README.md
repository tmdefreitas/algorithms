# BCB4002 - assignment2 #
##
Timothy DeFreitas

2/12/15 

## 1.  Hierarchical Clustering Algorithm
My code uses the Unweighted Pair Group Method with Arithmetic Mean (UPGMA) method to cluster the sequences. This method computes  distance between two clusters as the average distance between all pairs of elements from both sets. It is an iterative algorithm, with the two closest clusters joined into one at each step until every datapoint is in one cluster (Bottom-up agglomerative).

In my project, the underlying data are protein sequences that are read in from a file in FASTA format. The code then uses the Needleman-Wunsch Algorithm (further described in section 4 below), to compute an alignment score between the proteins. This score is used as the distance between two data points, and in my case, larger scores mean better matches (the opposite of euclidean methods).



## 2. How the Tree is Drawn
My tree is fairly similar in appearance to traditional phylogenic trees in Bioinformatic software. Starting from the top cluster, I draw a horizontal line to indicate the beginning of a cluster, and then split if there are children. The size of the offset of the split is determined mathematically by the size of the two child trees (see code). The horizontal distance to the next split is determined by the iteration step that the split occurred. Splits farther to the left in the tree correspond to clusters that were joined later in the UPGMA algorithm. The child trees are then drawn recursively.

 If no child trees exist, we are at a leaf node, and the text label of the data is drawn. 

## 3. Outside code
The Needleman-Wunsch algorithm requires a Position-Specific Scoring Matrix in order to calculate the alignment score for proteins. I found a java implementation of the BLOSUM-62 scoring matrix from the[ University of Texas](http://www.cs.utexas.edu/~mobios/cs329e/rosetta/src/Blosum.java) and adapted it to support my code. This section appears at the end of assignment2.pde and is clearly identified. All other classes and algorithms were developed from scratch.
 
## 4. Beyond the Requirements ##

###Biological Improvements

My code uses as input protein sequences, encoded in FASTA files. These files are the standard for text-based encoding of protein and nucleotide sequences, and allows my tree to be constructed from an arbitrary FASTA dataset. In addition I implmented the [Needleman-Wunsch algorithm](http://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm) which computes a global alignment score between two proteins.

The algorithm is significant because it uses principles of evolution to compare how similar two protein sequences are. Evolution is caused by random mutation, but it is widely understood that some sequences are conserved, and further, some Amino Acids (The letters of a protein sequence) are very similar. Thus, some mutations are more likely to survive than others. Similarly insertions and deletions are penalized using a gap penalty. Thus, the Needleman-Wunsch algorithm uses a Position-Specific Scoring Matrix (I use the BLOSUM-62 matrix, PAM is the other common one), which better distinguishes between "likely" mutations that simple levenshtein or manhattan distance. Thus the resulting hierarchy suggests which proteins are more closely related to eachother.

In my example file, I collected samples of hemoglobin subunit alpha from 10 different species. All mammals have some variant of the protein, but n the resulting tree, you can see that Humans are most similar to M. Mulatta (Rhesus Macaque), and more distantly to rabbits and rhinos. 

###Technical Improvements

Though the majority of my additions were biological, I also added mouseover highlighting to the tree, because this tree can be hard to read (and hard to distinguish between groups). Placing your mouse over the branching area of a cluster will highlight the cluster and all leaves in the cluster, so one can visually distinguish between the groups. 

You can also change the highlighter color or any of the other graphical properties by editing lines 13-19 of assignment2.pde.



## 5. Running the Program ##

The code for this program is contained entirely within **assignment2.pde**. Open the file in Processing.
The input for this program is *hemoglobin.fasta*, which contains the protein sequences of 10 example proteins that are used to construct the tree. Because of the way Processing is set up, you will have to manually change the location of the *hemoglobin.fasta* file on line 10 of assignment2.pde:

`final String DATA_FILE = "C:\\Users\\Tim\\Documents\\Biovis\\assignment2\\hemoglobin.fasta";`

Input the complete path to *hemoglobin.fasta* and run the sketch. The tree will then appear in the window (and a debugging version will appear in the console).