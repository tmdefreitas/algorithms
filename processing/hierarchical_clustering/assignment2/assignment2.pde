//Timothy DeFreitas
//BCB4002 - assignment2.pde
//2/12/15

import java.util.LinkedList;
import java.io.*;
import java.lang.*;

//Data source
final String DATA_FILE = "C:\\Users\\Tim\\Documents\\Biovis\\assignment2\\hemoglobin.fasta";

//Graphics Parameters
final int WIN_X = 600;
final int WIN_Y = 600;
color NORMAL_STROKE = color(0,0,0); //black
color BACKGROUND = color(255,255,255); // white
color HIGHLIGHT = color(255,255,0); //gold
int SCALE = 20;
int STROKE_WEIGHT = 2;

//Default Gap Penalty
final int DEFAULT_GAP = -5;

//Globals
Tree t;

//Processing setup function
void setup(){
  colorMode(RGB, 255,255,255);
  stroke(NORMAL_STROKE); // black
  strokeWeight(STROKE_WEIGHT);
  fill(NORMAL_STROKE);
  LinkedList<FastaSequence> fs = loadFile(DATA_FILE);
  t = new Tree(fs);
  System.out.println(t.toString());
  size(WIN_X, WIN_Y);
  background(BACKGROUND);
  
}

//Drawing loop
void draw(){
  fill(NORMAL_STROKE);
  stroke(NORMAL_STROKE);
  t.drawTree(0,WIN_X/2,SCALE);
}

/**
  Class to hold protein sequence. FASTA files have a header line and the 
  protein sequence stored as amino acid letter codes.
  */
class FastaSequence {
  
  String description;
  String sequence;
  
  FastaSequence(String description, String sequence){
    this.description = description;
    this.sequence = sequence.toUpperCase();
    
  }//FastaSequence()
  //Required to be overrited so that I can use the remove() method on Lists
  boolean equals(Object o){
    if (!(o instanceof FastaSequence)){
      return false;
    } else {
      FastaSequence fs = (FastaSequence) o;
      return this.description.equals(fs.description) && this.sequence.equals(fs.sequence);
    }
  }
  
  String getHeading(){
    return this.description;
  }
  
  int length(){
    return this.sequence.length();
  }
  
  char charAt(int index){
    return this.sequence.charAt(index);
  }
  
  String toString(){
    return description + "\n" + sequence;
    
  }
  /**
  
  Computes the alignment score using the Needleman-Wunsch algorithm: http://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm
  for a global alignment between the two sequences using the Blosum-62 matrix at the end of this file.
   */
  
  int alignmentScore(FastaSequence s, int gapPenalty){
    //Array to compute score
    int[][] array = new int[this.length()+1][s.length()+1];
    
    //Initialize first row to zeros
    for (int i = 0; i < this.length() + 1; i++){
      array[i][0] = gapPenalty*i; 
    }
    
    //Initialize first column to zeros
    for (int j = 0; j < s.length()+1; j++){
      array[0][j] = gapPenalty*j;
    }
    
    for (int i = 1; i < this.length() + 1; i++) {
      for (int j = 1; j < s.length()+1; j++){
        //Calculate possible scores for this box
        int matchScore = array[i-1][j-1] + Blosum.getDistance(this.charAt(i-1), s.charAt(j-1));
        int delete = array[i-1][j] + gapPenalty;
        int insert = array[i][j-1] + gapPenalty;
        
        //Select the best value 
        array[i][j] = max(matchScore, delete, insert);
      }
    }
    
    //Return the bottom corner, which is the global alignment score for these sequences
    return array[this.length()][s.length()];
    
  }
  
  //Default version
  int alignmentScore(FastaSequence s){
    return this.alignmentScore(s, DEFAULT_GAP);
  }
    
    
}//class FastaSequence

//Loading method for FastaSequences
LinkedList<FastaSequence> loadFile(String fromFile){
  LinkedList<FastaSequence> l = new LinkedList<FastaSequence>();
  BufferedReader br;

  try {
    br = new BufferedReader(new FileReader(fromFile));
    StringBuilder sb = new StringBuilder();
    String line = br.readLine();
    String header = null;
    
    
    while (line!=null) {
      //We hit a new sequence
      if (line.startsWith(">")){
        if (header != null) {
          //Log the previous sequence
          l.add(new FastaSequence(header, sb.toString()));
          
        }
        //Start a new sequence
        header = line.substring(1).trim();
        sb = new StringBuilder();      
        
      } else {
        sb.append(line.trim()); //Add part of the sequence
      }
      
      line = br.readLine();
    }
    
    //Add whatever remains of the final sequence
    if (header != null){ 
      l.add(new FastaSequence(header, sb.toString()));
    }
    
  } catch (Exception e) {
    println("Failed to load file exception:" + e);
  }//catch
  
  return l;

}

/**
  Node class for creating the phylogenetic tree.
*/
class TreeNode {
  private LinkedList<FastaSequence> sequences;
  private TreeNode left;
  private TreeNode right;
  
  //Iteration the nodes were merged
  private int age;
  
  //Distance between the two children
  private float similarity;
  
  //Drawing parameters (needed for mouseOver)
  private int x, y;
  
  //Default constructor (not used)
  TreeNode(LinkedList<FastaSequence> sequences, TreeNode left, TreeNode right){
    this.sequences = sequences;
    this.left = left;
    this.right = right;
  }
  
  //Simplified constructor, for initializing list of nodes
  TreeNode(FastaSequence fs){
   this.sequences = new LinkedList<FastaSequence>();
   this.sequences.add(fs);
   this.left = null;
   this.right = null;
   this.age = 1;
   this.similarity = -1.0;//Meaningless, since there is no left and right
   
    this.x = 0;
    this.y = 0;
  }
  //Most helpful Constructor, used by the UPGMA algorithm to join two (non-null) nodes
  TreeNode(TreeNode left, TreeNode right, int iteration){
    this.sequences = new LinkedList<FastaSequence>();
    this.sequences.addAll(left.sequences);
    this.sequences.addAll(right.sequences);
    this.left = left;
    this.right = right;
    this.similarity = left.upgma(right);
    this.age = iteration;
    this.x = 0;
    this.y = 0;
  }
  //Must be overwritten to use remove() method
  boolean equals(Object o){
    if (!(o instanceof TreeNode)){
      return false;
    } else {
      TreeNode other = (TreeNode) o;
      boolean result = true;
      
      if (this.right == null) {
        result = result && other.right == null;
      } else {
        result = result && this.right.equals(other.right);
      }
      
      if (this.left == null) {
        result = result && other.left == null;
      } else {
        result = result && this.left.equals(other.left);
      }
      if (this.sequences == null) {
        result = result && other.sequences==null;
      } else {
        result = result && this.sequences.equals(other.sequences);
      } 
 
      return result;
    }
  }
  
  //Calculates the unweighted pairwise group method with arithmetic mean (UPGMA) distance between two groups
  float upgma(TreeNode other) {
    float score = 0.0;
    int comparisons = this.sequences.size()*other.sequences.size();
    for (FastaSequence fs_i: this.sequences){
      for (FastaSequence fs_j: other.sequences){
        score += fs_i.alignmentScore(fs_j)/comparisons; 
        
      }//fs_j
    }//fs_i   
    return score;
  }
  
  //Rendering algorithm
  void drawTree(int x, int y, int scale, color c){
    this.x = x;
    this.y = y;
    
    //If this node, or it's parent node is highlighted (color is already changed)
    //Use the highlight color
    if (this.mouseOver()){
      c = HIGHLIGHT;
    }
    fill(c); //set colors
    stroke(c);    
    textSize(scale);
    
    
    if (this.left==null & this.right==null) { 
      //This is a leaf node
      line(x, y, x + scale, y);
      //Text label
      text(this.sequences.get(0).getHeading(), x + 3*scale/2, y + scale/2);      

    } else {
      //This is an interior node
      //Space between labels
      int padding = scale/8;      
      
      
      line(x,y,x+scale,y); //Draw handle for this tree
      
      //Figure out how much room to move leave each side      
      int left_y_offset = y - scale/2 + padding;
      int right_y_offset = y + scale/2 + padding;    
      if (this.left.age > 1){ //The left tree definitely has a right sub tree
        left_y_offset = y - this.left.right.size()*scale + padding; //The size of the right side of the left subtree        
      }      
      if (this.right.age > 1) {//The right tree has a sub tree
        right_y_offset = y + this.right.left.size()*scale + padding;
      }
      
      //Now that we know the vertical offset, draw in the vertical bar
      line(x + scale,right_y_offset, x + scale,left_y_offset);
      
      //where to start drawing the subtrees (X)
      int left_x_offset = x + scale +(this.age - this.left.age - 1)*scale; 
      int right_x_offset = x + scale +(this.age - this.right.age -1)*scale;
      
      //Draw the intermediate horizontal bars
      line(x + scale, right_y_offset, right_x_offset,right_y_offset);
      line(x + scale, left_y_offset, left_x_offset, left_y_offset);
      
      //Finally, recursively draw the left and right subtrees
      this.left.drawTree(left_x_offset, left_y_offset, scale, c);
      this.right.drawTree(right_x_offset, right_y_offset,scale, c);
      
    }
    
  }//drawTree
  
  //Returns true if the mouse is over this tree
  boolean mouseOver(){
    if (mouseY < this.y + SCALE && mouseY > this.y - SCALE 
        && mouseX >= this.x && mouseX < this.x +SCALE){
       return true;
     }
    if(this.left == null && this.right == null){
      //Leaf node
      if (mouseY > this.y - SCALE/2 && mouseY < this.y + SCALE/2 &&
          mouseX > this.x + SCALE && mouseX < this.x + SCALE*4){
            return true;
          }
    }
    return false;
  }
    
  //toString() helpful for debugging
  String toString(){
    String [] children;
    StringBuilder sb = new StringBuilder("");
       
    if (this.left == null && this.right == null){
      sb.append("----" + this.sequences.get(0).getHeading() + "\n");
    }
    else if (this.right == null && this.left !=null){
      children = this.left.toString().split("\n");
      //Add one string, and a bunch of blank ones
      sb.append("+---");
      sb.append(children[0] + "\n");
      for (int i = 1; i < children.length; i ++){
        sb.append("    "); //four spaces
        sb.append(children[i]+"\n");
      }
    } else if (this.left == null && this.right != null){
      children = this.right.toString().split("\n");
      //Add one string, and a bunch of blank ones
      sb.append("+---");
      sb.append(children[0] + "\n");
      for (int i = 1; i < children.length; i ++){
        sb.append("    "); //four spaces
        sb.append(children[i]+"\n");
      }   
    } else {
      children = this.left.toString().split("\n");
      //Add one string, and a bunch of blank ones
      sb.append("+---");
      sb.append(children[0] + "\n");
      for (int i = 1; i < children.length; i ++){
        sb.append("|   "); //"|" plus three spaces
        sb.append(children[i]+"\n");
      }
      sb.append("|\n");
      children = this.right.toString().split("\n");
      sb.append("+---");
      sb.append(children[0] + "\n");
      for (int i = 1; i < children.length; i ++){
        sb.append("    "); //four spaces
        sb.append(children[i]+"\n");
      }       
            
    } //else
    
    return sb.toString();
  } //toString()
    
  /**
    Returns the relative size of the tree node, in terms of the number of sequences contained in it.
  */
  int size()  {
    return this.sequences.size();
  }
  /**
    Returns the age(iteration that the node was created) of the root node
  */
  int rootAge(){
    return this.age;
  }
  /**
    returns the distance between the two children of this node. Helpful for rendering purposes
  */
  float similarity(){
    return this.similarity;
  }
  
}//TreeNode

class Tree {
  TreeNode root;
  
  //Construct a Tree using the given FastaSequences
  Tree(LinkedList<FastaSequence> fseqs){
    ArrayList<TreeNode> forest = new ArrayList<TreeNode>();
    for (FastaSequence fs:fseqs){
      forest.add(new TreeNode(fs));
    }
    int iteration = 1;
    
    //Each iteration will combine two TreeNodes into one based on the best upgma score between nodes
    //UPGMA Algorithm
    while (forest.size() > 1){
      iteration++;//increment iteration count
      TreeNode best1 = forest.get(0);
      TreeNode best2 = forest.get(1);
      float best_upgma = best1.upgma(best2);
      
      //Calculate all pairwise cluster distances
      for (int i = 0; i < forest.size(); i++) {
        for (int j = i+1; j < forest.size(); j++) {
          float this_upgma = forest.get(i).upgma(forest.get(j));
          
          //Because we are using similarity scores, the largest one is the best
          if (this_upgma > best_upgma){
            best1 = forest.get(i);
            best2 = forest.get(j);
            
            best_upgma = this_upgma;
          }
        }
      }//endfor
      //Merge the two closest nodes, and replace them with the merged node.
      TreeNode newTreeNode = new TreeNode(best1, best2, iteration);
      forest.remove(best1);
      forest.remove(best2);
      forest.add(newTreeNode);
      
    }//endwhile
    
    //When the forest contains only one node, it necessarily contains all other nodes.
    root = forest.get(0);
    
  }
  
  /** 
    Render the tree
    
  */
  void drawTree(int x, int y, int scale) {
    this.root.drawTree(x,y,scale, NORMAL_STROKE);//Default color, the TreeNode will override if mouseOver()
  }
  
  //toString
  String toString(){
    return this.root.toString();
  }
  
}



//--------------Outside Code ---------------------------------------/
// Scoring matrix class, taken (with some comments omitted) from Univesity of Texas:
//http://www.cs.utexas.edu/~mobios/cs329e/rosetta/src/Blosum.java
//**I made some slight modifications to shorten the code/ Remove extra error handling 
// --> Must ensure properly .fasta files or IndexOutofBoundsExceptions will occur

final public static class Blosum{

    /*
     * Array representation of Blosum-62 matrix 
     * Refer to above matrix for corrseponding amino acids
     * i.e. score(A, R) corresponds to  matrix[0][1]=matrix[1][0]=-1
    */  
    private static final int[][] matrix = {
  { 4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0},
  {-1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3},
  {-2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3},
  {-2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3},
  { 0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1},
  {-1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2},
  {-1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2},
  { 0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3},
  {-2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3},
  {-1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3},
  {-1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1},
  {-1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2},
  {-1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1},
  {-2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1},
  {-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2},
  { 1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2},
  { 0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0},
  {-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3},
  {-2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1},
  { 0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4}};


  // quick and dirty equivalent of typesafe enum pattern, can also use HashMap
    // or even better, EnumMap in Java 5. 
    // This code is for Java 1.4.2, so we will stick to the simple implementation
    private static int getIndex(char a) {
      // check for upper and lowercase characters
      switch ((String.valueOf(a)).toUpperCase().charAt(0)) {
        case 'A': return 0;
        case 'R': return 1;
        case 'N': return 2;
        case 'D': return 3;
        case 'C': return 4;
        case 'Q': return 5;
        case 'E': return 6;
        case 'G': return 7;
        case 'H': return 8;
        case 'I': return 9;
        case 'L': return 10;
        case 'K': return 11;
        case 'M': return 12;
        case 'F': return 13;
        case 'P': return 14;
        case 'S': return 15;
        case 'T': return 16;
        case 'W': return 17;
        case 'Y': return 18;
        case 'V': return 19;
        default: return -1; //Invalid code
      }
    }
    
    private static int getDistance(int i, int j) {
      return matrix[i][j];
    }

    public static int getDistance(char a1, char a2) {
      // toUpper
      return getDistance(getIndex(a1), getIndex(a2));    
    }
}
//--------------End Blosum class----------------------------//
