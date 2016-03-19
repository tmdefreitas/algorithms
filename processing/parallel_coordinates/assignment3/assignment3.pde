//Timothy DeFreitas
//assignment3.pde
//2/16/15

import java.util.LinkedList;
import java.util.EnumMap;
import java.io.*;
import java.util.Arrays;

import java.awt.Frame;
import java.awt.BorderLayout;
import controlP5.*;

//Data source
final String DATA_FILE = "C:\\Users\\Tim\\Documents\\Biovis\\assignment3\\wdbc.data";

//Graphics Parameters
final int WIN_X = 1000;
final int WIN_Y = 800;
final int BINS = 5;

//Globals
ControlP5 cp5;
ControlFrame cf;
EnumMap<FEATURE, Double> minima;
EnumMap<FEATURE, Double> maxima;
LinkedList<CancerSample> l;



//Setup Processing
void setup(){
  frameRate(25);
  colorMode(RGB,255);
  size(WIN_X, WIN_Y);
  cp5 = new ControlP5(this);
  
  cf = addControlFrame("Controls", 300, 400);
  //Test loading of .data file
  l = loadFile(DATA_FILE);
  calculateRanges(l);
  
  println(minima);
  println(maxima);
  println(l.size());
  
  background(255);
}

//Drawing loop
void draw(){
  background(255,255,255);
  LinkedList<FEATURE> activeFeatures = cf.activeFeatures();
  if (activeFeatures.size() == 1){
    histogram(activeFeatures.get(0), l);
    
  }
  
  
  axes(activeFeatures);
  for (CancerSample cs: l){
    cs.render(activeFeatures);
  }
}

//Histogram for if one feature is selected
void histogram(FEATURE f, LinkedList<CancerSample> l){

  
  //Get x width for histogram
  double domain_min = minima.get(f);
  double domain_max = maxima.get(f);
  
  //Count samples in each bin
  int[] countM = new int[BINS];
  int[] countB = new int[BINS];
  for (CancerSample cs : l){
    int inBin = floor((float) (cs.getFeature(f)/(domain_max-domain_min)));
    if (cs.isMalignant()) {
      countM[inBin] += 1;
    } else {
      countB[inBin] += 1;
    }
    
  }// for cs
  
  int range_max = max(max(countM), max(countB));
  
  //----Draw axes---------
  color c = color(0,0,0); //Black axes
  stroke(c);
  fill(c);
  
  line(50, 50,50, WIN_Y - 50); //Y-axis
  line(50, WIN_Y-50, WIN_X - 50, WIN_Y-50); //X-axis
  
  //Y axis labels
  textAlign(CENTER);
  textSize(20);
  text("0", 25, WIN_Y-50); //Lowest value (always 0)
  text(Integer.toString(range_max), 25, 50); //highest value (maximum of the histogram);
  
  //X axis labels
  for (int i = 0; i <= BINS; i++){
    String label = Double.toString(i*(domain_max-domain_min)/BINS); //Equal width binning
    text(label, 50 + i*(WIN_X-100)/BINS,WIN_Y-25);
  }//for i
  
  //Chart title
  text("Histogram of " + f.toString(), WIN_X/2, 50);
  
  //Legend labels
  textAlign(LEFT);
  text("Malignant", WIN_X -125, 50);
  text("Benign", WIN_X - 125, 80);
  
  //------plot Malignant counts------
  c = color(255, 0,0);
  stroke(c);
  fill(c);
  //color legend
  rect(WIN_X-150, 30, 20,20);
 
  float bar_width = (WIN_Y-100)/(2*BINS);
  
  for (int m = 0; m < BINS; m++){
    //plot a rectangle to the left side
    int h = countM[m]*(WIN_Y - 100)/range_max;
    int x = 50+m*(WIN_X-100)/BINS;
    rect(x, WIN_Y - 50 - h, bar_width, h);
  }//for m
  
  
  
  //-------plot Benign counts---------
  c = color(0,0,255);
  stroke(c);
  fill(c);
  //color legend
  rect(WIN_X-150, 60,20,20);
  
  for (int b = 0; b < BINS; b++){
    //plot a rectangle to the left side
    int h = countB[b]*(WIN_Y - 100)/range_max;
    float x = 50 +bar_width + b*(WIN_X-100)/BINS;
    rect(x, WIN_Y - 50 - h, bar_width, h);
  }//for b  
  
  
  
}//Histogram


/**
  Class to hold breast cancer data entries.  
*/

class CancerSample {
  private String id;
  private boolean malignant;
  private EnumMap<FEATURE, Double> featureMap;
  
  //Constructs a cancer sample from a .data file
  
  CancerSample(String data) {
    String[] features = data.split(",");
    
    //Id is first entry
    this.id = features[0];
    
    if (features[1].equals("M")){
      this.malignant = true;
    } else {
      this.malignant = false;
    }
    
    int idx = 2;
    
    this.featureMap = new EnumMap<FEATURE, Double>(FEATURE.class);
    
    for (FEATURE f: FEATURE.values()){
      //Parse the next feature as 
      this.featureMap.put(f, Double.parseDouble(features[idx]));
      idx++;
    }
    
  }//CancerSample(String data)
  
  //Get the value for a feature
  public double getFeature(FEATURE f){
    return this.featureMap.get(f);
  }
  
  public boolean isMalignant(){
    return this.malignant;
  }
  
  public String getID(){
    return this.id;
  }
  
  //Render based on the current active featureset
  public void render(LinkedList<FEATURE> activeFeatures){
    if (activeFeatures.size() < 2) return;
    color c;
    if (this.isMalignant()){
      c = color(255,0,0);
    } else {
      c = color(0,0,255);
    }
    stroke(c);
    
    int xWidth = (WIN_X - 100) / (activeFeatures.size() - 1);
    double left_scale = maxima.get(activeFeatures.get(0)) - minima.get(activeFeatures.get(0)); 
    int xcoord = 50;
    
    //Draw lines for active features 
    for (int i = 1; i < activeFeatures.size(); i ++){
      
      double right_scale = maxima.get(activeFeatures.get(i)) - minima.get(activeFeatures.get(i)); 
      double left_y_offset = (WIN_Y - 75)*(this.featureMap.get(activeFeatures.get(i-1))- minima.get(activeFeatures.get(i-1)))/left_scale;
      double right_y_offset = (WIN_Y - 75)*(this.featureMap.get(activeFeatures.get(i)) - minima.get(activeFeatures.get(i)))/right_scale;

      int left_y = new Double(WIN_Y - 25 - left_y_offset).intValue();
      int right_y = new Double(WIN_Y - 25 - right_y_offset).intValue();
      if (left_y > WIN_Y) {
        println(this.featureMap);
      }
      
      
      line(xcoord, left_y, xcoord+xWidth, right_y);
      left_scale = right_scale;
      xcoord += xWidth;
      
    }
  }
  
}//CancerSample

//Loading method for CancerSamples
LinkedList<CancerSample> loadFile(String fromFile){
  LinkedList<CancerSample> l = new LinkedList<CancerSample>();
  
  BufferedReader br;
  
  try {
    br = new BufferedReader(new FileReader(fromFile));
    
    String line = br.readLine();
    
    while (line != null){
      l.add(new CancerSample(line.trim()));
      line = br.readLine();
    }
  } catch (Exception e) {
    println("Failed to load file, exception:" + e);
  }
  
  return l;
}

//Puts the ranges of each feature into global maxima/minima 
void calculateRanges(LinkedList<CancerSample> samples){
  maxima = new EnumMap<FEATURE, Double>(FEATURE.class);
  minima = new EnumMap<FEATURE, Double>(FEATURE.class);
  
  for (CancerSample cs: samples){
    for (FEATURE f: Arrays.copyOfRange(FEATURE.values(), 0,10)){
      Double thisValue = cs.getFeature(f);
      
      Double max = maxima.get(f);
      Double min = minima.get(f);
      
      if (max == null || max < thisValue){
        maxima.put(f, thisValue);
      }
      if (min==null || min > thisValue){
        minima.put(f, thisValue);
      }
    }
  }
      
}

//Draw the axes for the parallel coordinates
public void axes(LinkedList<FEATURE> activeFeatures){
  if (activeFeatures.size() >= 2){
    //Amount of spacing between each parrallel coordinate
    int xWidth = (WIN_X - 100) / (activeFeatures.size() - 1);
    
    int xcoord = 50;
    stroke(0,0,0);
    fill(0,0,0);
    textSize(20);
    textAlign(CENTER);
    
    for (FEATURE f : activeFeatures){
      //Draw vertical line
      line(xcoord, 50, xcoord, WIN_Y - 25);
      
      //Add range numbers
      String max = String.valueOf(maxima.get(f));
      String min = String.valueOf(minima.get(f));
      text(max, xcoord, 40);
      text(min, xcoord, WIN_Y - 10);
      
      //Add text label
      text(f.toString(), xcoord, 15);
      
      xcoord += xWidth;
      
    }//for
  }//if active features > 2
  
} //axes


/**
  External control frame, taken (and modified) from control frame example: http://www.sojamo.de/libraries/controlP5/examples/extra/ControlP5frame/ControlP5frame.pde
  
  Creates a new processing applet inside a new java frame with a controlp5 object
  */
ControlFrame addControlFrame(String theName, int theWidth, int theHeight) {
  Frame f = new Frame(theName);
  ControlFrame p = new ControlFrame(this, theWidth, theHeight);
  f.add(p);
  p.init();
  f.setTitle(theName);
  f.setSize(p.w, p.h);
  f.setLocation(100,100);
  f.setResizable(false);
  f.setVisible(true);
  return p;
}


class ControlFrame extends PApplet {
  int w, h;
  int abc = 100;
  ControlP5 cp5;
  CheckBox checkBox;
  Object parent;
  
  public void setup(){
    size(w,h);
    frameRate(25);
    cp5 = new ControlP5(this);
    checkBox = cp5.addCheckBox("Checkbox");
    checkBox.setColorForeground(color(120))
       .setColorActive(color(255))
       .setColorLabel(color(255))
       .setSize(40,40)
       .setItemsPerRow(2)
       .setSpacingColumn(60)
       .setSpacingRow(30);
       
    for (FEATURE f: Arrays.copyOfRange(FEATURE.values(), 0,10)){
      checkBox.addItem(f.toString(), f.ordinal());
    }
    
  }
    
  public void draw(){
    background(abc);
  }
  
  private ControlFrame(){
  }
  
  public ControlFrame(Object theParent, int theWidth, int theHeight){
    parent = theParent;
    w = theWidth;
    h = theHeight;
  }
  
  public ControlP5 control(){
    return cp5;
  }
  //Returns a list of the current features to visualize with parallel coordinates
  public LinkedList<FEATURE> activeFeatures(){
    LinkedList<FEATURE> onFeatures = new LinkedList<FEATURE>();
    for (FEATURE f: Arrays.copyOfRange(FEATURE.values(), 0,10)){
      if (this.checkBox.getState(f.toString())){
        onFeatures.add(f);
      }
    }//for
    
    return onFeatures;    
  }
    
  
 
}//ControlFrame

