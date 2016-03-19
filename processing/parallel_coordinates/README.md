## BCB4002 - Assignment 2
Timothy DeFreitas

2/15/15

## 1. Dataset

I used the [breast cancer diagnostic dataset](http://http://archive.ics.uci.edu/ml/datasets/Breast+Cancer+Wisconsin+(Diagnostic)), one of the datasets listed on Matt Ward's course website. Specifically, the data represents attributes measured from digital images of a fine needle aspirate of breast masses. There are a total of ten real-valued features in addition to a classification of benign or malignant, as well as an identification number for each sample. More information on each feature is provided in `wdbc.names`. The 569 data points themselves are provided in `wdbc.data`.

## 2. Transformations 
The data is not significantly altered. With parallel coordinates, each dimension is normalized linearly to the available screen space by first calculating the possible range. The other factor I altered is to color each line based on the classification; red --> malignant, blue --> benign.

## 3. Versatility
My visualization will work (in principle) for a dataset of any size, though extremely large datasets may be hard to interpret, as individual lines wil start to blur together. However, the coloring may still show the trends between malignant and benign samples, since a color gradient would still be evident.

## 4. Outside Code

All of the code for the parallel coordinate representation is original, but I borrowed some code from the ControlP5 tutorials to get the checkbox controls to work. I modified [this example](http://www.sojamo.de/libraries/controlP5/examples/controllers/ControlP5checkBox/ControlP5checkBox.pde) on checkboxes and used [this tutorial](http://www.sojamo.de/libraries/controlP5/examples/extra/ControlP5frame/ControlP5frame.pde) to put the controls in a separate window. 

## 5. Beyond the Requirements

###Biological Additions

The primary biological improvement is the use of color to designate a classification target, in this case whether a biomass sample was malignant or benign. A researcher can use my visualizations to see differences in the values for each feature between benign and malignant samples. I chose a dataset with numerical attributes so the parallel coordinate representation would more easily differentiate the cases. The technical additions also allow for biological experimentation, since a researcher can focus on one attribute specifically (e.g. size), or a particular set of attributes that may signify something relevant.

Unfortunately, this dataset has fairly skewed distributions on most attributes, resulting in some less-than-useful histograms. However, a more varied distribution would show more interesting results.

###Technical Additions

I made a number of additions to allow the user to explore the data. My code supports two different views: a histogram view for comparing the distributions of one feature between malignant and benign, and a parallel coordinate representation that plots any user selected features and resizes the chart accordingly. Malignant and benign are color coded, and the ranges are shown so the user can see patterns of malignant versus benign across any number of features.

I also implemented a separate window using the ControlP5 library with checkboxes for each attribute to allow the user to change the visible features on the fly. The user can easily see and change the features at the same time he is looking at the output window.

## 6. Running the Program

The code for this program is primarily written in **assignment3.pde**. Open this file in processing. As with project 2, this project requires a dataset to first be read in. Change the value of DATA_FILE on line 15 to match the path to **wdbc.data** on your filesystem:


`final String DATA_FILE = "C:\\Users\\Tim\\Documents\\Biovis\\assignment3\\wdbc.data";`

**This assignment uses the ControlP5 library.** Make sure the library is installed via the library Manager, and then run the sketch. 

Two windows will appear, one with a list of available features, and one output window.  By selecting one feature, a histogram of the distribution of values for that feature will appear in the output window.

By selecting two or more features, a parallel coordinate visualization will appear in the output window. Feel free to choose any combination of input features. Red lines indicate malignant masses, and blue lines are benign samples.

