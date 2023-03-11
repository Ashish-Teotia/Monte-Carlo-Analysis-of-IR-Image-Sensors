# Monte-Carlo-Analysis-of-IR-Image-Sensors
Mercury Cadmium Telluride (HgCdTe) has been the dominant material for fabricating high performance infrared (IR) detectors in the complete useful range of IR radiations for a last few decades due to its tunable band-gap and high sensitivity and for that reason it continues to be the material of choice for future generation of IR detectors also  . However, being fragile and defect prone material, controlling its properties over a large area is a challenging task. Therefore, maintaining a reasonable uniformity in an imaging array of HgCdTe photodiodes is difficult.
In this project based on data generation using programming ,visualization and interpretation we present a method of analysing the performance non-uniformity of HgCdTe photodiode arrays for infrared imaging applications using programming and visualization techniques. 

Monte Carlo Simulation, also known as the Monte Carlo Analysis or a multiple probability simulation, is a mathematical technique, which is used to estimate the possible outcomes of an uncertain event. 
Unlike a normal forecasting model, Monte Carlo Simulation in this project predicts a set of outcomes based on an estimated range of values versus a set of fixed input values of input parameters .These input parameters values are normally distributed values. In other words, a Monte Carlo Simulation builds a model of possible results by leveraging a probability distribution, such as a uniform or normal distribution, for any variable that has inherent uncertainty like the various input parameters .In this project, the variations in the electrical properties which are input parameters that govern the current conduction mechanisms of a photodiode are responsible for the pixel-to-pixel non-uniformity in a photodiode array, which is highly undesirable from image quality point of view .
In a typical Monte Carlo experiment, this exercise can be repeated thousands of times to produce a large number of likely outcomes.So therefore the outcomes received are then visualised to depict the uncertainity and efficiency like quantum efficiency and so on.

The various steps involved in this project via coding in C++ ,Python etc.
1.Data Generation
Careful and comprehensive data preparation ensures analysts trust, understand, and ask better questions of their data, making their analyses more accurate and meaningful. From more meaningful data analysis comes better insights and, of course, better outcomes. I have generated 1,000 samples of normally distributed values for each of these parameters with the average and the range given in the Table below. The sample space is fairly large for statistical treatment and then further performing the task.The code is available in Values Generated cpp file.
![image](https://user-images.githubusercontent.com/125439405/224496428-203d0a3d-ce9e-4b6c-b21f-634afff8f120.png)

2.Input material Parameter distribution
To confirm and to visualise the input parameters ,the values are visualised via a Python Code.The Python code along with the plots are in Input Parameters Plot Python file.

3.Correlation Analysis
Correlation Analysis is statistical method that is used to discover if there is a relationship between two variables/datasets, and how strong that relationship may be.
There can be positive correlation,negative correlation or zero correlation.The code and outputs for this is given in Correlation Analysis file.

4.MONTE CARLO ANALYSIS FOR OUTPUT PARAMETERS(Code is available in code.cpp)
After generating records for input parameters,visualising the input parameters,checking statistical correlation of all combinations of input parameters,now Monte Carlo analysis is performed requires assigning multiple values to an uncertain variable to achieve multiple results and then averaging the results to obtain an estimate. So hence assigning values to all the input parameters and getting specific output parameters values for fairly large iterations the probability of various outcomes and more importantly quantum efficiency can be observed.

5.Output plots
Visualize the output plots as given in output file.











